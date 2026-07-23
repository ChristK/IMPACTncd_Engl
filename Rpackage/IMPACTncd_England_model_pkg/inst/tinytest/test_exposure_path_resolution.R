# Tests for user-defined exposure table locations.
#
# Exposure `file_name` entries in exposure_definitions may be:
#   1. a bare name  -> stock ./inputs/exposure_distributions/<name>
#   2. a relative path with "/" -> resolved against the sim_design YAML dir
#   3. an absolute path -> used verbatim
# Design$resolve_exposure_path() implements this; Design$xps_table_path()
# routes the reads that bypass Exposure objects (PARF, smok_relapse) through it.

library(tinytest)
library(IMPACTncdEngland)

# A minimal, self-contained but valid sim_prm (avoids depending on any shipped
# config, which may lag behind the current required_params). A single
# dependency-free disease keeps reorder_diseases() from needing igraph.
cfg <- list(
  locality = "England",
  clusternumber = 1L,
  clusternumber_export = 1L,
  logs = FALSE,
  cols_for_output = c("year", "age", "sex"),
  strata_for_output = c("year", "age", "sex"),
  exposures_for_output = list("bmi"),
  exposure_definitions = list(
    list(
      name = "bmi", file_name = "bmi", var_name = "bmi",
      rank_var = "rank_bmi", distribution = "continuous"
    )
  ),
  n = 10L,
  num_chunks = 1L,
  init_year_long = 2013L,
  sim_horizon_max = 2030L,
  ageL = 30L,
  ageH = 90L,
  diseases = list(
    chd = list(
      name = "chd",
      meta = list(incidence = list(influenced_by_disease_name = list()))
    )
  ),
  maxlag = 5L,
  smoking_relapse_limit = 3L,
  stochastic = TRUE,
  kismet = TRUE,
  jumpiness = 1.0,
  statin_adherence = 0.9,
  bpmed_adherence = 0.9,
  export_xps = TRUE,
  simsmok_calibration = FALSE,
  output_dir = file.path(tempdir(), "out_xps_test"),
  synthpop_dir = file.path(tempdir(), "sp_xps_test"),
  validation = FALSE,
  iteration_n_max = 10L
)

# Place the YAML in a known temp "scenario" directory so we can assert that
# YAML-relative table paths resolve against it (not against getwd()).
scen_dir <- file.path(tempdir(), "scen_xps_test")
dir.create(scen_dir, showWarnings = FALSE, recursive = TRUE)
yaml_path <- file.path(scen_dir, "sim_design_test.yaml")
yaml::write_yaml(cfg, yaml_path)

# Build from a clean working dir so no ./inputs is present -> skeleton mode
# (no data load), keeping the test hermetic regardless of where it runs.
clean_wd <- file.path(tempdir(), "clean_wd_xps_test")
dir.create(clean_wd, showWarnings = FALSE, recursive = TRUE)
old_wd <- setwd(clean_wd)
design <- tryCatch(
  Design$new(yaml_path),
  finally = setwd(old_wd)
)
setwd(old_wd)

priv <- design$.__enclos_env__$private
resolve <- priv$resolve_exposure_path

# --- design_dir captured from the YAML's location -------------------------
expect_equal(
  normalizePath(priv$design_dir),
  normalizePath(scen_dir)
)

# --- Case 1: bare name -> stock inputs dir (relative, backward compatible) --
expect_equal(
  resolve("bmi"),
  file.path("./inputs/exposure_distributions", "bmi")
)

# --- Case 2: relative path with "/" -> resolved against the YAML dir --------
expect_equal(
  resolve("tables/bmi_local"),
  normalizePath(file.path(scen_dir, "tables/bmi_local"), mustWork = FALSE)
)
# a leading "./" also contains "/" and anchors at the YAML dir
expect_equal(
  resolve("./bmi_local"),
  normalizePath(file.path(scen_dir, "./bmi_local"), mustWork = FALSE)
)

# --- Case 3: absolute path -> verbatim -------------------------------------
expect_equal(resolve("/mnt/data/bmi"), "/mnt/data/bmi")

# --- xps_table_path falls back to the resolver when no Exposure matches -----
# (skeleton mode: self$exposures is NULL, so smok_relapse -> stock path)
expect_equal(
  design$xps_table_path("smok_relapse"),
  file.path("./inputs/exposure_distributions", "smok_relapse")
)

# --- Exposure honours a pre-resolved file_path verbatim --------------------
# Build a tiny on-disk parquet dataset so Exposure$new()'s column pre-read
# succeeds, then confirm the resolved path is used unchanged.
if (requireNamespace("arrow", quietly = TRUE)) {
  ds_dir <- file.path(scen_dir, "tbl_custom")
  arrow::write_dataset(
    data.frame(year = 1:3, mu = c(0.1, 0.2, 0.3)),
    ds_dir,
    format = "parquet"
  )
  e <- Exposure$new(
    name = "custom",
    file_name = "ignored_when_file_path_supplied",
    file_path = ds_dir,
    var_name = "custom",
    rank_var = "rank_custom",
    distribution = "continuous"
  )
  expect_equal(normalizePath(e$file_path), normalizePath(ds_dir))
  expect_true(all(c("year", "mu") %in% e$get_col_names()))

  # Without file_path, a bare file_name still combines with base_path.
  e2 <- Exposure$new(
    name = "custom2",
    file_name = ds_dir, # absolute -> used as-is by Exposure's own fallback
    var_name = "custom2",
    rank_var = "rank_custom2",
    distribution = "continuous"
  )
  expect_equal(normalizePath(e2$file_path), normalizePath(ds_dir))
}
