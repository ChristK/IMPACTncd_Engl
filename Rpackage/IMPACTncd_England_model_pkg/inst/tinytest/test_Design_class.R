library(tinytest)
library(yaml)
library(IMPACTncdEngland)

# Build Designs from a clean working directory so no ./inputs is present: the
# constructor then stays in "skeleton mode" (no data load), keeping these
# sim_prm-processing tests hermetic and fast wherever the suite is launched.
clean_wd <- file.path(tempdir(), "test_design_clean_wd")
dir.create(clean_wd, showWarnings = FALSE, recursive = TRUE)
old_wd <- setwd(clean_wd)

# A complete, valid sim_prm list satisfying Design's required_params.
create_design_list <- function(overrides = list()) {
    base <- list(
        locality = "England",
        clusternumber = 1,
        clusternumber_export = 1,
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
        n = 10,
        init_year_long = 2010,
        sim_horizon_max = 2030,
        ageL = 35,
        ageH = 80,
        diseases = list(
            chd = list(name = "chd", meta = list(
                incidence = list(influenced_by_disease_name = character(0))
            ))
        ),
        maxlag = 5,
        smoking_relapse_limit = 3,
        stochastic = TRUE,
        kismet = 42,
        jumpiness = 0.2,
        statin_adherence = 0.85,
        bpmed_adherence = 0.9,
        export_xps = TRUE,
        simsmok_calibration = FALSE,
        output_dir = tempdir(),
        synthpop_dir = tempdir(),
        validation = FALSE,
        iteration_n_max = 2,
        num_chunks = 1
    )
    modifyList(base, overrides)
}

# ---- Initialization: derived time fields ----
dl <- create_design_list()
design <- Design$new(dl)

expect_true("init_year" %in% names(design$sim_prm))
expect_equal(design$sim_prm$init_year, 10) # 2010 - 2000
expect_equal(design$sim_prm$sim_horizon_max, 20) # 2030 - 2010

# ---- Default values set ----
expect_true(design$sim_prm$national_qimd)
expect_equal(design$sim_prm$init_year_fromGUI, 10)
expect_equal(design$sim_prm$sim_horizon_fromGUI, 20)

# ---- Directory paths made absolute (or Docker /outputs, /synthpop) ----
expect_true(grepl("^/", design$sim_prm$output_dir))
expect_true(grepl("^/", design$sim_prm$synthpop_dir))

# ---- Disease present ----
expect_true("chd" %in% names(design$sim_prm$diseases))

# ---- Topological disease ordering: a dependency precedes its dependent ----
dep <- create_design_list(list(diseases = list(
    b = list(name = "b", meta = list(
        incidence = list(influenced_by_disease_name = "a")
    )),
    a = list(name = "a", meta = list(
        incidence = list(influenced_by_disease_name = character(0))
    ))
)))
dep_design <- Design$new(dep)
ord <- names(dep_design$sim_prm$diseases)
expect_true(all(c("a", "b") %in% ord))
expect_true(which(ord == "a") < which(ord == "b"))

# ---- xps_table_path falls back to the stock path with no Exposure loaded ----
# (skeleton mode: self$exposures is empty, so smok_relapse -> stock inputs dir)
expect_equal(
    design$xps_table_path("smok_relapse"),
    file.path("./inputs/exposure_distributions", "smok_relapse")
)

# ---- The shipped default config is valid and loadable from a path ----
design_path <- system.file(
    "config/default_sim_design.yaml",
    package = "IMPACTncdEngland"
)
expect_true(nzchar(design_path))
design_yaml <- Design$new(design_path)
expect_true(inherits(design_yaml, "Design"))
# design_dir is captured from the YAML's own directory
expect_equal(
    normalizePath(design_yaml$.__enclos_env__$private$design_dir),
    normalizePath(dirname(design_path))
)
# init_year derived from the shipped init_year_long (2013)
expect_equal(design_yaml$sim_prm$init_year, 13)

# ---- Save to disk round trip ----
yaml_out <- tempfile(fileext = ".yaml")
design$save_to_disk(yaml_out)
roundtrip <- yaml::read_yaml(yaml_out)
expect_true("locality" %in% names(roundtrip))
expect_equal(roundtrip$locality, "England")
unlink(yaml_out)

setwd(old_wd)
