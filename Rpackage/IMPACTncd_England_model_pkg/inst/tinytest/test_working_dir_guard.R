# Tests: working-directory guard and PARF/RR registry anchoring (GitHub #11).
#
# The bug: the model resolves ./inputs and ./simulation against getwd() at
# CALL time, but Disease/ExposureEffect fix their caches at CONSTRUCTION.
# RStudio evaluates .Rmd chunks in the document's folder, so running a
# vignette in place changed the working directory between Simulation$new()
# and $run(). gen_parf_files() then looked for the PARF registry under the
# wrong folder (creating an empty one there), judged the cache invalid and
# del_parf_file() deleted the REAL PARF cache, before the run failed on a
# relative read.
#
# Fix: Design$initialize() records the root (getwd()); Disease/ExposureEffect
# build their registry paths from it (anchor), and methods that touch
# ./inputs or ./simulation stop if getwd() has changed (guard).
#
# tinytest runs each file with the working directory set to the test file's
# folder, so every root below is an explicit tempfile() directory.

library(tinytest)
library(IMPACTncdEngland)

same_dir <- IMPACTncdEngland:::.same_dir
assert_wd <- IMPACTncdEngland:::.assert_wd
old_wd <- getwd()

mk <- function(prefix) {
  d <- tempfile(prefix)
  dir.create(d, recursive = TRUE)
  normalizePath(d, winslash = "/")
}
root <- mk("wdguard_root_")  # the "repository root" the model is built in
elsewhere <- mk("wdguard_elsewhere_") # e.g. the vignettes/ folder

# ---- 1. .same_dir ------------------------------------------------------------
expect_true(same_dir(root, root))
expect_true(same_dir(paste0(root, "/"), root), info = "trailing slash")
expect_true(same_dir("/", "/"), info = "the filesystem root is kept as '/'")
expect_false(same_dir(root, elsewhere))
lnk <- file.path(tempdir(), paste0("wdguard_link_", basename(root)))
if (isTRUE(suppressWarnings(file.symlink(root, lnk)))) {
  expect_true(same_dir(lnk, root), info = "a symlink to the root is the root")
}

# ---- 2. .assert_wd -----------------------------------------------------------
setwd(root)
expect_true(assert_wd(root, "x"))
expect_true(assert_wd(NULL, "x"), info = "objects without a root are not checked")
setwd(elsewhere)
expect_error(assert_wd(root, "run()"), pattern = "working directory changed")
err <- tryCatch(assert_wd(root, "run()"), error = conditionMessage)
expect_true(grepl("^run\\(\\):", err), info = "message names the operation")
expect_true(grepl(root, err, fixed = TRUE) && grepl(elsewhere, err, fixed = TRUE),
            info = "message names both folders")
expect_true(grepl(sprintf('setwd("%s")', root), err, fixed = TRUE),
            info = "message gives the remedy")

# ---- 3. a skeleton Design and Simulation built in `root` ---------------------
setwd(root)
dl <- list(
  locality = "England", clusternumber = 1, clusternumber_export = 1,
  logs = FALSE, cols_for_output = c("year", "age", "sex"),
  strata_for_output = c("year", "age", "sex"),
  exposures_for_output = list("bmi"),
  exposure_definitions = list(list(
    name = "bmi", file_name = "bmi", var_name = "bmi",
    rank_var = "rank_bmi", distribution = "continuous"
  )),
  n = 10, init_year_long = 2010, sim_horizon_max = 2030, ageL = 35, ageH = 80,
  diseases = list(chd = list(name = "chd", meta = list(
    incidence = list(influenced_by_disease_name = character(0))
  ))),
  maxlag = 5, smoking_relapse_limit = 3, stochastic = TRUE, kismet = 42,
  jumpiness = 0.2, statin_adherence = 0.85, bpmed_adherence = 0.9,
  export_xps = FALSE, export_PARF = FALSE, simsmok_calibration = FALSE,
  output_dir = file.path(root, "out"), synthpop_dir = file.path(root, "sp"),
  validation = FALSE, iteration_n_max = 2, num_chunks = 1
)
design <- suppressMessages(suppressWarnings(Design$new(dl)))
expect_true(same_dir(design$get_root_dir(), root),
            info = "Design records the working directory as the root")
expect_true(same_dir(design$clone(deep = TRUE)$get_root_dir(), root),
            info = "the root survives clone(deep = TRUE)")
sim <- suppressMessages(suppressWarnings(Simulation$new(design)))

# From the root nothing is refused...
expect_true(design$assert_working_dir("x"))

# ... from elsewhere every guarded method stops BEFORE doing anything.
setwd(elsewhere)
expect_error(Simulation$new(design), pattern = "working directory changed",
             info = "Simulation$new(<Design>) checks the Design's root")
expect_error(design$load_data(force = TRUE), pattern = "working directory changed")
expect_error(design$load_exposure_effects(), pattern = "working directory changed")
expect_error(design$load_exposures(), pattern = "working directory changed")
expect_error(design$load_diseases(), pattern = "working directory changed")
expect_error(sim$run(1L, multicore = FALSE, "sc0"), pattern = "^run\\(\\): working directory changed")
expect_error(sim$del_parfs(), pattern = "working directory changed")
expect_error(sim$del_RR_cache(), pattern = "working directory changed")
expect_error(sim$validate(), pattern = "working directory changed")
expect_error(sim$reconstruct_large_files(), pattern = "working directory changed")
expect_error(sim$del_large_files(), pattern = "working directory changed")
expect_error(sim$split_large_files(), pattern = "working directory changed")
expect_error(sim$zenodo_check_inputs(), pattern = "working directory changed")
expect_error(sim$zenodo_download_all(), pattern = "working directory changed")
expect_error(sim$calibrate_incd_ftlt(1L), pattern = "working directory changed")
# Methods that use only output_dir/synthpop_dir keep working from anywhere.
expect_silent(sim$del_summaries())
expect_silent(sim$del_logs())
expect_false(dir.exists(file.path(elsewhere, "simulation")),
             info = "nothing was created in the wrong folder")

# ---- 4. structure: the guard is the FIRST statement of each method -----------
guarded <- c("run", "validate", "del_parfs", "del_RR_cache", "split_large_files",
             "reconstruct_large_files", "del_large_files", "zenodo_connect",
             "zenodo_check_inputs", "zenodo_download_inputs",
             "zenodo_download_PARFs_RRs", "zenodo_download_all",
             "zenodo_upload_inputs", "zenodo_upload_simulation",
             "zenodo_resync_simulation", "zenodo_upload_all",
             "zenodo_create_archives", "zenodo_create_manifest",
             "calibrate_incd_ftlt")
for (m in guarded) {
  first <- body(Simulation$public_methods[[m]])[[2]]
  expect_true(
    identical(deparse(first[[1]]), "private$assert_wd"),
    info = paste0(m, "() must call private$assert_wd() first")
  )
}
for (m in c("load_data", "load_exposure_effects", "load_exposures", "load_diseases")) {
  expect_true(
    grepl("assert_working_dir", deparse(body(Design$public_methods[[m]])[[2]])[1], fixed = TRUE),
    info = paste0("Design$", m, "() must check the working directory first")
  )
}

# A scenario function that calls setwd() inside run() is invisible to a
# method-entry check: run_sim() re-checks right after each hook.
rs <- paste(deparse(body(Simulation$private_methods$run_sim)), collapse = "\n")
p1 <- regexpr("private$primary_prevention_scn(sp)", rs, fixed = TRUE)
c1 <- regexpr('private$assert_wd("the primary prevention scenario function")', rs, fixed = TRUE)
p2 <- regexpr("private$secondary_prevention_scn(sp)", rs, fixed = TRUE)
c2 <- regexpr('private$assert_wd("the secondary prevention scenario function")', rs, fixed = TRUE)
expect_true(p1 > 0 && c1 > p1, info = "re-check after the primary scenario hook")
expect_true(p2 > 0 && c2 > p2, info = "re-check after the secondary scenario hook")

# ---- 5. anchor: a REAL Disease keeps its cache from another folder -----------
setwd(root)
dir.create(file.path(root, "simulation"), showWarnings = FALSE)
dir.create(file.path(root, "inputs", "disease_burden"), recursive = TRUE, showWarnings = FALSE)
writeLines("x", file.path(root, "inputs", "disease_burden", "testdis_incd.fst"))
dis <- Disease$new(
  name = "testdis", friendly_name = "Test disease",
  meta = list(incidence = list(type = 2L), mortality = list(type = 1L)),
  design_ = design, RR = list()
)
p <- dis$.__enclos_env__$private
expect_true(same_dir(p$root_dir, root), info = "Disease takes the Design's root")
expect_true(startsWith(p$parf_dir, root))
dir.create(p$parf_filenam, recursive = TRUE)
writeLines("cached PARF", file.path(p$parf_filenam, "part-0.parquet"))
p$update_disease_registry(TRUE) # a valid record, as shipped in git

setwd(elsewhere) # e.g. RStudio running a vignette chunk
expect_true(startsWith(p$get_generation_record_path(), root),
            info = "the record path no longer follows getwd()")
expect_identical(
  p$to_rel_path(file.path(root, "inputs", "RR", "bmi~chd.csvy")),
  "inputs/RR/bmi~chd.csvy",
  info = "to_rel_path() strips the root, not getwd()"
)
# Exactly the first block of gen_parf_files():
p$UpdateDiseaseSnapshotIfInvalid(TRUE, function() dis$del_parf_file())
expect_true(dir.exists(p$parf_filenam),
            info = "the PARF cache survives a changed working directory")
expect_false(dir.exists(file.path(elsewhere, "simulation")),
             info = "no stray simulation/parf_generation_record.csv elsewhere")

# ---- 6. anchor: ExposureEffect's registry helpers use the root too -----------
ee_env <- new.env(parent = asNamespace("IMPACTncdEngland"))
ee_priv <- new.env(parent = emptyenv())
ee_priv$root_dir <- root
for (nm in c("root_or_wd", "get_generation_record_path", "to_rel_path")) {
  f <- ExposureEffect$private_methods[[nm]]
  environment(f) <- ee_env
  assign(nm, f, envir = ee_priv)
}
ee_env$private <- ee_priv
expect_true(startsWith(ee_priv$get_generation_record_path(), root))
expect_identical(ee_priv$to_rel_path(file.path(root, "simulation", "rr", "x.fst")),
                 "simulation/rr/x.fst")
ee_priv$root_dir <- NULL
expect_true(startsWith(ee_priv$get_generation_record_path(), getwd()),
            info = "an object without a root falls back to getwd()")

setwd(old_wd)
unlink(c(root, elsewhere, lnk), recursive = TRUE)
