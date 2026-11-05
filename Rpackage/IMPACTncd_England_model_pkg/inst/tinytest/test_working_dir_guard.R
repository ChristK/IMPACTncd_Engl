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
# Every use of the new API goes through get0()/safe() below, so on a build
# WITHOUT the fix each check FAILS (with an explanation) instead of the file
# aborting at its first line -- the file is its own negative control.
#
# tinytest runs each file with the working directory set to the test file's
# folder, so every root below is an explicit tempfile() directory.

library(tinytest)
library(IMPACTncdEngland)

ns <- asNamespace("IMPACTncdEngland")
internal <- function(nm) get0(nm, envir = ns, inherits = FALSE)
canon_dir <- internal(".canon_dir")
same_dir <- internal(".same_dir")
assert_wd <- internal(".assert_wd")
# Evaluate, returning NA (with the error attached) instead of erroring.
safe <- function(expr) {
  tryCatch(expr, error = function(e) structure(NA, err = conditionMessage(e)))
}
errmsg <- function(expr) tryCatch({ force(expr); "<no error>" }, error = conditionMessage)
old_wd <- getwd()
# Simulation$new() may set data.table/fst/arrow to 1 thread (OpenMP BLAS
# guard): save the session state now, restore it at the end.
state0 <- list(dt = data.table::getDTthreads(), fst = fst::threads_fst(),
               arrow = arrow::cpu_count())

mk <- function(prefix) {
  d <- tempfile(prefix)
  dir.create(d, recursive = TRUE)
  normalizePath(d, winslash = "/")
}
root <- mk("wdguard_root_")           # the "repository root" the model is built in
elsewhere <- mk("wdguard_elsewhere_") # e.g. the vignettes/ folder

# ---- 1. path canonicalisation --------------------------------------------------
expect_true(is.function(canon_dir) && is.function(same_dir) && is.function(assert_wd),
            info = "the guard helpers exist")
if (is.function(canon_dir)) {
  # Non-existent paths pass through normalizePath() unchanged on Unix, so these
  # exercise canon_dir's own trailing-slash rule (not realpath()'s).
  expect_identical(canon_dir("/nonexistent_wdg/a//"), "/nonexistent_wdg/a",
                   info = "trailing slashes are dropped")
  expect_identical(canon_dir("/"), "/", info = "the filesystem root stays '/'")
  if (.Platform$OS.type != "windows") {
    expect_identical(canon_dir("C:/"), "C:/", info = "a drive root keeps its '/'")
  }
  expect_true(same_dir(root, root))
  expect_false(same_dir(root, elsewhere))
  lnk <- file.path(tempdir(), paste0("wdguard_link_", basename(root)))
  if (isTRUE(suppressWarnings(file.symlink(root, lnk)))) {
    expect_true(same_dir(lnk, root), info = "a symlink to the root is the root")
  }
}

# ---- 2. .assert_wd and its two messages ------------------------------------------
if (is.function(assert_wd)) {
  setwd(root)
  expect_true(assert_wd(root, "x"))
  expect_true(assert_wd(NULL, "x"), info = "objects without a root are not checked")
  setwd(elsewhere)
  err <- errmsg(assert_wd(root, "run()"))
  expect_true(grepl("^run\\(\\): working directory changed", err),
              info = "names the operation")
  expect_true(grepl(root, err, fixed = TRUE) && grepl(elsewhere, err, fixed = TRUE),
              info = "names both folders")
  expect_true(grepl(sprintf('setwd("%s")', root), err, fixed = TRUE) &&
                grepl("console", err, fixed = TRUE),
              info = "remedy: setwd() at the console")
  expect_true(grepl(".Rmd", err, fixed = TRUE) &&
                grepl("Evaluate chunks in directory", err, fixed = TRUE),
              info = "explains the RStudio .Rmd case, where a setwd() in a chunk is undone")
  err <- errmsg(assert_wd(root, "the primary prevention scenario function", "scenario"))
  expect_true(grepl("^the primary prevention scenario function changed the working directory", err),
              info = "scenario cause: says the scenario function moved it")
  expect_true(grepl("on.exit(setwd(old)", err, fixed = TRUE),
              info = "scenario cause: remedy is to restore it, not to setwd() and retry")
  setwd(root)
}

# ---- 3. a skeleton Design and Simulation built in `root` ------------------------
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
rd <- safe(design$get_root_dir())
expect_true(isTRUE(safe(same_dir(rd, root))),
            info = "Design records the working directory as the root")
expect_true(isTRUE(safe(same_dir(design$clone(deep = TRUE)$get_root_dir(), root))),
            info = "the root survives clone(deep = TRUE)")
sim <- suppressMessages(suppressWarnings(Simulation$new(design)))

# Relative output_dir / synthpop_dir are resolved against the root, so the
# unguarded del_*/export_* methods (which use only these) cannot act on
# another folder after a setwd(). normalizePath() leaves a not-yet-existing
# relative path unchanged, which is how './outputs' used to stay relative.
dl_rel <- modifyList(dl, list(output_dir = "out_rel", synthpop_dir = "sp_rel"))
design_rel <- suppressMessages(suppressWarnings(Design$new(dl_rel)))
expect_true(isTRUE(safe(same_dir(design_rel$sim_prm$output_dir, file.path(root, "out_rel")))),
            info = "a relative output_dir is anchored on the root")
expect_true(isTRUE(safe(same_dir(design_rel$sim_prm$synthpop_dir, file.path(root, "sp_rel")))),
            info = "a relative synthpop_dir is anchored on the root")

# An empty path, the root itself, or a folder containing it is refused:
# del_outputs()/del_synthpops() delete everything inside it. ("" used to
# resolve to the root, so del_outputs() would have emptied the checkout.)
# ".." after a folder that does not exist yet ("ghost/..") resolves to the
# root once the folder is created, so it must be refused up front too.
for (bad in list("", NA_character_, ".", root, dirname(root), "ghost/..",
                 "ghost/../..", file.path(root, "ghost2", ".."))) {
  m <- errmsg(suppressMessages(Design$new(modifyList(dl, list(output_dir = bad)))))
  expect_true(grepl("non-empty path|model's root folder or contains it", m),
              info = paste0("output_dir = ", deparse(bad), " must be refused; got: ", substr(m, 1, 100)))
  m <- errmsg(suppressMessages(Design$new(modifyList(dl, list(synthpop_dir = bad)))))
  expect_true(grepl("non-empty path|model's root folder or contains it", m),
              info = paste0("synthpop_dir = ", deparse(bad), " must be refused"))
}
# A folder inside the repository is fine (only the root and its ancestors are not).
expect_true(is.character(safe(suppressMessages(
  Design$new(modifyList(dl, list(output_dir = "outputs/run1"))))$sim_prm$output_dir)),
  info = "a sub-folder of the root is accepted")
# update_output_path()/update_synthpop_path() resolve the same way.
invisible(safe(suppressMessages(sim$update_output_path("rel_out"))))
expect_true(isTRUE(safe(same_dir(sim$design$sim_prm$output_dir, file.path(root, "rel_out")))),
            info = "update_output_path() anchors a relative path on the root")
expect_true(grepl("non-empty path|model's root folder", errmsg(sim$update_synthpop_path(""))),
            info = "update_synthpop_path('') is refused")
expect_true(grepl("model's root folder", errmsg(sim$update_output_path(root))),
            info = "update_output_path(<root>) is refused")
expect_true(grepl("model's root folder", errmsg(sim$update_output_path("ghost3/.."))),
            info = "update_output_path('ghost3/..') is refused")
expect_true(isTRUE(safe(same_dir(
  suppressMessages(Design$new(modifyList(dl, list(output_dir = "ghost/../outputs"))))$sim_prm$output_dir,
  file.path(root, "outputs")))), info = "'ghost/../outputs' resolves lexically to <root>/outputs")
invisible(safe(suppressMessages(sim$update_output_path(file.path(root, "out")))))
# carry_over_lifecourse_files_from is resolved against the root too, so one
# call from another folder cannot copy from that tree into this one.
for (b in c(root, elsewhere)) {
  dir.create(file.path(b, "prev_out", "lifecourse"), recursive = TRUE, showWarnings = FALSE)
}
writeLines("root", file.path(root, "prev_out", "lifecourse", "A_lifecourse.txt"))
writeLines("elsewhere", file.path(elsewhere, "prev_out", "lifecourse", "B_lifecourse.txt"))
setwd(elsewhere)
invisible(safe(suppressMessages(sim$update_output_path(
  file.path(root, "carry_dst"), carry_over_lifecourse_files_from = "prev_out"))))
setwd(root)
copied <- list.files(file.path(root, "carry_dst", "lifecourse"))
expect_true("A_lifecourse.txt" %in% copied && !"B_lifecourse.txt" %in% copied,
            info = "carry-over copies from <root>/prev_out, not from the current folder")
invisible(safe(suppressMessages(sim$update_output_path(file.path(root, "out")))))

expect_true(isTRUE(safe(design$assert_working_dir("x"))), info = "no refusal from the root")

# ---- 4. from elsewhere, every guarded method stops BEFORE doing anything ---------
refused <- function(expr, what) {
  m <- errmsg(expr)
  expect_true(grepl("working directory changed", m, fixed = TRUE),
              info = paste0(what, " must refuse from another folder; got: ", substr(m, 1, 120)))
}
setwd(elsewhere)
refused(Simulation$new(design), "Simulation$new(<Design>)")
refused(design$load_data(force = TRUE), "Design$load_data()")
refused(design$load_exposure_effects(), "Design$load_exposure_effects()")
refused(design$load_exposures(), "Design$load_exposures()")
refused(design$load_diseases(), "Design$load_diseases()")
refused(sim$run(1L, multicore = FALSE, "sc0"), "run()")
refused(sim$del_parfs(), "del_parfs()")
refused(sim$del_RR_cache(), "del_RR_cache()")
refused(sim$validate(), "validate()")
refused(sim$reconstruct_large_files(), "reconstruct_large_files()")
refused(sim$del_large_files(), "del_large_files()")
refused(sim$split_large_files(), "split_large_files()")
refused(sim$zenodo_check_inputs(), "zenodo_check_inputs()")
refused(sim$zenodo_download_all(), "zenodo_download_all()")
refused(sim$calibrate_incd_ftlt(1L), "calibrate_incd_ftlt()")
refused(SynthPop$new(1L, design), "SynthPop$new()")
expect_true(grepl("^run\\(\\): working directory changed", errmsg(sim$run(1L, FALSE, "sc0"))),
            info = "run() is refused by the guard itself, first")
# Methods that use only output_dir/synthpop_dir keep working from anywhere.
expect_silent(sim$del_summaries())
expect_silent(sim$del_logs())
expect_false(dir.exists(file.path(elsewhere, "simulation")),
             info = "nothing was created in the wrong folder")

# ---- 5. structure: the guard is the FIRST statement of each method ---------------
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
# SynthPop$new() and Simulation$new(<Design>) check too (behaviourally above).
expect_true(grepl("assert_working_dir", paste(deparse(body(SynthPop$public_methods$initialize)), collapse = " "), fixed = TRUE),
            info = "SynthPop$new() checks the working directory")
expect_true(grepl("assert_working_dir", paste(deparse(body(Simulation$public_methods$initialize)), collapse = " "), fixed = TRUE),
            info = "Simulation$new(<Design>) checks the working directory")
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
c1 <- regexpr("private$assert_wd(\"the primary prevention scenario function\"", rs, fixed = TRUE)
p2 <- regexpr("private$secondary_prevention_scn(sp)", rs, fixed = TRUE)
c2 <- regexpr("private$assert_wd(\"the secondary prevention scenario function\"", rs, fixed = TRUE)
expect_true(p1 > 0 && c1 > p1, info = "re-check after the primary scenario hook")
expect_true(p2 > 0 && c2 > p2, info = "re-check after the secondary scenario hook")

# ---- 6. anchor: a REAL Disease keeps its cache from another folder ---------------
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
expect_true(isTRUE(safe(same_dir(p$root_dir, root))), info = "Disease takes the Design's root")
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
# gen_parf_files() itself refuses: its PARF generation reads relative inputs.
# (Give the object one stub RR so it gets past the no-risk-factor early exit.)
# Make the snapshot STALE first, so a guard placed after the snapshot block
# would delete the cache -- otherwise '...deletes nothing first' cannot fail.
p$rr <- list(list(get_rr_file_path = function() NA_character_, name = "bmi",
                  lag = 1L, distribution = "lognormal"))
setwd(root)
unlink(p$get_generation_record_path())
setwd(elsewhere)
refused(dis$gen_parf_files(design, list(dis)), "Disease$gen_parf_files()")
expect_true(dir.exists(p$parf_filenam), info = "...and deletes nothing first")
# ... and structurally: the check comes before the snapshot/deletion block.
gpf <- vapply(as.list(body(Disease$public_methods$gen_parf_files)),
              function(x) paste(deparse(x), collapse = " "), character(1))
i_guard <- which(grepl(".assert_wd(", gpf, fixed = TRUE))[1]
i_snap <- which(grepl("UpdateDiseaseSnapshotIfInvalid", gpf, fixed = TRUE))[1]
expect_true(isTRUE(i_guard < i_snap),
            info = "gen_parf_files(): the working-directory check precedes the snapshot block")

# ---- 7. anchor: ExposureEffect's registry helpers use the root too ---------------
ee_env <- new.env(parent = ns)
ee_priv <- new.env(parent = emptyenv())
ee_priv$root_dir <- root
ee_ok <- TRUE
for (nm in c("root_or_wd", "get_generation_record_path", "to_rel_path")) {
  f <- ExposureEffect$private_methods[[nm]]
  if (is.null(f)) { ee_ok <- FALSE; next }
  environment(f) <- ee_env
  assign(nm, f, envir = ee_priv)
}
ee_env$private <- ee_priv
expect_true(ee_ok, info = "ExposureEffect has the anchored registry helpers")
if (ee_ok) {
  expect_true(startsWith(ee_priv$get_generation_record_path(), root))
  expect_identical(ee_priv$to_rel_path(file.path(root, "simulation", "rr", "x.fst")),
                   "simulation/rr/x.fst")
  ee_priv$root_dir <- NULL
  expect_true(startsWith(ee_priv$get_generation_record_path(), getwd()),
              info = "an object without a root falls back to getwd()")
}

setwd(old_wd)
unlink(c(root, elsewhere), recursive = TRUE)
if (exists("lnk")) unlink(lnk)
data.table::setDTthreads(state0$dt)
fst::threads_fst(state0$fst)
arrow::set_cpu_count(state0$arrow)
