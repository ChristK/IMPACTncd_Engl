#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# verify_init_prvl_crn.R
#
# End-to-end check, on the REAL model, that a scenario acting only from a
# FUTURE year leaves the simulated past untouched -- GitHub issue #8
# ("Initial prevalence differs between baseline and _future scenarios before
# sc_year").
#
# The model relies on common random numbers (CRN): a policy scenario and the
# baseline share every random draw that the policy does not causally touch,
# so that a scenario-minus-baseline contrast carries no Monte Carlo noise. A
# scenario that changes nothing before `sc_year` must therefore reproduce the
# baseline lifecourse EXACTLY, simulant by simulant, in every year before
# `sc_year`. Issue #8 reports that it does not: the number of prevalent cases
# per stratum is preserved, but WHICH simulants are prevalent at init_year
# changes.
#
# Five scenarios are run on ONE shared synthetic population. They act from
# `sc_year`, which is the design's init_year + 2 (2015 for the stock
# init_year_long of 2013, so that two simulated years precede it) unless
# --sc-year overrides it:
#
#   sc0               the baseline.
#   noop_future       a primary-prevention scenario that does nothing.
#                     NEGATIVE CONTROL: it must equal sc0 in every row and
#                     column of every year. If it does not, CRN is broken
#                     independently of #8.
#   xps_future        from sc_year, lower BMI, SBP and total cholesterol in a
#                     RANDOM HALF of the simulant-years (analogue of the
#                     reporter's "Modelled" arm; parf_deprivation =
#                     "original", the default). The half is drawn inside the
#                     scenario with rbinom() and NO seed of its own: run_sim()
#                     seeds every scenario from its name and restores the
#                     base-R and dqrng states afterwards, so this arm also
#                     checks that random draws made BY a scenario cannot
#                     shift any draw it shares with sc0.
#   dimd_future       from sc_year, move everyone to the least deprived decile
#                     and quintile (analogue of the reporter's "Non_Modelled"
#                     arm; parf_deprivation = "scenario": the cached PARF
#                     table is joined at the forced dimd).
#   dimd_future_orig  the same forced deprivation as dimd_future, but with the
#                     default parf_deprivation = "original": the PARF join and
#                     the initial-prevalence look-up use the pre-scenario
#                     dimd, which with_original_dimd() swaps in and back out
#                     again -- a code path neither other arm takes.
#
# Each non-baseline lifecourse is joined to sc0's on (pid, year) and compared,
# separately for rows before sc_year ("pre") and from sc_year on ("post"), on
# EVERY column the two lifecourses share (outcomes, weights, strata, markers),
# with per-column differing counts and the earliest differing year. The checks
# below decide the exit status (0 only if every check passes):
#
#   A. noop_future == sc0 everywhere: every column, every year, row set
#      included (CRN negative control).
#   B. xps_future, dimd_future and dimd_future_orig == sc0 in every column of
#      every row before sc_year, row set included. This is the property issue
#      #8 says is violated, so on UNFIXED code these checks are EXPECTED to
#      fail.
#   C. the same three DIFFER from sc0 from sc_year on (positive control: the
#      scenarios genuinely act, so B cannot pass vacuously). C looks only at
#      the OUTCOME columns (*_prvl, *_dgns, cms*, all_cause_mrtl) and at the
#      row set (a death that moved), not at every column: its job is to show
#      that the simulated health outcomes respond. A column that can change
#      without any outcome moving -- dimd, were the restore of the original
#      deprivation to break (D's job), or wt_esp, which only follows the row
#      set -- must not be able to satisfy it.
#   D. age, sex and dimd agree on every matched row. Outputs are stratified by
#      the ORIGINAL deprivation, so the dimd scenarios must not leak into them.
#   E. no PARF table was regenerated during the run (regeneration takes hours
#      and would mean the cached tables were not the ones used).
#
# Usage, from the repository root:
#   Rscript auxil/verify_init_prvl_crn.R <label> --base-dir=<path>
#           [--n=20000] [--sc-year=<short year>] [--rebuild] [--compare-only]
#
#   <label>            tags the output and results directories, e.g. "before"
#                      and "after" a fix. Letters, digits, '_' and '-' only.
#   --base-dir=<path>  the output root (see below). If omitted, the
#                      environment variable IMPACTNCD_VERIFY_DIR is used; one
#                      of the two is required, and the root's parent
#                      directory must already exist.
#   --n=<int>          population size of each synthpop (default 20000). The
#                      synthpop is cached per design checksum, so "before" and
#                      "after" runs with the same n reuse the same population.
#   --sc-year=<yy>     the year the scenarios act from, as a SHORT year (15 =
#                      2015). Default: the design's init_year + 2. At least one
#                      simulated year must precede it, and it must not lie
#                      beyond sim_horizon_max. With --compare-only the default
#                      is the sc_year recorded in results_<label>/
#                      provenance.yaml, and an override contradicting it is
#                      refused (the scenarios were simulated with that year).
#   --rebuild          source global.R first. global.R rebuilds IMPACTncdEngland
#                      from the working tree and reinstalls it into the shared
#                      user library (R_LIBS_USER) whenever the source changed
#                      since the last install -- which swaps the package under
#                      any other simulation running from that library, hence
#                      opt-in. Without it the INSTALLED build is loaded as is
#                      (a warning is printed if the source has changed since
#                      global.R last installed it). The mode used, and the md5
#                      of the installed .rdb, are recorded in provenance.yaml.
#   --compare-only     skip the simulation; redo only the comparison on an
#                      existing output_<label> tree. An arm that tree lacks
#                      (e.g. dimd_future_orig, in a tree written before that
#                      arm existed) is skipped with a note, not failed.
#
# Everything is written under the output root. The script refuses to create or
# delete anything outside it, and a full run deletes only inside a root that
# holds the marker file `.verify_init_prvl_crn`, which the script writes when
# it creates the root or finds it empty. (To reuse a non-empty root written by
# an earlier version of this script, create that file in it by hand.)
#   output_<label>/   model output (lifecourse, logs); recreated on each run.
#   synthpop/         SHARED across labels, so every label simulates the
#                     identical synthetic population.
#   results_<label>/  comparison CSVs, the design YAML used, provenance.
#
# The design is inputs/sim_design.yaml with a handful of overrides (see
# DESIGN_OVERRIDES); inputs/sim_design.yaml itself is never edited. The fields
# the cached PARF checksum depends on (init_year_long, sim_horizon_max, ageL,
# ageH, apply_RR_to_mrtl2, model_trends_in_residual_incd, the RR files) are
# deliberately left alone, and a pre-flight check stops the script before any
# simulation if a PARF table would have to be regenerated.
#
# Uses at most 2 cores: 2 forked workers, each with BLAS/OpenMP pinned to one
# thread (the script re-launches itself with OMP_NUM_THREADS=1 etc.). Runtime
# at n = 20000 with 8 CPUs available: about 1.2 min per scenario (both mc
# iterations together) when measured in October 2026; see timings.csv in
# results_<label>/ for each run's own figures.
#
# Examples (the output root is illustrative):
#   Rscript auxil/verify_init_prvl_crn.R after \
#     --base-dir=/mnt/storage_fast/verify_issue8
#   IMPACTNCD_VERIFY_DIR=/mnt/storage_fast/verify_issue8 \
#     Rscript auxil/verify_init_prvl_crn.R after --rebuild
#   Rscript auxil/verify_init_prvl_crn.R before --compare-only \
#     --base-dir=/mnt/storage_fast/verify_issue8
# ---------------------------------------------------------------------------

# ===========================================================================
# 0. Arguments, constants and guards
# ===========================================================================
args <- commandArgs(trailingOnly = TRUE)
usage <- paste0(
  "Usage: Rscript auxil/verify_init_prvl_crn.R <label> --base-dir=<path>\n",
  "         [--n=20000] [--sc-year=<short year>] [--rebuild] ",
  "[--compare-only]\n",
  "  The output root may instead be given by the environment variable\n",
  "  IMPACTNCD_VERIFY_DIR; its parent directory must already exist.\n",
  "  See the header of auxil/verify_init_prvl_crn.R for details."
)
pos_args <- args[!startsWith(args, "--")]
opt_args <- args[startsWith(args, "--")]
if (length(pos_args) != 1L || !grepl("^[A-Za-z0-9_-]+$", pos_args[[1L]])) {
  stop(usage, call. = FALSE)
}
bad_opts <- opt_args[!grepl(paste0(
  "^--(compare-only|rebuild)$|^--(n|sc-year)=[0-9]+$|^--base-dir=.+$"
), opt_args)]
if (length(bad_opts) > 0L) {
  stop("Unknown option(s): ", paste(bad_opts, collapse = " "), "\n", usage,
       call. = FALSE)
}
opt_names <- sub("=.*$", "", opt_args)
if (anyDuplicated(opt_names)) {
  stop("Option(s) given more than once: ",
       paste(unique(opt_names[duplicated(opt_names)]), collapse = " "), "\n",
       usage, call. = FALSE)
}
# Value of `--<name>=<value>`, or NULL when the option is absent.
opt_value <- function(name) {
  hit <- grep(paste0("^--", name, "="), opt_args, value = TRUE)
  if (length(hit)) sub(paste0("^--", name, "="), "", hit[[1L]]) else NULL
}
if (all(c("--rebuild", "--compare-only") %in% opt_args)) {
  stop("--rebuild and --compare-only exclude each other: a comparison-only ",
       "run never loads the model package.\n", usage, call. = FALSE)
}

# Output root: --base-dir, else IMPACTNCD_VERIFY_DIR. Checked here, before
# anything is created and before the re-launch below, so that a missing or
# misspelt root fails at once with the usage message rather than with an
# error from deep inside normalizePath(). An absolute path is kept exactly as
# given, because it is written into the design YAML (output_dir,
# synthpop_dir); a relative one is anchored at the working directory.
base_arg <- opt_value("base-dir")
if (is.null(base_arg)) base_arg <- Sys.getenv("IMPACTNCD_VERIFY_DIR")
if (!nzchar(base_arg)) {
  stop("No output root given: pass --base-dir=<path> or set the environment ",
       "variable IMPACTNCD_VERIFY_DIR.\n", usage, call. = FALSE)
}
base_arg <- sub("/+$", "", path.expand(base_arg))
if (!nzchar(base_arg) || basename(base_arg) %in% c(".", "..") ||
    !dir.exists(dirname(base_arg)) ||
    (file.exists(base_arg) && !dir.exists(base_arg))) {
  stop("Unusable output root '", base_arg, "': its parent directory must ",
       "exist, and the root itself must be a directory (or not exist yet), ",
       "not '/', '.' or '..'.\n", usage, call. = FALSE)
}
# Pin BLAS and OpenMP to ONE thread per process. Left unset, each forked
# worker's synthpop generation ran ~7 threads (measured ~690% CPU per worker
# on a machine with 8 CPUs available), so `clusternumber: 2` alone does not
# cap the run at 2 cores. These libraries read the variables when they are
# loaded, i.e. before this line runs, so the script re-launches itself with
# them set.
THREAD_VARS <- c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
                 MKL_NUM_THREADS = "1", FLEXIBLAS_NUM_THREADS = "1")
if (!all(Sys.getenv(names(THREAD_VARS)) == THREAD_VARS)) {
  self_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(self_path) != 1L) {
    stop("Run this file with Rscript, or set ",
         paste0(names(THREAD_VARS), "=1", collapse = " "), " yourself.",
         call. = FALSE)
  }
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(self_path), shQuote(args)),
    env = paste0(names(THREAD_VARS), "=", THREAD_VARS)
  )
  quit(save = "no", status = status)
}

LABEL        <- pos_args[[1L]]
COMPARE_ONLY <- "--compare-only" %in% opt_args
REBUILD      <- "--rebuild" %in% opt_args
n_opt        <- opt_value("n")
N_SIM        <- if (!is.null(n_opt)) as.integer(n_opt) else 20000L
stopifnot(!is.na(N_SIM), N_SIM > 0L)

REPO <- normalizePath(getwd())
if (!file.exists(file.path(REPO, "global.R")) ||
    !dir.exists(file.path(REPO, "Rpackage", "IMPACTncd_England_model_pkg"))) {
  stop("Run this script from the IMPACTncd_Engl repository root.", call. = FALSE)
}
PKG_SRC <- file.path(REPO, "Rpackage", "IMPACTncd_England_model_pkg")

BASE_DIR <- if (startsWith(base_arg, "/")) base_arg else
  file.path(normalizePath(dirname(base_arg), mustWork = TRUE),
            basename(base_arg))
OUT_DIR  <- file.path(BASE_DIR, paste0("output_", LABEL))
SP_DIR   <- file.path(BASE_DIR, "synthpop")
RES_DIR  <- file.path(BASE_DIR, paste0("results_", LABEL))
# A full run deletes inside BASE_DIR only when this file is present (see the
# header): it marks a root that this script created, or found empty.
MARKER   <- file.path(BASE_DIR, ".verify_init_prvl_crn")

MC          <- 1:2  # Monte Carlo iterations (integer, as run() requires)
SCENARIOS   <- c("noop_future", "xps_future", "dimd_future", "dimd_future_orig")
ACTIVE_SCN  <- c("xps_future", "dimd_future", "dimd_future_orig") # they act
PARF_DEPRIVATION <- c(
  noop_future = "original", xps_future = "original",
  dimd_future = "scenario", dimd_future_orig = "original"
)
XPS_SHARE   <- 0.5  # chance that xps_future changes a row from sc_year on

# Lifecourse columns that make up the OUTCOMES (check C), and columns that
# must agree on matched rows (check D). Checks A and B compare every column.
CMP_PATTERN <- "_prvl$|_dgns$|^all_cause_mrtl$|^cms"
SANITY_COLS <- c("age", "sex", "dimd")
# Strata of the initial-prevalence tables (inputs/disease_burden/*_prvl),
# used only by the init-year diagnostic table.
PRVL_STRATA <- c("year", "age", "sex", "dimd", "ethnicity", "sha")

# Design overrides. Nothing here may touch a field that enters the PARF
# checksum (see Disease_class.R, private$chksum), or the cached PARF tables
# in simulation/parf would no longer match and would be regenerated.
DESIGN_OVERRIDES <- list(
  n                    = N_SIM,
  num_chunks           = 1L,
  clusternumber        = 2L,
  clusternumber_export = 1L,
  logs                 = TRUE,
  keep_lifecourse      = TRUE, # informational; the lifecourse is always written
  export_xps           = FALSE,
  output_dir           = OUT_DIR,
  synthpop_dir         = SP_DIR
)
PARF_CHECKSUM_FIELDS <- c(
  "init_year_long", "sim_horizon_max", "ageL", "ageH",
  "apply_RR_to_mrtl2", "model_trends_in_residual_incd"
)
stopifnot(!any(names(DESIGN_OVERRIDES) %in% PARF_CHECKSUM_FIELDS))

# ---- the year the scenarios act from (sp$pop and the lifecourse store
# short years: 15 = 2015). Default: the design's init_year + 2, so that two
# simulated years precede it. A comparison must use the year the scenarios
# were SIMULATED with, so --compare-only reads it from the run's
# provenance.yaml (and the run's own design) whenever they exist.
design_years_path <- file.path(RES_DIR, "sim_design.yaml")
if (!COMPARE_ONLY || !file.exists(design_years_path)) {
  design_years_path <- file.path(REPO, "inputs", "sim_design.yaml")
}
design_years <- yaml::read_yaml(design_years_path)
INIT_YEAR_SHORT <- as.integer(design_years$init_year_long) - 2000L
LAST_YEAR_SHORT <- as.integer(design_years$sim_horizon_max) - 2000L
if (length(INIT_YEAR_SHORT) != 1L || length(LAST_YEAR_SHORT) != 1L ||
    anyNA(c(INIT_YEAR_SHORT, LAST_YEAR_SHORT))) {
  stop("Cannot read init_year_long and sim_horizon_max from ",
       design_years_path, call. = FALSE)
}
prov_sc_year <- NULL
prov_path <- file.path(RES_DIR, "provenance.yaml")
if (COMPARE_ONLY && file.exists(prov_path)) {
  prov_sc_year <- as.integer(yaml::read_yaml(prov_path)$sc_year_short)
  if (length(prov_sc_year) != 1L || is.na(prov_sc_year)) prov_sc_year <- NULL
}
sc_opt <- opt_value("sc-year")
if (!is.null(sc_opt)) {
  sc_opt <- as.integer(sc_opt)
  if (sc_opt >= 100L) {
    stop("--sc-year takes a SHORT year (e.g. 15 for 2015), not ", sc_opt, ".",
         call. = FALSE)
  }
  if (!is.null(prov_sc_year) && sc_opt != prov_sc_year) {
    stop(sprintf(paste0(
      "--sc-year=%d contradicts %s, which records that the scenarios in %s ",
      "were simulated with sc_year %d. Comparing at another year would ",
      "report differences the scenarios were designed to produce (or miss ",
      "ones they should not). Drop --sc-year, or rerun the simulation."),
      sc_opt, prov_path, OUT_DIR, 2000L + prov_sc_year), call. = FALSE)
  }
}
if (!is.null(sc_opt)) {
  SC_YEAR_SHORT  <- sc_opt
  SC_YEAR_SOURCE <- "--sc-year"
} else if (!is.null(prov_sc_year)) {
  SC_YEAR_SHORT  <- prov_sc_year
  SC_YEAR_SOURCE <- paste("recorded in", prov_path)
} else {
  SC_YEAR_SHORT  <- INIT_YEAR_SHORT + 2L
  SC_YEAR_SOURCE <- paste("init_year + 2, from", design_years_path)
}
if (SC_YEAR_SHORT <= INIT_YEAR_SHORT) {
  stop(sprintf(paste0(
    "sc_year %d leaves no simulated year before it: the simulation starts ",
    "in %d (init_year_long in %s), so there is nothing for checks A-B to ",
    "compare. Pass --sc-year=<short year> of at least %d."),
    2000L + SC_YEAR_SHORT, 2000L + INIT_YEAR_SHORT, design_years_path,
    INIT_YEAR_SHORT + 1L), call. = FALSE)
}
if (SC_YEAR_SHORT > LAST_YEAR_SHORT) {
  stop(sprintf(paste0(
    "sc_year %d lies beyond the last simulated year %d (sim_horizon_max in ",
    "%s): the scenarios would never act and check C could not pass. Pass ",
    "--sc-year=<short year> of at most %d."),
    2000L + SC_YEAR_SHORT, 2000L + LAST_YEAR_SHORT, design_years_path,
    LAST_YEAR_SHORT), call. = FALSE)
}

# ---- output root. A comparison only reads output_<label> and (re)writes
# files in results_<label>, so it needs an existing tree and nothing else. A
# full run deletes output_<label> and results_<label> (sections 1 and 2), so
# it insists on the marker file: a root it creates, or finds empty, is
# marked; any other root is refused, so that a --base-dir pointed at a shared
# directory cannot cost someone else's output_<label>.
if (COMPARE_ONLY) {
  if (!dir.exists(OUT_DIR)) {
    stop("--compare-only: no output tree at ", OUT_DIR, ". Run without ",
         "--compare-only first, or check --base-dir and the label.",
         call. = FALSE)
  }
} else {
  base_fresh <- !dir.exists(BASE_DIR) ||
    length(list.files(BASE_DIR, all.files = TRUE, no.. = TRUE)) == 0L
  dir.create(BASE_DIR, showWarnings = FALSE) # its parent exists (checked)
  if (base_fresh) {
    writeLines(c(
      "Output root of auxil/verify_init_prvl_crn.R (IMPACTncd_Engl).",
      "Full runs delete and recreate output_<label>/ and results_<label>/ here."
    ), MARKER)
  }
  if (!file.exists(MARKER)) {
    stop("Refusing a full run in ", BASE_DIR, ": it is not empty and does ",
         "not hold the marker file ", basename(MARKER), ", so it was not ",
         "created by this script, which would delete ", basename(OUT_DIR),
         "/ and ", basename(RES_DIR), "/ in it. Choose a new or empty ",
         "directory or, if this root WAS written by an earlier version of ",
         "this script, create ", MARKER, " to adopt it.", call. = FALSE)
  }
}

# Refuse to create or delete anything outside BASE_DIR. Both sides are
# compared fully resolved, because /mnt/storage_fast may be a symlink (here it
# resolves to /data); the parent is resolved rather than the path itself so
# that the answer does not depend on whether the path exists yet.
assert_inside_base <- function(path) {
  base_real <- normalizePath(BASE_DIR, mustWork = TRUE)
  p <- file.path(normalizePath(dirname(path), mustWork = FALSE), basename(path))
  if (!startsWith(p, paste0(base_real, "/"))) {
    stop("Refusing to touch a path outside ", BASE_DIR, ": ", p, call. = FALSE)
  }
  invisible(p)
}

# PASS/FAIL reporter, as in auxil/verify_two_agegrps_on_real_output.R.
ok <- TRUE
check_log <- character()
say <- function(pass, msg) {
  ok <<- ok && pass
  line <- sprintf("[%s] %s", if (pass) " OK " else "FAIL", msg)
  check_log <<- c(check_log, line)
  cat(line, "\n", sep = "")
}
# A check that could not be run: logged, but neither passes nor fails.
skip <- function(msg) {
  line <- sprintf("[SKIP] %s", msg)
  check_log <<- c(check_log, line)
  cat(line, "\n", sep = "")
}

# Phase timings.
timings <- list()
timed <- function(phase, expr) {
  t0 <- Sys.time()
  force(expr)
  t1 <- Sys.time()
  timings[[length(timings) + 1L]] <<- list(
    phase = phase,
    start = format(t0, "%Y-%m-%d %H:%M:%S"),
    end = format(t1, "%Y-%m-%d %H:%M:%S"),
    minutes = round(as.numeric(difftime(t1, t0, units = "mins")), 2)
  )
  cat(sprintf("  ... %s took %.1f min\n", phase,
              as.numeric(difftime(t1, t0, units = "mins"))))
  invisible(NULL)
}

PKG_MODE <- if (COMPARE_ONLY) {
  "not loaded (--compare-only)"
} else if (REBUILD) {
  "--rebuild: sourced global.R (rebuilds and reinstalls if the source changed)"
} else {
  "installed build, loaded as is (no rebuild; pass --rebuild to rebuild)"
}
T_START <- Sys.time()
cat(sprintf(
  "verify_init_prvl_crn.R  label = %s  %smc = %s  sc_year = %d  %s\n",
  LABEL, if (COMPARE_ONLY) "" else sprintf("n = %d  ", N_SIM),
  paste(MC, collapse = ","), 2000L + SC_YEAR_SHORT,
  if (COMPARE_ONLY) "(compare only)" else ""
))
cat("  output root: ", BASE_DIR, "\n",
    "  sc_year:     ", SC_YEAR_SOURCE, "\n",
    "  package:     ", PKG_MODE, "\n", sep = "")

# ===========================================================================
# 1. Environment
# ===========================================================================
# TRUE when the package source is unchanged since global.R last installed it
# (the snapshot that global.R's installLocalPackageIfChanged() compares
# against, re-checked here read-only), FALSE when it has changed, NA when no
# snapshot exists.
source_matches_install <- function() {
  snap <- file.path(REPO, "Rpackage",
                    ".IMPACTncd_England_model_pkg_snapshot.rds")
  if (!file.exists(snap)) return(NA)
  ch <- utils::changedFiles(readRDS(snap))
  !any(nzchar(unlist(ch[c("added", "deleted", "changed")])))
}
SOURCE_MATCHES_INSTALL <- NA
if (!COMPARE_ONLY && REBUILD) {
  # Initialises the environment and rebuilds the local package if its source
  # changed since the last install, so the run uses the CURRENT source. This
  # reinstalls into the shared user library, hence opt-in (see the header).
  timed("source(global.R)", source(file.path(REPO, "global.R")))
  SOURCE_MATCHES_INSTALL <- source_matches_install()
} else if (!COMPARE_ONLY) {
  # Only what global.R does that a RUN needs, without its rebuild/reinstall:
  #  * attach IMPACTncdEngland, and gamlss.dist: the design names exposure
  #    quantile functions (e.g. qZABNB) that are resolved by name at run
  #    time, so gamlss.dist must be on the search path;
  #  * single-threaded data.table, fst and arrow, not restored after a fork
  #    (data.table and arrow are re-set to 2 threads below, as in --rebuild);
  #  * an initialised .Random.seed: run_sim() saves and restores it around
  #    every scenario and fails if it does not exist yet.
  # The rest of global.R serves installation (CRAN mirror, download timeout,
  # creating the user library, installing missing dependencies, the rebuild
  # itself) or sets what is already the default under Rscript
  # (datatable.verbose and datatable.showProgress FALSE; rgl is not used).
  # The working directory is the repository root in both modes (checked
  # above), as global.R and the design's relative ./inputs paths assume.
  if (!nzchar(system.file(package = "IMPACTncdEngland"))) {
    stop("IMPACTncdEngland is not installed in ",
         paste(.libPaths(), collapse = ", "),
         ". Run once with --rebuild to build it from the working tree.",
         call. = FALSE)
  }
  timed("library(IMPACTncdEngland)", suppressPackageStartupMessages({
    library(IMPACTncdEngland)
    library(gamlss.dist)
  }))
  data.table::setDTthreads(threads = 1L, restore_after_fork = FALSE)
  fst::threads_fst(nr_of_threads = 1L, reset_after_fork = FALSE)
  arrow::set_cpu_count(1L)
  invisible(runif(1))
  SOURCE_MATCHES_INSTALL <- source_matches_install()
  if (isFALSE(SOURCE_MATCHES_INSTALL)) {
    message(
      "WARNING: the package source has changed since global.R last ",
      "installed IMPACTncdEngland, so this run verifies the INSTALLED build, ",
      "not the working tree. Rerun with --rebuild to test the current source."
    )
  } else if (is.na(SOURCE_MATCHES_INSTALL)) {
    message("NOTE: no global.R install snapshot found, so whether the ",
            "installed build matches the working tree is unknown.")
  }
}
if (!COMPARE_ONLY) {
  cat("  installed:   ", system.file(package = "IMPACTncdEngland"), " (",
      as.character(packageVersion("IMPACTncdEngland")), "); source ",
      if (isTRUE(SOURCE_MATCHES_INSTALL)) "unchanged since"
      else if (isFALSE(SOURCE_MATCHES_INSTALL)) "CHANGED since"
      else "unknown relative to",
      " the last global.R install\n", sep = "")
}
suppressPackageStartupMessages({
  library(data.table)
  library(CKutils)
  library(yaml)
})
setDTthreads(2L)
arrow::set_cpu_count(2L)

assert_inside_base(RES_DIR)
if (!COMPARE_ONLY && dir.exists(RES_DIR)) unlink(RES_DIR, recursive = TRUE)
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# ===========================================================================
# 2. Simulation (skipped with --compare-only)
# ===========================================================================

# Scenario functions are built with bquote() so that sc_year is a LITERAL in
# the function body: update_primary_prevention_scn() re-parents the function
# onto the Simulation object's environment, so it must not rely on variables
# from this script. data.table::set() keeps them free of any data.table
# NSE/awareness concerns in that environment.
mk_noop_future <- function() {
  function(sp) NULL
}

mk_xps_future <- function(sc_year, share) {
  eval(bquote(function(sp) {
    sc_year <- .(sc_year)
    pop <- sp$pop
    # Only a random share of the rows from sc_year on, drawn with NO seed of
    # our own: run_sim() seeds each primary-prevention scenario from its name
    # (and mc_aggr) and afterwards restores both the base-R and the dqrng
    # state, so these draws exercise a scenario's own RNG use without being
    # able to shift any draw it shares with sc0. One draw per row of sp$pop,
    # as `sp$pop[year >= sc_year & rbinom(.N, 1L, share) == 1L, ...]` would
    # make, and only then restricted to year >= sc_year.
    hit <- stats::rbinom(nrow(pop), 1L, .(share)) == 1L
    idx <- which(pop[["year"]] >= sc_year & hit)
    # Lower three exposures that feed many diseases' relative risks, clamped
    # to the minimum values declared in exposure_definitions.
    shifts <- list(
      bmi_curr_xps   = function(x) pmax(x * 0.9, 14),
      sbp_curr_xps   = function(x) pmax(x - 10, 75),
      tchol_curr_xps = function(x) pmax(x * 0.9, 0.5)
    )
    miss <- setdiff(names(shifts), names(pop))
    if (length(miss) > 0L) {
      stop("xps_future: column(s) missing from sp$pop: ",
           paste(miss, collapse = ", "))
    }
    for (col in names(shifts)) {
      v <- shifts[[col]](pop[[col]][idx])
      if (is.integer(pop[[col]])) v <- as.integer(round(v)) # keep the type
      data.table::set(pop, idx, col, v)
    }
    invisible(NULL)
  }))
}

mk_dimd_future <- function(sc_year) {
  eval(bquote(function(sp) {
    sc_year <- .(sc_year)
    pop <- sp$pop
    idx <- which(pop[["year"]] >= sc_year)
    targets <- c(dimd = "10 least deprived", qimd = "5 least deprived")
    if (!"dimd" %in% names(pop)) stop("dimd_future: sp$pop has no dimd column")
    # Force dimd (and qimd, consistently, when present), preserving the
    # factor levels so the column type is unchanged.
    for (col in intersect(names(targets), names(pop))) {
      old <- pop[[col]]
      lvl <- targets[[col]]
      if (is.factor(old)) {
        if (!lvl %in% levels(old)) {
          stop("dimd_future: level '", lvl, "' not among the levels of ", col)
        }
        new <- factor(rep(lvl, length(idx)), levels = levels(old))
      } else {
        new <- rep(lvl, length(idx))
      }
      data.table::set(pop, idx, col, new)
    }
    invisible(NULL)
  }))
}

# Pre-flight: every disease that needs a cached PARF table must find one whose
# checksum matches AND whose generation record is valid; otherwise run() would
# start regenerating PARFs, which takes hours. Mirrors the conditions in
# Disease$gen_parf_files() and private$UpdateDiseaseSnapshotIfInvalid().
preflight_parf <- function(sim) {
  out <- rbindlist(lapply(sim$design$diseases, function(d) {
    priv <- d$.__enclos_env__$private
    type <- d$meta$incidence$type
    needs <- !(is.numeric(type) && type < 2L) && length(priv$rr) > 0L
    fn <- d$get_parf_filename()
    data.table(
      disease = d$name,
      needs_parf = needs,
      parf_table = basename(fn),
      exists = dir.exists(fn) || file.exists(fn),
      registry_valid = if (needs) isTRUE(priv$check_disease_registry_valid(TRUE))
                       else NA
    )
  }))
  fwrite(out, file.path(RES_DIR, "parf_preflight.csv"))
  bad <- out[needs_parf == TRUE & !(exists & registry_valid)]
  if (nrow(bad) > 0L) {
    print(bad)
    stop("PARF table(s) would be regenerated (hours). Aborting before any ",
         "simulation. See ", file.path(RES_DIR, "parf_preflight.csv"),
         call. = FALSE)
  }
  cat(sprintf("Pre-flight: %d cached PARF table(s) found and valid.\n",
              out[needs_parf == TRUE, .N]))
  out
}

md5 <- function(f) unname(tools::md5sum(f))

if (!COMPARE_ONLY) {
  # ---- design YAML: stock design + overrides, written under the results dir
  design <- yaml::read_yaml(file.path(REPO, "inputs", "sim_design.yaml"))
  design[names(DESIGN_OVERRIDES)] <- DESIGN_OVERRIDES
  design_path <- file.path(RES_DIR, "sim_design.yaml")
  yaml::write_yaml(design, design_path)
  # The written YAML must read back to exactly the list we meant to write
  # (guards against YAML 1.1 surprises such as an unquoted `n` key).
  if (!identical(yaml::read_yaml(design_path), design)) {
    print(all.equal(yaml::read_yaml(design_path), design))
    stop("The design YAML does not round-trip; refusing to run.", call. = FALSE)
  }

  # ---- fresh output directory (the synthpop directory is kept and shared)
  assert_inside_base(OUT_DIR)
  assert_inside_base(SP_DIR)
  if (dir.exists(OUT_DIR)) unlink(OUT_DIR, recursive = TRUE)
  dir.create(SP_DIR, recursive = TRUE, showWarnings = FALSE)

  IMPACTncd <- Simulation$new(design_path)
  stopifnot(identical(
    normalizePath(IMPACTncd$design$sim_prm$output_dir, mustWork = FALSE),
    normalizePath(OUT_DIR, mustWork = FALSE)
  ))
  parf_tbl <- preflight_parf(IMPACTncd)

  # Baseline first: every other scenario copies its population weights from
  # the sc0 lifecourse. sc0 keeps the constructor's default (no-op) scenario.
  timed("run sc0", IMPACTncd$run(MC, multicore = TRUE, "sc0"))

  # dimd_future and dimd_future_orig run the SAME scenario function; they
  # differ only in parf_deprivation (see PARF_DEPRIVATION and the header).
  scn_funs <- list(
    noop_future      = mk_noop_future(),
    xps_future       = mk_xps_future(SC_YEAR_SHORT, XPS_SHARE),
    dimd_future      = mk_dimd_future(SC_YEAR_SHORT),
    dimd_future_orig = mk_dimd_future(SC_YEAR_SHORT)
  )
  stopifnot(setequal(names(scn_funs), SCENARIOS),
            setequal(names(PARF_DEPRIVATION), SCENARIOS))
  for (scn in SCENARIOS) {
    IMPACTncd$update_primary_prevention_scn(
      scn_funs[[scn]], parf_deprivation = PARF_DEPRIVATION[[scn]]
    )
    stopifnot(identical(IMPACTncd$get_parf_deprivation(),
                        PARF_DEPRIVATION[[scn]]))
    timed(paste("run", scn), IMPACTncd$run(MC, multicore = TRUE, scn))
  }

  # ---- provenance of the simulation
  git <- function(...) {
    tryCatch(
      suppressWarnings(system2("git", c("-C", REPO, ...), stdout = TRUE,
                               stderr = TRUE)),
      error = function(e) paste("git failed:", conditionMessage(e))
    )
  }
  src_files <- sort(c(
    list.files(file.path(PKG_SRC, "R"), pattern = "\\.R$", full.names = TRUE),
    list.files(file.path(PKG_SRC, "src"), pattern = "\\.(cpp|h|hpp)$",
               full.names = TRUE)
  ))
  sp_files <- sort(list.files(SP_DIR, pattern = "^synthpop_.*\\.fst$",
                              full.names = TRUE))
  rdb <- system.file("R", "IMPACTncdEngland.rdb", package = "IMPACTncdEngland")
  # Record the exact source the run used, including any uncommitted fix.
  writeLines(
    git("diff", "--", "Rpackage/IMPACTncd_England_model_pkg/R",
        "Rpackage/IMPACTncd_England_model_pkg/src"),
    file.path(RES_DIR, "package_source_uncommitted.diff")
  )
  provenance <- list(
    label = LABEL,
    script = "auxil/verify_init_prvl_crn.R",
    started_at = format(T_START, "%Y-%m-%d %H:%M:%S %Z"),
    simulation_finished_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    host = unname(Sys.info()[["nodename"]]),
    r_version = R.version.string,
    thread_env = as.list(Sys.getenv(names(THREAD_VARS))),
    git_head = git("rev-parse", "HEAD"),
    git_branch = git("rev-parse", "--abbrev-ref", "HEAD"),
    git_status_rpackage = git("status", "--porcelain", "--",
                              "Rpackage/IMPACTncd_England_model_pkg"),
    md5_Disease_class_R = md5(file.path(PKG_SRC, "R", "Disease_class.R")),
    md5_Simulation_class_R = md5(file.path(PKG_SRC, "R", "Simulation_class.R")),
    md5_IMPACTncd_sim_cpp = md5(file.path(PKG_SRC, "src", "IMPACTncd_sim.cpp")),
    md5_all_package_source = digest::digest(
      paste(basename(src_files), md5(src_files), collapse = ";"),
      algo = "md5", serialize = FALSE
    ),
    package_versions = list(
      IMPACTncdEngland = as.character(packageVersion("IMPACTncdEngland")),
      CKutils = as.character(packageVersion("CKutils")),
      data.table = as.character(packageVersion("data.table")),
      dqrng = as.character(packageVersion("dqrng")),
      wrswoR = as.character(packageVersion("wrswoR")),
      arrow = as.character(packageVersion("arrow"))
    ),
    package_build = if (REBUILD) "rebuild" else "installed",
    package_build_detail = PKG_MODE,
    source_unchanged_since_last_global_R_install = SOURCE_MATCHES_INSTALL,
    installed_package_path = dirname(dirname(rdb)),
    md5_installed_rdb = if (nzchar(rdb)) md5(rdb) else NA_character_,
    base_dir = BASE_DIR,
    design_source = "inputs/sim_design.yaml",
    md5_design_source = md5(file.path(REPO, "inputs", "sim_design.yaml")),
    design_used = design_path,
    design_overrides = DESIGN_OVERRIDES,
    mc = MC,
    sc_year_short = SC_YEAR_SHORT,
    sc_year_source = SC_YEAR_SOURCE,
    scenarios = c("sc0", SCENARIOS),
    parf_deprivation = as.list(PARF_DEPRIVATION),
    xps_future_row_share = XPS_SHARE,
    parf_tables_used = parf_tbl[needs_parf == TRUE, parf_table],
    synthpop_files = setNames(as.list(md5(sp_files)), basename(sp_files))
  )
  yaml::write_yaml(provenance, file.path(RES_DIR, "provenance.yaml"))
  rm(IMPACTncd)
  invisible(gc())
}

# ===========================================================================
# 3. Comparison of each scenario with sc0, joined on (pid, year)
# ===========================================================================
cat("\n=== Comparing lifecourses with sc0 ===\n")

lc_path <- function(mc, scn) {
  file.path(OUT_DIR, "lifecourse", paste0("mc=", mc), paste0("scenario=", scn))
}

# Which arms the output tree holds. After a full run every arm must be there.
# With --compare-only an arm the tree lacks -- typically one added to this
# script after the tree was written, such as dimd_future_orig -- is skipped
# with a note rather than failed; the baseline itself is always required.
arm_missing_mc <- lapply(setNames(nm = c("sc0", SCENARIOS)), function(s) {
  MC[!dir.exists(lc_path(MC, s))]
})
if (length(arm_missing_mc[["sc0"]])) {
  stop("No sc0 lifecourse for mc ", paste(arm_missing_mc[["sc0"]],
       collapse = ","), " under ", file.path(OUT_DIR, "lifecourse"),
       call. = FALSE)
}
SKIPPED_SCN <- SCENARIOS[lengths(arm_missing_mc[SCENARIOS]) > 0L]
if (length(SKIPPED_SCN) && !COMPARE_ONLY) {
  stop("The run did not produce a lifecourse for: ",
       paste(SKIPPED_SCN, collapse = ", "), call. = FALSE)
}
CMP_SCN <- setdiff(SCENARIOS, SKIPPED_SCN)
if (!length(CMP_SCN)) {
  stop("No scenario lifecourse to compare with sc0 under ",
       file.path(OUT_DIR, "lifecourse"), call. = FALSE)
}
for (scn in SKIPPED_SCN) {
  cat(sprintf(
    "NOTE: %s has no lifecourse for mc %s in %s (a tree written before the arm existed?); its checks are skipped.\n",
    scn, paste(arm_missing_mc[[scn]], collapse = ","), OUT_DIR
  ))
}

read_lc <- function(mc, scn) {
  p <- lc_path(mc, scn)
  if (!dir.exists(p)) stop("Lifecourse not found: ", p, call. = FALSE)
  dt <- CKutils::read_parquet_dt(p)
  dt[, `:=`(pid = as.integer(pid), year = as.integer(year))]
  if (anyDuplicated(dt, by = c("pid", "year"))) {
    stop("Duplicate (pid, year) rows in ", p, call. = FALSE)
  }
  setkeyv(dt, c("pid", "year"))
  dt
}

# Element-wise "differs", treating NA == NA as equal and NA vs value as a
# difference. Factors are compared on their labels.
neq <- function(a, b) {
  if (is.factor(a)) a <- as.character(a)
  if (is.factor(b)) b <- as.character(b)
  na_a <- is.na(a)
  na_b <- is.na(b)
  out <- na_a != na_b
  both <- !na_a & !na_b
  out[both] <- a[both] != b[both]
  out
}

compare_pair <- function(scn, mc) {
  lc0 <- read_lc(mc, "sc0")
  lcs <- read_lc(mc, scn)
  key <- c("pid", "year")

  # Years: the lifecourse stores short years (13 = 2013); detect rather than
  # assume, and report long years throughout.
  yr_off <- if (max(lc0$year) < 1000L) 2000L else 0L
  sc_cut <- 2000L + SC_YEAR_SHORT - yr_off
  init_y <- min(lc0$year) # lifecourse is pruned to year >= init_year

  # EVERY column the two lifecourses share is compared (checks A and B), so
  # that a residual common-random-number channel cannot hide in a column
  # outside the outcome set (weights, strata, markers, ...). The outcome
  # subset is kept separately for the positive control (check C).
  cols0 <- setdiff(names(lc0), key)
  cols1 <- setdiff(names(lcs), key)
  cmp_cols <- intersect(cols0, cols1)
  out_cols <- grep(CMP_PATTERN, cmp_cols, value = TRUE)
  schema_diff <- c(
    if (length(setdiff(cols0, cols1))) paste0("only in sc0: ",
      paste(setdiff(cols0, cols1), collapse = ",")),
    if (length(setdiff(cols1, cols0))) paste0("only in ", scn, ": ",
      paste(setdiff(cols1, cols0), collapse = ","))
  )
  san_cols <- intersect(SANITY_COLS, cmp_cols)

  # Row-set mismatches (e.g. a simulant who dies a year earlier in one arm).
  only0 <- lc0[!lcs, on = key, key, with = FALSE]
  onlys <- lcs[!lc0, on = key, key, with = FALSE]

  m <- merge(
    lc0[, c(key, cmp_cols), with = FALSE],
    lcs[, c(key, cmp_cols), with = FALSE],
    by = key, suffixes = c(".sc0", ".scn")
  )
  d_cols <- lapply(cmp_cols, function(cc) {
    neq(m[[paste0(cc, ".sc0")]], m[[paste0(cc, ".scn")]])
  })
  names(d_cols) <- cmp_cols
  # TRUE for a row on which any of `cols` differs.
  any_of <- function(cols) {
    if (length(cols)) Reduce(`|`, d_cols[cols]) else logical(nrow(m))
  }
  any_diff <- any_of(cmp_cols)
  out_diff <- any_of(out_cols)
  san_diff <- any_of(san_cols)

  period_of <- function(y) fifelse(y < sc_cut, "pre", "post")
  d <- data.table(pid = m$pid, year = m$year, period = period_of(m$year),
                  dif = any_diff, out = out_diff, san = san_diff)

  # ---- by year
  by_year <- d[, .(rows_compared = .N, rows_any_diff = sum(dif),
                   share_any_diff = mean(dif)), keyby = year]
  mm_year <- rbind(
    only0[, .(rows_only_sc0 = .N, rows_only_scn = 0L), keyby = year],
    onlys[, .(rows_only_sc0 = 0L, rows_only_scn = .N), keyby = year]
  )[, .(rows_only_sc0 = sum(rows_only_sc0), rows_only_scn = sum(rows_only_scn)),
    keyby = year]
  by_year <- merge(by_year, mm_year, by = "year", all = TRUE)
  setnafill(by_year, fill = 0,
            cols = c("rows_compared", "rows_any_diff", "rows_only_sc0",
                     "rows_only_scn"))

  # ---- by year x column, then by period x column
  by_year_col <- rbindlist(lapply(cmp_cols, function(cc) {
    data.table(year = m$year, x = d_cols[[cc]])[
      , .(rows = .N, n_diff = sum(x)), keyby = year][, column := cc]
  }))
  by_col <- by_year_col[, .(
    rows = sum(rows), n_diff = sum(n_diff), share_diff = sum(n_diff) / sum(rows),
    earliest_diff_year = if (any(n_diff > 0L)) min(year[n_diff > 0L]) + yr_off
                         else NA_integer_
  ), keyby = .(period = period_of(year), column)]
  by_col[, outcome_column := column %chin% out_cols]

  # ---- by period
  summ <- d[, .(
    rows_compared = .N,
    rows_any_diff = sum(dif),
    share_any_diff = mean(dif),
    rows_outcome_diff = sum(out),
    share_outcome_diff = mean(out),
    pids_compared = uniqueN(pid),
    pids_any_diff = uniqueN(pid[dif]),
    earliest_diff_year = if (any(dif)) min(year[dif]) + yr_off else NA_integer_,
    rows_sanity_diff = sum(san)
  ), keyby = period]
  # rep(.., .N) rather than a bare scalar: on a zero-row table a scalar would
  # be recycled into one spurious all-NA row.
  mm_per <- rbind(
    only0[, .(side = rep("sc0", .N), period = period_of(year), year)],
    onlys[, .(side = rep("scn", .N), period = period_of(year), year)]
  )
  summ[, `:=`(
    rows_only_sc0 = vapply(period, function(p) mm_per[side == "sc0" & period == p, .N], integer(1)),
    rows_only_scn = vapply(period, function(p) mm_per[side == "scn" & period == p, .N], integer(1)),
    earliest_mismatch_year = vapply(period, function(p) {
      y <- mm_per[period == p, year]
      if (length(y)) min(y) + yr_off else NA_integer_
    }, integer(1)),
    top_columns = vapply(period, function(p) {
      tc <- by_col[period == p & n_diff > 0L][order(-n_diff)][seq_len(min(.N, 5L))]
      if (nrow(tc)) paste0(tc$column, "(", tc$n_diff, ")", collapse = " ") else ""
    }, character(1)),
    schema_diff = paste(schema_diff, collapse = "; ")
  )]

  # ---- init-year prevalence diagnostic: are the per-stratum COUNTS of
  # prevalent cases preserved while the identity of the cases changes?
  dis_cols <- grep("^cms", grep("_prvl$", cmp_cols, value = TRUE),
                   value = TRUE, invert = TRUE)
  strata <- intersect(PRVL_STRATA, intersect(names(lc0), names(lcs)))
  iy0 <- lc0[year == init_y]
  iys <- lcs[year == init_y]
  mi <- m[year == init_y]
  init_tab <- rbindlist(lapply(dis_cols, function(cc) {
    a <- mi[[paste0(cc, ".sc0")]]
    b <- mi[[paste0(cc, ".scn")]]
    c0 <- iy0[iy0[[cc]] > 0L, .N, keyby = strata]
    c1 <- iys[iys[[cc]] > 0L, .N, keyby = strata]
    cc01 <- merge(c0, c1, by = strata, all = TRUE, suffixes = c("_sc0", "_scn"))
    setnafill(cc01, fill = 0L, cols = c("N_sc0", "N_scn"))
    data.table(
      column = cc,
      prevalent_sc0 = sum(iy0[[cc]] > 0L),
      prevalent_scn = sum(iys[[cc]] > 0L),
      matched_rows = nrow(mi),
      prevalence_status_differs = sum((a > 0L) != (b > 0L)),
      value_differs = sum(neq(a, b)),
      strata_with_cases = nrow(cc01),
      strata_count_equal = sum(cc01$N_sc0 == cc01$N_scn)
    )
  }))

  tag <- function(x) {
    x[, `:=`(scenario = scn, mc = mc)]
    setcolorder(x, c("scenario", "mc"))
    x
  }
  by_year[, year := year + yr_off]
  by_year_col[, year := year + yr_off]
  init_tab[, init_year := init_y + yr_off]
  list(
    summary = tag(summ), by_year = tag(by_year), by_col = tag(by_col),
    by_year_col = tag(by_year_col), init = tag(init_tab)
  )
}

res <- list()
timed("comparison", {
  for (scn in CMP_SCN) for (mc in MC) {
    res[[paste(scn, mc)]] <- compare_pair(scn, mc)
  }
})
pick <- function(what) rbindlist(lapply(res, `[[`, what), use.names = TRUE)
summary_dt  <- pick("summary")
by_year     <- pick("by_year")
by_col      <- pick("by_col")
by_year_col <- pick("by_year_col")
init_dt     <- pick("init")
summary_dt[, period := factor(period, levels = c("pre", "post"))]
setorder(summary_dt, scenario, mc, period)

fwrite(summary_dt,  file.path(RES_DIR, "summary.csv"))
fwrite(by_year,     file.path(RES_DIR, "by_year.csv"))
fwrite(by_col,      file.path(RES_DIR, "by_column.csv"))
fwrite(by_year_col, file.path(RES_DIR, "by_year_column.csv"))
fwrite(init_dt,     file.path(RES_DIR, "init_year_prevalence.csv"))

# ---- concise console summary
cat(sprintf("\nPeriods: pre = years < %d, post = years >= %d\n\n",
            2000L + SC_YEAR_SHORT, 2000L + SC_YEAR_SHORT))
options(width = 150L)
cat("share_diff: any column; share_out: outcome columns only (check C)\n")
print(summary_dt[, .(
  scenario, mc, period, rows = rows_compared,
  share_diff = sprintf("%.4f", share_any_diff),
  share_out = sprintf("%.4f", share_outcome_diff), pids_diff = pids_any_diff,
  first_diff = earliest_diff_year, only_sc0 = rows_only_sc0,
  only_scn = rows_only_scn
)], row.names = FALSE, class = FALSE)
cat("\nMost-differing columns (n rows):\n")
invisible(summary_dt[nzchar(top_columns),
  cat(sprintf("  %-16s mc %d %-4s %s\n", scenario, mc, period, top_columns),
      sep = "")])

cat("\nInit-year prevalence, all diseases pooled (counts per stratum vs ",
    "identity of the cases):\n", sep = "")
print(init_dt[, .(
  prevalent_sc0 = sum(prevalent_sc0), prevalent_scn = sum(prevalent_scn),
  status_differs = sum(prevalence_status_differs),
  strata_count_equal = sprintf("%d/%d", sum(strata_count_equal),
                               sum(strata_with_cases))
), keyby = .(scenario, mc)], row.names = FALSE, class = FALSE)

# ===========================================================================
# 4. Checks
# ===========================================================================
cat("\n=== Checks ===\n")
# TRUE when every period in `s` is identical to sc0, row set included. An
# empty selection is never "clean", so a missing period cannot pass silently.
clean <- function(s) {
  nrow(s) > 0L &&
    s[, all(rows_any_diff == 0L & rows_only_sc0 == 0L & rows_only_scn == 0L)]
}
schema_ok <- summary_dt[, all(!nzchar(schema_diff))]
say(schema_ok, sprintf(
  "every scenario carries the same columns as sc0 (%d compared, pid and year excluded)",
  uniqueN(by_col$column)
))
for (scn in SKIPPED_SCN) {
  skip(sprintf("%s: no lifecourse for mc %s in this output tree; checks A-D skipped for it",
               scn, paste(arm_missing_mc[[scn]], collapse = ",")))
}
SC_LONG <- 2000L + SC_YEAR_SHORT

# A. negative control
if ("noop_future" %in% CMP_SCN) for (mc_ in MC) {
  s <- summary_dt[scenario == "noop_future" & mc == mc_]
  say(clean(s), sprintf(
    "A. noop_future == sc0 in every row, column and year (mc %d)%s", mc_,
    if (clean(s)) "" else if (nrow(s) == 0L) " -- no rows compared" else
      paste0(" -- differing: ", paste(unique(by_col[
        scenario == "noop_future" & mc == mc_ & n_diff > 0L, column]),
        collapse = ", "))
  ))
}

# B. the issue-#8 property. Every message part is guarded, so that a period
# with no rows (a lifecourse that starts at or after sc_year) gives a clear
# FAIL line rather than an error from a zero-length condition.
for (scn in intersect(ACTIVE_SCN, CMP_SCN)) for (mc_ in MC) {
  s <- summary_dt[scenario == scn & mc == mc_ & period == "pre"]
  if (nrow(s) != 1L) {
    say(FALSE, sprintf(
      "B. %-16s == sc0 before %d (mc %d): NO simulated rows before %d to compare",
      scn, SC_LONG, mc_, SC_LONG
    ))
    next
  }
  say(clean(s), sprintf(
    "B. %-16s == sc0 before %d (mc %d): %d of %d rows differ (%.1f%%), %d+%d unmatched rows%s",
    scn, SC_LONG, mc_, s$rows_any_diff, s$rows_compared,
    100 * s$share_any_diff, s$rows_only_sc0, s$rows_only_scn,
    if (isTRUE(!is.na(s$earliest_diff_year)))
      sprintf(", first in %d", s$earliest_diff_year) else ""
  ))
}

# C. positive control, on the OUTCOME columns and the row set only (see the
# header for why not every column).
for (scn in intersect(ACTIVE_SCN, CMP_SCN)) for (mc_ in MC) {
  s <- summary_dt[scenario == scn & mc == mc_ & period == "post"]
  if (nrow(s) != 1L) {
    say(FALSE, sprintf(
      "C. %-16s differs from sc0 from %d (mc %d): NO simulated rows from %d on to compare",
      scn, SC_LONG, mc_, SC_LONG
    ))
    next
  }
  acts <- isTRUE(s$rows_outcome_diff > 0L || s$rows_only_sc0 > 0L ||
                   s$rows_only_scn > 0L)
  say(acts, sprintf(
    "C. %-16s differs from sc0 from %d (mc %d): outcomes differ in %.1f%% of matched rows, %d+%d unmatched rows",
    scn, SC_LONG, mc_, 100 * s$share_outcome_diff, s$rows_only_sc0,
    s$rows_only_scn
  ))
}

# D. strata agree on matched rows
say(summary_dt[, all(rows_sanity_diff == 0L)],
    "D. age, sex and dimd agree with sc0 on every matched row")

# E. no PARF regeneration
log_files <- list.files(file.path(OUT_DIR, "logs"), full.names = TRUE)
if (length(log_files) > 0L) {
  regen <- vapply(log_files, function(f) {
    any(grepl("Generating new parf file|Saving parf file to|Saving parf cache",
              readLines(f, warn = FALSE)))
  }, logical(1))
  say(!any(regen), sprintf(
    "E. no PARF table regenerated (%d log file(s) scanned)%s", length(log_files),
    if (any(regen)) paste0(" -- see ", paste(basename(log_files[regen]),
                                             collapse = ", ")) else ""
  ))
} else {
  skip(paste("E. no log files found under", file.path(OUT_DIR, "logs")))
}

timings_dt <- rbindlist(timings)
fwrite(timings_dt, file.path(RES_DIR,
       if (COMPARE_ONLY) "timings_compare_only.csv" else "timings.csv"))
writeLines(c(
  check_log,
  if (ok) "ALL CHECKS PASSED" else "*** SOME CHECKS FAILED ***"
), file.path(RES_DIR, "checks.txt"))

cat(sprintf("\nTotal %.1f min. Results in %s\n",
            as.numeric(difftime(Sys.time(), T_START, units = "mins")), RES_DIR))
cat(sprintf("\n%s\n", if (ok) "ALL CHECKS PASSED" else "*** SOME CHECKS FAILED ***"))
quit(status = if (ok) 0L else 1L)
