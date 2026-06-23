#!/usr/bin/env Rscript
# =============================================================================
# test_new_user_pipeline.R
#
# Simulates the full pipeline a new user would follow:
#   1. Clone the repo from GitHub into a temp directory
#   2. Install the R package (via global.R)
#   3. Create a Simulation object
#   4. Download data from Zenodo
#   5. Run simulate_testing (baseline + scenario, export)
#
# Requirements:
#   - git must be on the PATH
#   - A Zenodo token must be available (env var or pass)
#   - Internet access (GitHub clone + Zenodo download)
#
# Usage:
#   Rscript testing/test_new_user_pipeline.R
#
# Environment variables:
#   ZENODO_TOKEN              — Zenodo API token (production)
#   ZENODO_SANDBOX_TOKEN      — Zenodo sandbox token (default)
#   IMPACTNCD_CONCEPT_DOI     — override concept DOI (optional)
#   IMPACTNCD_SANDBOX         — "TRUE" or "FALSE" (default TRUE)
#   IMPACTNCD_BRANCH          — branch to clone (default: current branch)
#   IMPACTNCD_KEEP_TMPDIR     — "TRUE" to keep temp dir after test (default FALSE)
# =============================================================================

cat("
============================================================
  IMPACTncd England — New User Pipeline Test
============================================================
\n")

# --- Configuration ----------------------------------------------------------

# Determine branch: env var > current git branch > fallback "main"
branch <- Sys.getenv("IMPACTNCD_BRANCH", unset = "")
if (!nzchar(branch)) {
  branch <- tryCatch(
    trimws(system("git rev-parse --abbrev-ref HEAD", intern = TRUE)),
    error = function(e) "main"
  )
}

# Default to PRODUCTION: this exercises the real new-user path (the published
# public record). Set IMPACTNCD_SANDBOX=TRUE to test against the sandbox.
use_sandbox <- as.logical(
  Sys.getenv("IMPACTNCD_SANDBOX", unset = "FALSE")
)

concept_doi <- Sys.getenv("IMPACTNCD_CONCEPT_DOI", unset = "")
if (!nzchar(concept_doi)) {
  concept_doi <- if (use_sandbox) "10.5072/zenodo.442996" else "10.5281/zenodo.20812409"
}

# A token is OPTIONAL: downloading published public data works anonymously.
# A token is only required to upload/publish.
token <- Sys.getenv(
  if (use_sandbox) "ZENODO_SANDBOX_TOKEN" else "ZENODO_TOKEN",
  unset = ""
)
if (!nzchar(token)) {
  token <- Sys.getenv("ZENODO_TOKEN", unset = "")
}
if (!nzchar(token)) {
  message(
    "No Zenodo token found; connecting anonymously (fine for downloading ",
    "published data)."
  )
  token <- NULL
}

keep_tmpdir <- as.logical(
  Sys.getenv("IMPACTNCD_KEEP_TMPDIR", unset = "FALSE")
)

repo_url <- "https://github.com/ChristK/IMPACTncd_Engl.git"

cat("  Repo:        ", repo_url, "\n")
cat("  Branch:      ", branch, "\n")
cat("  Sandbox:     ", use_sandbox, "\n")
cat("  Concept DOI: ", concept_doi, "\n")
cat("  Keep tmpdir: ", keep_tmpdir, "\n\n")

# =============================================================================
# STEP 1: Clone the repository
# =============================================================================

cat("=== Step 1: Cloning repository ===\n")

tmpdir <- file.path(tempdir(), paste0("IMPACTncd_test_", format(Sys.time(), "%Y%m%d_%H%M%S")))
cat("  Temp directory: ", tmpdir, "\n")

clone_cmd <- paste0(
  "git clone --branch ", shQuote(branch),
  " --single-branch --depth 1 ",
  shQuote(repo_url), " ", shQuote(tmpdir)
)

clone_status <- system(clone_cmd)
if (clone_status != 0L) {
  stop("git clone failed with exit code ", clone_status, call. = FALSE)
}

cat("  Clone complete.\n\n")

# Clean up on exit unless user wants to keep it
if (!keep_tmpdir) {
  on.exit({
    cat("\nCleaning up temp directory: ", tmpdir, "\n")
    unlink(tmpdir, recursive = TRUE)
  }, add = TRUE)
}

# Switch working directory to the cloned repo (as a new user would)
original_wd <- getwd()
on.exit(setwd(original_wd), add = TRUE)
setwd(tmpdir)
cat("  Working directory: ", getwd(), "\n\n")

# =============================================================================
# STEP 2: Install packages via global.R
# =============================================================================

cat("=== Step 2: Installing packages (global.R) ===\n")

source("./global.R")

cat("  Package installation complete.\n\n")

# =============================================================================
# STEP 3: Create Simulation object
# =============================================================================

cat("=== Step 3: Creating Simulation object ===\n")

library(IMPACTncdEngland)
IMPACTncd <- Simulation$new("testing/sim_design_testing.yaml")

cat("  Simulation object created (skeleton mode — no data yet).\n\n")

# =============================================================================
# STEP 4: Download data from Zenodo
# =============================================================================

cat("=== Step 4: Connecting to Zenodo and downloading data ===\n")

# Install Zenodo dependencies if not already present
for (pkg in c("zen4R", "httr2")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Installing ", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

IMPACTncd$zenodo_connect(
  token       = token,
  concept_doi = concept_doi,
  sandbox     = use_sandbox
)

IMPACTncd$zenodo_download_all()

cat("  Download and initialization complete.\n\n")

# =============================================================================
# STEP 5: Run simulate_testing pipeline
# =============================================================================

cat("=== Step 5: Running simulation (baseline sc0) ===\n")

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:2, multicore = TRUE, "sc0")

cat("  Baseline scenario complete.\n\n")

# --- Intervention scenario ---

cat("=== Step 5b: Running intervention scenario (sc1) ===\n")

IMPACTncd$update_primary_prevention_scn(
  function(synthpop) {
    synthpop$pop[year >= 25L, bmi_curr_xps := bmi_curr_xps * 0.8]
  }
)

IMPACTncd$update_secondary_prevention_scn(
  function(synthpop) {
    sc_year <- 23L
    change <- 0.2
    synthpop$pop[
      year >= sc_year,
      `:=`(
        prb_af_mrtl1 = prb_af_mrtl1 * (1 - change),
        prb_af_mrtl2 = prb_af_mrtl2 * (1 - change)
      )
    ]
  }
)

IMPACTncd$run(1:2, multicore = TRUE, "sc1")

cat("  Intervention scenario complete.\n\n")

# --- Export results ---

cat("=== Step 5c: Exporting summaries and tables ===\n")

IMPACTncd$export_summaries(multicore = TRUE)
IMPACTncd$export_tables()

cat("  Export complete.\n\n")

# =============================================================================
# Summary
# =============================================================================

cat("
============================================================
  New User Pipeline Test — PASSED
============================================================
  Branch:    ", branch, "
  Temp dir:  ", tmpdir, "
  Cleanup:   ", if (keep_tmpdir) "SKIPPED (IMPACTNCD_KEEP_TMPDIR=TRUE)" else "will remove on exit", "
============================================================
\n")
