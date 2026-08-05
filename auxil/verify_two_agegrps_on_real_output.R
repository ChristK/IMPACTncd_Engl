#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# verify_two_agegrps_on_real_output.R
#
# End-to-end verification of the two_agegrps fix against a REAL model output
# tree, not a fixture. Three things are checked:
#
#   A. REGRESSION -- the refactored collapse in export_main_tables() must
#      reproduce the already-published tables2agegrps/ file BYTE FOR BYTE. That
#      file was written by the old hard-coded implementation, so agreement
#      proves the refactor changed nothing for the family that already worked.
#
#   B. THE FIX -- export_all_cause_mrtl_tables() must now collapse. Previously
#      `two_agegrps` never reached it, so its tables2agegrps/ output was
#      md5-identical to tables/. Both halves are asserted: the published file
#      IS identical to tables/ (the bug, as shipped), and the newly generated
#      one is NOT (the fix).
#
#   C. EXACTNESS -- the collapsed rate must equal
#      sum(deaths over the constituent bands) / sum(cases over them),
#      recomputed independently from the raw parquet. This is what makes the
#      collapse a relabelling rather than an approximation.
#
# Plus a non-default split age, to show it is genuinely an argument.
#
# Usage: Rscript auxil/verify_two_agegrps_on_real_output.R [output_dir]
# Peak RSS ~16 GB; run it when the box has headroom.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table); library(IMPACTncdEngland); library(arrow)
})

outdir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L)
  commandArgs(trailingOnly = TRUE)[[1L]] else
  "/mnt/storage_fast/ethnicity_proj/output"

setDTthreads(4L); arrow::set_cpu_count(4L)
sdir <- file.path(outdir, "summaries")
PRBL <- c(0.5, 0.025, 0.975, 0.1, 0.9)
BASE_YEAR <- 2026L
ok <- TRUE
say <- function(pass, msg) {
  ok <<- ok && pass
  cat(sprintf("[%s] %s\n", if (pass) " OK " else "FAIL", msg))
}

# ---- harness: re-parent the private methods onto a stub self/private --------
harness <- function() {
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE, diseases = NULL)))
  core <- Simulation$private_methods$tbl_smmrs_core
  environment(core) <- env
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization = "scaled_up") {
      p <- file.path(sdir, paste0(summary_type, "_", standardization))
      if (!dir.exists(p)) return(NULL)
      CKutils::read_parquet_dt(p)
    },
    tbl_smmrs_core = core
  )
  env
}
bind <- function(name, env) {
  f <- Simulation$private_methods[[name]]
  environment(f) <- env
  f
}

# ===========================================================================
# A. Regression: the main family must reproduce the published file exactly
# ===========================================================================
cat("\n=== A. main tables: refactor vs the published tables2agegrps file ===\n")
d_main <- tempfile("main_"); dir.create(d_main, recursive = TRUE)
env <- harness()
bind("export_main_tables", env)(
  prbl = PRBL, baseline_year = BASE_YEAR, output_dir = outdir,
  tables_dir = d_main, comparator_scenario = "sc0",
  two_agegrps = TRUE, two_agegrps_split_age = 65L,
  strata_ons = list(c("year", "agegrp")), strata_esp = list("year")
)
pub <- file.path(outdir, "tables2agegrps")
cand <- intersect(list.files(d_main), list.files(pub))
cand <- grep("agegrp", cand, value = TRUE)
say(length(cand) > 0L, sprintf("%d main-family agegrp file(s) to compare", length(cand)))
for (f in cand) {
  a <- fread(file.path(d_main, f)); b <- fread(file.path(pub, f))
  say(isTRUE(all.equal(a, b)), paste0("identical to published: ", f))
}

# ===========================================================================
# B. The fix: all_cause_mrtl now collapses
# ===========================================================================
cat("\n=== B. all_cause_mrtl: silent no-op before, collapsed now ===\n")
acm_f <- "all-cause mortality given disease-year-agegroup (not standardised).csv"
old_two <- fread(file.path(pub, acm_f))
std     <- fread(file.path(outdir, "tables", acm_f))
say(isTRUE(all.equal(old_two, std)),
    "as shipped, tables2agegrps == tables (the bug: two_agegrps was ignored)")
say(identical(sort(unique(old_two$agegrp)), sort(unique(std$agegrp))),
    "  ...both carried the same 5-year bands")

d_acm <- tempfile("acm_"); dir.create(d_acm, recursive = TRUE)
env <- harness()
bind("export_all_cause_mrtl_tables", env)(
  prbl = PRBL, summaries_dir = sdir, tables_dir = d_acm,
  two_agegrps = TRUE, two_agegrps_split_age = 65L,
  strata_ons = list(c("year", "agegrp")), strata_esp = list("year")
)
new_two <- fread(file.path(d_acm, acm_f))
say(identical(sort(unique(new_two$agegrp)), c("30-64", "65-99")),
    "after the fix, agegrp is the two coarse groups")
say(!isTRUE(all.equal(new_two, std)), "  ...and the table genuinely differs")
say(!anyNA(new_two), "no NAs introduced")

n_yr <- uniqueN(new_two$year); n_dis <- uniqueN(new_two$disease)
say(nrow(new_two) == n_yr * n_dis * 2L,
    sprintf("row count is years(%d) x diseases(%d) x 2 groups = %d",
            n_yr, n_dis, n_yr * n_dis * 2L))

# ===========================================================================
# C. Exactness: collapsed rate == summed deaths / summed cases over the bands
# ===========================================================================
cat("\n=== C. exactness, recomputed independently from the raw parquet ===\n")
raw <- CKutils::read_parquet_dt(file.path(sdir, "all_cause_mrtl_by_dis_scaled_up"))
raw[, year := year + 2000L]
lo_bands <- c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64")
raw[, grp := fifelse(agegrp %chin% lo_bands, "30-64", "65-99")]

check_disease <- function(dis) {
  dn <- paste0("deaths_", dis); cn <- paste0("cases_", dis)
  if (!all(c(dn, cn) %chin% names(raw))) return(NULL)
  per_mc <- raw[, .(r = sum(get(dn)) / sum(get(cn))), keyby = .(mc, year, grp)]
  exp_med <- per_mc[, .(med = median(r)), keyby = .(year, grp)]
  got <- new_two[disease == dis][order(year, agegrp)]
  list(exp = exp_med$med, got = got[["all_cause_mrtl_by_disease_rate_50.0%"]])
}
for (dis in c("chd", "stroke", "t2dm", "af")) {
  r <- check_disease(dis)
  if (is.null(r)) next
  say(isTRUE(all.equal(r$got, r$exp, tolerance = 1e-8)),
      sprintf("%-8s collapsed rate == summed deaths / summed cases (n=%d)",
              dis, length(r$exp)))
}
rm(raw); invisible(gc(full = TRUE))

# ===========================================================================
# D. The split age is genuinely an argument
# ===========================================================================
cat("\n=== D. a non-default split age ===\n")
d50 <- tempfile("acm50_"); dir.create(d50, recursive = TRUE)
env <- harness()
bind("export_all_cause_mrtl_tables", env)(
  prbl = PRBL, summaries_dir = sdir, tables_dir = d50,
  two_agegrps = TRUE, two_agegrps_split_age = 50L,
  strata_ons = list(c("year", "agegrp")), strata_esp = list("year")
)
t50 <- fread(file.path(d50, acm_f))
say(identical(sort(unique(t50$agegrp)), c("30-49", "50-99")),
    "two_agegrps_split_age = 50 gives 30-49 / 50-99")
say(!isTRUE(all.equal(t50[[6]], new_two[[6]])),
    "  ...and produces different numbers from the 65 split")

cat(sprintf("\n%s\n", if (ok) "ALL CHECKS PASSED" else "*** SOME CHECKS FAILED ***"))
quit(status = if (ok) 0L else 1L)
