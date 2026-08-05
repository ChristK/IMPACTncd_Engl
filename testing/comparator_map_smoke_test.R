# =============================================================================
# End-to-end smoke test of MAP MODE through the real export_tables() path.
#
# The tinytest suite mocks the R6 object, so it never reaches
# private$check_comparator_map() or private$probe_summary_scenarios(). This runs
# the genuine article on the small testing output (sc0 + sc1), where the map
# c(sc1 = "sc0") is EQUIVALENT to the scalar "sc0" -- so the numbers must match
# exactly and the only difference must be the extra `comparator` column plus the
# `comparators.csv` manifest.
#
# It also checks the validators fire on real data (unknown scenario -> error).
#
# Run from the project root:  Rscript testing/comparator_map_smoke_test.R
# =============================================================================

source("./global.R")
suppressMessages(library(data.table))

ok <- function(label, cond) {
  cat(sprintf("%-64s %s\n", label, if (isTRUE(cond)) "PASS" else "**FAIL**"))
  if (!isTRUE(cond)) assign(".failed", TRUE, envir = globalenv())
}

YAML <- "testing/sim_design_testing.yaml"
IMPACTncd <- Simulation$new(YAML)
out <- IMPACTncd$design$sim_prm$output_dir
tdir <- file.path(out, "tables")
keep <- file.path(out, "tables_scalar_reference")

# ---- 1. scalar reference ----------------------------------------------------
cat("\n=== scalar export ===\n")
IMPACTncd$export_tables(multicore = FALSE)
unlink(keep, recursive = TRUE)
file.rename(tdir, keep)
ok("scalar run wrote no comparators.csv",
   !file.exists(file.path(keep, "comparators.csv")))

# ---- 2. equivalent map ------------------------------------------------------
cat("\n=== map export: c(sc1 = 'sc0') ===\n")
IMPACTncd$export_tables(multicore = FALSE,
                        comparator_scenario = c(sc1 = "sc0"))

ok("map run wrote comparators.csv", file.exists(file.path(tdir, "comparators.csv")))
man <- fread(file.path(tdir, "comparators.csv"))
ok("manifest is the resolved map",
   identical(man$scenario, "sc1") && identical(man$comparator, "sc0"))

sf <- setdiff(sort(list.files(keep)), "comparators.csv")
mf <- setdiff(sort(list.files(tdir)), "comparators.csv")
ok("same file list either way", identical(sf, mf))

# Contrast tables must gain `comparator`; everything else must not.
contrast_pat <- "^(cases prevented|case-years prevented|deaths prevented|net QALYs|net costs|cost-effectiveness|equity )"
n_gain <- 0L
n_same <- 0L
bad <- character(0)
for (f in intersect(sf, mf)) {
  o <- fread(file.path(keep, f))
  n <- fread(file.path(tdir, f))
  has <- "comparator" %in% names(n)
  is_contrast <- grepl(contrast_pat, f)
  if (has != is_contrast) bad <- c(bad, f)
  if (has) {
    n_gain <- n_gain + 1L
    if (!all(n$comparator == "sc0")) bad <- c(bad, paste(f, "(bad value)"))
    n2 <- copy(n)[, comparator := NULL]
    setcolorder(n2, names(o))
    if (!isTRUE(all.equal(o, n2, check.attributes = FALSE)))
      bad <- c(bad, paste(f, "(numbers moved)"))
  } else {
    n_same <- n_same + 1L
    if (!isTRUE(all.equal(o, n, check.attributes = FALSE)))
      bad <- c(bad, paste(f, "(non-contrast table changed)"))
  }
}
cat(sprintf("  contrast tables gaining `comparator`: %d\n", n_gain))
cat(sprintf("  other tables unchanged              : %d\n", n_same))
ok("comparator column exactly on the contrast tables, numbers unmoved",
   length(bad) == 0L)
if (length(bad)) cat("   offenders:\n     ", paste(head(bad, 10), collapse = "\n     "), "\n")

# ---- 3. validators fire on real data ----------------------------------------
cat("\n=== validation ===\n")
ok("unknown scenario name is a hard error",
   inherits(try(IMPACTncd$export_tables(
     multicore = FALSE, comparator_scenario = c(ghost = "sc0")),
     silent = TRUE), "try-error"))
ok("self-comparison is a hard error",
   inherits(try(IMPACTncd$export_tables(
     multicore = FALSE, comparator_scenario = c(sc1 = "sc1")),
     silent = TRUE), "try-error"))
ok("two unnamed elements is a hard error",
   inherits(try(IMPACTncd$export_tables(
     multicore = FALSE, comparator_scenario = c("sc0", "sc1")),
     silent = TRUE), "try-error"))

# ---- restore the scalar tables ----------------------------------------------
unlink(tdir, recursive = TRUE)
file.rename(keep, tdir)
cat("\nrestored the scalar tables to ", tdir, "\n", sep = "")

cat("\n", if (isTRUE(get0(".failed", ifnotfound = FALSE)))
  "SMOKE TEST FAILED\n" else "smoke test passed\n", sep = "")
