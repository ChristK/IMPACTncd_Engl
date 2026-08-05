# =============================================================================
# Prove the comparator-map change leaves the SCALAR path byte-identical.
#
# `comparator_scenario = "sc0"` (an unnamed length-1 character) must produce
# exactly the tables it produced before the change: same file list, same bytes,
# no `comparator` column, no `comparators.csv`. This is the claim that matters,
# because these tables are already published -- the tinytest suite can only
# check a fixture, this checks 300+ real tables built from 100 MC iterations.
#
# Compares:
#   tables_pre_cmpmap  <- built by the equity-fix build (commit 7922334)
#   tables             <- rebuilt by the comparator-map build, scalar path
#
# Run from the project root, AFTER re-running genvasc_las_reexport_tables.R:
#   Rscript scenarios/polygenic/genvasc_las_scalar_identity_check.R
# =============================================================================

suppressMessages(library(data.table))

ROOT <- "/mnt/storage_fast4/genvasc_las/output"
OLD  <- file.path(ROOT, "tables_pre_cmpmap")
NEW  <- file.path(ROOT, "tables")
stopifnot(dir.exists(OLD), dir.exists(NEW))

# Other scripts write into the same directory (genvasc_cea.R, the deck builder).
# They are not export_tables() outputs, so a re-export neither produces nor
# should be judged on them.
NOT_EXPORT_TABLES <- c("genvasc_cea_iterations.csv", "genvasc_cea_summary.csv",
                       "genvasc_deck.html", "genvasc_presentation_data.json")

old_f <- setdiff(sort(list.files(OLD)), NOT_EXPORT_TABLES)
new_f <- setdiff(sort(list.files(NEW)), NOT_EXPORT_TABLES)

cat("files before:", length(old_f), " after:", length(new_f), "\n")
extra <- setdiff(new_f, old_f)
gone  <- setdiff(old_f, new_f)
cat("only after :", if (length(extra)) paste(extra, collapse = ", ") else "(none)", "\n")
cat("only before:", if (length(gone))  paste(gone,  collapse = ", ") else "(none)", "\n")

# `comparators.csv` is written in MAP MODE ONLY. Its presence here would mean
# the syntactic gate leaked.
cat("\ncomparators.csv present (must be FALSE in scalar mode):",
    file.exists(file.path(NEW, "comparators.csv")), "\n")

# --- byte-for-byte -----------------------------------------------------------
common <- intersect(old_f, new_f)
md5_old <- tools::md5sum(file.path(OLD, common))
md5_new <- tools::md5sum(file.path(NEW, common))
differing <- common[unname(md5_old) != unname(md5_new)]

cat("\n=== byte-for-byte over", length(common), "files ===\n")
cat("identical:", length(common) - length(differing),
    "  differing:", length(differing),
    if (length(differing) == 0L) " (as expected)" else " <- UNEXPECTED", "\n")

if (length(differing)) {
  cat("\nfiles that differ:\n")
  for (f in head(differing, 20)) cat("  ", f, "\n")
  # Where they differ, say HOW -- a new column is a different failure from a
  # changed number, and only one of them is a schema leak.
  cat("\nfirst differing file, column diff:\n")
  o <- fread(file.path(OLD, differing[1L]))
  n <- fread(file.path(NEW, differing[1L]))
  cat("  cols only in new:", paste(setdiff(names(n), names(o)), collapse = ", "), "\n")
  cat("  cols only in old:", paste(setdiff(names(o), names(n)), collapse = ", "), "\n")
  cat("  nrow old/new    :", nrow(o), "/", nrow(n), "\n")
  if (identical(names(o), names(n)) && nrow(o) == nrow(n)) {
    num <- names(o)[vapply(o, is.numeric, logical(1))]
    for (cl in num) {
      d <- max(abs(o[[cl]] - n[[cl]]), na.rm = TRUE)
      if (isTRUE(d > 0)) cat("    ", cl, " max abs diff ", d, "\n", sep = "")
    }
  }
}

cat("\n", if (length(differing) == 0L && length(extra) == 0L && length(gone) == 0L)
  "SCALAR PATH IS BYTE-IDENTICAL" else "SCALAR PATH CHANGED -- INVESTIGATE", "\n", sep = "")
