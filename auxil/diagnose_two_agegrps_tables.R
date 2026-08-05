#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# diagnose_two_agegrps_tables.R
#
# Why does export_tables(two_agegrps = TRUE) leave fewer files in
# <output_dir>/tables2agegrps than two_agegrps = FALSE leaves in
# <output_dir>/tables?
#
# Two independent things are going on, and only the first is about file COUNT:
#
#   (1) `strata`. export_tables() merges the caller's `strata` with defaults
#       PER FAMILY (build_strata_config): a family the caller names is REPLACED
#       wholesale, a family they omit takes the default. So the file count is
#       driven by the strata list passed to each call, not by two_agegrps.
#       two_agegrps only changes which DEFAULTS apply -- and only for families
#       the caller did not name.
#
#   (2) `two_agegrps` reaches only SOME families. export_tables_hlpr() passes it
#       to the "main" and "equity" tasks only; export_all_cause_mrtl_tables(),
#       export_disease_characteristics_tables() and export_xps_tables() do not
#       take the argument at all. So tables2agegrps/ is a MIXTURE: some families
#       carry the collapsed 30-64 / 65-99 bands, others keep the 5-year bands.
#
# This script reports both: the file-set difference attributed to strata, and
# the actual agegrp levels found in each family in each directory.
#
# Usage: Rscript auxil/diagnose_two_agegrps_tables.R [output_dir]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))

outdir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L)
  commandArgs(trailingOnly = TRUE)[[1L]] else
  "/mnt/storage_fast/ethnicity_proj/output"

d_std <- file.path(outdir, "tables")
d_two <- file.path(outdir, "tables2agegrps")
stopifnot(dir.exists(d_std), dir.exists(d_two))

f_std <- list.files(d_std, pattern = "\\.csv$")
f_two <- list.files(d_two, pattern = "\\.csv$")

cat(sprintf("tables/          %d files\ntables2agegrps/  %d files\n\n",
            length(f_std), length(f_two)))

# ---- (1) which strata combinations are missing -----------------------------
# Filenames are "<description>-<strata joined by ->  (<standardisation>).csv"
# (the description itself contains " by ", and agegrp is rewritten to agegroup
# in the all_cause_mrtl family). Pull out the strata part.
strata_of <- function(f) {
  x <- sub("\\s*\\((not standardised|[^)]*standardised)\\)\\.csv$", "", f)
  x <- sub("^.*? by ", "", x)                 # main / xps / disease_char
  x <- sub("^all-cause mortality given disease-", "", x)  # all_cause_mrtl
  sub(" popdenom$", "", x)
}

missing <- setdiff(f_std, f_two)
extra   <- setdiff(f_two, f_std)

cat("=== file-set difference ===\n")
cat(sprintf("only in tables/         : %d\n", length(missing)))
cat(sprintf("only in tables2agegrps/ : %d\n\n", length(extra)))

cat("Missing files, grouped by the strata combination they came from:\n")
m <- data.table(file = missing, strata = strata_of(missing))
print(m[, .(n_files = .N), keyby = strata])
if (length(extra) > 0L) {
  cat("\nPresent ONLY in tables2agegrps/:\n")
  print(data.table(file = extra, strata = strata_of(extra))[, .(n = .N), keyby = strata])
}

# ---- (2) did two_agegrps actually reach each family? -----------------------
# Read one file per family that carries an agegrp/agegroup column in BOTH
# directories and compare the levels actually present.
cat("\n=== agegrp levels actually written, by family ===\n")
common <- intersect(f_std, f_two)
has_age <- grep("agegrp|agegroup", common, value = TRUE)

family_of <- function(f) {
  x <- sub("\\s*\\((not standardised|[^)]*standardised)\\)\\.csv$", "", f)
  x <- sub(" popdenom$", "", x)
  if (grepl("^all-cause mortality given disease-", x)) "all_cause_mrtl (task 2)"
  else if (grepl("^disease characteristics", x))      "disease_char (task 3)"
  else if (grepl("^exposures", x))                    "xps (task 4)"
  else                                                "main (task 1)"
}

levels_in <- function(dir, f) {
  d <- fread(file.path(dir, f), nrows = 500000L)
  col <- intersect(c("agegrp", "agegroup"), names(d))
  if (length(col) == 0L) return(NA_character_)
  paste(sort(unique(as.character(d[[col[[1L]]]]))), collapse = ", ")
}

res <- rbindlist(lapply(has_age, function(f) {
  data.table(family = family_of(f), file = f,
             std = levels_in(d_std, f), two = levels_in(d_two, f))
}))
# One representative file per family is enough; they share the collapse.
rep <- res[, .SD[1L], keyby = family]
for (i in seq_len(nrow(rep))) {
  cat(sprintf("\n%s\n  e.g. %s\n", rep$family[i], rep$file[i]))
  cat(sprintf("    tables/         : %s\n", rep$std[i]))
  cat(sprintf("    tables2agegrps/ : %s\n", rep$two[i]))
  cat(sprintf("    -> two_agegrps %s\n",
              if (identical(rep$std[i], rep$two[i]))
                "had NO EFFECT on this family" else "collapsed the bands"))
}

cat("\n=== per-family verdict across every comparable file ===\n")
print(res[, .(files = .N,
              collapsed = sum(std != two),
              unchanged = sum(std == two)), keyby = family])
