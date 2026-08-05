# Tests for .tbl_col_prefixes(), the metric -> output-column-prefix map used by
# private$tbl_smmrs_core() to name every quantile column it writes.
#
# The bug these pin: "ftlt_change_relative" mapped to "ftlt_rate_" -- byte for
# byte the prefix of the "ftlt" LEVEL table. Two things went wrong as a result.
# The published column claimed to be a rate while holding a ratio to the
# baseline year (a 2013 median of 0.72 means "72% of the 2026 value", not a case
# fatality rate of 0.72). And because the name contains neither "change" nor
# "rel", any downstream code selecting change columns by pattern skipped those
# six tables silently rather than erroring -- the failure mode that hides.
#
# The assertions below are deliberately about the CONVENTION, not just the one
# entry, so that a future family added by copy-paste trips the same wire.

library(tinytest)
library(IMPACTncdEngland)

prefixes <- IMPACTncdEngland:::.tbl_col_prefixes()

# --- the specific regression -------------------------------------------------
expect_equal(unname(prefixes[["ftlt_change_relative"]]), "ftlt_change_relative_")
expect_true(prefixes[["ftlt_change_relative"]] != prefixes[["ftlt"]])
expect_true(prefixes[["ftlt_change_relative"]] != prefixes[["ftlt_change_absolute"]])

# --- structural sanity -------------------------------------------------------
expect_true(is.character(prefixes))
expect_true(!is.null(names(prefixes)))
expect_equal(anyDuplicated(names(prefixes)), 0L)          # keys unique
expect_true(all(nzchar(prefixes)))
expect_true(all(grepl("_$", prefixes)))                   # all end in "_"

change_keys <- grep("_change_(relative|absolute)$", names(prefixes), value = TRUE)
level_keys  <- setdiff(names(prefixes), change_keys)
expect_true(length(change_keys) >= 12L)

# --- a change metric must ADVERTISE that it is a change ----------------------
# This is the property that makes `grep("change", names(dt))` a safe way to pick
# out change columns downstream. It is exactly what ftlt_rate_ violated.
for (k in change_keys) {
  expect_true(grepl("change", prefixes[[k]], fixed = TRUE), info = paste("change key:", k))
}

# --- and a level metric must NOT claim to be one -----------------------------
for (k in level_keys) {
  expect_false(grepl("change", prefixes[[k]], fixed = TRUE), info = paste("level key:", k))
}

# --- each change variant is distinct from its own level ----------------------
# A family's relative and absolute tables are separate FILES, so a collision is
# not a duplicate-column error -- it is a silently mislabelled table. Nothing
# else catches it.
for (k in change_keys) {
  fam <- sub("_change_(relative|absolute)$", "", k)
  if (fam %in% names(prefixes)) {
    expect_true(prefixes[[k]] != prefixes[[fam]],
                info = paste(k, "collides with its level", fam))
  }
}

# --- and relative differs from absolute within a family ----------------------
fams <- unique(sub("_change_(relative|absolute)$", "", change_keys))
for (f in fams) {
  rel <- paste0(f, "_change_relative")
  abs <- paste0(f, "_change_absolute")
  if (all(c(rel, abs) %in% names(prefixes))) {
    expect_true(prefixes[[rel]] != prefixes[[abs]],
                info = paste("relative/absolute collision in", f))
  }
}

# --- every family that has one change variant has both -----------------------
for (f in fams) {
  expect_true(paste0(f, "_change_relative") %in% names(prefixes), info = f)
  expect_true(paste0(f, "_change_absolute") %in% names(prefixes), info = f)
}
