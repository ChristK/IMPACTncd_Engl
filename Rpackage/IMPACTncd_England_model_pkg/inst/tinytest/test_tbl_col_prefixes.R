# Tests for .tbl_col_prefixes(), the metric -> output-column-prefix map used by
# private$tbl_smmrs_core() to name every quantile column it writes.
#
# Two bugs are pinned here, both of the same species: a column name asserting
# something the values do not support.
#
#  1. "ftlt_change_relative" mapped to "ftlt_rate_" -- byte for byte the prefix
#     of the "ftlt" LEVEL table. The published column claimed to be a rate while
#     holding a ratio to the baseline year (a 2013 median of 0.72 means "72% of
#     the 2026 value", not a case fatality rate of 0.72). And because the name
#     contains neither "change" nor "rel", any downstream code selecting change
#     columns by pattern skipped those six tables silently rather than erroring
#     -- the failure mode that hides.
#
#  2. "prvl_change_relative" and "incd_change_relative" both mapped to
#     "prct_change_relative_", asserting a percentage where the value is a ratio
#     (1.0 at baseline). Taken at face value it understates every effect by a
#     factor of 100 without ever looking wrong.
#
# The assertions below are deliberately about the CONVENTION, not just the one
# entry, so that a future family added by copy-paste trips the same wire.

library(tinytest)
library(IMPACTncdEngland)

prefixes <- IMPACTncdEngland:::.tbl_col_prefixes()

# --- regression 1: ftlt_change_relative reused its own LEVEL prefix ----------
expect_equal(unname(prefixes[["ftlt_change_relative"]]), "ftlt_change_relative_")
expect_true(prefixes[["ftlt_change_relative"]] != prefixes[["ftlt"]])
expect_true(prefixes[["ftlt_change_relative"]] != prefixes[["ftlt_change_absolute"]])

# --- regression 2: prvl/incd relative change claimed to be a percentage ------
# Both mapped to "prct_change_relative_". The value is `value / i.value` -- a
# RATIO to the baseline year (1.0 at baseline, 0.9 for a 10% reduction), not a
# percentage. Reading "prct_" literally turns a 10% fall into a 0.9% fall: an
# order-of-magnitude error that looks perfectly plausible in a results table,
# which is precisely why it survived. The shared prefix also meant the column
# could not say which family it came from.
expect_equal(unname(prefixes[["prvl_change_relative"]]), "prvl_change_relative_")
expect_equal(unname(prefixes[["incd_change_relative"]]), "incd_change_relative_")
expect_true(prefixes[["prvl_change_relative"]] != prefixes[["incd_change_relative"]])

# No prefix may claim to be a percentage while holding a ratio. Nothing in this
# map is a percentage, so the token must not appear at all.
expect_false(any(grepl("prct", prefixes, fixed = TRUE)))
expect_false(any(grepl("percent", prefixes, fixed = TRUE)))

# --- every RELATIVE prefix is unique across families -------------------------
# Deliberately relative-only. "abs_change_" is still shared by prvl and incd,
# and that is fine: an absolute change is exactly what the name claims, so it is
# generic rather than wrong. Renaming it would break consumers for no
# correctness gain -- so this invariant must NOT be widened to all change keys.
rel_keys <- grep("_change_relative$", names(prefixes), value = TRUE)
rel_prefixes <- unname(prefixes[rel_keys])
expect_equal(anyDuplicated(rel_prefixes), 0L)
expect_equal(length(rel_prefixes), length(unique(rel_prefixes)))

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
