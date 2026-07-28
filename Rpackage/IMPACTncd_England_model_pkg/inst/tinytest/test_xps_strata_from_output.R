# Tests for xps_strata_from_output(), which derives the exposure-table
# stratification from the design's `strata_for_output`.
#
# Regression: the mapping (agegrp -> agegrp20, dimd -> qimd) is many-to-one, so
# a design naming BOTH `dimd` and `qimd` produced c(..., "qimd", "qimd"), and
# data.table's groupingsets() rejects a duplicated `by` with "Argument 'by' must
# have unique column names for grouping" -- failing every parallel task of
# run() with an opaque message.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

f <- IMPACTncdEngland:::xps_strata_from_output

# --- The stock design ------------------------------------------------------
expect_equal(f(c("scenario", "year", "agegrp", "sex", "dimd")),
             c("agegrp20", "sex", "qimd"),
             info = "agegrp -> agegrp20, dimd -> qimd, scenario/year dropped")

# --- Both sides of a mapped pair: the regression --------------------------
expect_equal(f(c("scenario", "year", "agegrp", "sex", "dimd", "qimd")),
             c("agegrp20", "sex", "qimd"),
             info = "dimd AND qimd collapse to a single qimd")
expect_equal(f(c("scenario", "year", "agegrp", "agegrp20", "sex")),
             c("agegrp20", "sex"),
             info = "agegrp AND agegrp20 collapse to a single agegrp20")
expect_equal(f(c("year", "dimd", "qimd", "agegrp", "agegrp20")),
             c("qimd", "agegrp20"),
             info = "both pairs collapsed at once")

# The result must be usable as a groupingsets() `by`, which is the thing that
# actually broke.
by <- c("year", f(c("scenario", "year", "agegrp", "sex", "dimd", "qimd")))
expect_false(anyDuplicated(by) > 0L, info = "no duplicated grouping columns")
dt <- data.table(year = 1:4, agegrp20 = 1:4, sex = 1:4, qimd = 1:4,
                 x = as.numeric(1:4))
expect_silent(groupingsets(dt, j = lapply(.SD, mean), by = by,
                           .SDcols = "x", sets = list(by)))

# --- Order and pass-through ------------------------------------------------
expect_equal(f(c("year", "sex", "agegrp", "dimd")),
             c("sex", "agegrp20", "qimd"),
             info = "first-occurrence order is preserved")
expect_equal(f(c("year", "ethnicity", "sha")), c("ethnicity", "sha"),
             info = "unmapped variables pass through unchanged")
expect_equal(f(c("scenario", "year")), character(0),
             info = "scenario and year alone leave nothing to stratify by")
expect_equal(f(character(0)), character(0))
