# Tests for .promote_integer_measures().
#
# Every rate in this file is built by dividing one summed column by another and
# writing the quotient back with `:=`. data.table assigns IN PLACE, so a double
# written into an integer column is TRUNCATED: a case-fatality rate of 0.034
# silently becomes 0, with nothing but a warning to show for it. A whole table
# of zeros looks like a real result.
#
# tbl_smmrs_core() has always promoted its summed columns first;
# export_all_cause_mrtl_tables() did not, and its three division sites were
# protected only by the stock summaries happening to store counts as `double`.
# That is a property of whoever wrote the parquet, not of this code. Both now
# call the same helper.
#
# The decisive test is section 3: the SAME fixture, once integer and once
# double, must produce identical tables.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

promote <- IMPACTncdEngland:::.promote_integer_measures

# ===========================================================================
# 1. The helper: promotes measures, never keys, and reports what it did
# ===========================================================================
d <- data.table(year = 1:3, grp = c("a", "b", "c"),
                n_int = 1:3, n_dbl = c(1.5, 2.5, 3.5), lab = c("x", "y", "z"))
expect_true(promote(d, c("year", "grp")), info = "reports TRUE when it promoted")
expect_true(is.double(d$n_int), info = "an integer measure is promoted")
expect_true(is.double(d$n_dbl), info = "a double measure is left alone")
expect_true(is.integer(d$year), info = "an integer KEY is NOT promoted")
expect_true(is.character(d$lab), info = "a non-numeric measure is untouched")

expect_false(promote(d, c("year", "grp")),
             info = "FALSE on a second pass: nothing left to promote")
expect_false(promote(data.table(a = 1L), "a"),
             info = "FALSE when every column is a key")
expect_false(promote(data.table(a = 1L, b = 2.0), "a"),
             info = "FALSE when no measure is integer")

# Values survive the promotion exactly.
d2 <- data.table(k = 1L, v = 7L)
promote(d2, "k")
expect_equal(d2$v, 7, info = "promotion preserves the value")

# ===========================================================================
# 2. The truncation this prevents, demonstrated
# ===========================================================================
# It matters WHICH form of `:=` is used, and the difference is easy to miss:
#
#   DT[, col := <double>]        replaces the whole column ("plonk"), so the
#                                column is PROMOTED to double. Harmless.
#   DT[i, col := <double>]       a join-update writes into a SUBSET of the
#                                existing column, so it must coerce to that
#                                column's type. An integer column TRUNCATES.
#
# Both division sites in export_all_cause_mrtl_tables() are the second form,
# which is exactly why the guard is needed there.
den <- data.table(k = 1:2, denom = c(400L, 300L))

plain <- data.table(k = 1:2, value = c(4L, 6L), denom = c(400L, 300L))
plain[, value := value / denom]
expect_equal(plain$value, c(0.01, 0.02),
             info = "plain := promotes the column; no truncation")

bad <- data.table(k = 1:2, value = c(4L, 6L))
suppressWarnings(bad[den, on = "k", value := value / denom])
expect_true(is.integer(bad$value),
            info = "join-update keeps the column integer...")
expect_equal(bad$value, c(0L, 0L),
             info = "...so unguarded, a genuine 0.01 / 0.02 rate becomes zero")

good <- data.table(k = 1:2, value = c(4L, 6L))
promote(good, "k")
good[den, on = "k", value := value / denom]
expect_equal(good$value, c(0.01, 0.02),
             info = "guarded: the same join-update keeps the rate")

# ===========================================================================
# 3. End to end: an integer summary must equal a double one
# ===========================================================================
BANDS <- c("30-34", "65-69", "70-74")

mk <- function(as_int) {
  d <- CJ(mc = 1:3, scenario = "sc0", year = 19:20, agegrp = BANDS,
          sex = c("men", "women"), sorted = FALSE)
  k <- d$mc + 10L * (d$year - 19L) + 100L * match(d$agegrp, BANDS) +
    1000L * match(d$sex, c("men", "women"))
  deaths <- 4L + k %% 17L
  cases <- 400L + k %% 97L
  d[, deaths_chd := if (as_int) deaths else as.numeric(deaths)]
  d[, cases_chd  := if (as_int) cases  else as.numeric(cases)]
  d[]
}
mk_prvl <- function(as_int) {
  d <- CJ(mc = 1:3, scenario = "sc0", year = 19:20, agegrp = BANDS,
          sex = c("men", "women"), sorted = FALSE)
  p <- 5000L + d$mc + 100L * match(d$agegrp, BANDS)
  d[, popsize := if (as_int) p else as.numeric(p)]
  d[]
}

run_acm <- function(as_int) {
  dir <- tempfile("intp_"); dir.create(dir, recursive = TRUE)
  src <- list(all_cause_mrtl_by_dis = mk(as_int), prvl = mk_prvl(as_int))
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  env$private <- list(read_summary_dataset = function(summary_type, standardization) {
    s <- src[[summary_type]]
    if (is.null(s) || identical(standardization, "esp")) NULL else copy(s)
  })
  f <- Simulation$private_methods$export_all_cause_mrtl_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975, 0.1, 0.9), summaries_dir = dir, tables_dir = dir,
    two_agegrps = FALSE, strata_ons = list(c("year", "agegrp")),
    strata_esp = list("year"))
  dir
}

d_int <- run_acm(TRUE)
d_dbl <- run_acm(FALSE)
files <- sort(list.files(d_int, pattern = "\\.csv$"))
expect_equal(files, sort(list.files(d_dbl, pattern = "\\.csv$")),
             info = "both runs write the same files")
expect_true(length(files) >= 2L,
            info = "both the disease-denominator and popdenom tables are covered")

MED <- "all_cause_mrtl_by_disease_rate_50.0%"
for (f in files) {
  a <- fread(file.path(d_int, f))
  b <- fread(file.path(d_dbl, f))
  expect_equal(a, b, info = paste0("integer summary == double summary: ", f))
  # ...and the rates are real, not a table of zeros that happens to match.
  expect_true(any(a[[MED]] > 0),
              info = paste0("rates are non-zero (guard actually fired): ", f))
}
