# Tests for the two-age-group collapse and its configurable split age.
#
# `export_tables(two_agegrps = TRUE)` reports `agegrp` as two coarse groups and
# writes to `tables2agegrps/`. Two things are pinned here:
#
#   A. The collapse REACHES every family whose summaries carry `agegrp`. It used
#      to be plumbed only into the `main` and `equity` tasks, so
#      export_all_cause_mrtl_tables() silently produced 5-year bands inside
#      `tables2agegrps/` -- files byte-identical to their `tables/` namesakes.
#      A silent no-op, which is exactly the failure mode worth a regression test.
#      `disease_char` has no `agegrp` column and `xps` uses `agegrp20` (20-year
#      bands), so both are correctly out of scope.
#
#   B. The split age is an ARGUMENT (`two_agegrps_split_age`, default 65L =
#      the historical 30-64 / 65-99), the labels are derived from the DATA
#      rather than hard-coded 30..99, and a split that would cut a band in half
#      is rejected in the parent rather than silently mis-binning.
#
# The load-bearing correctness property throughout: the collapse is a
# RELABELLING done before the group-by, so a collapsed count must equal the sum
# of the constituent 5-year bands. Every "it worked" assertion below is that
# identity, not a spot check on a label.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

# Internal helpers, reached the same way test_qimd_derivation.R reaches
# add_qimd_from_dimd().
two_agegrp_map <- IMPACTncdEngland:::two_agegrp_map
collapse_two_agegrps <- IMPACTncdEngland:::collapse_two_agegrps
validate_two_agegrps_split_age <- IMPACTncdEngland:::validate_two_agegrps_split_age

BANDS <- c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
           "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99")

# ===========================================================================
# 1. two_agegrp_map(): labels, boundaries, and loud failure
# ===========================================================================
m <- two_agegrp_map(BANDS, 65L)
expect_equal(unname(m[c("30-34", "60-64")]), c("30-64", "30-64"),
             info = "bands below the split take the lower label")
expect_equal(unname(m[c("65-69", "95-99")]), c("65-99", "65-99"),
             info = "bands from the split up take the upper label")
expect_equal(attr(m, "levels_ordered"), c("30-64", "65-99"),
             info = "default split reproduces the historical labels")

# The labels come from the data, not a hard-coded 30..99.
m2 <- two_agegrp_map(c("20-24", "25-29", "30-34", "35-39"), 30L)
expect_equal(attr(m2, "levels_ordered"), c("20-29", "30-39"),
             info = "bounds are read off the bands actually present")

# Every legal split, and only those, is accepted.
for (s in c(35L, 50L, 65L, 95L)) {
  expect_silent(two_agegrp_map(BANDS, s))
}
expect_error(two_agegrp_map(BANDS, 62L), pattern = "falls INSIDE",
             info = "a split mid-band is refused, not silently rounded")
# Degenerate splits leave one group empty -- always a mistake under an argument
# called `two_agegrps`, and the empty group's label would be undefined.
expect_error(two_agegrp_map(BANDS, 100L), pattern = "EMPTY",
             info = "a split past the top band is refused")
expect_error(two_agegrp_map(BANDS, 30L), pattern = "EMPTY",
             info = "a split at the lowest band start is refused")

# The split need NOT start a band that is actually present: what matters is
# that it falls between bands. A locality with an empty 5-year band around the
# split must still be splittable there.
sparse <- two_agegrp_map(c("30-34", "70-74"), 65L)
expect_equal(attr(sparse, "levels_ordered"), c("30-64", "65-74"),
             info = "a gap at the split is fine; labels come from the bin edges")

# An open-ended top band keeps its form rather than inventing an upper bound.
expect_equal(attr(two_agegrp_map(c("30-34", "65-69", "90+"), 65L),
                  "levels_ordered"),
             c("30-64", "65+"),
             info = "an open-ended band yields an open-ended upper label")

# Unparseable labels fail loudly (cf. deprivation_rank).
expect_error(two_agegrp_map(c("30-34", "All"), 65L), pattern = "unparseable",
             info = "a non-band label errors instead of landing in a group")

# ===========================================================================
# 2. collapse_two_agegrps(): in place, idempotent, and a safe no-op
# ===========================================================================
d <- data.table(agegrp = BANDS, n = 1L)
expect_true(collapse_two_agegrps(d, 65L))
expect_equal(sort(unique(d$agegrp)), c("30-64", "65-99"))
expect_equal(d[, sum(n), keyby = agegrp]$V1, c(7L, 7L),
             info = "seven 5-year bands land in each group")
before <- copy(d)
collapse_two_agegrps(d, 65L)
expect_equal(d, before, info = "idempotent: re-collapsing changes nothing")

expect_false(collapse_two_agegrps(data.table(x = 1), 65L),
             info = "no-op (FALSE) when there is no agegrp column")
expect_false(collapse_two_agegrps(NULL, 65L), info = "no-op on NULL")

# A factor column keeps its type and gains ordered levels.
df <- data.table(agegrp = factor(BANDS, levels = BANDS))
collapse_two_agegrps(df, 65L)
expect_true(is.factor(df$agegrp), info = "a factor stays a factor")
expect_equal(levels(df$agegrp), c("30-64", "65-99"),
             info = "and its levels are the two groups, in order")

# ===========================================================================
# 3. validate_two_agegrps_split_age(): structural, data-free
# ===========================================================================
expect_equal(validate_two_agegrps_split_age(65), 65L)
expect_equal(validate_two_agegrps_split_age(65L), 65L)
for (bad in list("65", c(60, 65), NA_integer_, -5L, 0L, 62.5, NULL, Inf)) {
  expect_error(validate_two_agegrps_split_age(bad),
               pattern = "single positive whole number")
}
# The default in the public signature is the historical split.
expect_equal(eval(formals(Simulation$public_methods$export_tables)$two_agegrps_split_age),
             65L, info = "export_tables() default split age is 65")

# ===========================================================================
# 4. The collapse reaches all_cause_mrtl -- the family it used to miss
# ===========================================================================
# The production quantiles. Worth using verbatim: `scales::percent()` infers its
# accuracy from the whole vector, so a shorter `prbl` would name the median
# column "...50%" instead of the published "...50.0%".
PRBL <- c(0.5, 0.025, 0.975, 0.1, 0.9)
MED <- "all_cause_mrtl_by_disease_rate_50.0%"

# Counts are stored as DOUBLE in the summaries (verified against
# all_cause_mrtl_by_dis_scaled_up). That matters: export_all_cause_mrtl_tables()
# has no integer -> numeric promotion of its own (unlike tbl_smmrs_core), so an
# integer fixture would truncate every rate to 0 and make section 5 vacuous.
mk_acm <- function() {
  d <- CJ(mc = 1:3, scenario = "sc0", year = 19:20, agegrp = BANDS,
          sex = c("men", "women"), dimd = c("1 most deprived", "10 least deprived"),
          sorted = FALSE)
  k <- d$mc + 10L * (d$year - 19L) + 100L * match(d$agegrp, BANDS) +
    1000L * match(d$sex, c("men", "women"))
  d[, deaths_chd := as.numeric(4L + k %% 17L)]
  d[, cases_chd  := as.numeric(400L + k %% 97L)]
  d[]
}
mk_prvl <- function() {
  d <- CJ(mc = 1:3, scenario = "sc0", year = 19:20, agegrp = BANDS,
          sex = c("men", "women"), dimd = c("1 most deprived", "10 least deprived"),
          sorted = FALSE)
  d[, popsize := as.numeric(5000L + d$mc + 100L * match(d$agegrp, BANDS))]
  d[]
}

run_acm <- function(two_agegrps, split_age = 65L, strata = list(c("year", "agegrp"))) {
  dir <- tempfile("acm_"); dir.create(dir, recursive = TRUE)
  src <- list(all_cause_mrtl_by_dis = mk_acm(), prvl = mk_prvl())
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  env$private <- list(read_summary_dataset = function(summary_type, standardization) {
    s <- src[[summary_type]]
    if (is.null(s) || identical(standardization, "esp")) NULL else copy(s)
  })
  f <- Simulation$private_methods$export_all_cause_mrtl_tables
  environment(f) <- env
  f(prbl = PRBL, summaries_dir = dir, tables_dir = dir,
    two_agegrps = two_agegrps, two_agegrps_split_age = split_age,
    strata_ons = strata, strata_esp = list("year"))
  dir
}

acm_file <- "all-cause mortality given disease-year-agegroup (not standardised).csv"
a_std <- fread(file.path(run_acm(FALSE), acm_file))
a_two <- fread(file.path(run_acm(TRUE), acm_file))

expect_equal(sort(unique(a_std$agegrp)), BANDS,
             info = "two_agegrps = FALSE keeps the 5-year bands")
expect_equal(sort(unique(a_two$agegrp)), c("30-64", "65-99"),
             info = "REGRESSION: two_agegrps now reaches all_cause_mrtl")
expect_false(isTRUE(all.equal(a_std, a_two)),
             info = "...and the tables genuinely differ (it was a silent no-op)")

# A non-default split is honoured end to end.
a_50 <- fread(file.path(run_acm(TRUE, 50L), acm_file))
expect_equal(sort(unique(a_50$agegrp)), c("30-49", "50-99"),
             info = "two_agegrps_split_age is plumbed through to the family")

# ===========================================================================
# 5. Correctness: collapsing == aggregating the constituent bands
# ===========================================================================
# The all-cause-mortality table is a RATE (deaths/cases), so the identity to
# check is summed-numerator / summed-denominator, computed here directly from
# the fixture and compared with what the collapsed export actually wrote.
raw <- mk_acm()[, year := year + 2000L]
raw[, grp := fifelse(match(agegrp, BANDS) <= 7L, "30-64", "65-99")]
expect_by_hand <- raw[, .(rate = sum(deaths_chd) / sum(cases_chd)),
                      keyby = .(mc, year, grp)][
                        , .(med = median(rate)), keyby = .(year, grp)]
got <- a_two[disease == "chd"][order(year, agegrp)]
expect_true(MED %in% names(got), info = "the median column is named as published")
expect_equal(got[[MED]], expect_by_hand$med, tolerance = 1e-6,
             info = "collapsed rate == summed deaths / summed cases over bands")

# The same identity must NOT hold for the uncollapsed table -- otherwise the
# check above would pass even if the collapse had done nothing.
expect_false(isTRUE(all.equal(
  a_std[disease == "chd"][order(year, agegrp)][[MED]], expect_by_hand$med)),
  info = "the 5-year-band table is a different quantity, as it must be")

# ===========================================================================
# 6. The CEA family is plumbed too
# ===========================================================================
expect_true("two_agegrps" %in%
              names(formals(Simulation$private_methods$export_cea_tables)),
            info = "export_cea_tables takes two_agegrps")
expect_true("two_agegrps_split_age" %in%
              names(formals(Simulation$private_methods$export_cea_tables)),
            info = "...and the split age")
for (f in c("export_main_tables", "export_all_cause_mrtl_tables",
            "export_cea_tables", "export_equity_tables")) {
  fm <- formals(Simulation$private_methods[[f]])
  expect_equal(eval(fm$two_agegrps_split_age), 65L,
               info = paste0(f, ": split-age default agrees with export_tables()"))
}
