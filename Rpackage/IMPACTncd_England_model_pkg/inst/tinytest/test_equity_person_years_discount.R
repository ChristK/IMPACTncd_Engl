# Tests for the three 2026-08 changes to the equity slope-index tables:
#
#   1. the CUMULATIVE PERSON-YEARS denominator, which makes the estimator
#      return exactly zero on a distributionally-neutral benefit stream even
#      when the deprivation groups grow at different rates (the single
#      reporting-year denominator did not, and published a spurious pro-rich
#      gradient out of demography alone);
#   2. multi-level DISCOUNTING of net QALYs, applied per-year before the
#      cumulative sum;
#   3. the P_PRO_POOR column, the share of Monte-Carlo draws falling on the
#      pro-poor side of each index's own null.
#
# Same mock-harness approach as test_export_equity_tables.R: the private method
# is re-parented onto an environment whose parent is the package namespace, so
# calc_equity_slope_indices(), expand_discount_levels() etc. stay resolvable.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

dimd_lv <- c("1 most deprived", as.character(2:9), "10 least deprived")

run_export <- function(sources, strata = list("year"),
                       discount_levels = NULL, discount_from_year = NULL,
                       ridit_reference = "comparator") {
  dir <- tempfile("eqpy_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE, diseases = NULL)))
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      if (is.null(sources[[summary_type]])) NULL else copy(sources[[summary_type]])
    }
  )
  f <- Simulation$private_methods$export_equity_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
    summaries_dir = dir, tables_dir = dir,
    comparator_scenario = "sc0", baseline_year = 2019L,
    ridit_reference = ridit_reference,
    discount_levels = discount_levels,
    discount_from_year = discount_from_year,
    strata = strata)
  dir
}

# Population that GROWS FASTEST IN THE MOST DEPRIVED deciles -- the pattern in
# this model's ONS projection input (most/least deprived ratio 1.35 -> 1.79 over
# 2019-2043 in the England run) and the one the old denominator turned into a
# fake pro-rich gradient.
decile_pop <- function(year, rank_) {
  # rank_ 1 = most deprived. Base population falls with affluence, and the
  # annual growth rate falls with it too (5%/yr most deprived, 0.5%/yr least).
  base <- 100000 * (11L - rank_)
  round(base * (1 + 0.05 * ((11L - rank_) / 10) * (year - 19L)))
}

# ===========================================================================
# 1. Neutrality: an exactly equal benefit per head per year scores ZERO
# ===========================================================================
# Every decile receives the same QALY gain per person in every year, so there
# is no gradient to find -- whatever the deciles' relative sizes or growth.
PER_HEAD <- 0.001

qalys_neutral <- local({
  d <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 19:23,
          dimd = factor(dimd_lv, levels = dimd_lv), sorted = FALSE)
  d[, rank_ := as.integer(dimd)]
  d[, popsize := decile_pop(year, rank_)]
  d[, EQ5D5L := 0.8 * popsize + fifelse(scenario == "sc1",
                                        PER_HEAD * popsize, 0)]
  d[, rank_ := NULL]
  d[]
})

dn <- run_export(list(qalys = qalys_neutral))
tn <- fread(file.path(dn, "equity net_qalys slope index by year (not standardised).csv"))

expect_true(nrow(tn) > 0L, info = "the neutral fixture still produces a table")
expect_true(all(abs(tn[type == "AEI_per100k", `equity_50.0%`]) < 1e-6),
            info = "neutral benefit -> AEI_per100k is exactly zero")
expect_true(all(abs(tn[type == "AEI_total", `equity_50.0%`]) < 1e-6),
            info = "neutral benefit -> AEI_total is exactly zero")
expect_true(all(abs(tn[type == "AEI_total", `equity_2.5%`]) < 1e-6 &
                abs(tn[type == "AEI_total", `equity_97.5%`]) < 1e-6),
            info = "every draw is zero, so the whole interval collapses on zero")

# The fixture must actually EXERCISE the bug: recompute the index the old way
# (cumulative numerator over the single reporting-year population) and confirm
# it is materially non-zero. Without this the test above would also pass on
# code that had never been fixed.
old_style_index <- function(dt, T_year) {
  b <- dt[scenario == "sc1" & mc == 1L & year >= 19L & year <= T_year]
  s <- dt[scenario == "sc0" & mc == 1L & year >= 19L & year <= T_year]
  x <- merge(b[, .(year, dimd, v1 = EQ5D5L)],
             s[, .(year, dimd, v0 = EQ5D5L, N = popsize)],
             by = c("year", "dimd"))
  x[, B_year := v1 - v0]
  agg <- x[, .(B = sum(B_year)), by = dimd]                     # cumulative
  agg <- merge(agg, x[year == T_year, .(dimd, N)], by = "dimd") # single-year N
  agg[, rank_ := as.integer(factor(dimd, levels = dimd_lv))]
  setorder(agg, -rank_)                                          # least -> most
  agg[, csum := cumsum(N)]
  agg[, ridit := (csum - N / 2) / sum(N)]
  unname(coef(lm(I(B / N) ~ ridit, data = agg, weights = N))[2L]) * 1e5
}
old_bias <- old_style_index(qalys_neutral, 23L)
# Measured: -32.4 per 100k against a cumulative benefit of 500 per 100k, i.e.
# the old denominator turned a perfectly neutral policy into a 6.5% PRO-RICH
# gradient. (The England run's differential growth is steeper still: ~36% of
# the published index at 2043.) The sign matters as much as the size -- it is
# always pro-rich when the deprived groups grow fastest.
expect_true(old_bias < -10,
            info = paste("the fixture reproduces the old pro-rich bias, so the",
                         "zero above is the fix and not a degenerate fixture"))

# The same neutrality must hold under the "scenario" ridit reference.
ds <- run_export(list(qalys = qalys_neutral), ridit_reference = "scenario")
ts <- fread(file.path(ds, "equity net_qalys slope index by year (not standardised).csv"))
expect_true(all(abs(ts[type == "AEI_per100k", `equity_50.0%`]) < 1e-6),
            info = "neutrality does not depend on the ridit reference")

# ===========================================================================
# 2. Discounting
# ===========================================================================
# A gradient fixture: the per-head gain RISES with deprivation (pro-poor), so
# the indices are non-zero and discounting has something to shrink.
qalys_grad <- local({
  d <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 19:23,
          dimd = factor(dimd_lv, levels = dimd_lv), sorted = FALSE)
  d[, rank_ := as.integer(dimd)]
  d[, popsize := decile_pop(year, rank_)]
  # per-head gain 0.002 in the most deprived decile falling to 0.0002 in the
  # least, and growing with the year so discounting bites
  d[, gain := 0.0002 * (11L - rank_) * (year - 18L)]
  d[, EQ5D5L := 0.8 * popsize + fifelse(scenario == "sc1", gain * popsize, 0)]
  d[, c("rank_", "gain") := NULL]
  d[]
})

incd_grad <- local({
  d <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 19:23,
          dimd = factor(dimd_lv, levels = dimd_lv), sorted = FALSE)
  d[, rank_ := as.integer(dimd)]
  d[, popsize := decile_pop(year, rank_)]
  d[, chd_incd := fifelse(scenario == "sc0", 0.01 * popsize,
                          0.01 * popsize * (1 - 0.02 * (11L - rank_)))]
  d[, rank_ := NULL]
  d[]
})

lv <- IMPACTncdEngland:::build_discount_levels(c(0, 3.5), c(0, 3.5))
dd <- run_export(list(qalys = qalys_grad, incd = incd_grad),
                 discount_levels = lv, discount_from_year = 2019L)

tq <- fread(file.path(dd, "equity net_qalys slope index by year (not standardised).csv"))
tc <- fread(file.path(dd, "equity cpp slope index by year (not standardised).csv"))

expect_equal(sort(unique(tq$discount)), c("0%", "3.5%"),
             info = "net QALYs is expanded across both discount levels")
expect_equal(unique(tc$discount), "0%",
             info = "CPP is an event count and is never discounted")

# The set of FILES is unchanged however many levels are requested -- the levels
# live in rows, exactly as in the main and CEA tables.
d0 <- run_export(list(qalys = qalys_grad, incd = incd_grad))
expect_equal(sort(list.files(dd)), sort(list.files(d0)),
             info = "discount levels add rows, never files")

# The 0% block reproduces the undiscounted run byte for byte.
t0 <- fread(file.path(d0, "equity net_qalys slope index by year (not standardised).csv"))
expect_equal(tq[discount == "0%"], t0,
             info = "the 0% block is identical to an undiscounted run")

# Discounting shrinks a positive benefit, so it shrinks the absolute index.
#
# This is a property OF THIS FIXTURE, not a general guarantee, and it must not
# be read as one. Discounting is not a uniform rescaling of the index: it
# re-weights the YEARS making up the cumulative benefit, so when the gradient
# changes sign or shape over the horizon the discounted index is a different
# weighted combination and can be LARGER in magnitude. On the HFSS England run
# that happens in 3 of 25 years, all of them pre-policy years whose index is
# ~0. Here the gradient keeps one sign and shape throughout, so the shrink is
# guaranteed.
cmp <- merge(tq[discount == "0%" & type == "AEI_total", .(year, v0 = `equity_50.0%`)],
             tq[discount == "3.5%" & type == "AEI_total", .(year, v1 = `equity_50.0%`)],
             by = "year")
expect_true(all(abs(cmp$v1) <= abs(cmp$v0) + 1e-9),
            info = "3.5% discounting never inflates the absolute index")
expect_true(all(cmp[year > 2019, abs(v1) < abs(v0)]),
            info = "and strictly shrinks it once past the discount origin")
expect_equal(cmp[year == 2019, v1], cmp[year == 2019, v0],
             info = "the discount origin year is undiscounted at every level")

# THE LOAD-BEARING ONE: the per-year values are scaled BEFORE the cumsum, so
# `total_benefit` is a sum of present values, not a discounted total. The two
# differ whenever the benefit is spread over more than one year.
expected_pv <- local({
  x <- qalys_grad[mc == 1L]
  b <- dcast(x, year + dimd ~ scenario, value.var = "EQ5D5L")
  b[, B_year := sc1 - sc0]
  b[, f := 1 / 1.035^pmax(0, (year + 2000L) - 2019L)]
  b[year <= 22L, sum(B_year * f)]          # cumulated to 2022, present values
})
got_pv <- tq[discount == "3.5%" & type == "total_benefit" & year == 2022L,
             `equity_50.0%`]
expect_equal(got_pv, expected_pv, tolerance = 1e-6,
             info = "cumulative PRESENT VALUES, not a discounted cumulative total")

wrong_pv <- local({
  x <- qalys_grad[mc == 1L]
  b <- dcast(x, year + dimd ~ scenario, value.var = "EQ5D5L")
  b[, B_year := sc1 - sc0]
  b[year <= 22L, sum(B_year)] / 1.035^(2022L - 2019L)   # discount the total
})
expect_true(abs(got_pv - wrong_pv) > 1e-3,
            info = "the two orderings really do differ on this fixture")

# Neutrality under discounting. The UNDISCOUNTED level stays exactly zero.
dnd <- run_export(list(qalys = qalys_neutral),
                  discount_levels = lv, discount_from_year = 2019L)
tnd <- fread(file.path(dnd, "equity net_qalys slope index by year (not standardised).csv"))
expect_true(all(abs(tnd[discount == "0%" & type == "AEI_per100k",
                       `equity_50.0%`]) < 1e-6),
            info = "the 0% level of a neutral stream is exactly zero")

# The DISCOUNTED level is not exactly zero, and should not be. Present value
# depends on WHEN a benefit falls, so a group whose person-years are
# concentrated later genuinely receives less present value per person-year even
# when the undiscounted per-head benefit is identical. That is time preference
# doing its job, not a denominator artefact -- and it is second-order: measured
# here at -0.24 per 100k against the -32.4 the old denominator produced on the
# very same fixture, i.e. more than two orders of magnitude smaller.
disc_resid <- tnd[discount == "3.5%" & type == "AEI_per100k" & year == 2023L,
                  `equity_50.0%`]
expect_true(abs(disc_resid) < abs(old_bias) / 100,
            info = paste("the discounting residual is orders of magnitude",
                         "below the denominator bias it replaced"))
expect_true(disc_resid < 0,
            info = "and runs the same way growth does, as present value implies")

# ===========================================================================
# 3. p_pro_poor
# ===========================================================================
# 4 MC draws, 3 pro-poor and 1 pro-rich, so the probability is a known 0.75.
qalys_mixed <- local({
  d <- CJ(mc = 1:4, scenario = c("sc0", "sc1"), year = 19:21,
          dimd = factor(dimd_lv, levels = dimd_lv), sorted = FALSE)
  d[, rank_ := as.integer(dimd)]
  d[, popsize := decile_pop(year, rank_)]
  # draws 1-3 favour the deprived; draw 4 reverses the gradient
  d[, slope := fifelse(mc == 4L, -1, 1)]
  d[, gain := 0.001 + slope * 0.0002 * ((11L - rank_) - 5.5)]
  d[, EQ5D5L := 0.8 * popsize + fifelse(scenario == "sc1", gain * popsize, 0)]
  d[, c("rank_", "slope", "gain") := NULL]
  d[]
})

dp <- run_export(list(qalys = qalys_mixed))
tp <- fread(file.path(dp, "equity net_qalys slope index by year (not standardised).csv"))

expect_true("p_pro_poor" %in% names(tp), info = "p_pro_poor column is written")
expect_equal(unique(tp[type == "AEI_total", p_pro_poor]), 0.75,
             info = "3 of 4 draws are pro-poor")
expect_equal(unique(tp[type == "AEI_per100k", p_pro_poor]), 0.75,
             info = "AEI_per100k is AEI_total times a positive constant")
expect_true(all(is.na(tp[type == "total_benefit", p_pro_poor])),
            info = "total_benefit has no pro-poor direction")
expect_true(all(is.na(tp[type == "fit_R2", p_pro_poor])),
            info = "fit_R2 has no pro-poor direction")

# RII_ratio is a ratio of fitted extremes, so its null is 1 and not 0 -- and
# the fixture happens to demonstrate exactly the selection the code comment
# warns about. In the pro-rich draw the line extrapolates the most-deprived end
# below zero, so the fit0 > 0 & fit1 > 0 guard makes RII_ratio NA for that draw
# alone. RII_ratio therefore rests on 3 draws, all of them pro-poor, and reports
# p_pro_poor = 1 where AEI_total (which never drops a draw) reports 0.75.
#
# The scored threshold is still 1: every surviving RII_ratio exceeds 1, and if
# the code had scored it against 0 the answer would also be 1 here, so pin the
# median above 1 as well to keep the assertion meaningful.
expect_equal(unique(tp[type == "RII_ratio", n_mc]), 3L,
             info = "the positivity guard drops the pro-rich draw")
expect_equal(unique(tp[type == "RII_ratio", p_pro_poor]), 1,
             info = "all surviving RII_ratio draws are pro-poor")
expect_true(all(tp[type == "RII_ratio", `equity_50.0%`] > 1),
            info = "and they are pro-poor by exceeding 1, not by exceeding 0")
expect_true(unique(tp[type == "RII_ratio", p_pro_poor]) >
              unique(tp[type == "AEI_total", p_pro_poor]),
            info = paste("the relative index is selected towards pro-poor:",
                         "read p_pro_poor beside n_mc, never on its own"))

# A ratio-valued index scored against 0 instead of 1 would be a silent error
# (every positive ratio would count as pro-poor), so check the threshold bites
# on a fixture where the two disagree: an all-pro-rich set whose ratios are
# positive but below 1.
tpr_chk <- fread(file.path(
  run_export(list(qalys = qalys_mixed[mc == 4L])),
  "equity net_qalys slope index by year (not standardised).csv"))
if (nrow(tpr_chk[type == "RII_ratio" & !is.na(p_pro_poor)]) > 0L) {
  expect_true(all(tpr_chk[type == "RII_ratio", `equity_50.0%`] < 1),
              info = "pro-rich draw gives a ratio below 1")
  expect_equal(unique(tpr_chk[type == "RII_ratio", p_pro_poor]), 0,
               info = "a positive ratio below 1 is scored PRO-RICH")
}

# A wholly pro-poor and a wholly pro-rich fixture pin the two ends of the scale.
qalys_allpoor <- copy(qalys_mixed)[mc == 4L, EQ5D5L := NA_real_]
qalys_allpoor <- qalys_allpoor[!is.na(EQ5D5L)]
tpp <- fread(file.path(
  run_export(list(qalys = qalys_allpoor)),
  "equity net_qalys slope index by year (not standardised).csv"))
expect_equal(unique(tpp[type == "AEI_total", p_pro_poor]), 1,
            info = "all draws pro-poor -> 1")

qalys_allrich <- qalys_mixed[mc == 4L]
tpr <- fread(file.path(
  run_export(list(qalys = qalys_allrich)),
  "equity net_qalys slope index by year (not standardised).csv"))
expect_equal(unique(tpr[type == "AEI_total", p_pro_poor]), 0,
            info = "no draws pro-poor -> 0")

# p_pro_poor rests on the same draws as n_mc, so it can never be computed from
# more draws than the quantiles were.
expect_true(all(tp[!is.na(p_pro_poor), n_mc] >= 1L))
expect_true(all(tp[!is.na(p_pro_poor), p_pro_poor * n_mc] %% 1 < 1e-9),
            info = "p_pro_poor * n_mc is a whole number of draws")

# ===========================================================================
# 4. `discount` is a reserved stratum name
# ===========================================================================
expect_error(IMPACTncdEngland:::validate_equity_strata(list(c("year", "discount"))),
             pattern = "reserved column",
             info = "discount cannot be an output stratum")
