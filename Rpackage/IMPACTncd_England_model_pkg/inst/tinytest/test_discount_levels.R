# Tests for multi-level discounting in the table exports.
#
# `qaly_discount_rate` / `cost_discount_rate` take VECTORS (default c(0, 3.5)),
# paired element-wise into discount levels. Every level is reported in the
# qalys, net_qalys, costs, net_costs AND cost-effectiveness tables, tagged in a
# `discount` column - discounting is no longer confined to the CEA tables.
#
# Coverage note: before this file, no test exercised the qalys / costs /
# net_qalys / net_costs branches of tbl_smmrs_core() at all. These tests drive
# all four through the REAL code path on a synthetic summary whose discounted
# values are known in closed form, so the arithmetic is pinned and not merely
# the schema.
#
# Same mock harness as test_export_main_tables_working_copy.R: the methods'
# only couplings to the R6 object are `self$design$sim_prm$logs` and
# `private$read_summary_dataset()`.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

bdl <- IMPACTncdEngland:::build_discount_levels
dfac <- IMPACTncdEngland:::discount_factor
expand_lv <- IMPACTncdEngland:::expand_discount_levels

BASE <- 2019L
YEARS <- 19:21              # summaries store short years; promoted to 2019:2021
PRBL <- c(0.5, 0.025, 0.975)
LVLS <- bdl(c(0, 3.5), c(0, 3.5))

# ===========================================================================
# 1. build_discount_levels(): pairing, recycling, validation
# ===========================================================================
expect_equal(nrow(LVLS), 2L, info = "c(0, 3.5) gives two levels")
expect_equal(LVLS$label, c("0%", "3.5%"), info = "equal rates label as one rate")
expect_equal(bdl(3.5, 3.5)$label, "3.5%", info = "a scalar rate gives one level")

# Rates are PAIRED with each other (a length-1 rate recycles), never crossed:
# c(0, 1.5) x 3.5 is two levels, not four.
r <- bdl(c(0, 1.5), 3.5)
expect_equal(nrow(r), 2L, info = "a length-1 rate is recycled, not crossed")
expect_equal(r$cost, c(3.5, 3.5), info = "...against every QALY rate")
expect_equal(r$label, c("QALYs 0%/costs 3.5%", "QALYs 1.5%/costs 3.5%"),
             info = "differing rates spell both out")
expect_equal(r$qaly_label, c("0%", "1.5%"),
             info = "the QALY-only label carries just the QALY rate")

expect_equal(nrow(bdl(c(0, 0, 3.5), c(0, 0, 3.5))), 2L,
             info = "duplicate levels are de-duplicated")
expect_error(bdl(c(0, 2, 4), c(0, 2)), info = "incompatible lengths are an error")
expect_error(bdl(-100, 0), info = "a rate of -100% or below is an error")
expect_error(bdl(NA_real_, 0), info = "a non-finite rate is an error")

# The contract is EQUAL lengths or one length-1 rate -- NOT "any divisor length".
# The guard used to be `n %% length(q) != 0L`, which accepted 4-vs-2 and 6-vs-3
# and let rep_len() invent pairings nobody asked for: c(0,1.5,3.5,5) against
# c(0,3.5) recycled the costs to c(0,3.5,0,3.5) and published
# "QALYs 3.5%/costs 0%" and "QALYs 5%/costs 3.5%" -- discounting health but not
# money. It also made rejection non-monotonic: 4-vs-2 passed while 3-vs-2 failed.
expect_error(bdl(c(0, 1.5, 3.5, 5), c(0, 3.5)),
             info = "a divisor length is NOT recycling: 4-vs-2 is an error")
expect_error(bdl(c(0, 1.5, 3, 3.5, 5, 7), c(0, 3.5, 5)),
             info = "6-vs-3 is an error too")
expect_error(bdl(c(0, 3.5), c(0, 1.5, 3.5, 5)),
             info = "the same holds with the longer vector on the cost side")
expect_silent(bdl(c(0, 1.5, 3.5, 5), 3.5))
expect_silent(bdl(3.5, c(0, 1.5, 3.5, 5)))
expect_equal(nrow(bdl(c(0, 1.5, 3.5, 5), 3.5)), 4L,
             info = "a genuine length-1 rate still recycles against any length")

# ===========================================================================
# 2. discount_factor() / expand_discount_levels()
# ===========================================================================
expect_equal(dfac(2023L, 3.5, BASE), 1 / 1.035^4,
             info = "compounds from the base year")
expect_equal(dfac(2017L, 3.5, BASE), 1,
             info = "years before the base year are undiscounted")
expect_true(all(dfac(2019:2050, 0, BASE) == 1),
            info = "a 0% rate leaves every year unchanged")

# A QALYs-only table depends on the QALY rate alone, so levels sharing that
# rate collapse; the same levels stay distinct in a costs table.
shared <- bdl(c(3.5, 3.5), c(0, 5))
one <- data.table(year = 2019:2021, v = 1)
expect_equal(uniqueN(expand_lv(one, "v", shared, BASE, "qaly")$discount), 1L,
             info = "levels sharing a QALY rate collapse in a QALYs table")
expect_equal(uniqueN(expand_lv(one, "v", shared, BASE, "cost")$discount), 2L,
             info = "...and stay distinct in a costs table")
expect_error(expand_lv(data.table(v = 1), "v", LVLS, BASE, "qaly"),
             info = "discounting without a `year` column is a clear error")

# ===========================================================================
# 3. tbl_smmrs_core(): the four money/QALY branches
# ===========================================================================
# Flat 100 QALYs and 100 of each cost per year per scenario, so the discounted
# cumulative total is a geometric sum we can write down exactly.
make_qalys <- function() {
  d <- CJ(mc = 1:2, scenario = c("sc0", "sc1"), year = YEARS, sorted = FALSE)
  d[, EQ5D5L := fifelse(scenario == "sc1", 110, 100)]
  d[, popsize := 1000][]
}
make_costs <- function() {
  d <- CJ(mc = 1:2, scenario = c("sc0", "sc1"), year = YEARS, sorted = FALSE)
  d[, healthcare_cost := fifelse(scenario == "sc1", 90, 100)]
  d[, total_cost := healthcare_cost]
  d[, popsize := 1000][]
}

run_core <- function(tt, what, levels = LVLS) {
  dir <- tempfile("disc_"); dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE, diseases = NULL)))
  env$private <- list()
  f <- Simulation$private_methods$tbl_smmrs_core
  environment(f) <- env
  tt <- copy(tt)
  tt[, year := year + 2000L]   # export_main_tables() promotes before dispatch
  f(tt = tt, what = what, population = "ons", strata = list("year"),
    prbl = PRBL, baseline_year = BASE, comparator_scenario = "sc0",
    comparison_starting_year = BASE, tables_dir = dir,
    discount_levels = levels, discount_from_year = BASE)
  dir
}
rd1 <- function(dir) fread(file = file.path(dir, list.files(dir, pattern = "\\.csv$")[1]))

# Quantile columns are named <metric prefix><percent> (qalys_50%, costs_50%,
# value_50% ...), so resolve the median column by pattern instead of hardcoding
# one metric's name.
med_col <- function(t) {
  cn <- grep("50%$", names(t), value = TRUE)
  stopifnot(length(cn) == 1L)
  cn
}
# Median series of `t`, in year order.
pick <- function(t) t[order(year)][[med_col(t)]]

# ---- qalys ----------------------------------------------------------------
q <- rd1(run_core(make_qalys(), "qalys"))
expect_true("discount" %in% names(q), info = "qalys table gains a `discount` column")
expect_equal(sort(unique(q$discount)), c("0%", "3.5%"),
             info = "qalys table carries both levels")

qc <- q[type == "QALYs_cuml" & scenario == "sc0"]
# 0% cumulative over 2019..2021 = 100, 200, 300.
expect_equal(pick(qc[discount == "0%"]), c(100, 200, 300),
             info = "the 0% level reproduces the undiscounted cumulative QALYs")
# 3.5%: 2019 undiscounted, then 1/1.035, 1/1.035^2 - cumulated.
expect_equal(pick(qc[discount == "3.5%"]), cumsum(100 / 1.035^(0:2)),
             info = "the 3.5% level discounts each year's flow before cumulating")

# ---- costs ----------------------------------------------------------------
cst <- rd1(run_core(make_costs(), "costs"))
expect_equal(sort(unique(cst$discount)), c("0%", "3.5%"),
             info = "costs table carries both levels")
cc <- cst[type == "costs_cuml" & scenario == "sc0" & costs_type == "total_cost"]
expect_equal(pick(cc[discount == "0%"]), c(100, 200, 300),
             info = "the 0% level reproduces the undiscounted cumulative costs")
# The costs branch quantiles with rounding = TRUE, so compare to the rounded
# closed form rather than the raw one.
expect_equal(pick(cc[discount == "3.5%"]), round(cumsum(100 / 1.035^(0:2))),
             info = "the 3.5% level discounts costs before cumulating")

# ---- net_qalys / net_costs ------------------------------------------------
nq <- rd1(run_core(make_qalys(), "net_qalys"))
expect_equal(sort(unique(nq$discount)), c("0%", "3.5%"),
             info = "net_qalys table carries both levels")
nqc <- nq[type == "net_QALYs_cuml"]
expect_equal(pick(nqc[discount == "0%"]), c(10, 20, 30),
             info = "undiscounted net QALYs are the raw sc1 - sc0 difference")
# Discounting is linear, so differencing then discounting == discounting then
# differencing: the net series is the gross one scaled by the same factors.
expect_equal(pick(nqc[discount == "3.5%"]), cumsum(10 / 1.035^(0:2)),
             info = "net QALYs discount identically to the gross series")

nc <- rd1(run_core(make_costs(), "net_costs"))
expect_equal(sort(unique(nc$discount)), c("0%", "3.5%"),
             info = "net_costs table carries both levels")
ncc <- nc[type == "net_costs_cuml" & costs_type == "total_cost"]
expect_equal(pick(ncc[discount == "0%"]), c(-10, -20, -30),
             info = "undiscounted net costs are the raw saving")

# ---- the discount level must NOT reach the file name ----------------------
d_q <- run_core(make_qalys(), "qalys")
expect_equal(list.files(d_q, pattern = "\\.csv$"),
             "QALYs by year (not standardised).csv",
             info = "levels live in a column, so the file set is unchanged")

# ---- a single level reproduces the classic one-block table ----------------
q1 <- rd1(run_core(make_qalys(), "qalys", levels = bdl(0, 0)))
expect_equal(unique(q1$discount), "0%",
             info = "one level -> one block of rows")
expect_equal(pick(q1[type == "QALYs_cuml" & scenario == "sc0"]), c(100, 200, 300),
             info = "...matching the undiscounted figures exactly")

# ===========================================================================
# 4. CEA: discount levels x wtp are CROSSED
# ===========================================================================
WTP <- c(20000, 30000)
NMB <- paste0("NMB_at_wtp_", format(WTP, scientific = FALSE, trim = TRUE))

run_cea <- function(levels = LVLS) {
  dir <- tempfile("ceadisc_"); dir.create(dir, recursive = TRUE)
  src <- list(qalys = make_qalys(), costs = make_costs())
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE, diseases = NULL)))
  env$private <- list(read_summary_dataset = function(summary_type, standardization) {
    if (is.null(src[[summary_type]])) NULL else copy(src[[summary_type]])
  })
  f <- Simulation$private_methods$export_cea_tables
  environment(f) <- env
  f(prbl = PRBL, summaries_dir = dir, tables_dir = dir,
    comparator_scenario = "sc0", baseline_year = BASE, wtp = WTP,
    discount_levels = levels, discount_from_year = NULL,
    custom_costs_in_healthcare = NULL, strata = list("year"))
  dir
}

cea <- fread(file = file.path(run_cea(),
  "cost-effectiveness by year (healthcare-EQ5D5L) (not standardised).csv"))
cea[, type := as.character(type)]

expect_equal(sort(unique(cea$discount)), c("0%", "3.5%"),
             info = "CEA tables carry both discount levels")

# Every threshold at every level: crossed, not paired.
got <- unique(cea[type %chin% NMB, .(discount, type)])
expect_equal(nrow(got), length(NMB) * 2L,
             info = "all (discount, wtp) combinations are written")
expect_equal(uniqueN(cea[type %chin% NMB, .N, by = .(discount, type)]$N), 1L,
             info = "every (discount, wtp) cell spans the same rows")

# ICER is wtp-independent: one row per (discount, scenario, year).
icer <- cea[type == "ICER"]
expect_equal(uniqueN(icer[, .(discount, scenario, year)]), nrow(icer),
             info = "ICER appears once per discount level, not once per wtp")

# NMB rises with wtp within a level. Quantile-safe: wtp * Q - C is increasing
# in wtp for every draw when Q > 0, hence for every quantile of it.
VCOL <- med_col(cea)
wide <- dcast(cea, discount + year ~ type, value.var = VCOL)
expect_true(all(wide[[NMB[1]]] <= wide[[NMB[2]]] + 1e-6),
            info = "NMB increases with wtp within each discount level")

# Discounting reaches the NMB, not just the raw columns.
tw <- dcast(cea[type == NMB[2]], year ~ discount, value.var = VCOL)
expect_false(isTRUE(all.equal(tw[year > BASE][["0%"]], tw[year > BASE][["3.5%"]])),
             info = "NMB at a given wtp differs between discount levels")


# ---- Two-digit year shorthand -------------------------------------------
# The summaries store short years, so callers naturally pass `19` for 2019.
# `baseline_year_for_change_outputs` had always been promoted; when multi-level
# discounting landed, `discount_from_year` was not. Both now go through
# promote_short_year(), so the promotion cannot be added to one and forgotten
# for the other.
psy <- IMPACTncdEngland:::promote_short_year

expect_equal(psy(19), 2019, info = "two-digit year is promoted by 2000")
expect_equal(psy(BASE), BASE, info = "a full year is left alone")
expect_equal(psy(100), 2100, info = "100 is shorthand (boundary is inclusive)")
expect_equal(psy(101), 101, info = "101 is not shorthand (boundary is exclusive)")
expect_equal(psy(0), 2000, info = "0 is shorthand for 2000")
expect_true(is.null(psy(NULL)), info = "NULL passes through for the is.null branch")
expect_equal(psy(NA_integer_), NA_integer_, info = "NA is not silently promoted")
expect_equal(psy(c(19, 20)), c(19, 20), info = "only a scalar is promoted")

# Why it matters: an unpromoted 19 is read as the year 19 AD, so the exponent
# becomes ~2000 and the discount factor underflows. Nothing errors -- the tables
# are simply published as zeros, which is why this is a test and not a comment.
expect_true(dfac(2043L, 3.5, 19) < 1e-30,
            info = "unpromoted shorthand collapses the discount factor to ~0")
expect_equal(dfac(2043L, 3.5, psy(19)), 1 / 1.035^24,
             info = "promoted shorthand gives the correct present value")
expect_equal(dfac(2043L, 3.5, psy(19)), dfac(2043L, 3.5, BASE),
             info = "19 and 2019 are equivalent once promoted")
