# Tests for the cost perspectives of private$export_cea_tables().
#
# Three perspectives are written, one file each per stratum and QALY scale:
#   societal              total_cost (+ all user *_costs)
#   healthcare            healthcare_cost (+ custom_costs_in_healthcare)
#   healthcare_socialcare healthcare_cost + socialcare_cost (+ the same custom)
#
# The point of the third is NICE's reference case (NHS + Personal Social
# Services). The tests below pin that it is genuinely a THIRD quantity -- not a
# relabelled copy of `healthcare` -- and that it is SKIPPED rather than silently
# degraded when `socialcare_cost` is absent, which is the failure mode that
# would otherwise publish misleading numbers under a name claiming to include
# social care.
#
# Same mock harness as test_export_equity_tables.R: the method's only couplings
# to the R6 object are `self$design$sim_prm$logs` and
# `private$read_summary_dataset()`.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

make_costs <- function(with_socialcare = TRUE) {
  d <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 19:21, sorted = FALSE)
  d[, healthcare_cost   := 1000 * mc + fifelse(scenario == "sc1", -50, 0)]
  d[, socialcare_cost   :=  400 * mc + fifelse(scenario == "sc1", -30, 0)]
  d[, informalcare_cost :=  200 * mc]
  # Productivity enters as a NEGATIVE cost, and a health intervention moves it:
  # sc1 keeps more people in work, so its economic output is larger. Without a
  # scenario term this column would cancel in the increment and the societal
  # perspective would collapse onto healthcare_socialcare.
  d[, indirect_cost := -5000 * mc + fifelse(scenario == "sc1", -500 * mc, 0)]
  d[, total_cost := healthcare_cost + socialcare_cost + informalcare_cost +
       indirect_cost]
  d[, prog_costs := fifelse(scenario == "sc1", 25 * mc, 0)]  # user-defined
  if (!with_socialcare) d[, socialcare_cost := NULL]
  d[]
}

make_qalys <- function() {
  d <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 19:21, sorted = FALSE)
  d[, EQ5D5L := 500 * mc + fifelse(scenario == "sc1", 3, 0)]
  d[]
}

run_cea <- function(sources, custom_in_hc = NULL, logs = FALSE) {
  dir <- tempfile("ceatest_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = logs, diseases = NULL)))
  env$private <- list(read_summary_dataset = function(summary_type, standardization) {
    if (is.null(sources[[summary_type]])) NULL else copy(sources[[summary_type]])
  })
  f <- Simulation$private_methods$export_cea_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975, 0.1, 0.9), summaries_dir = dir, tables_dir = dir,
    comparator_scenario = "sc0", baseline_year = 2019L,
    wtp = c(20000, 30000), qaly_discount_rate = 3.5, cost_discount_rate = 3.5,
    discount_from_year = NULL, custom_costs_in_healthcare = custom_in_hc,
    strata = list("year"))
  dir
}

persp_of <- function(dir) {
  f <- list.files(dir, pattern = "^cost-effectiveness by ")
  sort(sub(".*\\(([^-]+)-EQ5D5L\\).*", "\\1", f))
}

# ===========================================================================
# 0. The willingness-to-pay default is pinned
# ===========================================================================
# The thresholds appear in the OUTPUT COLUMN NAMES (NMB_at_wtp_25000), so the
# default is part of the published schema and must not drift unnoticed.
expect_equal(eval(formals(Simulation$public_methods$export_tables)$wtp),
             c(25000, 35000),
             info = "export_tables() default WTP thresholds")
expect_equal(eval(formals(Simulation$private_methods$export_cea_tables)$wtp),
             c(25000, 35000),
             info = "export_cea_tables() default agrees with its caller")

# ===========================================================================
# 1. All three perspectives are written
# ===========================================================================
d1 <- run_cea(list(costs = make_costs(), qalys = make_qalys()))
expect_equal(persp_of(d1),
             c("healthcare", "healthcare_socialcare", "societal"),
             info = "three perspectives, one file each")

rd <- function(dir, persp) {
  fread(file.path(dir, paste0(
    "cost-effectiveness by year (", persp, "-EQ5D5L) (not standardised).csv")))
}
hc  <- rd(d1, "healthcare")
hcs <- rd(d1, "healthcare_socialcare")
soc <- rd(d1, "societal")

expect_equal(names(hcs), names(hc),
             info = "the new perspective has the same schema as the others")

# ===========================================================================
# 2. It is a genuinely THIRD quantity
# ===========================================================================
dc <- function(t) t[type == "dCosts_cuml"][order(year), `value_50.0%`]
expect_false(isTRUE(all.equal(dc(hcs), dc(hc))),
             info = "healthcare_socialcare != healthcare")
expect_false(isTRUE(all.equal(dc(hcs), dc(soc))),
             info = "healthcare_socialcare != societal")

# sc1 saves on BOTH healthcare and social care, so the combined saving must be
# the larger (more negative incremental cost) of the two.
expect_true(all(dc(hcs) < dc(hc)),
            info = "adding social care savings increases the total saving")

# societal nets out productivity (indirect_cost is negative and dominates), so
# it must not be reconstructible from the two payer perspectives.
expect_true(all(abs(dc(soc)) > abs(dc(hcs))),
            info = "societal is dominated by the productivity term")

# ===========================================================================
# 3. Custom cost columns reach BOTH narrower perspectives
# ===========================================================================
d2 <- run_cea(list(costs = make_costs(), qalys = make_qalys()),
              custom_in_hc = "prog_costs")
hc2  <- rd(d2, "healthcare")
hcs2 <- rd(d2, "healthcare_socialcare")
expect_true(all(dc(hc2) > dc(hc)),
            info = "a programme cost raises the healthcare incremental cost")
expect_true(all(dc(hcs2) > dc(hcs)),
            info = "...and the healthcare_socialcare one by the same reasoning")
# The two must move by the SAME amount: they share custom_costs_in_healthcare.
expect_equal(dc(hc2) - dc(hc), dc(hcs2) - dc(hcs),
             info = "both narrower perspectives take the same custom columns")

# ===========================================================================
# 4. Missing socialcare_cost SKIPS the perspective, never degrades it
# ===========================================================================
d3 <- run_cea(list(costs = make_costs(with_socialcare = FALSE),
                   qalys = make_qalys()))
expect_equal(persp_of(d3), c("healthcare", "societal"),
             info = "no socialcare_cost -> the perspective is skipped")
expect_false(file.exists(file.path(d3, paste0(
  "cost-effectiveness by year (healthcare_socialcare-EQ5D5L) ",
  "(not standardised).csv"))),
  info = "and emphatically NOT written as a copy of the healthcare table")
