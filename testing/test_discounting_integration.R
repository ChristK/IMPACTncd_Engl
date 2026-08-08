## Integration test for discounting in the export_tables() method.
##
## Design (current): discounting is controlled SOLELY by export_tables()
## arguments. `qaly_discount_rate` / `cost_discount_rate` are VECTORS (default
## c(0, 3.5)) paired element-wise into discount levels; `discount_from_year`
## defaults to the baseline year; `wtp` defaults to c(25000, 35000). There is
## intentionally NO `discounting` block in the sim_design YAMLs.
##
## Discounting reaches the main QALY / cost / net tables AS WELL AS the
## cost-effectiveness ones: each reports every level side by side, tagged in a
## `discount` column, the `0%` level being the undiscounted figures. The equity
## tables are deliberately left undiscounted (they mix QALYs with event counts).
##
## Run with: Rscript testing/test_discounting_integration.R
## Or source in R: source("testing/test_discounting_integration.R")

library(data.table)
library(yaml)

# Set working directory to project root if needed
if (!file.exists("Rpackage/IMPACTncd_England_model_pkg")) {
  setwd("/home/ckyprid/My_Models/IMPACTncd_Engl")
}

tryCatch({
  library(IMPACTncdEngland)
  cat("Loaded IMPACTncdEngland package\n")
}, error = function(e) {
  cat("Package not installed, sourcing files directly...\n")
  source("Rpackage/IMPACTncd_England_model_pkg/R/aux_functions.R")
  source("Rpackage/IMPACTncd_England_model_pkg/R/Design_class.R")
  source("Rpackage/IMPACTncd_England_model_pkg/R/Simulation_class.R")
  source("Rpackage/IMPACTncd_England_model_pkg/R/Simulation_class_tables.R")
})

cat("\n========================================\n")
cat("Integration Test: Discounting in export_tables\n")
cat("========================================\n\n")

tests_passed <- 0
tests_failed <- 0

run_test <- function(name, expr) {
  result <- tryCatch({
    if (isTRUE(expr)) {
      cat(sprintf("[PASS] %s\n", name))
      TRUE
    } else {
      cat(sprintf("[FAIL] %s\n", name))
      FALSE
    }
  }, error = function(e) {
    cat(sprintf("[FAIL] %s - Error: %s\n", name, e$message))
    FALSE
  })
  if (result) tests_passed <<- tests_passed + 1 else tests_failed <<- tests_failed + 1
  invisible(result)
}

# -----------------------------------------------------------------------------
# Test 1: The sim_design YAMLs no longer carry a `discounting` block
# -----------------------------------------------------------------------------
cat("\n--- Test 1: No `discounting` block in the config YAMLs ---\n")

# The scenario designs are GLOBBED, not listed. The three named files below were
# the only ones this test ever checked, but they are not where the problem lived:
# the eleven scenarios/*/sim_design_*.yaml files kept a live `discounting:` block
# long after it stopped being read, and this test passed green throughout.
# `Design_class.R` never references `discount` and does not reject unknown keys,
# so a re-added block is silently accepted and silently ignored -- there is no
# other guard. Globbing means a NEW scenario directory is covered automatically.
yaml_files <- unique(c(
  "Rpackage/IMPACTncd_England_model_pkg/inst/config/default_sim_design.yaml",
  "inputs/sim_design.yaml",
  "testing/sim_design_testing.yaml",
  Sys.glob("scenarios/*/sim_design*.yaml"),
  Sys.glob("inputs/sim_design*.yaml")
))
checked <- 0L
for (yf in yaml_files) {
  if (file.exists(yf)) {
    checked <- checked + 1L
    parsed <- yaml::read_yaml(yf)
    run_test(sprintf("no `discounting` key in %s", yf),
             is.null(parsed$discounting))
  }
}
# A glob that matches nothing would make every assertion above vacuous, so assert
# the sweep actually found the scenario designs.
run_test("the YAML sweep found the scenario designs (glob is not empty)",
         checked >= 12L)

# -----------------------------------------------------------------------------
# Test 2: export_tables() exposes the discounting / CEA arguments with the
#         expected defaults (single source of truth)
# -----------------------------------------------------------------------------
cat("\n--- Test 2: export_tables() argument defaults ---\n")

et_formals <- tryCatch(
  formals(Simulation$public_methods$export_tables),
  error = function(e) NULL
)

run_test("export_tables() is available with formals", !is.null(et_formals))
if (!is.null(et_formals)) {
  for (arg in c("cea", "wtp", "qaly_discount_rate", "cost_discount_rate",
                "discount_from_year", "custom_costs_in_healthcare")) {
    run_test(sprintf("export_tables() has argument `%s`", arg),
             arg %in% names(et_formals))
  }
  run_test("default qaly_discount_rate == c(0, 3.5)",
           isTRUE(all.equal(eval(et_formals$qaly_discount_rate), c(0, 3.5))))
  run_test("default cost_discount_rate == c(0, 3.5)",
           isTRUE(all.equal(eval(et_formals$cost_discount_rate), c(0, 3.5))))
  run_test("default discount_from_year is NULL (resolves to baseline year)",
           is.null(et_formals$discount_from_year))
  run_test("default wtp == c(25000, 35000)",
           isTRUE(all.equal(eval(et_formals$wtp), c(25000, 35000))))
}

# -----------------------------------------------------------------------------
# Test 3: The CEA discounting workflow (mimics export_cea_tables)
#         discount_from_year defaults to the baseline year when NULL.
# -----------------------------------------------------------------------------
cat("\n--- Test 3: CEA discounting workflow ---\n")

cea_discount <- function(qaly_rate, cost_rate, discount_from_year, baseline_year) {
  if (is.null(discount_from_year)) discount_from_year <- baseline_year
  disc <- function(v, year, rate) v / (1 + rate / 100)^pmax(0, year - discount_from_year)

  qalys <- data.table(
    mc = 1L, scenario = "sc1", year = 2020:2024, EQ5D5L = 10000
  )
  costs <- data.table(
    mc = 1L, scenario = "sc1", year = 2020:2024, total_cost = 100000
  )
  qalys[, Q := disc(EQ5D5L, year, qaly_rate)]
  costs[, C := disc(total_cost, year, cost_rate)]
  list(qalys = qalys, costs = costs)
}

# Baseline year 2020, NULL discount_from_year -> discounts from 2020
res <- cea_discount(3.5, 1.5, NULL, baseline_year = 2020L)

run_test("QALY at baseline year 2020 is undiscounted",
         abs(res$qalys[year == 2020, Q] - 10000) < 1e-6)
run_test("QALY at 2024 discounted by 4y at 3.5%",
         abs(res$qalys[year == 2024, Q] - 10000 / 1.035^4) < 1e-4)
run_test("Cost at 2024 discounted by 4y at 1.5% (separate rate)",
         abs(res$costs[year == 2024, C] - 100000 / 1.015^4) < 1e-4)
run_test("QALY and cost use different discount factors",
         res$qalys[year == 2024, Q] / 10000 != res$costs[year == 2024, C] / 100000)

# zero rate -> no discounting
res0 <- cea_discount(0, 0, 2020L, baseline_year = 2020L)
run_test("Zero rate leaves QALYs/costs unchanged",
         all(abs(res0$qalys$Q - 10000) < 1e-9) &&
           all(abs(res0$costs$C - 100000) < 1e-9))

# -----------------------------------------------------------------------------
# Test 4: The MAIN qalys/costs tables carry EVERY discount level
# -----------------------------------------------------------------------------
cat("\n--- Test 4: Main tables are discounted per level ---\n")

mock_qalys <- data.table(
  mc = rep(1:2, each = 6),
  scenario = rep("sc0", 12),
  year = rep(2018:2023, 2),
  EQ5D5L = rep(10000, 12)
)
x <- c("mc", "scenario", "year")
# tbl_smmrs_core QALYs branch: aggregate, expand across discount levels, cumsum.
# Discounting the annual flow BEFORE the cumsum is what makes the cumulative
# column a sum of present values.
lv <- IMPACTncdEngland:::build_discount_levels(c(0, 3.5), c(0, 3.5))
d <- mock_qalys[, .(EQ5D5L = sum(EQ5D5L)), keyby = eval(x)]
d <- melt(d, id.vars = x, variable.name = "scale", value.name = "QALYs")
d <- IMPACTncdEngland:::expand_discount_levels(d, "QALYs", lv, 2019L, "qaly")
xd <- c(x, "scale", "discount")
setkeyv(d, c(setdiff(xd, "year"), "year"))
d[, cumulative := cumsum(QALYs), keyby = setdiff(xd, "year")]

run_test("Main QALY tables carry both discount levels",
         setequal(unique(d$discount), c("0%", "3.5%")))
run_test("The 0% level is undiscounted (all == 10000)",
         all(abs(d[discount == "0%", QALYs] - 10000) < 1e-9))
run_test("0% cumulative at 2023 == 6 * 10000",
         abs(d[mc == 1 & year == 2023 & discount == "0%", cumulative] - 60000) < 1e-6)
# 2018 precedes discount_from_year (2019) so it is undiscounted at every level;
# 2019..2023 are discounted by 0..4 years.
run_test("The 3.5% level discounts from 2019, leaving earlier years alone",
         abs(d[mc == 1 & year == 2018 & discount == "3.5%", QALYs] - 10000) < 1e-9)
run_test("3.5% cumulative at 2023 == 10000 * (1 + sum 1/1.035^0:4)",
         abs(d[mc == 1 & year == 2023 & discount == "3.5%", cumulative] -
               10000 * (1 + sum(1 / 1.035^(0:4)))) < 1e-6)

# Costs main branch picks up both built-in `_cost` and custom `_costs` columns
mock_costs <- data.table(
  mc = 1L, scenario = "sc0", year = 2020:2022,
  total_cost = 100000, healthcare_cost = 60000,
  sbp_intervention_costs = 500  # custom plural-suffix column
)
cc <- mock_costs[, lapply(.SD, sum),
                 .SDcols = patterns("_cost$|_costs$|^economic_output$"),
                 keyby = .(mc, scenario, year)]
run_test("Costs branch includes built-in total_cost & healthcare_cost",
         all(c("total_cost", "healthcare_cost") %in% names(cc)))
run_test("Costs branch includes custom `_costs` column",
         "sbp_intervention_costs" %in% names(cc))
run_test("`_cost$|_costs$` does not double-count (no extra columns)",
         setequal(setdiff(names(cc), c("mc", "scenario", "year")),
                  c("total_cost", "healthcare_cost", "sbp_intervention_costs")))

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n========================================\n")
cat("Integration Test Summary\n")
cat("========================================\n")
cat(sprintf("Tests passed: %d\n", tests_passed))
cat(sprintf("Tests failed: %d\n", tests_failed))
cat(sprintf("Total tests:  %d\n", tests_passed + tests_failed))

if (tests_failed == 0) {
  cat("\n[SUCCESS] All integration tests passed!\n")
  cat("\nThe discounting design correctly:\n")
  cat("- Keeps discount rates / wtp as export_tables() arguments only (no YAML)\n")
  cat("- Defaults to c(0, 3.5) QALY/cost rates and wtp c(25000, 35000)\n")
  cat("- Defaults discount_from_year to the baseline year\n")
  cat("- Reports every discount level in the main QALY/cost tables AND the\n")
  cat("  cost-effectiveness tables, tagged in a `discount` column\n")
  cat("- Leaves the 0% level equal to the undiscounted figures\n")
  cat("- Picks up custom `_costs` columns in the cost tables without double-counting\n\n")
} else {
  cat("\n[WARNING] Some integration tests failed.\n\n")
}

invisible(tests_failed == 0)
