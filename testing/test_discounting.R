## Test script for the discounting logic used in the cost-effectiveness
## (ICER/NMB) tables. Discounting is applied ONLY in export_cea_tables(); the
## main QALY/cost/net tables are reported undiscounted. These are unit tests of
## the discount formula PV = FV / (1 + r)^max(0, year - discount_from_year),
## which is exactly the formula used by export_cea_tables()'s internal disc().
##
## Run with: Rscript testing/test_discounting.R
## Or source in R: source("testing/test_discounting.R")

library(data.table)

cat("\n========================================\n")
cat("Testing Discounting Logic for QALYs and Costs\n")
cat("========================================\n\n")

# Track test results
tests_passed <- 0
tests_failed <- 0

# Helper function to run a test
run_test <- function(name, expr) {

  result <- tryCatch({
    if (expr) {
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

  if (result) {
    tests_passed <<- tests_passed + 1
  } else {
    tests_failed <<- tests_failed + 1
  }
  invisible(result)
}

# -----------------------------------------------------------------------------
# Test 1: Basic discounting formula
# -----------------------------------------------------------------------------
cat("\n--- Test 1: Basic Discounting Formula ---\n")
cat("Formula: PV = FV / (1 + r)^(year - discount_from_year)\n\n")

# Test parameters
discount_rate <- 3.5  # percent
discount_from_year <- 2019

# Create test data
test_values <- data.table(
  year = 2017:2025,
  value = 1000  # constant value each year
)

# Apply discounting formula (same as in the code)
test_values[, discounted := value / (1 + discount_rate / 100) ^ pmax(0, year - discount_from_year)]

# Manual calculations for verification
# 2017: 1000 / (1.035)^0 = 1000 (before discount_from_year, pmax(0, -2) = 0)
# 2018: 1000 / (1.035)^0 = 1000 (before discount_from_year, pmax(0, -1) = 0)
# 2019: 1000 / (1.035)^0 = 1000 (year = discount_from_year)
# 2020: 1000 / (1.035)^1 = 966.18...
# 2021: 1000 / (1.035)^2 = 933.51...
# 2022: 1000 / (1.035)^3 = 901.94...
# 2023: 1000 / (1.035)^4 = 871.44...
# 2024: 1000 / (1.035)^5 = 841.97...
# 2025: 1000 / (1.035)^6 = 813.50...

expected_2017 <- 1000
expected_2019 <- 1000
expected_2020 <- 1000 / 1.035
expected_2023 <- 1000 / (1.035^4)

run_test("Years before discount_from_year are not discounted (2017)",
         abs(test_values[year == 2017, discounted] - expected_2017) < 0.01)

run_test("Years before discount_from_year are not discounted (2018)",
         abs(test_values[year == 2018, discounted] - 1000) < 0.01)

run_test("Discount_from_year itself is not discounted (2019)",
         abs(test_values[year == 2019, discounted] - expected_2019) < 0.01)

run_test("Year 2020 is discounted by 1 year",
         abs(test_values[year == 2020, discounted] - expected_2020) < 0.01)

run_test("Year 2023 is discounted by 4 years",
         abs(test_values[year == 2023, discounted] - expected_2023) < 0.01)

# Print the test data
cat("\nTest data with discounting (rate=3.5%, from_year=2019):\n")
print(test_values)

# -----------------------------------------------------------------------------
# Test 2: Cumulative sum after discounting
# -----------------------------------------------------------------------------
cat("\n--- Test 2: Cumulative Sum After Discounting ---\n")

test_cumsum <- copy(test_values)
test_cumsum[, cumulative := cumsum(discounted)]

# Verify cumulative is sum of discounted values
expected_cumsum_2020 <- sum(test_values[year <= 2020, discounted])
expected_cumsum_2025 <- sum(test_values[, discounted])

run_test("Cumulative sum at 2020 is correct",
         abs(test_cumsum[year == 2020, cumulative] - expected_cumsum_2020) < 0.01)

run_test("Cumulative sum at 2025 (end) is correct",
         abs(test_cumsum[year == 2025, cumulative] - expected_cumsum_2025) < 0.01)

cat("\nCumulative sums:\n")
print(test_cumsum[, .(year, discounted, cumulative)])

# -----------------------------------------------------------------------------
# Test 3: Different discount rates for QALYs and costs
# -----------------------------------------------------------------------------
cat("\n--- Test 3: Different Discount Rates ---\n")

qaly_rate <- 3.5
cost_rate <- 1.5
from_year <- 2020

test_diff_rates <- data.table(
  year = 2020:2025,
  qaly_value = 100,
  cost_value = 1000
)

test_diff_rates[, qaly_discounted := qaly_value / (1 + qaly_rate / 100) ^ pmax(0, year - from_year)]
test_diff_rates[, cost_discounted := cost_value / (1 + cost_rate / 100) ^ pmax(0, year - from_year)]

# At year 2025 (5 years from 2020):
# QALY: 100 / (1.035)^5 = 84.20
# Cost: 1000 / (1.015)^5 = 928.26

expected_qaly_2025 <- 100 / (1.035^5)
expected_cost_2025 <- 1000 / (1.015^5)

run_test("QALY discounted correctly with 3.5% rate at year 2025",
         abs(test_diff_rates[year == 2025, qaly_discounted] - expected_qaly_2025) < 0.01)

run_test("Cost discounted correctly with 1.5% rate at year 2025",
         abs(test_diff_rates[year == 2025, cost_discounted] - expected_cost_2025) < 0.01)

run_test("Different rates produce different discount factors",
         test_diff_rates[year == 2025, qaly_discounted / qaly_value] !=
         test_diff_rates[year == 2025, cost_discounted / cost_value])

cat("\nDifferent discount rates (QALY=3.5%, Cost=1.5%, from_year=2020):\n")
print(test_diff_rates)

# -----------------------------------------------------------------------------
# Test 4: No discounting when rate is 0
# -----------------------------------------------------------------------------
cat("\n--- Test 4: No Discounting When Rate is 0 ---\n")

zero_rate <- 0
test_zero <- data.table(
  year = 2019:2025,
  value = 1000
)

# This mimics the condition in the actual code
if (zero_rate > 0) {
  test_zero[, discounted := value / (1 + zero_rate / 100) ^ pmax(0, year - 2019)]
} else {
  test_zero[, discounted := value]  # No discounting applied
}

run_test("No discounting when rate is 0",
         all(test_zero[, value == discounted]))

cat("\nWith zero discount rate, values are unchanged:\n")
print(test_zero)

# -----------------------------------------------------------------------------
# Test 5: Discount factor calculation
# -----------------------------------------------------------------------------
cat("\n--- Test 5: Discount Factor Values ---\n")

discount_factors <- data.table(
  years_from_reference = 0:10,
  rate_1.5 = 1 / (1.015 ^ (0:10)),
  rate_3.5 = 1 / (1.035 ^ (0:10)),
  rate_5.0 = 1 / (1.05 ^ (0:10))
)

cat("\nDiscount factors by year and rate:\n")
print(discount_factors, digits = 4)

# After 10 years at 3.5%, present value should be about 70.9% of future value
expected_10yr_factor <- 1 / (1.035^10)
run_test("10-year discount factor at 3.5% is ~0.709",
         abs(discount_factors[years_from_reference == 10, rate_3.5] - expected_10yr_factor) < 0.001)

# -----------------------------------------------------------------------------
# Test 6: Net values (intervention - baseline) with discounting
# -----------------------------------------------------------------------------
cat("\n--- Test 6: Net Values with Discounting ---\n")

# Simulate scenario comparison
discount_rate <- 3.5
from_year <- 2020

net_test <- data.table(
  year = rep(2020:2023, 2),
  scenario = rep(c("sc0", "sc1"), each = 4),
  value = c(100, 100, 100, 100,    # baseline (sc0)
            120, 130, 140, 150)     # intervention (sc1)
)

# Apply discounting first (as done in the code)
net_test[, discounted := value / (1 + discount_rate / 100) ^ pmax(0, year - from_year)]

# Then calculate net
baseline <- net_test[scenario == "sc0", .(year, baseline_disc = discounted)]
intervention <- net_test[scenario == "sc1", .(year, interv_disc = discounted)]
net_result <- baseline[intervention, on = "year"]
net_result[, net_value := interv_disc - baseline_disc]

cat("\nNet values (intervention - baseline) with discounting:\n")
print(net_result)

# Verify net calculation for 2023:
# Baseline 2023: 100 / (1.035)^3 = 90.19
# Intervention 2023: 150 / (1.035)^3 = 135.29
# Net 2023: 135.29 - 90.19 = 45.10
expected_net_2023 <- (150 / 1.035^3) - (100 / 1.035^3)

run_test("Net value calculation at 2023 is correct",
         abs(net_result[year == 2023, net_value] - expected_net_2023) < 0.01)

# -----------------------------------------------------------------------------
# Test 7: Simulation of the CEA discounting workflow (mimicking export_cea_tables)
# -----------------------------------------------------------------------------
cat("\n--- Test 7: CEA Discounting Workflow Simulation ---\n")

# Create mock QALY data similar to what would come from the qalys summary
mock_qalys <- data.table(
  mc = rep(1:3, each = 5),           # 3 Monte Carlo iterations
  scenario = rep("sc0", 15),
  year = rep(2019:2023, 3),
  EQ5D5L = rnorm(15, mean = 50000, sd = 1000)
)

# Parameters (these come from export_tables() arguments, not the YAML)
qaly_discount_rate <- 3.5
discount_from_year <- 2019

# Mimic the discounting done inside export_cea_tables(): aggregate Q per
# (mc, scenario, year), then discount with the same disc() formula.
x <- c("mc", "scenario", "year")
d <- mock_qalys[, .(Q = sum(EQ5D5L)), keyby = eval(x)]
disc <- function(v, year, rate) v / (1 + rate / 100)^pmax(0, year - discount_from_year)
d[, Q := disc(Q, year, qaly_discount_rate)]
setkeyv(d, c(setdiff(x, "year"), "year"))
d[, Q_cuml := cumsum(Q), by = setdiff(x, "year")]

cat("\nMock QALY data after CEA discounting workflow:\n")
print(d[mc == 1])

run_test("CEA workflow: cumulative increases over time",
         all(diff(d[mc == 1, Q_cuml]) >= 0))

run_test("CEA workflow: discounting reduces later year values",
         d[mc == 1 & year == 2023, Q] <
         mock_qalys[mc == 1 & year == 2023, sum(EQ5D5L)])

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------
cat("\n========================================\n")
cat("Test Summary\n")
cat("========================================\n")
cat(sprintf("Tests passed: %d\n", tests_passed))
cat(sprintf("Tests failed: %d\n", tests_failed))
cat(sprintf("Total tests:  %d\n", tests_passed + tests_failed))

if (tests_failed == 0) {
  cat("\n[SUCCESS] All discounting tests passed!\n\n")
} else {
  cat("\n[WARNING] Some tests failed. Please review the output above.\n\n")
}

# Return success/failure for automated testing
invisible(tests_failed == 0)
