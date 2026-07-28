# Tests for add_qimd_from_dimd(): deriving IMD quintiles from IMD deciles so
# ANY table family can be stratified by `qimd` without it having been named in
# `strata_for_output`.
#
# The correctness claim is that the derivation is EXACT, not approximate: the
# helper only relabels, and the caller's existing group-by does the collapsing,
# so summing two deciles reproduces exactly what aggregating by quintile at
# summary time would have produced. Rates follow because they are computed as
# summed-numerator / summed-denominator AFTER the group-by. The end-to-end test
# below asserts precisely that through the real export_main_tables() code path.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

add_qimd <- IMPACTncdEngland:::add_qimd_from_dimd
dimd_lv <- c("1 most deprived", as.character(2:9), "10 least deprived")
qimd_lv <- c("1 most deprived", as.character(2:4), "5 least deprived")

# ===========================================================================
# 1. The mapping itself
# ===========================================================================
d <- data.table(dimd = factor(dimd_lv, levels = dimd_lv))
expect_true(add_qimd(d))
expect_equal(levels(d$qimd), qimd_lv)
expect_equal(as.character(d$qimd),
             qimd_lv[c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)],
             info = "deciles 1,2 -> quintile 1, ... 9,10 -> quintile 5")
expect_equal(as.integer(d$qimd), c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L),
             info = "rank order preserved: 1 = most deprived at both scales")

# Works from a character dimd column too.
d2 <- data.table(dimd = dimd_lv)
expect_true(add_qimd(d2))
expect_equal(as.character(d2$qimd), as.character(d$qimd))

# ===========================================================================
# 2. No-ops
# ===========================================================================
d3 <- data.table(dimd = dimd_lv, qimd = "preexisting")
expect_false(add_qimd(d3))
expect_equal(unique(d3$qimd), "preexisting",
             info = "a native qimd column is never overwritten")
d4 <- data.table(sex = c("men", "women"))
expect_false(add_qimd(d4), info = "no dimd to derive from")
expect_false("qimd" %in% names(d4))
expect_false(add_qimd(NULL), info = "NULL summary is tolerated")

# Unrecognised deprivation labels must stop, not silently mis-map.
expect_error(add_qimd(data.table(dimd = c("1 most deprived", "banana"))),
             pattern = "unexpected `dimd` level")

# ===========================================================================
# 3. End-to-end: main tables by qimd, derived == natively aggregated
# ===========================================================================
# A prevalence summary at decile granularity. Counts and popsize vary by
# decile so the quintile collapse is a genuine re-aggregation, and the
# prevalence RATE is a non-linear function of the two (sum/sum), which is what
# makes this a real test rather than a restatement of the mapping.
make_prvl <- function() {
  d <- CJ(mc = 1:3, scenario = "sc0", year = 19:20,
          dimd = factor(dimd_lv, levels = dimd_lv),
          agegrp = c("50-54", "55-59"), sex = c("men", "women"),
          sorted = FALSE)
  k <- 11L - as.integer(d$dimd)
  d[, chd_prvl := as.numeric(k * 100L + mc)]
  d[, popsize := as.numeric(k * 1000L + 500L)]   # deliberately uneven
  d[]
}

run_main <- function(src) {
  dir <- tempfile("qimdtest_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  core <- Simulation$private_methods$tbl_smmrs_core
  environment(core) <- env
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      if (summary_type == "prvl" && standardization == "scaled_up") {
        copy(src)
      } else {
        NULL
      }
    },
    tbl_smmrs_core = core
  )
  f <- Simulation$private_methods$export_main_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
    baseline_year = 2019L, output_dir = dir, tables_dir = dir,
    comparator_scenario = "sc0", two_agegrps = FALSE,
    strata_ons = list(c("year", "qimd")), strata_esp = list())
  dir
}

# (a) summaries carry dimd only -> qimd is derived
derived_dir <- run_main(make_prvl())
fn <- "prevalence by year-qimd (not standardised).csv"
expect_true(file.exists(file.path(derived_dir, fn)),
            info = "main tables can be stratified by qimd with dimd-only summaries")
derived <- fread(file.path(derived_dir, fn))
expect_equal(sort(unique(derived$qimd)), sort(qimd_lv),
             info = "all five quintiles present")

# (b) summaries carry a native qimd, pre-aggregated at summary time
native <- make_prvl()
native[, qimd := factor(qimd_lv[c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L)][
  as.integer(dimd)], levels = qimd_lv)]
native[, dimd := NULL]
native <- native[, .(chd_prvl = sum(chd_prvl), popsize = sum(popsize)),
                 by = .(mc, scenario, year, qimd, agegrp, sex)]
native_dir <- run_main(native)
expect_equal(fread(file.path(native_dir, fn)), derived,
             info = "derived qimd == qimd aggregated at summary time (rates too)")

# (c) the quintile table is a real re-aggregation, not a relabelled decile table
run_main_dimd <- function() {
  dir <- tempfile("qimdtest_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  core <- Simulation$private_methods$tbl_smmrs_core
  environment(core) <- env
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      if (summary_type == "prvl" && standardization == "scaled_up") {
        make_prvl()
      } else {
        NULL
      }
    },
    tbl_smmrs_core = core
  )
  f <- Simulation$private_methods$export_main_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
    baseline_year = 2019L, output_dir = dir, tables_dir = dir,
    comparator_scenario = "sc0", two_agegrps = FALSE,
    strata_ons = list(c("year", "dimd")), strata_esp = list())
  dir
}
decile <- fread(file.path(run_main_dimd(),
                          "prevalence by year-dimd (not standardised).csv"))
expect_equal(nrow(decile), 2L * nrow(derived),
             info = "10 deciles vs 5 quintiles, same years and diseases")
# The quintile rate is the popsize-weighted blend of its two decile rates, so
# it must differ from either of them (the fixture's rates vary by decile).
q1 <- derived[qimd == "1 most deprived" & year == 2019L, `prvl_rate_50.0%`]
d12 <- decile[dimd %in% c("1 most deprived", "2") & year == 2019L,
              `prvl_rate_50.0%`]
expect_true(length(q1) == 1L && length(d12) == 2L)
expect_true(all(abs(q1 - d12) > 1e-9),
            info = "quintile rate is a genuine re-aggregation, not a copy")
expect_true(q1 > min(d12) && q1 < max(d12),
            info = "and it lies between the two decile rates it blends")
