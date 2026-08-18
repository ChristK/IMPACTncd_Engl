# Tests for add_agegrp_from_age(): deriving 5-year age bands from single years
# of age, so a strata entry naming `agegrp` keeps working against summaries
# written by `export_summaries(single_year_of_age = TRUE)`.
#
# The bug this guards against: that flag replaces `agegrp` with `age` in every
# summary, while build_strata_config() merges the caller's `strata` over the
# defaults BY NAME -- so a caller who overrode `ons` but not `mrtl_ons` got
# `object 'agegrp' not found` from `keyby = eval(outstrata)`, raised inside a
# forked worker after gigabytes had been loaded. Section 4 reproduces exactly
# that call and asserts it now succeeds; section 5 asserts the numbers are the
# ones an `agegrp` run would have published, not merely *some* numbers.
#
# The correctness claim is the same one add_qimd_from_dimd() makes: the
# derivation is EXACT, not approximate, because the helper only relabels and
# the caller's existing group-by does the collapsing. Rates follow because they
# are formed as summed-numerator / summed-denominator AFTER the group-by.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

add_agegrp <- IMPACTncdEngland:::add_agegrp_from_age
add_qimd <- IMPACTncdEngland:::add_qimd_from_dimd
dimd_lv <- c("1 most deprived", as.character(2:9), "10 least deprived")

# ===========================================================================
# 1. The mapping itself
# ===========================================================================
d <- data.table(age = 30:99)
expect_true(add_agegrp(d))
expect_true(is.character(d$agegrp),
            info = "character, matching what the summaries carry on disk")
expect_equal(sort(unique(d$agegrp)),
             c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64",
               "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99"),
             info = "the stock 30-99 design's 14 five-year bands")
expect_equal(d[age == 30L, agegrp], "30-34")
expect_equal(d[age == 34L, agegrp], "30-34")
expect_equal(d[age == 35L, agegrp], "35-39", info = "band edges are closed")
expect_equal(d[age == 99L, agegrp], "95-99")
expect_equal(unique(d[, .N, by = agegrp]$N), 5L,
             info = "every band holds exactly five single years")

# The labels are IDENTICAL to the ones the model itself writes, because both go
# through the same CKutils::to_agegrp() call -- not merely similar-looking.
ref <- data.table(age = 30:99)
CKutils::to_agegrp(ref, 5L, 99L, "age", "agegrp", to_factor = FALSE)
expect_equal(d$agegrp, ref$agegrp,
             info = "same derivation the lifecourse used to build the column")

# Band edges follow the DATA, not a hard-coded 30: a design with ageL = 20
# must not be relabelled onto 30-based bands.
d20 <- data.table(age = 20:89)
expect_true(add_agegrp(d20))
expect_equal(min(d20$agegrp), "20-24")
expect_equal(max(d20$agegrp), "85-89")
# ... and an age range that does not start on a band boundary still starts its
# first band on one (to_agegrp() floors to a multiple of the band width).
d32 <- data.table(age = 32:49)
expect_true(add_agegrp(d32))
expect_equal(d32[age == 32L, agegrp], "30-34")

# ===========================================================================
# 2. No-ops
# ===========================================================================
d3 <- data.table(age = 30:34, agegrp = "preexisting")
expect_false(add_agegrp(d3))
expect_equal(unique(d3$agegrp), "preexisting",
             info = "a native agegrp column is never overwritten")
d4 <- data.table(sex = c("men", "women"))
expect_false(add_agegrp(d4), info = "no age to derive from")
expect_false("agegrp" %in% names(d4))
expect_false(add_agegrp(NULL), info = "NULL summary is tolerated")

# dis_characteristics carries NEITHER column: the call site is unconditional,
# so this must stay a silent no-op rather than an error.
expect_false(add_agegrp(data.table(mc = 1L, year = 19L, sex = "men", dimd = "2")))

# An empty summary still gains the column, so a strata entry naming `agegrp`
# yields an empty table rather than the missing-column error.
d5 <- data.table(age = integer(0))
expect_true(add_agegrp(d5))
expect_true("agegrp" %in% names(d5))
expect_equal(nrow(d5), 0L)

# ===========================================================================
# 3. It composes with the two-age-group collapse
# ===========================================================================
# `two_agegrps` relabels `agegrp`, so the derivation must run FIRST -- and the
# split must land between derived bands exactly as between native ones.
d6 <- data.table(age = 30:99)
add_agegrp(d6)
expect_true(IMPACTncdEngland:::collapse_two_agegrps(d6, 65L))
expect_equal(sort(unique(d6$agegrp)), c("30-64", "65-99"))
expect_equal(d6[age == 64L, agegrp], "30-64")
expect_equal(d6[age == 65L, agegrp], "65-99")
# A split inside a derived band is refused, just as for a native one.
d7 <- data.table(age = 30:99)
add_agegrp(d7)
expect_error(IMPACTncdEngland:::collapse_two_agegrps(d7, 62L),
             pattern = "falls INSIDE")

# ===========================================================================
# 4. Regression: the failing call, end to end
# ===========================================================================
# A single-year-of-age all-cause-mortality summary, with the DEFAULT mrtl_ons
# strata -- i.e. the exact combination that produced
#   Error in { : task 2 failed - "object 'agegrp' not found"
#
# Deaths and cases both vary with age WITHIN a band, so collapsing is a real
# re-aggregation and the published rate is a non-linear function of the two
# sums -- a broken collapse could not accidentally match.
acm_age <- CJ(mc = 1:3, scenario = "sc0", year = 19:20,
              age = c(50:54, 65:69), sex = c("men", "women"),
              dimd = factor(dimd_lv[c(1, 5, 10)], levels = dimd_lv),
              sorted = FALSE)
acm_age[, cases_chd := as.numeric(age * 10L + mc)]
acm_age[, deaths_chd := as.numeric(age %/% 7L + mc)]
acm_age[, cases_stroke := as.numeric(age * 3L + 40L)]
acm_age[, deaths_stroke := as.numeric(age %/% 11L + 1L)]

# The population denominator, on the same rows as the numerator: under a
# derived `agegrp` the two are JOINED on the strata, so both sides must gain
# the column -- and both sides must be banded identically for the native
# comparison to be like-for-like.
prvl_age <- unique(acm_age[, .(mc, scenario, year, age, sex, dimd)])
prvl_age[, popsize := 1000 + 7 * .I]

# The natively banded counterparts: exactly the same numbers, summed into
# 5-year bands at "summary time" instead of at table time.
band <- function(d, keys) {
  d <- copy(d)
  add_agegrp(d)
  d[, lapply(.SD, sum), .SDcols = setdiff(names(d), c(keys, "age", "agegrp")),
    by = c(setdiff(keys, "age"), "agegrp")]
}
acm_keys <- c("mc", "scenario", "year", "age", "sex", "dimd")
acm_band <- band(acm_age, acm_keys)
prvl_band <- band(prvl_age, acm_keys)

run_acm <- function(acm, prvl, strata_ons, strata_esp = list(),
                    two_agegrps = FALSE) {
  dir <- tempfile("agegrptest_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      if (summary_type == "all_cause_mrtl_by_dis") copy(acm)
      else if (summary_type == "prvl" && standardization == "scaled_up") copy(prvl)
      else NULL
    }
  )
  f <- Simulation$private_methods$export_all_cause_mrtl_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975),
    summaries_dir = dir, tables_dir = dir,
    two_agegrps = two_agegrps, two_agegrps_split_age = 65L,
    strata_ons = strata_ons, strata_esp = strata_esp)
  dir
}

# The default mrtl_ons list, verbatim from build_strata_config().
default_mrtl_ons <- list("year", c("year", "sex"), c("year", "agegrp"),
                         c("year", "agegrp", "sex"),
                         c("year", "agegrp", "sex", "dimd"))

by_age_dir <- run_acm(acm_age, prvl_age, default_mrtl_ons)
expect_true(file.exists(file.path(
  by_age_dir, "all-cause mortality given disease-year-agegroup (not standardised).csv")),
  info = "agegrp strata work against a single-year-of-age summary")
expect_true(file.exists(file.path(
  by_age_dir,
  "all-cause mortality given disease-year-agegroup-sex-dimd (not standardised).csv")),
  info = "and so do the deeper agegrp strata")
# The published filenames are unchanged: this family spells `agegrp` as
# `agegroup` in the suffix, and deriving rather than substituting `age` is what
# keeps that true.
expect_false(any(grepl("-year-age (", list.files(by_age_dir), fixed = TRUE)),
             info = "no `-year-age` files appear")

# ===========================================================================
# 5. Exactness: derived == natively banded, table for table
# ===========================================================================
native_dir <- run_acm(acm_band, prvl_band, default_mrtl_ons)
expect_equal(sort(list.files(by_age_dir)), sort(list.files(native_dir)),
             info = "same set of files")
# The loop below is the whole of the exactness check, so assert it is not empty:
# 5 strata x (disease-denominator + popdenom) = 10 files.
expect_equal(length(list.files(native_dir, pattern = "\\.csv$")), 10L,
             info = "every stratification produced its two tables")
for (fn in list.files(native_dir, pattern = "\\.csv$")) {
  expect_equal(fread(file = file.path(by_age_dir, fn)),
               fread(file = file.path(native_dir, fn)),
               info = paste0("derived agegrp reproduces the banded summary: ", fn))
}

# The comparison is only meaningful if the rates actually vary across the
# cells being compared -- otherwise a broken collapse would still match.
# The median column is located by pattern rather than spelled out:
# scales::percent() picks its accuracy from the whole `prbl` vector, so the
# name is "50%" here but "50.0%" under the production default. Hard-coding
# either one makes these assertions pass vacuously against a NULL column.
median_col <- function(dt) {
  nm <- grep("_rate_50(\\.0)?%$", names(dt), value = TRUE)
  expect_equal(length(nm), 1L, info = "exactly one median rate column")
  nm[[1L]]
}
chk <- fread(file = file.path(
  by_age_dir, "all-cause mortality given disease-year-agegroup (not standardised).csv"))
expect_true(length(unique(chk[[median_col(chk)]])) > 1L,
            info = "fixture produces varying rates, so the test can fail")
expect_equal(sort(unique(chk$agegrp)), c("50-54", "65-69"),
             info = "derived bands reach the output, not single years")

# The population-denominator tables exercise the pairing of the numerator with
# `prvl` on the derived column: `prvl` is aggregated on the SAME strata before
# the join, so deriving on only one side fails in that group-by with the very
# error this helper exists to prevent.
pd <- fread(file = file.path(
  by_age_dir,
  "all-cause mortality given disease-year-agegroup popdenom (not standardised).csv"))
expect_true(nrow(pd) > 0L)
expect_false(anyNA(pd[[median_col(pd)]]),
             info = "numerator and denominator join on the derived agegrp")
expect_true(all(pd[[median_col(pd)]] > 0),
            info = "and the joined rate is a real number, not a zero-fill")

# ===========================================================================
# 6. It reaches the other families that default to agegrp strata
# ===========================================================================
# `ons` carries the same agegrp defaults and is reused verbatim by the CEA
# task, so the main tables must cope too. (Note this family keeps `agegrp` in
# its filenames; only the by-disease, xps and equity families spell it
# `agegroup`.)
prvl_main_age <- CJ(mc = 1:3, scenario = "sc0", year = 19:20,
                    age = c(50:54, 65:69), sex = c("men", "women"), sorted = FALSE)
prvl_main_age[, chd_prvl := as.numeric(age * 2L + mc)]
prvl_main_age[, popsize := as.numeric(age * 100L + 500L)]
prvl_main_band <- band(prvl_main_age, c("mc", "scenario", "year", "age", "sex"))

run_main <- function(src) {
  dir <- tempfile("agegrpmain_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  core <- Simulation$private_methods$tbl_smmrs_core
  environment(core) <- env
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      if (summary_type == "prvl" && standardization == "scaled_up") copy(src) else NULL
    },
    tbl_smmrs_core = core
  )
  f <- Simulation$private_methods$export_main_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975),
    baseline_year = 2019L, output_dir = dir, tables_dir = dir,
    comparator_scenario = "sc0", two_agegrps = FALSE,
    strata_ons = list(c("year", "agegrp")), strata_esp = list())
  dir
}

fn_main <- "prevalence by year-agegrp (not standardised).csv"
main_derived <- run_main(prvl_main_age)
main_native <- run_main(prvl_main_band)
expect_true(file.exists(file.path(main_derived, fn_main)),
            info = "main tables: agegrp strata survive single-year-of-age summaries")
expect_equal(sort(list.files(main_derived)), sort(list.files(main_native)))
expect_true(length(list.files(main_native, pattern = "\\.csv$")) >= 4L,
            info = "the comparison loop below is not empty")
for (fn in list.files(main_native, pattern = "\\.csv$")) {
  expect_equal(fread(file = file.path(main_derived, fn)),
               fread(file = file.path(main_native, fn)),
               info = paste0("main tables: derived == natively banded: ", fn))
}

# ===========================================================================
# 7. Ordering: the derivation must precede the two_agegrps collapse
# ===========================================================================
# collapse_two_agegrps() returns FALSE and does NOTHING when there is no
# `agegrp` column, so a family that collapses before it derives publishes the
# full 5-year bands into tables2agegrps/ with no error and no warning. Section
# 3 composes the two calls by hand and cannot catch that; these drive the real
# export paths, which is where the ordering actually lives.
acm_2ag <- run_acm(acm_age, prvl_age, list(c("year", "agegrp")),
                   two_agegrps = TRUE)
t2 <- fread(file = file.path(
  acm_2ag, "all-cause mortality given disease-year-agegroup (not standardised).csv"))
expect_equal(sort(unique(t2$agegrp)), c("50-64", "65-69"),
             info = "all_cause_mrtl: derived bands are then collapsed to two")
expect_equal(t2, fread(file = file.path(
  run_acm(acm_band, prvl_band, list(c("year", "agegrp")), two_agegrps = TRUE),
  "all-cause mortality given disease-year-agegroup (not standardised).csv")),
  info = "all_cause_mrtl: two_agegrps output matches the natively banded run")

run_main_2ag <- function(src) {
  dir <- tempfile("agegrpmain2_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  core <- Simulation$private_methods$tbl_smmrs_core
  environment(core) <- env
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      if (summary_type == "prvl" && standardization == "scaled_up") copy(src) else NULL
    },
    tbl_smmrs_core = core
  )
  f <- Simulation$private_methods$export_main_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975),
    baseline_year = 2019L, output_dir = dir, tables_dir = dir,
    comparator_scenario = "sc0", two_agegrps = TRUE, two_agegrps_split_age = 65L,
    strata_ons = list(c("year", "agegrp")), strata_esp = list())
  dir
}
m2 <- fread(file = file.path(run_main_2ag(prvl_main_age), fn_main))
expect_equal(sort(unique(m2$agegrp)), c("50-64", "65-69"),
             info = "main tables: derived bands are then collapsed to two")
expect_equal(m2, fread(file = file.path(run_main_2ag(prvl_main_band), fn_main)),
             info = "main tables: two_agegrps output matches the natively banded run")

# ===========================================================================
# 8. The equity family
# ===========================================================================
# `agegrp` is an ordinary output stratum here, and this is the family where a
# late derivation is silent rather than fatal -- collapse_two_agegrps() simply
# does nothing. It is also the family where the consequence is worst: the index
# is a fitted MODEL, so a two-group and a 5-year-band fit are different
# estimands, not two views of one number.
make_incd_age <- function() {
  d <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 19:21,
          dimd = factor(dimd_lv, levels = dimd_lv),
          age = c(50:54, 65:69), sex = c("men", "women"), sorted = FALSE)
  K <- length(dimd_lv)
  d[, base := (K + 1L - as.integer(dimd)) * 100 + mc + age]
  d[, chd_incd := fifelse(scenario == "sc0", base,
                          base * (1 - 0.02 * (K + 1L - as.integer(dimd))))]
  d[, popsize := 10000 * (K + 1L - as.integer(dimd))]
  d[, base := NULL]
  d[]
}
incd_age <- make_incd_age()
incd_band <- band(incd_age, c("mc", "scenario", "year", "dimd", "age", "sex"))

run_equity <- function(src, strata, two_agegrps = FALSE) {
  dir <- tempfile("agegrpeq_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE, diseases = NULL)))
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      if (summary_type == "incd") copy(src) else NULL
    }
  )
  f <- Simulation$private_methods$export_equity_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975),
    summaries_dir = dir, tables_dir = dir, comparator_scenario = "sc0",
    baseline_year = 2019L, ridit_reference = "comparator", strata = strata,
    two_agegrps = two_agegrps, two_agegrps_split_age = 65L)
  dir
}

fn_eq <- "equity cpp slope index by year-agegroup (not standardised).csv"
eq_derived <- fread(file = file.path(
  run_equity(incd_age, list(c("year", "agegrp"))), fn_eq))
expect_equal(sort(unique(eq_derived$agegrp)), c("50-54", "65-69"),
             info = "equity: agegrp is available as a stratum on age summaries")
expect_equal(eq_derived, fread(file = file.path(
  run_equity(incd_band, list(c("year", "agegrp"))), fn_eq)),
  info = "equity: derived == natively banded")

# ... and under two_agegrps the collapse must actually happen.
eq_2ag <- fread(file = file.path(
  run_equity(incd_age, list(c("year", "agegrp")), two_agegrps = TRUE), fn_eq))
expect_equal(sort(unique(eq_2ag$agegrp)), c("50-64", "65-69"),
             info = "equity: derivation precedes the two_agegrps collapse")
expect_false(isTRUE(all.equal(eq_2ag, eq_derived)),
             info = "and the collapse is not a silent no-op")
expect_equal(eq_2ag, fread(file = file.path(
  run_equity(incd_band, list(c("year", "agegrp")), two_agegrps = TRUE), fn_eq)),
  info = "equity: two_agegrps output matches the natively banded run")

# ===========================================================================
# 9. The CEA family (its strata default to `ons`, which names agegrp)
# ===========================================================================
make_cea <- function(what) {
  d <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 19:20,
          age = c(50:54, 65:69), sex = c("men", "women"), sorted = FALSE)
  d[, popsize := as.numeric(age * 100L + 500L)]
  if (what == "qalys") {
    d[, EQ5D5L := fifelse(scenario == "sc0", 0.7 + age / 1000, 0.72 + age / 1000)]
  } else {
    d[, healthcare_cost := as.numeric(age * 20L + mc)]
    d[, socialcare_cost := as.numeric(age * 5L)]
    d[, total_cost := healthcare_cost + socialcare_cost]
  }
  d[]
}
cea_keys <- c("mc", "scenario", "year", "age", "sex")
run_cea <- function(qalys, costs, two_agegrps = FALSE) {
  dir <- tempfile("agegrpcea_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      if (summary_type == "qalys") copy(qalys)
      else if (summary_type == "costs") copy(costs)
      else NULL
    }
  )
  f <- Simulation$private_methods$export_cea_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975), summaries_dir = dir, tables_dir = dir,
    comparator_scenario = "sc0", baseline_year = 2019L, wtp = 25000,
    two_agegrps = two_agegrps, two_agegrps_split_age = 65L,
    strata = list(c("year", "agegrp")))
  dir
}
cea_derived <- run_cea(make_cea("qalys"), make_cea("costs"))
cea_files <- list.files(cea_derived, pattern = "\\.csv$")
expect_true(length(cea_files) > 0L,
            info = "CEA tables are written against single-year-of-age summaries")
cea_native <- run_cea(band(make_cea("qalys"), cea_keys),
                      band(make_cea("costs"), cea_keys))
for (fn in cea_files) {
  expect_equal(fread(file = file.path(cea_derived, fn)),
               fread(file = file.path(cea_native, fn)),
               info = paste0("CEA: derived == natively banded: ", fn))
}
cea2 <- fread(file = file.path(run_cea(make_cea("qalys"), make_cea("costs"),
                                       two_agegrps = TRUE), cea_files[[1L]]))
expect_equal(sort(unique(cea2$agegrp)), c("50-64", "65-69"),
             info = "CEA: derivation precedes the two_agegrps collapse")

# ===========================================================================
# 10. Ages above the default top edge are banded, not wrapped
# ===========================================================================
# to_agegrp() derives its LABELS from max_age but its age vector from
# max(actual, max_age), so without the guard an age above 99 recycles the label
# vector and lands on the FIRST band -- age 100 published as "30-34".
d99 <- data.table(age = 30:104)
expect_silent(add_agegrp(d99))
expect_equal(d99[age == 100L, agegrp], "100-104")
expect_equal(d99[age == 99L, agegrp], "95-99",
             info = "and the bands at or below 99 are untouched")
