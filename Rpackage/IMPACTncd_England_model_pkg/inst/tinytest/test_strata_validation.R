# Parent-side validation of the `strata` argument of export_tables():
#   1. build_strata_config() rejects an element name that means nothing.
#   2. private$check_strata_columns() rejects a stratum the summaries cannot
#      supply, before any task is dispatched.
#
# Both exist for the same reason: `strata` is merged over the built-in defaults
# BY NAME, so a mistake in it does not announce itself -- the family you meant
# to change quietly keeps its defaults, and you learn about it from whatever
# those defaults do next. For `mrtl` instead of `mrtl_ons` on a
# single-year-of-age run, that was `object 'agegrp' not found` raised inside a
# forked worker minutes in.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

bsc <- Simulation$private_methods$build_strata_config
build <- function(user_strata, two_agegrps = FALSE) {
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  f <- bsc
  environment(f) <- env
  f(user_strata, two_agegrps)
}

VALID <- c("ons", "esp", "mrtl_ons", "mrtl_esp", "disease_char",
           "xps_ons", "xps_esp", "equity")

# ===========================================================================
# 1. Unknown element names
# ===========================================================================
expect_equal(sort(names(build(NULL))), sort(VALID),
             info = "the eight documented names are exactly what defaults holds")

# The near-miss that motivated this: `mrtl` is not a name, `mrtl_ons` is.
expect_error(build(list(ons = list("year"), mrtl = list("year"))),
             pattern = "unknown `strata` element name")
e <- tryCatch(build(list(mrtl = list("year"))), error = conditionMessage)
expect_true(grepl("Did you mean", e), info = "a near-miss is named, not just listed")
expect_true(grepl("mrtl_ons", e), info = "and the suggestion is the right one")

# A name with no plausible neighbour still errors, just without a suggestion.
e2 <- tryCatch(build(list(bananas = list("year"))), error = conditionMessage)
expect_true(grepl("unknown `strata` element name", e2))
expect_false(grepl("Did you mean", e2))

# An unnamed element cannot be dispatched to any family either.
expect_error(build(list(list("year"))), pattern = "must be a NAMED list")
expect_error(build(list(ons = list("year"), list("year"))),
             pattern = "must be a NAMED list")

# Valid names still merge, and only the named element is replaced.
r <- build(list(ons = list(c("year", "age"))))
expect_equal(r$ons, list(c("year", "age")))
expect_equal(r$mrtl_ons, build(NULL)$mrtl_ons,
             info = "an element you did not name keeps its default in full")
expect_silent(build(NULL))
for (nm in setdiff(VALID, "equity")) {
  expect_silent(build(setNames(list(list("year")), nm)))
}

# ===========================================================================
# 2. Requested columns must exist in the summaries
# ===========================================================================
# A fake summaries tree: one parquet per dataset, columns only -- the check
# reads schemas, never rows.
make_tree <- function(cols, datasets = c("prvl", "incd", "mrtl", "dis_mrtl",
                                         "qalys", "costs",
                                         "all_cause_mrtl_by_dis"),
                      std = "scaled_up") {
  root <- tempfile("stratacols_")
  for (d in datasets) {
    p <- file.path(root, "summaries", paste0(d, "_", std))
    dir.create(p, recursive = TRUE)
    dt <- as.data.table(setNames(
      lapply(cols, function(x) if (x == "age") 50L else "a"), cols))
    dt[, value := 1]
    arrow::write_parquet(dt, file.path(p, "0.parquet"))
  }
  root
}

check <- function(root, strata_cfg) {
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  env$private <- list(output_dir = function(x = "") file.path(root, x))
  f <- Simulation$private_methods$check_strata_columns
  environment(f) <- env
  env$private$check_strata_columns <- f
  f(strata_cfg)
}

banded <- make_tree(c("mc", "scenario", "year", "agegrp", "sex", "dimd"))
by_age <- make_tree(c("mc", "scenario", "year", "age", "sex", "dimd"))

# --- what must pass --------------------------------------------------------
expect_true(check(banded, build(NULL)),
            info = "the shipped defaults pass against a banded run")
expect_true(check(by_age, build(NULL)),
            info = "and against a single-year run, because agegrp is derived")
expect_true(check(by_age, build(list(ons = list(c("year", "age")),
                                     mrtl_ons = list(c("year", "age"))))),
            info = "`age` is legal when the summaries carry it")
expect_true(check(banded, build(list(ons = list(c("year", "qimd"))))),
            info = "`qimd` is legal when the summaries carry `dimd`")

# --- what must fail, in the parent -----------------------------------------
# The motivating reverse case: `age` against a banded run. Single years cannot
# be recovered from bands, so this can only fail -- the question is where.
e3 <- tryCatch(check(banded, build(list(ons = list(c("year", "age"))))),
               error = conditionMessage)
expect_true(grepl("do not carry", e3))
expect_true(grepl("'age'", e3), info = "the offending column is named")
expect_true(grepl("strata\\$ons", e3), info = "and so is the family")
expect_true(grepl("single_year_of_age", e3),
            info = "the remedy is a different export_summaries() call, so say so")

# A stratum the design never asked for.
nodimd <- make_tree(c("mc", "scenario", "year", "agegrp", "sex"))
e4 <- tryCatch(check(nodimd, build(list(mrtl_ons = list(c("year", "dimd"))))),
               error = conditionMessage)
expect_true(grepl("'dimd'", e4))
expect_true(grepl("strata\\$mrtl_ons", e4))
expect_false(grepl("single_year_of_age", e4),
             info = "the age hint is not appended to unrelated failures")
# ... and qimd is not available either, since it is derived FROM dimd.
expect_error(check(nodimd, build(list(ons = list(c("year", "qimd"))))),
             pattern = "'qimd'")

# A plain typo in a column name.
expect_error(check(banded, build(list(ons = list(c("year", "agegroup"))))),
             pattern = "'agegroup'")

# --- conservatism: never fire where the worker would have coped ------------
# A family with no datasets on disk is not checked at all.
empty <- tempfile("stratacols_empty_")
dir.create(file.path(empty, "summaries"), recursive = TRUE)
expect_true(check(empty, build(list(ons = list(c("year", "banana"))))),
            info = "no summaries yet -> nothing to validate against")
# A family whose own summary is missing is skipped even when others are present.
only_acm <- make_tree(c("mc", "scenario", "year", "agegrp", "sex", "dimd"),
                      datasets = "all_cause_mrtl_by_dis")
expect_true(check(only_acm, build(list(ons = list(c("year", "nonesuch"))))),
            info = "ons reads none of the datasets present, so it is not judged")
expect_error(check(only_acm, build(list(mrtl_ons = list(c("year", "nonesuch"))))),
             pattern = "'nonesuch'")
# The ESP element is judged against the _esp datasets, not the scaled_up ones.
esp_only <- make_tree(c("mc", "scenario", "year", "sex", "dimd"), std = "esp")
expect_error(check(esp_only, build(list(esp = list(c("year", "agegrp"))))),
             pattern = "'agegrp'")
expect_true(check(esp_only, build(list(ons = list(c("year", "banana"))))),
            info = "scaled_up datasets are absent, so `ons` is not judged")

# The families this check deliberately does not police.
expect_true(check(banded, build(list(xps_ons = list(c("year", "agegrp20"))))),
            info = "xps tables are not built from summaries/")
expect_true(check(banded, build(list(equity = list(c("year", "ethnicity"))))),
            info = "equity keeps its own warn-and-skip, pinned by its own tests")

# ===========================================================================
# 3. The derivation rules used here match the ones the workers apply
# ===========================================================================
drv <- IMPACTncdEngland:::.derivable_strata_columns
expect_equal(sort(drv(c("year", "dimd"))), c("dimd", "qimd", "year"))
expect_equal(sort(drv(c("year", "age"))), c("age", "agegrp", "year"))
expect_equal(sort(drv(c("year", "age", "dimd"))),
             c("age", "agegrp", "dimd", "qimd", "year"))
expect_equal(drv(c("year", "sex")), c("year", "sex"),
             info = "nothing is invented out of thin air")
expect_equal(sort(drv(c("year", "agegrp"))), c("agegrp", "year"),
             info = "bands do not yield single years -- the derivation is one-way")
