# End-to-end tests for private$export_equity_tables(): which files it writes,
# what gradient axis each is fit over, the dimd -> qimd derivation, and the
# warnings raised when a requested gradient is unavailable.
#
# The method is exercised directly against an in-memory fixture rather than a
# simulation run: its only couplings to the R6 object are
# `self$design$sim_prm$logs` and `private$read_summary_dataset()`, so it is
# re-parented onto a mock environment whose parent is the package namespace
# (keeping calc_equity_slope_indices(), safe_fquantile_byid() etc. resolvable).

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

dimd_lv <- c("1 most deprived", as.character(2:9), "10 least deprived")
qimd_lv <- c("1 most deprived", as.character(2:4), "5 least deprived")

# --- Fixture: incd summary, 3 MC x 2 scenarios x 3 years x 10 deciles -------
# Deprivation-patterned incidence with a pro-poor intervention effect, so the
# slope index is well determined and non-zero.
make_incd <- function(dep_col = "dimd", levs = dimd_lv) {
  d <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 19:21,
          dep = factor(levs, levels = levs),
          agegrp = c("50-54", "55-59"), sex = c("men", "women"),
          sorted = FALSE)
  K <- length(levs)
  # baseline events rise with deprivation; sc1 removes a share of them that is
  # itself larger in the more deprived groups (pro-poor)
  d[, base := (K + 1L - as.integer(dep)) * 100 + mc]
  d[, chd_incd := fifelse(scenario == "sc0", base,
                          base * (1 - 0.02 * (K + 1L - as.integer(dep))))]
  d[, popsize := 10000 * (K + 1L - as.integer(dep))]
  d[, base := NULL]
  setnames(d, "dep", dep_col)
  d[]
}

# --- Mock harness -----------------------------------------------------------
run_export <- function(strata, sources = list(incd = make_incd()), logs = FALSE) {
  # tempfile() is guaranteed unique within a session; do NOT derive the name
  # from the RNG, or two runs can share a directory and read each other's files.
  dir <- tempfile("eqtest_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = logs)))
  env$private <- list(
    # the method mutates tt by reference, so hand out a fresh copy each call
    read_summary_dataset = function(summary_type, standardization) {
      if (is.null(sources[[summary_type]])) NULL else copy(sources[[summary_type]])
    }
  )
  f <- Simulation$private_methods$export_equity_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),   # the export_tables() default
    summaries_dir = dir, tables_dir = dir,
    comparator_scenario = "sc0", baseline_year = 2019L,
    ridit_reference = "comparator", strata = strata)
  dir
}

files_in <- function(dir) sort(basename(list.files(dir, pattern = "\\.csv$")))

# Subsetting a data.table with `==` leaves an auto-index behind, which
# all.equal.data.table reports as a difference. Strip it before comparing.
noidx <- function(x) { setindex(x, NULL); x[] }

# ===========================================================================
# 1. Defaults are unchanged: implicit dimd, historical filenames
# ===========================================================================
d1 <- run_export(list("year", c("year", "sex")))
expect_equal(files_in(d1),
             c("equity cpp slope index by year (not standardised).csv",
               "equity cpp slope index by year-sex (not standardised).csv"),
             info = "implicit gradient keeps the historical filenames")

t1 <- fread(file.path(d1, "equity cpp slope index by year (not standardised).csv"))
expect_true("gradient" %in% names(t1), info = "gradient column is written")
expect_equal(unique(t1$gradient), "dimd")
expect_equal(names(t1)[1:5],
             c("gradient", "scenario", "year", "disease", "type"),
             info = "gradient leads the column order")
expect_equal(sort(unique(t1$type)),
             c("AEI_per100k", "AEI_total", "REI_rel", "RII_ratio"))
expect_true(all(t1[type == "AEI_per100k", `equity_50.0%`] > 0),
            info = "pro-poor fixture gives a positive absolute index")
expect_false("sex" %in% names(t1))
expect_true("sex" %in% names(fread(file.path(
  d1, "equity cpp slope index by year-sex (not standardised).csv"))))

# ===========================================================================
# 2. Case 2 from the report: an explicit qimd gradient actually works
# ===========================================================================
d2 <- run_export(list(c("year", "qimd")))
expect_equal(files_in(d2),
             "equity cpp slope index by year-qimd (not standardised).csv",
             info = "the qimd token is echoed in the filename")
t2 <- fread(file.path(d2, files_in(d2)))
expect_equal(unique(t2$gradient), "qimd")

# It must NOT be the old silent duplicate of the year-only table: the fixture's
# gradient is non-linear across deciles, so the quintile index really differs.
expect_false(isTRUE(all.equal(t1[, -"gradient"], t2[, -"gradient"])),
             info = "qimd gradient is genuinely refit, not a relabelled copy")
# ...but it is the same order of magnitude and the same (pro-poor) sign.
expect_true(all(t2[type == "AEI_per100k", `equity_50.0%`] > 0))
expect_equal(nrow(t1), nrow(t2),
             info = "same output rows: the gradient is consumed, not stratified")

# ===========================================================================
# 3. Both tokens in one entry -> one table per gradient
# ===========================================================================
d3 <- run_export(list(c("year", "dimd", "qimd"), c("year", "sex", "qimd")))
expect_equal(files_in(d3),
             c("equity cpp slope index by year-dimd (not standardised).csv",
               "equity cpp slope index by year-qimd (not standardised).csv",
               "equity cpp slope index by year-sex-qimd (not standardised).csv"))
expect_equal(unique(fread(file.path(d3, files_in(d3)[1]))$gradient), "dimd")
expect_equal(unique(fread(file.path(d3, files_in(d3)[3]))$gradient), "qimd")
# The explicit dimd table matches the implicit one bar the filename.
expect_equal(fread(file.path(d3, files_in(d3)[1])), noidx(t1),
             info = "explicit dimd == implicit dimd")

# ===========================================================================
# 4. The derived qimd equals a natively-aggregated qimd summary, exactly
# ===========================================================================
native <- make_incd()
native[, qimd := factor(qimd_lv[c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L)][
  as.integer(dimd)], levels = qimd_lv)]
native[, dimd := NULL]
native <- native[, .(chd_incd = sum(chd_incd), popsize = sum(popsize)),
                 by = .(mc, scenario, year, qimd, agegrp, sex)]
d4 <- run_export(list(c("year", "qimd")), sources = list(incd = native))
expect_equal(fread(file.path(d4, files_in(d4))), noidx(t2),
             info = "qimd derived from dimd == qimd aggregated at summary time")

# ===========================================================================
# 5. Unavailable gradients warn (not a logs-gated message) and skip
# ===========================================================================
# qimd-only summaries: deciles cannot be recovered, so a dimd request warns.
qonly <- copy(native)
d5 <- NULL
expect_warning(d5 <- run_export(list(c("year", "dimd"), c("year", "qimd")),
                                sources = list(incd = qonly)),
               pattern = "gradient 'dimd' was requested")
expect_equal(files_in(d5),
             "equity cpp slope index by year-qimd (not standardised).csv",
             info = "the available gradient is still written")

# No deprivation column at all -> warn once per metric, write nothing.
nodep <- make_incd()
nodep[, dimd := NULL]
nodep <- nodep[, .(chd_incd = sum(chd_incd), popsize = sum(popsize)),
               by = .(mc, scenario, year, agegrp, sex)]
d6 <- NULL
expect_warning(d6 <- run_export(list("year"), sources = list(incd = nodep)),
               pattern = "no deprivation column")
expect_equal(length(files_in(d6)), 0L)

# Missing comparator scenario -> warn, write nothing.
nocmp <- make_incd()[scenario != "sc0"]
d7 <- NULL
expect_warning(d7 <- run_export(list("year"), sources = list(incd = nocmp)),
               pattern = "comparator scenario 'sc0'")
expect_equal(length(files_in(d7)), 0L)

# ===========================================================================
# 6. Any variable the summaries carry works as an output stratum
# ===========================================================================
d8 <- run_export(list(c("year", "agegrp"), c("year", "agegrp", "sex", "qimd")))
expect_equal(files_in(d8),
             c("equity cpp slope index by year-agegroup (not standardised).csv",
               "equity cpp slope index by year-agegroup-sex-qimd (not standardised).csv"),
             info = "agegrp is spelled agegroup in the filename, as elsewhere")
t8 <- fread(file.path(d8, files_in(d8)[1]))
expect_true("agegrp" %in% names(t8))
expect_equal(sort(unique(t8$agegrp)), c("50-54", "55-59"))
expect_equal(nrow(t8), 2L * nrow(noidx(t1)),
             info = "one gradient fit per age group")
expect_equal(unique(t8$gradient), "dimd", info = "no token -> dimd, as ever")

t8b <- fread(file.path(d8, files_in(d8)[2]))
expect_equal(unique(t8b$gradient), "qimd")
expect_true(all(c("agegrp", "sex") %in% names(t8b)))
expect_equal(names(t8b)[1:6],
             c("gradient", "scenario", "year", "agegrp", "sex", "disease"))

# The defining invariant: stratifying by X must equal filtering to a level of X
# and not stratifying. (This is what silently failed before -- the extra
# variable was dropped, so the "stratified" table was the unstratified one.)
d9 <- run_export(list("year"),
                 sources = list(incd = make_incd()[agegrp == "50-54"]))
expect_equal(t8[agegrp == "50-54", -"agegrp"],
             noidx(fread(file.path(d9, files_in(d9)))),
             info = "stratify-by-agegrp == filter-to-agegrp then fit")

# A stratum the summaries do not carry warns and skips (it must not silently
# collapse to a mislabelled copy of a coarser table).
d10 <- NULL
expect_warning(d10 <- run_export(list(c("year", "ethnicity"), "year")),
               pattern = "stratification variable")
expect_equal(files_in(d10),
             "equity cpp slope index by year (not standardised).csv",
             info = "the usable stratum is still written")

# Reserved columns are rejected up front.
expect_error(run_export(list(c("year", "mc"))), pattern = "reserved column")
expect_error(run_export(list(c("year", "scenario"))), pattern = "reserved column")
expect_error(run_export(list(c("agegrp", "sex"))), pattern = "must include")

# ===========================================================================
# 7. build_strata_config(): the wiring from export_tables(strata = ...)
# ===========================================================================
# Validation must fire here, i.e. up front in export_tables(), not deep inside
# a parallel worker after the other tables have already been built.
bsc <- Simulation$private_methods$build_strata_config

expect_equal(bsc(NULL, FALSE)$equity, list("year", c("year", "sex")),
             info = "default equity strata are unchanged")
expect_equal(bsc(NULL, TRUE)$equity, list("year", c("year", "sex")),
             info = "two_agegrps defaults are unchanged too")
expect_equal(bsc(list(equity = list(c("year", "qimd"))), FALSE)$equity,
             list(c("year", "qimd")),
             info = "a user equity list overrides the default")
expect_true("dimd" %in% unlist(bsc(list(equity = list(c("year", "qimd"))),
                                   FALSE)$ons),
            info = "overriding equity leaves the other strata lists alone")
expect_silent(bsc(list(equity = list(c("year", "agegrp", "qimd"))), FALSE))
expect_error(bsc(list(equity = list(c("year", "dimd"), "sex")), FALSE),
             pattern = "must include")
expect_error(bsc(list(equity = list(c("year", "mc"))), FALSE),
             pattern = "reserved column")
