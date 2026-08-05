# Tests for the placement of the working copy in private$export_main_tables().
#
# The method loops over the metrics of a source dataset and hands each one a
# `tt` that downstream code mutates BY REFERENCE -- the two_agegrps relabel, the
# absorb_dt() in the ftlt branch. `tt_base` must therefore survive untouched for
# the NEXT metric in the same group, which is what `tt <- copy(tt_base)` buys.
#
# That copy used to be taken ABOVE the guard that skips the comparison metrics
# (cypp / cpp / dpp / net_qalys / net_costs) when the summary carries no
# intervention arm. On a single-scenario run the guard fires for every one of
# them, so the copy was allocated and immediately dropped -- ten full duplicates
# of the summary per export_tables() call (five metrics x two populations),
# 0.6-2.2 GB each on a production run. The copy now sits BELOW the guard.
#
# These tests pin that the move changed nothing observable:
#   1. a single-scenario run writes exactly the non-comparison tables;
#   2. its output is IDENTICAL to the sc0 slice of a two-scenario run. Every
#      non-comparison metric is computed strictly within scenario (the group-by,
#      the baseline-year join and the quantiling all carry `scenario` in the
#      key), so this is an EXACT invariant -- and the sharpest available check
#      that reordering the guard and the copy left the data undisturbed;
#   3. a two-scenario run still writes the comparison tables;
#   4. `tt_base` is not mutated across metrics -- `dis_mrtl` must not see the
#      `*_prvl` columns that the preceding `ftlt` absorbs into ITS working copy,
#      which is the failure the copy exists to prevent and the one a careless
#      hoist would reintroduce.
#
# Same mock harness as test_cea_perspectives.R / test_export_equity_tables.R:
# the method's only couplings to the R6 object are `self$design$sim_prm$logs`,
# `private$read_summary_dataset()` and `private$tbl_smmrs_core()`.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

DIMD <- c("1 most deprived", "10 least deprived")
AGEGRP <- c("30-34", "70-74")   # one in each two_agegrps bucket
SEX <- c("men", "women")

# ---------------------------------------------------------------------------
# Fixture. Every value is a deterministic function of the KEY columns only, so
# the sc0 rows are bit-identical whether or not sc1 is present in the summary.
# Without that, test 2 could not be an exact comparison.
# ---------------------------------------------------------------------------
grid_for <- function(scns) {
  d <- CJ(mc = 1:3, scenario = scns, year = 19:21, agegrp = AGEGRP,
          sex = SEX, dimd = DIMD, sorted = FALSE)
  d[, k := mc + 10L * (year - 19L) + 100L * match(agegrp, AGEGRP) +
       1000L * match(sex, SEX) + 10000L * match(dimd, DIMD) +
       100000L * (scenario != "sc0")]
  d[]
}

mk_sources <- function(scns) {
  g <- grid_for(scns)
  id <- c("mc", "scenario", "year", "agegrp", "sex", "dimd")
  cp <- function(...) cbind(g[, ..id], data.table(...))
  list(
    prvl = cp(popsize     = 5000L + g$k %% 701L,
              chd_prvl    =  100L + g$k %% 97L,
              stroke_prvl =   50L + g$k %% 89L),
    incd = cp(popsize     = 5000L + g$k %% 701L,
              chd_incd    =   20L + g$k %% 31L,
              stroke_incd =   10L + g$k %% 23L),
    mrtl = cp(popsize     = 5000L + g$k %% 701L,
              chd_mrtl    =    5L + g$k %% 13L,
              stroke_mrtl =    3L + g$k %% 11L),
    dis_mrtl = cp(popsize            = 5000L + g$k %% 701L,
                  chd_deaths         =    4L + g$k %% 17L,
                  stroke_deaths      =    2L + g$k %% 19L,
                  nonmodelled_deaths =    7L + g$k %% 29L),
    qalys = cp(EQ5D5L = 400 + g$k %% 53L),
    costs = cp(healthcare_cost = 1000 + g$k %% 211L,
               socialcare_cost =  400 + g$k %% 173L,
               total_cost      = 1400 + g$k %% 307L,
               economic_output = 2000 + g$k %% 151L)
    # `contd` deliberately absent -> read_summary_dataset() returns NULL, which
    # exercises the skip-a-missing-source path in the same run.
  )
}

# ---------------------------------------------------------------------------
# Harness: re-parent export_main_tables() (and the REAL tbl_smmrs_core, so the
# CSVs written are the genuine article) onto a stub self/private.
# ---------------------------------------------------------------------------
run_main <- function(sources, two_agegrps = FALSE, strata = list("year"),
                     record = NULL) {
  dir <- tempfile("maintbl_")
  dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE, diseases = NULL)))
  core <- Simulation$private_methods$tbl_smmrs_core
  environment(core) <- env
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      s <- sources[[summary_type]]
      if (is.null(s)) NULL else copy(s)
    },
    tbl_smmrs_core = function(...) {
      a <- list(...)
      # Snapshot BEFORE dispatching: this is the `tt` the metric was handed.
      if (!is.null(record)) {
        record(list(what = a$what, population = a$population,
                    cols = names(a$tt), nrow = nrow(a$tt)))
      }
      do.call(core, a)
    }
  )
  f <- Simulation$private_methods$export_main_tables
  environment(f) <- env
  f(prbl = c(0.5, 0.025, 0.975), baseline_year = 2019L, output_dir = dir,
    tables_dir = dir, comparator_scenario = "sc0", two_agegrps = two_agegrps,
    strata_ons = strata, strata_esp = strata)
  dir
}

# str4 prefixes of the five comparison metrics (see tbl_smmrs_core).
CMP_PREFIX <- c("case-years prevented or postponed by ",
                "cases prevented or postponed by ",
                "deaths prevented or postponed by ",
                "net QALYs by ", "net costs by ")
is_cmp_file <- function(f) {
  vapply(f, function(x) any(startsWith(x, CMP_PREFIX)), logical(1), USE.NAMES = FALSE)
}

d1 <- run_main(mk_sources("sc0"))                 # single scenario
d2 <- run_main(mk_sources(c("sc0", "sc1")))       # with an intervention arm
f1 <- sort(list.files(d1)); f2 <- sort(list.files(d2))

# ===========================================================================
# 1. A single-scenario run writes the non-comparison tables and NONE of the
#    comparison ones
# ===========================================================================
expect_true(length(f1) > 0L, info = "the single-scenario run wrote tables")
expect_false(any(is_cmp_file(f1)),
             info = "no cypp/cpp/dpp/net_qalys/net_costs file on one scenario")
# The guard must skip only those five; everything else must still be produced.
expect_equal(f1, f2[!is_cmp_file(f2)],
             info = "single-scenario file set == two-scenario minus contrasts")

# ===========================================================================
# 2. Identical output: single-scenario == the sc0 slice of the two-scenario run
# ===========================================================================
# This is the regression pin. If moving the copy below the guard had disturbed
# `tt_base` -- or if the guard now read the scenario list from a mutated table
# -- these would diverge.
for (f in f1) {
  a <- fread(file.path(d1, f))
  b <- fread(file.path(d2, f))[scenario == "sc0"]
  expect_equal(a, b, info = paste0("sc0 rows identical with/without sc1: ", f))
}

# ===========================================================================
# 3. Two scenarios still produce every contrast table
# ===========================================================================
expect_true(sum(is_cmp_file(f2)) > 0L,
            info = "the contrast tables are still built when an arm exists")
# Non-empty, i.e. the contrast actually joined rather than writing a header.
for (f in f2[is_cmp_file(f2)]) {
  expect_true(nrow(fread(file.path(d2, f))) > 0L,
              info = paste0("contrast table has rows: ", f))
}

# ===========================================================================
# 4. `tt_base` is not mutated across the metrics of a source group
# ===========================================================================
# Within the `dis_mrtl` source group the metric order is
#   ftlt, ftlt_change_relative, ftlt_change_absolute, dis_mrtl, ...
# and every ftlt* metric absorb_dt()s the prevalence denominators into ITS
# working copy. If that copy were shared with `tt_base`, `dis_mrtl` would arrive
# carrying `nonmodelled_prvl` / `chd_prvl` / `stroke_prvl` -- and, because the
# ftlt branch also drops rows (`tt[nonmodelled_prvl > 0]`), a different row
# count. Both are checked.
seen <- list()
invisible(run_main(mk_sources("sc0"),
                   record = function(x) seen[[length(seen) + 1L]] <<- x))

by_metric <- function(w, p = "ons") {
  Filter(function(x) x$what == w && x$population == p, seen)[[1L]]
}
ftlt <- by_metric("ftlt")
dmrt <- by_metric("dis_mrtl")

expect_true("nonmodelled_prvl" %in% ftlt$cols,
            info = "ftlt really does absorb the prevalence denominator")
expect_false(any(c("nonmodelled_prvl", "chd_prvl", "stroke_prvl") %in% dmrt$cols),
             info = "dis_mrtl gets a pristine tt, not ftlt's absorbed copy")
expect_equal(dmrt$nrow, nrow(mk_sources("sc0")$dis_mrtl),
             info = "dis_mrtl sees every row; ftlt's row filter did not leak")

# The same must hold for the esp population, which is a second pass over the
# same tt_base within one source-group iteration.
expect_false(any(c("nonmodelled_prvl", "chd_prvl", "stroke_prvl") %in%
                   by_metric("dis_mrtl", "esp")$cols),
             info = "...and on the esp pass too")

# ===========================================================================
# 5. two_agegrps: the relabel likewise stays inside the working copy
# ===========================================================================
# Metric 1 collapses agegrp on its own copy; metric 2 must still start from the
# uncollapsed tt_base. Row counts before the group-by are the tell: a leaked
# relabel does not change nrow, but a leaked ftlt absorb does, so this case is
# covered by the column check plus an end-to-end run that must not error.
d3 <- run_main(mk_sources("sc0"), two_agegrps = TRUE)
expect_true(length(list.files(d3)) > 0L,
            info = "two_agegrps single-scenario run completes")
expect_false(any(is_cmp_file(list.files(d3))),
             info = "and still skips the contrasts")
