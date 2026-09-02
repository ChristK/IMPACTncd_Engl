# Exposure DISTRIBUTIONS, and the pack-years denominator.
#
# Two changes to export_xps(), both about what a mean cannot tell you.
#
# A. `export_xps_quantiles: yes` adds weighted quantiles of every exposure
#    (`xps20_qnt` / `xps5_qnt`) next to the means, so a scenario that shifts a
#    whole distribution can be told apart from one that only moves its mean.
#    Off by default: it roughly sextuples export_xps(), which runs once per mc
#    chunk per scenario.
#
# B. Pack-years, smoking duration and cigarettes/day are set to NA for
#    never-smokers before aggregation, so the published means are CONDITIONAL
#    ON EVER HAVING SMOKED. That is a defensible denominator, but it makes a
#    scenario that moves people into never-smoking look as though smoking got
#    worse -- the remaining ever-smokers are a more heavily-smoking selection,
#    so their mean pack-years RISES while the population total falls.
#    `*_pop_curr_xps` companions count never-smokers as the zero they are, so
#    both readings are available. The test that matters is the one asserting
#    the two genuinely disagree in that direction.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

wq <- IMPACTncdEngland::weighted_quantile

# ===========================================================================
# 1. weighted_quantile()
# ===========================================================================
x <- c(1, 2, 3, 4, 5)
p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
expect_equal(wq(x, rep(1, 5), p),
             stats::quantile(x, p, type = 1, names = FALSE),
             info = "equal weights reproduce the type-1 (inverse CDF) quantile")

# Replication invariance: integer weights must equal repeating the values. This
# is the property that makes it a *population* quantile rather than a sample one.
set.seed(1)
v <- rnorm(200)
cnt <- sample(1:5, 200, TRUE)
expect_equal(wq(v, cnt, p),
             stats::quantile(rep(v, cnt), p, type = 1, names = FALSE),
             info = "weights == counts is the same as replicating the values")

# Weighting actually moves the answer (otherwise the above proves nothing).
expect_true(wq(x, c(1, 1, 1, 50, 50), 0.5) > wq(x, rep(1, 5), 0.5),
            info = "mass on the upper values pulls the median up")

# No interpolation: a discrete exposure must not gain values it cannot take.
expect_true(all(wq(0:7, rep(1, 8), c(0.3, 0.6, 0.9)) %in% 0:7),
            info = "results are always observed values")

# Edges and degenerate input.
expect_equal(wq(x, rep(1, 5), 1 - 1e-12), max(x))
expect_equal(wq(x, rep(1, 5), 1e-12), min(x))
expect_equal(wq(7, 1, c(0.1, 0.9)), c(7, 7))
expect_true(is.na(wq(x, rep(0, 5), 0.5)), info = "no positive weight -> NA")
expect_true(is.na(wq(rep(NA_real_, 5), rep(1, 5), 0.5)))
expect_equal(length(wq(x, rep(1, 5), p)), length(p))
# A missing weight drops that observation rather than poisoning the result.
expect_equal(wq(c(1, 2, NA, 4), c(1, 1, 1, 1), 0.5),
             wq(c(1, 2, 4), c(1, 1, 1), 0.5))

# ===========================================================================
# 1b. The hoisted core is the SAME function, not a second implementation
# ===========================================================================
# export_xps() calls .weighted_quantile_pre() directly, having hoisted the
# weight cleaning out of its per-exposure loop (~20% of the cost, since the
# weights do not depend on the exposure). weighted_quantile() is a thin wrapper
# over the same core -- so the only thing that can drift is the hoisting, and
# that is what these fuzz cases check.
pre <- IMPACTncdEngland:::.weighted_quantile_pre
hoist <- function(x, w, probs) {
  wsel <- !is.na(w) & w > 0
  if (all(wsel)) { wsel <- NULL; W <- sum(w) } else { W <- sum(w[wsel]) }
  pre(x, w, wsel, W, probs)
}

set.seed(99)
for (i in 1:300) {
  n <- sample(1:60, 1L)
  xx <- rnorm(n)
  ww <- runif(n, 0, 2)
  # inject the awkward cases the hoisting has to survive
  if (i %% 3L == 0L) xx[sample(n, max(1L, n %/% 4L))] <- NA_real_
  if (i %% 4L == 0L) ww[sample(n, max(1L, n %/% 5L))] <- 0
  if (i %% 5L == 0L) ww[sample(n, max(1L, n %/% 6L))] <- NA_real_
  if (i %% 7L == 0L) xx[] <- NA_real_          # nothing usable
  if (i %% 11L == 0L) ww[] <- 0                # no positive weight
  if (i %% 13L == 0L) xx <- round(xx)          # heavy ties
  expect_identical(hoist(xx, ww, p), wq(xx, ww, p),
                   info = paste("hoisted == public, fuzz case", i))
}

# The two paths through the wrapper (all weights usable vs not) must agree
# where they should: a zero weight contributes nothing, so dropping it cannot
# change the answer.
xa <- c(1, 2, 3, 4); wa <- c(1, 1, 1, 1)
expect_equal(wq(c(xa, 99), c(wa, 0), p), wq(xa, wa, p),
             info = "a zero-weight observation is inert")

# ===========================================================================
# 2. The design flag: off unless asked for, and validated
# ===========================================================================
# Guarded so `export_xps_quantiles` cannot be a typo that silently does nothing.
dc <- paste(deparse(Design$public_methods$initialize), collapse = "\n")
expect_true(grepl("export_xps_quantiles", dc, fixed = TRUE),
            info = "defaulted in Design, so old YAMLs keep loading")
expect_true(grepl("xps_quantile_probs", dc, fixed = TRUE))

xps_txt <- paste(deparse(Simulation$private_methods$export_xps), collapse = "\n")
expect_true(grepl("export_xps_quantiles", xps_txt, fixed = TRUE),
            info = "export_xps consults the flag")
expect_true(grepl("weighted_quantile", xps_txt, fixed = TRUE),
            info = "...and uses the weighted definition, not stats::quantile")
expect_true(grepl("xps20_qnt", xps_txt, fixed = TRUE))
expect_true(grepl("xps5_qnt", xps_txt, fixed = TRUE))

# ===========================================================================
# 3. Pack-years: the ever-smoker mean and the population mean disagree,
#    in the direction that caused the confusion
# ===========================================================================
# Two worlds with the SAME ever-smokers, differing only in how many
# never-smokers there are -- i.e. exactly a "fewer people start smoking"
# scenario, with no change to any smoker's behaviour.
mk_world <- function(n_never) {
  ever <- data.table(smok_status_curr_xps = "4",
                     smok_packyrs_curr_xps = c(20, 30, 40, 50),
                     wt = 1)
  never <- data.table(smok_status_curr_xps = "1",
                      smok_packyrs_curr_xps = 0, wt = 1)[rep(1L, n_never)]
  rbind(ever, never)
}
baseline <- mk_world(4L)
scenario <- mk_world(20L)   # many more never-smokers, same smokers

# What export_xps publishes for `smok_packyrs_curr_xps`: NA for never-smokers,
# then weighted.mean(na.rm = TRUE).
ever_mean <- function(d) {
  d <- copy(d)
  d[smok_status_curr_xps == "1", smok_packyrs_curr_xps := NA]
  d[, weighted.mean(smok_packyrs_curr_xps, wt, na.rm = TRUE)]
}
# What the `_pop_` companion publishes: never-smokers count as 0.
pop_mean <- function(d) d[, weighted.mean(smok_packyrs_curr_xps, wt)]

expect_equal(ever_mean(baseline), ever_mean(scenario),
             info = "the ever-smoker mean is BLIND to the scenario's effect")
expect_true(pop_mean(scenario) < pop_mean(baseline),
            info = "the population mean shows the fall the scenario produced")
expect_true(pop_mean(baseline) / pop_mean(scenario) > 2,
            info = "and the gap is large, not marginal")

# The companions are built before the NAs are introduced, are picked up by the
# same `_curr_xps$` grep, and do not leak into the lifecourse parquet.
expect_true(grepl("_pop_curr_xps", xps_txt, fixed = TRUE),
            info = "export_xps creates the population companions")
expect_true(grepl("smok_packyrs_pop_curr_xps", xps_txt, fixed = TRUE))
expect_true(grepl("smok_quit_yrs", xps_txt, fixed = TRUE))
# quit-years deliberately has NO companion: a population mean over
# never-smokers and current smokers alike would be meaningless.
expect_false(grepl("smok_quit_yrs_pop_curr_xps", xps_txt, fixed = TRUE),
             info = "no companion for smok_quit_yrs, on purpose")

# ===========================================================================
# 4. The quantile TABLES, end to end through export_xps_tables()
# ===========================================================================
# A synthetic xps20_qnt tree in the layout export_xps() writes, including the
# "All" marginals groupingsets produces.
PRBL <- c(0.5, 0.025, 0.975, 0.1, 0.9)
PROBS <- c(0.1, 0.5, 0.9)

mk_qnt <- function() {
  full <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 2019:2020,
             agegrp20 = c("30-49", "50-69"), sex = c("men", "women"),
             qimd = c("1 most deprived", "5 least deprived"),
             xps = c("bmi_curr_xps", "sbp_curr_xps"), prob = PROBS,
             sorted = FALSE)
  full[, value := 25 + 10 * prob + mc * 0.1 +
         2 * (scenario == "sc1") + (sex == "women")]
  # The year-only marginal, as groupingsets emits it: a POPULATION-WEIGHTED
  # aggregate over unequally-sized cells, so it does NOT equal the unweighted
  # mean of the cell values. That distinction is the whole point -- if the
  # marginal were the plain mean of its cells, a table built from the
  # unfiltered dataset would coincide with the correct one and the check below
  # could not detect a missing filter.
  full[, .w := 1 + as.integer(factor(paste(agegrp20, sex, qimd)))]
  marg <- full[, .(value = weighted.mean(value, .w)),
               keyby = .(mc, scenario, year, xps, prob)]
  marg[, `:=`(agegrp20 = "All", sex = "All", qimd = "All")]
  full[, .w := NULL]
  d <- rbind(full, marg, use.names = TRUE)
  # FACTORS, not character. export_xps() writes these strata columns from
  # sp$pop, where they are factors, so parquet stores them dictionary-encoded
  # and read_parquet_dt() brings them back as factors. A character fixture is
  # what let a real bug through: detect_filterable_vars() tested only
  # is.character(), so in production it found nothing, no filtering happened,
  # and every stratification's file got the whole dataset.
  for (cc in c("agegrp20", "sex", "qimd")) set(d, NULL, cc, factor(d[[cc]]))
  d[]
}
QNT <- mk_qnt()

run_xps_tab <- function(strata_ons, esp_strata = list("year")) {
  out <- tempfile("xpsq_"); dir.create(out, recursive = TRUE)
  tab <- tempfile("tab_"); dir.create(tab, recursive = TRUE)
  dir.create(file.path(out, "xps", "xps20_qnt"), recursive = TRUE)
  # Written with FULL years, as export_xps() does -- these datasets are the one
  # family that is NOT in short form, and adding 2000 again published years
  # 4013-4043.
  arrow::write_dataset(QNT, file.path(out, "xps", "xps20_qnt"),
                       format = "parquet")
  # The ESP variant too, so the standardisation LABEL path is exercised. It
  # has no agegrp20 (xps5 standardises over age), which is exactly why `age`
  # has to be added to the label vocabulary unconditionally.
  dir.create(file.path(out, "xps", "xps5_qnt"), recursive = TRUE)
  arrow::write_dataset(QNT[, !"agegrp20"],
                       file.path(out, "xps", "xps5_qnt"), format = "parquet")
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(
    logs = FALSE, exposures_for_output = character(0)
  )))
  env$private <- list()
  f <- Simulation$private_methods$export_xps_tables
  environment(f) <- env
  f(prbl = PRBL, output_dir = out, tables_dir = tab,
    strata_ons = strata_ons, strata_esp = esp_strata)
  tab
}

tab <- run_xps_tab(list("year"))
fn <- "exposure quantiles by year (not standardised).csv"
expect_true(fn %in% basename(list.files(tab)),
            info = "the quantile family is published as CSV")
q <- fread(file.path(tab, fn))

# ---- the year must NOT be promoted a second time -----------------------
expect_equal(sort(unique(q$year)), c(2019L, 2020L),
             info = "xps datasets already carry full years; +2000 gave 4019")

# ---- the 'All' marginal filter must actually fire ----------------------
# A `by year` table is the year-only marginal, so it must contain ONLY the
# rows groupingsets emitted with every other stratum as "All" -- not those
# averaged together with the agegrp20 x sex x qimd cells.
marg <- QNT[agegrp20 == "All" & sex == "All" & qimd == "All",
            .(v = mean(value)), keyby = .(year, scenario, xps, prob)]
setkey(marg, year, scenario, xps, prob)
gq0 <- q[order(year, scenario, exposure, prob)]
setkeyv(gq0, c("year", "scenario", "exposure", "prob"))
expect_equal(gq0[["xps_qnt_50.0%"]], marg$v, tolerance = 1e-10,
             info = "the year table is the marginal, not a mean over the cells")
# ...and the mixed-in version must give a DIFFERENT answer, or the assertion
# above would pass on unfiltered code.
mixed <- QNT[, .(v = mean(value)), keyby = .(year, scenario, xps, prob)]
setkey(mixed, year, scenario, xps, prob)
expect_false(isTRUE(all.equal(gq0[["xps_qnt_50.0%"]], mixed$v)),
             info = "filtered and unfiltered genuinely differ")

# ---- a stratification the data cannot support gets no file -------------
tab_ag <- run_xps_tab(list("year", c("year", "agegrp20")))
expect_false("exposure quantiles by year-agegroup (not standardised).csv" %in%
               basename(list.files(tab_ag)),
             info = "no file for a stratification whose cells are all 'All'")

# ---- the ESP label must name what was standardised over ---------------
# An empty `( standardised)` was the third symptom of the same root cause.
expect_false(any(grepl("\\( standardised\\)", list.files(tab))),
             info = "no file named '( standardised)' with nothing named")

# The ESP tables are age-standardised by construction, so `age` must appear in
# their label even though the dataset has no age column to be absent from the
# strata.
esp <- run_xps_tab(list("year"), esp_strata = list("year", c("year", "sex", "qimd")))
esp_files <- grep("^exposure quantiles", list.files(esp), value = TRUE)
expect_true(any(grepl("\\(age-sex-qimd standardised\\)", esp_files)),
            info = "a year-only ESP table is standardised over age, sex and qimd")
expect_true(any(grepl("\\(age standardised\\)", esp_files)),
            info = "a fully-stratified ESP table is still age-standardised")
expect_false(any(grepl("\\( standardised\\)", esp_files)),
             info = "and never labelled with nothing")

# `age` is stated, not derived from the strata. The xps5 dataset has no age
# axis, so age can never be one of the strata -- and a stray `agegrp20` in
# strata$xps_esp must NOT be able to remove it from the label, which is what a
# setdiff against the strata allowed.
# The ESP dataset has no `agegrp20` to stratify by, so that entry is skipped
# with a warning rather than raising `object 'agegrp20' not found` from inside
# a forked worker after the whole dataset has been read.
expect_warning(
  esp2 <- run_xps_tab(list("year"),
                      esp_strata = list("year", c("year", "agegrp20"))),
  pattern = "no column")
esp2_files <- grep("^exposure quantiles", list.files(esp2), value = TRUE)
expect_true(length(esp2_files) > 0L,
            info = "the stratifications that ARE supported are still written")
expect_false(any(grepl("agegroup", esp2_files)),
             info = "and the unsupported one produced no file")
expect_true(all(grepl("not standardised|\\(age", esp2_files)),
            info = "every ESP label still begins with age")
expect_false(any(grepl("\\(sex-qimd standardised\\)", esp2_files)),
             info = "a stray agegrp20 stratum cannot strip age from the label")
expect_true(all(c("year", "scenario", "prob", "exposure") %in% names(q)))
expect_true("xps_qnt_50.0%" %in% names(q),
            info = "value columns use the xps_qnt_ prefix, not xps_mean_")
expect_false(any(grepl("xps_mean_", names(q))),
             info = "and are not confusable with the means family")
expect_equal(sort(unique(q$prob)), PROBS,
             info = "`prob` is the exposure quantile, kept as rows")
expect_equal(sort(unique(q$exposure)), c("bmi_curr_xps", "sbp_curr_xps"))

# The published median is the MC median of the per-iteration exposure quantile.
bh <- QNT[agegrp20 == "All" & sex == "All" & qimd == "All",
          .(med = median(value)), keyby = .(year, scenario, xps, prob)]
setkey(bh, year, scenario, xps, prob)
gq <- q[order(year, scenario, exposure, prob)]
setkeyv(gq, c("year", "scenario", "exposure", "prob"))
expect_equal(gq[["xps_qnt_50.0%"]], bh$med, tolerance = 1e-10,
             info = "MC median of the per-iteration quantile, per exposure")

# Quantiles must be monotone in `prob` -- the single most obvious way for a
# mis-shaped long format to be wrong.
expect_true(gq[, all(diff(`xps_qnt_50.0%`) > 0),
               by = .(year, scenario, exposure)][, all(V1)],
            info = "the published quantiles increase with prob")

# A run without the flag wrote no quantile directory: that must be silent, not
# an error, since it is the default.
out2 <- tempfile("noq_"); dir.create(out2, recursive = TRUE)
tab2 <- tempfile("tab2_"); dir.create(tab2, recursive = TRUE)
env <- new.env(parent = asNamespace("IMPACTncdEngland"))
env$self <- list(design = list(sim_prm = list(
  logs = FALSE, exposures_for_output = character(0)
)))
env$private <- list()
f <- Simulation$private_methods$export_xps_tables
environment(f) <- env
expect_silent(f(prbl = PRBL, output_dir = out2, tables_dir = tab2,
                strata_ons = list("year"), strata_esp = list("year")))
expect_equal(length(list.files(tab2)), 0L)

# ===========================================================================
# 5. export_xps() itself, end to end
# ===========================================================================
# Sections 2-3 pin the wiring by inspection; this runs the real function, so
# the new groupingsets call and the `_pop_` companions are actually executed.
mk_sp <- function(n_never = 200L) {
  set.seed(7)
  n <- 1200L
  pop <- data.table(
    pid = seq_len(n),
    year = rep(19:20, each = n / 2L),
    age = sample(35:80, n, TRUE),
    sex = factor(sample(c("men", "women"), n, TRUE)),
    qimd = factor(sample(c("1 most deprived", "5 least deprived"), n, TRUE)),
    all_cause_mrtl = 0L,
    wt = runif(n, 0.8, 1.2),
    wt_esp = runif(n, 0.8, 1.2),
    bmi_curr_xps = rnorm(n, 27, 4),
    # a 0/1 indicator, like ets / statin_px / smok_never / smok_active: its
    # quantiles can only be 0 or 1, so they must be skipped
    ets_curr_xps = as.double(sample(0:1, n, TRUE)),
    # and a constant, which is equally uninformative
    const_curr_xps = 1
  )
  # never-smokers get zero smoking history, ever-smokers a real one
  pop[, smok_status_curr_xps := factor(
    c(rep("1", n_never), rep("4", n - n_never)), levels = c("1", "2", "3", "4"))]
  pop[, smok_cig_curr_xps := fifelse(smok_status_curr_xps == "1", 0,
                                     round(runif(n, 5, 30)))]
  pop[, smok_dur_curr_xps := fifelse(smok_status_curr_xps == "1", 0,
                                     round(runif(n, 5, 40)))]
  pop[, smok_packyrs_curr_xps := round(smok_cig_curr_xps * smok_dur_curr_xps / 20)]
  pop[, smok_quit_yrs_curr_xps := 0]
  list(pop = pop, mc = 1L, mc_aggr = 1L)
}

run_export_xps <- function(quantiles, n_never = 200L,
                           q_strata = list("year", c("year", "qimd"))) {
  out <- tempfile("xpsrun_"); dir.create(out, recursive = TRUE)
  sp <- mk_sp(n_never)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(
    logs = FALSE, ageL = 30L, ageH = 99L, init_year = 19L,
    exposures_for_output = character(0),
    # As in production: agegrp -> agegrp20 and dimd -> qimd
    # (xps_strata_from_output()).
    strata_for_output = c("year", "agegrp", "sex", "dimd"),
    export_xps_quantiles = quantiles,
    xps_quantile_probs = PROBS,
    xps_quantile_strata = q_strata
  )))
  env$private <- list(output_dir = function(x = "") file.path(out, x))
  f <- Simulation$private_methods$export_xps
  environment(f) <- env
  f(sp, "sc0")
  list(dir = out, pop = sp$pop)
}

# --- default: means only, no quantile directory ---------------------------
r_off <- run_export_xps(quantiles = FALSE)
expect_true(dir.exists(file.path(r_off$dir, "xps", "xps20")),
            info = "the means are written as before")
expect_false(dir.exists(file.path(r_off$dir, "xps", "xps20_qnt")),
             info = "and nothing extra when the flag is off")

# The `_pop_` companions are export-only and must not survive into sp$pop
# (which goes on to become the lifecourse parquet).
expect_false(any(grepl("_pop_curr_xps", names(r_off$pop))),
             info = "companions are cleaned up with agegrp20/smok_never")
expect_false("agegrp20" %in% names(r_off$pop))

# --- flag on: both quantile datasets appear -------------------------------
r_on <- run_export_xps(quantiles = TRUE)
expect_true(dir.exists(file.path(r_on$dir, "xps", "xps20_qnt")))
expect_true(dir.exists(file.path(r_on$dir, "xps", "xps5_qnt")))

qd <- CKutils::read_parquet_dt(file.path(r_on$dir, "xps", "xps20_qnt"))
expect_true(all(c("xps", "prob", "value", "year", "mc", "scenario") %in% names(qd)))
expect_equal(sort(unique(qd$prob)), PROBS)
expect_true("All" %in% qd$sex,
            info = "the groupingsets 'All' marginals are labelled as in the means")
expect_true(is.numeric(qd$value),
            info = "`value` stays numeric -- the All-relabel must skip it")

# ---- only the configured stratifications are computed -------------------
# This is the speed lever: cost is ~linear in the number of stratifications,
# so the quantiles deliberately use fewer than the means. With the default
# `list("year", c("year", "qimd"))`, `sex` and `agegrp20` are never grouped by.
# as.character(): a strata column round-trips through parquet as a
# dictionary, so it comes back a factor -- with "All" added as a level by
# data.table's set(), exactly as in the means datasets.
expect_equal(sort(unique(as.character(qd$sex))), "All",
             info = "sex is not a quantile stratification by default")
expect_equal(sort(unique(as.character(qd$agegrp20))), "All",
             info = "nor is agegrp20")
expect_true(length(unique(qd$qimd)) > 1L,
            info = "but qimd is, so the deprivation shift is visible")

# Asking for more gets more...
qd_more <- CKutils::read_parquet_dt(file.path(
  run_export_xps(TRUE, q_strata = list("year", c("year", "sex")))$dir,
  "xps", "xps20_qnt"))
expect_true(length(unique(qd_more$sex)) > 1L,
            info = "xps_quantile_strata is honoured")
expect_equal(sort(unique(as.character(qd_more$qimd))), "All",
             info = "...and only what was asked for is computed")

# ...and the ESP population drops agegrp20, since it standardises over age, so
# an agegrp20 entry must not survive into xps5_qnt as a phantom grouping.
q5 <- CKutils::read_parquet_dt(file.path(
  run_export_xps(TRUE, q_strata = list(c("year", "agegrp20")))$dir,
  "xps", "xps5_qnt"))
expect_false("agegrp20" %in% names(q5),
             info = "xps5 has no agegrp20 to group by")
expect_true(nrow(q5) > 0L,
            info = "the entry degrades to year-only rather than vanishing")

# ---- indicator and constant columns are skipped -------------------------
# Their quantiles could only ever return 0/1 (or the constant), so they are
# noise; the mean already carries everything they have to say.
expect_true("bmi_curr_xps" %in% qd$xps)
expect_false("ets_curr_xps" %in% qd$xps,
             info = "a 0/1 indicator gets no quantiles")
expect_false("const_curr_xps" %in% qd$xps,
             info = "nor does a constant")
# ...but they keep their MEANS, which are meaningful (a prevalence).
m20 <- CKutils::read_parquet_dt(file.path(r_on$dir, "xps", "xps20"))
expect_true("ets_curr_xps" %in% names(m20),
            info = "indicators are still summarised by their mean")

# The population companions reach the quantile output, and their medians differ
# from the ever-smoker ones in the expected direction.
expect_true("smok_packyrs_pop_curr_xps" %in% qd$xps,
            info = "the companion is picked up by the same _curr_xps grep")
med_all <- qd[sex == "All" & qimd == "All" & agegrp20 == "All" & prob == 0.5]
ever <- med_all[xps == "smok_packyrs_curr_xps", value]
popm <- med_all[xps == "smok_packyrs_pop_curr_xps", value]
expect_true(length(ever) > 0L && length(popm) > 0L)
expect_true(all(popm <= ever),
            info = "including never-smokers as zeros cannot raise the median")

# And the real quantiles are monotone in prob.
mono <- qd[sex == "All" & qimd == "All" & agegrp20 == "All" &
             xps == "bmi_curr_xps"][order(year, prob)]
expect_true(mono[, all(diff(value) > 0), by = year][, all(V1)],
            info = "BMI quantiles increase with prob in the real output")

# --- the ever-smoker/population divergence, through the real code ---------
# Same smokers, many more never-smokers: the conditional mean must be ~flat
# while the population mean falls.
m_base <- CKutils::read_parquet_dt(
  file.path(run_export_xps(TRUE, n_never = 100L)$dir, "xps", "xps20"))
m_scn <- CKutils::read_parquet_dt(
  file.path(run_export_xps(TRUE, n_never = 800L)$dir, "xps", "xps20"))
pick <- function(d, col) d[sex == "All" & qimd == "All" & agegrp20 == "All"][[col]]
expect_true(
  mean(pick(m_scn, "smok_packyrs_pop_curr_xps")) <
    mean(pick(m_base, "smok_packyrs_pop_curr_xps")),
  info = "population pack-years fall when fewer people smoke")
# na.rm here is the point, not a convenience: the ever-smoker mean is
# UNDEFINED (NaN) in a stratum with no ever-smokers, whereas the population
# mean is defined everywhere. That asymmetry is another reason not to read the
# conditional column as a population figure.
expect_true(anyNA(pick(m_scn, "smok_packyrs_curr_xps")),
            info = "no ever-smokers in a stratum -> the conditional mean is NaN")
expect_false(anyNA(pick(m_scn, "smok_packyrs_pop_curr_xps")),
             info = "...while the population mean is defined for every stratum")
expect_true(
  mean(pick(m_scn, "smok_packyrs_curr_xps"), na.rm = TRUE) >
    mean(pick(m_base, "smok_packyrs_pop_curr_xps"), na.rm = TRUE),
  info = "the ever-smoker mean sits above the population mean, as it must")


# ===========================================================================
# 6. `xps_quantile_strata` is validated in the PARENT
# ===========================================================================
# Same slot and same reason as export_tables()'s strata checks: a typo must
# fail at Design construction, not from inside a forked worker after a
# scenario has already run. The vocabulary is the xps one -- `agegrp20` and
# `qimd`, NOT `agegrp`/`dimd` -- which is exactly the slip worth catching.
# Built from a clean working directory so Design stays in skeleton mode (no
# data load), as test_Design_class.R does.
clean_wd <- file.path(tempdir(), "test_xpsq_clean_wd")
dir.create(clean_wd, showWarnings = FALSE, recursive = TRUE)
old_wd <- setwd(clean_wd)
on.exit(setwd(old_wd), add = TRUE)

mk_dl <- function(overrides = list()) {
  base <- list(
    locality = "England", clusternumber = 1, clusternumber_export = 1,
    logs = FALSE,
    cols_for_output = c("year", "age", "agegrp", "sex", "dimd"),
    strata_for_output = c("year", "agegrp", "sex", "dimd"),
    exposures_for_output = list("bmi"),
    exposure_definitions = list(list(
      name = "bmi", file_name = "bmi", var_name = "bmi",
      rank_var = "rank_bmi", distribution = "continuous")),
    n = 10, init_year_long = 2010, sim_horizon_max = 2030,
    ageL = 35, ageH = 80,
    diseases = list(chd = list(name = "chd", meta = list(
      incidence = list(influenced_by_disease_name = character(0))))),
    maxlag = 5, smoking_relapse_limit = 3, stochastic = TRUE, kismet = 42,
    jumpiness = 0.2, statin_adherence = 0.85, bpmed_adherence = 0.9,
    export_xps = TRUE, simsmok_calibration = FALSE,
    output_dir = tempdir(), synthpop_dir = tempdir(),
    validation = FALSE, iteration_n_max = 2, num_chunks = 1)
  modifyList(base, overrides)
}
qstrata <- function(overrides = list()) {
  suppressMessages(Design$new(mk_dl(overrides))$sim_prm$xps_quantile_strata)
}

# ---- the default: the cheap two, not the means' eight -------------------
expect_equal(qstrata(), list("year", c("year", "qimd")),
             info = "default is year + year-qimd (measures at parity with the means)")
expect_equal(suppressMessages(Design$new(mk_dl()))$sim_prm$export_xps_quantiles,
             FALSE, info = "and it is off unless asked for")
expect_equal(suppressMessages(Design$new(mk_dl()))$sim_prm$xps_quantile_probs,
             c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))

# ---- a DEFAULT must never be able to fail -------------------------------
# `qimd` only exists in the exposure output if `dimd` is in strata_for_output.
# A design without it must degrade, not error on a value the user never wrote.
expect_silent(dg <- qstrata(list(
  cols_for_output = c("year", "age", "sex"),
  strata_for_output = c("year", "age", "sex"))))
expect_equal(dg, list("year"),
             info = "no dimd -> the default degrades to year alone")

# ---- a value the USER wrote is validated strictly -----------------------
expect_equal(qstrata(list(xps_quantile_strata = list("year", c("year", "sex")))),
             list("year", c("year", "sex")),
             info = "an explicit list is honoured verbatim")
# `year` defines the rows, so every entry needs it.
expect_error(qstrata(list(xps_quantile_strata = list(c("sex", "qimd")))),
             pattern = "must include")
expect_error(qstrata(list(xps_quantile_strata = list(character(0)))),
             pattern = "must include")
# The xps vocabulary is agegrp20/qimd, NOT agegrp/dimd -- exactly the slip
# worth catching in the parent rather than inside a fork.
e_dimd <- tryCatch(qstrata(list(xps_quantile_strata = list(c("year", "dimd")))),
                   error = conditionMessage)
expect_true(grepl("cannot supply", e_dimd))
expect_true(grepl("agegrp20", e_dimd),
            info = "and the message names the vocabulary that does work")
expect_error(qstrata(list(xps_quantile_strata = list(c("year", "agegrp")))),
             pattern = "cannot supply")
expect_error(qstrata(list(xps_quantile_strata = list(c("year", "bananas")))),
             pattern = "cannot supply")

# Duplicates collapse rather than doubling the work.
expect_equal(qstrata(list(xps_quantile_strata = list("year", "year",
                                                     c("year", "qimd")))),
             list("year", c("year", "qimd")))

# probs are validated too.
expect_error(suppressMessages(Design$new(mk_dl(list(xps_quantile_probs = c(0, 0.5))))),
             pattern = "strictly between 0 and 1")
expect_error(suppressMessages(Design$new(mk_dl(list(xps_quantile_probs = numeric(0))))),
             pattern = "strictly between 0 and 1")

setwd(old_wd)
