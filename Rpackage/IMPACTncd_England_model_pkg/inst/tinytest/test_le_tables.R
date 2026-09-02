# The life-expectancy / healthy-life-expectancy table family.
#
# `export_summaries()` has always written `le`, `le60`, `hle_1st_cond` and
# `hle_cmsmm1.5`, but no table family consumed them, so a run produced the
# summaries and no CSV -- "life expectancy tables are missing" was a missing
# feature, not a fault.
#
# THE LOAD-BEARING CORRECTNESS PROPERTY
# -------------------------------------
# These values are NOT counts and NOT summed-numerator/summed-denominator
# rates. Each summary is built by
#     SELECT <strata>, SUM(w) AS popsize, SUM(age * w) / SUM(w) AS LE
# so the value is a weight-averaged age at an event and `popsize` is the
# weighted count of those events -- the very denominator the average was formed
# with. Collapsing over strata is therefore exact, but ONLY as a
# `weighted.mean(value, popsize)`; a plain mean or a sum would be wrong. Every
# "it worked" assertion below is that identity computed by hand from the
# fixture, plus its mirror image: the same assertion must FAIL for an
# unweighted mean, or it would not detect the mistake it exists to catch.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

PRBL <- c(0.5, 0.025, 0.975, 0.1, 0.9)   # the production quantiles
MED  <- "le_mean_50.0%"

DIMD_LV <- c("1 most deprived", as.character(2:9), "10 least deprived")
MEASURES <- c("le", "le60", "hle_1st_cond", "hle_cmsmm1.5")

# A fixture with a real gradient in BOTH the value and the weight: if the two
# were uncorrelated, weighted and unweighted means would coincide and section 3
# could not tell them apart.
mk_le <- function(measure, offset) {
  d <- CJ(mc = 1:3, scenario = c("sc0", "sc1"), year = 19:20,
          sex = c("men", "women"), dimd = DIMD_LV, sorted = FALSE)
  i <- match(d$dimd, DIMD_LV)
  d[, popsize := as.numeric(2000L - 100L * i + 10L * mc)]   # falls with decile
  d[, value := 70 + 0.8 * i + 0.5 * (sex == "women") +
       0.1 * (year - 19L) + 0.05 * mc + offset]             # rises with decile
  setnames(d, "value", if (startsWith(measure, "le")) "LE" else "HLE")
  d[]
}

SRC <- list(
  le            = mk_le("le", 0),
  le60          = mk_le("le60", 2),
  hle_1st_cond  = mk_le("hle_1st_cond", -20),
  `hle_cmsmm1.5` = mk_le("hle_cmsmm1.5", -5)
)

run_le <- function(strata_ons = list("year"),
                   strata_esp = list("year"),
                   src = SRC,
                   esp = FALSE) {
  dir <- tempfile("le_"); dir.create(dir, recursive = TRUE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(logs = FALSE)))
  env$private <- list(
    read_summary_dataset = function(summary_type, standardization) {
      if (!esp && identical(standardization, "esp")) return(NULL)
      s <- src[[summary_type]]
      if (is.null(s)) NULL else copy(s)
    }
  )
  f <- Simulation$private_methods$export_le_tables
  environment(f) <- env
  f(prbl = PRBL, summaries_dir = dir, tables_dir = dir,
    strata_ons = strata_ons, strata_esp = strata_esp)
  dir
}

# ===========================================================================
# 1. Files, names and shape
# ===========================================================================
d1 <- run_le(strata_ons = list("year", c("year", "sex")))
expect_equal(sort(basename(list.files(d1))),
             sort(c("life expectancy by year (not standardised).csv",
                    "life expectancy by year-sex (not standardised).csv")),
             info = "one file per stratification, none per measure")

le1 <- fread(file.path(d1, "life expectancy by year (not standardised).csv"))
expect_true(all(c("year", "scenario", "type", "n_events", MED) %in% names(le1)))
expect_equal(sort(unique(le1$type)), sort(MEASURES),
             info = "all four measures share one file, tagged by `type`")
expect_equal(sort(unique(le1$scenario)), c("sc0", "sc1"))
expect_equal(sort(unique(le1$year)), c(2019L, 2020L),
             info = "short years are promoted to full years")
expect_equal(nrow(le1), 2L * 2L * length(MEASURES),
             info = "year x scenario x type, with mc collapsed into quantiles")

# The ESP variant names what the weights standardised for.
d2 <- run_le(strata_esp = list("year", c("year", "sex")), esp = TRUE)
expect_true(all(c("life expectancy by year (age-sex-dimd standardised).csv",
                  "life expectancy by year-sex (age-dimd standardised).csv")
                %in% basename(list.files(d2))),
            info = "a stratum is not reported as standardised over itself")

# ===========================================================================
# 2. Correctness: the collapse is the popsize-weighted mean, exactly
# ===========================================================================
by_hand <- SRC$le[, .(v = weighted.mean(LE, popsize)),
                  keyby = .(mc, year, scenario)][
                    , .(med = median(v)), keyby = .(year, scenario)]
got <- le1[type == "le"][order(scenario, year)]
by_hand <- by_hand[order(scenario, year)]
expect_equal(got[[MED]], by_hand$med, tolerance = 1e-10,
             info = "collapsed LE == popsize-weighted mean over sex x dimd")

# Every measure, not just the first.
for (m in MEASURES) {
  vcol <- if (startsWith(m, "le")) "LE" else "HLE"
  bh <- SRC[[m]][, .(v = weighted.mean(get(vcol), popsize)),
                 keyby = .(mc, year, scenario)][
                   , .(med = median(v)), keyby = .(year, scenario)][
                     order(scenario, year)]
  expect_equal(le1[type == m][order(scenario, year)][[MED]], bh$med,
               tolerance = 1e-10, info = paste("weighted collapse exact for", m))
}

# ===========================================================================
# 3. Mirror image: an UNWEIGHTED mean must give a different answer
# ===========================================================================
# Without this, section 2 would pass on code that had dropped the weights in a
# fixture where they happened not to matter.
unwtd <- SRC$le[, .(v = mean(LE)), keyby = .(mc, year, scenario)][
  , .(med = median(v)), keyby = .(year, scenario)][order(scenario, year)]
expect_false(isTRUE(all.equal(got[[MED]], unwtd$med)),
             info = "the fixture correlates value with weight, so it can tell")

# ===========================================================================
# 4. n_events: the weighted event count behind each figure
# ===========================================================================
n_hand <- SRC$le[, .(n = sum(popsize)), keyby = .(mc, year, scenario)][
  , .(n = median(n)), keyby = .(year, scenario)][order(scenario, year)]
expect_equal(got$n_events, n_hand$n, tolerance = 1e-10,
             info = "n_events is the median weighted event count over mc")

# ===========================================================================
# 5. qimd strata on a dimd-only fixture, and their exactness
# ===========================================================================
dq <- run_le(strata_ons = list(c("year", "qimd")))
leq <- fread(file.path(dq, "life expectancy by year-qimd (not standardised).csv"))
expect_equal(length(unique(leq$qimd)), 5L,
             info = "qimd is derived from dimd without being in the summaries")

qmap <- data.table(dimd = factor(DIMD_LV, levels = DIMD_LV),
                   qimd = rep(c("1 most deprived", "2", "3", "4",
                                "5 least deprived"), each = 2L))
bh <- merge(SRC$le, qmap, by = "dimd")[
  , .(v = weighted.mean(LE, popsize)), keyby = .(mc, year, qimd, scenario)][
    , .(med = median(v)), keyby = .(year, qimd, scenario)]
gq <- leq[type == "le"]
setkey(bh, year, qimd, scenario); setkey(gq, year, qimd, scenario)
expect_equal(gq[[MED]], bh$med, tolerance = 1e-10,
             info = "quintile collapse == weighted mean over the decile pair")

# ===========================================================================
# 6. A summary that was never written is skipped, not fatal
# ===========================================================================
# `le60` is not produced when the design's age range does not span 60.
partial <- SRC[setdiff(names(SRC), "le60")]
expect_silent(dp <- run_le(src = partial))
lep <- fread(file.path(dp, "life expectancy by year (not standardised).csv"))
expect_equal(sort(unique(lep$type)), sort(setdiff(MEASURES, "le60")),
             info = "the other three measures are still published")
# ...and nothing at all is still not an error.
expect_silent(de <- run_le(src = list()))
expect_equal(length(list.files(de)), 0L)

# ===========================================================================
# 7. Wiring: the family is dispatched, and can be switched off
# ===========================================================================
et <- Simulation$public_methods$export_tables
expect_true("life_expectancy" %in% names(formals(et)))
expect_true(eval(formals(et)$life_expectancy), info = "on by default")
et_txt <- paste(deparse(body(et)), collapse = "\n")
expect_true(grepl('type = "le"', et_txt, fixed = TRUE),
            info = "the task is added to the task list")
expect_true(grepl("strata_cfg$le_ons", et_txt, fixed = TRUE),
            info = "...with its own strata element")
hlpr_txt <- paste(deparse(body(Simulation$private_methods$export_tables_hlpr)),
                  collapse = "\n")
expect_true(grepl("export_le_tables", hlpr_txt, fixed = TRUE),
            info = "and the dispatcher knows how to run it")
