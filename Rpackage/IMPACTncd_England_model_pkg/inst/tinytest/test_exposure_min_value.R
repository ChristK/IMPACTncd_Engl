## Exposure bounds are INCLUSIVE (2026-11-18): draws are truncated to
## [min_value, max_value] with minq = P(Y < min_value) and maxq = P(Y <= max_value),
## by the distribution's support (.exposure_cdf_bounds(), Exposure `family_type`).
## Until then minq was P(Y <= min_value): a point mass AT min_value was never
## drawn, so `min_value: 0` on a count gave no non-drinkers and no one eating no
## fruit or veg, and the rescaled rank range moved every simulant up the
## distribution. Ported from IMPACTncd_Engl2026 (fixed there 2026-10-11).
##
##   1  .exposure_family_type(): qfun name -> gamlss.dist family type (every
##      stock qfun); a declared value wins; Exposure$new() stops without a type;
##      Design$load_exposures() passes a YAML family_type on
##   2  the bounds by family: discrete (min 0 -> 0 WITHOUT calling pfun; min 1 ->
##      F(0), the old minq; fractional min / max handed to pfun as 2 / 80),
##      continuous (F(min), the old formula), mixed (ZAGA: 0 at min 0); a stale
##      stored minq is reported once per object and recalculated every time; a
##      pfun that is not a CDF (gamlss.dist 6.1.11 pZABNB) is refused, after
##      rows a non-vectorised pfun flipped (6.1.11 pGG) are redone one by one
##   3  the draw, for each family type: a point mass at min_value is drawn with
##      its probability (discrete ZINBI, mixed ZAGA); continuous draws stay in
##      [min, max]; CONTROL: the old minq draws no zero
##   4  the PARF path uses the SAME bounds as the synthetic population:
##      .refresh_xps_bounds() == add_qmin_qmax(), with no random numbers and no
##      message; every bounded draw in Disease$gen_sp_forPARF() refreshes its
##      table between reading and joining it; the PARF population also gets the
##      synthpop's age-12 rescale, quit-years guard and cigarettes 0 -> 1
##   5  ordinal threshold detection ignores minq / maxq (active_days was 1..8)
##   6  the shipped designs: smoking counts min_value + offset >= 1 (option B,
##      the previous bounds), fruit / veg / alcohol 0, every bounded exposure
##      typed, smok_cig_ex on CKutils' ZABNB, and Design refuses a smoking count
##      that could draw 0 or a pZABNB pfun. Set IMPACTNCD_PROJECT_ROOT to also check every
##      tracked design YAML of a checkout.
##   7  the caches key on the PARF exposures' definitions (canonicalised) and the
##      bounds rule; a superseded PARF directory is removed -- only this
##      disease's -- while every .qs is kept
##   8  the smoking engine: an ex-smoker's quit years >= 1 (synthpop guard after
##      the age-12 rescale; simsmok() / simsmok_sc() stop instead of reading
##      pr_relapse before the matrix)
library(tinytest)
library(IMPACTncdEngland)
library(data.table)
suppressPackageStartupMessages(library(gamlss.dist)) # pZAGA / qZAGA for the mixed case

ns <- asNamespace("IMPACTncdEngland")
ft <- get(".exposure_family_type", ns)
cdf_bounds <- get(".exposure_cdf_bounds", ns)
refresh <- get(".refresh_xps_bounds", ns)
even_idx <- get(".even_sample_idx", ns)
`%||%` <- function(x, y) if (is.null(x)) y else x
# (deparse() writes `x := v` as `:=`(x, v), and the patterns below are written that way)
ws <- function(f) gsub("[[:space:]]", "", paste(deparse(f, width.cutoff = 500L), collapse = ""))
# tinytest's expect_silent() checks warnings and errors only (a message() passes
# it), so messages are counted
n_msg <- function(expr) {
  n <- 0L
  withCallingHandlers(expr, message = function(m) { n <<- n + 1L; invokeRestart("muffleMessage") })
  n
}

## ---- 1. the family type --------------------------------------------------------
# every qfun of the stock design, and the testing design's extra ones
stock <- c(fqMN4 = "Discrete", fqDPO = "Discrete", fqNBI = "Discrete", fqZINBI = "Discrete",
           fqZABNB = "Discrete", qZABNB = "Discrete", fqZISICHEL = "Discrete", fqDEL = "Discrete",
           fqBCPEo = "Continuous", fqBCT = "Continuous", qGB1 = "Continuous",
           qGA = "Continuous", qBCTo = "Continuous", qGG = "Continuous")
for (q in names(stock)) {
  expect_identical(ft(q), stock[[q]], info = paste0("1: ", q, " -> ", stock[[q]]))
}
expect_identical(ft("my_qBCPEo"), "Continuous", info = "1: a prefixed qfun")
expect_identical(ft("qZAGA"), "Mixed", info = "1: a zero-adjusted continuous family is Mixed")
expect_null(ft("tt_qnone"), info = "1: no gamlss family in the name -> NULL")
expect_null(ft("qNOSUCHFAMILY"), info = "1: a family-shaped name gamlss.dist lacks -> NULL")
expect_null(ft(NULL), info = "1: no qfun -> NULL")
expect_identical(ft("tt_qnone", "discrete"), "Discrete", info = "1: a declared type wins, any case")
expect_identical(ft("fqZINBI", "MIXED"), "Mixed", info = "1: ...even over the inferred one")
expect_error(ft("fqZINBI", "count"), "family_type must be one of", info = "1: an unknown declared type stops")

## ---- fixtures: tiny parameter tables on disk (Exposure reads their column names)
td <- file.path(tempdir(), "tt_exposure_min_value"); unlink(td, recursive = TRUE)
dir.create(td, recursive = TRUE)
mk <- function(nm, tbl) {
  p <- file.path(td, nm)
  CKutils::write_parquet_dt(tbl, p, partition = "year")
  p
}
zinbi <- data.table(year = 1L, age = 1:8, mu = c(0.5, 1, 2, 5, 5, 10, 20, 40),
                    sigma = c(0.5, 1, 2, 1, 0.3, 1, 0.8, 1.5), nu = c(0.1, 0.6, 0.3, 0.3, 0.05, 0.2, 0.4, 0.5))
bcpe <- data.table(year = 1L, age = 1:4, mu = c(22, 27, 31, 40), sigma = c(0.1, 0.15, 0.2, 0.18),
                   nu = c(-1, -0.5, 0, 0.5), tau = c(1.5, 2, 2.5, 3))
zaga <- data.table(year = 1L, age = 1:4, mu = c(2, 10, 30, 60), sigma = c(0.5, 0.8, 1, 1.2),
                   nu = c(0.05, 0.3, 0.5, 0.8))
p_zinbi <- mk("zinbi", zinbi); p_bcpe <- mk("bcpe", bcpe); p_zaga <- mk("zaga", zaga)

# a ZINBI pfun that records the q it is handed and refuses a negative one (as
# CKutils' fast pfuns up to 0.1.35 do), so "min 0 never calls pfun below 0" is
# observable
.tt_mv <- new.env()
.tt_mv$q <- numeric(0)
assign("tt_mv_fpZINBI", function(q, mu, sigma, nu) {
  if (any(q < 0)) stop("tt_mv_fpZINBI: negative q")
  .tt_mv$q <- c(.tt_mv$q, unique(q))
  CKutils::fpZINBI(q, mu, sigma, nu)
}, envir = globalenv())
# a "CDF" that is 0 from q = 1 on, as gamlss.dist 6.1.11's pZABNB is
assign("tt_broken_p", function(q, mu, sigma, nu) ifelse(q < 1, 0.3, 0), envir = globalenv())

xps <- function(path, qfun, pfun = NULL, min, max, qparams, family_type = NULL) {
  Exposure$new(name = "zz", file_name = "zz", var_name = "zz", rank_var = "rank_zz",
               distribution = "continuous", qfun = qfun, pfun = pfun, min_value = min,
               max_value = max, qparams = qparams, family_type = family_type, file_path = path)
}
bounds <- function(x, tbl) {
  t <- copy(tbl)
  suppressMessages(x$add_qmin_qmax(t))
  t
}
pr <- c("mu", "sigma", "nu")

## ---- 1b. Exposure$new(), Design$load_exposures() and the type -----------------
expect_error(xps(p_zinbi, "tt_qnone", "tt_mv_fpZINBI", 0, 250, pr),
             "family_type", info = "1b: bounds on a qfun of unknown family stop, naming family_type")
x <- tryCatch(xps(p_zinbi, "tt_qnone", "tt_mv_fpZINBI", 0, 250, pr, family_type = "discrete"),
              error = function(e) NULL)
expect_false(is.null(x), info = "1b: ...and declaring family_type lets the Exposure be built")
expect_identical(x$family_type, "Discrete", info = "1b: ...as the declared type")
expect_silent(xps(p_zinbi, "tt_qnone", "tt_mv_fpZINBI", NULL, NULL, pr))
expect_identical(xps(p_zinbi, "fqZINBI", NULL, 0, 250, NULL)$pfun, "fpZINBI",
                 info = "1b: pfun derivation unchanged (fqX -> fpX)")
expect_identical(xps(p_zinbi, "qZABNB", NULL, 0, 80, NULL)$pfun, "pZABNB", info = "1b: qX -> pX")
expect_identical(xps(p_bcpe, "my_qBCPEo", NULL, 14, 70, NULL)$pfun, "my_pBCPEo", info = "1b: prefix_qX -> prefix_pX")

# Design$load_exposures() run on a stub self / private (no YAML, no inputs)
le_with <- function(defs) {
  le <- Design$public_methods$load_exposures
  stub_self <- new.env()
  stub_self$assert_working_dir <- function(...) invisible(TRUE)
  stub_self$sim_prm <- list(exposure_definitions = defs)
  environment(le) <- list2env(list(self = stub_self,
                                   private = list(resolve_exposure_path = function(x) x,
                                                  reorder_exposures = function(x) x)), parent = ns)
  suppressMessages(le())
  stub_self
}
def <- function(name, qfun, min, max, ...) {
  c(list(name = name, file_name = p_zinbi, var_name = name, rank_var = paste0("rank_", name),
         distribution = "continuous", qfun = qfun, min_value = min, max_value = max, qparams = pr), list(...))
}
got <- tryCatch(vapply(le_with(list(def("zz", "tt_qnone", 0, 250, pfun = "tt_mv_fpZINBI", family_type = "discrete"),
                                    def("yy", "fqZINBI", 0, 250, family_type = "mixed")))$exposures,
                       function(e) e$family_type %||% NA_character_, ""),
                error = function(e) conditionMessage(e))
expect_identical(got, c(zz = "Discrete", yy = "Mixed"),
                 info = "1b: Design$load_exposures() passes the YAML's family_type (declared wins over inferred)")

## ---- 2. the bounds, by family ----------------------------------------------------
F0 <- zinbi[, CKutils::fpZINBI(rep(0, .N), mu, sigma, nu)]
F2 <- zinbi[, CKutils::fpZINBI(rep(2, .N), mu, sigma, nu)]
F80 <- zinbi[, CKutils::fpZINBI(rep(80, .N), mu, sigma, nu)]
F250 <- zinbi[, CKutils::fpZINBI(rep(250, .N), mu, sigma, nu)]
expect_true(all(F0 > 0.05), info = "2: (non-vacuity) the fixture has a real mass at 0")

.tt_mv$q <- numeric(0)
t0 <- tryCatch(bounds(xps(p_zinbi, "fqZINBI", "tt_mv_fpZINBI", 0, 250, pr), zinbi), error = function(e) NULL)
expect_false(is.null(t0), info = "2: min 0 never hands pfun a q below 0 (the recording pfun would stop)")
expect_identical(t0$minq, rep(0, nrow(zinbi)), info = "2: discrete, min 0 -> minq = P(Y < 0) = 0")
expect_identical(t0$maxq, F250, info = "2: discrete, maxq = P(Y <= 250)")
expect_identical(.tt_mv$q, 250, info = "2: ...and pfun was called for the maximum only: never with a q below 0")
t1 <- bounds(xps(p_zinbi, "fqZINBI", "tt_mv_fpZINBI", 1, 250, pr), zinbi)
expect_identical(t1$minq, F0, info = "2: discrete, min 1 -> minq = F(0): exactly what min 0 gave before (the smoking counts, option B)")
.tt_mv$q <- numeric(0)
t25 <- bounds(xps(p_zinbi, "fqZINBI", "tt_mv_fpZINBI", 2.5, 80.7, pr), zinbi)
expect_identical(t25$minq, F2, info = "2: discrete, min 2.5 -> P(Y < 2.5) = F(2)")
expect_identical(t25$maxq, F80, info = "2: discrete, max 80.7 -> F(80)")
expect_identical(.tt_mv$q, c(2, 80),
                 info = "2: ...pfun is handed ceiling(2.5) - 1 = 2 and floor(80.7) = 80 (CKutils' pfuns floor q anyway; a custom one may not)")
expect_true(all(t0$minq < F0), info = "2: CONTROL: the old minq for min 0, P(Y <= 0), differs from the new one")

prb <- c("mu", "sigma", "nu", "tau")
tb <- bounds(xps(p_bcpe, "fqBCPEo", NULL, 14, 70, prb), bcpe)
expect_identical(tb$minq, bcpe[, CKutils::fpBCPEo(rep(14, .N), mu, sigma, nu, tau)],
                 info = "2: continuous -> F(min), the formula as it was (no continuous output moves)")
expect_identical(tb$maxq, bcpe[, CKutils::fpBCPEo(rep(70, .N), mu, sigma, nu, tau)])
expect_true(all(tb$minq > 0), info = "2: (non-vacuity) the continuous minq is a real truncation, not 0")

tz <- bounds(xps(p_zaga, "qZAGA", NULL, 0, 250, pr), zaga)
expect_identical(tz$minq, rep(0, nrow(zaga)), info = "2: mixed, min 0 -> 0: ZAGA's mass at 0 is kept")
expect_true(all(zaga[, pZAGA(rep(0, .N), mu, sigma, nu)] > 0.04), info = "2: (non-vacuity) ...which F(0) would have cut")
tz5 <- bounds(xps(p_zaga, "qZAGA", NULL, 0.5, 250, pr), zaga)
expect_identical(tz5$minq, zaga[, pZAGA(rep(0.5, .N), mu, sigma, nu)], info = "2: mixed, min above the mass -> F(min)")

expect_error(cdf_bounds(list(mu = 1), 1L, "fpZINBI", NULL, 0, 1, "zz"), "family_type",
             info = "2: the helper refuses a missing type")
expect_error(cdf_bounds(list(mu = 1), 1L, "fpZINBI", "count", 0, 1, "zz"), "family_type",
             info = "2: ...and an unknown one")
# a pfun that is not a CDF (maxq below minq, or 0) is refused, not drawn from: gamlss.dist 6.1.11's pZABNB made
# every smok_cig_ex draw 0
xb <- xps(p_zinbi, "tt_qnone", "tt_broken_p", 1, 80, pr, family_type = "discrete")
expect_error(xb$add_qmin_qmax(copy(zinbi)), "not a valid CDF", info = "2: a broken pfun stops (maxq < minq)")
xb0 <- xps(p_zinbi, "tt_qnone", "tt_broken_p", 0, 80, pr, family_type = "discrete")
expect_error(xb0$add_qmin_qmax(copy(zinbi)), "not a valid CDF", info = "2: ...also when minq is 0 (maxq 0)")
expect_error(cdf_bounds(as.list(zinbi[, .(mu, sigma, nu)]), nrow(zinbi), "pZABNB_like", "Discrete", 1, 80, "zz") ,
             "object 'pZABNB_like' not found", info = "2: (a missing pfun still errors as itself)")
m_broken <- tryCatch(xb$add_qmin_qmax(copy(zinbi)), error = conditionMessage)
expect_false(grepl("ZABNB", m_broken), info = "2: the message names the ZABNB remedy only for a ZABNB pfun")
# a pfun that is valid row by row but not vectorised over a sign parameter (gamlss.dist 6.1.11's pGG hands pGA a
# vector lower.tail, so rows whose nu sign differs from the first row's get 1 - F): those rows are redone one at a time
assign("tt_signvec_p", function(q, mu, sigma, nu) {
  out <- pnorm(q, mu, sigma)
  flip <- sign(nu) != sign(nu[1L])
  out[flip] <- 1 - out[flip]
  out
}, envir = globalenv())
sv <- data.table(mu = c(30, 40, 50, 60, 70), sigma = 5, nu = c(-1, -1, 1, -0.5, 2))
b_sv <- cdf_bounds(as.list(sv), nrow(sv), "tt_signvec_p", "Continuous", 0, 200, "zz")
expect_equal(b_sv$minq, pnorm(0, sv$mu, 5), tolerance = 0, info = "2: rows a non-vectorised pfun flipped are recomputed row by row (minq)")
expect_equal(b_sv$maxq, pnorm(200, sv$mu, 5), tolerance = 0, info = "2: ...and maxq; no error")
gg <- data.table(mu = c(30, 28, 35, 31, 26), sigma = c(0.5, 0.6, 0.4, 0.5, 0.7), nu = c(-0.4, -0.2, 0.6, -0.5, 1.1))
b_gg <- tryCatch(cdf_bounds(as.list(gg), nrow(gg), "pGG", "Continuous", 0, 2000, "hfss_food_weight"), error = function(e) NULL)
expect_false(is.null(b_gg), info = "2: real gamlss.dist pGG with mixed-sign nu (the testing design's hfss_food_weight) does not stop")
row_gg <- vapply(seq_len(nrow(gg)), function(i) pGG(2000, gg$mu[i], gg$sigma[i], gg$nu[i]), 0)
expect_equal(b_gg$maxq, row_gg, tolerance = 0, info = "2: ...and its maxq equals pGG row by row")

# a stale stored minq (the shipped count tables hold the old one) is recomputed and reported once per object
x0 <- xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr)
stale <- copy(zinbi)[, `:=`(minq = F0, maxq = F250)]
expect_identical(n_msg(x0$add_qmin_qmax(stale)), 1L, info = "2: a table holding the old minq is reported...")
expect_identical(stale$minq, rep(0, nrow(zinbi)), info = "2: ...and recalculated")
stale2 <- copy(zinbi)[, `:=`(minq = F0, maxq = F250)]
expect_identical(n_msg(x0$add_qmin_qmax(stale2)), 0L, info = "2: ...not reported again by the same object")
expect_identical(stale2$minq, rep(0, nrow(zinbi)), info = "2: ...although recalculated every time")
stale3 <- copy(zinbi)[, `:=`(minq = F0, maxq = F250)]
expect_identical(n_msg(xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr)$add_qmin_qmax(stale3)), 1L,
                 info = "2: a NEW object reports again (the flag is per object)")
fresh <- copy(zinbi)[, `:=`(minq = 0, maxq = F250)]
expect_identical(n_msg(xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr)$add_qmin_qmax(fresh)), 0L,
                 info = "2: a table already on the inclusive rule passes without a message")
expect_identical(fresh$minq, rep(0, nrow(zinbi)))

# the synthetic-population route still draws its check rows with sample.int(), exactly as before: the fix leaves the
# ambient base-R stream where it was (same consumption whether or not the bounds are stale)
big <- zinbi[rep(seq_len(.N), 2500L)][, `:=`(minq = F0[age], maxq = F250[age])] # 20,000 rows -> 2 check rows
set.seed(42L); suppressMessages(xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr)$add_qmin_qmax(copy(big)))
s_fix <- .Random.seed
set.seed(42L); invisible(sample.int(nrow(big), 2L))
expect_identical(s_fix, .Random.seed, info = "2: add_qmin_qmax() consumes the base-R stream as one sample.int(), as before")

## ---- 3. the draw, by family type -------------------------------------------------
N <- 100000L
rk <- (seq_len(N) - 0.5) / N
dsg <- list(sim_prm = list(logs = FALSE, keep_simulants_rn = FALSE))
# generate_continuous() on a population carrying the table's columns and the bounds
draw_direct <- function(x, params, minq_override = NULL) {
  pop <- data.table(rank_zz = rk)
  for (p in names(params)) set(pop, j = p, value = params[[p]])
  tb <- as.data.table(params)
  suppressMessages(x$add_qmin_qmax(tb))
  pop[, `:=`(minq = if (is.null(minq_override)) tb$minq else minq_override, maxq = tb$maxq)]
  x$.__enclos_env__$private$generate_continuous(pop)
  pop$zz
}
# generate(): the table looked up from the parquet partition (year 1, the given age's row)
draw_live <- function(x, age_) {
  pop <- data.table(year = 1L, age = as.integer(age_), rank_zz = rk)
  suppressMessages(x$generate(pop, dsg))
  pop$zz
}

# -- discrete (ZINBI, age 4 row: mu 5, sigma 1, nu 0.3)
par <- as.list(zinbi[age == 4L, .(mu, sigma, nu)])
p0 <- CKutils::fpZINBI(0, par$mu, par$sigma, par$nu)
pmax_ <- CKutils::fpZINBI(250, par$mu, par$sigma, par$nu)
y0d <- draw_direct(xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr), par)
y0l <- draw_live(xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr), 4L)
expect_true(abs(mean(y0d == 0) - p0 / pmax_) < 1e-3, info = "3 discrete: min_value 0 draws zeros at their own share")
expect_identical(y0l, y0d, info = "3 discrete: generate() draws what generate_continuous() draws")
y1 <- draw_live(xps(p_zinbi, "fqZINBI", NULL, 1, 250, pr), 4L)
expect_equal(min(y1), 1, info = "3 discrete: min_value 1 draws no zero (the smoking counts, option B)")
yold <- draw_direct(xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr), par, minq_override = p0)
expect_identical(sum(yold == 0), 0L, info = "3 discrete: CONTROL: the old minq, P(Y <= 0), draws no zero -- the check above can fail")
expect_true(mean(yold) > mean(y0d) / (1 - p0) * 0.99,
            info = "3 discrete: CONTROL: ...and inflates the mean by ~1/(1 - P0): every simulant moved up")

# -- mixed (ZAGA, age 2 row: mass nu = 0.3 at 0)
parz <- as.list(zaga[age == 2L, .(mu, sigma, nu)])
pz0 <- parz$nu
pzmax <- pZAGA(250, parz$mu, parz$sigma, parz$nu)
z0 <- draw_live(xps(p_zaga, "qZAGA", NULL, 0, 250, pr), 2L)
expect_true(abs(mean(z0 == 0) - pz0 / pzmax) < 1e-3, info = "3 mixed: the mass at 0 is drawn with its probability")
zold <- draw_direct(xps(p_zaga, "qZAGA", NULL, 0, 250, pr), parz, minq_override = pZAGA(0, parz$mu, parz$sigma, parz$nu))
expect_identical(sum(zold == 0), 0L, info = "3 mixed: CONTROL: the old minq draws no zero")

# -- continuous (BCPEo, age 2 row): no point mass; the draws stay in [min, max] and follow the truncated CDF
parb <- as.list(bcpe[age == 2L, .(mu, sigma, nu, tau)])
b <- draw_live(xps(p_bcpe, "fqBCPEo", NULL, 14, 70, prb), 2L)
Fb <- function(q) CKutils::fpBCPEo(q, parb$mu, parb$sigma, parb$nu, parb$tau)
expect_true(min(b) >= 14 - 1e-8 && max(b) <= 70 + 1e-8, info = "3 continuous: draws within [min, max]")
expect_true(abs(mean(b <= 25) - (Fb(25) - Fb(14)) / (Fb(70) - Fb(14))) < 1e-3,
            info = "3 continuous: ...with the truncated distribution's probabilities")

## ---- 4. the PARF path and the synthetic population use the same bounds ------------
x_sp <- xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr)
x_pf <- xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr)
d_pf <- list(exposures = list(zz = x_pf))
for (case in c("stale", "fresh")) {
  tb0 <- copy(big)
  if (case == "fresh") tb0[, minq := 0]
  t_sp <- copy(tb0); suppressMessages(x_sp$add_qmin_qmax(t_sp))
  t_pf <- copy(tb0); refresh(t_pf, d_pf, "zz")
  expect_identical(t_pf$minq, t_sp$minq, info = paste0("4 (", case, " table): PARF minq == synthpop minq"))
  expect_identical(t_pf$maxq, t_sp$maxq, info = paste0("4 (", case, " table): PARF maxq == synthpop maxq"))
}
expect_true(all(t_pf$minq == 0), info = "4: (non-vacuity) both are the inclusive bound")
expect_identical(n_msg(refresh(copy(big), list(exposures = list(zz = xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr))), "zz")), 0L,
                 info = "4: the PARF route recalculates a stale stock table silently (one message per forked task otherwise)")
x_q <- xps(p_zinbi, "fqZINBI", NULL, 0, 250, pr)
refresh(copy(big), list(exposures = list(zz = x_q)), "zz")
expect_identical(n_msg(x_q$add_qmin_qmax(copy(big))), 1L,
                 info = "4: a quiet PARF refresh does not use up the object's stale-table report (the synthpop path still reports)")
set.seed(7L); s0 <- .Random.seed
refresh(copy(big), d_pf, "zz")
expect_identical(.Random.seed, s0, info = "4: the PARF route consumes no random numbers (evenly spaced check rows)")
expect_identical(even_idx(20000L), c(1L, 20000L), info = "4: 0.01% of the rows, evenly spaced, first and last")
expect_identical(even_idx(0L), integer(0))
t_miss <- copy(big)
expect_warning(refresh(t_miss, d_pf, "not_defined"), "not in exposure_definitions",
               info = "4: an exposure the design lacks keeps its stored bounds, with a warning")
expect_identical(t_miss$minq, big$minq)
# only continuous exposures carry quantile bounds; an ordinal (binary, custom) table is left as it is
pa_nb <- data.table(year = 1L, age = 1:2, pa0 = c(0.2, 0.3), pa1 = c(0.5, 0.6), pa2 = c(0.8, 0.9))
p_pa_nb <- mk("pa_nb", pa_nb)
x_ord <- Exposure$new(name = "pa", file_name = "pa", var_name = "pa", rank_var = "rank_pa", distribution = "ordinal",
                      file_path = p_pa_nb)
pa_t <- copy(pa_nb)
refresh(pa_t, list(exposures = list(pa = x_ord)), "pa")
expect_identical(names(pa_t), names(pa_nb), info = "4: .refresh_xps_bounds() adds no bounds to an ordinal table")

# every bounded draw of the PARF generator refreshes its table between reading and joining it (13 draws read minq)
parf_fn <- Disease$private_methods$gen_sp_forPARF %||% Disease$public_methods$gen_sp_forPARF
src <- paste(deparse(parf_fn, width.cutoff = 500L), collapse = "\n")
parf_xps <- c("fruit", "veg", "smok_status", "smok_quit_yrs", "smok_dur_ex", "smok_dur_curr", "smok_cig_curr",
              "smok_cig_ex", "alcohol", "bmi", "sbp", "tchol", "hdl_to_tchol")
for (nm in parf_xps) {
  at_read <- regexpr(paste0("xps_table_path(\"", nm, "\")"), src, fixed = TRUE)
  rest <- substring(src, at_read)
  at_join <- regexpr("(lookup|absorb)_dt\\(ff, tbl", rest)
  expect_true(at_read > 0L && at_join > 0L &&
                grepl(paste0(".refresh_xps_bounds(tbl, design_, \"", nm, "\")"), substr(rest, 1L, at_join),
                      fixed = TRUE),
              info = paste0("4: gen_sp_forPARF refreshes ", nm, "'s bounds after reading its table, before the join"))
}
n_minq_draws <- lengths(regmatches(src, gregexpr("minq + ", src, fixed = TRUE)))
n_refresh <- lengths(regmatches(src, gregexpr(".refresh_xps_bounds(", src, fixed = TRUE)))
expect_identical(n_refresh, n_minq_draws, info = "4: as many bounds refreshes as minq draws (a new PARF draw needs one)")
expect_identical(n_minq_draws, length(parf_xps), info = "4: the 13 known draws")
# .PARF_EXPOSURES (what the PARF key hashes) is exactly the set of tables gen_sp_forPARF reads
read_names <- unique(regmatches(src, gregexpr("(?<=xps_table_path\\(\")[A-Za-z0-9_]+(?=\"\\))", src, perl = TRUE))[[1L]])
expect_identical(sort(read_names), sort(get(".PARF_EXPOSURES", ns)),
                 info = "4: .PARF_EXPOSURES lists exactly the exposures gen_sp_forPARF reads")
# the PARF population gets the synthpop's smoking post-processing (pre-existing gaps, fixed 2026-11-19)
s_parf <- ws(parf_fn)
expect_true(grepl("smok_dur_curr_xps+smok_quit_yrs_curr_xps>.age12", s_parf, fixed = TRUE),
            info = "4: PARF smoking histories are rescaled to start at age 12, as in the synthpop")
expect_true(grepl("smok_status_curr_xps%in%2:3&smok_quit_yrs_curr_xps<1L,`:=`(smok_quit_yrs_curr_xps,1L)", s_parf, fixed = TRUE),
            info = "4: ...and an ex-smoker's quit years are at least 1")
expect_true(grepl("smok_cig_curr_xps==0L&smok_status_curr_xps!=\"1\",`:=`(smok_cig_curr_xps,1L)", s_parf, fixed = TRUE),
            info = "4: ...and every ever-smoker smokes at least 1 cigarette a day, as in the synthpop")

## ---- 5. ordinal threshold detection ignores minq / maxq (active_days) -------------
pa <- data.table(year = 1L, age = 1:2, pa0 = c(0.2, 0.3), pa1 = c(0.5, 0.6), pa2 = c(0.8, 0.9), minq = 0, maxq = 1)
p_pa <- mk("pa", pa)
ord <- function(thr = NULL) {
  e <- Exposure$new(name = "pa", file_name = "pa", var_name = "pa", rank_var = "rank_pa", distribution = "ordinal",
                    thresholds = thr, offset = -1, file_path = p_pa)
  e$cast_to_int <- TRUE
  pop <- data.table(year = 1L, age = rep(1:2, each = 1000L), rank_pa = rep((seq_len(1000L) - 0.5) / 1000, 2L))
  e$generate(pop, dsg)
  list(e = e, y = pop$pa)
}
auto <- ord()
expect_identical(auto$e$thresholds, c("pa0", "pa1", "pa2"), info = "5: auto-detected thresholds exclude minq / maxq")
expect_identical(range(auto$y), c(0L, 3L), info = "5: offset -1 gives 0..3 (active_days 0..7)")
expect_identical(auto$y, ord(c("pa0", "pa1", "pa2"))$y, info = "5: ...the same as declaring the thresholds")
bug <- ord(c("maxq", "minq", "pa0", "pa1", "pa2"))$y
expect_identical(bug, auto$y + 1L, info = "5: CONTROL: with minq / maxq taken as thresholds every value is +1 (the bug)")

## ---- 6. the shipped designs --------------------------------------------------------
ymls <- c(default = system.file("config", "default_sim_design.yaml", package = "IMPACTncdEngland"))
root <- Sys.getenv("IMPACTNCD_PROJECT_ROOT", "")
if (nzchar(root) && dir.exists(root)) {
  more <- c(file.path(root, "inputs", "sim_design.yaml"), file.path(root, "testing", "sim_design_testing.yaml"),
            Sys.glob(file.path(root, "scenarios", "*", "sim_design_*.yaml")))
  more <- more[file.exists(more)]
  more <- more[vapply(more, function(f) !is.null(yaml::read_yaml(f)$exposure_definitions), NA)]
  names(more) <- sub(paste0("^", normalizePath(root), "/"), "", normalizePath(more))
  ymls <- c(ymls, more)
}
expect_true(file.exists(ymls[["default"]]), info = "6: the package's default design is installed")
smok_counts <- get(".SMOKING_COUNT_EXPOSURES", ns)
check_smok <- get(".check_legacy_exposure_bounds", ns)
for (k in names(ymls)) {
  defs <- yaml::read_yaml(ymls[[k]])$exposure_definitions
  nm <- vapply(defs, function(e) e$name, "")
  names(defs) <- nm
  mv <- vapply(defs, function(e) if (is.null(e$min_value)) NA_real_ else as.numeric(e$min_value), 0)
  off <- vapply(defs, function(e) if (is.null(e$offset)) 0 else as.numeric(e$offset), 0)
  expect_true(all(smok_counts %in% nm), info = paste0("6 (", k, "): the five smoking counts are there"))
  expect_true(all(mv[smok_counts] + off[smok_counts] >= 1),
              info = paste0("6 (", k, "): every smoking count's lowest value is >= 1 (simsmok() reads pr_relapse at quit years - 1)"))
  expect_true(check_smok(unname(defs)), info = paste0("6 (", k, "): ...which Design's check accepts"))
  expect_true(all(mv[intersect(c("fruit", "veg", "alcohol"), nm)] == 0),
              info = paste0("6 (", k, "): fruit, veg and alcohol keep their zeros"))
  expect_identical(defs[["smok_cig_ex"]]$qfun, "fqZABNB",
                   info = paste0("6 (", k, "): smok_cig_ex on CKutils' ZABNB (gamlss.dist 6.1.11 pZABNB returns 0 for q >= 1)"))
  bounded <- Filter(function(e) !is.null(e$min_value) && !is.null(e$max_value), defs)
  typ <- vapply(bounded, function(e) ft(e$qfun, e$family_type) %||% NA_character_, "")
  expect_false(anyNA(typ), info = paste0("6 (", k, "): every bounded exposure has a family type: ",
                                         paste(names(bounded)[is.na(typ)], collapse = ", ")))
}
# a design written before 2026-11-18 (smoking counts at min_value 0) is refused, not silently re-read as
# "draw the zeros"; offset 1 (the newer model's convention) is accepted; other exposures are not affected
expect_error(check_smok(list(list(name = "smok_quit_yrs", min_value = 0L, max_value = 90L))), "min_value: 1",
             info = "6: a smoking count at min_value 0 is refused, naming the remedy")
expect_error(check_smok(list(list(name = "smok_cig_ex", qfun = "qZABNB", min_value = 1L, max_value = 80L))), "qfun: fqZABNB",
             info = "6: a qfun whose pfun is gamlss.dist's broken pZABNB is refused, naming fqZABNB")
legacy <- tryCatch(check_smok(list(list(name = "smok_quit_yrs", qfun = "fqDPO", min_value = 0L, max_value = 90L),
                                   list(name = "smok_cig_ex", qfun = "qZABNB", min_value = 0L, max_value = 80L))),
                   error = conditionMessage)
expect_true(grepl("min_value: 1", legacy) && grepl("qfun: fqZABNB", legacy) && grepl("smok_quit_yrs, smok_cig_ex", legacy),
            info = "6: a pre-2026-11-18 design gets BOTH edits in one message")
expect_true(check_smok(list(list(name = "smok_cig_ex", qfun = "qZABNB", pfun = "fpZABNB", min_value = 1L, max_value = 80L))),
            info = "6: qZABNB with an explicit CKutils pfun is accepted")
expect_error(check_smok(list(list(name = "smok_cig_ex", min_value = 0.5))), "smok_cig_ex", info = "6: ...any of the five")
expect_true(check_smok(list(list(name = "smok_quit_yrs", min_value = 0L, offset = 1L))), info = "6: min 0 + offset 1 is accepted")
expect_true(check_smok(list(list(name = "alcohol", min_value = 0L))), info = "6: other exposures may draw 0")
expect_error(le_with(list(def("smok_dur_ex", "fqZINBI", 0, 90))), "INCLUSIVE",
             info = "6: Design$load_exposures() runs the check")

## ---- 7. the caches key on the PARF exposures' definitions and the bounds rule -------
key <- get(".exposure_bounds_key", ns)
defs0 <- yaml::read_yaml(ymls[["default"]])$exposure_definitions
d1 <- list(sim_prm = list(exposure_definitions = defs0))
defs2 <- defs0
i_alc <- which(vapply(defs2, function(e) e$name, "") == "alcohol")
defs2[[i_alc]]$max_value <- 200L
d2 <- list(sim_prm = list(exposure_definitions = defs2))
defs3 <- c(defs0, list(list(name = "hfss_energy", file_name = "hfss_energy", distribution = "continuous",
                            qfun = "qBCTo", min_value = 0L, max_value = 3000L)))
d3 <- list(sim_prm = list(exposure_definitions = defs3))
expect_identical(key(d1)$exposure_bounds_rule, get(".EXPOSURE_BOUNDS_RULE", ns), info = "7: the PARF key carries the rule")
expect_false(identical(key(d1), key(d2)), info = "7: ...and a PARF exposure's bound moves it")
expect_identical(key(d1), key(d3),
                 info = "7: an exposure the PARF never draws does not (the testing design shares the stock PARFs)")
defs4 <- lapply(rev(defs0), function(e) { e <- e[rev(names(e))]; lapply(e, function(v) if (is.integer(v)) as.double(v) else v) })
expect_identical(key(d1), key(list(sim_prm = list(exposure_definitions = defs4))),
                 info = "7: reordered entries / keys and 1 vs 1.0 keep the same PARF key")
uc <- SynthPop$private_methods$get_unique_characteristics
sim_prm0 <- yaml::read_yaml(ymls[["default"]])
sim_prm2 <- sim_prm0
sim_prm2$exposure_definitions <- defs2 # (not modifyList(): it merges unnamed lists into nothing)
expect_false(identical(uc(list(sim_prm = sim_prm0))$exposure_definitions,
                       uc(list(sim_prm = sim_prm2))$exposure_definitions),
             info = "7: the synthpop checksum's characteristics carry the exposure definitions")
# the synthpop checksum VALUE moves with the bounds rule: gen_checksum() run with its helpers stubbed
chk <- function(rule) {
  f <- SynthPop$private_methods$gen_checksum
  environment(f) <- list2env(list(
    private = list(get_unique_characteristics = function(d) list(n = 1L, exposure_definitions = "h"),
                   get_unique_LSOAs = function(d) c("E01000001", "E01000002")),
    .EXPOSURE_BOUNDS_RULE = rule), parent = ns)
  f(list())
}
expect_false(identical(chk("rule_A"), chk("rule_B")), info = "7: the synthpop checksum value moves with the bounds rule")
expect_identical(chk("rule_A"), chk("rule_A"), info = "7: ...and is deterministic for a given rule")
expect_true(grepl("rr_distributions,.exposure_bounds_key(design_)", ws(Disease$public_methods$initialize), fixed = TRUE),
            info = "7: the PARF directory key carries the exposure bounds key")
expect_true(grepl(".exposure_bounds_key(design_))),\".qs\"", ws(Disease$public_methods$gen_parf_files), fixed = TRUE),
            info = "7: the PARF intermediate .qs key carries it too")
# del_parf_file(invert = TRUE) removes this disease's superseded PARF directories and nothing else: not the current
# one, not another disease's, not a disease whose name extends this one ("chd" vs "chd_extra"; the old "^PARF_chd_"
# prefix matched it), and no .qs (a design whose directory another one evicted rebuilds it from its own .qs)
pd <- file.path(td, "parf"); dir.create(pd)
h <- function(ch) strrep(ch, 32L)
cur <- file.path(pd, paste0("PARF_chd_", h("a"))); old <- file.path(pd, paste0("PARF_chd_", h("b")))
ext <- file.path(pd, paste0("PARF_chd_extra_", h("c")))
invisible(lapply(c(cur, old, ext), dir.create))
qs_keep <- file.path(pd, paste0("PARF_chd_", h("d"), ".qs")); qs_old <- file.path(pd, paste0("PARF_chd_", h("e"), ".qs"))
qs_ckd <- file.path(pd, paste0("PARF_ckd_", h("f"), ".qs")); qs_ext <- file.path(pd, paste0("PARF_chd_extra_", h("1"), ".qs"))
invisible(file.create(c(qs_keep, qs_old, qs_ckd, qs_ext)))
dpf <- Disease$public_methods$del_parf_file
environment(dpf) <- list2env(list(self = list(name = "chd"), private = list(parf_dir = pd, parf_filenam = cur)),
                             parent = ns)
dpf(invert = TRUE)
expect_identical(dir.exists(c(cur, old, ext)), c(TRUE, FALSE, TRUE),
                 info = "7: the superseded PARF directory goes; the current one and 'chd_extra''s stay")
expect_identical(file.exists(c(qs_keep, qs_old, qs_ckd, qs_ext)), c(TRUE, TRUE, TRUE, TRUE),
                 info = "7: every .qs stays")

## ---- 8. the smoking engine: an ex-smoker's quit years are >= 1 (pre-existing, fixed 2026-11-19) -----------
# simsmok() reads pr_relapse(row, quit_yrs - 1) with no bounds check: quit years 0 read BEFORE the matrix
# (undefined behaviour that made identical runs differ). It now stops instead.
m_rel <- matrix(0.5, nrow = 10L, ncol = 15L)
smk <- function(q0) data.table(smok_status = c(3L, 3L), prb_smok_incid = 0, prb_smok_cess = 0, rankstat_simsmok = 0.9,
                               pid_mrk = c(TRUE, FALSE), sex = c(1L, 1L), qimd = c(1L, 1L),
                               smok_quit_yrs = c(q0, 0L), smok_dur = c(5L, 5L))
expect_error(simsmok(smk(0L), m_rel, 3L), "quit years must be >= 1", info = "8: simsmok() refuses an ex-smoker with quit years 0")
ok <- smk(1L)
simsmok(ok, m_rel, 3L)
expect_identical(ok$smok_quit_yrs, c(1L, 2L), info = "8: CONTROL: quit years 1 runs (no relapse at rn 0.9 > 0.5: 1 -> 2)")
smk_sc <- function(q0) data.table(smok_status_sc = c(3L, 3L), prb_smok_incid_sc = 0, prb_smok_cess_sc = 0,
                                  rankstat_simsmok = 0.9, pid_mrk_sc = c(TRUE, FALSE), sex = c(1L, 1L), qimd_sc = c(1L, 1L),
                                  smok_quit_yrs_sc = c(q0, 0L), smok_dur_sc = c(5L, 5L))
expect_error(simsmok_sc(smk_sc(0L), m_rel, 3L, 1:2), "quit years must be >= 1", info = "8: so does simsmok_sc()")
ok_sc <- smk_sc(1L)
simsmok_sc(ok_sc, m_rel, 3L, 1:2)
expect_identical(ok_sc$smok_quit_yrs_sc, c(1L, 2L), info = "8: CONTROL: simsmok_sc() with quit years 1 runs")
# the synthpop restores the invariant right after the age-12 rescale, whose as.integer() truncation turns a quit
# spell of 1 into 0 (only first rows: later rows are simsmok()'s), and before simsmok() is called
gs <- ws(SynthPop$private_methods$gen_synthpop)
at_rescale <- regexpr("smok_quit_yrs=as.integer(smok_quit_yrs/((smok_dur+smok_quit_yrs)/(age-12L)))", gs, fixed = TRUE)
at_guard <- regexpr("dtb[(pid_mrk)&smok_status%in%2:3&smok_quit_yrs<1L,`:=`(smok_quit_yrs,1L)]", gs, fixed = TRUE)
at_simsmok <- regexpr("simsmok(dtb,", gs, fixed = TRUE)
expect_true(at_rescale > 0L && at_guard > at_rescale && at_simsmok > at_guard,
            info = "8: gen_synthpop() guards quit years >= 1 after the age-12 rescale and before simsmok()")
# the rescale + guard arithmetic over a grid: no ex-smoker ends at quit years 0, and the history fits from age 12
g <- CJ(age = 20:99, smok_dur = 0:90, smok_quit_yrs = 1:90)
g[age - smok_quit_yrs < 12L, smok_quit_yrs := age - 12L]
g[age - smok_dur < 12L, smok_dur := age - 12L]
g[age - smok_dur - smok_quit_yrs < 12L, `:=`(
  smok_dur = as.integer(smok_dur / ((smok_dur + smok_quit_yrs) / (age - 12L))),
  smok_quit_yrs = as.integer(smok_quit_yrs / ((smok_dur + smok_quit_yrs) / (age - 12L))))]
n_zero <- g[smok_quit_yrs < 1L, .N]
expect_true(n_zero > 0L, info = "8: (non-vacuity) the rescale alone does produce quit years 0")
g[smok_quit_yrs < 1L, smok_quit_yrs := 1L]
expect_true(g[, all(smok_quit_yrs >= 1L & smok_dur + smok_quit_yrs <= age - 12L)],
            info = "8: ...and the guard restores >= 1 without breaking the age-12 budget")

unlink(td, recursive = TRUE)
rm("tt_mv_fpZINBI", "tt_broken_p", "tt_signvec_p", envir = globalenv())
