# Common random numbers for initial disease prevalence (GitHub issue #8).
#
# THE BUG THIS PINS
# -----------------
# `Disease$set_init_prvl()` reallocates initial prevalent cases within each
# prevalence stratum by weighted sampling without replacement. It used to call
# `wrswoR::sample_int_expj()` once per stratum (`by = strata`), all strata on
# ONE base-R stream, visited in order of first appearance. `sample_int_expj()`
# consumes a data-dependent number of variates, so any change to one stratum
# shifted the stream for every stratum after it. Synthpop pids are ordered by
# age in the initial year, so the strata of future entrants come before every
# `init_year` stratum: a scenario that touched only entrants from `sc_year` on
# re-drew WHICH simulants were prevalent in `init_year`. On the real model
# (auxil/verify_init_prvl_crn.R, n = 20000, sc_year 2015) 0.2-44% of pre-2015
# person-years differed from the baseline under an exposure-only scenario and
# ~75% under a forced-deprivation one; a no-op scenario matched exactly.
# The duration draw had the same flaw (uniforms indexed among the PREVALENT
# rows only): reverting just that draw re-opens ~68% of the forced-deprivation
# discrepancy on its own.
#
# The reallocation now goes through select_wrswor_by(): one exponential variate
# per ROW, key = variate / weight, keep the `size` smallest keys per stratum
# (Efraimidis & Spirakis 2006). Pinned here:
#   1-6. its contract, loud failures and edge cases;
#   2.   the CRN property itself (other strata perturbed -> stratum unchanged);
#   7.   its distribution equals successive sampling (sample.int(prob =));
#   8.   a NEGATIVE CONTROL: the replaced per-stratum sample_int_expj() pattern
#        fails the very property section 2 asserts, so section 2 is not
#        vacuous;
#   9.   a tripwire that set_init_prvl() still uses it, and still draws its
#        duration and diagnosis random numbers per row;
#  10.   that set() coerces a double duration into the integer *_prvl column
#        exactly as the `:=` it replaced did.
# The duration and diagnosis draws are inline in set_init_prvl() and need the
# model data, so beyond the section 9 tripwire they are covered only by the
# production verification: auxil/verify_init_prvl_crn.R.

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

select_wrswor_by <- IMPACTncdEngland:::select_wrswor_by

# A small multi-stratum table: 4 strata of 6, 4, 5 and 6 rows, sizes 2, 0, 5
# and 2 (stratum "c" has size == n).
mk <- function() {
  data.table(
    s = rep(c("a", "b", "c", "d"), c(6L, 4L, 5L, 6L)),
    ncases = rep(c(2L, 0L, 5L, 2L), c(6L, 4L, 5L, 6L)),
    w = c(0.5, 1, 1, 2, 3, 0.7, 1, 1, 1, 1, 1, 2, 3, 4, 5,
          1.5, 0.6, 2, 1, 0.8, 1.2)
  )
}
rx_of <- function(n, seed = 1L) {
  set.seed(seed)
  stats::rexp(n)
}

# ===========================================================================
# 1. Basic contract
# ===========================================================================
dt <- mk()
rx <- rx_of(nrow(dt))
sel <- select_wrswor_by(dt, by = "s", size = "ncases", weight = "w", rexp = rx)
expect_true(is.integer(sel), info = "returns integer row numbers")
expect_identical(sel, sort(sel), info = "row numbers are ascending")
expect_identical(
  as.vector(table(factor(dt$s[sel], levels = c("a", "b", "c", "d")))),
  c(2L, 0L, 5L, 2L),
  info = "exactly `size` rows per stratum"
)
expect_true(all(11:15 %in% sel), info = "size == n selects the whole stratum")
expect_identical(nrow(dt), 21L, info = "dt is not modified (no rows added)")
expect_identical(names(dt), c("s", "ncases", "w"),
                 info = "dt is not modified (no columns added)")

# The selection is the `size` smallest rexp / weight within each stratum.
key <- rx / dt$w
exp_a <- which(dt$s == "a")[order(key[dt$s == "a"])][1:2]
expect_identical(sel[dt$s[sel] == "a"], sort(exp_a),
                 info = "smallest keys are selected")

# ===========================================================================
# 2. Common random numbers
# ===========================================================================
# Changing another stratum's weights, sizes or membership leaves a stratum
# untouched; this is the property issue #8 needs. Checked both for a stratum
# that comes AFTER the perturbed ones in row order ("d", as the init_year
# population comes after the entrants in sp$pop) and for one BEFORE them ("a").
dt1 <- mk()
dt1[s == "a", w := w * c(3, 0.2, 1, 0.5, 2, 1)]
dt1[s == "b", ncases := 3L] # four new candidate rows ahead of "d"
dt1[s == "c", ncases := 4L]
sel1 <- select_wrswor_by(dt1, "s", "ncases", "w", rx)
expect_identical(sel1[dt1$s[sel1] == "d"], sel[dt$s[sel] == "d"],
                 info = "CRN: earlier strata perturbed, later stratum d unchanged")
expect_false(
  identical(sel1[dt1$s[sel1] == "a"], sel[dt$s[sel] == "a"]) &&
    identical(sel1[dt1$s[sel1] == "b"], sel[dt$s[sel] == "b"]),
  info = "non-vacuity: the perturbation does change the perturbed strata"
)
dt2 <- mk()
dt2[s == "c", w := w * c(5, 0.1, 1, 3, 0.2)]
dt2[s == "b", ncases := 3L]
dt2[s == "c", ncases := 4L]
sel2 <- select_wrswor_by(dt2, "s", "ncases", "w", rx)
expect_identical(sel2[dt2$s[sel2] == "a"], sel[dt$s[sel] == "a"],
                 info = "CRN: other strata perturbed, stratum a unchanged")
dt3 <- mk()
dt3[c(7L, 12L), s := "a2"] # two rows leave their strata for a new one
dt3[s == "a2", ncases := 1L]
dt3[s == "c", ncases := 4L]
sel3 <- select_wrswor_by(dt3, "s", "ncases", "w", rx)
expect_identical(sel3[dt3$s[sel3] == "a"], sel[dt$s[sel] == "a"],
                 info = "CRN: stratum membership elsewhere changed, stratum a unchanged")

# Invariance to row order: permuting rows together with their variates
# selects the same rows.
set.seed(3L)
p <- sample.int(nrow(dt))
selp <- select_wrswor_by(dt[p], "s", "ncases", "w", rx[p])
expect_identical(sort(p[selp]), sel, info = "row-order invariant")

# Several `by` columns, factors and NA as a stratum value.
dtm <- data.table(
  y = rep(13:14, each = 6L),
  g = factor(rep(c("m", "f", NA), 4L)),
  ncases = 1L,
  w = 1
)
selm <- select_wrswor_by(dtm, c("y", "g"), "ncases", "w", rx_of(nrow(dtm), 4L))
expect_identical(length(selm), 6L,
                 info = "one per (y, g) stratum incl. NA strata")
expect_identical(nrow(unique(dtm[selm, .(y, g)])), 6L,
                 info = "NA stratum treated as one group")

# ===========================================================================
# 3. Zero weights: sample.int() semantics
# ===========================================================================
dz <- data.table(s = 1L, ncases = 2L, w = c(0, 1, 0, 1, 2))
for (sd in 1:50) {
  sz <- select_wrswor_by(dz, "s", "ncases", "w", rx_of(5L, sd))
  if (any(dz$w[sz] == 0)) break
}
expect_true(all(dz$w[sz] > 0),
            info = "zero-weight rows never selected while positives remain")
expect_error(
  select_wrswor_by(data.table(s = 1L, ncases = 2L, w = c(0, 0, 1)),
                   "s", "ncases", "w", rx_of(3L)),
  pattern = "too few positive weights"
)
# A zero variate on a zero weight would be 0/0 = NaN; it must stay unselectable.
expect_identical(
  select_wrswor_by(data.table(s = 1L, ncases = 1L, w = c(0, 1)),
                   "s", "ncases", "w", c(0, 5)),
  2L,
  info = "0/0 key handled as Inf"
)

# ===========================================================================
# 4. Loud failures (sample_int_expj() silently misbehaves on all of these)
# ===========================================================================
dt <- mk()
rx <- rx_of(nrow(dt))
bad_w <- function(v) {
  d <- mk()
  d[1L, w := v]
  select_wrswor_by(d, "s", "ncases", "w", rx)
}
expect_error(bad_w(NA_real_), pattern = "finite, non-negative and not NA")
expect_error(bad_w(NaN), pattern = "finite, non-negative and not NA")
expect_error(bad_w(-1), pattern = "finite, non-negative and not NA")
expect_error(bad_w(Inf), pattern = "finite, non-negative and not NA")
# ... but weights of rows in size-0 strata are never looked at.
d_ok <- mk()
d_ok[s == "b", w := NA_real_]
expect_silent(select_wrswor_by(d_ok, "s", "ncases", "w", rx))

bad_rx <- function(v) {
  r <- rx
  r[1L] <- v
  select_wrswor_by(mk(), "s", "ncases", "w", r)
}
expect_error(bad_rx(NA_real_), pattern = "'rexp' must be finite")
expect_error(bad_rx(-0.1), pattern = "'rexp' must be finite")
expect_error(bad_rx(Inf), pattern = "'rexp' must be finite")
expect_error(select_wrswor_by(mk(), "s", "ncases", "w", rx[-1L]),
             pattern = "one variate per row")

bad_size <- function(v, rows = 1L) {
  d <- mk()
  d[, ncases := as.numeric(ncases)]
  d[rows, ncases := v]
  select_wrswor_by(d, "s", "ncases", "w", rx)
}
expect_error(bad_size(NA_real_), pattern = "'size' must not be NA")
expect_error(bad_size(-1), pattern = "non-negative")
expect_error(bad_size(2.7), pattern = "finite whole numbers")
expect_error(bad_size(Inf), pattern = "finite whole numbers")
expect_error(bad_size(3), pattern = "constant within each group")
expect_error(bad_size(7, rows = 1:6),
             pattern = "cannot take a sample larger than the population")
# Whole-number doubles are accepted.
expect_identical(bad_size(2), select_wrswor_by(mk(), "s", "ncases", "w", rx))

expect_error(select_wrswor_by(as.data.frame(mk()), "s", "ncases", "w", rx),
             pattern = "must be a data.table")
expect_error(select_wrswor_by(mk(), "nope", "ncases", "w", rx),
             pattern = "not found")
expect_error(select_wrswor_by(mk(), "s", "ncases", "nope", rx),
             pattern = "not found")
expect_error(select_wrswor_by(mk(), character(0), "ncases", "w", rx),
             pattern = "non-empty")

# ===========================================================================
# 5. Empty inputs
# ===========================================================================
expect_identical(select_wrswor_by(mk()[0L], "s", "ncases", "w", numeric(0)),
                 integer(0))
d0 <- mk()
d0[, ncases := 0L]
expect_identical(select_wrswor_by(d0, "s", "ncases", "w", rx), integer(0))

# ===========================================================================
# 6. data.table numeric rounding is restored
# ===========================================================================
old <- getNumericRounding()
setNumericRounding(2L)
sel_r <- select_wrswor_by(mk(), "s", "ncases", "w", rx)
expect_identical(getNumericRounding(), 2L,
                 info = "user's rounding setting restored")
setNumericRounding(old)
expect_identical(sel_r, select_wrswor_by(mk(), "s", "ncases", "w", rx),
                 info = "result independent of the rounding setting")
# Two keys that differ only in the low mantissa bits: under rounding they
# would tie and row order would pick row 1; on all 64 bits row 2 is smaller.
setNumericRounding(2L)
expect_identical(
  select_wrswor_by(data.table(s = 1L, ncases = 1L, w = c(1, 1)),
                   "s", "ncases", "w", c(1 + 2^-45, 1)),
  2L,
  info = "keys compared on all 64 bits under a user's numeric rounding"
)
setNumericRounding(old)

# ===========================================================================
# 7. Distribution: inclusion probabilities, many strata in one call
# ===========================================================================
# Exact successive-sampling inclusion probabilities for w = (0.5, 1, 1, 2.5, 4),
# size 2, from enumerating all ordered draws -- the law of
# sample.int(5, 2, prob = w) and of wrswoR::sample_int_expj().
w5 <- c(0.5, 1, 1, 2.5, 4)
incl_exact <- local({
  n <- length(w5)
  incl <- numeric(n)
  for (i in seq_len(n)) for (j in seq_len(n)[-i]) {
    p <- w5[i] / sum(w5) * w5[j] / (sum(w5) - w5[i])
    incl[c(i, j)] <- incl[c(i, j)] + p
  }
  incl
})
R <- 2e4L
dd <- data.table(g = rep(seq_len(R), each = 5L), pos = rep(1:5, R),
                 w = rep(w5, R), ncases = 2L)
sd <- select_wrswor_by(dd, "g", "ncases", "w", rx_of(nrow(dd), 99L))
incl_mc <- tabulate(dd$pos[sd], nbins = 5L) / R
z <- (incl_mc - incl_exact) / sqrt(incl_exact * (1 - incl_exact) / R)
expect_true(all(abs(z) < 4.5),
            info = "MC inclusion probabilities within 4.5 SE of exact")

# ===========================================================================
# 8. Negative control: the replaced per-stratum sample_int_expj() pattern
# ===========================================================================
# Forty strata of eight rows, three cases each. Perturbation: stratum 1 (first
# in row order, like a future entrant stratum) needs ONE more case -- what a
# scenario moving an entrant into it does to `ncases`. Strata 2..40 are not
# touched, so a CRN-safe sampler must select exactly the same rows in them.
set.seed(7L)
nc0 <- data.table(
  s = rep(1:40, each = 8L),
  row = seq_len(320L),
  w = round(stats::runif(320L, 0.5, 3), 2),
  ncases = 3L
)
nc1 <- copy(nc0)[s == 1L, ncases := 4L]
untouched <- function(d, rows) sort(rows[d$s[rows] != 1L])

# Five seeds each. A shared stream can re-synchronise after an offset, so for
# a single seed the old pattern occasionally re-draws nothing (~8% of seeds on
# this fixture); "any of five" makes the control robust to the seed choice.
seeds <- 1:5

# The code set_init_prvl() used before the fix (one base-R stream, seeded once
# per disease, shared across strata in order of first appearance). wrswoR is
# no longer called by the package, so do not let its absence fail the suite.
if (requireNamespace("wrswoR", quietly = TRUE)) {
  old_way <- function(d, seed) {
    set.seed(seed)
    d[ncases > 0L][,
      .(row = row[wrswoR::sample_int_expj(unique(.N), unique(ncases), w)]),
      by = s
    ]$row
  }
  redrawn <- vapply(seeds, function(sd) {
    !identical(untouched(nc0, old_way(nc0, sd)), untouched(nc1, old_way(nc1, sd)))
  }, logical(1))
  expect_true(
    any(redrawn),
    info = "negative control: old per-stratum sample_int_expj() re-draws untouched strata"
  )
}

# The fix, with the same per-row variates in both runs.
kept <- vapply(seeds, function(sd) {
  rx40 <- rx_of(nrow(nc0), sd)
  identical(
    untouched(nc0, select_wrswor_by(nc0, "s", "ncases", "w", rx40)),
    untouched(nc1, select_wrswor_by(nc1, "s", "ncases", "w", rx40))
  )
}, logical(1))
expect_true(all(kept),
            info = "select_wrswor_by() leaves untouched strata exactly as they were")

# ===========================================================================
# 9. Tripwire: set_init_prvl() must keep its CRN-safe draws
# ===========================================================================
# deparse() drops comments, and whitespace is squashed, so only a change to
# the code itself trips these. They pin the inline draws no unit test reaches:
# a revert of any of them re-opens issue #8 while sections 1-8 stay green.
squash <- function(x) gsub("[[:space:]]+", "", paste(x, collapse = ""))
body_txt <- squash(
  deparse(IMPACTncdEngland::Disease$public_methods$set_init_prvl)
)
expect_false(grepl("sample_int_expj", body_txt, fixed = TRUE),
             info = "set_init_prvl() no longer calls sample_int_expj()")
for (pattern in c(
  "select_wrswor_by(",                     # reallocation by per-row keys
  "rexp = dqrexp(nrow(sp$pop))",           # ... one variate per row
  "u_dur <- dqrunif(nrow(sp$pop))",        # duration: one uniform per row
  "list(p = u_dur[rows])",                 # ... used by explicit row number
  "dqset.seed(private$seed, stream = sp$mc * 10 + 4L)" # diagnosis stream
)) {
  expect_true(grepl(squash(pattern), body_txt, fixed = TRUE),
              info = paste("set_init_prvl() still contains:", pattern))
}

# ===========================================================================
# 10. set() coerces like the `:=` it replaced (duration write-back)
# ===========================================================================
# The q*() duration functions return whole numbers as doubles; written into
# the integer *_prvl column they must arrive as the same integers, silently.
d_set <- data.table(p = c(0L, 1L, 1L, 0L))
d_col <- copy(d_set)
rows <- which(d_set$p == 1L)
expect_silent(set(d_set, i = rows, j = "p", value = 2L + c(3, 7)))
d_col[rows, p := 2L + c(3, 7)]
expect_identical(d_set, d_col, info = "set() and := give the same integer column")
expect_identical(d_set$p, c(0L, 5L, 9L, 0L))
