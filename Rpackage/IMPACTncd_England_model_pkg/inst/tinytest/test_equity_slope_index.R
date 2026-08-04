# Tests for the equity slope-index computation (calc_equity_slope_indices),
# the statistical core of export_equity_tables(). These verify that the
# closed-form population-weighted slope is identical to a weighted lm() fit and
# that the four reported indices behave as specified (incl. the NA edge cases).

library(tinytest)
library(IMPACTncdEngland)
library(data.table)

calc <- IMPACTncdEngland:::calc_equity_slope_indices

# --- Fixture: 10 DIMD deciles, UNEQUAL populations, a known linear gradient ---
set.seed(1L)
K <- 10L
rank <- 1:K                                  # 1 = most deprived .. 10 = least
N <- as.numeric(round(runif(K, 3e6, 4e6)))   # deliberately unequal decile sizes
totN <- sum(N)

# Manual ridit (population-weighted midpoint rank), least -> most deprived.
o <- order(-rank)
csum <- cumsum(N[o])
ridit_s <- (csum - N[o] / 2) / totN
ridit <- numeric(K); ridit[o] <- ridit_s

# per-capita benefit truly linear in the ridit (slope 5, pro-poor) + tiny noise
y <- 2 + 5 * ridit + rnorm(K, 0, 0.05)
B <- y * N

d <- data.table(rank = rank, N = N, B = B)
res <- calc(d, by = character(0))

# Reference: weighted OLS on the SAME ridit.
fit <- lm(y ~ ridit, weights = N)
beta_lm <- unname(coef(fit)["ridit"])
p0 <- unname(predict(fit, newdata = data.frame(ridit = 0)))  # least deprived
p1 <- unname(predict(fit, newdata = data.frame(ridit = 1)))  # most deprived

# 1. Closed-form weighted slope == lm() weighted slope (this is the SII).
expect_equal(res$AEI_per100k / 1e5, beta_lm,
             info = "closed-form slope equals weighted lm() coefficient")

# 2. AEI_total is the slope rescaled to the whole reference population.
expect_equal(res$AEI_total, beta_lm * totN,
             info = "AEI_total == slope * sum(N)")

# 3. Relative index REI_rel == slope / population-weighted mean benefit.
expect_equal(res$REI_rel, beta_lm / (sum(B) / totN),
             info = "REI_rel == SII / weighted mean")

# 4. RII_ratio == fitted(most deprived) / fitted(least deprived).
expect_equal(res$RII_ratio, p1 / p0,
             info = "RII_ratio == fitted(r=1)/fitted(r=0)")

# 5. Ridit midpoints match the manual cumulative-population formula. Recompute
#    them the way the helper does and compare (fit slope must match the manual).
expect_true(abs(beta_lm - 5) < 0.2,
            info = "recovered slope is close to the true gradient (5)")

# --- Sign-crossing gradient: RII_ratio must be NA ---
# NOTE this fixture is symmetric about r = 1/2, and the population-weighted mean
# ridit is exactly 1/2, so its weighted mean benefit is exactly zero (to
# floating point, -2.1e-17). REI_rel is therefore NOT defined here either: the
# old `ybar != 0` guard let a value that is numerically zero through and
# published beta/-2.1e-17 ~ -2.9e17 as a "finite" relative index.
y_cross <- -3 + 6 * ridit           # crosses zero within [0, 1]
res_cross <- calc(data.table(rank = rank, N = N, B = y_cross * N),
                  by = character(0))
expect_true(is.na(res_cross$RII_ratio),
            info = "RII_ratio is NA when the fitted line crosses zero")
expect_true(is.na(res_cross$REI_rel),
            info = "REI_rel is NA when the weighted mean benefit is zero")

# A line that crosses zero but has a genuinely POSITIVE mean keeps REI_rel:
# it is the sign of the mean that matters, not the crossing.
res_cross_pos <- calc(data.table(rank = rank, N = N, B = (-1 + 6 * ridit) * N),
                      by = character(0))
expect_true(is.na(res_cross_pos$RII_ratio),
            info = "RII_ratio still NA when the fitted line crosses zero")
expect_true(is.finite(res_cross_pos$REI_rel),
            info = "REI_rel is defined when the weighted mean benefit is > 0")

# --- Net-harm gradient: the relative indices must be WITHHELD, not inverted ---
# Every group loses, but the most deprived lose LEAST, so the gradient is
# pro-poor and the absolute index must stay positive. beta/ybar with ybar < 0
# flips sign, and fit1/fit0 with both ends negative inverts the ordering, so a
# pro-poor result would be published as pro-rich. Both are therefore NA.
res_harm <- calc(data.table(rank = rank, N = N, B = (-8 + 6 * ridit) * N),
                 by = character(0))
expect_true(res_harm$AEI_per100k > 0,
            info = "net-harm, most deprived lose least -> absolute index pro-poor")
expect_true(is.na(res_harm$REI_rel),
            info = "REI_rel withheld when the mean benefit is negative")
expect_true(is.na(res_harm$RII_ratio),
            info = "RII_ratio withheld when both fitted ends are negative")

# --- total_benefit and fit_R2 ---
expect_equal(res$total_benefit, sum(B),
             info = "total_benefit == sum of the group benefits")
expect_equal(res$fit_R2, summary(fit)$r.squared,
             info = "fit_R2 == weighted lm() R-squared")
res_curve <- calc(
  data.table(rank = rank, N = N, B = (2 + 5 * ridit^3) * N), by = character(0))
expect_true(res_curve$fit_R2 < res$fit_R2,
            info = "fit_R2 falls when the gradient departs from linearity")
# calc_equity_slope_indices() has no N > 0 precondition of its own, so a
# zero-population group must not divide by zero.
expect_true(is.na(calc(data.table(rank = rank, N = 0, B = 0),
                       by = character(0))$fit_R2),
            info = "fit_R2 is NA (not NaN/Inf) when every group has N == 0")

# --- Flat gradient: slope 0, so AEI == 0, REI_rel == 0, RII_ratio == 1 ---
res_flat <- calc(data.table(rank = rank, N = N, B = 7 * N), by = character(0))
expect_equal(res_flat$AEI_per100k, 0, info = "flat gradient -> zero absolute index")
expect_equal(res_flat$REI_rel, 0, info = "flat gradient -> zero relative slope")
expect_equal(res_flat$RII_ratio, 1, info = "flat gradient -> RII ratio of 1")

# --- Grouping: two groups with different true slopes are fit independently ---
dg <- rbind(
  data.table(rank = rank, N = N, B = (2 + 5 * ridit) * N, grp = "a"),
  data.table(rank = rank, N = N, B = (1 + 3 * ridit) * N, grp = "b")
)
resg <- calc(dg, by = "grp")
expect_equal(nrow(resg), 2L, info = "one row per group")
expect_equal(resg[grp == "a", AEI_per100k] / 1e5, 5, info = "group a slope == 5")
expect_equal(resg[grp == "b", AEI_per100k] / 1e5, 3, info = "group b slope == 3")

# --- Degenerate group (single decile) -> indices NA, no error ---
res_one <- calc(data.table(rank = 1L, N = 1e6, B = 5e5), by = character(0))
expect_true(is.na(res_one$AEI_total),
            info = "single-decile group yields NA (no gradient to fit)")

# --- Sign convention: pro-rich gradient (benefit falls with deprivation) < 0 ---
y_prorich <- 10 - 5 * ridit
res_pr <- calc(data.table(rank = rank, N = N, B = y_prorich * N),
               by = character(0))
expect_true(res_pr$AEI_per100k < 0,
            info = "pro-rich gradient gives a negative absolute index")


# ===========================================================================
# Gradient-axis selection (dimd vs qimd)
# ===========================================================================
validate_eq <- IMPACTncdEngland:::validate_equity_strata
build_plans <- IMPACTncdEngland:::build_equity_plans
dep_rank    <- IMPACTncdEngland:::deprivation_rank

# --- validate_equity_strata: structural rules only -------------------------
# Any variable the summaries carry can be an output stratum; whether it really
# exists is checked against the data in export_equity_tables(), not here.
expect_silent(validate_eq(list("year", c("year", "sex"))))
expect_silent(validate_eq(list(c("year", "dimd"), c("year", "qimd"))))
expect_silent(validate_eq(list(c("year", "sex", "dimd", "qimd"))))
expect_silent(validate_eq(list(c("year", "agegrp", "sex", "qimd"))))
expect_silent(validate_eq(list(c("year", "ethnicity"), c("year", "sha"))))

expect_error(validate_eq(list(c("sex", "dimd"))),
             pattern = "must include",
             info = "an equity stratum without year is rejected")
expect_error(validate_eq(list(c("year", "sex", "sex"))),
             pattern = "duplicated")
expect_error(validate_eq(list(c("year", "mc"))), pattern = "reserved column",
             info = "mc is the Monte-Carlo axis, not a stratum")
expect_error(validate_eq(list(c("year", "scenario"))),
             pattern = "reserved column")
expect_error(validate_eq(list(character(0))), pattern = "non-empty")
expect_error(validate_eq(list(1:2)), pattern = "non-empty")
expect_error(validate_eq("year"), pattern = "must be a list")

# --- build_equity_plans: one plan per output table -------------------------
# No gradient token -> dimd, and the historical filename suffix is preserved.
p_implicit <- build_plans(list("year", c("year", "sex")))
expect_equal(length(p_implicit), 2L)
expect_equal(vapply(p_implicit, `[[`, character(1), "gradient"),
             c("dimd", "dimd"), info = "no token falls back to dimd")
expect_equal(vapply(p_implicit, `[[`, character(1), "suffix"),
             c("year", "year-sex"),
             info = "implicit dimd keeps the historical filename suffix")

# One token -> that gradient, echoed in the suffix.
p_q <- build_plans(list(c("year", "qimd"), c("year", "sex", "qimd")))
expect_equal(vapply(p_q, `[[`, character(1), "gradient"), c("qimd", "qimd"))
expect_equal(vapply(p_q, `[[`, character(1), "suffix"),
             c("year-qimd", "year-sex-qimd"))
expect_equal(p_q[[2L]]$out_vars, c("year", "sex"),
             info = "the gradient token is not an output stratum")

# Both tokens in one entry -> one table per gradient.
p_both <- build_plans(list(c("year", "dimd", "qimd")))
expect_equal(length(p_both), 2L, info = "both tokens -> two tables")
expect_equal(vapply(p_both, `[[`, character(1), "gradient"), c("dimd", "qimd"))
expect_equal(vapply(p_both, `[[`, character(1), "suffix"),
             c("year-dimd", "year-qimd"),
             info = "the two tables get distinct filenames")

# Token order in the entry does not matter; gradients come out canonical.
expect_equal(build_plans(list(c("qimd", "year", "dimd"))),
             build_plans(list(c("year", "dimd", "qimd"))))

# --- deprivation_rank: 1 = most deprived, and no silent mis-mapping ---------
dimd_lv <- c("1 most deprived", as.character(2:9), "10 least deprived")
qimd_lv <- c("1 most deprived", as.character(2:4), "5 least deprived")
expect_equal(dep_rank(factor(dimd_lv, levels = dimd_lv), "dimd"), 1:10)
expect_equal(dep_rank(qimd_lv, "qimd"), 1:5)
expect_equal(dep_rank(rev(dimd_lv), "dimd"), 10:1,
             info = "rank follows the label, not the row order")

# The trap this guards: quintile labels share 4 of 5 labels with decile labels,
# so a column named `dimd` holding quintiles maps to 1,2,3,4,NA -- a silently
# mis-ordered gradient. It must stop, not return a plausible wrong number.
expect_error(dep_rank(qimd_lv, "dimd"), pattern = "unexpected `dimd` level",
             info = "quintile labels in a dimd column are rejected")
expect_error(dep_rank(c("1 most deprived", "banana"), "dimd"),
             pattern = "unexpected `dimd` level")

# --- The decile -> quintile collapse preserves a linear gradient exactly ---
# Collapsing adjacent groups keeps the population-weighted ridit midpoint
# (the algebra cancels for ANY population split), so a benefit that is exactly
# linear in the ridit yields the IDENTICAL slope index over 5 quintiles as over
# 10 deciles. This is what makes `qimd` a legitimate coarser gradient, and what
# makes deriving qimd from dimd by summing exact rather than approximate.
dimd_to_qimd <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L)
y_lin <- 2 + 5 * ridit                       # no noise: exactly linear
dec <- data.table(rank = rank, N = N, B = y_lin * N)
# Both rank scales run 1 = most deprived, so deciles 1,2 -> quintile 1 directly.
qui <- dec[, .(N = sum(N), B = sum(B)), keyby = .(rank = dimd_to_qimd[rank])]
expect_equal(nrow(qui), 5L)
expect_equal(calc(dec, by = character(0))$AEI_per100k,
             calc(qui, by = character(0))$AEI_per100k,
             info = "quintile SII == decile SII for an exactly linear gradient")
expect_equal(calc(dec, by = character(0))$AEI_total,
             calc(qui, by = character(0))$AEI_total,
             info = "quintile AEI_total == decile AEI_total (same reference pop)")
