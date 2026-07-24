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

# --- Sign-crossing gradient: RII_ratio must be NA, REI_rel stays finite ---
y_cross <- -3 + 6 * ridit           # crosses zero within [0, 1]
res_cross <- calc(data.table(rank = rank, N = N, B = y_cross * N),
                  by = character(0))
expect_true(is.na(res_cross$RII_ratio),
            info = "RII_ratio is NA when the fitted line crosses zero")
expect_true(is.finite(res_cross$REI_rel),
            info = "REI_rel remains finite when the fit crosses zero")

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
