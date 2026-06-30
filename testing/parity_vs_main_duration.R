# Parity test of the ported stochastic forward-duration draw vs the original
# (main) engine. Because the get_dur_forward linear predictors are identical
# between main and the branch, duration parity reduces to whether the ZANBI
# scalar functions agree. This compares, over the REAL asthma duration parameter
# space (mu/sigma/nu reconstructed from the fitted coefficients across
# age/sex/dIMD/year), the branch's CKutils fqZANBI_scalar/fpZANBI_scalar against
# main's my_qZANBI_scalar/my_pZANBI_scalar (transcribed in the .cpp), and against
# gamlss::qZANBI as a canonical reference.
suppressMessages({library(Rcpp); library(data.table)})
sourceCpp("testing/parity_vs_main_duration.cpp")

co <- yaml::read_yaml("inputs/disease_burden/asthma_dur_forward.yaml")
stopifnot(!is.null(co$mu), !is.null(co$sigma), !is.null(co$nu))
lp <- function(b, age, year, sex, dimd)
  b[["intercept"]] + b[["log(age)"]] * log(age) + b[["log(year)"]] * log(year) +
  b[["sex"]][sex] + b[["dimd"]][dimd]

grid <- CJ(age = seq(20, 90, 2), sex = 1:2, dimd = 1:10, year = c(2000, 2020, 2043))
grid[, mu    := exp(lp(co$mu,    age, year, sex, dimd))]
grid[, sigma := exp(lp(co$sigma, age, year, sex, dimd))]
grid[, nu    := plogis(lp(co$nu, age, year, sex, dimd))]   # antilogit link for nu

P <- seq(0.0005, 0.9995, length.out = 50)
test <- grid[rep(seq_len(.N), each = length(P))]
test[, p := rep(P, nrow(grid))]
test[, q_m := q_main(p, mu, sigma, nu)]
test[, q_b := q_branch(p, mu, sigma, nu)]

cat(sprintf("Asthma parameter sets: %d (age x sex x dIMD x year); probs/set: %d; comparisons: %d\n",
            nrow(grid), length(P), nrow(test)))
cat(sprintf("mu in [%.2f, %.2f], sigma in [%.2f, %.2f], nu in [%.3f, %.3f]\n",
            min(grid$mu), max(grid$mu), min(grid$sigma), max(grid$sigma), min(grid$nu), max(grid$nu)))

nm <- test[q_m == q_b, .N]
cat(sprintf("\n[QUANTILE] main my_qZANBI_scalar  vs  branch fqZANBI_scalar:\n"))
cat(sprintf("  exact integer match: %d / %d  (%.4f%%)   max |diff in drawn years|: %d\n",
            nm, nrow(test), 100 * nm / nrow(test), test[, max(abs(q_m - q_b))]))
mm <- test[q_m != q_b]
if (nrow(mm)) { cat("  --- largest discrepancies ---\n"); print(head(mm[, .(p, mu, sigma, nu, q_m, q_b)][order(-abs(q_m - q_b))], 12)) }

qcap <- pmin(test$q_m, 1000L)
test[, pc_m := p_main(qcap, mu, sigma, nu)]
test[, pc_b := p_branch(qcap, mu, sigma, nu)]
cat(sprintf("\n[CDF] main my_pZANBI_scalar  vs  branch fpZANBI_scalar:  max |diff|: %.3e\n",
            test[, max(abs(pc_m - pc_b))]))

ref <- if (requireNamespace("gamlss.dist", quietly = TRUE)) "gamlss.dist" else
       if (requireNamespace("gamlss", quietly = TRUE)) "gamlss" else NA_character_
if (!is.na(ref)) {
  qref <- get("qZANBI", asNamespace(ref))
  u <- unique(test[, .(mu, sigma, nu, p)])
  u[, q_ref := as.integer(qref(p, mu = mu, sigma = sigma, nu = nu))]
  u[, q_b := q_branch(p, mu, sigma, nu)]
  cat(sprintf("\n[REF] branch fqZANBI_scalar  vs  %s::qZANBI (canonical):\n", ref))
  cat(sprintf("  exact match: %d / %d  (%.4f%%)   max |diff|: %d\n",
              u[q_ref == q_b, .N], nrow(u), 100 * u[q_ref == q_b, .N] / nrow(u), u[, max(abs(q_ref - q_b))]))
}

verdict <- if (test[, max(abs(q_m - q_b))] == 0) "PASS - branch draw is identical to main over the full asthma parameter space" else
           "REVIEW - discrepancies found (see above)"
cat(sprintf("\nRESULT: %s\nPARITY_DONE\n", verdict))
