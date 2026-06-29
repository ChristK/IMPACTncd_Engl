// Parity harness for the stochastic forward-duration port.
//
// Question: does the branch's ported duration draw -- which now uses CKutils'
// header-only scalar ZANBI (fqZANBI_scalar / fpZANBI_scalar) -- reproduce the
// ORIGINAL engine's draw, which used main's my_qZANBI_scalar / my_pZANBI_scalar?
// The get_dur_forward linear predictors (mu/sigma/nu) are byte-identical between
// main and the branch (copied verbatim), so duration parity reduces to whether
// these two ZANBI implementations agree over the asthma parameter space.
//
// main's scalar functions below are transcribed verbatim (lower_tail=TRUE,
// log_p=FALSE paths, which is all get_dur_forward uses) from
//   origin/main:Rpackage/IMPACTncd_Engl_model_pkg/src/NBI_distribution.cpp
//
// Build via Rcpp::sourceCpp(); LinkingTo CKutils supplies <distr_ZANBI.h>.

// [[Rcpp::depends(CKutils)]]
#include <Rcpp.h>
#include <cmath>
#include <distr_ZANBI.h>   // branch side: CKutils fqZANBI_scalar / fpZANBI_scalar
using namespace Rcpp;

// ----- main reference (verbatim, lower_tail=true / log_p=false) ---------------
static double main_pNBI_scalar(int q, double mu, double sigma) {
  return (sigma < 1e-4) ? R::ppois(q, mu, 1, 0)
                        : R::pnbinom_mu(q, 1.0 / sigma, mu, 1, 0);
}
static int main_qNBI_scalar(double p, double mu, double sigma) {
  double q = (sigma < 1e-4) ? R::qpois(p, mu, 1, 0)
                            : R::qnbinom_mu(p, 1.0 / sigma, mu, 1, 0);
  return (int) q;
}
static int main_qZANBI_scalar(double p, double mu, double sigma, double nu) {
  double p_ = p;                              // lower_tail=true, log_p=false
  p_ = (p_ - nu) / (1.0 - nu) - 1e-10;
  double cdf0 = main_pNBI_scalar(0, mu, sigma);
  p_ = cdf0 * (1.0 - p_) + p_;
  if (p_ < 0.0) p_ = 0.0;
  return main_qNBI_scalar(p_, mu, sigma);
}
static double main_pZANBI_scalar(int q, double mu, double sigma, double nu) {
  double cdf0 = main_pNBI_scalar(0, mu, sigma);
  double cdf1 = main_pNBI_scalar(q, mu, sigma);
  double cdf3 = nu + ((1.0 - nu) * (cdf1 - cdf0) / (1.0 - cdf0));
  return q == 0 ? nu : cdf3;
}

// ----- exposed element-wise comparators --------------------------------------
// [[Rcpp::export]]
IntegerVector q_main(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector nu) {
  int n = p.size(); IntegerVector o(n);
  for (int i = 0; i < n; ++i) o[i] = main_qZANBI_scalar(p[i], mu[i], sigma[i], nu[i]);
  return o;
}
// [[Rcpp::export]]
IntegerVector q_branch(NumericVector p, NumericVector mu, NumericVector sigma, NumericVector nu) {
  int n = p.size(); IntegerVector o(n);
  for (int i = 0; i < n; ++i) o[i] = fqZANBI_scalar(p[i], mu[i], sigma[i], nu[i], true, false);
  return o;
}
// [[Rcpp::export]]
NumericVector p_main(IntegerVector q, NumericVector mu, NumericVector sigma, NumericVector nu) {
  int n = q.size(); NumericVector o(n);
  for (int i = 0; i < n; ++i) o[i] = main_pZANBI_scalar(q[i], mu[i], sigma[i], nu[i]);
  return o;
}
// [[Rcpp::export]]
NumericVector p_branch(IntegerVector q, NumericVector mu, NumericVector sigma, NumericVector nu) {
  int n = q.size(); NumericVector o(n);
  for (int i = 0; i < n; ++i) o[i] = fpZANBI_scalar(q[i], mu[i], sigma[i], nu[i], true, false);
  return o;
}
