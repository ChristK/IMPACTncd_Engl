/* IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by Chris
 Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
 Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
 funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngl: Health Outcomes
 Research Simulation Environment.  The views expressed are those of the
 authors and not necessarily those of the NHS, the NIHR or the Department of
 Health.

 Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos

 IMPACTncdEngl is free software; you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the Free Software
 Foundation; either version 3 of the License, or (at your option) any later
 version. This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 details. You should have received a copy of the GNU General Public License
 along with this program; if not, see <http://www.gnu.org/licenses/> or write
 to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 Boston, MA 02110-1301 USA. */

#include <Rcpp.h>
#include <math.h>
#include <Rmath.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#include <omp.h>

using namespace Rcpp;
//dBCT

double my_dBCCG_scalar(const double& x,
                       const double& mu = 1.0,
                       const double& sigma = 0.1,
                       const double& nu = 1.0,
                       const bool& log_ = false)
{
  // const double pi = 3.14159265358979323846264338328;
  const double pi = 4*atan(1); // calculate pi const
  double z;
  if (nu != 0)  z = (std::pow(x/mu, nu) - 1.0)/(nu * sigma);
    else z = log(x/mu)/sigma;
  double loglik = nu * log(x/mu) - log(sigma) - (z * z)/2.0 - log(x) -
        (log(2.0 * pi))/2.0 - log(R::pnorm(1.0/(sigma * std::abs(nu)), 0.0, 1.0, 1, 0));
  if (!log_) loglik = exp(loglik);
  return loglik;
}


//' @export
// [[Rcpp::export]]
NumericVector my_dBCT(const NumericVector& x,
                        const NumericVector& mu,
                        const NumericVector& sigma,
                        const NumericVector& nu,
                        const NumericVector& tau,
                        const bool& log_ = false,
                        const int& n_cpu = 1)
{
  if (x.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length() || nu.length() != tau.length()
  ) stop("Distribution parameters must be of same length");
  // TODO recycle if unequal lengths

  const int n = x.length();
  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (x[i] < 0.0) stop("x must be >=0");
    if (mu[i] <= 0.0) stop("mu must be positive");
    if (sigma[i] <= 0.0) stop("sigma must be positive");
    if (tau[i] <= 0.0 ) stop("tau must be positive");

    double z;
    if (nu[i] != 0) z = (std::pow(x[i]/mu[i], nu[i])-1.0)/(nu[i]*sigma[i]);
    else z = log(x[i]/mu[i])/sigma[i];

    out[i] = (nu[i]-1.0)*log(x[i])-nu[i]*log(mu[i])-log(sigma[i]);
    double fTz = R::lgammafn((tau[i]+1.0)/2.0)-R::lgammafn(tau[i]/2.0)-0.5*log(tau[i])-R::lgammafn(0.5) -
      ((tau[i]+1.0)/2.0) * log(1.0+(z*z)/tau[i]);
    out[i] += fTz-log(R::pt(1.0/(sigma[i]*std::abs(nu[i])),tau[i], 1, 0));
    if (tau[i]>1000000) out[i] = my_dBCCG_scalar(x[i],mu[i],sigma[i],nu[i],true);
  }

  if (!log_) out = exp(out);

  if (any(is_na(out))) warning("NaNs or NAs were produced");
  return out;

}




//pBCT
// double my_pBCT_scalar(const double& q,
//                       const double& mu,
//                       const double& sigma,
//                       const double& nu,
//                       const double& tau)
// {
//   // TODO test input within accepted range
//   double z;
//   if (nu != 0) z = ((std::pow(q/mu, nu) - 1)/(nu * sigma));
//   else z = log(q/mu)/sigma;
//   double FYy1 = R::pt(z, tau, true, false);
//   double FYy2 = 0;
//   if (nu > 0) FYy2 = R::pt(-1/(sigma * abs(nu)), tau, true, false);
//   double FYy3 = R::pt(1/(sigma * abs(nu)), tau, true, false);
//   return (FYy1 - FYy2)/FYy3;
// }



//' @export
// [[Rcpp::export]]
NumericVector my_pBCT(const NumericVector& q,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const NumericVector& nu,
                      const NumericVector& tau,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& n_cpu = 1)
{
  if (q.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length() || nu.length() != tau.length()
  ) stop("Distribution parameters must be of same length");

  const int n = q.length();
  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (q[i] < 0.0) stop("x must be >=0");
    if (mu[i] < 0.0) stop("mu must be positive");
    if (sigma[i] < 0.0) stop("sigma must be positive");
    if (tau[i] < 0.0 ) stop("tau must be positive");

    double z;
    if (nu[i] != 0) z = (std::pow(q[i]/mu[i], nu[i]) - 1.0)/(nu[i] * sigma[i]);
    else z = log(q[i]/mu[i])/sigma[i];
    double FYy1 = R::pt(z, tau[i], true, false);
    double FYy2 = 0.0;
    if (nu[i] > 0) FYy2 = R::pt(-1.0/(sigma[i] * std::abs(nu[i])), tau[i], true, false);
    double FYy3 = R::pt(1.0/(sigma[i] * std::abs(nu[i])), tau[i], true, false);
    out[i] = (FYy1 - FYy2)/FYy3;
  }
  if (!lower_tail) out = 1.0 - out;
  if (log_p) out = log(out);
  if (any(is_na(out))) warning("NaNs or NAs were produced");
  return out;
}

// qBCT
// double my_qBCT_scalar(const double& p,
//                       const double& mu,
//                       const double& sigma,
//                       const double& nu,
//                       const double& tau)
// {
//   // TODO test input within accepted range
//   double z;
//   if (nu <= 0) z = R::qt(p * R::pt(1/(sigma * abs(nu)), tau, true, false), tau, true, false);
//   else z = R::qt(1 - (1 - p) * R::pt(1/(sigma * abs(nu)), tau, true, false), tau, true, false);
//   double ya;
//   if (nu != 0) ya = mu * pow((nu * sigma * z + 1), (1/nu));
//   else ya = mu * exp(sigma * z);
//   return ya;
// }


//' @export
// [[Rcpp::export]]
NumericVector my_qBCT(      NumericVector  p,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const NumericVector& nu,
                      const NumericVector& tau,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& n_cpu = 1)
{

  if (p.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length() || nu.length() != tau.length()
  ) stop("Distribution parameters must be of same length");

  const int n = p.length();
  NumericVector out(n);
  if (log_p) p = exp(p);
  if(!lower_tail) p = 1.0 - p;

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (p[i] <= 0.0 || p[i] >= 1.0) stop("p must be between 0 and 1");
    if (mu[i] < 0.0) stop("mu must be positive");
    if (sigma[i] < 0.0) stop("sigma must be positive");
    if (tau[i] < 0.0 ) stop("tau must be positive");

    double z;
    if (nu[i] <= 0.0) z = R::qt(p[i] * R::pt(1.0/(sigma[i] * std::abs(nu[i])), tau[i], true, false), tau[i], true, false);
    else z = R::qt(1.0 - (1.0 - p[i]) * R::pt(1.0/(sigma[i] * std::abs(nu[i])), tau[i], true, false), tau[i], true, false);
    if (nu[i] != 0) out[i] = mu[i] * pow((nu[i] * sigma[i] * z + 1.0), (1.0/nu[i]));
    else out[i] = mu[i] * exp(sigma[i] * z);
  }
  if (any(is_na(out))) warning("NaNs or NAs were produced");
  return out;
}

// double my_qBCT_trunc_scalar(const double& p,
//                             const double& mu,
//                             const double& sigma,
//                             const double& nu,
//                             const double& tau,
//                             const double& lower_lim,
//                             const double& upper_lim) {
//   double pp1 = my_pBCT_scalar(lower_lim, mu, sigma, nu, tau);
//   double pp2 = my_pBCT_scalar(upper_lim, mu, sigma, nu, tau);
//   return my_qBCT_scalar((p * (pp2 - pp1) + pp1),  mu, sigma, nu, tau);
// }
//
// //' @export
// // [[Rcpp::export]]
// NumericVector my_qBCT_trunc(const NumericVector& p,
//                             const NumericVector& mu,
//                             const NumericVector& sigma,
//                             const NumericVector& nu,
//                             const NumericVector& tau,
//                             const double& lower_lim,
//                             const double& upper_lim,
//                             const int& n_cpu) {
//   // TODO check all have same length and are > 0 etc.
//   const int n = p.length();
//   NumericVector out(n);
//
//   omp_set_num_threads(n_cpu); // Use n_cpu threads for all
//
//   // consecutive parallel regions
// #pragma omp parallel for default(shared)
//   for (int i = 0; i < n; i++)
//   {
//     out[i] = my_qBCT_trunc_scalar(p[i], mu[i], sigma[i], nu[i], tau[i], lower_lim, upper_lim);
//   }
//
//   return out;
// }


/*** R

# N <- 1e5
# n_cpu <- 1L
# x <- sample(1e2, N, TRUE)
# q <- sample(1e2, N, TRUE)
# p <- runif(N, 0, 1)
# mu <- runif(N, 0, 30)
# sigma <- runif(N, 0, 5)
# nu <- runif(N, 0, 1)
# tau <- runif(N, 1e5, 1e7)
# log <- FALSE
# lower_tail <- TRUE
# max_value <- 10000L
#
# library(microbenchmark)
# library(ggplot2)
# library(gamlss)

# all.equal(my_dBCT(x, mu, sigma, nu, tau, log), gamlss.dist::dBCT(x, mu, sigma, nu, tau, log))
# autoplot(microbenchmark(my_dBCT(x, mu, sigma, nu, tau, log), gamlss.dist::dBCT(x, mu, sigma, nu, tau, log), times = 10))

# all.equal(my_pBCT(q, mu, sigma, nu, tau, lower_tail, log), gamlss.dist::pBCT(q, mu, sigma, nu, tau, lower_tail, log))
# autoplot(microbenchmark(my_pBCT(q, mu, sigma, nu, tau, lower_tail, log), gamlss.dist::pBCT(q, mu, sigma, nu, tau, lower_tail, log), times = 10))

# all.equal(my_qBCT(p, mu, sigma, nu, tau, lower_tail, log), gamlss.dist::qBCT(p, mu, sigma, nu, tau, lower_tail, log))
# autoplot(microbenchmark(my_qBCT(p, mu, sigma, nu, tau, lower_tail, log, 20), gamlss.dist::qBCT(p, mu, sigma, nu, tau, lower_tail, log), times = 10))

*/
