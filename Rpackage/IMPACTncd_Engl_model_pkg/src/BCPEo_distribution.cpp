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

// dBCPEo -----
double my_dBCPEo_hlp_f_T (const double& t, const double& tau,
                                 const bool& log_ = false)
{
  double log_c;
  double c;
  double log_lik;

  log_c = 0.5 * (-(2.0/tau) * log(2.0) + R::lgammafn(1/tau) - R::lgammafn(3/tau));
  c = exp(log_c);
  log_lik = log(tau) - log_c - (0.5 * std::pow(std::abs(t/c),tau)) -
    (1.0 + (1.0/tau)) * log(2.0) - R::lgammafn(1.0/tau);

  if (log_) return log_lik;
  else return exp(log_lik);
}

double my_dBCPEo_hlp_F_T (const double& t, const double& tau)
{
  double log_c;
  double c;
  double s;
  double F_s;
  double cdf;

    log_c = 0.5 * (-(2.0/tau) * log(2.0) + R::lgammafn(1.0/tau) - R::lgammafn(3.0/tau));
    c = exp(log_c);
    s = 0.5 * std::pow(std::abs(t/c), tau);
    F_s = R::pgamma(s, 1.0/tau, 1.0, 1, 0);
    cdf = 0.5 * (1.0 + F_s * R::sign(t));

      return cdf;
}

//' @export
// [[Rcpp::export]]
NumericVector my_dBCPEo(const NumericVector& x,
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
  NumericVector loglik(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (x[i] < 0.0) stop("x must be >=0");
    if (mu[i] < 0.0) stop("mu must be positive");
    if (sigma[i] < 0.0) stop("sigma must be positive");
    if (tau[i] < 0.0 ) stop("tau must be positive");

    double z;
    if (nu[i] != 0) z = ((std::pow(x[i]/mu[i], nu[i]) - 1.0)/(nu[i] * sigma[i]));
      else z = log(x[i]/mu[i])/sigma[i];

    double logfZ =  my_dBCPEo_hlp_f_T(z, tau[i], true) -
      log(my_dBCPEo_hlp_F_T(1.0/(sigma[i] * std::abs(nu[i])), tau[i]));
    double logder = (nu[i] - 1.0) * log(x[i]) - nu[i] * log(mu[i]) - log(sigma[i]);

    loglik[i] = logder + logfZ;
  }


  if (!log_) loglik = exp(loglik);

  if (any(is_na(loglik))) warning("NaNs or NAs were produced");
  return loglik;
}

// pBCPEo ------
double my_pBCPEo_hlp_F_T(const double& t, const double& tau) {
  double log_c = 0.5 * (-(2.0/tau) * log(2.0) + R::lgammafn(1.0/tau) - R::lgammafn(3.0/tau));
  double c = exp(log_c);
  double s = 0.5 * (std::pow((std::abs(t/c)), tau));
  double F_s = R::pgamma(s, 1/tau, 1.0, 1, 0);
  double cdf = 0.5 * (1.0 + F_s * R::sign(t));
  return cdf;
}



//' @export
// [[Rcpp::export]]
NumericVector my_pBCPEo(const NumericVector& q,
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
 // TODO recycle if unequal lengths


  const int n = q.length();
  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (q[i] < 0.0) stop("q must be positive");
    if (mu[i] < 0.0) stop("mu must be positive");
    if (sigma[i] < 0.0) stop("sigma must be positive");
    if (tau[i] < 0.0 ) stop("tau must be positive");

    double z;
    if (nu[i] == 0.0) z = log(q[i]/mu[i])/sigma[i];
    else z = ((std::pow((q[i]/mu[i]),nu[i]) - 1.0)/(nu[i] * sigma[i]));

    double FYy1 = my_pBCPEo_hlp_F_T(z, tau[i]);

    double FYy2 = 0.0;
    if (nu[i] > 0.0)
    {
      FYy2 = my_pBCPEo_hlp_F_T(-1.0/(sigma[i] * std::abs(nu[i])), tau[i]);
    }
    double FYy3 = my_pBCPEo_hlp_F_T(1.0/(sigma[i] * std::abs(nu[i])), tau[i]);
    out[i] = (FYy1 - FYy2)/FYy3;
  }

  if (!lower_tail) out = 1.0 - out;
  if (log_p) out = log(out);

  if (any(is_na(out))) warning("NaNs or NAs were produced");

  return out;
}

// qBCPEo ------

double my_qBCPEo_hlp_q_T(const double& p, const double& tau) {
  double log_c = 0.5 * (-(2.0/tau) * log(2.0) + std::lgamma(1.0/tau) - std::lgamma(3.0/tau));
  double c = exp(log_c);
  double s = R::qgamma((2.0 * p - 1.0) * R::sign(p - 0.5), (1.0/tau), 1.0, true, false);
  return R::sign(p - 0.5) * std::pow((2.0 * s),(1.0/tau)) * c;
}


// double my_qBCPEo_scalar(const double& p,
//                         const double& mu,
//                         const double& sigma,
//                         const double& nu,
//                         const double& tau,
//                         const bool& lower_tail = true,
//                         const bool& log_p = false)
// {
//   double za = 0.0;
//   if (nu < 0)
//   {
//     za = my_qBCPEo_hlp_q_T(p * my_pBCPEo_hlp_F_T(1.0/(sigma * std::abs(nu)), tau), tau);
//   } else if (nu == 0)
//   {
//     za = my_qBCPEo_hlp_q_T(p, tau);
//   } else // if nu > 0
//   {
//     za = my_qBCPEo_hlp_q_T((1.0 - (1.0 - p) * my_pBCPEo_hlp_F_T(1.0/(sigma * std::abs(nu)), tau)), tau);
//   }
//   double ya = 0.0;
//   if (nu == 0.0) ya = mu * exp(sigma * za); else ya = mu * std::pow((nu * sigma * za + 1.0),(1.0/nu));
//
//   return ya;
// }



//' @export
// [[Rcpp::export]]
NumericVector my_qBCPEo(      NumericVector  p,
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
  // TODO recycle if unequal lengths

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

    double za;
    if (nu[i] < 0.0)
    {
      za = my_qBCPEo_hlp_q_T(p[i] * my_pBCPEo_hlp_F_T(1.0/(sigma[i] * std::abs(nu[i])), tau[i]), tau[i]);
    } else if (nu[i] == 0)
    {
      za = my_qBCPEo_hlp_q_T(p[i], tau[i]);
    } else // if nu > 0
    {
      za = my_qBCPEo_hlp_q_T((1.0 - (1.0 - p[i]) * my_pBCPEo_hlp_F_T(1.0/(sigma[i] * std::abs(nu[i])), tau[i])), tau[i]);
    }

    if (nu[i] == 0.0) out[i] = mu[i] * exp(sigma[i] * za);
    else out[i] = mu[i] * std::pow((nu[i] * sigma[i] * za + 1.0),(1.0/nu[i]));
    }
  if (any(is_na(out))) warning("NaNs or NAs were produced");
  return out;
}


// double my_qBCPEo_trunc_scalar(const double& p,
//                               const double& mu,
//                               const double& sigma,
//                               const double& nu,
//                               const double& tau,
//                               const double& lower_lim,
//                               const double& upper_lim) {
//   double pp1 = my_pBCPEo_scalar(lower_lim, mu, sigma, nu, tau);
//   double pp2 = my_pBCPEo_scalar(upper_lim, mu, sigma, nu, tau);
//   return my_qBCPEo_scalar((p * (pp2 - pp1) + pp1),  mu, sigma, nu, tau);
// }

// //' @export
// // [[Rcpp::export]]
// NumericVector my_qBCPEo_trunc(const NumericVector& p,
//                               const NumericVector& mu,
//                               const NumericVector& sigma,
//                               const NumericVector& nu,
//                               const NumericVector& tau,
//                               const double& lower_lim,
//                               const double& upper_lim,
//                               const int& n_cpu) {
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
//     out[i] = my_qBCPEo_trunc_scalar(p[i], mu[i], sigma[i], nu[i], tau[i], lower_lim, upper_lim);
//   }
//
//   return out;
// }
//

/*** R

# N <- 1e5
# n_cpu <- 1L
# x <- sample(1e2, N, TRUE)
# q <- sample(1e2, N, TRUE)
# p <- runif(N, 0, 1)
# mu <- runif(N, 0, 30)
# sigma <- runif(N, 0, 5)
# nu <- runif(N, 0, 1)
# tau <- runif(N, 0, 10)
# log <- FALSE
# lower_tail <- TRUE
# max_value <- 10000L
#
# library(microbenchmark)
# library(ggplot2)
# library(gamlss)

# all.equal(my_dBCPEo(x, mu, sigma, nu, tau, log), gamlss.dist::dBCPEo(x, mu, sigma, nu, tau, log))
# autoplot(microbenchmark(my_dBCPEo(x, mu, sigma, nu, tau, log), gamlss.dist::dBCPEo(x, mu, sigma, nu, tau, log), times = 10))

# all.equal(my_pBCPEo(q, mu, sigma, nu, tau, lower_tail, log), gamlss.dist::pBCPEo(q, mu, sigma, nu, tau, lower_tail, log))
# autoplot(microbenchmark(my_pBCPEo(q, mu, sigma, nu, tau, lower_tail, log), gamlss.dist::pBCPEo(q, mu, sigma, nu, tau, lower_tail, log), times = 10))

# all.equal(my_qBCPEo(p, mu, sigma, nu, tau, lower_tail, log), gamlss.dist::qBCPEo(p, mu, sigma, nu, tau, lower_tail, log))
# autoplot(microbenchmark(my_qBCPEo(p, mu, sigma, nu, tau, lower_tail, log), gamlss.dist::qBCPEo(p, mu, sigma, nu, tau, lower_tail, log), times = 10))



*/
