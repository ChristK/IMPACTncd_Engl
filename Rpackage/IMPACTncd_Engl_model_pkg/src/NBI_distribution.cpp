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
#include "NBI_distribution.h"

using namespace Rcpp;


// dNBI ----
double my_dNBI_scalar(const int& x,
                      const double& mu,
                      const double& sigma,
                      const bool& log_p,
                      const bool& check) {

  if (check)
  {
    if (mu    <= 0.0) stop("mu must be greater than 0");
    if (sigma <= 0.0) stop("sigma must be greater than 0");
    if (x      < 0.0) stop("x must be >=0");
  }

  double fy =  (sigma < 1e-4) ?
  R::dpois(x, mu, log_p) :
  R::dnbinom_mu(x, 1/sigma, mu, log_p);

  return fy;
}


//' @export
// [[Rcpp::export]]
NumericVector my_dNBI(const IntegerVector& x,
                        const NumericVector& mu,
                        const NumericVector& sigma,
                        const bool& log_p,
                        const int& n_cpu)
  {
  if (x.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  const int n = x.length();
  for (int i = 0; i < n; i++)
  {
    if (x[i] < 0.0) stop("x must be >=0");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
  }

 NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_dNBI_scalar(x[i], mu[i], sigma[i], log_p, false);
  }

  return out;
}

// pNBI ----
double my_pNBI_scalar(const int& q,
                      const double& mu,
                      const double& sigma,
                      const bool& lower_tail,
                      const bool& log_p,
                      const bool& check)
  {

  if (check)
  {
    if (mu    <= 0.0) stop("mu must be greater than 0");
    if (sigma <= 0.0) stop("sigma must be greater than 0");
    if (q      < 0) stop("q must be >=0");
  }

    double cdf = (sigma < 1e-4) ?
    R::ppois(q, mu, lower_tail, log_p) :
    R::pnbinom_mu(q, 1/sigma, mu, lower_tail, log_p);

    return cdf;
  }



//' @export
// [[Rcpp::export]]
NumericVector my_pNBI(const IntegerVector& q,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const bool& lower_tail,
                      const bool& log_p,
                      const int& n_cpu)
  {
  if (q.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  const int n = q.length();
  for (int i = 0; i < n; i++)
  {
    if (q[i] < 0) stop("q must be >=0");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
  }

  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
  #pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_pNBI_scalar(q[i], mu[i], sigma[i], lower_tail, log_p, false);
  }

    return out;
  }

// qNBI ----

int my_qNBI_scalar(const double& p,
                    const double& mu,
                    const double& sigma,
                    const bool& lower_tail,
                    const bool& log_p,
                    const bool& check)
{
  if (check)
  {
    if (mu    <= 0) stop("mu must be greater than 0");
    if (sigma <= 0) stop("sigma must be greater than 0");
    if (p < 0.0 || p > 1.0) stop("p must be >=0 and <=1");
  }

  int q = (sigma < 1e-4) ?
  R::qpois(p, mu, lower_tail, log_p) :
  R::qnbinom_mu(p, 1/sigma, mu, lower_tail, log_p);

  return q;
}

//' @export
// [[Rcpp::export]]
IntegerVector my_qNBI(const NumericVector& p,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const bool& lower_tail,
                      const bool& log_p,
                      const int& n_cpu)
{
  if (p.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  const int n = p.length();
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0) stop("p must be >=0 and <=1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
  }

  IntegerVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qNBI_scalar(p[i], mu[i], sigma[i], lower_tail, log_p, false);
  }

  return out;
}

// qZANBI ----

int my_qZANBI_scalar(const double& p,
                     const double& mu,
                     const double& sigma,
                     const double& nu,
                     const bool& lower_tail,
                     const bool& log_p,
                     const bool& check)
{
  if (check)
  {
    if (p < 0.0 || p > 1.0) stop("p must be >=0 and <=1");
    if (mu    <= 0) stop("mu must be greater than 0");
    if (sigma <= 0) stop("sigma must be greater than 0");
    if (nu    <= 0 || nu > 1.0) stop("nu must be >=0 and <=1");
  }

  double p_ = p;
  if (log_p) p_ = exp(p_);
  if (!lower_tail) p_ = 1.0 - p_;

  p_ = (p_ - nu)/(1 - nu) - 1e-10;
  double cdf0 = my_pNBI_scalar(0.0, mu, sigma, true, false, false);
  p_ = cdf0 * (1 - p_) + p_;
  if (p_ < 0.0) p_ = 0.0;

  return my_qNBI_scalar(p_, mu, sigma, true, false, false);
}


//' @export
// [[Rcpp::export]]
IntegerVector my_qZANBI(const NumericVector& p,
                        const NumericVector& mu,
                        const NumericVector& sigma,
                        const NumericVector& nu,
                        const bool& lower_tail,
                        const bool& log_p,
                        const int& n_cpu)
{
  if (p.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != nu.length()) stop("Distribution parameters must be of same length");
  const int n = p.length();
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0) stop("p must be >=0 and <=1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    if (nu[i] <= 0.0 || nu[i] > 1.0) stop("nu must be >=0 and <=1");
  }

  if (n_cpu == 1)
  {
    IntegerVector out(n);
    for (int i = 0; i < n; i++)
    {
      out[i] = my_qZANBI_scalar(p[i], mu[i], sigma[i], nu[i], lower_tail, log_p, false);
    }
    return out;
  }
  else // if parallel
  {
    std::vector<int> out(n);
    std::vector<int> mu_(mu.begin(), mu.end());
    std::vector<int> sigma_(sigma.begin(), sigma.end());
    std::vector<int> nu_(nu.begin(), nu.end());

    omp_set_num_threads(n_cpu); // Use n_cpu threads for all

    // consecutive parallel regions
#pragma omp parallel for simd default(shared)
    for (int i = 0; i < n; i++)
    {
      out[i] = my_qZANBI_scalar(p[i], mu_[i], sigma_[i], nu_[i], lower_tail, log_p, false);
    }

    //
    return Rcpp::IntegerVector(out.begin(), out.end());
  }
}


// pZANBI ----
//' @export
// [[Rcpp::export]]
double my_pZANBI_scalar(const int& q,
                     const double& mu,
                     const double& sigma,
                     const double& nu,
                     const bool& lower_tail,
                     const bool& log_p,
                     const bool& check)
{
  if (check)
  {
    if (q < 0) stop("q must be >=0 and <=1");
    if (mu    <= 0.0) stop("mu must be greater than 0");
    if (sigma <= 0.0) stop("sigma must be greater than 0");
    if (nu    <= 0.0 || nu > 1.0) stop("nu must be >=0 and <=1");
  }

  double cdf0 = my_pNBI_scalar(0, mu, sigma);
  double cdf1 = my_pNBI_scalar(q, mu, sigma);
  double cdf3 = nu + ((1.0 - nu) * (cdf1 - cdf0)/(1 - cdf0));
  double cdf  = q == 0 ? nu : cdf3;
  if (!lower_tail) cdf = 1.0 - cdf;
  if (log_p) cdf = log(cdf);

  return cdf;
}


/*** R
# set.seed(42L)
# N <- 1e5
# n_cpu <- 1L
# x <- sample(1e4, N, TRUE)
# q <- sample(1e2, N, TRUE)
# p <- runif(N, 0, 0.999)
# mu <- runif(N, 0.01, 100)
# sigma <- c(runif(N/2, 1e-9, 1e-4), runif(N/2, 1e-4, 1e4))
# table(sigma < 1e-4)
# nu <- runif(N)
# log_p <- FALSE
# lower_tail <- TRUE
# library(microbenchmark)
# library(ggplot2)
# library(gamlss)

# my_qZANBI_scalar(p[1], mu[1], sigma[1], nu[1], lower_tail, log_p)
# all.equal(my_dNBI(x, mu, sigma, log_p, n_cpu), dNBI(x, mu, sigma, log_p))
# all.equal(my_pNBI(q, mu, sigma, lower_tail, log_p, n_cpu), pNBI(q, mu, sigma, lower_tail, log_p))
# all.equal(my_qNBI(p, mu, sigma, lower_tail, log_p, n_cpu), qNBI(p, mu, sigma, lower_tail, log_p))
# all.equal(my_qZANBI(p, mu, sigma, nu, lower_tail, log_p, n_cpu), qZANBI(p, mu, sigma, nu, lower_tail, log_p))

# autoplot(microbenchmark(my_dNBI(x, mu, sigma, log_p, n_cpu), dNBI(x, mu, sigma, log_p), times = 10))
# autoplot(microbenchmark(my_pNBI(q, mu, sigma, lower_tail, log_p, n_cpu), pNBI(q, mu, sigma, lower_tail, log_p), times = 10))
# autoplot(microbenchmark(my_qNBI(p, mu, sigma, lower_tail, log_p, n_cpu), qNBI(p, mu, sigma, lower_tail, log_p), times = 10))
# autoplot(microbenchmark(my_qZANBI(p, mu, sigma, nu, lower_tail, log_p, n_cpu), qZANBI(p, mu, sigma, nu, lower_tail, log_p), times = 10))

# n_cpu <- 10L
# autoplot(microbenchmark(my_dNBI(x, mu, sigma, log_p, n_cpu), dNBI(x, mu, sigma, log_p), times = 10))
# autoplot(microbenchmark(my_pNBI(q, mu, sigma, lower_tail, log_p, n_cpu), pNBI(q, mu, sigma, lower_tail, log_p), times = 10))
# autoplot(microbenchmark(my_qNBI(p, mu, sigma, lower_tail, log_p, n_cpu), qNBI(p, mu, sigma, lower_tail, log_p), times = 10))
# autoplot(microbenchmark(my_qZANBI(p, mu, sigma, nu, lower_tail, log_p, n_cpu), qZANBI(p, mu, sigma, nu, lower_tail, log_p), times = 10))


*/
