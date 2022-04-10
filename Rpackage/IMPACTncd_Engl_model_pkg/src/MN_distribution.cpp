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
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#include <omp.h>

using namespace Rcpp;


// qMN4 ----
int my_qMN4_scalar(const double& p,
                   const double& mu = 1.0,
                   const double& sigma = 1.0,
                   const double& nu = 1.0,
                   const bool&   lower_tail = true,
                   const bool&   log_p = false)
{
  if (mu    <= 0) stop("mu must be greater than 0");
  if (sigma <= 0) stop("sigma must be greater than 0");
  if (nu    <= 0) stop("nu must be greater than 0");
  if (p < 0.0 || p > 1.0) stop("q must be >=0 and <=1");

  double p_ = p;
  if (log_p) p_ = exp(p); // not in R code
  if (!lower_tail) p_ = 1.0 - p; // not in R code

  int q = 1;
  double tt = 1.0 + mu + sigma + nu;

  if (p_ >= (mu/tt))                q = 2;
  if (p_ >= ((mu + sigma)/tt))      q = 3;
  if (p_ >= ((mu + sigma + nu)/tt)) q = 4;
  return q;
}


//' @export
// [[Rcpp::export]]
IntegerVector my_qMN4(const NumericVector& p,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const NumericVector& nu,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& n_cpu = 1)
{
  if (p.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != nu.length()) stop("Distribution parameters must be of same length");

  const int n = p.length();
  IntegerVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qMN4_scalar(p[i], mu[i], sigma[i], nu[i], lower_tail, log_p);
  }

  return out;
}


/*** R

# N <- 1e6
# n_cpu <- 1L
# x <- sample(1e4, N, TRUE)
# q <- sample(1e2, N, TRUE)
# p <- runif(N)
# mu <- runif(N)
# sigma <- runif(N)
# nu <- runif(N)
# tau <- runif(N)
# log <- FALSE
# lower_tail <- TRUE
# max_value <- 10000L
# library(microbenchmark)
# library(ggplot2)
# library(gamlss)
# all.equal(my_qMN4(p, mu, sigma, nu, lower_tail, log, n_cpu), qMN4(p, mu, sigma, nu, lower_tail, log))
# autoplot(microbenchmark(my_qMN4(p, mu, sigma, nu, lower_tail, log, n_cpu),
#                         qMN4(p, mu, sigma, nu, lower_tail, log), times = 10))

*/
