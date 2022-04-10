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


// dBNB ----
double my_dBNB_scalar(const int& x,
                      const double& mu = 1.0,
                      const double& sigma = 1.0,
                      const double& nu = 1.0,
                      const bool& log = false) {
    // if (mu    <= 0.0) stop("mu must be greater than 0");
    // if (sigma <= 0.0) stop("sigma must be greater than 0");
    // if (nu    <= 0.0) stop("nu must be greater than 0");
    // if (x      < 0.0) stop("x must be >=0");
    double m = (1.0/sigma) + 1.0;
    double n = (mu * nu)/sigma;
    double k = 1.0/nu;
    double logL = R::lbeta(x + n, m + k) - R::lbeta(n, m) - R::lgammafn(x + 1) -
     R::lgammafn(k) + R::lgammafn(x + k);
    double lik = logL;
    if (!log) lik = exp(logL);
    return lik;
}


//' @export
// [[Rcpp::export]]
NumericVector my_dBNB(const NumericVector& x,
                        const NumericVector& mu,
                        const NumericVector& sigma,
                        const NumericVector& nu,
                        const bool& log = false,
                        const int& n_cpu = 1)
  {
  if (x.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != nu.length()) stop("Distribution parameters must be of same length");
  const int n = x.length();
  for (int i = 0; i < n; i++)
  {
    if (x[i] < 0.0) stop("x must be >=0");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    if (nu[i] <= 0.0) stop("nu must be greater than 0");
  }

 NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_dBNB_scalar(x[i], mu[i], sigma[i], nu[i], log);
  }

  return out;
}

// pBNB ----
double my_pBNB_scalar(const int& q,
                      const double& mu = 1.0,
                      const double& sigma = 1.0,
                      const double& nu = 1.0,
                      const bool& lower_tail = true,
                      const bool& log_p = false)
  {
    // if (mu    <= 0.0) stop("mu must be greater than 0");
    // if (sigma <= 0.0) stop("sigma must be greater than 0");
    // if (nu    <= 0.0) stop("nu must be greater than 0");
    // if (q      < 0) stop("q must be >=0");

    double cdf = 0.0;
    for(int i = 0; i <= q; i++)
    {
      cdf += my_dBNB_scalar(i, mu, sigma, nu);
    }
    if (!lower_tail) cdf = 1.0 - cdf;
    if (log_p) cdf = log(cdf);
    return cdf;
  }

//' @export
// [[Rcpp::export]]
NumericVector my_pBNB(const IntegerVector& q,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const NumericVector& nu,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& n_cpu = 1)
  {
  if (q.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != nu.length()) stop("Distribution parameters must be of same length");
  const int n = q.length();
  for (int i = 0; i < n; i++)
  {
    if (q[i] < 0) stop("q must be >=0");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    if (nu[i] <= 0.0) stop("nu must be greater than 0");
  }

  NumericVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
  #pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_pBNB_scalar(q[i], mu[i], sigma[i], nu[i], lower_tail, log_p);
  }

    return out;
  }

// qBNB ----
// slow
int my_qBNB_scalar2(const double& p,
                   const double& mu = 1.0,
                   const double& sigma = 1.0,
                   const double& nu = 1.0,
                   const bool& lower_tail = true,
                   const bool& log_p = false,
                   const int& max_value = 10000)
{
  // if (mu    <= 0) stop("mu must be greater than 0");
  // if (sigma <= 0) stop("sigma must be greater than 0");
  // if (nu    <= 0) stop("nu must be greater than 0");
  // if (p < 0.0 || p > 1.0001) stop("p must be >=0 and <=1"); //I don't like this but it comes from original function

  double p_ = p;
  if (log_p) p_ = exp(p_);
  if (!lower_tail) p_ = 1.0 - p_;

  int QQQ = 0;
  double cumpro = 0.0;
  if (p_ + 1e-09 >= 1.0) QQQ = R_PosInf;
  else
    {
    for(int j = 0; j <= max_value; j++)
      {
      cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);
      QQQ = j;
      if (p_ <= cumpro) break;
      }
    }
  return QQQ;
}

// fast
int my_qBNB_scalar(const double& p,
                    const double& mu = 1.0,
                    const double& sigma = 1.0,
                    const double& nu = 1.0,
                    const bool& lower_tail = true,
                    const bool& log_p = false)
{
  // if (mu    <= 0) stop("mu must be greater than 0");
  // if (sigma <= 0) stop("sigma must be greater than 0");
  // if (nu    <= 0) stop("nu must be greater than 0");
  // if (p < 0.0 || p > 1.0001) stop("p must be >=0 and <=1"); //I don't like this but it comes from original function

  double p_ = p;
  if (log_p) p_ = exp(p_);
  if (!lower_tail) p_ = 1.0 - p_;

  int QQQ = 0;
  double cumpro = 0.0;
  int j = 0;

  if (p_ + 1e-09 >= 1.0) QQQ = R_PosInf;
  else
  {
    // implement a divide & conquer algo
    j = mu * (p + 0.5) + 1; // guess a starting value
    cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);
    if (p <= cumpro)
    {
      while (j > 0 && p <= cumpro)
      {
        j /= 2;
        cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);
      }
      for (int k = j; k <= (j * 2 + 1); k++)
      {
        cumpro = my_pBNB_scalar(k, mu, sigma, nu, true, false);
        if (p <= cumpro)
        {
          QQQ = k;
          break;
        }
      }
    }
    else // if p > cumpro
    {
      while (j < INT_MAX && p > cumpro)
      {
        j *= 2;
        cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);
      }
      j /= 2;
      cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);
      if ((j*2) > (j+1000))
      {
        while (j < INT_MAX && p > cumpro)
        {
          j += 1000;
          cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);
        }
        if (j >= 1000) j -= 1000;
        cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);
      }
      while (j < INT_MAX && p > cumpro)
      {
        j += 100;
        cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);
      }
      if (j >= 100) j -= 100;
      cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);

      while (j < INT_MAX && p > cumpro)
      {
        j += 10;
        cumpro = my_pBNB_scalar(j, mu, sigma, nu, true, false);
      }
      if (j >= 10) j -= 10;
      for (int k = j; k <= INT_MAX; k++)
      {
        cumpro = my_pBNB_scalar(k, mu, sigma, nu, true, false);
        if (p <= cumpro)
        {
          QQQ = k;
          break;
        }
      }
    }
  }
  return QQQ;
}


IntegerVector my_qBNB2(const NumericVector& p,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const NumericVector& nu,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& max_value = 10000,
                      const int& n_cpu = 1)
{
  if (p.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != nu.length()) stop("Distribution parameters must be of same length");
  const int n = p.length();
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0001) stop("p must be >=0 and <=1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    if (nu[i] <= 0.0) stop("nu must be greater than 0");
  }

  IntegerVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qBNB_scalar2(p[i], mu[i], sigma[i], nu[i], lower_tail, log_p, max_value);
  }

  return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector my_qBNB(const NumericVector& p,
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
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0001) stop("p must be >=0 and <=1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    if (nu[i] <= 0.0) stop("nu must be greater than 0");
  }

  IntegerVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qBNB_scalar(p[i], mu[i], sigma[i], nu[i], lower_tail, log_p);
  }

  return out;
}


// qZIBNB ----
int my_qZIBNB_scalar(const double& p,
                     const double& mu = 1.0,
                     const double& sigma = 1.0,
                     const double& nu = 1.0,
                     const double& tau = 0.1,
                     const bool& lower_tail = true,
                     const bool& log_p = false)
{
  // if (mu    <= 0) stop("mu must be greater than 0");
  // if (sigma <= 0) stop("sigma must be greater than 0");
  // if (nu    <= 0) stop("nu must be greater than 0");
  // if (tau <= 0.0 || tau >= 1.0) stop("tau must be >0 and <1");
  // if (p < 0.0 || p > 1.0001) stop("p must be >=0 and <=1"); //I don't like this but it comes from original function

  double p_ = p;
  if (log_p) p_ = exp(p_);
  if (!lower_tail) p_ = 1.0 - p_;

  p_ = (p_ - tau)/(1.0 - tau) - (1e-07);
  if (p_ <= 0) p_ = 0.0;
  return my_qBNB_scalar(p_, mu, sigma, nu, true, false);
}


//' @export
// [[Rcpp::export]]
IntegerVector my_qZIBNB(const NumericVector& p,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const NumericVector& nu,
                      const NumericVector& tau,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& n_cpu = 1)
{
  if (p.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != nu.length()) stop("Distribution parameters must be of same length");
  if (nu.length() != tau.length()) stop("Distribution parameters must be of same length");
  const int n = p.length();
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0001) stop("p must be >=0 and <=1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    if (nu[i] <= 0.0) stop("nu must be greater than 0");
    if (tau[i] <= 0.0 || tau[i] >= 1.0) stop("tau must be >0 and <1");
  }

  IntegerVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qZIBNB_scalar(p[i], mu[i], sigma[i], nu[i], tau[i], lower_tail, log_p);
  }

  return out;
}

// qZABNB ----
int my_qZABNB_scalar(const double& p,
                     const double& mu = 1.0,
                     const double& sigma = 1.0,
                     const double& nu = 1.0,
                     const double& tau = 0.1,
                     const bool& lower_tail = true,
                     const bool& log_p = false)
{
  // if (mu    <= 0) stop("mu must be greater than 0");
  // if (sigma <= 0) stop("sigma must be greater than 0");
  // if (nu    <= 0) stop("nu must be greater than 0");
  // if (tau <= 0.0 || tau >= 1.0) stop("tau must be >0 and <1");
  // if (p < 0.0 || p > 1.0) stop("p must be >=0 and <=1"); //I don't like this but it comes from original function


  double p_ = p;
  if (log_p) p_ = exp(p_);
  if (!lower_tail) p_ = 1.0 - p_;

  p_ = (p_ - tau)/(1.0 - tau) - (1e-010);
  double cdf0 = my_pBNB_scalar(0, mu, sigma, nu);
  p_ = cdf0 * (1.0 - p_) + p_;
  if (p_ < 0.0) p_ = 0.0;

  return my_qBNB_scalar(p_, mu, sigma, nu, true, false);
}


//' @export
// [[Rcpp::export]]
IntegerVector my_qZABNB(const NumericVector& p,
                        const NumericVector& mu,
                        const NumericVector& sigma,
                        const NumericVector& nu,
                        const NumericVector& tau,
                        const bool& lower_tail = true,
                        const bool& log_p = false,
                        const int& n_cpu = 1)
{
  if (p.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != nu.length()) stop("Distribution parameters must be of same length");
  if (nu.length() != tau.length()) stop("Distribution parameters must be of same length");
  const int n = p.length();
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0) stop("p must be >=0 and <=1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    if (nu[i] <= 0.0) stop("nu must be greater than 0");
    if (tau[i] <= 0.0 || tau[i] >= 1.0) stop("tau must be >0 and <1");
  }

  IntegerVector out(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    out[i] = my_qZABNB_scalar(p[i], mu[i], sigma[i], nu[i], tau[i], lower_tail, log_p);
  }

  return out;
}


/*** R
# N <- 1e3
# n_cpu <- 1L
# x <- sample(1e4, N, TRUE)
# q <- sample(1e2, N, TRUE)
# p <- runif(N, 0, 0.999)
# mu <- runif(N, 0.01, 100)
# sigma <- runif(N)
# nu <- runif(N)
# tau <- runif(N)
# log <- FALSE
# lower_tail <- TRUE
# max_value <- 10000L
# library(microbenchmark)
# library(ggplot2)
# library(gamlss)

# all.equal(my_dBNB(x, mu, sigma, nu, log, n_cpu), dBNB(x, mu, sigma, nu, log))
# all.equal(my_pBNB(q, mu, sigma, nu, lower_tail, log, n_cpu), pBNB(q, mu, sigma, nu, lower_tail, log))
# all.equal(my_qBNB(p, mu, sigma, nu, lower_tail, log, n_cpu), qBNB(p, mu, sigma, nu, lower_tail, log, max_value))
# all.equal(my_qBNB(p, mu, sigma, nu, lower_tail, log, n_cpu), my_qBNB2(p, mu, sigma, nu, lower_tail, log, max_value, n_cpu))
# all.equal(my_qZIBNB(p, mu, sigma, nu, tau, lower_tail, log, n_cpu), qZIBNB(p, mu, sigma, nu, tau, lower_tail, log, max_value))
# all.equal(my_qZABNB(p, mu, sigma, nu, tau, lower_tail, log, n_cpu), qZABNB(p, mu, sigma, nu, tau, lower_tail, log, max_value))

# autoplot(microbenchmark(my_dBNB(x, mu, sigma, nu, log, n_cpu), dBNB(x, mu, sigma, nu, log), times = 10))
# autoplot(microbenchmark(my_pBNB(q, mu, sigma, nu, lower_tail, log, n_cpu), pBNB(q, mu, sigma, nu, lower_tail, log), times = 10))
# autoplot(microbenchmark(my_qBNB(p, mu, sigma, nu, lower_tail, log, n_cpu), qBNB(p, mu, sigma, nu, lower_tail, log, max_value), times = 10))
# autoplot(microbenchmark(my_qBNB(p, mu, sigma, nu, lower_tail, log, n_cpu), my_qBNB2(p, mu, sigma, nu, lower_tail, log, max_value, n_cpu), times = 10))
# autoplot(microbenchmark(my_qZIBNB(p, mu, sigma, nu, tau, lower_tail, log, n_cpu), qZIBNB(p, mu, sigma, nu, tau, lower_tail, log, max_value), times = 10))
# autoplot(microbenchmark(my_qZABNB(p, mu, sigma, nu, tau, lower_tail, log, n_cpu), qZABNB(p, mu, sigma, nu, tau, lower_tail, log, max_value), times = 10))


# my_qBNB_scalar(0.999, 99.64497, 26.9589, 12.45386, lower_tail, log) # 3040
# autoplot(microbenchmark(my_qZIBNB(p, mu, sigma, nu, tau, lower_tail, log, 1), my_qZIBNB(p, mu, sigma, nu, tau, lower_tail, log, 20), times = 10))


*/
