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
#include <algorithm>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
using namespace Rcpp;

// dPO (necessary for dDL)
double my_dPO_scalar(const int &x, const double &mu = 1.0, const double &sigma = 1.0, const bool &log_ = false)
{
  // if (x < 0) stop("x must be >=0");
  // if (mu <= 0.0) stop("mu must be greater than 0");
  // if (sigma <= 0.0) stop("sigma must be greater than 0");
  double fy;
  if (sigma > 1e-04)
  {
    fy = R::dnbinom_mu(x, 1.0 / sigma, mu, log_);
  }
  else
  {
    fy = R::dpois(x, mu, log_);
  }
  return fy;
}

// from gamlss.dist
NumericVector my_tofydel2(const NumericVector &y, const NumericVector &mu,
                          const NumericVector &sigma,
                          const NumericVector &nu)
{
  int maxyp1 = max(y) + 2;
  int ny = y.length();
  NumericVector ans(ny);
  double tofY[maxyp1];
  int iy, i, j;
  double sumT, dum;
  for (i = 0; i < ny; i++)
  {
    iy = y[i] + 1;
    tofY[0] = mu[i] * nu[i] + mu[i] * (1 - nu[i]) / (1 + mu[i] * sigma[i] * (1 - nu[i]));
    sumT = 0;
    for (j = 1; j < iy; j++)
    {
      dum = 1 + 1 / (mu[i] * sigma[i] * (1 - nu[i]));
      tofY[j] = (j + mu[i] * nu[i] + 1 / (sigma[i] * (1 - nu[i])) - (mu[i] * nu[i] * j) / tofY[j - 1]) / dum;
      sumT += log(tofY[j - 1]);
    }
    ans[i] = sumT;
  }
  return ans;
}

double my_tofydel2_scalar(const int &y, const double &mu,
                          const double &sigma,
                          const double &nu)
{
  double tofY[y + 2];
  int iy, j;
  double sumT, dum;
  iy = y + 1;
  tofY[0] = mu * nu + mu * (1.0 - nu) / (1.0 + mu * sigma * (1.0 - nu));
  sumT = 0;
  for (j = 1; j < iy; j++)
  {
    dum = 1 + 1 / (mu * sigma * (1.0 - nu));
    tofY[j] = (j + mu * nu + 1.0 / (sigma * (1.0 - nu)) - (mu * nu * j) / tofY[j - 1]) / dum;
    sumT += log(tofY[j - 1]);
  }
  return sumT;
}

// dDEL

//' @export
// [[Rcpp::export]]
NumericVector my_dDEL(const IntegerVector &x,
                      const NumericVector &mu,
                      const NumericVector &sigma,
                      const NumericVector &nu,
                      const bool &log_ = false,
                      const int &n_cpu = 1)
{
  if (x.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length())
    stop("Distribution parameters must be of same length");
  // TODO recycle if unequal lengths

  const int n = x.length();
  NumericVector logfy(n);
  NumericVector logpy0(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all

// consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (x[i] < 0.0)
      stop("x must be >=0");
    if (mu[i] <= 0.0)
      stop("mu must be greater than 0");
    if (sigma[i] <= 0.0)
      stop("sigma must be greater than 0");
    if (nu[i] <= 0.0 || nu[i] >= 1.0)
      stop("nu must be between 0 and 1");

    if (sigma[i] < 1e-04)
      logfy[i] = R::dpois(x[i], mu[i], (int)log_);
    else
    {
      logpy0[i] = -mu[i] * nu[i] - (1.0 / sigma[i]) * (log(1.0 + mu[i] * sigma[i] *
                                                                     (1.0 - nu[i])));
      double S = my_tofydel2_scalar(x[i], mu[i], sigma[i], nu[i]);
      logfy[i] = logpy0[i] - R::lgammafn(x[i] + 1) + S;
      if (!log_)
        logfy[i] = exp(logfy[i]);
    }
  }

  // If nu is far away from 0 (i.e. 100 or -100) bessel_k produces very small numbers
  // that lead to NaNs in the output
  if (any(is_na(logfy)))
    warning("NaNs or NAs were produced");
  return logfy; // despite the name, not logged if log_ = false
}

// pDEL
double my_dDEL_scalar(const int &x,
                      const double &mu,
                      const double &sigma,
                      const double &nu,
                      const bool &log_ = false)
{
  double logfy = 0.0;
  if (sigma < 1e-04)
    logfy = R::dpois(x, mu, (int)log_);
  else
  {
    double logpy0 = -mu * nu - (1.0 / sigma) * (log(1.0 + mu * sigma *
                                                              (1.0 - nu)));
    double S = my_tofydel2_scalar(x, mu, sigma, nu);
    logfy = logpy0 - R::lgammafn(x + 1) + S;
    if (!log_)
      logfy = exp(logfy);
  }
  // If nu is far away from 0 (i.e. 100 or -100) bessel_k produces very small numbers
  // that lead to NaNs in the output
  // if (is_na(logfy)) warning("NaNs or NAs were produced");
  return logfy; // despite the name, not logged if log_ = false
}

double my_pDEL_hlp_fn(const int &q,
                      const double &mu,
                      const double &sigma,
                      const double &nu)
{
  double ans = 0.0;
  for (int i = 0; i <= q; i++)
  {
    ans += my_dDEL_scalar(i, mu, sigma, nu);
  }

  return ans;
}

//' @export
// [[Rcpp::export]]
NumericVector my_pDEL(const IntegerVector &q,
                      const NumericVector &mu,
                      const NumericVector &sigma,
                      const NumericVector &nu,
                      const bool &lower_tail = true,
                      const bool &log_p = false,
                      const int &n_cpu = 1)
{
  if (q.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length())
    stop("Distribution parameters must be of same length");
  // TODO recycle if unequal lengths

  const int n = q.length();
  NumericVector cdf(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all
                              // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (q[i] < 0)
      stop("q must be >=0");
    if (mu[i] <= 0.0)
      stop("mu must be greater than 0");
    if (sigma[i] <= 0.0)
      stop("sigma must be greater than 0");
    if (nu[i] <= 0.0 || nu[i] >= 1)
      stop("nu must be between 0 and 1");
    cdf[i] = my_pDEL_hlp_fn(q[i], mu[i], sigma[i], nu[i]);
  }

  if (!lower_tail)
    cdf = 1.0 - cdf;
  if (log_p)
    cdf = log(cdf);

  if (any(is_na(cdf)))
    warning("NaNs or NAs were produced");
  return cdf;
}

double my_pDEL_scalar(const int &q,
                      const double &mu,
                      const double &sigma,
                      const double &nu,
                      const bool &lower_tail = true,
                      const bool &log_p = false)
{
  double cdf = my_pDEL_hlp_fn(q, mu, sigma, nu);
  if (!lower_tail)
    cdf = 1.0 - cdf;
  if (log_p)
    cdf = log(cdf);

  return cdf;
}

// qDEL
// fast! and not restricted by max_value
//' @export
// [[Rcpp::export]]
IntegerVector my_qDEL(NumericVector p,
                      const NumericVector &mu,
                      const NumericVector &sigma,
                      const NumericVector &nu,
                      const bool &lower_tail = true,
                      const bool &log_p = false,
                      const int &n_cpu = 1)
{
  // TODO implement recycling with rep_len()
  if (p.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length())
    stop("Distribution parameters must be of same length");

  const int n = p.length();
  IntegerVector QQQ(n);
  double cumpro = 0.0;
  int j = 0;

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all
// consecutive parallel regions
#pragma omp parallel for private(cumpro, j) shared(QQQ)
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0001)
      stop("p must be between 0 and 1");
    if (mu[i] <= 0.0)
      stop("mu must be greater than 0");
    if (sigma[i] <= 0.0)
      stop("sigma must be greater than 0");
    if (nu[i] <= 0.0 || nu[i] >= 1)
      stop("nu must be between 0 and 1");

    if (log_p)
      p[i] = exp(p[i]);
    if (!lower_tail)
      p[i] = 1.0 - p[i];

    if (p[i] + 1e-09 >= 1.0)
      QQQ[i] = R_PosInf; // I don't like this
    else
    {
      // implement a divide & conquer algo
      j = mu[i] * (p[i] + 0.5) + 1; // guess a starting value
      cumpro = my_pDEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
      if (p[i] <= cumpro)
      {
        while (j > 0 && p[i] <= cumpro)
        {
          j /= 2;
          cumpro = my_pDEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        }
        for (int k = j; k <= (j * 2 + 1); k++)
        {
          cumpro = my_pDEL_scalar(k, mu[i], sigma[i], nu[i], true, false);
          if (p[i] <= cumpro)
          {
            QQQ[i] = k;
            break;
          }
        }
      }
      else // if p[i] > cumpro
      {
        while (j < INT_MAX && p[i] > cumpro)
        {
          j *= 2;
          cumpro = my_pDEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        }
        j /= 2;
        cumpro = my_pDEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        if ((j * 2) > (j + 1000))
        {
          while (j < INT_MAX && p[i] > cumpro)
          {
            j += 1000;
            cumpro = my_pDEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
          }
          if (j >= 1000)
            j -= 1000;
          cumpro = my_pDEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        }

        while (j < INT_MAX && p[i] > cumpro)
        {
          j += 100;
          cumpro = my_pDEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        }
        if (j >= 100)
          j -= 100;
        cumpro = my_pDEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        while (j < INT_MAX && p[i] > cumpro)
        {
          j += 10;
          cumpro = my_pDEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        }
        if (j >= 10)
          j -= 10;
        for (int k = j; k <= INT_MAX; k++)
        {
          cumpro = my_pDEL_scalar(k, mu[i], sigma[i], nu[i], true, false);
          if (p[i] <= cumpro)
          {
            QQQ[i] = k;
            break;
          }
        }
      }
    }
  }
  if (any(is_na(QQQ)))
    warning("NaNs or NAs were produced");
  return QQQ;
}

// # N <- 1e4
// # n_cpu <- 1L
// # x <- sample(1e2, N, TRUE)
// # q <- sample(1e2, N, TRUE)
// # p <- runif(N, 0, 0.99)
// # mu <- runif(N, 0.8, 6)
// # sigma <- runif(N, 0, 21)
// # sigma[1] <- 0.000001
// # nu <- runif(N, 0, 1)
// # log <- FALSE
// # lower_tail <- TRUE
// # max_value <- 10000L
// #
// # library(microbenchmark)
// # library(ggplot2)
// # library(gamlss)
//
// # all.equal(my_dDEL(x, mu, sigma, nu, log, n_cpu), gamlss.dist::dDEL(x, mu, sigma, nu, log))
// # all.equal(my_dDEL(x, mu, sigma, nu, !log, n_cpu), gamlss.dist::dDEL(x, mu, sigma, nu, !log))
// # all.equal(my_dDEL(x, mu, sigma, nu, log, 20), gamlss.dist::dDEL(x, mu, sigma, nu, log))
//
// # autoplot(microbenchmark(my_dDEL(x, mu, sigma, nu, log, 20), gamlss.dist::dDEL(x, mu, sigma, nu, log), times = 10))
//
// # all.equal(my_pDEL(q, mu, sigma, nu, lower_tail, log, n_cpu), gamlss.dist::pDEL(q, mu, sigma, nu, lower_tail, log))
// # all.equal(my_pDEL(q, mu, sigma, nu, !lower_tail, log, n_cpu), gamlss.dist::pDEL(q, mu, sigma, nu, !lower_tail, log))
// # all.equal(my_pDEL(q, mu, sigma, nu, lower_tail, !log, n_cpu), gamlss.dist::pDEL(q, mu, sigma, nu, lower_tail, !log))
// # all.equal(my_pDEL(q, mu, sigma, nu, lower_tail, log, 20), gamlss.dist::pDEL(q, mu, sigma, nu, lower_tail, log))
//
// # autoplot(microbenchmark(my_pDEL(q, mu, sigma, nu, lower_tail, log, 20), gamlss.dist::pDEL(q, mu, sigma, nu, lower_tail, log), times = 10))
//
// # all.equal(my_qDEL(p, mu, sigma, nu, lower_tail, log, n_cpu), gamlss.dist::qDEL(p, mu, sigma, nu, lower_tail, log, max_value))
// # all.equal(my_qDEL(p, mu, sigma, nu, !lower_tail, log, n_cpu), gamlss.dist::qDEL(p, mu, sigma, nu, !lower_tail, log, max_value)) # FIXME
// # all.equal(my_qDEL(p, mu, sigma, nu, lower_tail, log, n_cpu), gamlss.dist::qDEL(p, mu, sigma, nu, lower_tail, log, max_value))
// # all.equal(my_qDEL(p, mu, sigma, nu, lower_tail, log, 20), gamlss.dist::qDEL(p, mu, sigma, nu, lower_tail, log, max_value))
//
// # autoplot(microbenchmark(my_qDEL(p, mu, sigma, nu, lower_tail, log, 20), gamlss.dist::qDEL(p, mu, sigma, nu, lower_tail, log, max_value), times = 10))
