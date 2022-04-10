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

// TODO parallel

// from gamlss.dist
NumericVector my_dDPOgetC5_C (const NumericVector& mu, const NumericVector& sigma,
                           const int& lmu, const int& ly) {
  double   sumC, mus, lmus, lsig2, invs, ls;
  NumericVector ylogofy(ly), lga(ly), ym(ly), ans(lmu);

  int i,j;
  for (j=0 ; j < ly ; j++) {
    ylogofy[j] = j * ((j==0)? 1 : log(j));
    lga[j] = R::lgammafn(j + 1);
    ym[j] = (j - ylogofy[j]);
  }
  for (i=0; i < lmu; i++){
    sumC = 0;
    mus = mu[i] / sigma[i];
    lsig2 = -0.5 * log(sigma[i]);
    lmus = log(mu[i]) / sigma[i] - 1;
    invs = 1 / sigma[i];
    ls = lsig2 - mus;
    for (j=0 ; j < ly ; j++){
      sumC += exp(ls - lga[j] + ylogofy[j] + j * lmus + invs * ym[j]);
    }
    ans[i] = pow(sumC,-1);
  }
  return ans;
}

double my_dDPOgetC5_C_scalar (const double& mu, const double& sigma,
                              const int& lmu, const int& ly) {
  double   sumC, mus, lmus, lsig2, invs, ls, ans(lmu);
  NumericVector ylogofy(ly), lga(ly), ym(ly);

  int i,j;
  for (j=0 ; j < ly ; j++) {
    ylogofy[j] = j * ((j==0)? 1 : log(j));
    lga[j] = R::lgammafn(j + 1);
    ym[j] = (j - ylogofy[j]);
  }
  for (i=0; i < lmu; i++){
    sumC = 0;
    mus = mu / sigma;
    lsig2 = -0.5 * log(sigma);
    lmus = log(mu) / sigma - 1;
    invs = 1 / sigma;
    ls = lsig2 - mus;
    for (j=0 ; j < ly ; j++){
      sumC += exp(ls - lga[j] + ylogofy[j] + j * lmus + invs * ym[j]);
    }
    ans = pow(sumC,-1);
  }
  return ans;
}

// this is fast
//' @export
// [[Rcpp::export]]
NumericVector my_get_C(const IntegerVector& x,
                       const NumericVector& mu,
                       const NumericVector& sigma)
{
  int maxV = std::max(Rcpp::max(x) * 3, 500);
  int lmu   = std::max(std::max(x.length(), mu.length()), sigma.length());
  return log(my_dDPOgetC5_C(mu, sigma, lmu, maxV + 1));
}




// dDPO

//' @export
// [[Rcpp::export]]
NumericVector my_dDPO(const IntegerVector& x,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const bool& log_ = false,
                      const int& n_cpu = 1)
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

  NumericVector logofx(n, 1.0); // filled with 1.0
  NumericVector lh(n);
  NumericVector theC =     my_get_C(x, mu, sigma);
  for (int i = 0; i < n; i++)
  {
    if (x[i] > 0.0) logofx[i] = log(x[i]);
    lh[i] = -0.5 * log(sigma[i]) - (mu[i]/sigma[i]) - R::lgammafn(x[i] + 1) + x[i] *
      logofx[i] - x[i] + (x[i] * log(mu[i]))/sigma[i] + x[i]/sigma[i] - (x[i] *
      logofx[i])/sigma[i] + theC[i];
  }
  if (!log_) lh = exp(lh);

  return lh;
}

// pDPO

//' @export
// [[Rcpp::export]]
NumericVector my_pDPO(const IntegerVector& q,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const bool& lower_tail = true,
                      const bool& log_p = false)
{
  if (q.length() != mu.length()) stop("Distribution parameters must be of same length (q mu)");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length (mu sigma)");
  const int n = q.length();
  for (int i = 0; i < n; i++)
  {
    if (q[i] < 0) stop("q must be >=0");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
  }


  int maxV = std::max(Rcpp::max(q) * 3, 500);

  NumericVector den(n);
  for (int i = 0; i < n; i++)
  {
    den[i] = my_dDPOgetC5_C_scalar(mu[i], sigma[i], 1, q[i] + 1 );
  }

  NumericVector num = my_dDPOgetC5_C(mu, sigma, n, maxV + 1);

  NumericVector cdf = num/den;

  if (!lower_tail) cdf = 1.0 - cdf;
  if (log_p) cdf = log(cdf);

  return cdf;
}

double my_pDPO_scalar(const int& q,
                      const double& mu,
                      const double& sigma,
                      const bool& lower_tail = true,
                      const bool& log_p = false)
{

    // if (q < 0) stop("q must be >=0");
    //
    // if (mu<= 0.0) stop("mu must be greater than 0");
    //
    // if (sigma <= 0.0) stop("sigma must be greater than 0");


  int maxV = std::max(q * 3, 500);

  double den = my_dDPOgetC5_C_scalar(mu, sigma, 1, q + 1);

  double num = my_dDPOgetC5_C_scalar(mu, sigma, 1, maxV + 1);

  double cdf = num/den;

  if (!lower_tail) cdf = 1.0 - cdf;
  if (log_p) cdf = log(cdf);

  return cdf;
}

// qDPO

// slow and restricted by max_value
IntegerVector my_qDPO2(      NumericVector p,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& max_value = 10000)
{
  if (p.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  const int n = p.length();
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0001) stop("p must be between 0 and 1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
  }

  if (log_p) p = exp(p);
  if (!lower_tail) p = 1.0 - p;

  IntegerVector QQQ(n);

  for (int i = 0; i < n; i++)
  {
    if (p[i] + 1e-09 >= 1.0) QQQ[i] = R_PosInf;
    else
      {
      for (int j = 0; j <= max_value; j++)
        {
        double cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
            if (p[i] <= cumpro) {
              QQQ[i] = j;
              break;
            }
        }
      }
  }
  return QQQ;
}


// fast! and not restricted by max_value
//' @export
// [[Rcpp::export]]
IntegerVector my_qDPO(      NumericVector p,
                            const NumericVector& mu,
                            const NumericVector& sigma,
                            const bool& lower_tail = true,
                            const bool& log_p = false,
                            const int& max_value = 0,
                            const int& n_cpu = 1)
{
  // TODO implement recycling with rep_len()
  if (p.length() != mu.length()) stop("Distribution parameters must be of same length");
  if (sigma.length() != mu.length()) stop("Distribution parameters must be of same length");
  const int n = p.length();

  if (log_p) p = exp(p);
  if (!lower_tail) p = 1.0 - p;

  IntegerVector QQQ(n);
  double cumpro = 0.0;
  int j = max_value;
  omp_set_num_threads(n_cpu); // Use n_cpu threads for all
  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0001) stop("p must be between 0 and 1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
  }

  for (int i = 0; i < n; i++)
  {
    if (p[i] + 1e-09 >= 1.0) QQQ[i] = R_PosInf; // I don't like this
    else
    {
      // implement a divide & conquer algo
      if (max_value == 0) {
        j = mu[i] * (p[i] + 0.5) + 1; // guess a starting value
      } // else j = max_value;
      cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
      if (p[i] <= cumpro)
      {
        while (j > 0 && p[i] <= cumpro)
        {
          j /= 2;
          cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
        }
        for (int k = j; k <= (j * 2 + 1); k++)
        {
          cumpro = my_pDPO_scalar(k, mu[i], sigma[i], true, false);
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
          cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
        }
        j /= 2;
        cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
        if ((j*2) > (j+1000))
        {
          while (j < INT_MAX && p[i] > cumpro)
          {
            j += 1000;
            cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
          }
          if (j >= 1000) j -= 1000;
          cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
        }

        while (j < INT_MAX && p[i] > cumpro)
        {
          j += 100;
          cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
        }
        if (j >= 100) j -= 100;
        cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
        while (j < INT_MAX && p[i] > cumpro)
        {
          j += 10;
          cumpro = my_pDPO_scalar(j, mu[i], sigma[i], true, false);
        }
        if (j >= 10) j -= 10;
        for (int k = j; k <= INT_MAX; k++)
        {
          cumpro = my_pDPO_scalar(k, mu[i], sigma[i], true, false);
          if (p[i] <= cumpro)
          {
            QQQ[i] = k;
            break;
          }
        }
      }
    }
  }
  return QQQ;
}

/*** R
# N <- 1e2
# n_cpu <- 1L
# x <- sample(1e2, N, TRUE)
# q <- sample(1e2, N, TRUE)
# p <- runif(N, 0, 0.999)
# mu <- runif(N, 0, 1e2)
# sigma <- runif(N, 0, 1e2)
# nu <- runif(N, 0, 1e2)
# tau <- runif(N, 0, 1)
# log <- FALSE
# lower_tail <- TRUE
# max_value <- 10000L
#
# library(microbenchmark)
# library(ggplot2)
# library(gamlss)

# all.equal(my_get_C(x, mu, sigma), gamlss.dist::get_C(x, mu, sigma))
# autoplot(microbenchmark(my_get_C(x, mu, sigma), gamlss.dist::get_C(x, mu, sigma), times = 10))

# all.equal(my_dDPO(x, mu, sigma, log, n_cpu), gamlss.dist::dDPO(x, mu, sigma, log))
# autoplot(microbenchmark(my_dDPO(x, mu, sigma, log), gamlss.dist::dDPO(x, mu, sigma, log), times = 10))

# all.equal(my_pDPO(q, mu, sigma, lower_tail, log), gamlss.dist::pDPO(q, mu, sigma, lower_tail, log))
# autoplot(microbenchmark(my_pDPO(q, mu, sigma, lower_tail, log), gamlss.dist::pDPO(q, mu, sigma, lower_tail, log), times = 10))

# all.equal(my_qDPO(p, mu, sigma, lower_tail, log, n_cpu), gamlss.dist::qDPO(p, mu, sigma, lower_tail, log, max_value))
# all.equal(my_qDPO(p, mu, sigma, lower_tail, log, 2), gamlss.dist::qDPO(p, mu, sigma, lower_tail, log, max_value))

# autoplot(microbenchmark(my_qDPO(p, mu, sigma, lower_tail, log), gamlss.dist::qDPO(p, mu, sigma, lower_tail, log, max_value), times = 10))
# autoplot(microbenchmark(my_qDPO(p, mu, sigma, lower_tail, log, 20), my_qDPO(p, mu, sigma, lower_tail, log, 1), times = 10))

# all.equal(my_qDPO2(p, mu, sigma, lower_tail, log, max_value), my_qDPO(p, mu, sigma, lower_tail, log))
# autoplot(microbenchmark(my_qDPO2(p, mu, sigma, lower_tail, log, max_value), my_qDPO(p, mu, sigma, lower_tail, log), times = 10))

# xx <- my_qDPO(p, mu, sigma, lower_tail, log)
# plot(mu*p, xx)
*/
