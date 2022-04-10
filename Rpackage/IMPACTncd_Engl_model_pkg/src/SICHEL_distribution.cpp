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

// TODO parallelise

// from gamlss.dist
NumericVector my_cdfSICHEL(const IntegerVector& y,
                           const NumericVector& mu,
                           const NumericVector& sigma,
                           const NumericVector& nu) {
  int ny = y.length();
  NumericVector ans = no_init(ny);
  int i, j, lyp1;
  //double alpha[*ny], cvec[*ny];
  double sumT, lbes, cvec, alpha;
  for (i = 0; i < ny; i++)
  {
    lyp1 = y[i]+1;
    sumT = 0.0;
    cvec = exp(log(R::bessel_k(1/sigma[i], nu[i]+1, 1)) - log(R::bessel_k(1/sigma[i], nu[i], 1)) );
    alpha = sqrt(1 + 2*sigma[i]*mu[i]/cvec)/sigma[i];
    lbes = log(R::bessel_k(alpha, nu[i]+1, 1)) - log(R::bessel_k(alpha, nu[i], 1));
    double tynew[lyp1];
    double lpnew[lyp1];
    tynew[0] = (mu[i]/cvec)*pow(1+2*sigma[i]*mu[i]/cvec,-0.5)*exp(lbes);
    lpnew[0] = -nu[i]*log(sigma[i]*alpha) +  log(R::bessel_k(alpha,nu[i],1)) - log(R::bessel_k(1/sigma[i],nu[i],1));
    for (j = 1; j < lyp1; j++)
    {
      tynew[j] = (cvec*sigma[i]*(2*(j + nu[i])/mu[i])+(1/tynew[j-1]))*pow(mu[i]/(sigma[i]*alpha*cvec),2);
      lpnew[j] = lpnew[j-1] + log(tynew[j-1]) - log(j);
    }
    for (j=0 ; j< lyp1; j++) sumT += exp(lpnew[j]);
    ans[i] = sumT;
  }
  return ans;
}

double my_cdfSICHEL_scalar(const int& y,
                           const double& mu,
                           const double& sigma,
                           const double& nu) {
  int j, lyp1;
  //double alpha[*ny], cvec[*ny];
  double sumT, lbes, cvec, alpha;
    lyp1 = y + 1;
    sumT = 0.0;
    cvec = exp(log(R::bessel_k(1/sigma, nu+1, 1)) - log(R::bessel_k(1/sigma, nu, 1)) );
    alpha = sqrt(1 + 2*sigma*mu/cvec)/sigma;
    lbes = log(R::bessel_k(alpha, nu+1, 1)) - log(R::bessel_k(alpha, nu, 1));
    double tynew[lyp1];
    double lpnew[lyp1];
    tynew[0] = (mu/cvec)*pow(1+2*sigma*mu/cvec,-0.5)*exp(lbes);
    lpnew[0] = -nu*log(sigma*alpha) +  log(R::bessel_k(alpha,nu,1)) - log(R::bessel_k(1/sigma,nu,1));
    for (j = 1; j < lyp1; j++)
    {
      tynew[j] = (cvec*sigma*(2*(j + nu)/mu)+(1/tynew[j-1]))*pow(mu/(sigma*alpha*cvec),2);
      lpnew[j] = lpnew[j-1] + log(tynew[j-1]) - log(j);
    }
    for (j=0 ; j< lyp1; j++) sumT += exp(lpnew[j]);
  return sumT;
}


// from gamlss.dist
NumericVector my_tofySICHEL2 (const NumericVector& y, const NumericVector& mu,
                              const NumericVector& sigma,
                              const NumericVector& nu, const NumericVector& lbes,
                              const NumericVector& cvec)
{
  int maxyp1 = max(y) + 2;
  int ny = y.length();
  NumericVector ans(ny);
  double tofY[maxyp1];
  int iy, i, j;
  double alpha, sumT;
  for (i = 0; i < ny; i++)
  {
    iy = y[i]+1;
    tofY[0] = (mu[i]/cvec[i])*pow(1+2*sigma[i]*mu[i]/cvec[i],-0.5)*exp(lbes[i]);
    alpha = sqrt(1 + 2*sigma[i]*mu[i]/cvec[i])/sigma[i];
    sumT = 0.0;
    for (j = 1; j < iy; j++)
    {
      tofY[j] = (cvec[i]*sigma[i]*(2*(j+nu[i])/mu[i]) + (1/tofY[j-1])) * pow(mu[i]/(sigma[i]*alpha*cvec[i]),2);
      sumT += log(tofY[j-1]);
    }
    ans[i] = sumT;
  }
  return ans;
}

double my_tofySICHEL2_scalar (const int& y, const double& mu,
                              const double& sigma,
                              const double& nu, const double& lbes,
                              const double& cvec)
{
  double tofY[y + 2];
  int iy = y+1;
  tofY[0] = (mu/cvec)*pow(1+2*sigma*mu/cvec,-0.5)*exp(lbes);
  double alpha = sqrt(1 + 2*sigma*mu/cvec)/sigma;
  double sumT = 0.0;
  for (int j = 1; j < iy; j++)
  {
    tofY[j] = (cvec*sigma*(2*(j+nu)/mu) + (1/tofY[j-1])) * pow(mu/(sigma*alpha*cvec),2);
    sumT += log(tofY[j-1]);
  }
  return sumT;
}

// dNBI (necessary for dSICHEL)
double my_dNBI_scalar (const int& x, const double& mu = 1.0, const double& sigma = 1.0, const bool& log_ = false)
{
  if (x < 0) stop("x must be >=0");
  if (mu <= 0.0) stop("mu must be greater than 0");
  if (sigma <= 0.0) stop("sigma must be greater than 0");
  double fy;
  if (sigma > 1e-04)
  {
    fy = R::dnbinom_mu(x, 1.0/sigma, mu, log_);
  }
  else
  {
    fy = R::dpois(x, mu, log_);
  }
  return fy;
}

// dSICHEL
//' @export
// [[Rcpp::export]]
NumericVector my_dSICHEL(const IntegerVector& x,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const NumericVector& nu,
                      const bool& log_ = false,
                      const int& n_cpu = 1)
{
  // OMP failed. R::bessel_k could be the cause as is probably not thread safe
  if (x.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length()
  ) stop("Distribution parameters must be of same length");
  // TODO recycle if unequal lengths

  const int n = x.length();
  NumericVector logfy = no_init(n);
  double cvec = 0.0;
  double alpha = 0.0;
  double lbes = 0.0;
  double sumlty = 0.0;

  for (int i = 0; i < n; i++)
  {
    if (x[i] < 0.0) stop("x must be >=0");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    cvec = exp(log(R::bessel_k((1.0/sigma[i]), nu[i] + 1.0, 1.0)) - log(R::bessel_k((1.0/sigma[i]),
                                  nu[i], 1.0)));
    alpha = sqrt(1.0 + 2.0 * sigma[i] * (mu[i]/cvec))/sigma[i];
    lbes = log(R::bessel_k(alpha, nu[i] + 1.0, 1.0)) - log(R::bessel_k(alpha,
                              nu[i], 1.0));
    sumlty = my_tofySICHEL2_scalar(x[i], mu[i], sigma[i], nu[i], lbes, cvec);
    logfy[i] = -R::lgammafn(x[i] + 1.0) - nu[i] * log(sigma[i] * alpha) + sumlty +
      log(R::bessel_k(alpha, nu[i], 1.0)) - log(R::bessel_k((1.0/sigma[i]), nu[i], 1.0));

    if (!log_) logfy[i] = exp(logfy[i]);
    if (sigma[i] > 10000.0 && nu[i] > 0.0) logfy[i] = my_dNBI_scalar(x[i], mu[i], 1.0/nu[i], log_);
  }

  // If nu is far away from 0 (i.e. 100 or -100) bessel_k produces very small numbers
  // that lead to NaNs in the output
  if (any(is_na(logfy))) warning("NaNs or NAs were produced");
  return logfy; // despite the name, not logged if log_ = false
}

// pSICHEL

//' @export
// [[Rcpp::export]]
NumericVector my_pSICHEL(const IntegerVector& q,
                      const NumericVector& mu,
                      const NumericVector& sigma,
                      const NumericVector& nu,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& n_cpu = 1)
{
  if (q.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length()
  ) stop("Distribution parameters must be of same length");

  const int n = q.length();
  // NumericVector cdf(n);
 omp_set_num_threads(n_cpu); // Use n_cpu threads for all

  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (q[i] < 0) stop("q must be >=0");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    // cdf[i] = my_cdfSICHEL_scalar(q[i], mu[i], sigma[i], nu[i]);
  }
  NumericVector cdf = my_cdfSICHEL(q, mu, sigma, nu);

  if (!lower_tail) cdf = 1.0 - cdf;
  if (log_p) cdf = log(cdf);
  if (any(is_na(cdf))) warning("NaNs or NAs were produced");

  return cdf;
}


double my_pSICHEL_scalar(int q,
                         const double& mu,
                         const double& sigma,
                         const double& nu,
                         const bool& lower_tail = true,
                         const bool& log_p = false)
{
    // if (q < 0) stop("q must be >=0");
    // if (mu <= 0.0) stop("mu must be greater than 0");
    // if (sigma <= 0.0) stop("sigma must be greater than 0");

  double cdf = my_cdfSICHEL_scalar(q, mu, sigma, nu);

  if (!lower_tail) cdf = 1.0 - cdf;
  if (log_p) cdf = log(cdf);

  return cdf;
}


// double my_pSICHEL_scalar(const int& q,
//                       const double& mu,
//                       const double& sigma,
//                       const bool& lower_tail = true,
//                       const bool& log_p = false)
// {
//
//   if (q < 0) stop("q must be >=0");
//
//   if (mu<= 0.0) stop("mu must be greater than 0");
//
//   if (sigma <= 0.0) stop("sigma must be greater than 0");
//
//
//   int maxV = std::max(q * 3, 500);
//
//   double den = my_dSICHELgetC5_C_scalar(mu, sigma, 1, q + 1 );
//
//   double num = my_dSICHELgetC5_C_scalar(mu, sigma, 1, maxV + 1);
//
//   double cdf = num/den;
//
//   if (!lower_tail) cdf = 1.0 - cdf;
//   if (log_p) cdf = log(cdf);
//
//   return cdf;
// }

// qSICHEL

// fast! and not restrivted by max_value

//' @export
// [[Rcpp::export]]
IntegerVector my_qSICHEL(      NumericVector p,
                            const NumericVector& mu,
                            const NumericVector& sigma,
                            const NumericVector& nu,
                            const bool& lower_tail = true,
                            const bool& log_p = false,
                            const int& n_cpu = 1)
{
  // TODO implement recycling with rep_len()
  if (p.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length()
  ) stop("Distribution parameters must be of same length");

  const int n = p.length();
  IntegerVector QQQ(n);
  double cumpro = 0.0;
  int j = 0;
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0001) stop("p must be between 0 and 1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");

    if (log_p) p[i] = exp(p[i]);
    if (!lower_tail) p[i] = 1.0 - p[i];

    if (p[i] + 1e-09 >= 1.0) QQQ[i] = R_PosInf; // I don't like this
    else
    {
      // implement a divide & conquer algo
      j = mu[i] * (p[i] + 0.5) + 1; // guess a starting value
      cumpro = my_pSICHEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
      if (p[i] <= cumpro)
      {
        while (j > 0 && p[i] <= cumpro)
        {
          j /= 2;
          cumpro = my_pSICHEL_scalar(j, mu[i], sigma[i], nu[i],  true, false);
        }
        for (int k = j; k <= (j * 2 + 1); k++)
        {
          cumpro = my_pSICHEL_scalar(k, mu[i], sigma[i], nu[i], true, false);
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
          cumpro = my_pSICHEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        }
        j /= 2;
        cumpro = my_pSICHEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        if ((j*2) > (j+1000))
        {
          while (j < INT_MAX && p[i] > cumpro)
          {
            j += 1000;
            cumpro = my_pSICHEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
          }
          if (j >= 1000) j -= 1000;
          cumpro = my_pSICHEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        }

        while (j < INT_MAX && p[i] > cumpro)
        {
          j += 100;
          cumpro = my_pSICHEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        }
        if (j >= 100) j -= 100;
        cumpro = my_pSICHEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        while (j < INT_MAX && p[i] > cumpro)
        {
          j += 10;
          cumpro = my_pSICHEL_scalar(j, mu[i], sigma[i], nu[i], true, false);
        }
        if (j >= 10) j -= 10;
        for (int k = j; k <= INT_MAX; k++)
        {
          cumpro = my_pSICHEL_scalar(k, mu[i], sigma[i], nu[i], true, false);
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
// fast! and not restricted by max_value
//' @export
// [[Rcpp::export]]
IntegerVector my_qZISICHEL(          NumericVector  p,
                               const NumericVector& mu,
                               const NumericVector& sigma,
                               const NumericVector& nu,
                               const NumericVector& tau,
                               const bool& lower_tail = true,
                               const bool& log_p = false,
                               const int& n_cpu = 1)
{
  // TODO implement recycling with rep_len()
  if (p.length() != mu.length() || mu.length() != sigma.length() ||
      sigma.length() != nu.length() || nu.length() != tau.length()
  ) stop("Distribution parameters must be of same length");

  const int n = p.length();
  NumericVector pnew(n); // = no_init(n);

 omp_set_num_threads(n_cpu); // Use n_cpu threads for all
  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (p[i] < 0.0 || p[i] > 1.0001) stop("p must be between 0 and 1");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    if (tau[i] <= 0.0 || tau[i] >= 1.0) stop("tau must be between 0 and 1");

    if (log_p) p[i] = exp(p[i]);
    if (!lower_tail) p[i] = 1.0 - p[i];

    pnew[i] = (p[i] - tau[i])/(1.0 - tau[i]) - (1e-07);

    if (pnew[i] < 0.0) pnew[i] = 0.0;
  }

  return my_qSICHEL(pnew, mu, sigma, nu);
}


// fast! and not restricted by max_value
//' @export
// [[Rcpp::export]]
NumericVector my_pZISICHEL(          NumericVector  q,
  const NumericVector& mu,
  const NumericVector& sigma,
  const NumericVector& nu,
  const NumericVector& tau,
  const bool& lower_tail = true,
  const bool& log_p = false,
  const int& n_cpu = 1)
{
  // TODO implement recycling with rep_len()
  if (q.length() != mu.length() || mu.length() != sigma.length() ||
    sigma.length() != nu.length() || nu.length() != tau.length()
  ) stop("Distribution parameters must be of same length");

  const int n = q.length();
  NumericVector cdf = no_init(n);

  omp_set_num_threads(n_cpu); // Use n_cpu threads for all
  // consecutive parallel regions
#pragma omp parallel for default(shared)
  for (int i = 0; i < n; i++)
  {
    if (q[i] < 0.0) stop("q must be >=0");
    if (mu[i] <= 0.0) stop("mu must be greater than 0");
    if (sigma[i] <= 0.0) stop("sigma must be greater than 0");
    if (tau[i] <= 0.0 || tau[i] >= 1.0) stop("tau must be between 0 and 1");

    cdf[i] = my_pSICHEL_scalar(q[i], mu[i], sigma[i], nu[i]);
    cdf[i] = tau[i] + (1.0 - tau[i]) * cdf[i];

    if (!lower_tail) cdf[i] = 1.0 - cdf[i];
    if (log_p) cdf[i] = log(cdf[i]);
  }

  return cdf;
}


// /*** R
// # N <- 1e4
// # n_cpu <- 1L
// # x <- sample(1e2, N, TRUE)
// # q <- sample(1e2, N, TRUE)
// # p <- runif(N, 0, 0.999)
// # mu <- runif(N, 0, 1e2)
// # sigma <- runif(N, 0, 1e3) # gamlss bug when sigma > 1e4
// # nu <- runif(N, -5, 5) # large values return NAN
// # tau <- runif(N, 0, 1)
// # log <- FALSE
// # lower_tail <- TRUE
// # max_value <- 10000L
// # library(microbenchmark)
// # library(ggplot2)
// # library(gamlss)
// # all.equal(my_dSICHEL(x, mu, sigma, nu, log, 1L), gamlss.dist::dSICHEL(x, mu, sigma, nu, log))
// # all.equal(my_dSICHEL(x, mu, sigma, nu, !log, 1L), gamlss.dist::dSICHEL(x, mu, sigma, nu, !log))
// # autoplot(microbenchmark(my_dSICHEL(x, mu, sigma, nu, log, 20L), gamlss.dist::dSICHEL(x, mu, sigma, nu, log), times = 10))
// # all.equal(my_pSICHEL(q, mu, sigma, nu, lower_tail, log, 1), gamlss.dist::pSICHEL(q, mu, sigma, nu, lower_tail, log))
// # all.equal(my_pSICHEL(q, mu, sigma, nu, lower_tail, !log, 1), gamlss.dist::pSICHEL(q, mu, sigma, nu, lower_tail, !log))
// # all.equal(my_pSICHEL(q, mu, sigma, nu, !lower_tail, log, 1), gamlss.dist::pSICHEL(q, mu, sigma, nu, !lower_tail, log))
// # all.equal(my_pSICHEL(q, mu, sigma, nu, lower_tail, log, 20), my_pSICHEL(q, mu, sigma, nu, lower_tail, log, 1))
// # autoplot(microbenchmark(my_pSICHEL(q, mu, sigma, nu, lower_tail, log, 20), gamlss.dist::pSICHEL(q, mu, sigma, nu, lower_tail, log), times = 10))
// # all.equal(my_qSICHEL(p, mu, sigma, nu, lower_tail, log), gamlss.dist::qSICHEL(p, mu, sigma, nu, lower_tail, log, max_value))
// # all.equal(my_qSICHEL(p, mu, sigma, nu, !lower_tail, log), gamlss.dist::qSICHEL(p, mu, sigma, nu, !lower_tail, log, max_value)) # FIXME
// # all.equal(my_qSICHEL(p, mu, sigma, nu, lower_tail, !log), gamlss.dist::qSICHEL(p, mu, sigma, nu, lower_tail, !log, max_value))
// # autoplot(microbenchmark(my_qSICHEL(p, mu, sigma, nu, lower_tail, log), gamlss.dist::qSICHEL(p, mu, sigma, nu, lower_tail, log, max_value), times = 10))
// # all.equal(my_qZISICHEL(p, mu, sigma, nu, tau, lower_tail, log), gamlss.dist::qZISICHEL(p, mu, sigma, nu, tau, lower_tail, log, max_value))
// # autoplot(microbenchmark(my_qZISICHEL(p, mu, sigma, nu, tau, lower_tail, log), gamlss.dist::qZISICHEL(p, mu, sigma, nu, tau, lower_tail, log, max_value), times = 10))
// # xx <- my_qSICHEL(p, mu, sigma, lower_tail, log)
// # plot(mu*p, xx)
// # for (i in 1:N) {
// #  if (!identical(my_dSICHEL(x[i], mu[i], sigma[i], nu[i], log, n_cpu),
// #            gamlss.dist::dSICHEL(x[i], mu[i], sigma[i], nu[i], log)))   print(i)
// # }
// # dSICHEL(x = c(1, 1), mu = c(1, 1), sigma = c(1e5, 1e5), nu = c(-1, 1))
// # ifelse((sigma > 10000) & (nu > 0),
// #       dNBI(x, mu = mu[which((sigma > 10000) & (nu > 0))],
// #            sigma = 1/nu[which((sigma > 10000) & (nu > 0))],
// #            log = log),
// #       gamlss.dist::dSICHEL(x[which(!((sigma > 10000) & (nu > 0)))],
// #                            mu[which(!((sigma > 10000) & (nu > 0)))],
// #                            sigma[which(!((sigma > 10000) & (nu > 0)))],
// #                            nu[which(!((sigma > 10000) & (nu > 0)))],
// #                            log))
// # ifelse((sigma > 10000) & (nu > 0),
// #       dNBI(x, mu = mu[which((sigma > 10000) & (nu > 0))], sigma = 1/nu[which((sigma > 10000) & (nu > 0))], log = log), NA)
// */
