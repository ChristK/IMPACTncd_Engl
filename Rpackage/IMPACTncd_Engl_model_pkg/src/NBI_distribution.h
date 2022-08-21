#ifndef NBI_distribution_H
#define NBI_distribution_H


double my_dNBI_scalar(const int& x,
                      const double& mu = 1.0,
                      const double& sigma = 1.0,
                      const bool& log_p = false,
                      const bool& check = true);

Rcpp::NumericVector my_dNBI(const Rcpp::IntegerVector& x,
                      const Rcpp::NumericVector& mu,
                      const Rcpp::NumericVector& sigma,
                      const bool& log_p = false,
                      const int& n_cpu = 1);

double my_pNBI_scalar(const int& q,
                      const double& mu = 1.0,
                      const double& sigma = 1.0,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const bool& check = true);


Rcpp::NumericVector my_pNBI(const Rcpp::IntegerVector& q,
                      const Rcpp::NumericVector& mu,
                      const Rcpp::NumericVector& sigma,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& n_cpu = 1);

int my_qNBI_scalar(const double& p,
                   const double& mu = 1.0,
                   const double& sigma = 1.0,
                   const bool& lower_tail = true,
                   const bool& log_p = false,
                   const bool& check = true);

Rcpp::IntegerVector my_qNBI(const Rcpp::NumericVector& p,
                      const Rcpp::NumericVector& mu,
                      const Rcpp::NumericVector& sigma,
                      const bool& lower_tail = true,
                      const bool& log_p = false,
                      const int& n_cpu = 1);

int my_qZANBI_scalar(const double& p,
                     const double& mu = 1.0,
                     const double& sigma = 1.0,
                     const double& nu = 0.3,
                     const bool& lower_tail = true,
                     const bool& log_p = false,
                     const bool& check = true);

Rcpp::IntegerVector my_qZANBI(const Rcpp::NumericVector& p,
                        const Rcpp::NumericVector& mu,
                        const Rcpp::NumericVector& sigma,
                        const Rcpp::NumericVector& nu,
                        const bool& lower_tail = true,
                        const bool& log_p = false,
                        const int& n_cpu = 1);

double my_pZANBI_scalar(const int& q,
                     const double& mu = 1.0,
                     const double& sigma = 1.0,
                     const double& nu = 0.3,
                     const bool& lower_tail = true,
                     const bool& log_p = false,
                     const bool& check = true);


#endif
