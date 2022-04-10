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
#include <Rmath.h>
#include <algorithm>
using namespace Rcpp;

// code based on https://root-forum.cern.ch/t/two-weighted-sample-tests-kolmogorov-smirnov-cramer-von-mises-anderson-darling-tests/33378
// from jakub.trusina@fjfi.cvut.cz
// Original code has homogeneity tests for weighted data samples - Kolmogorov-Smirnov, Cramér-von Mises, Anderson-Darling test
// I will only use code for the Anderson-Darling statistic, without calculating the p value
// The p-value will be calculated in R with non parametric bootstrap

// if a sample is unweighted pass a vector of ones

//' @export
// [[Rcpp::export]]
double wtd_ADstat(const NumericVector& a, const NumericVector& wa,
                         const NumericVector& b, const NumericVector& wb)
{
  // a: sample 1 needs to be sorted
  // wa: sample 1 weights
  // b: sample2 needs to be sorted
  // wb: sample2 weights
  int na = a.length();
  int nb = b.length();
  double ADstat = 0.0;
    // Require at least two points in each graph
  // if (na <= 2 || nb <= 2) {
  //   stop("Sets must have more than 2 points");
  // }
  // if (na != wa.length() || nb != wb.length()) {
  //   stop("Weight vectors need to be equal length with their respective sample version");
  // }
  // Constants needed
  int ia = 0;
  int ib = 0;
  double eventsa = 0;
  double eventsb = 0;
  double sq_eventsa = 0;
  double sq_eventsb = 0;

  // Calculating effective entries size
  for (int i = 0; i < na; i++) {
    eventsa += wa[i];
    sq_eventsa += std::pow(wa[i], 2.0);
  }
  for (int j = 0; j < nb; j++) {
    eventsb += wb[j];
    sq_eventsb += std::pow(wb[j], 2.0);
  }
  double effna = std::pow(eventsa, 2.0) / sq_eventsa;
  double effnb = std::pow(eventsb, 2.0) / sq_eventsb;
  double effn = effna + effnb;

  // Auxiliary variables
  double x;
  double sumwa = 0;
  double sumwb = 0;
  double FA = 0;
  double FB = 0;
  double H = 0;
  bool enda = false;
  bool endb = false;

  // Main loop over point sets
  while ((!enda) && (!endb)) {
    if (enda) {
      x = b[ib];
    }
    else if (endb) {
      x = a[ia];
    }
    else {
      x = std::min(a[ia], b[ib]);
    }
    while (a[ia] == x && ia < na) {
      sumwa += wa[ia];
      ia++;
    }
    while (b[ib] == x && ib < nb) {
      sumwb += wb[ib];
      ib++;
    }
    FA += sumwa / eventsa;
    FB += sumwb / eventsb;
    H = (effna * FA + effnb * FB) / effn;
    if ((H > 0) && (H < 1)) {
      ADstat += std::pow(FA - FB, 2.0) / ((1 - H) * H) * (sumwa / eventsa * effna + sumwb / eventsb * effnb) / effn;
    }

    // reseting sumwa, sumwb
    sumwa = 0;
    sumwb = 0;

    // set last point to infinity
    if (ia == na) {
      enda = true;
    }
    if (ib == nb) {
      endb = true;
    }
  }

  // Computing AD test's p-value
  ADstat = effna * effnb / effn * ADstat;

  return ADstat; //Anderson-Darling test statistic
}

//' @export
// [[Rcpp::export]]
double wtd_KSstat(const NumericVector& a, const NumericVector& wa,
                  const NumericVector& b, const NumericVector& wb)
{
  // a: sample 1 need to be sorted
  // wa: sample 1 weights
  // b: sample2 need to e sorted
  // wb: sample2 weights
  int na = a.length();
  int nb = b.length();
  double KSstat = 0.0;
  // Require at least two points in each graph
  // if (na <= 2 || nb <= 2) {
  //   stop("Sets must have more than 2 points");
  // }
  // if (na != wa.length() || nb != wb.length()) {
  //   stop("Weight vectors need to be equal length with their respective sample version");
  // }
  // Constants needed
  int ia = 0;
  int ib = 0;
  double eventsa = 0;
  double eventsb = 0;
  double sq_eventsa = 0;
  double sq_eventsb = 0;

  // Calculating effective entries size
  for (int i = 0; i < na; i++) {
    eventsa += wa[i];
    sq_eventsa += std::pow(wa[i], 2.0);
  }
  for (int j = 0; j < nb; j++) {
    eventsb += wb[j];
    sq_eventsb += std::pow(wb[j], 2.0);
  }
  double effna = std::pow(eventsa, 2.0) / sq_eventsa;
  double effnb = std::pow(eventsb, 2.0) / sq_eventsb;
  double effn = effna + effnb;

  // Auxiliary variables
  double x;
  double sumwa = 0;
  double sumwb = 0;
  double FA = 0;
  double FB = 0;
  double H = 0;
  bool enda = false;
  bool endb = false;

  // Main loop over point sets
  while ((!enda) && (!endb)) {
    if (enda) {
      x = b[ib];
    }
    else if (endb) {
      x = a[ia];
    }
    else {
      x = std::min(a[ia], b[ib]);
    }
    while (a[ia] == x && ia < na) {
      sumwa += wa[ia];
      ia++;
    }
    while (b[ib] == x && ib < nb) {
      sumwb += wb[ib];
      ib++;
    }
    FA += sumwa / eventsa;
    FB += sumwb / eventsb;
    H = (effna * FA + effnb * FB) / effn;
    KSstat = std::max(KSstat, std::abs(FA - FB));

    // reseting sumwa, sumwb
    sumwa = 0;
    sumwb = 0;

    // set last point to infinity
    if (ia == na) {
      enda = true;
    }
    if (ib == nb) {
      endb = true;
    }
  }

  // Computing KS test's p-value
  KSstat = effna * effnb / effn * KSstat;

  return KSstat; //Kolmogorov-Smirnov Max Dist
}
