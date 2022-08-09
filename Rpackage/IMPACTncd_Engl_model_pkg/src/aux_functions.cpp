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
#include "aux_functions.h"
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
IntegerVector carry_forward(IntegerVector& x,
                            const LogicalVector& pid_mrk,
                            const int& y,
                            const bool& byref)
{
  const int n = x.size();
  if (byref) // Alters x by reference
  {
    for (int i = 0; i < n; i++)
    {
      if (!pid_mrk[i] && x[i - 1] == y) x[i] = y;
    }
    return x;
  }
  else // Returns a new IntegerVector
  {
    IntegerVector out = clone(x);
    for (int i = 0; i < n; i++)
    {
      if (!pid_mrk[i] && out[i - 1] == y) out[i] = y;
    }
    return out;
  }

}

//' @export
// [[Rcpp::export]]
IntegerVector carry_forward_incr(IntegerVector& x, const LogicalVector& pid_mrk,
                                 const bool& recur, const int& y,
                                 const bool& byref) {
  // byref = true changes input x inplace
  // recur = false means that the value constantly increasing until it meets a new pid
  // recur = true means that the process restarts as soon as it finds a value < y

  const int n = x.size();
  if (recur)
  {
    if (byref) // Alters x by reference
    {
      for (int i = 0; i < n; i++)
      {
        if (!pid_mrk[i] && x[i] >= y && x[i - 1] >= y) x[i] = x[i - 1] + 1;
      }
      return x;
    } // end byref = true
    else // Returns a new IntegerVector
    {
      IntegerVector out = clone(x);

      for (int i = 0; i < n; i++)
      {
        if (!pid_mrk[i] && out[i] >= y && out[i - 1] >= y) out[i] = out[i - 1] + 1;
      }
      return out;
    } // end byref = false
  } // end if  recur
  else
  { // if not recur
    if (byref) // Alters x by reference
    {
      for (int i = 0; i < n; i++)
      {
        if (!pid_mrk[i] && x[i - 1] >= y) x[i] = x[i - 1] + 1;
      }
      return x;
    } // end byref = true
    else // Returns a new IntegerVector
    {
      IntegerVector out = clone(x);

      for (int i = 0; i < n; i++)
      {
        if (!pid_mrk[i] && out[i - 1] >= y) out[i] = out[i - 1] + 1;
      }
      return out;
    } // end byref = false
  } // end if not recur

  }



//' @export
// [[Rcpp::export]]
IntegerVector carry_backward(const IntegerVector& x, const LogicalVector& pid_mrk,
                             const int& y) {
  const int n = x.size();
  IntegerVector out = clone(x);
  for (int i = n - 1; i > 0; i--) // Go backwards but stop from one row before the last
  {
    if (!pid_mrk[i] && out[i] > y) out[i - 1] = out[i] - 1;
    if (i < (n - 1) && pid_mrk[i] && !pid_mrk[i + 1] && out[i + 1] > y && out[i] == y) out[i] = out[i + 1] - 1;
    if (out[i] < 0) out[i] = 0;
  }
  return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector carry_backward_decr(const IntegerVector& x, const LogicalVector& pid_mrk) {
  const int n = x.size();
  IntegerVector out = clone(x);
  for (int i = n - 1; i > 0; i--) // Go backwards but stop from one row before the last
  {
    if (!pid_mrk[i] && out[i] > 0) out[i - 1] = out[i] - 1;
    if (i < (n - 1) && pid_mrk[i] && !pid_mrk[i + 1] && out[i + 1] > 0 && out[i] == 0) out[i] = out[i + 1] - 1;
    if (out[i] < 0) out[i] = 0;
  }
  return out;
}

//' @export
// [[Rcpp::export]]
LogicalVector mk_new_simulant_markers(const IntegerVector& pid)
{
  // pid should be sorted and same length as x
  const int n = pid.size();
  LogicalVector new_simulant_markers(n);
  new_simulant_markers[0] = true;
  int previous_pid = pid[0];
  // Loop with no conditional branches in the body (therefore branch predictor should get it right almost every time) and minimal memory access by retaining previous_pid.
  for (int i = 1; i < n; i++)
  {
    new_simulant_markers[i] = pid[i] != previous_pid;
    previous_pid = pid[i];
  }
  return new_simulant_markers;
}

//' @export
// [[Rcpp::export]]
LogicalVector identify_longdead(const IntegerVector& x, const LogicalVector& pid) {
  const int n = x.size();
  LogicalVector out(n);
  for (int i = 0; i < n; i++)
  {
    if (!pid[i] && x[i-1] != 0) out[i] = true; else out[i] = false;
  }
  return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector identify_invitees(const IntegerVector& elig,
                                const IntegerVector& prev_inv,
                                const NumericVector& prb,
                                const IntegerVector& freq,
                                const LogicalVector& pid
) {
  // prev_inv only relevant for serial ensembles. Signifies whether invited in previous scenario
  const int n = elig.size();
  int counter = 1000; // to ensure > max(freq)
  // int counter1 = 1;
  IntegerVector out(n, 0);
  for (int i = 0; i < n; i++)
  {
    if (pid[i]) counter = 1000; //how many years since previous invitation?
    if (prev_inv[i] == 1) counter = 0;
    else if (counter < 1000) counter++; // counter == 2, 2 year has passed
    if (elig[i] == 1 && counter >= freq[i]) out[i] = R::rbinom(1.0, prb[i]);
    if (out[i] == 1) counter = 0;
    // {
    //   counter1 = 1;
    //   while (counter1 <= i && counter1 < freq[i] && !pid[i-counter1] && out[i] == 1)
    //   {
    //     if (out[i-counter1] == 1) out[i] = 0; // invited recently
    //     counter1++;
    //   }
    // }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
IntegerVector hc_effect(const IntegerVector& x,
                             const double& prb_of_continuation,
                             const LogicalVector& pid)
  {
  const int n = x.size();
  IntegerVector out = clone(x);
  for (int i = 0; i < n; i++)
  {
    if (!pid[i] && out[i-1] == 1) out[i] = R::rbinom(1, prb_of_continuation);
  }
  return out;
}


//' @export
// [[Rcpp::export]]
NumericVector fbound(const NumericVector &x, NumericVector &a, NumericVector &b) {
  const int n = x.size();
  NumericVector out(n);
  for(int i = 0; i < n; i++) {
    if (NumericVector::is_na(x[i])) out[i] = NA_REAL;
    else
    {
      if (a[i] > b[i]) {double c = a[i]; a[i] = b[i]; b[i] = c;}; // ensure a < b

      if (x[i] < a[i]) out[i] = a[i];
      else if (x[i] > b[i]) out[i] = b[i];
      else out[i] = x[i];
    }
  }
  return out;
}

//' @export
// [[Rcpp::export]]
double antilogit(const double& x)
{
  return exp(x)/(1 + exp(x));
}
