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
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector shift_bypidNum(const NumericVector& x, const int& lag,
                            const double& replace, const IntegerVector& id) {
  // id should be sorted and same length as x
  int n = x.size();
  NumericVector out(n);
    for (int i = 0; i < lag; i++)
    {
      out[i] = replace; // TODO FIXME!!
    }

    for (int i = lag; i < n; i++)
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    return out;
  // } else {
  //   stop("This function does not work because number of ids is too small!");
  // }
}

//' @export
// [[Rcpp::export]]
IntegerVector shift_bypidInt(const IntegerVector& x, const int& lag,
                            const int& replace, const IntegerVector& id) {
  // id should be sorted and same length as x
  int n = x.size();
  IntegerVector out(n);
  for (int i = 0; i < lag; i++)
  {
    out[i] = replace; // TODO FIXME!!
  }

    for (int i = lag; i < n; i++)
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    // SEXP cl = x.attr("class");
    // if (!Rf_isNull(cl) && Rf_length(cl) > 1) {
    //   out.attr("names") = VECTOR_ELT(dm, 1);
    // }
    if (x.hasAttribute("class")) {
      CharacterVector tt = x.attr("class");
      if (tt[0] == "factor") { // only works if factor is the first class
        out.attr("levels") = x.attr("levels");
        out.attr("class") = "factor";
      }
    }

    return out;
}


//' @export
// [[Rcpp::export]]
IntegerVector shift_bypidBool(const LogicalVector& x, const int& lag,
                             const bool& replace, const IntegerVector& id) {
  // id should be sorted and same length as x
  int n = x.size();
  IntegerVector out(n);
  for (int i = 0; i < lag; i++)
  {
    out[i] = replace; // TODO FIXME!!
  }

    for (int i = lag; i < n; i++)
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    return out;
}

//' @export
// [[Rcpp::export]]
StringVector  shift_bypidStr(const CharacterVector& x, const int& lag,
                             const std::string& replace, const IntegerVector& id) {
  // id should be sorted and same length as x
  int n = x.size();
  StringVector out(n);
  for (int i = 0; i < lag; i++)
  {
    out[i] = replace; // TODO FIXME!!
  }

    for (int i = lag; i < n; i++)
    {
      if (id[i] == id[i-lag]) {
        out[i] = x[i-lag];
      } else {
        out[i] = replace;
      }
    }
    return out;
}



