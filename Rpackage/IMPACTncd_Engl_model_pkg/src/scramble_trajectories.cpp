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
// #define BOOST_PENDING_INTEGER_LOG2_HPP // to prevent warning from BH. TODO remove this and next line when BH version >= 1.70
// #include <boost/integer/integer_log2.hpp>
// [[Rcpp::depends(dqrng)]]
// for an R package use LinkingTo: dqrng and remove above
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::plugins(openmp)]]
#include <omp.h>
#include <dqrng.h>

using namespace Rcpp;


// helper function that makes a double jump within x - jump and x + jump
// special care to avoid values outside (0, 1)
// seed and RNG type are defined in R side with dqset.seed(1234L); dqRNGkind("pcg64")
double fscramble_hlp (const double& x, const double& jump)
{
  // dqrng::dqRNGkind("pcg64");
  // dqrng::dqset_seed(IntegerVector::create(42));
  // if (jumpiness <= 0.0) stop("Jumpiness should be > 0"); //not necessary as I test the exported functions
  double out = dqrng::runif(x - jump, x + jump);
  if ((out <= 0.0) || (out >= 1.0)) out = x;

  return out;
}

//' @export
// [[Rcpp::export]]
NumericVector fscramble_trajectories (NumericVector& x, const LogicalVector& pid,
                                      const double& jumpiness = 1.0,
                                      const bool& inplace = true)
{
  // pid should be sorted and same length as x
  if (jumpiness <= 0.0) stop("Jumpiness should be positive");
  int n = x.size();
  NumericVector jump = dqrng::dqrexp(n, 50.0/jumpiness);  // mean jump = 0.02 with defaults

  if (inplace)
  {
    for (int i = 0; i < n; ++i)
      if (!pid[i]) x[i] = fscramble_hlp(x[i-1], jump[i]);
    return NULL;
  }
  else
  {
    NumericVector out(n);
    for (int i = 0; i < n; ++i)
    {
      if (pid[i]) out[i] = x[i];
        else out[i] = fscramble_hlp(out[i-1], jump[i]);
    }
  return out;
  }
}


// NumericVector fscramble_trajectories_inplace (NumericVector& x, const LogicalVector& pid, const double& jumpiness = 1.0) {
//   // pid should be sorted and same length as x
//   if (jumpiness <= 0.0) stop("Jumpiness should be > 0");
//     for (int i = 0; i < x.length(); ++i)
//     {
//       if (!pid[i]) x[i] = fscramble_hlp(x[i-1], jumpiness);
//     }
//     return NULL;
// }

/*** R
# dqset.seed(42)
# fscramble_hlp(0.8, 0.03)
# fscramble_hlp(0.8, 0.03)
*/
