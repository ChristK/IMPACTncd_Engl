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

#ifndef aux_functions_H
#define aux_functions_H

Rcpp::IntegerVector carry_forward(Rcpp::IntegerVector& x,
                            const Rcpp::LogicalVector& pid_mrk,
                            const int& y,
                            const bool& byref = false);

Rcpp::IntegerVector carry_forward_incr(Rcpp::IntegerVector& x,
                                       const Rcpp::LogicalVector& pid_mrk,
                                 const bool& recur,
                                 const int& y = 1,
                                 const bool& byref = false);


Rcpp::IntegerVector carry_backward(const Rcpp::IntegerVector& x,
                                   const Rcpp::LogicalVector& pid_mrk,
                             const int& y = 0);

Rcpp::IntegerVector carry_backward_decr(const Rcpp::IntegerVector& x,
                                        const Rcpp::LogicalVector& pid_mrk);


Rcpp::LogicalVector mk_new_simulant_markers(const Rcpp::IntegerVector& pid);


Rcpp::LogicalVector identify_longdead(const Rcpp::IntegerVector& x,
                                      const Rcpp::LogicalVector& pid);


Rcpp::IntegerVector identify_invitees(const Rcpp::IntegerVector& elig,
                                const Rcpp::IntegerVector& prev_inv,
                                const Rcpp::NumericVector& prb,
                                const Rcpp::IntegerVector& freq,
                                const Rcpp::LogicalVector& pid);

Rcpp::IntegerVector hc_effect(const Rcpp::IntegerVector& x,
                             const double& prb_of_continuation,
                             const Rcpp::LogicalVector& pid);

Rcpp::NumericVector fbound(const Rcpp::NumericVector& x,
                           Rcpp::NumericVector& a,
                           Rcpp::NumericVector& b);

double antilogit(const double& x);


#endif
