## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

#' IMPACTncdEnglmisc: A package with functions for the IMPACTncdEngl model.
#'
#' Description of your package.
#
# @section Foo functions:
# The foo functions ...
#
#' @docType package
#' @author Chris Kypridemos
#' @import Rcpp R6 data.table CKutils
#' @importFrom R6 R6Class
#' @importFrom Rcpp evalCpp
#' @importFrom fst read_fst
#' @importFrom dqrng dqrunif dqsample
#' @importFrom mc2d qpert
#' @importFrom cowplot ggsave2
#' @importFrom gamlss fitDist fitDistPred predictAll
#' @importFrom stats as.formula na.omit qunif var weighted.mean loess predict qbinom
#' @importFrom stats quantile rbinom rpois runif sigma
#' @importFrom utils tail
#' @importFrom digest digest2int
#' @importFrom yaml read_yaml write_yaml
#' @importFrom foreach foreach
#' @useDynLib IMPACTncdEnglmisc
#' @name IMPACTncdEnglmisc
NULL

