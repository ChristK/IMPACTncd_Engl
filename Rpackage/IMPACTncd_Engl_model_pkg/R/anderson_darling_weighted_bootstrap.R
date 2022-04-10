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

# loosely based on nsRFA::ADbootstrap.test
#' @export
wtd_ADtest <- function(x, cod, wt, Nsim = 500) {
  dt1 <- data.table(x1 = x, cod1 = cod, wt1 = wt)
  dt2 <- split(dt1, by = "cod1", keep.by = FALSE)
  invisible(lapply(dt2, setkeyv, "x1"))
  A2kN <- wtd_ADstat(dt2[[1]]$x1, dt2[[1]]$wt1,
                     dt2[[2]]$x1, dt2[[2]]$wt1)

  A2kNs <- vector("numeric", Nsim)
  for (i in 1:Nsim) {
    dt2 <- split(dt1[dqsample(.N, sum(wt1), TRUE)], by = "cod1", keep.by = FALSE)
    invisible(lapply(dt2, setkeyv, "x1"))
    A2kNs[i] <-  wtd_ADstat(dt2[[1]]$x1, dt2[[1]]$wt1,
                            dt2[[2]]$x1, dt2[[2]]$wt1)
  }
  ecdfA2kNs <- ecdf(A2kNs)
  probabilita <- ecdfA2kNs(A2kN) # of 2 samples from same population. If < 0.05 then different
  output <- c(A2kN, probabilita)
  names(output) <- c("wtdADstat", "P")
  output
}
