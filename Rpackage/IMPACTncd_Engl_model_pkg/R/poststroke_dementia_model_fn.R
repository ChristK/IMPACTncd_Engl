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

#' @export
poststroke_dementia_model <-
  function(
    scenario_nam,
    mc,
    dt,
    design_,
    diagnosis_prb = 1,
    timing = TRUE)
  {
    message("Loading poststroke dementia model...")
    if (timing) ptm <- proc.time()

    # This is a special condition. It only depends on stroke and no other risk
    # factors. To fit this in the existing C++ framework I modified C++ code to
    # allow disease dependency only on incidence and not prevalence (incidence
    # type 5). The multiplier is 1 for the incident stroke year and 0 for every
    # other year.

    # The probability is based on Eithne's syst review
    if (!nzchar(scenario_nam)) { # first run for scenario ""
      tt <- RR$stroke_dementia$get_rr(mc, design_, drop = TRUE)
      setnames(tt, "stroke_rr", "prb_poststroke_dementia_incd_nostroke")
      absorb_dt(dt, tt)
      setnafill(dt, "c", 1, cols = "prb_poststroke_dementia_incd_nostroke")


      set(dt, NULL, "poststroke_dementia_incd_stroke_mltp", 1.0) # C++ converts this to 0 for non stroke incident cases

      set(dt, NULL, "prb_poststroke_dementia_dgn", diagnosis_prb)
      set(dt, NULL, "poststroke_dementia_dgn", 0L)
    } else {
      colnam <- "poststroke_dementia_prvl_sc"
      set(dt, NULL, colnam, 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, (colnam) := poststroke_dementia_prvl]
      dt[year == design_$sim_prm$init_year_fromGUI & poststroke_dementia_prvl > 1L,
         (colnam) := poststroke_dementia_prvl]

      colnam <- "prb_poststroke_dementia_incd_nostroke_sc"
      dt[, (colnam) := prb_poststroke_dementia_incd_nostroke]
      # tt <- RR$stroke_dementia$get_rr(mc, design_, drop = TRUE)
      # setnames(tt, "stroke_rr", colnam)
      # absorb_dt(dt, tt)
      # setnafill(dt, "c", 1, cols = colnam)

      colnam <- "prb_poststroke_dementia_dgn_sc"
      dt[, (colnam) := prb_poststroke_dementia_dgn]

      colnam <- "poststroke_dementia_dgn_sc"
      set(dt, NULL, colnam, 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, (colnam) := poststroke_dementia_dgn]
      dt[year == design_$sim_prm$init_year_fromGUI & poststroke_dementia_dgn > 1L,
         (colnam) := poststroke_dementia_dgn]

    }

  }
