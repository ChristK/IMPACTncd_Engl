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
htn_model <-
  function(
    scenario_nam,
    mc_iter, # Not used currently, for consistency
    dt,
    design_,
    diagnosis_prb = 0.4,
    postHC_diagnosis_prb = 0.9, # Assumes almost perfect diagnosis of HTN after HC
    timing = TRUE) {
    message("Estimating HTN...")
    if (timing) ptm <- proc.time()

    if (!nzchar(scenario_nam)) { # first run for scenario ""
    set(dt, NULL, "prb_htn_incd", 0L)
    dt[sbp_curr_xps > 140 | bpmed_curr_xps > 0, prb_htn_incd := 1L]


    # Prevalence
    set(dt, NULL, "htn_prvl", 0L)
    dt[year == design_$sim_prm$init_year &
        (sbp_curr_xps > 140 | bpmed_curr_xps > 0),
      htn_prvl := 1L + rnbinom(.N, age, 0.8)] # Duration of disease is set arbitrarily (not used elsewhere)
    dt[htn_prvl > 0L & htn_prvl > (age - 20L), htn_prvl := age - 20L]

    # Diagnosis
    # hist(rnbinom(1e4, 1, diagnosis_prb), 100) # the distribution of the number of years until diagnosis
    set(dt, NULL, "prb_htn_dgn", diagnosis_prb)
    dt[bpmed_curr_xps > 0, prb_htn_dgn := 1]
    # dt[bpmed_curr_xps == 0, prb_htn_dgn := 0]

    set(dt, NULL, "htn_dgn", 0L)
    dt[htn_prvl > 0, htn_dgn := clamp(htn_prvl - 20L, 0, 100)]


    # dt[dead == F, prop_if(htn_prvl > 0 & htn_dgn > 0)]
    # dt[year == 15 & dead == F, prop_if(htn_prvl > 0 & htn_dgn == 0)]
    #
    # dt[htn_prvl > 0 & year == 15 & dead == F, prop_if(htn_dgn > 0)]
    } else { # subsequent runs for sc1 etc.
      scenario_nam <- "_sc"
      atte_colnam <- paste0("attendees", scenario_nam)
      colnam <- paste0("htn_prvl", scenario_nam)
      colnam_dgn <- paste0("htn_dgn", scenario_nam)
      colnam_bio <- paste0("sbp", scenario_nam)
      colnam_tr <- paste0("bpmed_px", scenario_nam)

      #incidence (not affected by HC)
      set(dt, NULL, paste0("prb_htn_incd", scenario_nam), dt$prb_htn_incd)

      # Prevalence
      set(dt, NULL, colnam, 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, (colnam) := htn_prvl]
      dt[year == design_$sim_prm$init_year_fromGUI & htn_prvl > 1L, (colnam) := htn_prvl]

      # Diagnosis
      set(dt, NULL, colnam_dgn, 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, (colnam_dgn) := htn_dgn]
      dt[year == design_$sim_prm$init_year_fromGUI & htn_dgn > 1L, (colnam_dgn) := htn_dgn]

      dt[, (paste0("prb_htn_dgn", scenario_nam)) := prb_htn_dgn]
      dt[htn_dgn == 0L & attendees_sc == 1L, paste0("prb_htn_dgn", scenario_nam) := postHC_diagnosis_prb]

    }
    if (timing) print(proc.time() - ptm)
}


