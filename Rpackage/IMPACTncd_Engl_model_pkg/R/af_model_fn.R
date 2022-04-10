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
af_model <- function(
  scenario_nam,
  mc_iter, # Not currently used. For consistency
  dt,
  design_,
  diagnosis_prb = NA,
  postHC_diagnosis_prb = 0.1, # currently not checked in a HC. TODO link to GUI
  timing = TRUE) {
  message("Estimating AF...")
  if (timing) ptm <- proc.time()

  if (!nzchar(scenario_nam)) { # first run for scenario ""
    # Incidence
    dt[, prb_af_incd := clamp(af_prvl_curr_xps)]

    # Prevalence
    set(dt, NULL, "af_prvl", 0L)
    dt[year == design_$sim_prm$init_year &
        af_prvl_curr_xps > 0, af_prvl := af_prvl_curr_xps]

    # Diagnosis
    set(dt, NULL, "af_dgn", 0L)
    dt[year == design_$sim_prm$init_year, af_dgn := af_dgn_curr_xps]

    # dt[, prb_af_dgn := prop_if(af_prvl_curr_xps > 0 & af_dgn_curr_xps > 0)/prop_if(af_prvl_curr_xps > 0), keyby = .(age, sex, qimd)]
    dt[, prb_af_dgn := fifelse(af_dgn_curr_xps == 0, 0, 1)]
  } else { # subsequent runs for sc1 etc.
    # l <- fromGUI_scenario_parms(scenario_nam, parameters_dt) # TODO link with GUI
    scenario_nam <- "_sc"
    atte_colnam <- "attendees_sc"
    colnam      <- "af_prvl_sc"
    colnam_dgn  <- "af_dgn_sc"

    dt[, af_prvl_sc := af_prvl]

    # dt[year < design_$sim_prm$init_year_fromGUI, af_prvl_sc := af_prvl]
    # dt[year == design_$sim_prm$init_year_fromGUI & af_prvl > 1L, af_prvl_sc := af_prvl]

    dt[, af_dgn_sc := af_dgn]
    # dt[year < design_$sim_prm$init_year_fromGUI, af_dgn_sc := af_dgn]
    # dt[year == design_$sim_prm$init_year_fromGUI & af_dgn > 1L, af_dgn_sc := af_dgn]


    dt[, prb_af_incd_sc := prb_af_incd]
    dt[, prb_af_dgn_sc  := prb_af_dgn]

    dt[attendees_sc == 1 & af_prvl > 0 & af_dgn == 0,
       af_dgn_sc := rbinom(.N, 1, postHC_diagnosis_prb)] # TODO link with GUI
    dt[, af_dgn_sc := carry_forward(af_dgn_sc, pid_mrk, 1L)] # This replaces years of dgn with 1. Not a problem since it gets fixed when I run Peter's code
    dt[af_dgn_sc == 1L, prb_af_dgn_sc := 1]
  }

  # TODO consider linking with SBP from Okin et al Hypertension. 2015;66:368–373.Although link with SBP reduction in scenarios possibly better

  if (timing) print(proc.time() - ptm)
  # NOTE: other approaches with eval(quote()) have no speed benefit.
  # NOTE: there is a small benefit with ie
  # eval(str2lang("dt[attendees_sc1 == 1 & af_prvl > 0 & af_dgn == 0,.N]"))
  # but it does not justify the additional complexity
}

