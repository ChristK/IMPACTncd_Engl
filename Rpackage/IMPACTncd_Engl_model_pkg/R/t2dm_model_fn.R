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
t2dm_model <-
  function(
    scenario_nam,
    mc,
    dt,
    design_,
    diagnosis_prb = 0.6, # probability to diagnosis every year
    postHC_diagnosis_prb = 0.9, # 90% of T2DM after a HC is diagnosed
    timing = TRUE) {
    message("Estimating T2DM...")
    if (timing) ptm <- proc.time()

    if (!nzchar(scenario_nam)) { # first run for scenario ""
    # Prevalence
    dt[, t2dm_prvl := t2dm_prvl_curr_xps]

    # Diagnosis
    # hist(rnbinom(1e4, 1, diagnosis_prb), 100) # the distribution of the number of years until diagnosis
    set(dt, NULL, "prb_t2dm_dgn", diagnosis_prb)
    setnafill(dt, "c", 0L, cols = "t2dm_dgn")

    } else { # subsequent runs for sc1 etc.

      scenario_nam             <- "_sc"
      atte_colnam              <- "attendees_sc"
      colnam_prb_incd_nocvd    <- "prb_t2dm_incd_nocvd_sc"
      colnam_incd_cvd_mltp <- "t2dm_incd_cvd_mltp_sc" # suffix is an exception here
      colnam_dgn               <- "t2dm_dgn_sc"

      # incidence
      dt[, (c(colnam_prb_incd_nocvd, colnam_incd_cvd_mltp)) :=
          QDiabetes_vec_inputs(age, sex, cst_prvl,
            # fifelse(bpmed_px_sc == 1L | bpmed_curr_xps == 1L, 1L, 0L),
            bpmed_curr_xps,
            # above disabled as produced disproportionally high diabetes cases
            # and it is not necesarily causally related. Instead I explicitly
            # model the statin effect on diabetes, below

            bmi_sc, ethnicity, fam_t2dm, smoke_cat_sc, tds)]



      # OR for statins on t2dm from Sattar N, Preiss D, Murray HM, Welsh P,
      # Buckley BM, de Craen AJ, et al. Statins and risk of incident diabetes: a
      # collaborative meta-analysis of randomised statin trials. The Lancet
      # 2010;375:735–42. # They report OR, I assume RR is approximately equal.
      # Mean observation time was 4 years so I lag for 4 years
      tt <- RR$statins_t2dm$get_rr(mc, design_, drop = TRUE)

      # Lagged exposures
      exps_tolag <- c("statin_px_sc", "statin_px_curr_xps")
      exps_nam   <- paste0(exps_tolag, "_lagged")
      for (i in seq_along(exps_tolag)) {
        set(dt, NULL, exps_nam[i],
            dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$statin_t2dm_lag, pid, 0L)])
      }

      dt[statin_px_sc_lagged == 1L & statin_px_curr_xps_lagged == 0L,
         (c("prb_t2dm_incd_nocvd_sc", "t2dm_incd_cvd_mltp_sc")) :=
           .(tt * prb_t2dm_incd_nocvd_sc, tt * t2dm_incd_cvd_mltp_sc)]
      # no link to adherence because "Meta-regression showed that risk of
      # development of diabetes with statins was highest in trials with older
      # participants, but neither baseline body-mass index nor change in
      # LDL-cholesterol concentrations accounted for residual variation in
      # risk."
      dt[, (exps_nam) := NULL]
      dt[, (colnam_incd_cvd_mltp) := t2dm_incd_cvd_mltp_sc / prb_t2dm_incd_nocvd_sc]

      # prevalence
      set(dt, NULL, "t2dm_prvl_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, t2dm_prvl_sc := t2dm_prvl]
      dt[year == design_$sim_prm$init_year_fromGUI & t2dm_prvl > 0L, t2dm_prvl_sc := t2dm_prvl]

      # Diagnosis
      set(dt, NULL, colnam_dgn, 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, (colnam_dgn) := t2dm_dgn]
      dt[year == design_$sim_prm$init_year_fromGUI & t2dm_dgn > 1L, (colnam_dgn) := t2dm_dgn]


      dt[, (paste0("prb_t2dm_dgn", scenario_nam)) := prb_t2dm_dgn]
      dt[t2dm_dgn == 0L & attendees_sc == 1L, paste0("prb_t2dm_dgn", scenario_nam) := postHC_diagnosis_prb]
    }

    if (timing) print(proc.time() - ptm)
  }
