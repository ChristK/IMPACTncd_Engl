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
stroke_model <-
  function(
    scenario_nam,
    mc,
    dt,
    design_,
    diagnosis_prb = 1,
    timing = TRUE) {

    message("Loading stroke (I60-I69) model...")
    if (timing) ptm <- proc.time()

    if (!nzchar(scenario_nam)) { # first run for scenario ""
    # Lagged exposures
    exps_tolag <-
      grep(
        "^act|^fru|^veg|^smok_st|^ets|^alco|^bmi|^sbp|^tchol_c|^t2dm_prvl_c|^af_prvl_c|^af_dgn_c",
        names(dt),
        value = TRUE
      )
    exps_nam <-  gsub("_curr_xps$", "_lagged", exps_tolag)
    for (i in seq_along(exps_tolag)) {
      set(dt, NULL, exps_nam[i], dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$cvd_lag, pid)])
    }
    } else {
      exps_tolag <-
        c(paste0(c("active_days_", "alcohol_", "bmi_",
          "sbp_", "tchol_", "smok_status_", "t2dm_prvl_"), "sc"),
          "ets_curr_xps", "fruit_curr_xps", "veg_curr_xps", "af_prvl_curr_xps", "af_dgn_curr_xps")
      # TODO AF_prvl_sc & af_dgn_sc influence stroke
      # TODO expand C++ code to allow influence from multiple disease and diagnosis

      exps_nam <-  gsub("_curr_xps$|_sc$", "_lagged", exps_tolag)
      for (i in seq_along(exps_tolag)) {
        set(dt, NULL, exps_nam[i], dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$cvd_lag, pid)])
      }
    }
    ## RR for AF from (Christiansen et al. 2016). For no risk factor and many limitations
    tt <-
      RR$af_stroke$get_rr(mc, design_, drop = TRUE)[, i.af_prvl_lagged := 1L]
    dt[tt, on = .(age, af_prvl_lagged >= i.af_prvl_lagged), af_rr := i.af_rr]
    setnafill(dt, "c", 1, cols = "af_rr")
    # dt[, summary(af_rr)]

    # From ATRIA study, assume diagnosed cases on warfarin reduce risk by ~66%
    # Further assume that ~80% of those who should be on warfarin are on it
    # TODO find prevalence of warfarin and newer anticoag
    # dt[af_dgn_lagged > 0L &
    #     rbinom(.N, 1, 0.8) == 1L, af_rr := 1 + ((af_rr - 1) * 0.34)]
    # TODO effect of health check


    # RR for tobacco from Ezzati M, Henley SJ, Thun MJ, Lopez AD. Role of Smoking in Global and Regional
    # Cardiovascular Mortality. Circulation. 2005 Jul 26;112(4):489–97.
    # Table 1 Model B

    #ex-smokers
    # Cigarette smoking as a risk factor for stroke: The Framingham study
    # "Stroke risk decreased significantly by two years and was at the
    # level of nonsmokers by five years after cessation of cigarette smoking"
    #cat("smoking RR\n")
    tt <- RR$tobacco_stroke$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "smok_status", "smok_status_lagged")
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "tobacco_rr")
    # dt[, summary(tobacco_rr)]

    # Calculate PAF of ETS for stroke
    # RR from Oono IP, Mackay DF, Pell JP. Meta-analysis of the association between secondhand smoke exposure and stroke.
    # J Public Health 2011;33:496–502. doi:10.1093/pubmed/fdr025
    tt <-
      data.table(
        ets_lagged = 1L,
        smok_status_lagged = as.character(1:3),
        ets_rr = RR$ets_stroke$get_rr(mc, design_, drop = TRUE)
      )
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "ets_rr")
    # dt[, prop_if(ets_rr > 1), keyby = smok_status_lagged]

    # Calculate RR for stroke. Optimal SBP level at 115mmHg and RR(HR) of dying from
    # stroke was taken from "Age-specific relevance of usual blood pressure to
    # vascular mortality: a meta-analysis of individual data for one million adults in 61 prospective studies.
    # The Lancet. 2002 Dec 14;360(9349):1903–1913"
    # Figure 3
    #cat("sbp RR\n")
    tt <- RR$sbp_stroke$get_rr(mc, design_, drop = TRUE)
    absorb_dt(dt, tt)
    dt[, sbp_rr := clamp(sbp_rr ^ ((
      RR$sbptmred_stroke$get_rr(mc, design_, drop = TRUE) - sbp_lagged
    ) / 20), 1, 20)]
    setnafill(dt, "c", 1, cols = "sbp_rr")
    # dt[, summary(sbp_rr)]


    # Calculate RR for stroke. Optimal chol level at 3.8 mmol/L and RR(HR) of
    # dying from stroke was taken from "Blood cholesterol and
    # vascular mortality by age, sex, and blood pressure: a meta-analysis of
    # individual data from 61 prospective studies
    # with 55.000 vascular deaths. The Lancet. 2007;370:1829–39.
    # Figure 4 (for total stroke. I used only significant HR's).
    #cat("chol RR\n")
    tt <- RR$tchol_stroke$get_rr(mc, design_, drop = TRUE)
    absorb_dt(dt, tt)
    dt[, tchol_rr := clamp(tchol_rr ^ (3.8 - tchol_lagged), 1, 20)]
    setnafill(dt, "c", 1, cols = "tchol_rr")
    # dt[, summary(tchol_rr)]


    # RR for BMI from "The Emerging Risk Factors Collaboration.
    # Separate and combined associations of body-mass index and abdominal adiposity
    # with cardiovascular disease: collaborative analysis of 58 prospective studies.
    # The Lancet 2011;377:1085–95. doi:10.1016/S0140-6736(11)60105-0
    # Table 1 (Adjusted for age, sex, smoking status, systolic blood pressure,
    # history of diabetes, and total and HDL cholesterol)
    # BMI not significant for ischaemic stroke but other obesity metrics are.
    #!! NEED TO decide if I want to use it
    #cat("bmi RR\n")
    tt <- RR$bmi_stroke$get_rr(mc, design_, drop = TRUE)
    absorb_dt(dt, tt)
    dt[, bmi_rr := clamp(bmi_rr ^ ((bmi_lagged - 20) / 4.56), 1, 20)]
    setnafill(dt, "c", 1, cols = "bmi_rr")
    # dt[, summary(bmi_rr)]

    # RR for diabetes from The Emerging Risk Factors Collaboration.
    # Diabetes mellitus, fasting blood glucose concentration,
    # and risk of vascular disease: a collaborative
    # meta-analysis of 102 prospective studies. The Lancet 2010;375:2215–22
    # figure 2 (HRs were adjusted for age, smoking status, body-mass index,
    # and  systolic blood pressure)
    #cat("diab RR\n")
    tt <- RR$t2dm_stroke$get_rr(mc, design_, drop = TRUE)
    tt[, i.t2dm_prvl_lagged := 1L]
    dt[tt, on = .(age, t2dm_prvl_lagged >= i.t2dm_prvl_lagged), t2dm_parf_rr := i.t2dm_rr]
    setnafill(dt, "c", 1, cols = "t2dm_parf_rr")
    # dt[, summary(t2dm_parf_rr)]

    # multiplier for the risk of stroke for diabetic
    dt[tt, on = .(age), stroke_incd_t2dm_mltp := i.t2dm_rr]
    setnafill(dt, "c", 1, cols = "stroke_incd_t2dm_mltp")

    # dt[, summary(prb_stroke_incd_t2dm_mltp)]

    # Calculate RR for stroke. From Dauchet L, Amouyel P, Dallongeville J.
    # Fruit and vegetable consumption and risk of stroke A meta-analysis of cohort
    # studies. Neurology. 2005;65:1193–7.
    # To avoid negative PAF an optimal level of F&V has to be set arbitrarily. I set it to 7
    #cat("fv RR\n")
    dt[,
      fv_rr := clamp(RR$fv_stroke$get_rr(mc, design_, drop = TRUE) ^
          ((fruit_lagged + veg_lagged - 80 * 7) / 80),
        1,
        20)] # no benefit after 7 portions
    setnafill(dt, "c", 1, cols = "fv_rr")
    # dt[, summary(fv_rr)]

    # RR for PA 1. WHO | Comparative Quantification of Health Risks [Internet].
    # WHO [cited 2014 Jan 30];Available from: http://www.who.int/publications/cra/en/
    # Table 10.20 (with adjustment for measurement error)
    tt <- RR$pa_stroke$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "active_days", "active_days_lagged")
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "pa_rr")
    # dt[, summary(pa_rr)]

    # From GBD 2017
    # Highest alcohol intake 72g/d
    tt <- RR$alcohol_stroke$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "alcohol", "alcohol_lagged")
    dt[alcohol_lagged > max(tt$alcohol_lagged), alcohol_lagged := max(tt$alcohol_lagged)]
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "alcohol_rr")
    # dt[, summary(alcohol_rr)]

    # ethnicity not important

    nam <- grep("_rr$", names(dt), value = TRUE)
    invisible(dt[, lapply(.SD, clamp, 0, 20, TRUE), .SDcols = nam]) # No rr > 20
    # dt[, lapply(.SD, max), .SDcols = nam]
        dt[, (exps_nam) := NULL]

    # Estimate PARF ------------------------------------------------------------
    #cat("Estimating stroke PAF...\n")
    if (!"p0_stroke" %in% names(dt)) {
      strokeparf <-
        dt[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH) &
            stroke_prvl == 0 & year == design_$sim_prm$init_year,
          .(parf = 1 - 1 / (
            sum(
              tobacco_rr * ets_rr * sbp_rr * tchol_rr * bmi_rr *
                t2dm_parf_rr * fv_rr * pa_rr * af_rr * alcohol_rr
            ) / .N
          )),
          keyby = .(age, sex, qimd)]
      # strokeparf[, parf := clamp(predict(loess(parf ~ age, span = 0.5))), by = .(sex, qimd)]
      # strokeparf[, {
      #   plot(age, parf, main = paste0(.BY[[1]],"-", .BY[[2]]), ylim = c(0, 1))
      #   lines(age, parf2)
      # }
      # , keyby = .(sex, qimd)]

      absorb_dt(strokeparf,
        get_disease_epi_mc(mc, "stroke", "i", "v", design_$sim_prm$stochastic))
      strokeparf[, p0_stroke := incidence * (1 - parf)]
      # strokeparf[, summary(p0_stroke)]
      strokeparf[is.na(p0_stroke), p0_stroke := incidence]
      strokeparf[, c("incidence", "parf") := NULL]
      absorb_dt(dt, strokeparf)
      setnafill(dt, "c", 0, cols = "p0_stroke")
      rm(strokeparf)
    }
    if (!nzchar(scenario_nam)) { # first run for scenario ""
    # Estimate stroke incidence prbl -------------------------------
    #cat("Estimating stroke incidence without diabetes...\n\n")
    set(dt, NULL, "prb_stroke_incd_not2dm", 0)
    dt[, prb_stroke_incd_not2dm :=
        p0_stroke * tobacco_rr * ets_rr * sbp_rr * tchol_rr * bmi_rr *
        fv_rr * pa_rr * af_rr * alcohol_rr] # Remeber no t2dm as this will be a multiplier
    dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

    # Assume a probability of diagnosis ----
    # hist(rnbinom(1e4, 1, diagnosis_prb), 100) # the distribution of the number of years until diagnosis
    set(dt, NULL, "prb_stroke_dgn", diagnosis_prb)
    set(dt, NULL, "stroke_dgn", 0L)
    dt[stroke_prvl > 0, stroke_dgn := clamp(stroke_prvl - 5L, 0, 100)]

    # Estimate case fatality ----
    absorb_dt(dt, get_disease_epi_mc(mc, "stroke", "f", "v", design_$sim_prm$stochastic))
    setnames(dt, "fatality", "prb_stroke_mrtl")

    } else {

      set(dt, NULL, "stroke_prvl_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, stroke_prvl_sc := stroke_prvl]
      dt[year == design_$sim_prm$init_year_fromGUI & stroke_prvl > 1L, stroke_prvl_sc := stroke_prvl]

      # Estimate stroke incidence prbl
      #cat("Estimating stroke incidence without diabetes...\n\n")
      colnam <-"prb_stroke_incd_not2dm_sc"
      set(dt, NULL, colnam, 0)
      dt[, (colnam) :=
          p0_stroke * tobacco_rr * ets_rr * sbp_rr * tchol_rr * bmi_rr *
          fv_rr * pa_rr * af_rr * alcohol_rr] # Remember no t2dm as this will be a multiplier
      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis
      dt[, prb_stroke_dgn_sc := prb_stroke_dgn]

      set(dt, NULL, "stroke_dgn_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, stroke_dgn_sc := stroke_dgn]
      dt[year == design_$sim_prm$init_year_fromGUI & stroke_dgn > 1L, stroke_dgn_sc := stroke_dgn]

      # Estimate case fatality
      if (!"prb_stroke_mrtl_sc" %in% names(dt))
        set(dt, NULL, "prb_stroke_mrtl_sc", dt$prb_stroke_mrtl)
    }

    if (timing) print(proc.time() - ptm)
  }

