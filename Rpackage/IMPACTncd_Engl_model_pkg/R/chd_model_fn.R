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
chd_model <-
  function(scenario_nam,
    mc,
    dt,
    design_,
    diagnosis_prb = 1,
    timing = TRUE) {
    message("Loading CHD (I20-I25) model...")
    if (timing)
      ptm <- proc.time()

    if (!nzchar(scenario_nam)) {
      # first run for scenario ""
      # Lagged exposures
      exps_tolag <-
        grep(
          "^act|^fru|^veg|^smok_st|^ets|^alco|^bmi|^sbp|^tchol_c|^t2dm_prvl_c",
          names(dt),
          value = TRUE
        )
      exps_nam <-  gsub("_curr_xps$", "_lagged", exps_tolag)
      for (i in seq_along(exps_tolag)) {
        set(dt, NULL, exps_nam[i],
          dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$cvd_lag, pid)])
      }

    } else {
      exps_tolag <-
        c(paste0(
          c(
            "active_days_",
            "alcohol_",
            "bmi_",
            "sbp_",
            "tchol_",
            "smok_status_",
            "t2dm_prvl_"
          ),
          "sc"
        ),
          "ets_curr_xps",
          "fruit_curr_xps",
          "veg_curr_xps")

      exps_nam <-
        gsub("_curr_xps$|_sc$",
          "_lagged",
          exps_tolag)
      for (i in seq_along(exps_tolag)) {
        set(dt, NULL, exps_nam[i], dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$cvd_lag, pid)])
      }
    }

    # RR for tobacco from Ezzati M, Henley SJ, Thun MJ, Lopez AD. Role of
    # Smoking in Global and Regional Cardiovascular Mortality. Circulation. 2005
    # Jul 26;112(4):489–97. Table 1 Model B

    #RR for ex-smokers from Huxley RR, Woodward M. Cigarette smoking as a risk
    #factor for coronary heart disease in women compared with men: a systematic
    #review and meta-analysis of prospective cohort studies. The Lancet. 2011
    #Oct 14;378(9799):1297–305. Appendix webfigure 8 cat("smoking RR\n")
    tt <- RR$tobacco_chd$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "smok_status", "smok_status_lagged")
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "tobacco_rr")
    # dt[, summary(tobacco_rr)]

    # RR for ETS He J, Vupputuri S, Allen K, Prerost MR, Hughes J, Whelton PK.
    # Passive Smoking and the Risk of Coronary Heart Disease — A Meta-Analysis
    # of Epidemiologic Studies. New England Journal of Medicine.
    # 1999;340(12):920–6. Table 3. Adjusted RR
    tt <-
      data.table(
        ets_lagged = 1L,
        smok_status_lagged = as.character(1:3),
        ets_rr = RR$ets_chd$get_rr(mc, design_, drop = TRUE)
      )
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "ets_rr")
    # dt[, prop_if(ets_rr > 1), keyby = smok_status_lagged]
    # dt[, summary(ets_rr)]

    # RR for SBP from Optimal SBP level at 115mmHg and RR(HR) of dying from CHD
    # was taken from "Age-specific relevance of usual blood pressure to vascular
    # mortality: a meta-analysis of individual data for one million adults in 61
    # prospective studies. The Lancet. 2002 Dec 14;360(9349):1903–1913" Figure 5
    #cat("sbp RR\n")
    tt <- RR$sbp_chd$get_rr(mc, design_, drop = TRUE)
    absorb_dt(dt, tt)
    dt[, sbp_rr := clamp(sbp_rr ^ ((
      RR$sbptmred_chd$get_rr(mc, design_, drop = TRUE) - sbp_lagged
    ) / 20), 1, 20)]
    setnafill(dt, "c", 1, cols = "sbp_rr")
    # dt[, summary(sbp_rr)]


    # RR for Chol from "Blood cholesterol and vascular mortality by age, sex,
    # and blood pressure: a meta-analysis of individual data from 61 prospective
    # studies with 55.000 vascular deaths. The Lancet. 2007 Dec
    # 7;370(9602):1829–39. Appendix Webtable 6  fully adjusted
    #cat("chol RR\n")
    tt <- RR$tchol_chd$get_rr(mc, design_, drop = TRUE)
    absorb_dt(dt, tt)
    dt[, tchol_rr := clamp(tchol_rr ^ (3.8 - tchol_lagged), 1, 20)]
    setnafill(dt, "c", 1, cols = "tchol_rr")
    # dt[, summary(tchol_rr)]


    # RR for BMI from "The Emerging Risk Factors Collaboration. Separate and
    # combined associations of body-mass index and abdominal adiposity with
    # cardiovascular disease: collaborative analysis of 58 prospective studies.
    # The Lancet 2011;377:1085–95. doi:10.1016/S0140-6736(11)60105-0 Table 1
    # (Adjusted for age, sex, smoking status, systolic blood pressure, history
    # of diabetes, and total and HDL cholesterol) and figure 2 for age specific
    # gradient
    #cat("bmi RR\n")
    tt <- RR$bmi_chd$get_rr(mc, design_, drop = TRUE)
    absorb_dt(dt, tt)
    dt[, bmi_rr := clamp(bmi_rr ^ ((bmi_lagged - 20) / 4.56), 1, 20)]
    setnafill(dt, "c", 1, cols = "bmi_rr")
    # dt[, summary(bmi_rr)]

    # RR for diab from The Emerging Risk Factors Collaboration. Diabetes
    # mellitus, fasting blood glucose concentration, and risk of vascular
    # disease: a collaborative meta-analysis of 102 prospective studies. The
    # Lancet 2010;375:2215–22. doi:10.1016/S0140-6736(10)60484-9 figure 2 (HRs
    # were adjusted for age, smoking status, body-mass index, and  systolic
    # blood pressure)
    #cat("diab RR\n")
    tt <- RR$t2dm_chd$get_rr(mc, design_, drop = TRUE)
    tt[, i.t2dm_prvl_lagged := 1L]
    dt[tt, on = .(age, t2dm_prvl_lagged >= i.t2dm_prvl_lagged), t2dm_parf_rr := i.t2dm_rr]
    setnafill(dt, "c", 1, cols = "t2dm_parf_rr")
    # dt[, summary(t2dm_parf_rr)]

    # multiplier for the risk of CHD for diabetic
    dt[tt, on = .(age), chd_incd_t2dm_mltp := i.t2dm_rr]
    setnafill(dt, "c", 1, cols = "chd_incd_t2dm_mltp")

    # RR for F&V from From Dauchet L, Amouyel P, Hercberg S, Dallongeville J.
    # Fruit and Vegetable Consumption and Risk of Coronary Heart Disease: A
    # Meta-Analysis of Cohort Studies. J Nutr. 2006 Oct 1;136(10):2588–93. To
    # avoid negative PAF an optimal level of F&V has to be set arbitrarily. I
    # set it to 7 when convert porftvg from categorical to numeric I create
    # bias. eg 1=less than 1 portion
    #cat("fv RR\n")
    dt[, fv_rr := clamp(RR$fv_chd$get_rr(mc, design_, drop = TRUE) ^
        ((fruit_lagged + veg_lagged - 80 * 7) / 80),
      1,
      30)] # no benefit after 7 portions
    setnafill(dt, "c", 1, cols = "fv_rr")
    # dt[, summary(fv_rr)]

    # RR for PA 1. WHO | Comparative Quantification of Health Risks [Internet].
    # WHO [cited 2014 Jan 30];Available from: http://www.who.int/publications/cra/en/
    # Table 10.19 (with adjustment for measurement error)
    #cat("pa RR\n")
    tt <- RR$pa_chd$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "active_days", "active_days_lagged")
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "pa_rr")
    # dt[, summary(pa_rr)]

    # RR for ethnicity 1. George J, Mathur R, Shah AD, Pujades-Rodriguez M,
    # Denaxas S, Smeeth L, et al. Ethnicity and the first diagnosis of a wide
    # range of cardiovascular diseases: Associations in a linked electronic
    # health record cohort of 1 million patients. PLOS ONE. 2017 Jun
    # 9;12(6):e0178945. Fig 3
    #cat("ethnicity RR\n")
    tt <- RR$ethnicity_chd$get_rr(mc, design_, drop = TRUE)
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "ethnicity_rr")
    # dt[, summary(ethnicity_rr)]

    # From GBD 2017
    # Highest alcohol intake 72g/d
    tt <- RR$alcohol_chd$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "alcohol", "alcohol_lagged")
    dt[alcohol_lagged > max(tt$alcohol_lagged),
      alcohol_lagged := max(tt$alcohol_lagged)]
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "alcohol_rr")
    # dt[, summary(alcohol_rr)]

    nam <- grep("_rr$", names(dt), value = TRUE)
    invisible(dt[, lapply(.SD, clamp, 0, 20, TRUE), .SDcols = nam]) # No rr > 20
    # dt[, lapply(.SD, max), .SDcols = nam]

    dt[, (exps_nam) := NULL]

    # Estimate PARF ------------------------------------------------------------
    #cat("Estimating CHD PAF...\n")
    # TODO this formula is accurate and meaningful only when RR >= 1. Otherwise
    # may estimate -ve PARF. I need to rescale all RR to >= 1 (i.e. ethnicity,
    # alcohol)
    if (!"p0_chd" %in% names(dt)) {
      chdparf <-
        dt[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH) &
            chd_prvl == 0 & year == design_$sim_prm$init_year,
          .(parf = 1 - 1 / (
            sum(
              tobacco_rr * ets_rr * sbp_rr * tchol_rr * bmi_rr *
                t2dm_parf_rr * fv_rr * pa_rr * ethnicity_rr * alcohol_rr
            ) / .N
          )),
          keyby = .(age, sex, qimd)]
      # chdparf[, parf := clamp(predict(loess(parf ~ age, span = 0.5))), by = .(sex, qimd)]
      # chdparf[, {
      #   plot(age, parf, main = paste0(.BY[[1]],"-", .BY[[2]]), ylim = c(0, 1))
      #   lines(age, parf2)
      # }
      # , keyby = .(sex, qimd)]

      absorb_dt(chdparf,
        get_disease_epi_mc(mc, "chd", "i", "v", design_$sim_prm$stochastic))
      chdparf[, p0_chd := incidence * (1 - parf)]
      # chdparf[, summary(p0_chd)]
      chdparf[is.na(p0_chd), p0_chd := incidence]
      chdparf[, c("incidence", "parf") := NULL]
      absorb_dt(dt, chdparf)
      setnafill(dt, "c", 0, cols = "p0_chd")
      rm(chdparf)
    }
    if (!nzchar(scenario_nam)) {
      # first run for scenario ""
      # Estimate CHD incidence prbl -------------------------------
      #cat("Estimating CHD incidence without diabetes...\n\n")
      set(dt, NULL, "prb_chd_incd_not2dm", 0)
      dt[, prb_chd_incd_not2dm :=
          p0_chd * tobacco_rr * ets_rr * sbp_rr * tchol_rr * bmi_rr *
          fv_rr * pa_rr * ethnicity_rr * alcohol_rr] # Remember no t2dm as this will be a multiplier
      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis ----
      # hist(rnbinom(1e4, 1, diagnosis_prb), 100) # the distribution of the number
      # of years until diagnosis
      set(dt, NULL, "prb_chd_dgn", diagnosis_prb)

      set(dt, NULL, "chd_dgn", 0L)
      dt[chd_prvl > 0, chd_dgn := clamp(chd_prvl - 5L, 0, 100)]

      # Estimate case fatality ----
      absorb_dt(dt,
        get_disease_epi_mc(mc, "chd", "f", "v", design_$sim_prm$stochastic))
      setnames(dt, "fatality", "prb_chd_mrtl")

    } else {
      set(dt, NULL, "chd_prvl_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, chd_prvl_sc := chd_prvl]
      dt[year == design_$sim_prm$init_year_fromGUI &
          chd_prvl > 1L, chd_prvl_sc := chd_prvl]

      # Estimate CHD incidence prbl
      #cat("Estimating CHD incidence without diabetes...\n\n")
      colnam <- "prb_chd_incd_not2dm_sc"
      set(dt, NULL, colnam, 0)
      dt[, (colnam) :=
          p0_chd * tobacco_rr * ets_rr * sbp_rr * tchol_rr * bmi_rr *
          fv_rr * pa_rr * ethnicity_rr * alcohol_rr] # Remember no t2dm as this will be a multiplier
      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis
      dt[, prb_chd_dgn_sc := prb_chd_dgn]

      set(dt, NULL, "chd_dgn_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, chd_dgn_sc := chd_dgn]
      dt[year == design_$sim_prm$init_year_fromGUI &
          chd_dgn > 1L, chd_dgn_sc := chd_dgn]

      # Estimate case fatality
      if (!"prb_chd_mrtl_sc" %in% names(dt))
        set(dt, NULL, "prb_chd_mrtl_sc", dt$prb_chd_mrtl)
    }

    if (timing)
      print(proc.time() - ptm)
  }
