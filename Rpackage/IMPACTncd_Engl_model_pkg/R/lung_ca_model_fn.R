## IMPACTncdEngl model: a decision support tool for primary prevention of NCDs
## Copyright (C) 2019 Chris Kypridemos

## IMPACTncdEngl model is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>
## or write to the Free Software Foundation, Inc., 51 Franklin Street,
## Fifth Floor, Boston, MA 02110-1301  USA.

#' @export
lung_ca_model <-
  function(scenario_nam,
    mc,
    dt,
    design_,
    diagnosis_prb = 0.7,
    timing = TRUE) {
    # TODO also associated with t2dm. However the current C++ framework does not allow
    # dependency on 2 diseases (copd and t2dm). A work around could be to use
    # group copd and t2dm together as I do with cvd but then only one multiplier
    # would be allowed

    message("Loading lung cancer (C33-C34) model...")
    if (timing)
      ptm <- proc.time()

    set(dt, NULL, "copd_forparf", 0L)
    if (!"p0_lung_ca" %in% names(dt)) {
      # I need copd prevalence estimate for PARF. I treat copd appropriately for all other calculations
      dt[year == design_$sim_prm$init_year, copd_forparf := fifelse(copd_prvl == 0L, 0L, 1L)]
    }

    if (!nzchar(scenario_nam)) {
      # first run for scenario ""
      # Lagged exposures
      exps_tolag <- grep("^smok_.*_curr_xps$|^ets|^fruit|^bmi",
        names(dt), value = TRUE)

      exps_nam <-  gsub("_curr_xps$", "_lagged", exps_tolag)
      exps_nam_calag <-
        grep("^smok_|^bmi_",
          exps_nam,
          value = TRUE,
          invert = TRUE)
      exps_nam_6lag <-
        grep("^smok_|^bmi_",
          exps_nam,
          value = TRUE,
          invert = FALSE) # Fixed lag for PLCO
      exps_tolag_calag <-
        grep("^smok_|^bmi_",
          exps_tolag,
          value = TRUE,
          invert = TRUE)
      exps_tolag_6lag <-
        grep("^smok_|^bmi_",
          exps_tolag,
          value = TRUE,
          invert = FALSE) # Fixed lag for PLCO

      for (i in seq_along(exps_nam_calag)) {
        set(dt, NULL, exps_nam_calag[i], dt[, shift_bypid(get(exps_tolag_calag[i]),
          design_$lags_mc$cancer_lag,
          pid)])
      }
      for (i in seq_along(exps_nam_6lag)) {
        set(dt, NULL, exps_nam_6lag[i], dt[, shift_bypid(get(exps_tolag_6lag[i]), design_$lags_mc$plco_lag, pid)])
      }
    } else {
      exps_tolag <-
        c(
          grep("^smok_.*_sc$", names(dt), value = TRUE),
          "ets_curr_xps",
          "fruit_curr_xps",
          "bmi_curr_xps"
        )
      # NOTE bmi effect of HC is intentionally ignored here. Otherwise would mean
      # that reducing bmi increases lung cancer which is unlikely (see A. Renhan
      # sys review for possible explanation)

      exps_nam <-
        gsub("_curr_xps$|_sc$",
          "_lagged",
          exps_tolag)

      exps_nam_calag <-
        grep("^smok_|^bmi_",
          exps_nam,
          value = TRUE,
          invert = TRUE)
      exps_nam_6lag <-
        grep("^smok_|^bmi_",
          exps_nam,
          value = TRUE,
          invert = FALSE) # Fixed lag for PLCO
      exps_tolag_calag <-
        grep("^smok_|^bmi_",
          exps_tolag,
          value = TRUE,
          invert = TRUE)
      exps_tolag_6lag <-
        grep("^smok_|^bmi_",
          exps_tolag,
          value = TRUE,
          invert = FALSE) # Fixed lag for PLCO

      for (i in seq_along(exps_nam_calag)) {
        set(dt, NULL, exps_nam_calag[i], dt[, shift_bypid(get(exps_tolag_calag[i]),
          design_$lags_mc$cancer_lag,
          pid)])
      }
      for (i in seq_along(exps_nam_6lag)) {
        set(dt, NULL, exps_nam_6lag[i], dt[, shift_bypid(get(exps_tolag_6lag[i]), design_$lags_mc$plco_lag, pid)])
      }
    }
    # RR for tobacco from Tammemägi MC, et al. Evaluation of the lung cancer
    # risks at which to screen ever- and never-smokers: screening rules applied
    # to the PLCO and NLST cohorts. PLoS Med 2014;11:e1001764. Table S1
    dt[year >= design_$sim_prm$init_year, c("lung_ca_forparf_rr",
      "lung_ca_nocopd_rr",
      "lung_ca_withcopd_rr") :=
        lung_ca_rr(
          smok_status_lagged,
          smok_cig_lagged,
          smok_dur_lagged,
          smok_quit_yrs_lagged,
          age - design_$lags_mc$plco_lag,
          # the age when the observation started. 6 years later is the event
          education,
          ethnicity,
          bmi_lagged,
          copd_forparf,
          history_of_ca,
          fam_lung_ca
        )]
    setnafill(
      dt,
      "c",
      1,
      cols = c(
        "lung_ca_forparf_rr",
        "lung_ca_nocopd_rr",
        "lung_ca_withcopd_rr"
      )
    )
    # dt[, summary(lung_ca_forparf_rr)]
    # dt[, summary(lung_ca_nocopd_rr)]
    # dt[, summary(lung_ca_withcopd_rr)]

    # RR for ETS from Kim CH, et al. Exposure to secondhand tobacco smoke and lung
    # cancer by histological type: A pooled analysis of the International Lung
    # Cancer Consortium (ILCCO). Int J Cancer 2014;135:1918-30.
    # doi:10.1002/ijc.28835
    # Very similar to Taylor R, et al. Meta-analysis of studies of passive smoking and lung cancer:
    # effects of study type and continent. Int J Epidemiol. 2007 Jan
    # 10;36(5):1048-59. Table 4 for Europe
    tt <-
      data.table(
        ets_lagged = 1L,
        smok_status_lagged = as.character(1:2),
        ets_rr = RR$ets_lung_ca$get_rr(mc, design_, drop = TRUE)
      )
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "ets_rr")
    # dt[, summary(ets_rr)]

    # RR for fruit from Vieira AR, et al. Fruits, vegetables and lung cancer risk: a
    # systematic review and meta-analysis. Ann Oncol 2016;27:81-96.
    # doi:10.1093/annonc/mdv381
    # Very similar to Wang Y, et al. Fruit and vegetable
    # consumption and risk of lung cancer: A dose-response meta-analysis of
    # prospective cohort studies. Lung Cancer 2015;88:124-30.
    # doi:10.1016/j.lungcan.2015.02.015
    dt[, fruit_rr := clamp(RR$fruit_lung_ca$get_rr(mc, design_, drop = TRUE) ^
        ((fruit_lagged - 80 * 5) / 80),
      1,
      20)] # no benefit after 5 portions
    setnafill(dt, "c", 1, cols = "fruit_rr")
    # dt[, summary(fruit_rr)]

    nam <- grep("_rr$", names(dt), value = TRUE)
    invisible(dt[, lapply(.SD, clamp, 0, 20, TRUE), .SDcols = nam]) # No rr > 20
    # dt[, lapply(.SD, max), .SDcols = nam]

    dt[, (exps_nam) := NULL]

    # Estimate PARF ------------------------------------------------------------
    #cat("Estimating lung cancer PAF...\n")
    if (!"p0_lung_ca" %in% names(dt)) {
      lung_caparf <-
        dt[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH) &
            lung_ca_prvl == 0 & year == design_$sim_prm$init_year,
          .(parf = 1 - 1 / (sum(lung_ca_forparf_rr * ets_rr * fruit_rr) / .N)),
          keyby = .(age, sex, qimd)]
      # lung_caparf[, parf := clamp(predict(loess(parf ~ age, span = 0.5))), by = .(sex, qimd)]
      # lung_caparf[, {
      #   plot(age, parf, main = paste0(.BY[[1]],"-", .BY[[2]]), ylim = c(0, 1))
      #   lines(age, parf2)
      # }
      # , keyby = .(sex, qimd)]

      absorb_dt(
        lung_caparf,
        get_disease_epi_mc(mc, "lung_ca", "i", "v", design_$sim_prm$stochastic)
      )
      lung_caparf[, p0_lung_ca := incidence * (1 - parf)]
      # lung_caparf[, summary(p0_lung_ca)]
      lung_caparf[is.na(p0_lung_ca), p0_lung_ca := incidence]
      lung_caparf[, c("incidence", "parf") := NULL]
      absorb_dt(dt, lung_caparf)
      setnafill(dt, "c", 0, cols = "p0_lung_ca")
      rm(lung_caparf)
    }

    if (!nzchar(scenario_nam)) {
      # Estimate lung cancer incidence prbl -------------------------------
      #cat("Estimating lung cancer incidence...\n\n")
      set(dt, NULL, "prb_lung_ca_incd_nocopd", 0)
      dt[, prb_lung_ca_incd_nocopd :=
          p0_lung_ca * lung_ca_nocopd_rr * ets_rr * fruit_rr] # Remember no copd as this will be a multiplier
      dt[, lung_ca_incd_copd_mltp := lung_ca_withcopd_rr / lung_ca_nocopd_rr]
      # TODO check that multiplier is OK to remain the same for all scenarios

      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis ----
      # probability to diagnosis every year
      # hist(rnbinom(1e4, 1, diagnosis_prb), 100) # the distribution of the number of years until diagnosis
      set(dt, NULL, "prb_lung_ca_dgn", diagnosis_prb)

      set(dt, NULL, "lung_ca_dgn", 0L)
      dt[lung_ca_prvl > 0, lung_ca_dgn := clamp(lung_ca_prvl - 5L, 0, 100)]

      # Estimate case fatality ----
      absorb_dt(dt,
        get_disease_epi_mc(mc, "lung_ca", "f", "v", design_$sim_prm$stochastic))
      setnames(dt, "fatality", "prb_lung_ca_mrtl")


    } else {
      set(dt, NULL, "lung_ca_prvl_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, lung_ca_prvl_sc := lung_ca_prvl]
      dt[year == design_$sim_prm$init_year_fromGUI &
          lung_ca_prvl > 1L, lung_ca_prvl_sc := lung_ca_prvl]

      # Estimate lung cancer incidence prbl
      #cat("Estimating lung cancer incidence without diabetes...\n\n")
      colnam <- "prb_lung_ca_incd_nocopd_sc"
      set(dt, NULL, colnam, 0)
      dt[, (colnam) :=
          p0_lung_ca * lung_ca_nocopd_rr * ets_rr * fruit_rr]
      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis
      dt[, prb_lung_ca_dgn_sc := prb_lung_ca_dgn]

      set(dt, NULL, "lung_ca_dgn_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, lung_ca_dgn_sc := lung_ca_dgn]
      dt[year == design_$sim_prm$init_year_fromGUI &
          lung_ca_dgn > 1L, lung_ca_dgn_sc := lung_ca_dgn]

      # Estimate case fatality
      if (!"prb_lung_ca_mrtl_sc" %in% names(dt))
        set(dt, NULL, "prb_lung_ca_mrtl_sc", dt$prb_lung_ca_mrtl)
    }
    dt[, copd_forparf := NULL]

    if (timing)
      print(proc.time() - ptm)
  }
