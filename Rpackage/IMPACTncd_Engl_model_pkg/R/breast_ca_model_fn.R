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
breast_ca_model <-
  function(scenario_nam,
    mc,
    dt,
    design_,
    diagnosis_prb = 0.7,
    timing = TRUE) {
    message("Loading breast cancer (C50) model...")
    if (timing)
      ptm <- proc.time()

    if (!nzchar(scenario_nam)) {
      # first run for scenario ""
      # Lagged exposures
      exps_tolag <-
        grep(
          "^smok_sta|^ets_|^active_|^bmi_|^alcohol_|^t2dm_prvl_curr_xps",
          names(dt),
          value = TRUE
        )
      exps_nam <-  gsub("_curr_xps$", "_lagged", exps_tolag)
      for (i in seq_along(exps_tolag)) {
        set(dt, NULL, exps_nam[i],
          dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$cancer_lag, pid)])
      }
    } else {
      exps_tolag <-
        c(paste0(
          c(
            "active_days_",
            "alcohol_",
            "bmi_",
            "smok_status_",
            "t2dm_prvl_"
          ),
          "sc"
        ),
          "ets_curr_xps")

      exps_nam <-
        gsub("_curr_xps$|_sc$",
          "_lagged",
          exps_tolag)
      for (i in seq_along(exps_tolag)) {
        set(dt, NULL, exps_nam[i], dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$cvd_lag, pid)])
      }

    }

    # RR for tobacco from Macacu A, Autier P, Boniol M, Boyle P. Active
    # and passive smoking and risk of breast cancer: a meta-analysis.
    # Breast Cancer Res Treat 2015;154:213–24.
    # Table 1
    tt <- RR$tobacco_breast_ca$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "smok_status", "smok_status_lagged")
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "tobacco_rr")
    dt[sex == "men", tobacco_rr := 1]
    # dt[, summary(tobacco_rr)]

    # RR for ETS GBD
    tt <-
      data.table(
        ets_lagged = 1L,
        sex = "women",
        smok_status_lagged = as.character(1:2),
        ets_rr = RR$ets_breast_ca$get_rr(mc, design_, drop = TRUE)
      )
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "ets_rr")

    # RR for alcohol from GBD 2016
    tt <-
      RR$alcohol_breast_ca$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "alcohol", "alcohol_lagged")
    dt[alcohol_lagged > max(tt$alcohol_lagged),
      alcohol_lagged := max(tt$alcohol_lagged)]
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "alcohol_rr")
    dt[sex == "men", alcohol_rr := 1]
    # dt[, summary(alcohol_rr)]

    # RR for PA from 1. Bull FC, Armstrong TP, Dixon T, Ham S, Neiman A, Pratt M.
    # Comparative quantification of health risks. Chapter 10: physical inactivity.
    # Geneva: World Health Organisation; 2004.
    # Available from: http://www.who.int/publications/cra/en/
    # table 10.21, estimates with adjustment for measurement error
    tt <-
      RR$pa_breast_ca$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "active_days", "active_days_lagged")
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "pa_rr")
    dt[sex == "men", pa_rr := 1]
    # dt[, summary(pa_rr)]

    # RR for BMI from GBD
    tt <- RR$bmi_breast_ca$get_rr(mc, design_, drop = TRUE)
    absorb_dt(dt, tt)
    dt[, bmi_rr := clamp(bmi_rr ^ ((bmi_lagged - 20) / 5), 1, 20)]
    setnafill(dt, "c", 1, cols = "bmi_rr")
    dt[sex == "men", bmi_rr := 1]
    # dt[, summary(bmi_rr)]

    # RR for t2dm from GBD
    set(dt, NULL, "t2dm_rr", 1)
    tt <- RR$t2dm_breast_ca$get_rr(mc, design_, drop = TRUE)
    dt[t2dm_prvl_lagged > 0, t2dm_rr := tt]
    dt[sex == "men", t2dm_rr := 1]

    nam <- grep("_rr$", names(dt), value = TRUE)
    invisible(dt[, lapply(.SD, clamp, 0, 20, TRUE), .SDcols = nam]) # No rr > 20
    # dt[, lapply(.SD, max), .SDcols = nam]

    dt[, (exps_nam) := NULL]

    # Estimate PARF ------------------------------------------------------------
    #cat("Estimating breast_ca PAF...\n")
    if (!"p0_breast_ca" %in% names(dt)) {
      breast_caparf <-
        dt[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH) &
            breast_ca_prvl == 0 & year == design_$sim_prm$init_year,
          .(parf = 1 - 1 / (
            sum(tobacco_rr * ets_rr * alcohol_rr * pa_rr * bmi_rr *
                t2dm_rr) / .N
          )),
          keyby = .(age, sex, qimd)]
      # breast_caparf[, parf := clamp(predict(loess(parf ~ age, span = 0.5))), by = .(sex, qimd)]
      # breast_caparf[, {
      #   plot(age, parf, main = paste0(.BY[[1]],"-", .BY[[2]]), ylim = c(0, 1))
      #   lines(age, parf)
      # },
      # keyby = .(sex, qimd)]

      absorb_dt(
        breast_caparf,
        get_disease_epi_mc(mc, "breast_ca", "i", "v", design_$sim_prm$stochastic)
      )
      breast_caparf[, p0_breast_ca := incidence * (1 - parf)]
      # breast_caparf[, summary(p0_breast_ca)]
      breast_caparf[is.na(p0_breast_ca), p0_breast_ca := incidence]
      breast_caparf[, c("incidence", "parf") := NULL]
      absorb_dt(dt, breast_caparf)
      setnafill(dt, "c", 0, cols = "p0_breast_ca")
      rm(breast_caparf)
    }

    if (!nzchar(scenario_nam)) {
      # Estimate breast_ca incidence prbl -------------------------------
      #cat("Estimating breast_ca incidence...\n\n")
      set(dt, NULL, "prb_breast_ca_incd_not2dm", 0)
      dt[, prb_breast_ca_incd_not2dm :=
          p0_breast_ca * tobacco_rr * ets_rr * alcohol_rr * pa_rr * bmi_rr]
      # Remember no t2dm as this will be a multiplier
      dt[, breast_ca_incd_t2dm_mltp :=
          RR$t2dm_breast_ca$get_rr(mc, design_, drop = TRUE)]
      dt[sex == "men", breast_ca_incd_t2dm_mltp := 1]

      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis ----
      # probability to diagnosis every year
      # hist(rnbinom(1e4, 1, diagnosis_prb), 100) # the distribution of the number of years until diagnosis
      set(dt, NULL, "prb_breast_ca_dgn", diagnosis_prb)

      set(dt, NULL, "breast_ca_dgn", 0L)
      dt[breast_ca_prvl > 0, breast_ca_dgn := clamp(breast_ca_prvl - 5L, 0, 100)]

      # Estimate case fatality ----
      absorb_dt(dt,
        get_disease_epi_mc(mc, "breast_ca", "f", "v", design_$sim_prm$stochastic))
      setnames(dt, "fatality", "prb_breast_ca_mrtl")

    } else {
      set(dt, NULL, "breast_ca_prvl_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, breast_ca_prvl_sc := breast_ca_prvl]
      dt[year == design_$sim_prm$init_year_fromGUI &
          breast_ca_prvl > 1L, breast_ca_prvl_sc := breast_ca_prvl]

      # Estimate breast cancer incidence prbl
      #cat("Estimating breast cancer incidence without diabetes...\n\n")
      colnam <- "prb_breast_ca_incd_not2dm_sc"
      set(dt, NULL, colnam, 0)
      dt[, (colnam) :=
          p0_breast_ca * tobacco_rr * ets_rr * alcohol_rr * pa_rr * bmi_rr] # Remeber no t2dm as this will be a multiplier
      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis
      dt[, prb_breast_ca_dgn_sc := prb_breast_ca_dgn]

      set(dt, NULL, "breast_ca_dgn_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, breast_ca_dgn_sc := breast_ca_dgn]
      dt[year == design_$sim_prm$init_year_fromGUI &
          breast_ca_dgn > 1L, breast_ca_dgn_sc := breast_ca_dgn]

      # Estimate case fatality
      if (!"prb_breast_ca_mrtl_sc" %in% names(dt))
        set(dt, NULL, "prb_breast_ca_mrtl_sc", dt$prb_breast_ca_mrtl)

    }
    if (timing)
      print(proc.time() - ptm)
  }
# dt[, lapply(.SD, anyNA)]
