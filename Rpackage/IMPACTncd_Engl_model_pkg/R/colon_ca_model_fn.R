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
colon_ca_model <-
  function(scenario_nam,
    mc,
    dt,
    design_,
    diagnosis_prb = 0.7,
    timing = TRUE) {
    message("Loading colon cancer (C18) model...")
    if (timing)
      ptm <- proc.time()

    if (!nzchar(scenario_nam)) {
      # first run for scenario ""

      # Calculate pack years
      dt[, smok_packyrs_curr_xps := smok_cig_curr_xps * smok_dur_curr_xps / 20]

      # Lagged exposures
      exps_tolag <-
        grep("^smok_pack|^active_|^bmi_|^alcohol_|^t2dm_prvl_curr_xps",
          names(dt),
          value = TRUE)
      exps_nam <-  gsub("_curr_xps$", "_lagged", exps_tolag)
      for (i in seq_along(exps_tolag)) {
        set(dt, NULL, exps_nam[i], dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$cancer_lag, pid)])
      }

    } else {
      # Calculate pack years
      dt[, smok_packyrs_sc := smok_cig_sc * smok_dur_sc / 20]

      exps_tolag <-
        c(paste0(
          c(
            "active_days_",
            "alcohol_",
            "bmi_",
            "smok_packyrs_",
            "t2dm_prvl_"
          ),
          "sc"
        ))

      exps_nam <-
        gsub("_curr_xps$|_sc$", "_lagged", exps_tolag)
      for (i in seq_along(exps_tolag)) {
        set(dt, NULL, exps_nam[i], dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$cvd_lag, pid)])
      }
    }

    # RR for tobacco from 1. Liang PS, Chen T-Y, Giovannucci E. Cigarette
    # smoking and colorectal cancer incidence and mortality: Systematic review
    # and meta-analysis. International Journal of Cancer [Internet]
    # 2009;124:2406–15. Available from:
    # https://onlinelibrary.wiley.com/doi/abs/10.1002/ijc.24191
    # Few studies have separate colon cancer RR thus results are non significant
    # I will use packyears for this and only consider incidence and not mortality
    # Table III and text for CI
    tt <- RR$packyears_colon_ca$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "packyears", "smok_packyrs_lagged")
    dt[smok_packyrs_lagged > max(tt$smok_packyrs_lagged),
      smok_packyrs_lagged := max(tt$smok_packyrs_lagged)]
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "packyears_rr")
    # dt[, summary(packyears_rr)]

    # RR for alcohol from GBD 2016
    tt <- RR$alcohol_colon_ca$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "alcohol", "alcohol_lagged")
    dt[alcohol_lagged > max(tt$alcohol_lagged),
      alcohol_lagged := max(tt$alcohol_lagged)]
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "alcohol_rr")
    # dt[, summary(alcohol_rr)]

    # RR for PA from 1. Bull FC, Armstrong TP, Dixon T, Ham S, Neiman A, Pratt M.
    # Comparative quantification of health risks. Chapter 10: physical inactivity.
    # Geneva: World Health Organisation; 2004.
    # Available from: http://www.who.int/publications/cra/en/
    # table 10.22, p111, estimates with adjustment for measurement error
    tt <- RR$pa_colon_ca$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "active_days", "active_days_lagged")
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "pa_rr")
    # dt[, summary(pa_rr)]

    # RR for BMI from GBD
    tt <- RR$bmi_colon_ca$get_rr(mc, design_, drop = TRUE)
    absorb_dt(dt, tt)
    dt[, bmi_rr := clamp(bmi_rr ^ ((bmi_lagged - 20) / 5), 1, 20)]
    setnafill(dt, "c", 1, cols = "bmi_rr")
    # dt[, summary(bmi_rr)]

    # RR for t2dm from GBD
    set(dt, NULL, "t2dm_rr", 1)
    tt <- RR$t2dm_colon_ca$get_rr(mc, design_, drop = TRUE)
    dt[t2dm_prvl_lagged > 0, t2dm_rr := tt]

    nam <- grep("_rr$", names(dt), value = TRUE)
    invisible(dt[, lapply(.SD, clamp, 0, 20, TRUE), .SDcols = nam]) # No rr > 20
    # dt[, lapply(.SD, max), .SDcols = nam]

    dt[, (exps_nam) := NULL]

    # Estimate PARF ------------------------------------------------------------
    #cat("Estimating colon_ca PAF...\n")
    if (!"p0_colon_ca" %in% names(dt)) {
      colon_caparf <-
        dt[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH) &
            colon_ca_prvl == 0 & year == design_$sim_prm$init_year,
          .(parf = 1 - 1 / (sum(
            packyears_rr * alcohol_rr * pa_rr * bmi_rr *
              t2dm_rr
          ) / .N)),
          keyby = .(age, sex, qimd)]
      # colon_caparf[, parf := clamp(predict(loess(parf ~ age, span = 0.5))), by = .(sex, qimd)]
      # colon_caparf[, {
      #   plot(age, parf, main = paste0(.BY[[1]],"-", .BY[[2]]), ylim = c(0, 1))
      #   lines(age, parf2)
      # }
      # , keyby = .(sex, qimd)]

      absorb_dt(
        colon_caparf,
        get_disease_epi_mc(mc, "colon_ca", "i", "v", design_$sim_prm$stochastic)
      )
      colon_caparf[, p0_colon_ca := incidence * (1 - parf)]
      # colon_caparf[, summary(p0_colon_ca)]
      colon_caparf[is.na(p0_colon_ca), p0_colon_ca := incidence]
      colon_caparf[, c("incidence", "parf") := NULL]
      absorb_dt(dt, colon_caparf)
      setnafill(dt, "c", 0, cols = "p0_colon_ca")
      rm(colon_caparf)
    }

    if (!nzchar(scenario_nam)) {
      # first run for scenario ""

      # Estimate colon_ca incidence prbl -------------------------------
      #cat("Estimating colon_ca incidence...\n\n")
      set(dt, NULL, "prb_colon_ca_incd_not2dm", 0)
      dt[, prb_colon_ca_incd_not2dm :=
          p0_colon_ca * packyears_rr * alcohol_rr * pa_rr * bmi_rr]
      # Remember no t2dm as this will be a multiplier
      dt[, colon_ca_incd_t2dm_mltp := RR$t2dm_colon_ca$get_rr(mc, design_, drop = TRUE)]


      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis ----
      # probability to diagnosis every year
      # hist(rnbinom(1e4, 1, diagnosis_prb), 100) # the distribution of the number of years until diagnosis
      set(dt, NULL, "prb_colon_ca_dgn", diagnosis_prb)

      set(dt, NULL, "colon_ca_dgn", 0L)
      dt[colon_ca_prvl > 0, colon_ca_dgn := clamp(colon_ca_prvl - 5L, 0, 100)]

      # Estimate case fatality ----
      absorb_dt(dt,
        get_disease_epi_mc(mc, "colon_ca", "f", "v", design_$sim_prm$stochastic))
      setnames(dt, "fatality", "prb_colon_ca_mrtl")


      dt[, smok_packyrs_curr_xps := NULL]
    } else {
      set(dt, NULL, "colon_ca_prvl_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, colon_ca_prvl_sc := colon_ca_prvl]
      dt[year == design_$sim_prm$init_year_fromGUI &
          colon_ca_prvl > 1L,
        colon_ca_prvl_sc := colon_ca_prvl]

      # Estimate colon cancer incidence prbl
      #cat("Estimating colon cancer incidence without diabetes...\n\n")
      colnam <- "prb_colon_ca_incd_not2dm_sc"
      set(dt, NULL, colnam, 0)
      dt[, (colnam) :=
          p0_colon_ca * packyears_rr * alcohol_rr * pa_rr * bmi_rr]
      # Remember no t2dm as this will be a multiplier
      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis
      dt[, prb_colon_ca_dgn_sc := prb_colon_ca_dgn]

      set(dt, NULL, "colon_ca_dgn_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, colon_ca_dgn_sc := colon_ca_dgn]
      dt[year == design_$sim_prm$init_year_fromGUI &
          colon_ca_dgn >= 1L, colon_ca_dgn_sc := colon_ca_dgn]

      # Estimate case fatality
      if (!"prb_colon_ca_mrtl_sc" %in% names(dt))
        set(dt, NULL, "prb_colon_ca_mrtl_sc", dt$prb_colon_ca_mrtl)

      dt[, smok_packyrs_sc := NULL]
    }
    if (timing)
      print(proc.time() - ptm)
  }
# dt[, lapply(.SD, anyNA)]
