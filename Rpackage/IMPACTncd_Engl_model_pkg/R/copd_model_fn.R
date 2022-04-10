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
copd_model <-
  function(
    scenario_nam,
    mc,
    dt,
    design_,
    diagnosis_prb = 0.7,
    timing = TRUE) {
    message("Loading COPD (J40-J47) model...")
    if (timing)
      ptm <- proc.time()

    if (!nzchar(scenario_nam)) { # first run for scenario ""
      # Calculate pack years
      dt[, smok_packyrs_curr_xps := smok_cig_curr_xps * smok_dur_curr_xps / 20]

    # Lagged exposures
    exps_tolag <- grep("^smok_pack|^smok_st|^ets",
      names(dt), value = TRUE)
    exps_nam <-  gsub("_curr_xps$", "_lagged", exps_tolag)
    for (i in seq_along(exps_tolag)) {
      set(dt, NULL, exps_nam[i], dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$copd_lag, pid)])
    }
    } else {
      # Calculate pack years
      dt[, smok_packyrs_sc := smok_cig_sc * smok_dur_sc / 20]

      exps_tolag <-
        c(paste0(
          c(
            "smok_packyrs_", "smok_status_"
          ),
          "sc"
        ),
          "ets_curr_xps")

      exps_nam <-
        gsub("_curr_xps$|_sc$", "_lagged", exps_tolag)
      for (i in seq_along(exps_tolag)) {
        set(dt, NULL, exps_nam[i], dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$copd_lag, pid)])
      }
}
    # RR for tobacco from Forey BA, Thornton AJ, Lee PN. Systematic review with meta
    # -analysis of the epidemiological evidence relating smoking to COPD, chronic
    # bronchitis and emphysema. BMC Pulmonary Medicine. 2011 Jun 14;11(1):36.
    # A huge source of info. I will use RR for packyears. Note that duration is not
    # significant but intensity and packyears are with packyears indirectly capture
    # age effect. Most studies for packyears were about incidence rather than
    # mortality copd. The problem is that there is no differentiation between current
    # and ex-smokers. This dilutes the effect but currently I see no workaround.

    tt <- RR$packyears_copd$get_rr(mc, design_, drop = TRUE)
    setnames(tt, "packyears", "smok_packyrs_lagged")
    dt[smok_packyrs_lagged > max(tt$smok_packyrs_lagged),
      smok_packyrs_lagged := max(tt$smok_packyrs_lagged)]
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "packyears_rr")
    # dt[, summary(packyears_rr)]

    # RR for ETS from Fischer F, Kraemer A. Meta-analysis of the association between
    # second-hand smoke exposure and ischaemic heart diseases, COPD and stroke.
    # BMC Public Health. 2015 Dec;15(1):1202.
    tt <-
      data.table(
        ets_lagged = 1L,
        smok_status_lagged = as.character(1:2),
        ets_rr = RR$ets_copd$get_rr(mc, design_, drop = TRUE)
      )
    absorb_dt(dt, tt)
    setnafill(dt, "c", 1, cols = "ets_rr")

    nam <- grep("_rr$", names(dt), value = TRUE)
    invisible(dt[, lapply(.SD, clamp, 0, 20, TRUE), .SDcols = nam]) # No rr > 20
    # dt[, lapply(.SD, max), .SDcols = nam]

    dt[, (exps_nam) := NULL]

    # Estimate PARF ------------------------------------------------------------
    #cat("Estimating COPD PAF...\n")
    if (!"p0_copd" %in% names(dt)) {
      copdparf <-
        dt[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH) &
            copd_prvl == 0 & year == design_$sim_prm$init_year,
          .(parf = 1 - 1 / (sum(packyears_rr * ets_rr) / .N)),
          keyby = .(age, sex, qimd)]
      # copdparf[, parf := clamp(predict(loess(parf ~ age, span = 0.5))), by = .(sex, qimd)]
      # copdparf[, {
      #   plot(age, parf, main = paste0(.BY[[1]],"-", .BY[[2]]), ylim = c(0, 1))
      #   lines(age, parf2)
      # }
      # , keyby = .(sex, qimd)]

      absorb_dt(copdparf,
        get_disease_epi_mc(mc, "copd", "i", "v", design_$sim_prm$stochastic))
      copdparf[, p0_copd := incidence * (1 - parf)]
      # copdparf[, summary(p0_copd)]
      copdparf[is.na(p0_copd), p0_copd := incidence]
      copdparf[, c("incidence", "parf") := NULL]
      absorb_dt(dt, copdparf)
      setnafill(dt, "c", 0, cols = "p0_copd")
      rm(copdparf)
    }

    if (!nzchar(scenario_nam)) { # first run for scenario ""
    # Estimate COPD incidence prbl -------------------------------
    #cat("Estimating COPD incidence...\n\n")
    set(dt, NULL, "prb_copd_incd", 0)
    dt[, prb_copd_incd :=
        p0_copd * packyears_rr * ets_rr]
    dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

    # Assume a probability of diagnosis ----
    # probability to diagnosis every year
    # hist(rnbinom(1e4, 1, diagnosis_prb), 100) # the distribution of the number of years until diagnosis
    set(dt, NULL, "prb_copd_dgn", diagnosis_prb)

    set(dt, NULL, "copd_dgn", 0L)
    dt[copd_prvl > 0, copd_dgn := clamp(copd_prvl - 5L, 0, 100)]

    # Estimate case fatality ----
    absorb_dt(dt, get_disease_epi_mc(mc, "copd", "f", "v", design_$sim_prm$stochastic))
    setnames(dt, "fatality", "prb_copd_mrtl")

    dt[, smok_packyrs_curr_xps := NULL]
    } else {

      set(dt, NULL, "copd_prvl_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, copd_prvl_sc := copd_prvl]
      dt[year == design_$sim_prm$init_year_fromGUI & copd_prvl > 1L, copd_prvl_sc := copd_prvl]

      # Estimate COPD incidence prbl
      #cat("Estimating COPD incidence...\n\n")
      colnam <- "prb_copd_incd_sc"
      set(dt, NULL, colnam, 0)
      dt[, (colnam) :=
          p0_copd * packyears_rr * ets_rr]
      dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

      # Assume a probability of diagnosis
      dt[, prb_copd_dgn_sc := prb_copd_dgn]

      set(dt, NULL, "copd_dgn_sc", 0L)
      dt[year < design_$sim_prm$init_year_fromGUI, copd_dgn_sc := copd_dgn]
      dt[year == design_$sim_prm$init_year_fromGUI & copd_dgn > 1L, copd_dgn_sc := copd_dgn]

      # Estimate case fatality
      if (!"prb_copd_mrtl_sc" %in% names(dt))
        set(dt, NULL, "prb_copd_mrtl_sc", dt$prb_copd_mrtl)

      dt[, smok_packyrs_sc := NULL]
}
    if (timing) print(proc.time() - ptm)
  }

