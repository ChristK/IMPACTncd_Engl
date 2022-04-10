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
nonmodelled_model <- function(
  scenario_nam,
  mc, # Not currently used. For consistency
  dt,
  design_,
  timing = TRUE) {

  message("Estimating deaths from non-modelled causes...")
  if (timing) ptm <- proc.time()

  if (!nzchar(scenario_nam)) { # first run for scenario ""
    # Lagged exposures
    exps_tolag <-
      grep(
        "^act|^smok_st|^alco|^sbp|^t2dm_prvl_c",
        names(dt),
        value = TRUE
      )
    exps_nam <-  gsub("_curr_xps$", "_lagged", exps_tolag)
    for (i in seq_along(exps_tolag)) {
      set(dt, NULL, exps_nam[i],
        dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$nonmodelled, pid)])
    }

  } else {

    exps_tolag <-
      c(paste0(
        c(
          "active_days_", "alcohol_", "sbp_", "smok_status_",
          "t2dm_prvl_"
        ),
        "sc"
      ))

    exps_nam <-
      gsub("_curr_xps$|_sc$",
        "_lagged",
        exps_tolag)
    for (i in seq_along(exps_tolag)) {
      set(dt, NULL, exps_nam[i],
        dt[, shift_bypid(get(exps_tolag[i]), design_$lags_mc$nonmodelled_lag, pid)])
    }
  }

  # RR from  Stringhini S, et al. Socioeconomic status and the 25 × 25 risk
  # factors as determinants of premature mortality: a multicohort study and
  # meta-analysis of 1·7 million men and women. The Lancet 2017;389:1229–37.
  # figure 4
  tt <- RR$tobacco_nonmodelled$get_rr(mc, design_, drop = TRUE)
  set(dt, NULL, "tobacco_rr", 1)
  dt[smok_status_lagged == "4", tobacco_rr := tt]
  # dt[, summary(tobacco_rr)]

  tt <- RR$t2dm_nonmodelled$get_rr(mc, design_, drop = TRUE)
  set(dt, NULL, "t2dm_rr", 1)
  dt[t2dm_prvl_lagged > 0L, t2dm_rr := tt]
  dt[, nonmodelled_mrtl_t2dm_mltp := tt]
  # dt[, summary(t2dm_rr)]

  tt <- RR$pa_nonmodelled$get_rr(mc, design_, drop = TRUE)
  set(dt, NULL, "pa_rr", 1)
  dt[active_days_lagged < 2L, pa_rr := tt]
  # dt[, summary(pa_rr)]

  # One unit of alcohol (UK) is defined as 8 grams of pure alcohol.
  # The study uses cut-off 21 units for men and 14 for women per week.
  # In grams 24g for men and 16g for women
  tt <- RR$alcohol_nonmodelled$get_rr(mc, design_, drop = TRUE)
  set(dt, NULL, "alcohol_rr", 1)
  dt[(sex == "men" & alcohol_lagged > 24)|
      (sex == "women" & alcohol_lagged > 16), alcohol_rr := tt]
  # dt[, summary(alcohol_rr)]

  tt <- RR$sbp_nonmodelled$get_rr(mc, design_, drop = TRUE)
  set(dt, NULL, "sbp_rr", 1)
  dt[sbp_lagged > 140, sbp_rr := tt]
  # dt[, summary(sbp_rr)]

  dt[, (exps_nam) := NULL]

   # Estimate PARF ------------------------------------------------------------
  #cat("Estimating nonmodelled PAF...\n")
  if (!"p0_nonmodelled" %in% names(dt)) {
    nonmodelledparf <-
      dt[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH) & year == design_$sim_prm$init_year,
         .(parf = 1 - 1 / (
           sum(
             tobacco_rr * sbp_rr * pa_rr * alcohol_rr * t2dm_rr
           ) / .N
         )),
         # I should include year  to avoid double counting of exposure trends
         # that have been implicitly included in mortality forecasts. However I
         # don't have t2dm prevalence post 2013 so this is currently impossible
         keyby = .(age, sex, qimd)]
    # nonmodelledparf[, parf := clamp(predict(loess(parf ~ age, span = 0.5))), by = .(sex, qimd)]
    # nonmodelledparf[, {
    #   plot(age, parf, main = paste0(.BY[[1]],"-", .BY[[2]]), ylim = c(0, 1))
    #   lines(age, parf2)
    # }
    # , keyby = .(sex, qimd)]

    tt <- get_lifetable_all(mc, "nonmodelled", design_$sim_prm, "qx")[year == design_$sim_prm$init_year][, year := NULL]
    absorb_dt(nonmodelledparf, tt)
    nonmodelledparf[, p0_nonmodelled := qx_mc * (1 - parf)]
    # nonmodelledparf[, summary(p0_nonmodelled)]
    nonmodelledparf[is.na(p0_nonmodelled), p0_nonmodelled := qx_mc]
    nonmodelledparf[, c("qx_mc", "parf") := NULL]
    absorb_dt(dt, nonmodelledparf)
    setnafill(dt, "c", 0, cols = "p0_nonmodelled")
    rm(nonmodelledparf)
  }
  if (!nzchar(scenario_nam)) { # first run for scenario ""
    # Estimate nonmodelled mortality prbl -------------------------------
    # cat("Estimating nonmodelled incidence without diabetes...\n\n")
    set(dt, NULL, "prb_nonmodelled_mrtl_not2dm", 0)
    dt[, prb_nonmodelled_mrtl_not2dm :=
         p0_nonmodelled * tobacco_rr * sbp_rr * pa_rr * alcohol_rr] # Remember no t2dm as this will be a multiplier

    # Calibration
    # dt[, prb_nonmodelled_mrtl_not2dm := prb_nonmodelled_mrtl_not2dm * 1.15]
    # tt <- get_lifetable_all(mc, "nonmodelled", design_, "qx")
    # dt[tt, on = .NATURAL, prb_nonmodelled_mrtl_not2dm := qx_mc]
    # dt[, nonmodelled_mrtl_t2dm_mltp := 1]

  } else {
    # Estimate nonmodelled mortality prbl
    # cat("Estimating nonmodelled incidence without diabetes...\n\n")
    colnam <- "prb_nonmodelled_mrtl_not2dm_sc"
    set(dt, NULL, colnam, 0)

    if ("p0_nonmodelled_sc" %in% names(dt)) {
      dt[, (colnam) :=
          p0_nonmodelled_sc * tobacco_rr * sbp_rr * pa_rr * alcohol_rr] # Remember no t2dm as this will be a multiplier
    } else {
      dt[, (colnam) :=
          p0_nonmodelled * tobacco_rr * sbp_rr * pa_rr * alcohol_rr] # Remember no t2dm as this will be a multiplier
    }
  }

  dt[, (grep("_rr$", names(dt), value = TRUE)) := NULL]

  if (timing) print(proc.time() - ptm)
}



