source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = TRUE)

scenario_fn <- function(sp) NULL

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:200, multicore = TRUE, "")

# IMPACTncd$export_summaries(multicore = TRUE)
# source("./aux/CPRD_sim_validation_plots.R")

scenario_fn <- function(sp) {
  sc_year <- 22L # The year the change starts
  change <- 0.2 # positive means health improvement. Do not set to 0


  sp$pop[year >= sc_year, bmi_curr_xps := bmi_curr_xps * (1 - change)]
  sp$pop[year >= sc_year, sbp_curr_xps := sbp_curr_xps * (1 - change)]
  sp$pop[year >= sc_year, tchol_curr_xps := tchol_curr_xps * (1 - change)]
  sp$pop[year >= sc_year, alcohol_curr_xps := as.integer(round(alcohol_curr_xps * (1 - change)))]
  sp$pop[year >= sc_year, fruit_curr_xps := as.integer(round(fruit_curr_xps * (1 + change)))]
  sp$pop[year >= sc_year, veg_curr_xps := as.integer(round(veg_curr_xps * (1 + change)))]

  # (20%) of the population increases physical activity by 1 day and met by 20%
  tt <- sp$pop[year >= sc_year, .(unique(pid))]
  tt <- tt[as.logical(rbinom(.N, 1, abs(change))), V1]
  sp$pop[pid %in% tt & year >= sc_year, `:=` (
    active_days_curr_xps = clamp(x = active_days_curr_xps + as.integer(sign(change)), a = 0L, b = 7L, inplace = FALSE),
    met_curr_xps = as.integer(round(met_curr_xps * (1 + change)))
    )]

  # smoking
  if (change > 0) { # Note positive means better health

    sp$pop[year >= sc_year & smok_status_curr_xps == "4",
       hc_eff := rbinom(.N, 1L, change)]

    sp$pop[, (c("smok_status_curr_xps", "smok_quit_yrs_curr_xps", "smok_dur_curr_xps", "smok_cig_curr_xps")) :=
         simsmok_policy_impact_decr(
           smok_status_curr_xps,
           smok_quit_yrs_curr_xps,
           smok_dur_curr_xps,
           smok_cig_curr_xps,
           pid_mrk,
           hc_eff
         )]
  } else {
    # calculate policy effect with those quit smoking recently be more
    # likely to relapse
    tt <- sp$pop[year >= sc_year, .("ex"   = sum(smok_status_curr_xps == "3"),
                                "curr" = sum(smok_status_curr_xps == "4")), keyby = year]
    tt[, impacted := round(curr * (1 - change))] # Note change here is -ve so 1 - change is an increase

    # Make change to add up every year (for Vincy's SCC abstract)
    # tt[, impacted := round(curr * scenario_parms$sc_str_smk_change *
    #     (year - min(year) + 1L))]

    sp$pop[tt, `:=`(impacted = i.impacted,
                ex = i.ex), on = "year"]
    sp$pop[year >= sc_year & smok_status_curr_xps == "3",
       rid := 1:.N, by = year]
    sp$pop[, hc_eff := 0L]
    tt <- sp$pop[year >= sc_year & smok_status_curr_xps == "3",
             .(rid = sample_int_expj(first(ex), first(impacted),
                                     (smok_quit_yrs_curr_xps + 1L) ^
                                       -1)),
             keyby = year]
    sp$pop[tt, hc_eff := 1L, on = .(year, rid)]
    sp$pop[, c("impacted", "ex", "rid") := NULL]

    sp$pop[, (c("smok_status_curr_xps", "smok_quit_yrs_curr_xps", "smok_dur_curr_xps")) :=
         simsmok_policy_impact_incr(
           smok_status_curr_xps,
           smok_quit_yrs_curr_xps,
           smok_dur_curr_xps,
           pid_mrk,
           hc_eff
         )]
  }
  sp$pop[, hc_eff := NULL]

  sp$pop[year >= sc_year, smok_cig_curr_xps := as.integer(round(smok_cig_curr_xps * (1 - change)))]

  sp$pop[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]



  # ets need to be after smoking
  if (change > 0) {
    sp$pop[ets_curr_xps == 1L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change)),
           ets_curr_xps := 0L]
  } else {
    sp$pop[ets_curr_xps == 0L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change)),
           ets_curr_xps := 1L]
  }

}

IMPACTncd$
  run(1:200, multicore = TRUE, "sc1")$
  export_summaries(multicore = TRUE)

source("./aux/CPRD_sim_validation_plots.R")
source("./aux/process_out_for_HF.R")
# IMPACTncd$export_summaries(multicore = TRUE)


