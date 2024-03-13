source("global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

# IMPACTncd$update_output_path("/mnt/storage_fast/output/test")
# IMPACTncd$update_synthpop_path("/mnt/storage_fast/synthpop/test")
# g <- IMPACTncd$get_causal_structure(print_plot = TRUE)
# g <- IMPACTncd$get_causal_structure(processed = FALSE, print_plot = TRUE, focus = "chd")

#plot(igraph::make_ego_graph(g, order = 1, c("pain"), "in")[[1]])

scenario_fn_primary_prevention   <- function(sp) NULL
scenario_fn_secondary_prevention <- function(sp) NULL
# multicore - F
# baseline
IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:2, multicore = F, "sc0")

# IMPACTncd$export_summaries(multicore = TRUE)
# source("./auxil/CPRD_sim_validation_plots.R")

#' Apply Primary Prevention Strategies to Simulate Health Improvements
#'
#' This function models the impact of primary prevention strategies on various
#' risk factors in a population over time. It adjusts BMI, systolic blood
#' pressure, cholesterol, alcohol consumption, fruit and vegetable intake,
#' physical activity, smoking prevalence, smoking consumption, and passive
#' smoking prevalence based on specified changes.
#'
#' @details
#' The major risk factors included in our model in three cases have “healthy”
#' ranges according to government guidelines, body mass index (BMI), systolic
#' blood pressure and cholesterol.
#'
#' @param sp A data frame containing population information.
#' @param change A numeric value indicating the magnitude of health improvement.
#' A positive value means health improvement. Do not set to 0.
#' @param sc_year An integer specifying the year when the change starts.
#'
#' @return The input data frame `sp` with modified risk factors based on the
#' specified primary prevention strategies.
#'
#' @export
scenario_fn_primary_prevention <- function(sp) {

  # The major risk factors included in our model in three cases have “healthy”
  # ranges according to government guidelines, body mass index (BMI), systolic
  # blood pressure and cholesterol: for example the NHS advises an adults BMI in
  # the range 18.5-24.9. In this setting, the 20% improvement of risk factors is
  # relative to the distance from the middle of the range (21.7): someone with a
  # BMI of 31.7 is 10kg/m2 above the middle of the range, so a 20% improvement
  # is a 2kg/m2 reduction to 29.7. The NHS advised ranges for blood pressure and
  # cholesterol are 90/60mmHg - 120/80mmHg, for total cholesterol the
  # recommendation is below 5mmol/L, but above 1mmol/L for “good cholesterol”.
  # We model improvements in the other risk factors as follows: fruit and
  # vegetable intake is 20% higher than the base case for each person; alcohol
  # intake is 20% lower; 20% of the population become more active by one active
  # day per week; smoking prevalence drops by 20% (relative); smoking
  # consumption for those that still smoke drops by 20%; passive smoking
  # prevalence drops by 20% (relative).
  sc_year <- 23L # The year the change starts
  change <- 0.2 # positive means health improvement. Do not set to 0

  bmi_target <- mean(c(18.5, 24.9))
  sp$pop[year >= sc_year & bmi_curr_xps > bmi_target, bmi_curr_xps := bmi_curr_xps + (bmi_target - bmi_curr_xps) * change]

  sbp_target <- mean(c(90, 120))
  sp$pop[year >= sc_year & sbp_curr_xps > sbp_target, sbp_curr_xps := sbp_curr_xps + (sbp_target - sbp_curr_xps) * change]

  tchol_target <- mean(c(3, 5)) # Lower bound arbitrary but high enough to ensure HDL > 1mmol/L possible
  sp$pop[year >= sc_year & tchol_curr_xps > tchol_target, tchol_curr_xps := tchol_curr_xps + (tchol_target - tchol_curr_xps) * change]



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
    sp$pop[year >= sc_year & ets_curr_xps == 1L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change)),
           ets_curr_xps := 0L]
  } else {
    sp$pop[year >= sc_year & ets_curr_xps == 0L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change)),
           ets_curr_xps := 1L]
  }

}


IMPACTncd$
  run(1:2, multicore = F, "sc1")$
  export_summaries(multicore = F)

# IMPACTncd$export_summaries(multicore = TRUE)
source("./auxil/process_out_for_HF.R")
# if (IMPACTncd$design$sim_prm$validation) source("./auxil/CPRD_sim_validation_plots_CK.R")

