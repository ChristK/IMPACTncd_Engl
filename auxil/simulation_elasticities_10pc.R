source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = T)


scenario_fn_primary_prevention   <- function(sp) NULL
scenario_fn_secondary_prevention <- function(sp) NULL

#change_1pc <- 0.01
change_10pc <- 0.1

sc_year <- 23L # The year the change starts

n_runs <- 100 #Number of runs; check don't crash and then run 100


IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:n_runs, multicore = TRUE, "")


# Checking minimal theoretical exposures
# IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl() # 21.9502
# IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl() #  112.4769
# IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl() # 3.874455
# IMPACTncd$RR$`fruit~chd`$get_ideal_xps_lvl() # 3.959
# IMPACTncd$RR$`veg~chd`$get_ideal_xps_lvl() # 3.975
#
#


# 10pc ----
#
# BMI ----
scenario_fn_primary_prevention <- function(sp) {
  bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & bmi_curr_xps > bmi_target, bmi_curr_xps := bmi_curr_xps + (bmi_target - bmi_curr_xps) * change_10pc]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "bmi_10pc")




# SBP ----

scenario_fn_primary_prevention <- function(sp) {
  sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & sbp_curr_xps > sbp_target, sbp_curr_xps := sbp_curr_xps + (sbp_target - sbp_curr_xps) * change_10pc]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "sbp_10pc")



# Total chol ----
scenario_fn_primary_prevention <- function(sp) {
  tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & tchol_curr_xps > tchol_target, tchol_curr_xps := tchol_curr_xps + (tchol_target - tchol_curr_xps) * change_10pc]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "tchol_10pc")



# Alcohol ----
# alcohol_target <- 0
#
# scenario_fn_primary_prevention <- function(sp) {
#   sp$pop[year >= sc_year & alcohol_curr_xps > alcohol_target, alcohol_curr_xps := alcohol_curr_xps + (alcohol_target - alcohol_curr_xps) * change_10pc]
# }
#
# IMPACTncd$
#   run(1:n_runs, multicore = TRUE, "alcohol_10pc")


# frvg ----

scenario_fn_primary_prevention <- function(sp) {
  sp$pop[year >= sc_year, fruit_curr_xps := as.integer(round(fruit_curr_xps * (1 + change_10pc)))]
  sp$pop[year >= sc_year, veg_curr_xps := as.integer(round(veg_curr_xps * (1 + change_10pc)))]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "frvg_10pc")



# Physical activity ----

scenario_fn_primary_prevention <- function(sp) {
  tt <- sp$pop[year >= sc_year & active_days_curr_xps == 0, .(unique(pid))]
  tt <- tt[as.logical(rbinom(.N, 1, abs(change_10pc))), V1]
  sp$pop[pid %in% tt & year >= sc_year, `:=` (
    active_days_curr_xps = 1,
    met_curr_xps = as.integer(floor(1 * (3L + qbinom(runif(.N), 8, 3/11)) *
                                      (30 + qexp(runif(.N), 1/7)) / 100))
  )]
}
IMPACTncd$
  run(1:n_runs, multicore = TRUE, "pa_10pc")


# Smoking ----


# # only works for decreasing smokers
scenario_fn_primary_prevention <- function(sp) {
  sp$pop[year >= sc_year & smok_status_curr_xps == "4" ,
         hc_eff := rbinom(.N, 1L, change_10pc)]

  sp$pop[, (c("smok_status_curr_xps", "smok_quit_yrs_curr_xps", "smok_dur_curr_xps", "smok_cig_curr_xps")) :=
           simsmok_policy_impact_decr(
             smok_status_curr_xps,
             smok_quit_yrs_curr_xps,
             smok_dur_curr_xps,
             smok_cig_curr_xps,
             pid_mrk,
             hc_eff
           )]
  sp$pop[, hc_eff := NULL]

  sp$pop[year >= sc_year, smok_cig_curr_xps := as.integer(round(smok_cig_curr_xps * (1 - change_10pc)))]

  sp$pop[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]

  sp$pop[year >= sc_year & ets_curr_xps == 1L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change_10pc)),
         ets_curr_xps := 0L]
}


IMPACTncd$
  run(1:n_runs, multicore = TRUE, "smk_10pc")


# All ----
scenario_fn_primary_prevention <- function(sp) {
  # BMI
  bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & bmi_curr_xps > bmi_target, bmi_curr_xps := bmi_curr_xps + (bmi_target - bmi_curr_xps) * change_10pc]

  # SBP
  sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & sbp_curr_xps > sbp_target, sbp_curr_xps := sbp_curr_xps + (sbp_target - sbp_curr_xps) * change_10pc]

  # Total chol
  tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & tchol_curr_xps > tchol_target, tchol_curr_xps := tchol_curr_xps + (tchol_target - tchol_curr_xps) * change_10pc]

  # ALc - not doing for now
#  sp$pop[year >= sc_year & alcohol_curr_xps > alcohol_target, alcohol_curr_xps := alcohol_curr_xps + (alcohol_target - alcohol_curr_xps) * change_10pc]

  # Fr & vg
  sp$pop[year >= sc_year, fruit_curr_xps := as.integer(round(fruit_curr_xps * (1 + change_10pc)))]
  sp$pop[year >= sc_year, veg_curr_xps := as.integer(round(veg_curr_xps * (1 + change_10pc)))]

  # Physical activity
  tt <- sp$pop[year >= sc_year & active_days_curr_xps == 0, .(unique(pid))]
  tt <- tt[as.logical(rbinom(.N, 1, abs(change_10pc))), V1]
  sp$pop[pid %in% tt & year >= sc_year, `:=` (
    active_days_curr_xps = 1,
    met_curr_xps = as.integer(floor(1 * (3L + qbinom(runif(.N), 8, 3/11)) *
                                      (30 + qexp(runif(.N), 1/7)) / 100))
  )]

# Smoking
  sp$pop[year >= sc_year & smok_status_curr_xps == "4",
         hc_eff := rbinom(.N, 1L, change_10pc)]


  sp$pop[, (c("smok_status_curr_xps", "smok_quit_yrs_curr_xps", "smok_dur_curr_xps", "smok_cig_curr_xps")) :=
           simsmok_policy_impact_decr(
             smok_status_curr_xps,
             smok_quit_yrs_curr_xps,
             smok_dur_curr_xps,
             smok_cig_curr_xps,
             pid_mrk,
             hc_eff
           )]
  sp$pop[, hc_eff := NULL]

  sp$pop[year >= sc_year, smok_cig_curr_xps := as.integer(round(smok_cig_curr_xps * (1 - change_10pc)))]

  sp$pop[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]

  sp$pop[year >= sc_year & ets_curr_xps == 1L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change_10pc)),
         ets_curr_xps := 0L]
}



IMPACTncd$
  run(1:n_runs, multicore = TRUE, "all_10pc") $
    export_summaries(multicore = TRUE, type = c("dis_char", "prvl",
                     "incd", "dis_mrtl", "mrtl",
                     "allcause_mrtl_by_dis", "cms"))

