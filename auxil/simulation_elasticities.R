source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design_elasticities.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = T)

scenario_fn_primary_prevention   <- function(sp) NULL
scenario_fn_secondary_prevention <- function(sp) NULL

change1 <- 0.01
change2 <- 0.1

sc_year <- 23L # The year the change starts

n_runs <- 100 #Number of runs; check don't crash and then run 100


IMPACTncd$
  # del_logs()$
  # del_outputs()$
  run(1:n_runs, multicore = TRUE, "")



# BMI ----

scenario_fn_primary_prevention <- function(sp) {
  bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & bmi_curr_xps > bmi_target, bmi_curr_xps := bmi_curr_xps + (bmi_target - bmi_curr_xps) * change1]
  }

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "bmi1")

scenario_fn_primary_prevention <- function(sp) {
  bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & bmi_curr_xps > bmi_target, bmi_curr_xps := bmi_curr_xps + (bmi_target - bmi_curr_xps) * change2]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "bmi2")


scenario_fn_primary_prevention <- function(sp) {
  bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year, bmi_curr_xps := bmi_target] #doesn't matter what these are as below
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "bmi_parf")



# SBP ----

scenario_fn_primary_prevention <- function(sp) {
  sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & sbp_curr_xps > sbp_target, sbp_curr_xps := sbp_curr_xps + (sbp_target - sbp_curr_xps) * change1]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "sbp1")

scenario_fn_primary_prevention <- function(sp) {
  sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & sbp_curr_xps > sbp_target, sbp_curr_xps := sbp_curr_xps + (sbp_target - sbp_curr_xps) * change2]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "sbp2")


scenario_fn_primary_prevention <- function(sp) {
  sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year, sbp_curr_xps := sbp_target]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "sbp_parf")



# Total chol ----
scenario_fn_primary_prevention <- function(sp) {
  tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & tchol_curr_xps > tchol_target, tchol_curr_xps := tchol_curr_xps + (tchol_target - tchol_curr_xps) * change1]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "tchol1")

scenario_fn_primary_prevention <- function(sp) {
  tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year & tchol_curr_xps > tchol_target, tchol_curr_xps := tchol_curr_xps + (tchol_target - tchol_curr_xps) * change2]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "tchol2")


scenario_fn_primary_prevention <- function(sp) {
  tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  sp$pop[year >= sc_year, tchol_curr_xps := tchol_target]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "tchol_parf")



# Alcohol ----
alcohol_target <- 0

scenario_fn_primary_prevention <- function(sp) {
  sp$pop[year >= sc_year & alcohol_curr_xps > alcohol_target, alcohol_curr_xps := alcohol_curr_xps + (alcohol_target - alcohol_curr_xps) * change1]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "alcohol1")

scenario_fn_primary_prevention <- function(sp) {
  sp$pop[year >= sc_year & alcohol_curr_xps > alcohol_target, alcohol_curr_xps := alcohol_curr_xps + (alcohol_target - alcohol_curr_xps) * change2]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "alcohol2")


scenario_fn_primary_prevention <- function(sp) {
  sp$pop[year >= sc_year, alcohol_curr_xps := alcohol_target]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "alcohol_parf")



# frvg ----

scenario_fn_primary_prevention <- function(sp) {
     sp$pop[year >= sc_year, fruit_curr_xps := as.integer(round(fruit_curr_xps * (1 + change1)))]
     sp$pop[year >= sc_year, veg_curr_xps := as.integer(round(veg_curr_xps * (1 + change1)))]
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "frvg1")

scenario_fn_primary_prevention <- function(sp) {
  sp$pop[year >= sc_year, fruit_curr_xps := as.integer(round(fruit_curr_xps * (1 + change2)))]
  sp$pop[year >= sc_year, veg_curr_xps := as.integer(round(veg_curr_xps * (1 + change2)))]}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "frvg2")


scenario_fn_primary_prevention <- function(sp) {
  fruit_target <- IMPACTncd$RR$`fruit~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...
  veg_target <- IMPACTncd$RR$`veg~chd`$get_ideal_xps_lvl(sp$mc_aggr) #may need to be self$RR...

  sp$pop[year >= sc_year, fruit_curr_xps := fruit_target] #approx 3.959
  sp$pop[year >= sc_year, veg_curr_xps := veg_target] # approx 3.975
}

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "frvg_parf")



# Physical activity ----
active_days_target <- 7
met_curr_xps <- 5000

scenario_fn_primary_prevention <- function(sp) {
  tt <- sp$pop[year >= sc_year & active_days_curr_xps == 0, .(unique(pid))]
  tt <- tt[as.logical(rbinom(.N, 1, abs(change1))), V1]
  sp$pop[pid %in% tt & year >= sc_year, `:=` (
       active_days_curr_xps = 1,
       met_curr_xps = as.integer(floor(1 * (3L + qbinom(runif(.N), 8, 3/11)) *
                                         (30 + qexp(runif(.N), 1/7)) / 100))
       )]
}
IMPACTncd$
  run(1:n_runs, multicore = TRUE, "pa1")

scenario_fn_primary_prevention <- function(sp) {
  tt <- sp$pop[year >= sc_year & active_days_curr_xps == 0, .(unique(pid))]
  tt <- tt[as.logical(rbinom(.N, 1, abs(change1))), V1]
  sp$pop[pid %in% tt & year >= sc_year, `:=` (
    active_days_curr_xps = 1,
    met_curr_xps = as.integer(floor(1 * (3L + qbinom(runif(.N), 8, 3/11)) *
                                      (30 + qexp(runif(.N), 1/7)) / 100))
  )]
}
IMPACTncd$
  run(1:n_runs, multicore = TRUE, "pa2")


scenario_fn_primary_prevention <- function(sp) {
  sp$pop[year >= sc_year, `:=` (
    active_days_curr_xps = 7,
    met_curr_xps = 5000)]
}
IMPACTncd$
  run(1:n_runs, multicore = TRUE, "pa_parf")


# Smoking ----


# only works for decreasing smokers
scenario_fn_primary_prevention <- function(sp) {
    sp$pop[year >= sc_year & smok_status_curr_xps == "4",
       hc_eff := rbinom(.N, 1L, change1)]

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

       sp$pop[year >= sc_year, smok_cig_curr_xps := as.integer(round(smok_cig_curr_xps * (1 - change1)))]

       sp$pop[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]

       sp$pop[year >= sc_year & ets_curr_xps == 1L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change1)),
                   ets_curr_xps := 0L]
}


IMPACTncd$
  run(1:n_runs, multicore = TRUE, "smk1")

scenario_fn_primary_prevention <- function(sp) {
  sp$pop[year >= sc_year & smok_status_curr_xps == "4",
         hc_eff := rbinom(.N, 1L, change2)]

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

  sp$pop[year >= sc_year, smok_cig_curr_xps := as.integer(round(smok_cig_curr_xps * (1 - change2)))]

  sp$pop[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]

  sp$pop[year >= sc_year & ets_curr_xps == 1L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change2)),
         ets_curr_xps := 0L]
}


IMPACTncd$
  run(1:n_runs, multicore = TRUE, "smk2")$
  export_summaries(multicore = TRUE)


# #Parf
# scenario_fn_primary_prevention <- function(sp) {
#   sp$pop[year >= sc_year & smok_status_curr_xps == "4",
#          hc_eff := rbinom(.N, 1L, change)]
#
#   sp$pop[, (c("smok_status_curr_xps", "smok_quit_yrs_curr_xps", "smok_dur_curr_xps", "smok_cig_curr_xps")) :=
#            simsmok_policy_impact_decr(
#              smok_status_curr_xps,
#              smok_quit_yrs_curr_xps,
#              smok_dur_curr_xps,
#              smok_cig_curr_xps,
#              pid_mrk,
#              hc_eff
#            )]
#   sp$pop[, hc_eff := NULL]
#
#   sp$pop[year >= sc_year, smok_cig_curr_xps := as.integer(round(smok_cig_curr_xps * (1 - change)))]
#
#   sp$pop[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]
#
#   sp$pop[year >= sc_year & ets_curr_xps == 1L & smok_status_curr_xps != "4" & rbinom(.N, 1, abs(change)),
#          ets_curr_xps := 0L]
# }


# All ----

# scenario_fn_primary_prevention <- function(sp) {
#   sp$pop[year >= sc_year, bmi_curr_xps := bmi_target]
#   sp$pop[year >= sc_year, sbp_curr_xps := sbp_target]
#   sp$pop[year >= sc_year, tchol_curr_xps := tchol_target]
#   sp$pop[year >= sc_year, alcohol_curr_xps := alcohol_target]
#   sp$pop[year >= sc_year, fruit_curr_xps := fruit_target]
#   sp$pop[year >= sc_year, veg_curr_xps := veg_target]
#   sp$pop[year >= sc_year, `:=` (active_days_curr_xps = active_days_target,
#                                 met_curr_xps = met_target)]
#   sp$pop[year >= sc_year, `:=` (smok_status_curr_xps = smok_status_target,
#                                 smok_quit_yrs_curr_xps = smok_quit_yrs_target,
#                                 smok_dur_curr_xps = smok_dur_target,
#                                 smok_cig_curr_xps = smok_cig_target,
#                                 smok_packyrs_curr_xps = smok_packyrs_target)]
#   sp$pop[year >= sc_year, ets_curr_xps := ets_target]
# }
#
#
#
# IMPACTncd$
#   run(1:n_runs, multicore = TRUE, "all")$
#   export_summaries(multicore = TRUE, type = c("incd", "prvl", "mrtl"))
#







# # Incd not standardised----
# prbl <- c(0.5, 0.025, 0.975, 0.1, 0.9)
# tt <- fread("/mnt/storage_fast/output/hf_real_parf/summaries/incd_scaled_up.csv.gz"
# )[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
#
# outstrata <- c("mc", "year", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", scales::percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, "/mnt/storage_fast/output/hf_real_parf/tables/incidence by year (not standardised).csv")
#
#
#
# outstrata <- c("mc", "year", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", scales::percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, "/mnt/storage_fast/output/hf_real_parf/tables/incidence by year-dimd (not standardised).csv")



# IMPACTncd$export_summaries(multicore = TRUE)

# source("./auxil/CPRD_sim_validation_plots.R")

# IMPACTncd$export_summaries(multicore = TRUE)

# rr <- fread("/mnt/storage_fast/output/hf_real_parf/lifecourse/12_lifecourse.csv.gz")
# rr[, table(scenario)]
# fread("/mnt/storage_fast/output/hf_real_parf/summaries/prvl_scaled_up.csv.gz")[, table(scenario, mc)]
