# This script is the basis of the model that is presented in the manuscript 
# draft, currently under review at Nature Comms. Model run summer 2024.
# Processing of outputs for paper done in process_out_elasticites_NEW.R & 
# process_out_elasticities_summary_results.R
# SA of no impact on mortality run in simulation_elasticities_nomort.R (requires
# different sim_design.yaml file)

source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design_elast.yaml")
sc0_output_folder <- IMPACTncd$.__enclos_env__$self$design$sim_prm$output_dir
n_runs <- 100


IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:n_runs, multicore = TRUE, "sc0")#$
#export_summaries(multicore = TRUE)





IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    change_10pc <- 0.1
    sc_year <- 23L
    bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr) 
    sp$pop[year >= sc_year & bmi_curr_xps > bmi_target, 
           bmi_curr_xps := bmi_curr_xps + (bmi_target - bmi_curr_xps) * change_10pc]
    
  }
)

IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/bmi_10pc", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "bmi_10pc")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))



# BMI with effect on SBP & tot chol ----

# Effect on Tchol: 10.1210/clinem/dgaa673 Using the pharmacotherapy effect for mean
# reduction in tchol from a 1 unit reduction in BMI after 12 months: -8.7 (–11.7,–5.94)mg/DL
# which is -0.483mmol/L (-0.649, -0.330). For lifestyle interventions, effect on BMI
# not stat significant, but for 1kg weight reduction the effect sizes for pharmaco
# effects v similar to those for lifestyle. Only applying effect to those with 
#BMI >= 25

# Effect on SBP: https://doi.org/10.1111/jch.14661. Using effect on clinic 
# systolic blood pressure (SBP): educed by 5.79 mmHg (95% CI, 3.54–8.05) 
# after a mean body mass index (BMI) reduction of 2.27 kg/m2. Only applying 
# effect to those with BMI >= 25

IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    change_10pc <- 0.1
    sc_year <- 23L
    if (exists(".Random.seed", .GlobalEnv)) { # save rng state
      oldseed <- .GlobalEnv$.Random.seed
    } else {
      oldseed <- NULL
    }
    set.seed(545390654 + sp$mc_aggr) # http://www.cookbook-r.com/Numbers/Saving_the_state_of_the_random_number_generator/
    on.exit({if (!is.null(oldseed)) { # restore rng state
      .GlobalEnv$.Random.seed <- oldseed
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    }
    )
    bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr) 
    #If effect pushes sbp/tchol below theoretical minimum exposure, set it back to that
    sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr) 
    tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr) 
    
    sbp_sd  <- (8.05 - 3.54)/(2*(qnorm(0.975)))
    sbp_effect <- rnorm(1, mean = -5.79, sd= sbp_sd)/2.27 #change provided is for 2.27 decrease in BMI
    sbp_effect[sbp_effect > 0] <- 0 # make sure effect is nor reversed
    tchol_sd <- (11.7 / 38.67 - 5.94 / 38.67) / (2 * (qnorm(0.975)))
    tchol_effect <- rnorm(1, mean = -8.87 / 38.67, tchol_sd)
    tchol_effect[tchol_effect > 0] <- 0 # make sure effect is nor reversed
    sp$pop[year >= sc_year, 
           bmi_chnge :=  fifelse(bmi_curr_xps > bmi_target,
                                 (bmi_target - bmi_curr_xps) * change_10pc, 
                                 0)]
    sp$pop[bmi_curr_xps >= 25 & year >= sc_year, 
           `:=` (sbp_curr_xps = sbp_curr_xps - (bmi_chnge * sbp_effect), # NOTE both bmi_chnge & sbp_effect -ve
                 tchol_curr_xps = tchol_curr_xps - (bmi_chnge * tchol_effect)) # same here
    ]
    sp$pop[year >= sc_year,
           `:=` (bmi_curr_xps =  bmi_curr_xps + bmi_chnge, 
                 sbp_curr_xps = fifelse(sbp_curr_xps < sbp_target, sbp_target, sbp_curr_xps),
                 tchol_curr_xps = fifelse(tchol_curr_xps < tchol_target, tchol_target, tchol_curr_xps))]
    sp$pop[,   bmi_chnge := NULL]
  }
)


IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/bmi_10pc2", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "bmi_10pc2")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle")
  )






# SBP ----

IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    change_10pc <- 0.1
    sc_year <- 23L
    sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr) 
    sp$pop[year >= sc_year & sbp_curr_xps > sbp_target, 
           sbp_curr_xps := sbp_curr_xps + (sbp_target - sbp_curr_xps) * change_10pc]
  }
)


IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/sbp_10pc", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "sbp_10pc")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))



# Total chol ----
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    change_10pc <- 0.1
    sc_year <- 23L
    tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr) 
    sp$pop[year >= sc_year & tchol_curr_xps > tchol_target,
           tchol_curr_xps := tchol_curr_xps + (tchol_target - tchol_curr_xps) * change_10pc]
  }
)


IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/tchol_10pc", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "tchol_10pc")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))


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

IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    change_10pc <- 0.1
    sc_year <- 23L
    sp$pop[year >= sc_year, fruit_curr_xps := as.integer(round(fruit_curr_xps * (1 + change_10pc)))]
    sp$pop[year >= sc_year, veg_curr_xps := as.integer(round(veg_curr_xps * (1 + change_10pc)))]
  }
)

IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/frvg_10pc", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "frvg_10pc")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))


# Physical activity ----
# Doing 10% of people with <7 days activity in every year 
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    change_10pc <- 0.1
    sc_year <- 23L
    if (exists(".Random.seed", .GlobalEnv)) { # save rng state
      oldseed <- .GlobalEnv$.Random.seed
    } else {
      oldseed <- NULL
    }
    set.seed(545390654 + sp$mc)
    on.exit({if (!is.null(oldseed)) { # restore rng state
      .GlobalEnv$.Random.seed <- oldseed
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    }
    )
    tt <- sp$pop[year >= sc_year & active_days_curr_xps < 7, .(year, pid)]
    tt <- tt[as.logical(rbinom(.N, 1, abs(change_10pc))), .(year, pid)]
    sp$pop[tt, on = c("year", "pid"), alter := 1L]
    sp$pop[alter == 1L & year >= sc_year & active_days_curr_xps != 7, `:=` (
      active_days_curr_xps = active_days_curr_xps + 1L,
      met_curr_xps = as.integer(floor((active_days_curr_xps + 1L) * (3L + qbinom(runif(.N), 8, 3/11)) *
                                        (30 + qexp(runif(.N), 1/7)) / 100))
    )]
    sp$pop[, alter := NULL ]
  }
)


IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/pa_10pc", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "pa_10pc")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))

# Smoking ----


# # only works for decreasing smokers
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    change_10pc <- 0.1
    sc_year <- 23L
    if (exists(".Random.seed", .GlobalEnv)) { # save rng state
      oldseed <- .GlobalEnv$.Random.seed
    } else {
      oldseed <- NULL
    }
    set.seed(545390654 + sp$mc)
    on.exit({if (!is.null(oldseed)) { # restore rng state
      .GlobalEnv$.Random.seed <- oldseed
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    }
    )
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
)


IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/smk_10pc", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "smk_10pc")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))


# All ----
IMPACTncd$update_primary_prevention_scn(
  function(sp) { # BMI (without indirect effects)
    change_10pc <- 0.1
    sc_year <- 23L
    if (exists(".Random.seed", .GlobalEnv)) { # save rng state
      oldseed <- .GlobalEnv$.Random.seed
    } else {
      oldseed <- NULL
    }
    set.seed(545390654 + sp$mc)
    on.exit({if (!is.null(oldseed)) { # restore rng state
      .GlobalEnv$.Random.seed <- oldseed
    } else {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    }
    )

    bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr)
    sp$pop[year >= sc_year & bmi_curr_xps > bmi_target, bmi_curr_xps := bmi_curr_xps + (bmi_target - bmi_curr_xps) * change_10pc]
    
    # SBP
    sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr)
    sp$pop[year >= sc_year & sbp_curr_xps > sbp_target, sbp_curr_xps := sbp_curr_xps + (sbp_target - sbp_curr_xps) * change_10pc]
    
    # Total chol
    tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr)
    sp$pop[year >= sc_year & tchol_curr_xps > tchol_target, tchol_curr_xps := tchol_curr_xps + (tchol_target - tchol_curr_xps) * change_10pc]
    
    # ALc - not doing for now
    #  sp$pop[year >= sc_year & alcohol_curr_xps > alcohol_target, alcohol_curr_xps := alcohol_curr_xps + (alcohol_target - alcohol_curr_xps) * change_10pc]
    
    # Fr & vg
    sp$pop[year >= sc_year, fruit_curr_xps := as.integer(round(fruit_curr_xps * (1 + change_10pc)))]
    sp$pop[year >= sc_year, veg_curr_xps := as.integer(round(veg_curr_xps * (1 + change_10pc)))]
    
    # Physical activity
    tt <- sp$pop[year >= sc_year & active_days_curr_xps < 7, .(year, pid)]
    tt <- tt[as.logical(rbinom(.N, 1, abs(change_10pc))), .(year, pid)]
    sp$pop[tt, on = c("year", "pid"), alter := 1L]
    sp$pop[alter == 1L & year >= sc_year & active_days_curr_xps != 7, `:=` (
      active_days_curr_xps = active_days_curr_xps + 1L,
      met_curr_xps = as.integer(floor((active_days_curr_xps + 1L) * (3L + qbinom(runif(.N), 8, 3/11)) *
                                        (30 + qexp(runif(.N), 1/7)) / 100))
    )]
    sp$pop[, alter := NULL ]
    
    
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
)



IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/all_10pc", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "all_10pc")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))



# PARF ----

# bmi
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 23L
    bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr)
    sp$pop[year >= sc_year, bmi_curr_xps := bmi_target] #doesn't matter what these are as below
  }
)

IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/bmi_parf", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "bmi_parf")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))



# sbp
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 23L
    sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr)
    sp$pop[year >= sc_year, sbp_curr_xps := sbp_target]
  }
)


IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/sbp_parf", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "sbp_parf")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))

# tchol
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 23L
    tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr)
    sp$pop[year >= sc_year, tchol_curr_xps := tchol_target]
  }
)


IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/tchol_parf", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "tchol_parf")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))


#
# #alc
# # scenario_fn_primary_prevention <- function(sp) {
# #   sp$pop[year >= sc_year, alcohol_curr_xps := 0]
# # }
# #
# # IMPACTncd$
# #   run(1:n_runs, multicore = TRUE, "alcohol_parf")
#
#
#frvg

IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 23L
    fruit_target <- IMPACTncd$RR$`fruit~chd`$get_ideal_xps_lvl(sp$mc_aggr) * 80
    veg_target <- IMPACTncd$RR$`veg~chd`$get_ideal_xps_lvl(sp$mc_aggr) * 80
    
    sp$pop[year >= sc_year, fruit_curr_xps := fruit_target] #approx 3.959
    sp$pop[year >= sc_year, veg_curr_xps := veg_target] # approx 3.975
  }
)

IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/frvg_parf", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "frvg_parf")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))


# pa
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 23L
    sp$pop[year >= sc_year, `:=` (
      active_days_curr_xps = 7,
      met_curr_xps = 50)] #Is this correct? The max in the baseline is 54 and mean is 6... Should it be 50 because the unit is 100 mins/week?
  }
)


IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/pa_parf", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "pa_parf")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))

# smk
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 23L
    setkey(sp$pop, pid, year)
    
    
    sp$pop[, (c("smok_status_curr_xps", "smok_quit_yrs_curr_xps", "smok_dur_curr_xps")) :=
             simsmok_complete_cessation(
               smok_status_curr_xps,
               smok_quit_yrs_curr_xps,
               smok_dur_curr_xps,
               pid_mrk,
               year,
               sc_year
             )]
    simsmok_complete_cessation_cig(sp$pop)
    
    sp$pop[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]
    
    sp$pop[year >= sc_year,
           ets_curr_xps := 0L]
  }
)

IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/smk_parf", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "smk_parf")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))

# All ----

IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 23L
    setkey(sp$pop, pid, year)
    # BMI
    bmi_target <- IMPACTncd$RR$`bmi~chd`$get_ideal_xps_lvl(sp$mc_aggr)
    sp$pop[year >= sc_year, bmi_curr_xps := bmi_target]
    # SBP
    sbp_target <- IMPACTncd$RR$`sbp~chd`$get_ideal_xps_lvl(sp$mc_aggr)
    sp$pop[year >= sc_year, sbp_curr_xps := sbp_target]
    # Tchol
    tchol_target <- IMPACTncd$RR$`tchol~chd`$get_ideal_xps_lvl(sp$mc_aggr)
    sp$pop[year >= sc_year, tchol_curr_xps := tchol_target]
    # Alcohol
    # sp$pop[year >= sc_year, alcohol_curr_xps := 0]
    # Fruit
    fruit_target <- IMPACTncd$RR$`fruit~chd`$get_ideal_xps_lvl(sp$mc_aggr) * 80
    sp$pop[year >= sc_year, fruit_curr_xps := fruit_target] #approx 3.959
    # Veg
    veg_target <- IMPACTncd$RR$`veg~chd`$get_ideal_xps_lvl(sp$mc_aggr) * 80
    sp$pop[year >= sc_year, veg_curr_xps := veg_target] # approx 3.975
    # PA
    sp$pop[year >= sc_year, `:=` (
      active_days_curr_xps = 7,
      met_curr_xps = 50)] #Is this correct? The max in the baseline is 54 and mean is 6... Should it be 50 because the unit is 100 mins/week?
    # Smk
    sp$pop[, (c("smok_status_curr_xps", "smok_quit_yrs_curr_xps", "smok_dur_curr_xps")) :=
             simsmok_complete_cessation(
               smok_status_curr_xps,
               smok_quit_yrs_curr_xps,
               smok_dur_curr_xps,
               pid_mrk,
               year,
               sc_year
             )]
    simsmok_complete_cessation_cig(sp$pop)
    sp$pop[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps * smok_dur_curr_xps / 20))]
    sp$pop[year >= sc_year,
           ets_curr_xps := 0L]
  }
)


#
IMPACTncd$
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/all_parf", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "all_parf")$
  export_summaries(multicore = TRUE, type = c(
    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))



#source("./auxil/process_out_for_HF.R")
