source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design_nomort.yaml")

sc0_output_folder <- IMPACTncd$.__enclos_env__$self$design$sim_prm$output_dir
n_runs <- 100


# IMPACTncd$
#   del_logs()$
#   del_outputs()$
#   run(1:n_runs, multicore = TRUE, "sc0")#$
# #export_summaries(multicore = TRUE)
# 

# All 10pc ----
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
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/no_mort_all_10pc", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "all_10pc")$
  export_summaries(multicore = TRUE, type = c(
    "prvl"
    #    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
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
  update_output_path("/mnt/storage_fast4/output/hf_real_elasticities/no_mort_all_parf", carry_over_lifecourse_files_from = sc0_output_folder)$
  run(1:n_runs, multicore = TRUE, "all_parf")$
  export_summaries(multicore = TRUE, type = c(
     "prvl"
#    "dis_char", "prvl", "incd", "mrtl", "cms", "hle"
  ))
