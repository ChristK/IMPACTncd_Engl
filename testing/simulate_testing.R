# https://stackoverflow.com/questions/53622354/how-to-debug-line-by-line-rcpp-generated-code-in-windows
# R -d gdb -e "source('debug.r')"
# break simcpp
# run
source("./global.R")
IMPACTncd <- Simulation$new("testing/sim_design_testing.yaml")

# IMPACTncd$del_parfs()
# IMPACTncd$del_synthpops()

# g <- IMPACTncd$get_causal_structure(print_plot = TRUE)
# g <- IMPACTncd$get_causal_structure(processed = FALSE, print_plot = TRUE, focus = "chd")

# plot(igraph::make_ego_graph(g, order = 1, c("chd"), "in")[[1]])

IMPACTncd$
  del_logs()$
  del_outputs()$
  del_synthpops()$
  # del_parfs()$
  run(1:2, multicore = TRUE, "sc0")

# example of primary prevention scenario function
IMPACTncd$update_primary_prevention_scn(
  function(synthpop) {
    synthpop$pop[year >= 25, sbp_curr_xps := sbp_curr_xps * 0.9]
  }
)

IMPACTncd$
  run(1:2, multicore = TRUE, "sc1")


IMPACTncd$export_summaries(multicore = TRUE)

# IMPACTncd$export_tables(baseline_year_for_change_outputs = 2001L)

# source("./auxil/process_out.R")

print("Simulation has finished!")

# sp$pop colnames 
#  [1] "pid"                     "year"                    "age"                     "sex"                     "type"                   
#  [6] "rank_Fruit_vege"         "rankstat_Smoking_act"    "rankstat_Smoking_ex"     "rankstat_Med_HT"         "rankstat_Med_HL"        
# [11] "rankstat_Med_DM"         "rank_PA_days"            "rank_BMI"                "rank_HbA1c"              "rank_LDLc"              
# [16] "rank_SBP"                "rankstat_Smoking_number" "Fruit_vege_curr_xps"     "Smoking_curr_xps"        "Smoking_number_curr_xps"
# [21] "Med_HT_curr_xps"         "Med_HL_curr_xps"         "Med_DM_curr_xps"         "PA_days_curr_xps"        "BMI_curr_xps"           
# [26] "HbA1c_curr_xps"          "LDLc_curr_xps"           "SBP_curr_xps"            "pid_mrk"                 "wt_immrtl"              
# [31] "all_cause_mrtl"          "cms_score"               "cms_count" 
