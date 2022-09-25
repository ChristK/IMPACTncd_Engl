source("./global.R")

IMPACTncd <- Simulation$new("./aux/sim_design_parf.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = T)

scenario_fn <- function(sp) NULL

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:100, multicore = TRUE, "")



scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, bmi_curr_xps := 15]
}

IMPACTncd$
    run(1:100, multicore = TRUE, "bmi")


scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, sbp_curr_xps := 90]
}

IMPACTncd$
  run(1:100, multicore = TRUE, "sbp")

scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, tchol_curr_xps := 2]
}

IMPACTncd$
  run(1:100, multicore = TRUE, "tchol")

scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, alcohol_curr_xps := 0]
}

IMPACTncd$
  run(1:100, multicore = TRUE, "alc")

scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, alcohol_curr_xps := 0]
}

IMPACTncd$
  run(1:100, multicore = TRUE, "alc")


scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, fruit_curr_xps := 800L]
  sp$pop[year >= sc_year, veg_curr_xps := 800L]
  }

IMPACTncd$
  run(1:100, multicore = TRUE, "frvg")


scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, `:=` (
    active_days_curr_xps = 7,
    met_curr_xps = 5000)]
}

IMPACTncd$
  run(1:100, multicore = TRUE, "pa")

scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, `:=` (
    smok_status_curr_xps = 1,
    smok_quit_yrs_curr_xps = 0,
    smok_dur_curr_xps = 0,
    smok_cig_curr_xps = 0,
    smok_packyrs_curr_xps = 0)]
  sp$pop[year >= sc_year,
         ets_curr_xps := 0L]
}

IMPACTncd$
  run(1:100, multicore = TRUE, "smk")


scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, bmi_curr_xps := 15]
  sp$pop[year >= sc_year, sbp_curr_xps := 90]
  sp$pop[year >= sc_year, tchol_curr_xps := 2]
  sp$pop[year >= sc_year, alcohol_curr_xps := 0]
  sp$pop[year >= sc_year, alcohol_curr_xps := 0]
  sp$pop[year >= sc_year, fruit_curr_xps := 800L]
  sp$pop[year >= sc_year, veg_curr_xps := 800L]
  sp$pop[year >= sc_year, `:=` (active_days_curr_xps = 7,
                                met_curr_xps = 5000)]
  sp$pop[year >= sc_year, `:=` (smok_status_curr_xps = 1,
                                smok_quit_yrs_curr_xps = 0,
                                smok_dur_curr_xps = 0,
                                smok_cig_curr_xps = 0,
                                smok_packyrs_curr_xps = 0)]
  sp$pop[year >= sc_year, ets_curr_xps := 0L]
}



IMPACTncd$
  run(1:100, multicore = TRUE, "all")$
  export_summaries(multicore = TRUE)

# IMPACTncd$export_summaries(multicore = TRUE)

# source("./aux/CPRD_sim_validation_plots.R")

# IMPACTncd$export_summaries(multicore = TRUE)
