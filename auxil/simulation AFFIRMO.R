source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design_AFFIRMO.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = TRUE)

scenario_fn <- function(sp) NULL


IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:10, multicore = TRUE, "sc0")

scenario_fn <- function(sp) {
  sc_year <- 23L # The year the change starts
  change <- 0.2 # mmHg

  sp$pop[year >= sc_year, sbp_curr_xps := sbp_curr_xps - change]
}


IMPACTncd$
  run(1:10, multicore = TRUE, "sc1")$

scenario_fn <- function(sp) {
  sc_year <- 23L # The year the change starts
  change <- 20 # mmHg

  sp$pop[year >= sc_year, sbp_curr_xps := sbp_curr_xps - change]
}


IMPACTncd$
  run(1:10, multicore = TRUE, "sc3")$
  export_summaries(multicore = TRUE)



