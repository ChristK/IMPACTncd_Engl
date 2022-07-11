source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = T)

scenario_fn <- function(sp) NULL

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:2, multicore = TRUE, "")

scenario_fn <- function(sp) sp$pop[year >= 19, bmi_curr_xps := bmi_curr_xps * (1 - 0.1)]

IMPACTncd$
  run(1:2, multicore = TRUE, "sc1")$
  export_summaries(multicore = TRUE)

# IMPACTncd$export_summaries(multicore = TRUE)

source("./validation_internal/internal_validation_plots.R")
