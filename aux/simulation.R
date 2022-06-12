source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = T)

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:10, multicore = TRUE)

source("./validation_internal/internal_validation_plots.R")
