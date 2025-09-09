source("global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design_NotinghamLA.yaml")

n_runs <- 100L

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:n_runs, multicore = TRUE, "sc0")$
  export_summaries(multicore = TRUE)

# IMPACTncd$export_summaries(multicore = TRUE)
source("./auxil/process_out_for_NotinghamLA.R")

