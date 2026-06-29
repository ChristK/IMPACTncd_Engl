# Nottingham ICS area simulation.
#
# Ported from main's auxil/simulation_NotinghamICS.R into the scenarios/ layout
# and the current API. The old auxil/process_out_for_NotinghamICS.R (880 lines of
# CSV-era post-processing) is intentionally NOT ported: the in-package
# export_summaries()/export_tables() (parquet/DuckDB) supersede it.
source("./global.R")

IMPACTncd <- Simulation$new("scenarios/NotinghamICS/sim_design_NotinghamICS.yaml")

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:100, multicore = TRUE, "sc0")$
  export_summaries(multicore = TRUE)

IMPACTncd$export_tables(multicore = TRUE)

print("Nottingham ICS simulation has finished!")
