# Nottingham Local Authority area simulation.
#
# Ported from main's auxil/simulation_NotinghamLA.R into the scenarios/ layout
# and the current API. The old auxil/process_out_for_NotinghamLA.R (880 lines of
# CSV-era post-processing) is intentionally NOT ported: the in-package
# export_summaries()/export_tables() (parquet/DuckDB) supersede it.
source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design_NotinghamLA.yaml")

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:100, multicore = TRUE, "sc0")$
  export_summaries(multicore = TRUE)

IMPACTncd$export_tables(multicore = TRUE)

print("Nottingham LA simulation has finished!")
