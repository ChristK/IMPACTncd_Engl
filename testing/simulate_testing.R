# remove.packages("IMPACTncdEngland")
# # Also delete snapshot so installLocalPackageIfChanged() will reinstall
# unlink("./Rpackage/.IMPACTncd_England_model_pkg_snapshot.rds")
source("./global.R")
IMPACTncd <- Simulation$new("testing/sim_design_testing.yaml")

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:2, multicore = TRUE, "sc0")

# example of primary prevention scenario function
IMPACTncd$update_primary_prevention_scn(
  function(synthpop) {
    synthpop$pop[year >= 25, bmi_curr_xps := bmi_curr_xps * 0.8]
  }
)

print("Simulating sc1...")

IMPACTncd$
  run(1:2, multicore = TRUE, "sc1")

IMPACTncd$export_summaries(multicore = TRUE)
IMPACTncd$export_tables()

print("Simulation has finished!")