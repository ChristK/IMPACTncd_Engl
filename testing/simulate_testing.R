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
    synthpop$pop[year >= 25L, bmi_curr_xps := bmi_curr_xps * 0.8]
  }
)

IMPACTncd$update_secondary_prevention_scn(
  function(synthpop) {
    sc_year <- 23L # The year the change starts
    change <- 0.2 # 20% reduction in mortality

    synthpop$pop[
      year >= sc_year,
      `:=`(
        prb_af_mrtl1 = prb_af_mrtl1 * (1 - change),
        prb_af_mrtl2 = prb_af_mrtl2 * (1 - change)
      )
    ]
  }
)

print("Simulating sc1...")

IMPACTncd$
  run(1:2, multicore = TRUE, "sc1")

IMPACTncd$export_summaries(multicore = TRUE)
IMPACTncd$export_tables()

print("Simulation has finished!")