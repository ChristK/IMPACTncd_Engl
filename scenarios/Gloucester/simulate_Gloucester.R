# if (requireNamespace("IMPACTncdEngland", quietly = TRUE)) {
#   remove.packages("IMPACTncdEngland")
#   # Also delete snapshot so installLocalPackageIfChanged() will reinstall
#   unlink("./Rpackage/.IMPACTncd_England_model_pkg_snapshot.rds")
# }
source("./global.R")
IMPACTncd <- Simulation$new("scenarios/Gloucester/sim_design_Gloucester.yaml")

IMPACTncd$
  del_logs()$
  del_outputs()$
  # del_synthpops()$
  # del_parfs()$
  run(1:100, multicore = TRUE, "sc0")

IMPACTncd$export_summaries(
  multicore = TRUE,
  type = "prvl",
  single_year_of_age = TRUE
)
IMPACTncd$export_tables(multicore = TRUE,
  strata = list(
    ons = list(c("year", "age")),  # only prevalence by year × age
    esp = list()  
    )
  )
print("Simulation has finished!")
