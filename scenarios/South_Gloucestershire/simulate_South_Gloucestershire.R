# if (requireNamespace("IMPACTncdEngland", quietly = TRUE)) {
#   remove.packages("IMPACTncdEngland")
#   # Also delete snapshot so installLocalPackageIfChanged() will reinstall
#   unlink("./Rpackage/.IMPACTncd_England_model_pkg_snapshot.rds")
# }
source("./global.R")
IMPACTncd <- Simulation$new("scenarios/South_Gloucestershire/sim_design_South_Gloucestershire.yaml")

IMPACTncd$
  del_logs()$
  del_outputs()$
  # del_synthpops()$
  # del_parfs()$
  run(1:100, multicore = TRUE, "sc0")

IMPACTncd$export_summaries(
  multicore = TRUE,
  # type = c("prvl", "incd", "mrtl"),
  single_year_of_age = TRUE
)

IMPACTncd$export_tables(
  multicore = TRUE,
  cea = FALSE,
  equity = FALSE, # sc0-only run: nothing to contrast
  strata = list(
    ons = list(c("year", "age")),
    esp = list(),
    mrtl_ons = list(c("year", "age")), # ← was defaulting to agegrp
    mrtl_esp = list()
  )
)
print("Simulation has finished!")
