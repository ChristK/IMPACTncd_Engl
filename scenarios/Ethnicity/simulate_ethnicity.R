# if (requireNamespace("IMPACTncdEngland", quietly = TRUE)) {
#   remove.packages("IMPACTncdEngland")
#   # Also delete snapshot so installLocalPackageIfChanged() will reinstall
#   unlink("./Rpackage/.IMPACTncd_England_model_pkg_snapshot.rds")
# }
source("./global.R")
IMPACTncd <- Simulation$new("scenarios/Ethnicity/sim_design_ethnicity.yaml")

IMPACTncd$
  del_logs()$
  del_outputs()$
  # del_synthpops()$
  # del_parfs()$
  run(1:100, multicore = TRUE, "sc0")

IMPACTncd$export_summaries(multicore = TRUE)

IMPACTncd$export_tables(
  baseline_year_for_change_outputs = 2026L,
  two_agegrps = FALSE,
  multicore = TRUE,
  strata = list(
    ons = list(
      "year",
      c("year", "ethnicity"),
      c("year", "dimd"),
      c("year", "ethnicity", "dimd"),
      c("year", "agegrp"),
      c("year", "agegrp", "ethnicity"),
      c("year", "agegrp", "ethnicity", "dimd")
    ),
    esp = list(
      "year",
      c("year", "ethnicity"),
      c("year", "dimd"),
      c("year", "ethnicity", "dimd")
    ),
    mrtl_ons = list(
      "year",
      c("year", "ethnicity"),
      c("year", "dimd"),
      c("year", "agegrp"),
      c("year", "agegrp", "ethnicity"),
      c("year", "agegrp", "ethnicity", "dimd")
    ),
    mrtl_esp = list(
      "year",
      c("year", "ethnicity"),
      c("year", "dimd"),
      c("year", "ethnicity", "dimd")
    ),
    disease_char = list(
      "year",
      c("year", "ethnicity"),
      c("year", "dimd"),
      c("year", "ethnicity", "dimd")
    ),
    xps_ons = list(
      "year",
      c("year", "ethnicity"),
      c("year", "agegrp20"),
      c("year", "agegrp20", "ethnicity"),
      c("year", "agegrp20", "ethnicity", "qimd")
    ),
    xps_esp = list(
      "year",
      c("year", "ethnicity"),
      c("year", "qimd"),
      c("year", "ethnicity", "qimd")
    )
  )
)
IMPACTncd$export_tables(
  baseline_year_for_change_outputs = 2026L,
  two_agegrps = TRUE,
  multicore = TRUE,
  strata = list(
    ons = list(
      c("year", "agegrp"),
      c("year", "agegrp", "ethnicity"),
      c("year", "agegrp", "ethnicity", "dimd")
    ),
    esp = list(
      "year",
      c("year", "ethnicity"),
      c("year", "dimd"),
      c("year", "ethnicity", "dimd")
    ),
    mrtl_ons = list(
      c("year", "agegrp"),
      c("year", "agegrp", "ethnicity"),
      c("year", "agegrp", "ethnicity", "dimd")
    ),
    mrtl_esp = list(
      "year",
      c("year", "ethnicity"),
      c("year", "dimd"),
      c("year", "ethnicity", "dimd")
    ),
    disease_char = list(
      "year",
      c("year", "ethnicity"),
      c("year", "dimd"),
      c("year", "ethnicity", "dimd")
    ),
    xps_ons = list(
      c("year", "agegrp20"),
      c("year", "agegrp20", "ethnicity"),
      c("year", "agegrp20", "ethnicity", "qimd")
    ),
    xps_esp = list(
      "year",
      c("year", "ethnicity"),
      c("year", "qimd"),
      c("year", "ethnicity", "qimd")
    )
  )
)

print("Simulation has finished!")
