library(tinytest)
library(yaml)
library(IMPACTncdEngl)

design_path <- system.file("config/default_sim_design.yaml", package = "IMPACTncdEngl")
design_list <- yaml::read_yaml(file_path)

create_design_list <- function(overrides = list()) {
    base <- list(
        simulation_files_overwrite = TRUE,
        sTag = "test",
        bOverwriteFilesOnDeploy = TRUE,
        RootDirPath = "path",
        sToken = "tok",
        locality = "England",
        clusternumber = 1,
        logs = FALSE,
        scenarios = list(sc0 = list()),
        cols_for_output = c("year", "age", "sex"),
        strata_for_output = c("year", "age", "sex"),
        exposures = list(),
        n = 10,
        init_year_long = 2010,
        sim_horizon_max = 2030,
        ageL = 35,
        ageH = 80,
        diseases = list(
            chd = list(name = "chd", meta = list(incidence = list(influenced_by_disease_name = character(0))))
        ),
        maxlag = 5,
        smoking_relapse_limit = 3,
        stochastic = TRUE,
        kismet = 42,
        jumpiness = 0.2,
        statin_adherence = 0.85,
        bpmed_adherence = 0.9,
        decision_aid = FALSE,
        export_xps = TRUE,
        simsmok_calibration = FALSE,
        output_dir = tempdir(),
        synthpop_dir = tempdir(),
        validation = FALSE,
        iteration_n_max = 2,
        n_synthpop_aggregation = 1
    )
    modifyList(base, overrides)
}

# ---- Initialization checks ----
dl <- create_design_list()
design <- Design$new(design_path)

expect_true("init_year" %in% names(design$sim_prm))
expect_equal(design$sim_prm$init_year, 10) # 2010 - 2000
expect_equal(design$sim_prm$sim_horizon_max, 20) # 2030 - 2010

# ---- Test for default values set ----
expect_true(design$sim_prm$national_qimd)
expect_equal(design$sim_prm$init_year_fromGUI, 10)
expect_equal(design$sim_prm$sim_horizon_fromGUI, 20)

# ---- Directory paths normalized ----
expect_true(grepl("^/", design$sim_prm$output_dir))
expect_true(grepl("^/", design$sim_prm$synthpop_dir))

# ---- Disease ordering respected ----
d_order <- names(design$sim_prm$diseases)
expect_true("chd" %in% d_order)

# ---- Cycle detection with no cycles ----
cycles <- design$.__enclos_env__$private$detect_cycles(design$sim_prm)
expect_equal(length(cycles), 0)

# ---- Cycle detection with artificial loop ----
looped <- create_design_list()
looped$diseases <- list(
    A = list(name = "A", meta = list(incidence = list(influenced_by_disease_name = "B"))),
    B = list(name = "B", meta = list(incidence = list(influenced_by_disease_name = "A")))
)
looped_design <- Design$new(looped)
looped_cycles <- looped_design$.__enclos_env__$private$detect_cycles(looped_design$sim_prm)
expect_equal(length(looped_cycles) > 0, TRUE)

# ---- GUI update modifies correct fields ----
gui_input <- list(
    national_qimd_checkbox = FALSE,
    locality_select = "LADs",
    iteration_n_gui = 99,
    iteration_n_final_gui = 9,
    n_gui = 1000,
    n_synthpop_aggregation_gui = 2,
    n_primers_gui = 1,
    cancer_cure_gui = 0.7,
    jumpiness_gui = 0.3,
    statin_adherence_gui = 0.95,
    bpmed_adherence_gui = 0.85,
    decision_aid_gui = TRUE,
    logs_gui = TRUE,
    timeframe_slider = 2010:2030 # fallback for some GUIs
)

# patch fromGUI_timeframe() if needed
fromGUI_timeframe <- function(gui_input) c("init year" = 2015, "horizon" = 10)

design$update_fromGUI(gui_input)

expect_equal(design$sim_prm$national_qimd, FALSE)
expect_equal(design$sim_prm$locality, "LADs")
expect_equal(design$sim_prm$iteration_n, 99)
expect_equal(design$sim_prm$cancer_cure, 0.7)
expect_equal(design$sim_prm$decision_aid, TRUE)

# ---- Output columns updated for local QIMD ----
expect_true("nqimd" %in% design$sim_prm$cols_for_output)
expect_false("lqimd" %in% design$sim_prm$cols_for_output)

# ---- Save to disk round trip ----
yaml_out <- tempfile(fileext = ".yaml")
design$save_to_disk(yaml_out)
roundtrip <- yaml::read_yaml(yaml_out)
expect_true("sTag" %in% names(roundtrip))
unlink(yaml_out)