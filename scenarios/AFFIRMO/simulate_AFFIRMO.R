# AFFIRMO scenario suite: atrial-fibrillation interventions among prevalent-AF
# cases (case-fatality, risk factors, stroke risk, and combinations).
#
# Ported from main's auxil/simulation_AFFIRMO.R to the scenarios/ layout and the
# current API:
#   - old global scenario_fn_primary/secondary_prevention()  -> IMPACTncd$
#     update_primary_prevention_scn() / update_secondary_prevention_scn().
#   - reading the sc0 lifecourse changed from CSV (<mc>_lifecourse.csv.gz) to the
#     partitioned parquet output (lifecourse/mc=<m>/scenario=sc0/). Inside a
#     scenario function the enclosing environment is the Simulation object, so
#     `self$design$sim_prm$...` provides output_dir / horizon (NOT the global
#     `IMPACTncd`, whose binding is dropped when the function's environment is
#     reset).
#   - the trailing CSV-summary analysis block (DPP/CPP/CyPP on *_scaled_up.csv.gz)
#     is dropped; use export_tables()/export_summaries() output instead.
source("./global.R")

IMPACTncd <- Simulation$new("scenarios/AFFIRMO/sim_design_AFFIRMO.yaml")

n_runs <- 1:20

# --- sc0: baseline -----------------------------------------------------------
IMPACTncd$update_primary_prevention_scn(function(sp) NULL)
IMPACTncd$update_secondary_prevention_scn(function(sp) NULL)
IMPACTncd$
  del_logs()$
  del_outputs()$
  run(n_runs, multicore = TRUE, "sc0")

# --- sc1: 20% reduction in AF case fatality ----------------------------------
IMPACTncd$update_primary_prevention_scn(function(sp) NULL)
IMPACTncd$update_secondary_prevention_scn(function(sp) {
  sc_year <- 23L # year the change starts
  change <- 0.2  # 20% reduction in AF mortality
  sp$pop[year >= sc_year, `:=`(
    prb_af_mrtl1 = prb_af_mrtl1 * (1 - change),
    prb_af_mrtl2 = prb_af_mrtl2 * (1 - change)
  )]
  NULL
})
IMPACTncd$run(n_runs, multicore = TRUE, "sc1")

# --- sc2: 20% reduction in SBP/BMI/Tchol among those with prevalent AF -------
IMPACTncd$update_secondary_prevention_scn(function(sp) NULL)
IMPACTncd$update_primary_prevention_scn(function(sp) {
  sc_year <- 23L
  change <- 0.2
  # person-years of prevalent AF in baseline (sc0), expanded from each person's
  # first prevalent year to the horizon (sc0 mortality may differ under the scn).
  lc_dir <- file.path(self$design$sim_prm$output_dir, "lifecourse",
                      paste0("mc=", sp$mc_aggr), "scenario=sc0")
  fls <- list.files(lc_dir, pattern = "\\.parquet$", full.names = TRUE, recursive = TRUE)
  lc <- data.table::rbindlist(lapply(fls, function(f)
    data.table::as.data.table(arrow::read_parquet(f, col_select = c("pid", "year", "af_prvl")))))
  pids <- lc[year >= sc_year & af_prvl > 1, .(year = min(year)), keyby = pid]
  horizon <- self$design$sim_prm$init_year + self$design$sim_prm$sim_horizon_max
  pids <- pids[rep.int(seq_len(.N), times = horizon - year + 1L)]
  pids[, year := year + (0:(.N - 1L)), keyby = pid]
  sp$pop[pids, on = .(year, pid), `:=`(
    sbp_curr_xps   = sbp_curr_xps   * (1 - change),
    bmi_curr_xps   = bmi_curr_xps   * (1 - change),
    tchol_curr_xps = tchol_curr_xps * (1 - change)
  )]
  NULL
})
IMPACTncd$run(n_runs, multicore = TRUE, "sc2")

# --- sc3: 20% reduction in stroke risk among those with prevalent AF ---------
# (the AF->stroke multiplier is 1 for non-AF, so a year-wide change only affects
#  prevalent-AF person-years.)
IMPACTncd$update_primary_prevention_scn(function(sp) NULL)
IMPACTncd$update_secondary_prevention_scn(function(sp) {
  sc_year <- 23L
  change <- 0.2
  sp$pop[year >= sc_year,
         stroke_incd_af_prvl_mltp := ((stroke_incd_af_prvl_mltp - 1) * (1 - change)) + 1]
  NULL
})
IMPACTncd$run(n_runs, multicore = TRUE, "sc3")

# --- sc4: combined sc2 (risk factors) + sc3 (stroke risk) among prevalent AF -
IMPACTncd$update_primary_prevention_scn(function(sp) {
  sc_year <- 23L
  change <- 0.2
  lc_dir <- file.path(self$design$sim_prm$output_dir, "lifecourse",
                      paste0("mc=", sp$mc_aggr), "scenario=sc0")
  fls <- list.files(lc_dir, pattern = "\\.parquet$", full.names = TRUE, recursive = TRUE)
  lc <- data.table::rbindlist(lapply(fls, function(f)
    data.table::as.data.table(arrow::read_parquet(f, col_select = c("pid", "year", "af_prvl")))))
  pids <- lc[year >= sc_year & af_prvl > 1, .(year = min(year)), keyby = pid]
  horizon <- self$design$sim_prm$init_year + self$design$sim_prm$sim_horizon_max
  pids <- pids[rep.int(seq_len(.N), times = horizon - year + 1L)]
  pids[, year := year + (0:(.N - 1L)), keyby = pid]
  sp$pop[pids, on = .(year, pid), `:=`(
    sbp_curr_xps   = sbp_curr_xps   * (1 - change),
    bmi_curr_xps   = bmi_curr_xps   * (1 - change),
    tchol_curr_xps = tchol_curr_xps * (1 - change)
  )]
  NULL
})
IMPACTncd$update_secondary_prevention_scn(function(sp) {
  sc_year <- 23L
  change <- 0.2
  lc_dir <- file.path(self$design$sim_prm$output_dir, "lifecourse",
                      paste0("mc=", sp$mc_aggr), "scenario=sc0")
  fls <- list.files(lc_dir, pattern = "\\.parquet$", full.names = TRUE, recursive = TRUE)
  lc <- data.table::rbindlist(lapply(fls, function(f)
    data.table::as.data.table(arrow::read_parquet(f, col_select = c("pid", "year", "af_prvl")))))
  pids <- lc[year >= sc_year & af_prvl > 1, .(year = min(year)), keyby = pid]
  horizon <- self$design$sim_prm$init_year + self$design$sim_prm$sim_horizon_max
  pids <- pids[rep.int(seq_len(.N), times = horizon - year + 1L)]
  pids[, year := year + (0:(.N - 1L)), keyby = pid]
  sp$pop[pids, on = .(year, pid),
         stroke_incd_af_prvl_mltp := ((stroke_incd_af_prvl_mltp - 1) * (1 - change)) + 1]
  NULL
})
IMPACTncd$run(n_runs, multicore = TRUE, "sc4")

IMPACTncd$export_summaries(multicore = TRUE)
IMPACTncd$export_tables(multicore = TRUE)

print("AFFIRMO simulation has finished!")
