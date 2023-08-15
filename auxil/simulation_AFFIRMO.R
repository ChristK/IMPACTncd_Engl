source("./global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design_AFFIRMO.yaml")

# g <- IMPACTncd$get_causal_structure(processed = FALSE, print_plot = TRUE)
# g <- IMPACTncd$get_causal_structure(processed = FALSE, print_plot = FALSE)

# plot(igraph::make_ego_graph(g, order = 1, c("af"), "all")[[1]])
# plot(igraph::make_ego_graph(g, order = 1, c("af"), "in")[[1]])
# plot(igraph::make_ego_graph(g, order = 1, c("af"), "out")[[1]])

scenario_fn_primary_prevention   <- function(sp) NULL
scenario_fn_secondary_prevention <- function(sp) NULL

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:20, multicore = TRUE, "sc0")



# 20% reduction in AF case fatality
scenario_fn_primary_prevention   <- function(sp) NULL
scenario_fn_secondary_prevention <- function(sp) {
  sc_year <- 23L # The year the change starts
  change <- 0.2 # 20% reduction in mortality

  sp$pop[year >= sc_year, `:=` (prb_af_mrtl1 = prb_af_mrtl1 * (1 - change),
                                prb_af_mrtl2 = prb_af_mrtl2 * (1 - change))]
  NULL
}

IMPACTncd$
  run(1:20, multicore = TRUE, "sc1")


# 20% reduction in SBP/BMI/Tchol among those with prevalent AF
scenario_fn_secondary_prevention <- function(sp) NULL
scenario_fn_primary_prevention <- function(sp) {
  sc_year <- 23L # The year the change starts
  change <- 0.2 # 20% reduction

    # load lifecourse from sc0 and select person years with disease
  pth <- file.path(IMPACTncd$design$sim_prm$output_dir, "lifecourse", paste0(sp$mc_aggr, "_lifecourse.csv.gz"))
  lc <- fread(pth, select = c("pid", "scenario", "year", "af_prvl"))
  # because people may die in sc0 but survive longer in this sc the year in the
  # lc 0 is not accuraye. I will recreate the person-years to bypass mortality
  # in sc0
  pids <- lc[year >= sc_year & af_prvl > 1 & scenario == "sc0", .(pid, year)][, .(year = min(year)), keyby = pid] # incd & prvl cases
  pids <- pids[rep.int(1:.N, times = IMPACTncd$design$sim_prm$init_year + IMPACTncd$design$sim_prm$sim_horizon_max - year + 1L)]
  pids[, year := year + (0:(.N - 1L)), keyby = pid]


  sp$pop[pids, on = .(year, pid), `:=` (
    sbp_curr_xps = sbp_curr_xps * (1 - change),
    bmi_curr_xps = bmi_curr_xps * (1 - change),
    tchol_curr_xps = tchol_curr_xps * (1 - change)
  )]
  NULL
}

IMPACTncd$
  run(1:20, multicore = TRUE, "sc2")


# 20% reduction in stroke risk among those with prevalent AF
scenario_fn_primary_prevention   <- function(sp) NULL
scenario_fn_secondary_prevention <- function(sp) {
  sc_year <- 23L # The year the change starts
  change <- 0.2 # 20% reduction

  sp$pop[year >= sc_year, `:=` (
    stroke_incd_af_prvl_mltp = ((stroke_incd_af_prvl_mltp - 1) * (1 - change)) + 1
  )]
  NULL
}

IMPACTncd$
  run(1:20, multicore = TRUE, "sc3")


# 20% reduction in stroke risk & 20% reduction in SBP/BMI/Tchol among those with prevalent AF
scenario_fn_primary_prevention <- function(sp) {
  sc_year <- 23L # The year the change starts
  change <- 0.2 # 20% reduction
  
  # load lifecourse from sc0 and select person years with disease
  pth <- file.path(IMPACTncd$design$sim_prm$output_dir, "lifecourse", paste0(sp$mc_aggr, "_lifecourse.csv.gz"))
  lc <- fread(pth, select = c("pid", "scenario", "year", "af_prvl"))
  # because people may die in sc0 but survive longer in this sc the year in the
  # lc 0 is not accuraye. I will recreate the person-years to bypass mortality
  # in sc0
  pids <- lc[year >= sc_year & af_prvl > 1 & scenario == "sc0", .(pid, year)][, .(year = min(year)), keyby = pid] # incd & prvl cases
  pids <- pids[rep.int(1:.N, times = IMPACTncd$design$sim_prm$init_year + IMPACTncd$design$sim_prm$sim_horizon_max - year + 1L)]
  pids[, year := year + (0:(.N - 1L)), keyby = pid]


  sp$pop[pids, on = .(year, pid), `:=` (
    sbp_curr_xps = sbp_curr_xps * (1 - change),
    bmi_curr_xps = bmi_curr_xps * (1 - change),
    tchol_curr_xps = tchol_curr_xps * (1 - change)
  )]
  NULL
}

scenario_fn_secondary_prevention <- function(sp) {
  sc_year <- 23L # The year the change starts
  change <- 0.2 # 20% reduction
  
  # load lifecourse from sc0 and select person years with disease
  pth <- file.path(IMPACTncd$design$sim_prm$output_dir, "lifecourse", paste0(sp$mc_aggr, "_lifecourse.csv.gz"))
  lc <- fread(pth, select = c("pid", "scenario", "year", "af_prvl"))
  # because people may die in sc0 but survive longer in this sc the year in the
  # lc 0 is not accuraye. I will recreate the person-years to bypass mortality
  # in sc0
  pids <- lc[year >= sc_year & af_prvl > 1 & scenario == "sc0", .(pid, year)][, .(year = min(year)), keyby = pid] # incd & prvl cases
  pids <- pids[rep.int(1:.N, times = IMPACTncd$design$sim_prm$init_year + IMPACTncd$design$sim_prm$sim_horizon_max - year + 1L)]
  pids[, year := year + (0:(.N - 1L)), keyby = pid]


  sp$pop[pids, on = .(year, pid), `:=` (
    stroke_incd_af_prvl_mltp = ((stroke_incd_af_prvl_mltp - 1) * (1 - change)) + 1
  )]
  NULL
}

IMPACTncd$
  run(1:20, multicore = TRUE, "sc4")$
  export_summaries(multicore = TRUE)

# process summaries ----
simulationParameters <- read_yaml(base::normalizePath("./inputs/sim_design_AFFIRMO.yaml", mustWork=TRUE))
sSummariesSubDirPath <- file.path(simulationParameters$output_dir,"summaries/")
sTablesSubDirPath <- file.path(simulationParameters$output_dir,"tables/")
output_dir <- simulationParameters$output_dir

# DPP
tt <- fread(file.path(sSummariesSubDirPath, "mrtl_scaled_up.csv.gz"))
tt0 <- tt[scenario == "sc0"]
tt[tt0, on = c("mc", "year", "agegrp", "sex", "dimd"), `:=` (
  popsize0 = i.popsize, 
   all_cause_mrtl0 = i.all_cause_mrtl
)]
tt[, dpp := all_cause_mrtl0 - all_cause_mrtl]
tt[between(year, 23, 42), sum(dpp), keyby = .(scenario, mc)][, as.list(signif(quantile(V1, prob = c(0.5, 0.025, 0.975)), 2)), keyby = scenario]

# CPP Multimorbidity (CMS > 1.5)
tt <- fread(file.path(sSummariesSubDirPath, "prvl_scaled_up.csv.gz"))
tt0 <- tt[scenario == "sc0"]
tt[tt0, on = c("mc", "year", "agegrp", "sex", "dimd"), `:=` (
  cmsmm1.5_prvl0 = i.cmsmm1.5_prvl, 
   stroke_prvl0 = i.stroke_prvl
)]
tt[, mm_cpp := cmsmm1.5_prvl0 - cmsmm1.5_prvl]
tt[between(year, 23, 42), sum(mm_cpp), keyby = .(scenario, mc)][, as.list(signif(quantile(V1, prob = c(0.5, 0.025, 0.975)), 2)), keyby = scenario]

# CyPP Multimorbidity (CMS > 1.5)
tt <- fread(file.path(sSummariesSubDirPath, "prvl_scaled_up.csv.gz"))
tt0 <- tt[scenario == "sc0"]
tt[tt0, on = c("mc", "year", "agegrp", "sex", "dimd"), `:=` (
  cmsmm1.5_prvl0 = i.cmsmm1.5_prvl, 
   stroke_prvl0 = i.stroke_prvl
)]
tt[, mm_cpp := cmsmm1.5_prvl0 - cmsmm1.5_prvl]
tt[between(year, 23, 42), sum(mm_cpp), keyby = .(scenario, mc)][, as.list(signif(quantile(V1, prob = c(0.5, 0.025, 0.975)), 2)), keyby = scenario]

# CyPP Stroke
tt[, stroke_cpp := stroke_prvl0 - stroke_prvl]
tt[between(year, 23, 42), sum(stroke_cpp), keyby = .(scenario, mc)][, as.list(signif(quantile(V1, prob = c(0.5, 0.025, 0.975)), 2)), keyby = scenario]

# AF all cases (case years)
tt0[between(year, 23, 42), sum(af_prvl), keyby = .(scenario, mc)][, as.list(signif(quantile(V1, prob = c(0.5, 0.025, 0.975)), 2)), keyby = scenario]


# CPP Stroke
tt <- fread(file.path(sSummariesSubDirPath, "incd_scaled_up.csv.gz"))
tt0 <- tt[scenario == "sc0"]
tt[tt0, on = c("mc", "year", "agegrp", "sex", "dimd"), `:=` (
  cmsmm1.5_incd0 = i.cmsmm1.5_incd, 
   stroke_incd0 = i.stroke_incd
)]
tt[, stroke_cpp := stroke_incd0 - stroke_incd]
tt[between(year, 23, 42), sum(stroke_cpp), keyby = .(scenario, mc)][, as.list(signif(quantile(V1, prob = c(0.5, 0.025, 0.975)), 2)), keyby = scenario]

# AF new cases
tt0[between(year, 23, 42), sum(af_incd), keyby = .(scenario, mc)][, as.list(signif(quantile(V1, prob = c(0.5, 0.025, 0.975)), 2)), keyby = scenario]
