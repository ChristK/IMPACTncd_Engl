source("./global.R")

IMPACTncd <- Simulation$new("./aux/sim_design_parf.yaml")

# g <- IMPACTncd$get_causal_structure(print_plot = T)

scenario_fn <- function(sp) NULL

IMPACTncd$
  # del_logs()$
  # del_outputs()$
  run(1:200, multicore = TRUE, "")



scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, bmi_curr_xps := 15]
}

IMPACTncd$
    run(1:200, multicore = TRUE, "bmi")


scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, sbp_curr_xps := 90]
}

IMPACTncd$
  run(1:200, multicore = TRUE, "sbp")

scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, tchol_curr_xps := 2]
}

IMPACTncd$
  run(1:200, multicore = TRUE, "tchol")

scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, alcohol_curr_xps := 0]
}

IMPACTncd$
  run(1:200, multicore = TRUE, "alc")

scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, alcohol_curr_xps := 0]
}

IMPACTncd$
  run(1:200, multicore = TRUE, "alc")


scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, fruit_curr_xps := 800L]
  sp$pop[year >= sc_year, veg_curr_xps := 800L]
  }

IMPACTncd$
  run(1:200, multicore = TRUE, "frvg")


scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, `:=` (
    active_days_curr_xps = 7,
    met_curr_xps = 5000)]
}

IMPACTncd$
  run(1:200, multicore = TRUE, "pa")

scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, `:=` (
    smok_status_curr_xps = 1,
    smok_quit_yrs_curr_xps = 0,
    smok_dur_curr_xps = 0,
    smok_cig_curr_xps = 0,
    smok_packyrs_curr_xps = 0)]
  sp$pop[year >= sc_year,
         ets_curr_xps := 0L]
}

IMPACTncd$
  run(1:200, multicore = TRUE, "smk")


scenario_fn <- function(sp) {
  sc_year <- 13L # The year the change starts
  sp$pop[year >= sc_year, bmi_curr_xps := 15]
  sp$pop[year >= sc_year, sbp_curr_xps := 90]
  sp$pop[year >= sc_year, tchol_curr_xps := 2]
  sp$pop[year >= sc_year, alcohol_curr_xps := 0]
  sp$pop[year >= sc_year, alcohol_curr_xps := 0]
  sp$pop[year >= sc_year, fruit_curr_xps := 800L]
  sp$pop[year >= sc_year, veg_curr_xps := 800L]
  sp$pop[year >= sc_year, `:=` (active_days_curr_xps = 7,
                                met_curr_xps = 5000)]
  sp$pop[year >= sc_year, `:=` (smok_status_curr_xps = 1,
                                smok_quit_yrs_curr_xps = 0,
                                smok_dur_curr_xps = 0,
                                smok_cig_curr_xps = 0,
                                smok_packyrs_curr_xps = 0)]
  sp$pop[year >= sc_year, ets_curr_xps := 0L]
}



IMPACTncd$
  run(1:200, multicore = TRUE, "all")$
  export_summaries(multicore = TRUE, type = c("incd", "prvl", "mrtl"))

# Incd not standardised----

tt <- fread("/mnt/storage_fast/output/hf_real_parf/summaries/incd_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year (not standardised).csv")



outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-dimd (not standardised).csv")



# IMPACTncd$export_summaries(multicore = TRUE)

# source("./aux/CPRD_sim_validation_plots.R")

# IMPACTncd$export_summaries(multicore = TRUE)

# rr <- fread("/mnt/storage_fast/output/hf_real_parf/lifecourse/12_lifecourse.csv.gz")
# rr[, table(scenario)]
