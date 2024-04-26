source("global.R")

IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")

# IMPACTncd$update_output_path("/mnt/storage_fast/output/test")
# IMPACTncd$update_synthpop_path("/mnt/storage_fast/synthpop/test")
# g <- IMPACTncd$get_causal_structure(print_plot = TRUE)
# g <- IMPACTncd$get_causal_structure(processed = FALSE, print_plot = TRUE, focus = "chd")

#plot(igraph::make_ego_graph(g, order = 1, c("pain"), "in")[[1]])

n_runs <- 2L

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:n_runs, multicore = TRUE, "sc0")


IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    change_10pc <- 0.1
    sc_year <- 23L
    sp$pop[year >= sc_year,
           bmi_curr_xps := bmi_curr_xps * (1 - change_10pc)]
    }
  )

IMPACTncd$update_secondary_prevention_scn(
  function(sp) {
    change_10pc <- 0.1
    sc_year <- 23L
    sp$pop[year >= sc_year, `:=`(
      prb_af_mrtl1 = prb_af_mrtl1 * (1 - change_10pc),
      prb_af_mrtl2 = prb_af_mrtl2 * (1 - change_10pc)
    )]
  }
)

IMPACTncd$
  run(1:n_runs, multicore = TRUE, "sc1")


IMPACTncd$
  export_summaries(multicore = TRUE)

source("./auxil/process_out_for_HF.R")

