source("global.R")

IMPACTncd <- Simulation$new("./scenarios/sim_design_ineq.yaml")

# IMPACTncd$update_output_path("/mnt/storage_fast/output/test")
# IMPACTncd$update_synthpop_path("/mnt/storage_fast/synthpop/test")
# g <- IMPACTncd$get_causal_structure(print_plot = TRUE)
# g <- IMPACTncd$get_causal_structure(processed = FALSE, print_plot = TRUE, focus = "chd")

#plot(igraph::make_ego_graph(g, order = 1, c("pain"), "in")[[1]])

n_runs <- 10L

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:n_runs, multicore = TRUE, "sc0")

# Everyone goes to 3 ----
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sc_year <- 26L
    sp$pop[year >= sc_year, `:=`(dimd = 5, qimd = 3)]

    # check Rpackage/IMPACTncd_Engl_model_pkg/R/SynthPop_class.R
    tbl <-
      read_fst(
        "./inputs/exposure_distributions/income_table.fst",
        as.data.table = TRUE
      )
    nam <- intersect(names(dt), names(tbl))
    sp$pop[
      tbl,
      income := (rank_income > inc1) +
        (rank_income > inc2) +
        (rank_income > inc3) +
        (rank_income > inc4) +
        1L,
      on = nam
    ]
    sp$pop[,
      income := factor(
        income,
        levels = 1:5,
        labels = c("1 Highest", "2", "3", "4", "5 Lowest")
      )
    ]

    tbl <-
      read_fst(
        "./inputs/exposure_distributions/active_days_table.fst",
        as.data.table = TRUE
      )
    nam <- intersect(names(dt), names(tbl))
    sp$pop[
      tbl,
      active_days := (rank_pa > pa0) +
        (rank_pa > pa1) +
        (rank_pa > pa2) +
        (rank_pa > pa3) +
        (rank_pa > pa4) +
        (rank_pa > pa5) +
        (rank_pa > pa6),
      on = nam
    ]

                sp$pop[,
                  met := as.integer(floor(
                    active_days *
                      (3L + qbinom(rankstat_pa_met, 8, 3 / 11)) *
                      (30 + qexp(rankstat_pa_dur, 1 / 7)) /
                      100
                  ))
                ] # TODO make data driven

    tbl <-
      read_fst(
        "./inputs/exposure_distributions/frtpor_table.fst",
        as.data.table = TRUE
      )
    col_nam <-
      setdiff(names(tbl), intersect(names(dt), names(tbl)))
    lookup_dt(dt, tbl, check_lookup_tbl_validity = design_$sim_prm$logs)
    sp$pop[,
      fruit_curr_xps := my_qZISICHEL(rank_fruit, mu, sigma, nu, tau, n_cpu = 1L) * 80L
    ] # g/d
    

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

