# Direct verification that the year=init_year-1 buffer in set_init_prvl
# produces byte-identical init_prvl across scenarios that only differ at
# year == init_year.
#
# NOTE: this script depends on a diagnostic dump in
# Simulation_class.R$run_sim() that was REMOVED after the initial
# verification (engine kept clean). To re-run this verification, re-add
# the dump (gated on env var IMPACTNCD_DUMP_INIT_PRVL=1) immediately
# after the second `lapply` over diseases (set_rr / set_incd_prb / ...)
# and before the secondary_prevention_scn RNG isolation block. The dump
# should write sp$pop[year %in% c(iy-1, iy), c("pid", "year",
# "pid_mrk", "age", "bmi_curr_xps", "sbp_curr_xps", "dimd",
# grep("_prvl$", ..., value = TRUE))] to outputs/debug/init_prvl_
# <scenario>_mc<n>.fst. Original verification result (Apr 2026): 4 diffs
# out of 144,000 cells, all in obesity_prvl, all attributable to new
# entrants at year == init_year (pid_mrk == TRUE, no year=iy-1 row).

Sys.setenv(IMPACTNCD_DUMP_INIT_PRVL = "1")

source("./global.R")

library(data.table)
library(fst)

IMPACTncd <- Simulation$new("testing/sim_design_testing.yaml")

IMPACTncd$
  del_logs()$
  del_outputs()$
  run(1:2, multicore = TRUE, "sc0")

# sc1: edit at year == init_year (buffer should absorb)
IMPACTncd$update_primary_prevention_scn(
  function(synthpop) {
    iy <- self$design$sim_prm$init_year
    synthpop$pop[
      year == iy,
      `:=`(
        bmi_curr_xps = bmi_curr_xps * 0.5,
        sbp_curr_xps = sbp_curr_xps * 0.7,
        tchol_curr_xps = tchol_curr_xps * 0.6,
        dimd = "10 least deprived",
        qimd = "5 least deprived"
      )
    ]
  }
)
IMPACTncd$run(1:2, multicore = TRUE, "sc1_init_year_only")

# sc2: same edit at year == init_year - 1 (buffer should propagate)
IMPACTncd$update_primary_prevention_scn(
  function(synthpop) {
    py <- self$design$sim_prm$init_year - 1L
    synthpop$pop[
      year == py,
      `:=`(
        bmi_curr_xps = bmi_curr_xps * 0.5,
        sbp_curr_xps = sbp_curr_xps * 0.7,
        tchol_curr_xps = tchol_curr_xps * 0.6,
        dimd = "10 least deprived",
        qimd = "5 least deprived"
      )
    ]
  }
)
IMPACTncd$run(1:2, multicore = TRUE, "sc2_prev_year_only")

# === Compare pre-simcpp init_prvl dumps directly ===
out_dir <- IMPACTncd$design$sim_prm$output_dir
debug_dir <- file.path(out_dir, "debug")

load_dump <- function(scenario_nam, mc_) {
  f <- file.path(debug_dir, sprintf("init_prvl_%s_mc%d.fst", scenario_nam, mc_))
  read_fst(f, as.data.table = TRUE)
}

mc_iters <- c(1L, 2L)

cat("\n=== Pre-simcpp init_prvl: per-MC byte-identity check ===\n")
for (mc_ in mc_iters) {
  d0 <- load_dump("sc0", mc_)
  d1 <- load_dump("sc1_init_year_only", mc_)
  d2 <- load_dump("sc2_prev_year_only", mc_)

  cat(sprintf(
    "mc=%d  sc0 vs sc1: %s   sc0 vs sc2: %s\n",
    mc_,
    if (identical(d0, d1)) "IDENTICAL" else "DIFFERS",
    if (identical(d0, d2)) "IDENTICAL" else "DIFFERS"
  ))
}

cat("\n=== Per-disease difference counts (sc0 vs sc1) ===\n")
agg_diffs <- function(scenario_nam, label) {
  prvl_cols <- NULL
  out <- rbindlist(lapply(mc_iters, function(mc_) {
    d0 <- load_dump("sc0", mc_)
    dx <- load_dump(scenario_nam, mc_)
    cols <- setdiff(names(d0), "pid")
    if (is.null(prvl_cols)) prvl_cols <<- cols
    setkey(d0, pid)
    setkey(dx, pid)
    stopifnot(identical(d0$pid, dx$pid))
    data.table(
      mc = mc_,
      disease = cols,
      n_diff = vapply(
        cols,
        function(c) sum(d0[[c]] != dx[[c]]),
        integer(1)
      )
    )
  }))
  out[, .(total_diff = sum(n_diff)), keyby = disease][order(-total_diff)]
}

cat("\nsc1_init_year_only (buffer should absorb): top differences\n")
print(agg_diffs("sc1_init_year_only", "sc1")[total_diff > 0L])
print(head(agg_diffs("sc1_init_year_only", "sc1"), 20))

cat("\nsc2_prev_year_only (buffer should propagate): top differences\n")
print(head(
  agg_diffs("sc2_prev_year_only", "sc2")[order(-total_diff)],
  20
))

cat("\nVerification complete.\n")
