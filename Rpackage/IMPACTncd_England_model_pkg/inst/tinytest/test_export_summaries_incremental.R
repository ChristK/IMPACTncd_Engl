# Tests: export_summaries() exports the right iterations.
#
# Regression (2026-11-04, found in the vignette audit and reproduced on the
# real model): export_summaries()
#   (1) gave its worker the iteration's POSITION in mc_set (1..k) instead of
#       its mc NUMBER, so after adding iterations it re-exported the first k
#       and never the new ones, and any numbering not starting at 1 was
#       exported under the wrong mc; and
#   (2) skipped an iteration if its mc appeared in ONE marker summary (the
#       first requested type in a fixed precedence), whatever scenarios that
#       summary held -- so run("sc0")$export_summaries() followed by
#       run("sc1")$export_summaries() never summarised sc1, and a type never
#       exported before was skipped too.
#
# Method: the REAL export_summaries() body (and the real private helper
# summaries_complete_mc()) run against fake hive-partitioned output trees in
# tempdir(), with `self`/`private` replaced by minimal stubs. The stub worker
# export_summaries_hlpr() records each call (mc number and worker index) as a
# FILE, so the forked multicore path is observed as well.

library(data.table)

sim_gen <- IMPACTncdEngland::Simulation

# Build a fake output_dir: `lc` = list(mc -> scenarios) for the lifecourse;
# `smry` = list(marker_dir -> data.frame(mc, scenario)) for existing summaries.
make_tree <- function(lc, smry = list()) {
  out <- tempfile("exp_smry_")
  for (m in names(lc)) for (s in lc[[m]]) {
    d <- file.path(out, "lifecourse", paste0("mc=", m), paste0("scenario=", s))
    dir.create(d, recursive = TRUE)
    arrow::write_parquet(data.frame(pid = 1:3, year = 13L, wt = 1),
                         file.path(d, "1_lifecourse.parquet"))
  }
  for (mk in names(smry)) {
    d <- file.path(out, "summaries", mk)
    dir.create(d, recursive = TRUE)
    df <- smry[[mk]]
    for (m in unique(df$mc)) {
      arrow::write_parquet(
        data.frame(mc = as.integer(m), scenario = df$scenario[df$mc == m], value = 1),
        file.path(d, paste0(m, "_", mk, ".parquet"))
      )
    }
  }
  out
}

# Run the real export_summaries() on `out`; return the recorded worker calls.
run_export <- function(out, type = "le", multicore = FALSE, ncores = 1L) {
  rec <- tempfile("rec_")
  dir.create(rec)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(
    logs = FALSE, clusternumber_export = ncores, output_dir = out
  )))
  priv <- list(
    start_console_log = function(...) invisible(NULL),
    stop_console_log = function(...) invisible(NULL),
    time_mark = function(...) invisible(NULL),
    output_dir = function(x = "") file.path(out, x),
    execute_sql = function(con, sql, context = "") DBI::dbExecute(con, sql),
    query_sql = function(con, sql, context = "") DBI::dbGetQuery(con, sql),
    export_summaries_hlpr = function(mcaggr, worker_index = NA, ...) {
      file.create(file.path(rec, sprintf("mc%s_w%s", mcaggr, worker_index)))
      invisible(NULL)
    }
  )
  helper <- sim_gen$private_methods$summaries_complete_mc
  if (!is.null(helper)) {
    environment(helper) <- env
    priv$summaries_complete_mc <- helper
  }
  env$private <- priv
  f <- sim_gen$public_methods$export_summaries
  environment(f) <- env
  f(multicore = multicore, type = type)
  calls <- list.files(rec)
  data.table(
    mc = as.integer(sub("^mc(\\d+)_w.*$", "\\1", calls)),
    worker = suppressWarnings(as.integer(sub("^mc\\d+_w", "", calls)))
  )[order(mc)]
}

sc0 <- function(mcs) data.frame(mc = mcs, scenario = "sc0")

# ---- 1. fresh export, numbered from 1 --------------------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = "sc0", `4` = "sc0"))
got <- run_export(out)
expect_identical(got$mc, 1:4, info = "fresh: every iteration exported")
expect_identical(sort(got$worker), 1:4, info = "worker index is the position")

# ---- 2. adding iterations: the reported bug --------------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = "sc0", `4` = "sc0"),
                 list(le_scaled_up = sc0(1:2)))
got <- run_export(out)
expect_identical(got$mc, 3:4,
                 info = "only the NEW iterations 3:4 (the bug exported 1:2 again)")
expect_identical(sort(got$worker), 1:2,
                 info = "worker index stays positional for the staggered start")

# ---- 3. numbering that does not start at 1 --------------------------------
out <- make_tree(list(`101` = "sc0", `102` = "sc0", `103` = "sc0"))
expect_identical(run_export(out)$mc, 101:103,
                 info = "mc numbers are passed, not positions 1:3")

# ---- 4. a scenario run after the iterations were summarised ---------------
out <- make_tree(list(`1` = c("sc0", "sc1"), `2` = c("sc0", "sc1")),
                 list(le_scaled_up = sc0(1:2)))
expect_identical(run_export(out)$mc, 1:2,
                 info = "iterations whose summaries lack sc1 are re-exported")

# ---- 5. ... only where it is missing ---------------------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = c("sc0", "sc1")),
                 list(le_scaled_up = sc0(1:3)))
expect_identical(run_export(out)$mc, 3L,
                 info = "only the iteration that gained a scenario")

# ---- 6. a type never exported before --------------------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0"),
                 list(le_scaled_up = sc0(1:2)))
expect_identical(run_export(out, type = c("le", "contd"))$mc, 1:2,
                 info = "contd has no summaries yet, so every iteration is exported")

# ---- 7. everything already complete -> nothing exported -------------------
out <- make_tree(list(`1` = c("sc0", "sc1"), `2` = c("sc0", "sc1")),
                 list(le_scaled_up = data.frame(mc = rep(1:2, each = 2),
                                                scenario = c("sc0", "sc1")),
                      contd_scaled_up = data.frame(mc = rep(1:2, each = 2),
                                                   scenario = c("sc0", "sc1"))))
expect_identical(nrow(run_export(out, type = c("le", "contd"))), 0L,
                 info = "complete for every type and scenario: no worker called")

# ---- 8. unknown type is an error (it used to be ignored) -------------------
out <- make_tree(list(`1` = "sc0"))
expect_error(run_export(out, type = c("le", "prevalence")),
             pattern = "Unknown summary type")

# ---- 9. the forked multicore path passes mc numbers too -------------------
if (.Platform$OS.type != "windows") {
  out <- make_tree(list(`1` = "sc0", `2` = "sc0", `5` = "sc0", `6` = "sc0"),
                   list(le_scaled_up = sc0(c(1L, 2L))))
  got <- run_export(out, multicore = TRUE, ncores = 2L)
  expect_identical(got$mc, c(5L, 6L), info = "multicore: mc 5 and 6, not 1 and 2")
}

# ---- 10. an unreadable marker means re-export everything, with a warning ---
if (!is.null(sim_gen$private_methods$summaries_complete_mc)) {
  out <- make_tree(list(`1` = "sc0", `2` = "sc0"))
  d <- file.path(out, "summaries", "le_scaled_up")
  dir.create(d, recursive = TRUE)
  writeLines("not parquet", file.path(d, "1_le_scaled_up.parquet"))
  expect_warning(got <- run_export(out), pattern = "Could not read")
  expect_identical(got$mc, 1:2, info = "unreadable marker: export all")
}
