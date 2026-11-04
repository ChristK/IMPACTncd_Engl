# Tests: export_summaries() exports the right iterations, and refuses rather
# than lose data.
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
# Review of the fix (same day) added: (3) re-exporting an iteration whose
# summaries hold a scenario its lifecourse no longer has (pruned) must stop,
# not silently discard them; (4) `contd` without `_contd` columns must not
# force a re-export on every call; (5) a call in the other single_year_of_age
# shape must stop; (6) a caller's own sinks must survive the call.
#
# Method: the REAL export_summaries() body and its real private helpers run
# against fake hive-partitioned output trees in tempdir(), with `self` and
# the other private methods replaced by minimal stubs. The stub worker
# export_summaries_hlpr() records each call (mc number and worker index) as
# a FILE, so the forked multicore path is observed as well.

library(data.table)

sim_gen <- IMPACTncdEngland::Simulation

# Fake output_dir: `lc` = list(mc -> scenarios) for the lifecourse, with the
# extra lifecourse columns `lc_cols`; `smry` = list(marker_dir -> data.frame)
# of existing summaries (one file per mc).
make_tree <- function(lc, smry = list(), lc_cols = list()) {
  out <- tempfile("exp_smry_")
  for (m in names(lc)) for (s in lc[[m]]) {
    d <- file.path(out, "lifecourse", paste0("mc=", m), paste0("scenario=", s))
    dir.create(d, recursive = TRUE)
    arrow::write_parquet(
      do.call(data.frame, c(list(pid = 1:3, year = 13L, wt = 1), lc_cols)),
      file.path(d, "1_lifecourse.parquet")
    )
  }
  for (mk in names(smry)) {
    d <- file.path(out, "summaries", mk)
    dir.create(d, recursive = TRUE)
    df <- smry[[mk]]
    for (m in unique(df$mc)) {
      arrow::write_parquet(df[df$mc == m, , drop = FALSE],
                           file.path(d, paste0(m, "_", mk, ".parquet")))
    }
  }
  out
}

# Run the real export_summaries() on `out`, recording worker calls in `rec`.
run_export <- function(out, type = "le", multicore = FALSE, ncores = 1L,
                       single_year_of_age = FALSE, rec = tempfile("rec_")) {
  dir.create(rec, showWarnings = FALSE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(
    logs = FALSE, clusternumber_export = ncores, output_dir = out
  )))
  priv <- list(
    console_log = NULL,
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
  # (summaries_complete_mc: its predecessor, so the file also runs -- and
  # fails behaviourally -- against an older build used as a negative control)
  for (nm in c("summaries_status", "summaries_complete_mc",
               "check_summary_age_shape", "release_console_log_on_error")) {
    f <- sim_gen$private_methods[[nm]]
    if (!is.null(f)) {
      environment(f) <- env
      priv[[nm]] <- f
    }
  }
  env$private <- priv
  f <- sim_gen$public_methods$export_summaries
  environment(f) <- env
  f(multicore = multicore, type = type, single_year_of_age = single_year_of_age)
  recorded(rec)
}
recorded <- function(rec) {
  calls <- list.files(rec)
  data.table(
    mc = as.integer(sub("^mc(\\d+)_w.*$", "\\1", calls)),
    worker = suppressWarnings(as.integer(sub("^mc\\d+_w", "", calls)))
  )[order(mc)]
}

smry <- function(mcs, scenarios = "sc0", ...) {
  d <- expand.grid(scenario = scenarios, mc = as.integer(mcs),
                   stringsAsFactors = FALSE)
  data.frame(mc = d$mc, scenario = d$scenario, value = 1, ...)
}
with_contd <- list(x_contd = 1.5) # a lifecourse that has a `_contd` column

# ---- 1. fresh export, numbered from 1 --------------------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = "sc0", `4` = "sc0"))
got <- run_export(out)
expect_identical(got$mc, 1:4, info = "fresh: every iteration exported")
expect_identical(sort(got$worker), 1:4, info = "worker index is the position")

# ---- 2. adding iterations: the reported bug --------------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = "sc0", `4` = "sc0"),
                 list(le_scaled_up = smry(1:2)))
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
                 list(le_scaled_up = smry(1:2)))
expect_identical(run_export(out)$mc, 1:2,
                 info = "iterations whose summaries lack sc1 are re-exported")

# ---- 5. ... only where it is missing ---------------------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = c("sc0", "sc1")),
                 list(le_scaled_up = smry(1:3)))
expect_identical(run_export(out)$mc, 3L,
                 info = "only the iteration that gained a scenario")

# ---- 6. a type never exported before --------------------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0"),
                 list(le_scaled_up = smry(1:2)), lc_cols = with_contd)
expect_identical(run_export(out, type = c("le", "contd"))$mc, 1:2,
                 info = "contd has no summaries yet, so every iteration is exported")

# ---- 7. everything already complete -> nothing exported -------------------
out <- make_tree(list(`1` = c("sc0", "sc1"), `2` = c("sc0", "sc1")),
                 list(le_scaled_up = smry(1:2, c("sc0", "sc1")),
                      contd_scaled_up = smry(1:2, c("sc0", "sc1"))),
                 lc_cols = with_contd)
expect_identical(nrow(run_export(out, type = c("le", "contd"))), 0L,
                 info = "complete for every type and scenario: no worker called")

# ---- 8. unknown type is an error (it used to be ignored) -------------------
out <- make_tree(list(`1` = "sc0"))
expect_error(run_export(out, type = c("le", "prevalence")),
             pattern = "Unknown summary type")

# ---- 9. the forked multicore path passes mc numbers too -------------------
if (.Platform$OS.type != "windows") {
  # export_summaries(multicore = TRUE) sets 1-thread data.table/fst/arrow and
  # registers a doParallel backend in THIS process: restore them, so later
  # test files are not run single-threaded on a stale backend.
  dt_threads <- data.table::getDTthreads()
  fst_threads <- fst::threads_fst()
  arrow_cpus <- arrow::cpu_count()
  out <- make_tree(list(`1` = "sc0", `2` = "sc0", `5` = "sc0", `6` = "sc0"),
                   list(le_scaled_up = smry(c(1L, 2L))))
  got <- run_export(out, multicore = TRUE, ncores = 2L)
  expect_identical(got$mc, c(5L, 6L), info = "multicore: mc 5 and 6, not 1 and 2")
  data.table::setDTthreads(dt_threads)
  fst::threads_fst(fst_threads)
  arrow::set_cpu_count(arrow_cpus)
  foreach::registerDoSEQ()
}

# ---- 10. an unreadable marker means re-export everything, with a warning ---
if (!is.null(sim_gen$private_methods$summaries_status)) {
  out <- make_tree(list(`1` = "sc0", `2` = "sc0"))
  d <- file.path(out, "summaries", "le_scaled_up")
  dir.create(d, recursive = TRUE)
  writeLines("not parquet", file.path(d, "1_le_scaled_up.parquet"))
  expect_warning(got <- run_export(out), pattern = "Could not read")
  expect_identical(got$mc, 1:2, info = "unreadable marker: export all")
}

# ---- 11. a summarised scenario whose lifecourse is gone: refuse ------------
# sc0 was summarised for mc 1:2 and its lifecourse then deleted; sc1 has been
# run since. Re-exporting would rewrite mc 1:2 from sc1 alone and lose sc0.
out <- make_tree(list(`1` = "sc1", `2` = "sc1"),
                 list(le_scaled_up = smry(1:2)))
rec <- tempfile("rec_")
expect_error(run_export(out, rec = rec), pattern = "would discard existing summaries")
expect_identical(nrow(recorded(rec)), 0L, info = "refused before any export")
expect_identical(length(list.files(file.path(out, "summaries", "le_scaled_up"))), 2L,
                 info = "existing summaries untouched")
err <- tryCatch({ run_export(out); "" }, error = conditionMessage)
expect_true(grepl("sc0", err, fixed = TRUE) && grepl("del_summaries()", err, fixed = TRUE),
            info = "names the scenario at risk and the way out")
# ... but an iteration that needs no export is not blocked by it.
out <- make_tree(list(`1` = c("sc1"), `2` = c("sc0", "sc1")),
                 list(le_scaled_up = rbind(smry(1, c("sc0", "sc1")), smry(2))))
expect_identical(run_export(out)$mc, 2L,
                 info = "mc 1 is complete (sc0 orphaned but nothing to rewrite); mc 2 lacks sc1")

# ---- 12. contd requested without `_contd` columns converges ----------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0"),
                 list(le_scaled_up = smry(1:2)))
expect_message(got <- run_export(out, type = c("le", "contd")),
               pattern = "no `_contd` columns")
expect_identical(nrow(got), 0L,
                 info = "le is complete; contd has nothing to export: no re-export")

# ---- 13. the other single_year_of_age shape: refuse --------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = "sc0"),
                 list(prvl_scaled_up = smry(1:2, age = 40L)))
expect_error(run_export(out, type = "prvl", single_year_of_age = FALSE),
             pattern = "single_year_of_age = TRUE")
expect_identical(run_export(out, type = "prvl", single_year_of_age = TRUE)$mc, 3L,
                 info = "the matching shape exports the new iteration")

# ---- 14. a caller's own sink survives the call --------------------------------
out <- make_tree(list(`1` = "sc0"))
tf <- tempfile()
depth0 <- sink.number()
sink(tf)
invisible(run_export(out))
depth_after <- sink.number()
sink()
expect_identical(depth_after, depth0 + 1L,
                 info = "export_summaries() unwinds only the sinks it added")
