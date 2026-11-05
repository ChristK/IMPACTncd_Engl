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
# Round two (same day): the orphan check must read EVERY folder of a type,
# an unreadable file must stop the export, the age-shape check must cover
# every requested type (incl. cms), run even when nothing needs exporting,
# use the columns the call would write, and single_year_of_age must be TRUE
# or FALSE. The file saves and restores the session state it changes.
#
# Method: the REAL export_summaries() body and its real private helpers run
# against fake hive-partitioned output trees in tempdir(), with `self` and
# the other private methods replaced by minimal stubs. The stub worker
# export_summaries_hlpr() records each call (mc number and worker index) as
# a FILE, so the forked multicore path is observed as well.

library(data.table)

# export_summaries() sets data.table/fst/arrow threads in THIS session, and
# its multicore path registers a foreach backend and attaches packages:
# save the state at the very top, restore it at the end.
state0 <- list(dt = data.table::getDTthreads(), fst = fst::threads_fst(),
               arrow = arrow::cpu_count(), search = search())

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
                       single_year_of_age = FALSE, rec = tempfile("rec_"),
                       strata = c("scenario", "year", "agegrp", "sex", "dimd")) {
  dir.create(rec, showWarnings = FALSE)
  env <- new.env(parent = asNamespace("IMPACTncdEngland"))
  env$self <- list(design = list(sim_prm = list(
    logs = FALSE, clusternumber_export = ncores, output_dir = out,
    strata_for_output = strata
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
                      contd_scaled_up = smry(1:2, c("sc0", "sc1"), agegrp = "30-34")),
                 lc_cols = with_contd)
expect_identical(nrow(run_export(out, type = c("le", "contd"))), 0L,
                 info = "complete for every type and scenario: no worker called")

# ---- 8. unknown type is an error (it used to be ignored) -------------------
out <- make_tree(list(`1` = "sc0"))
expect_error(run_export(out, type = c("le", "prevalence")),
             pattern = "Unknown summary type")

# ---- 9. the forked multicore path passes mc numbers too -------------------
if (.Platform$OS.type != "windows") {
  out <- make_tree(list(`1` = "sc0", `2` = "sc0", `5` = "sc0", `6` = "sc0"),
                   list(le_scaled_up = smry(c(1L, 2L))))
  got <- run_export(out, multicore = TRUE, ncores = 2L)
  expect_identical(got$mc, c(5L, 6L), info = "multicore: mc 5 and 6, not 1 and 2")
}

# ---- 10. an unreadable summary file: refuse, naming it --------------------
# (Re-exporting everything used to discard the summaries of pruned scenarios
# that the unreadable file's siblings held.)
out <- make_tree(list(`1` = "sc0", `2` = "sc0"),
                 list(le_scaled_up = smry(2)))
bad_file <- file.path(out, "summaries", "le_scaled_up", "1_le_scaled_up.parquet")
writeLines("not parquet", bad_file)
rec <- tempfile("rec_")
err <- tryCatch({ run_export(out, rec = rec); "" }, error = conditionMessage)
expect_true(grepl("cannot read existing summary file", err, fixed = TRUE),
            info = "an unreadable summary file stops the export")
expect_true(grepl("1_le_scaled_up.parquet", err, fixed = TRUE),
            info = "...naming the file")
expect_identical(nrow(recorded(rec)), 0L, info = "...before any export")

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
expect_warning(got <- run_export(out, type = c("le", "contd")),
               pattern = "no `_contd` columns")
expect_identical(nrow(got), 0L,
                 info = "le is complete; contd has nothing to export: no re-export")

# ---- 13. the other single_year_of_age shape: refuse --------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = "sc0"),
                 list(prvl_scaled_up = smry(1:2, age = 40L)))
expect_error(run_export(out, type = "prvl", single_year_of_age = FALSE),
             pattern = "do not have the age columns")
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
unlink(tf)
expect_identical(depth_after, depth0 + 1L,
                 info = "export_summaries() unwinds only the sinks it added")

# ---- 15. the orphan check reads EVERY folder of the type ----------------------
# mc 2's marker file is missing (so mc 2 needs exporting) but mrtl_esp still
# holds mc 2's sc0, whose lifecourse is gone: re-exporting would discard it.
out <- make_tree(list(`1` = "sc1", `2` = "sc1"),
                 list(mrtl_scaled_up = smry(1, "sc1", agegrp = "30-34"),
                      mrtl_esp = smry(2, c("sc0", "sc1"), agegrp = "30-34")))
err <- tryCatch({ run_export(out, type = "mrtl"); "" }, error = conditionMessage)
expect_true(grepl("would discard existing summaries", err, fixed = TRUE) &&
              grepl("mrtl_esp", err, fixed = TRUE),
            info = "an orphan in a sibling folder (mrtl_esp) is caught and named")

# ---- 16. the refusal names the type and the narrower ways out -----------------
out <- make_tree(list(`1` = "sc1"),
                 list(mrtl_scaled_up = smry(1, "sc0", agegrp = "30-34")))
err <- tryCatch({ run_export(out, type = c("mrtl", "prvl")); "" }, error = conditionMessage)
expect_true(grepl("'mrtl'", err, fixed = TRUE), info = "names the affected type")
expect_true(grepl('type = c("prvl")', err, fixed = TRUE),
            info = "offers exporting only the unaffected type")
expect_true(grepl("Do not use del_summaries()", err, fixed = TRUE),
            info = "warns that del_summaries() would delete far more")

# ---- 17. every requested age-shaped type is checked, not just the first ------
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = "sc0"),
                 list(prvl_scaled_up = smry(1:2, agegrp = "30-34"),
                      incd_scaled_up = smry(1:2, age = 40L)))
err <- tryCatch({ run_export(out, type = c("prvl", "incd")); "" }, error = conditionMessage)
expect_true(grepl("incd_scaled_up", err, fixed = TRUE),
            info = "banded call: the single-year incd folder is refused")
err <- tryCatch({ run_export(out, type = c("incd", "prvl"), single_year_of_age = TRUE); "" },
                error = conditionMessage)
expect_true(grepl("prvl_scaled_up", err, fixed = TRUE),
            info = "single-year call: the banded prvl folder is refused, whatever the order of type")

# ---- 18. cms is age-shaped too ------------------------------------------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0"),
                 list(cms_count_scaled_up = smry(1, agegrp = "30-34")))
expect_error(run_export(out, type = "cms", single_year_of_age = TRUE),
             pattern = "cms_count_scaled_up")

# ---- 19. the other shape is refused even when nothing needs exporting --------
out <- make_tree(list(`1` = "sc0", `2` = "sc0"),
                 list(mrtl_scaled_up = smry(1:2, agegrp = "30-34")))
expect_error(run_export(out, type = "mrtl", single_year_of_age = TRUE),
             pattern = "do not have the age columns")

# ---- 20. single_year_of_age must be TRUE or FALSE -----------------------------
out <- make_tree(list(`1` = "sc0"))
expect_error(run_export(out, single_year_of_age = 1), pattern = "must be TRUE or FALSE")
expect_error(run_export(out, single_year_of_age = "TRUE"), pattern = "must be TRUE or FALSE")

# ---- 21. a design whose two modes write the same columns is never refused ----
# (age in strata_for_output, no agegrp: banded and single-year are identical)
out <- make_tree(list(`1` = "sc0", `2` = "sc0", `3` = "sc0"),
                 list(mrtl_scaled_up = smry(1:2, age = 40L)))
got21 <- tryCatch(
  run_export(out, type = "mrtl",
             strata = c("scenario", "year", "age", "sex", "dimd"))$mc,
  error = conditionMessage
)
expect_identical(got21, 3L,
                 info = "no false refusal: the banded call writes the same age column")

# ---- 22. completeness counts EVERY folder of the type ------------------------
# mc 2's mrtl_esp file is missing (e.g. deleted after it was found unreadable):
# mc 2 must be re-exported, not left with a permanent gap.
out <- make_tree(list(`1` = "sc0", `2` = "sc0"),
                 list(mrtl_scaled_up = smry(1:2, agegrp = "30-34"),
                      mrtl_esp = smry(1, agegrp = "30-34")))
expect_identical(run_export(out, type = "mrtl")$mc, 2L,
                 info = "an iteration missing from a sibling folder is re-exported")

# ---- 23. only the folders the model writes are read ---------------------------
out <- make_tree(list(`1` = "sc0"), list(le_scaled_up = smry(1)))
d <- file.path(out, "summaries", "le_scaled_up_old") # a user's backup copy
dir.create(d)
writeLines("not parquet", file.path(d, "1_le_scaled_up.parquet"))
got23 <- tryCatch(nrow(run_export(out)), error = conditionMessage)
expect_identical(got23, 0L, info = "a user's own folder is ignored: no false refusal")

# ---- 24. a design without `scenario` in strata_for_output ---------------------
# Its summaries have no scenario column: checked by mc only, never 'unreadable'.
out <- make_tree(list(`1` = "sc0", `2` = "sc0"),
                 list(le_scaled_up = data.frame(mc = 1L, value = 1)))
got24 <- tryCatch(run_export(out, strata = c("year", "agegrp", "sex"))$mc,
                  error = conditionMessage)
expect_identical(got24, 2L, info = "no scenario axis: only the missing mc 2 is exported")

# ---- 25. a changed strata_for_output is named as a possible cause ------------
out <- make_tree(list(`1` = "sc0", `2` = "sc0"),
                 list(prvl_scaled_up = smry(1, agegrp = "30-34")))
err <- tryCatch({ run_export(out, type = "prvl",
                             strata = c("scenario", "year", "sex", "dimd")); "" },
                error = conditionMessage)
expect_true(grepl("different strata_for_output", err, fixed = TRUE) &&
              grepl("'prvl_scaled_up' has agegrp", err, fixed = TRUE),
            info = "states the columns found and both possible causes")

# ---- 26. an event-filtered folder may lack a scenario: still complete --------
# hle_cmsmm1.5 keeps only onset rows; sc1 had none in mc 1, so that folder holds
# sc0 only. The iteration is complete (its file exists there), so it must not
# be re-exported on every call -- nor refused once sc0's lifecourse is pruned.
out <- make_tree(list(`1` = c("sc0", "sc1")),
                 list(hle_1st_cond_scaled_up = smry(1, c("sc0", "sc1")),
                      hle_cmsmm1.5_scaled_up = smry(1, "sc0")))
expect_identical(nrow(run_export(out, type = "hle")), 0L,
                 info = "a scenario legitimately absent from a filtered folder: complete")
out <- make_tree(list(`1` = "sc1"),
                 list(hle_1st_cond_scaled_up = smry(1, c("sc0", "sc1")),
                      hle_cmsmm1.5_scaled_up = smry(1, "sc0")))
got26 <- tryCatch(nrow(run_export(out, type = "hle")), error = conditionMessage)
expect_identical(got26, 0L, info = "...and no false refusal after sc0's lifecourse is pruned")

# ---- 27. single-year strata never write `age` twice ---------------------------
# (with `age` already in strata_for_output, c(strata, "age") wrote `age_1`)
hlpr_body <- paste(deparse(body(sim_gen$private_methods$export_summaries_hlpr)), collapse = " ")
expect_true(grepl('unique(c(strata_noagegrp, "age"))', hlpr_body, fixed = TRUE),
            info = "strata_age is unique()d")

# ---- restore the session state (see the top) ----------------------------------
data.table::setDTthreads(state0$dt)
fst::threads_fst(state0$fst)
arrow::set_cpu_count(state0$arrow)
foreach::registerDoSEQ()
for (p in setdiff(search(), state0$search)) {
  try(detach(p, character.only = TRUE), silent = TRUE)
}
unlink(list.files(tempdir(), pattern = "^(exp_smry_|rec_)", full.names = TRUE),
       recursive = TRUE)
