# =============================================================================
# Quantify a PRE-EXISTING bug in tbl_smmrs_core()'s cypp/cpp/dpp block.
#
# The contrast is an update-join whose target is `value` ITSELF:
#     d[d_sc0, on = c(setdiff(x, "scenario"), "variable"), value := i.value - value]
# data.table's x[i, :=] updates only the rows of `x` that `i` matches. An
# intervention row with NO comparator counterpart is therefore left UNCHANGED --
# it keeps its raw LEVEL (a prevalence or incidence count) and is published in a
# file called "cases prevented or postponed".
#
# That happens whenever the two arms do not have identical cell coverage, which
# is normal at fine strata: differential survival means a (year, agegrp, sex,
# dimd) cell can be occupied in one arm and empty in the other.
#
# net_qalys / net_costs are NOT affected: they write a NEW column, so an
# unmatched row gets NA, which cumsum() propagates and the quantile step drops.
#
# Run from the project root:
#   Rscript testing/diagnose_unmatched_comparator_cells.R [output_dir]
# =============================================================================

suppressMessages({
  library(data.table)
  library(CKutils)
})

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args)) args[[1L]] else "/mnt/storage_fast4/genvasc_las/output"
COMPARATOR <- "sc0"
BASELINE <- 2019L

# The default non-standardised strata, from build_strata_config().
STRATA <- list("year", c("year", "sex"), c("year", "agegrp"),
               c("year", "dimd"), c("year", "agegrp", "sex"),
               c("year", "agegrp", "sex", "dimd"))

# metric -> (summary, column pattern)
CFG <- list(
  cypp = list(src = "prvl", pat = "_prvl$|^popsize$"),
  cpp  = list(src = "incd", pat = "_incd$|^popsize$"),
  dpp  = list(src = "mrtl", pat = "_mrtl$|^popsize$")
)

tot <- data.table()
for (metric in names(CFG)) {
  cfg <- CFG[[metric]]
  p <- file.path(OUT, "summaries", paste0(cfg$src, "_scaled_up"))
  if (!dir.exists(p)) { cat("skip", metric, "- no", cfg$src, "summary\n"); next }
  tt <- read_parquet_dt(p)
  setDT(tt)
  tt[, year := year + 2000L]
  scns <- setdiff(unique(tt$scenario), COMPARATOR)

  for (st in STRATA) {
    if (!all(st %chin% names(tt))) next
    x <- c("mc", "scenario", st)
    d <- tt[, lapply(.SD, sum), .SDcols = patterns(cfg$pat), keyby = x]
    d <- melt(d, id.vars = x)
    cmp <- d[scenario == COMPARATOR & year >= BASELINE]
    di  <- d[scenario != COMPARATOR & year >= BASELINE]
    if (!nrow(di) || !nrow(cmp)) next
    on <- c(setdiff(x, "scenario"), "variable")
    u <- di[!cmp, on = on]
    tot <- rbind(tot, data.table(
      metric = metric, strata = paste(st, collapse = "+"),
      rows = nrow(di), unmatched = nrow(u),
      pct = round(100 * nrow(u) / nrow(di), 4),
      # what those rows would be published AS: the raw level
      max_raw = if (nrow(u)) max(u$value, na.rm = TRUE) else NA_real_,
      scenarios = if (nrow(u)) uniqueN(u$scenario) else 0L))
  }
  rm(tt)
  gc(verbose = FALSE)
}

cat("\n=== intervention rows with NO comparator cell ===\n")
cat("(these keep their RAW LEVEL and are published as prevented cases)\n\n")
print(tot)

bad <- tot[unmatched > 0]
cat("\nstrata affected:", nrow(bad), "of", nrow(tot), "\n")
if (nrow(bad)) {
  cat("worst:", bad[which.max(pct), paste0(metric, " by ", strata, " -- ",
      unmatched, " rows (", pct, "%), largest raw value published as a ",
      "prevented count: ", format(round(max_raw), big.mark = ","))], "\n")
  cat("\nNOTE the corruption does not stop at those cells: the value feeds\n",
      "cumsum() within (mc, scenario, strata, variable), so every LATER year\n",
      "in the same series carries the error too.\n", sep = "")
}
