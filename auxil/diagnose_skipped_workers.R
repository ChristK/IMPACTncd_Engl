#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# diagnose_skipped_workers.R
#
# In a SINGLE-SCENARIO run (only `sc0`), export_tables() writes no CEA and no
# equity tables -- every contrast needs an intervention arm. But both workers
# are still forked, and both READ their summary datasets BEFORE they discover
# there is nothing to contrast:
#
#   export_cea_tables()    Simulation_class_tables.R
#     L2390  qalys <- read_summary_dataset("qalys", "scaled_up")
#     L2391  costs <- read_summary_dataset("costs", "scaled_up")
#     L2423  non_comparator <- setdiff(unique(qalys$scenario), comparator)
#     L2424  if (length(non_comparator) == 0L) return(NULL)      <-- skip, too late
#
#   export_equity_tables()
#     L2667  tt <- read_summary_dataset(cfg$source, "scaled_up") # per metric:
#     L2733  non_comparator <- setdiff(scns, comparator_scenario) # incd, prvl,
#     L2734  if (length(non_comparator) == 0L) { rm(tt); next }   # mrtl, qalys
#
# This script measures what those two doomed workers actually allocate.
#
# Usage: Rscript auxil/diagnose_skipped_workers.R [output_dir]
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

outdir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1L)
  commandArgs(trailingOnly = TRUE)[[1L]] else
  "/mnt/storage_fast/ethnicity_proj/output"

setDTthreads(1L); arrow::set_cpu_count(1L)
sdir <- file.path(outdir, "summaries")

rss <- function() {
  as.numeric(gsub("\\D", "",
    grep("^VmRSS:", readLines("/proc/self/status"), value = TRUE))) / 2^20
}
hwm <- function() {
  as.numeric(gsub("\\D", "",
    grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE))) / 2^20
}
rd <- function(fam) CKutils::read_parquet_dt(file.path(sdir, paste0(fam, "_scaled_up")))

cat("baseline RSS: ", sprintf("%.2f GB\n\n", rss()))

# ---- CEA worker: holds qalys AND costs simultaneously, then bails -----------
cat("== export_cea_tables(): reads qalys + costs, THEN checks scenarios ==\n")
q <- rd("qalys");  cat(sprintf("  qalys loaded  RSS=%.2f GB\n", rss()))
co <- rd("costs"); cat(sprintf("  costs loaded  RSS=%.2f GB\n", rss()))
cat(sprintf("  scenarios present: %s -> nothing to contrast, returns NULL\n",
            paste(unique(q$scenario), collapse = ", ")))
cea_peak <- hwm()
cat(sprintf("  ** CEA worker wasted peak: %.2f GB **\n\n", cea_peak))
rm(q, co); invisible(gc(full = TRUE))

# ---- equity worker: 4 metrics, one dataset at a time, each bails ------------
cat("== export_equity_tables(): reads incd, prvl, mrtl, qalys in turn ==\n")
base <- hwm()
for (m in c(cpp = "incd", cypp = "prvl", dpp = "mrtl", net_qalys = "qalys")) {
  before <- hwm()
  tt <- rd(m)
  cat(sprintf("  %-6s (%-5s): rows=%-11s RSS=%.2f GB -> skipped\n",
              names(which(c(cpp = "incd", cypp = "prvl", dpp = "mrtl",
                            net_qalys = "qalys") == m))[1], m,
              format(nrow(tt), big.mark = ","), rss()))
  rm(tt); invisible(gc(full = TRUE))
}
cat(sprintf("  ** equity worker wasted peak: %.2f GB **\n\n", hwm()))

cat(sprintf("Combined waste from the two no-op workers: ~%.1f GB of the run's peak,\n",
            cea_peak + (hwm() - base)))
cat("plus their share of the 6 concurrent fork slots and their read time.\n")
cat("Both are avoidable with cea = FALSE, equity = FALSE on a 1-scenario run.\n")
