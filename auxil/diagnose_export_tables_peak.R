#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# diagnose_export_tables_peak.R
#
# Companion to diagnose_export_tables_memory.R (which sizes the *base* summary
# tables). This one measures the part that actually dominates peak RSS: the
# wide -> long `melt()` inside private$tbl_smmrs_core().
#
# tbl_smmrs_core() does, per (metric x strata) combination:
#   d <- tt[, lapply(.SD, sum), .SDcols = patterns(...), keyby = x]   # x incl. mc, scenario
#   d <- melt(d, id.vars = x)                                        # <-- x n_measure_cols
#   ... cumsum / join to baseline year / quantiles
# so the long table has  nrow(d) x n_measure_cols  rows, each carrying the FULL
# set of id columns. That is where tens of GB appear from a 1 GB parquet file.
#
# Peak RSS is measured for real (via /proc/self/status VmHWM) on ONE
# representative metric+strata, then extrapolated to the concurrent workers.
#
# Usage:
#   Rscript auxil/diagnose_export_tables_peak.R <output_dir> [family] [strata...]
# e.g.
#   Rscript auxil/diagnose_export_tables_peak.R \
#     /mnt/storage_fast/ethnicity_proj/output all_cause_mrtl_by_dis_scaled_up \
#     year agegrp ethnicity dimd
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

args   <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1L) args[[1L]] else
  "/mnt/storage_fast/ethnicity_proj/output"
family <- if (length(args) >= 2L) args[[2L]] else "all_cause_mrtl_by_dis_scaled_up"
strata <- if (length(args) >= 3L) args[-(1:2)] else
  c("year", "agegrp", "ethnicity", "dimd")

setDTthreads(1L)                    # mirror export_tables(multicore = TRUE)
arrow::set_cpu_count(1L)

hwm_gb <- function() {
  x <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("\\D", "", x)) / 2^20
}
rss_gb <- function() {
  x <- grep("^VmRSS:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("\\D", "", x)) / 2^20
}

path <- file.path(outdir, "summaries", family)
stopifnot(dir.exists(path))

cat(sprintf("family : %s\nstrata : %s\n\n", family,
            paste(c("mc", "scenario", strata), collapse = ", ")))

cat("[1] reading the whole family into one data.table ...\n")
t0 <- Sys.time()
tt <- CKutils::read_parquet_dt(path)
cat(sprintf("    rows=%s cols=%s  size=%.2f GB  RSS=%.1f GB  (%.0f s)\n",
            format(nrow(tt), big.mark = ","), ncol(tt),
            as.numeric(object.size(tt)) / 2^30, rss_gb(),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

x <- c("mc", "scenario", strata)
x <- intersect(x, names(tt))
measure_cols <- setdiff(names(tt), c(x, "sex", "age", "agegrp", "dimd",
                                     "ethnicity", "year", "mc", "scenario"))
cat(sprintf("\n[2] group-by %d id cols, %d measure cols\n",
            length(x), length(measure_cols)))

d <- tt[, lapply(.SD, sum), .SDcols = measure_cols, keyby = eval(x)]
cat(sprintf("    aggregated rows=%s  size=%.2f GB  RSS=%.1f GB\n",
            format(nrow(d), big.mark = ","),
            as.numeric(object.size(d)) / 2^30, rss_gb()))

cat("\n[3] melt() wide -> long (this is the expensive step)\n")
dl <- melt(d, id.vars = x)
cat(sprintf("    long rows=%s (= %s x %d)  size=%.2f GB  RSS=%.1f GB\n",
            format(nrow(dl), big.mark = ","), format(nrow(d), big.mark = ","),
            length(measure_cols),
            as.numeric(object.size(dl)) / 2^30, rss_gb()))

cat("\n[4] the *_change_* branch additionally holds a baseline-year copy + join\n")
d19 <- dl[year == min(year)][, year := NULL]
cat(sprintf("    baseline slice rows=%s  size=%.2f GB  RSS=%.1f GB\n",
            format(nrow(d19), big.mark = ","),
            as.numeric(object.size(d19)) / 2^30, rss_gb()))

cat(sprintf("\n=== PEAK for this ONE metric x ONE strata: %.1f GB (VmHWM) ===\n",
            hwm_gb()))
cat("export_tables() runs ~6 such workers CONCURRENTLY, and each iterates over\n")
cat("every metric x strata combination in turn, so this figure is the floor of\n")
cat("what a single worker needs, not the total.\n")
