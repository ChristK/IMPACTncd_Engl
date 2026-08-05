#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# diagnose_export_tables_memory.R
#
# Why does `Simulation$export_tables(multicore = TRUE)` kill the R session (and,
# on a shared box, take the login session down with it)?
#
# `export_tables(multicore = TRUE)` forks up to 6 workers -- one per task:
#   1 main   2 all_cause_mrtl   3 disease_char   4 xps   5 cea   6 equity
# Each worker calls private$read_summary_dataset(), which is
# `CKutils::read_parquet_dt(<dir>)` -- i.e. it materialises the WHOLE summary
# family (all MC iterations, all scenarios) into one data.table, then `copy()`s
# it once per metric inside export_main_tables().
#
# So peak RSS ~= sum over concurrently-running workers of
#      (dataset in RAM) x (1 base + 1 working copy) + group-by intermediates.
#
# This script measures the in-RAM size of each summary family WITHOUT running
# the export, so the peak can be predicted before it bites.
#
# Usage:
#   Rscript auxil/diagnose_export_tables_memory.R <output_dir> [--load-largest]
# e.g.
#   Rscript auxil/diagnose_export_tables_memory.R \
#       /mnt/storage_fast/ethnicity_proj/output
#
# --load-largest additionally reads the single largest family fully into RAM to
# measure the true parquet -> data.table expansion ratio (instead of assuming
# one). Costs a few GB of RAM for ~30 s.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1L]] else
  "/mnt/storage_fast/ethnicity_proj/output"
load_largest <- "--load-largest" %in% args

summaries_dir <- file.path(output_dir, "summaries")
stopifnot(dir.exists(summaries_dir))

# Which summary families each forked task touches. Mirrors
# export_main_tables() / export_all_cause_mrtl_tables() /
# export_disease_characteristics_tables() / export_xps_tables() /
# export_cea_tables() / export_equity_tables() in Simulation_class_tables.R.
task_families <- list(
  main = c("prvl", "incd", "dis_mrtl", "mrtl", "qalys", "costs"),
  all_cause_mrtl = "all_cause_mrtl_by_dis",
  disease_char = "dis_characteristics",
  cea = c("qalys", "costs"),
  equity = c("prvl", "incd", "mrtl", "qalys")
)

# ---- per-family footprint --------------------------------------------------
fams <- sort(basename(list.dirs(summaries_dir, recursive = FALSE)))

measure <- function(fam_dir) {
  p <- file.path(summaries_dir, fam_dir)
  files <- list.files(p, pattern = "\\.parquet$", full.names = TRUE)
  if (length(files) == 0L) return(NULL)
  ds <- arrow::open_dataset(p)
  nr <- nrow(ds)                       # metadata only, no data read
  sch <- ds$schema
  nc <- length(sch$names)
  # Every summary column is numeric/integer once in a data.table; assume 8 B.
  # Keys/factors are cheaper, so this is an upper-ish bound on the *base* table.
  est_ram <- nr * nc * 8 / 2^30
  data.table(
    family    = fam_dir,
    n_files   = length(files),
    disk_gb   = sum(file.size(files)) / 2^30,
    n_rows    = nr,
    n_cols    = nc,
    est_ram_gb = est_ram
  )
}

res <- rbindlist(lapply(fams, measure))
setorder(res, -est_ram_gb)

cat("\n=== Summary families on disk vs estimated in RAM ===\n")
print(res[, .(family, n_files, disk_gb = round(disk_gb, 2), n_rows,
              n_cols, est_ram_gb = round(est_ram_gb, 2))])

# ---- optional: true expansion ratio ---------------------------------------
ratio <- NA_real_
if (load_largest) {
  big <- res[1L, family]
  cat("\nLoading '", big, "' fully to measure the real expansion ratio ...\n",
      sep = "")
  gc(full = TRUE)
  before <- sum(gc()[, "used"] * c(56, 8)) / 2^30
  dt <- CKutils::read_parquet_dt(file.path(summaries_dir, big))
  actual <- as.numeric(object.size(dt)) / 2^30
  ratio <- actual / res[1L, est_ram_gb]
  cat(sprintf("  actual in-RAM: %.2f GB  (estimate %.2f GB, ratio %.2f)\n",
              actual, res[1L, est_ram_gb], ratio))
  rm(dt); gc(full = TRUE)
}

# ---- per-task peak ---------------------------------------------------------
# Standardisations: each family exists as *_scaled_up (ons) and *_esp.
# export_main_tables() holds ONE standardisation at a time, but keeps
# `tt_base` AND a `copy(tt_base)` per metric -> factor 2. The ftlt branch also
# holds prvl on top of dis_mrtl -> that is the worst single moment.
fam_ram <- function(base) {
  x <- res[family %chin% paste0(base, c("_scaled_up", "_esp")), est_ram_gb]
  if (length(x) == 0L) return(0)
  max(x)   # one standardisation resident at a time
}

cat("\n=== Predicted peak RSS per forked worker (GB) ===\n")
peak <- rbindlist(lapply(names(task_families), function(tk) {
  f <- task_families[[tk]]
  # worst single moment within a task: largest family x 2 (base + working copy)
  worst <- max(vapply(f, fam_ram, numeric(1)))
  # the dis_mrtl branch of `main` additionally holds prvl (case-fatality denom)
  extra <- if (tk == "main") fam_ram("prvl") else 0
  data.table(task = tk, peak_gb = round(worst * 2 + extra, 2))
}))
setorder(peak, -peak_gb)
print(peak)

cat(sprintf("\nAll 6 workers fork CONCURRENTLY -> combined peak ~ %.1f GB\n",
            sum(peak$peak_gb)))
cat("(plus the parent's own copy, which fork() shares copy-on-write until the\n",
    " first write touches a page -- data.table's `:=` touches them immediately.)\n",
    sep = "")

cat("\n=== Host headroom right now ===\n")
system("free -g", intern = FALSE)
