# =============================================================================
# Re-export the GENVASC LAs tables with the corrected equity slope-index code.
#
# Why: two defects were fixed in export_equity_tables() /
# calc_equity_slope_indices() (commit 7922334) --
#   1. `all_diseases_sum` summed the design's umbrella conditions (dm, ctdra,
#      cancer) alongside their own components, double-counting ~9.5% of events;
#   2. REI_rel / RII_ratio inverted their sign convention when the fitted
#      benefit was not positive.
# It also adds `total_benefit`, `fit_R2` and an `n_mc` column.
#
# This re-runs the WHOLE export_tables() step, exactly as genvasc_run.R:385
# does (`export_tables(multicore = TRUE)`, i.e. all defaults: cea = TRUE,
# equity = TRUE, 3.5% discounting, default strata). Running only the equity
# task is not an option the API offers, and matching the original call is what
# keeps the non-equity tables byte-identical rather than subtly re-derived
# under different arguments.
#
# The simulation is NOT re-run: this reads the existing summaries under
# `output_dir` and rewrites `output_dir/tables`.
#
# NOTE /mnt/storage_fast and /mnt/storage_fast4 are bind-mounts of the same
# filesystem (verified: same device and inode), so the YAML's
# `/mnt/storage_fast/genvasc_las/output` IS `/mnt/storage_fast4/genvasc_las/output`.
#
# The previous tables were copied to `output/tables_pre_equity_fix` first;
# genvasc_las_reexport_diff.R compares the two.
#
# Run from the project root:
#   Rscript scenarios/polygenic/genvasc_las_reexport_tables.R
# =============================================================================

source("./global.R")

GV_DESIGN <- "scenarios/polygenic/sim_design_genvasc_las.yaml"

IMPACTncd <- Simulation$new(GV_DESIGN)

out_dir <- IMPACTncd$design$sim_prm$output_dir
message("output_dir: ", out_dir)
stopifnot(dir.exists(file.path(out_dir, "summaries")))

message("=== export_tables() starting at ", format(Sys.time()), " ===")
IMPACTncd$export_tables(multicore = TRUE)
message("=== export_tables() finished at ", format(Sys.time()), " ===")
