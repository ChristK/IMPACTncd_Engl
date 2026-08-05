# =============================================================================
# Prove export_tables() now covers the GENVASC CEA case end to end, on the real
# GENVASC LAs summaries -- the comparator map, the programme start year, the
# payer perspectives and the WTP thresholds all at once.
#
# This is the configuration genvasc_run.R now uses. It writes into a SEPARATE
# tables directory so the published `tables/` is untouched.
#
# Run from the project root:
#   Rscript scenarios/polygenic/genvasc_las_cea_via_package.R
# =============================================================================

source("./global.R")
suppressMessages(library(data.table))

IMPACTncd <- Simulation$new("scenarios/polygenic/sim_design_genvasc_las.yaml")
out <- IMPACTncd$design$sim_prm$output_dir

GV_PROGRAMME_START <- 2023L
GV_COST_COLS <- c("gv_prs_costs", "gv_statin_costs", "gv_ahmed_costs",
                  "gv_labs_costs", "gv_clinical_costs")

scn <- unique(CKutils::read_parquet_dt(
  file.path(out, "summaries", "qalys_scaled_up"))$scenario)
prs <- grep("_prs$", scn, value = TRUE)
prs <- prs[sub("_prs$", "_qrisk", prs) %chin% scn]
cmp <- c("sc0", setNames(sub("_prs$", "_qrisk", prs), prs))

cat("scenarios :", paste(sort(scn), collapse = ", "), "\n")
cat("comparators:\n")
for (s in names(cmp)[nzchar(names(cmp))]) cat("   ", s, "vs", cmp[[s]], "\n")
cat("    everything else vs", unname(cmp[!nzchar(names(cmp))]), "\n\n")

# Keep the published tables intact: redirect to a sibling directory.
tdir <- file.path(out, "tables")
keep <- file.path(out, "tables_published_backup")
if (dir.exists(keep)) unlink(keep, recursive = TRUE)
file.rename(tdir, keep)

IMPACTncd$export_tables(
  multicore = TRUE,
  baseline_year_for_change_outputs = GV_PROGRAMME_START,
  comparator_scenario = cmp,
  custom_costs_in_healthcare = GV_COST_COLS
)

new_dir <- tdir
cat("\n=== comparators.csv manifest ===\n")
print(fread(file.path(new_dir, "comparators.csv")))

cat("\n=== CEA files written ===\n")
f <- list.files(new_dir, pattern = "^cost-effectiveness by year \\(")
print(sort(f))

rd <- function(p) fread(file.path(new_dir, paste0(
  "cost-effectiveness by year (", p, "-EQ5D5L) (not standardised).csv")))

cat("\n=== final-year incrementals, by perspective and contrast ===\n")
res <- rbindlist(lapply(c("healthcare", "healthcare_socialcare", "societal"),
  function(p) {
    t <- rd(p)
    t <- t[year == max(year) & type %chin% c("dQALYs_cuml", "dCosts_cuml", "ICER")]
    t[, .(perspective = p, scenario, comparator, type,
          median = signif(`value_50.0%`, 5))]
  }))
print(dcast(res, perspective + scenario + comparator ~ type,
            value.var = "median"))

cat("\n=== NMB column names (pin the WTP default) ===\n")
print(grep("^NMB", unique(rd("healthcare")$type), value = TRUE))

# Restore the published tables; keep the new run alongside for inspection.
probe <- file.path(out, "tables_via_package")
if (dir.exists(probe)) unlink(probe, recursive = TRUE)
file.rename(tdir, probe)
file.rename(keep, tdir)
cat("\npublished tables restored to ", tdir,
    "\nthis run left in            ", probe, "\n", sep = "")
