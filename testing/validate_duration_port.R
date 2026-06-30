# Functional validation of the forward-duration port: run a small baseline
# simulation and inspect the asthma duration distribution. If the fix works,
# asthma spells last a stochastic number of years (prvl takes values 2,3,4,...),
# not the 1-year truncation the broken engine produced.
source("./global.R")

IMPACTncd <- Simulation$new("testing/sim_design_testing.yaml")
IMPACTncd$del_logs()$del_outputs()$run(1L, multicore = FALSE, "sc0")

suppressMessages({library(arrow); library(data.table)})

fs <- list.files("outputs", pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
cat("=== parquet outputs written:", length(fs), "===\n")
lc <- fs[grepl("lifecourse", fs)]
if (!length(lc)) lc <- fs

# Find a file that actually carries an asthma prevalence/duration column
acol <- NULL; sample_cols <- NULL
for (f in head(lc, 200)) {
  nm <- names(arrow::open_dataset(f)$schema)
  sample_cols <- nm
  hit <- grep("^asthma.*prvl$|^prvl.*asthma|^asthma$", nm, value = TRUE)
  if (length(hit)) { acol <- hit[1]; break }
}

if (is.null(acol)) {
  cat("No asthma prevalence column found. Columns in a sample output:\n")
  print(sample_cols)
} else {
  cat("Asthma duration column:", acol, "\n")
  dt <- rbindlist(lapply(lc, function(f) {
    d <- tryCatch(as.data.table(read_parquet(f, col_select = acol)), error = function(e) NULL)
    d
  }), fill = TRUE)
  setnames(dt, acol, "asthma_prvl")
  cat("=== asthma duration (years) distribution among diseased person-years ===\n")
  print(dt[asthma_prvl > 0, .N, keyby = asthma_prvl][order(asthma_prvl)])
  pos <- dt[asthma_prvl > 0, .N]
  cat("max spell duration (years):", max(dt$asthma_prvl, na.rm = TRUE), "\n")
  cat("diseased person-years:", pos, "\n")
  cat("share with duration > 1 year:", if (pos) round(dt[asthma_prvl > 1, .N] / pos, 4) else NA, "\n")
  cat(if (pos && dt[asthma_prvl > 1, .N] > 0) "RESULT: PASS - asthma spells span multiple years (stochastic duration active)\n"
      else "RESULT: FAIL - asthma still truncated to 1 year\n")
}
cat("VALIDATION_SCRIPT_DONE\n")
