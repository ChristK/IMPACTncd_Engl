# Test script to compare export_tables() output with process_out_Bradford.R
# This verifies that the refactored Simulation$export_tables() produces
# identical results to the original standalone script.

library(data.table)

source("./global.R")

# Configuration
design <- Design$new("./testing/sim_design_testing.yaml")
output_dir <- design$sim_prm$output_dir
tables_dir <- file.path(output_dir, "tables")
tables_reference_dir <- file.path(output_dir, "tables_reference")

# Clean up any existing directories
if (dir.exists(tables_dir)) {
 unlink(tables_dir, recursive = TRUE)
}
if (dir.exists(tables_reference_dir)) {
 unlink(tables_reference_dir, recursive = TRUE)
}

# ============================================================================
# STEP 1: Generate reference tables using process_out_Bradford.R
# ============================================================================
message("\n=== STEP 1: Running process_out_Bradford.R ===\n")

# Create tables directory (process_out_Bradford.R doesn't create it)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

source("./testing/process_out.R")

# Rename to reference
file.rename(tables_dir, tables_reference_dir)
message("Reference tables saved to: ", tables_reference_dir)

# ============================================================================
# STEP 2: Generate tables using Simulation$export_tables()
# ============================================================================
message("\n=== STEP 2: Running Simulation$export_tables() ===\n")
IMPACTncd <- Simulation$new("./testing/sim_design_testing.yaml")
IMPACTncd$export_tables()
message("New tables saved to: ", tables_dir)

# ============================================================================
# STEP 3: Compare all CSV files
# ============================================================================
message("\n=== STEP 3: Comparing CSV files ===\n")

ref_files <- list.files(tables_reference_dir, pattern = "\\.csv$", full.names = FALSE)
new_files <- list.files(tables_dir, pattern = "\\.csv$", full.names = FALSE)

# Check for missing files
missing_in_new <- setdiff(ref_files, new_files)
extra_in_new <- setdiff(new_files, ref_files)

if (length(missing_in_new) > 0) {
 message("WARNING: Files missing in new output:")
 for (f in missing_in_new) message("
- ", f)
}

if (length(extra_in_new) > 0) {
 message("INFO: Extra files in new output:")
 for (f in extra_in_new) message("
+ ", f)
}

# Compare common files
common_files <- intersect(ref_files, new_files)
message("\nComparing ", length(common_files), " common files...\n")

results <- data.table(
 file = character(),
 identical = logical(),
 max_abs_diff = numeric(),
 max_rel_diff = numeric(),
 notes = character()
)

for (f in common_files) {
 ref_path <- file.path(tables_reference_dir, f)
 new_path <- file.path(tables_dir, f)

 ref_dt <- fread(ref_path)
 new_dt <- fread(new_path)

 # Check dimensions
 if (nrow(ref_dt) != nrow(new_dt) || ncol(ref_dt) != ncol(new_dt)) {
   results <- rbind(results, data.table(
     file = f,
     identical = FALSE,
     max_abs_diff = NA_real_,
     max_rel_diff = NA_real_,
     notes = sprintf("Dimension mismatch: ref(%d x %d) vs new(%d x %d)",
                     nrow(ref_dt), ncol(ref_dt), nrow(new_dt), ncol(new_dt))
   ))
   next
 }

 # Check column names
 if (!identical(names(ref_dt), names(new_dt))) {
   results <- rbind(results, data.table(
     file = f,
     identical = FALSE,
     max_abs_diff = NA_real_,
     max_rel_diff = NA_real_,
     notes = "Column names differ"
   ))
   next
 }

 # Sort both tables by key columns for comparison
 key_cols <- intersect(names(ref_dt), c("scenario", "year", "sex", "agegrp", "disease",
                                         "type", "scale", "costs_type", "exposure", "variable"))
 if (length(key_cols) > 0) {
   setkeyv(ref_dt, key_cols)
   setkeyv(new_dt, key_cols)
 }

 # Check for exact equality first
 if (identical(ref_dt, new_dt)) {
   results <- rbind(results, data.table(
     file = f,
     identical = TRUE,
     max_abs_diff = 0,
     max_rel_diff = 0,
     notes = "Exact match"
   ))
   next
 }

 # Check numeric columns for near-equality
 numeric_cols <- names(ref_dt)[sapply(ref_dt, is.numeric)]
 max_abs <- 0
 max_rel <- 0
 na_mismatch <- FALSE
 na_mismatch_col <- ""

 for (col in numeric_cols) {
   ref_vals <- ref_dt[[col]]
   new_vals <- new_dt[[col]]

   # Skip if both are all NA
   if (all(is.na(ref_vals)) && all(is.na(new_vals))) next

   # Check for NA pattern match
   if (!identical(is.na(ref_vals), is.na(new_vals))) {
     na_mismatch <- TRUE
     na_mismatch_col <- col
     break
   }

   # Calculate differences where both are not NA
   valid_idx <- !is.na(ref_vals) & !is.na(new_vals)
   if (any(valid_idx)) {
     abs_diff <- abs(ref_vals[valid_idx] - new_vals[valid_idx])
     max_abs <- max(max_abs, max(abs_diff))

     # Relative difference (avoid division by zero)
     denom <- pmax(abs(ref_vals[valid_idx]), 1e-10)
     rel_diff <- abs_diff / denom
     max_rel <- max(max_rel, max(rel_diff))
   }
 }

 if (na_mismatch) {
   results <- rbind(results, data.table(
     file = f,
     identical = FALSE,
     max_abs_diff = NA_real_,
     max_rel_diff = NA_real_,
     notes = paste("NA pattern differs in column:", na_mismatch_col)
   ))
   next
 }

 # Check non-numeric columns
 non_numeric_cols <- setdiff(names(ref_dt), numeric_cols)
 char_match <- TRUE
 char_diff_col <- ""
 for (col in non_numeric_cols) {
   if (!identical(ref_dt[[col]], new_dt[[col]])) {
     char_match <- FALSE
     char_diff_col <- col
     break
   }
 }

 if (!char_match) {
   results <- rbind(results, data.table(
     file = f,
     identical = FALSE,
     max_abs_diff = max_abs,
     max_rel_diff = max_rel,
     notes = paste("Non-numeric column differs:", char_diff_col)
   ))
 } else if (max_abs < 1e-10 && max_rel < 1e-10) {
   results <- rbind(results, data.table(
     file = f,
     identical = TRUE,
     max_abs_diff = max_abs,
     max_rel_diff = max_rel,
     notes = "Numerically identical (within tolerance)"
   ))
 } else {
   results <- rbind(results, data.table(
     file = f,
     identical = FALSE,
     max_abs_diff = max_abs,
     max_rel_diff = max_rel,
     notes = "Numeric differences found"
   ))
 }
}

# ============================================================================
# STEP 4: Print results
# ============================================================================
message("\n=== COMPARISON RESULTS ===\n")
print(results)

n_identical <- sum(results$identical)
n_different <- sum(!results$identical)

message("\n=== SUMMARY ===")
message("Reference files: ", length(ref_files))
message("New files: ", length(new_files))
message("Common files compared: ", nrow(results))
message("Identical: ", n_identical)
message("Different: ", n_different)
message("Missing in new: ", length(missing_in_new))
message("Extra in new: ", length(extra_in_new))

if (n_different > 0) {
 message("\nFiles with differences:")
 print(results[identical == FALSE])
}

if (n_identical == nrow(results) && length(missing_in_new) == 0) {
 message("\n*** SUCCESS: All outputs are identical! ***")
} else {
 message("\n*** DIFFERENCES FOUND - Review results above ***")
}
