# Test script for ZenodoAssetManager
# Customize this script before running

# Load the package (or source the file directly for testing)
# library(IMPACTncdEngland)
source("./Rpackage/IMPACTncd_England_model_pkg/R/ZenodoManager_class.R")

# ============================================================================
# KEYRING SETUP (run once to store your token securely)
# ============================================================================
#
# This script uses a FILE-BASED keyring backend, which works on servers without
# a desktop environment (where GNOME Keyring/KWallet aren't available).
#
# The keyring file is stored at: ~/.config/r-keyring/
# It is encrypted with a master password you choose.
#
# ============================================================================
# FIRST-TIME SETUP (run these commands interactively in R console):
# ============================================================================
#
# # 1. Install keyring if needed
# install.packages("keyring")
#
# # 2. Create the file backend and keyring
# kb <- keyring::backend_file$new()
# kb$keyring_create("zenodo_keyring")
# # Enter a master password when prompted (characters won't show - that's normal)
#
# # 3. Store your Zenodo sandbox token
# kb$set("zenodo_sandbox", username = "token", keyring = "zenodo_keyring")
# # Enter your Zenodo sandbox token when prompted
#
# # 4. Store your Zenodo production token (when ready)
# kb$set("zenodo", username = "token", keyring = "zenodo_keyring")
# # Enter your Zenodo production token when prompted
#
# ============================================================================
# Get your tokens from:
#   Sandbox:    https://sandbox.zenodo.org/account/settings/applications/tokens/new/
#   Production: https://zenodo.org/account/settings/applications/tokens/new/
#   Required scopes: deposit:write, deposit:actions
#
# To list stored keys:
#   kb <- keyring::backend_file$new()
#   kb$list(keyring = "zenodo_keyring")
#
# To delete a stored token:
#   kb <- keyring::backend_file$new()
#   kb$delete("zenodo_sandbox", username = "token", keyring = "zenodo_keyring")
#
# To list available keyrings:
#   kb <- keyring::backend_file$new()
#   kb$keyring_list()
#
# ============================================================================

# ============================================================================
# CONFIGURATION - CUSTOMIZE THESE
# ============================================================================

# Use sandbox (TRUE) or production (FALSE)
USE_SANDBOX <- TRUE

# Keyring name (file-based)
KEYRING_NAME <- "zenodo_keyring"

# Helper function to retrieve token from file-based keyring
get_zenodo_token <- function(sandbox = TRUE, keyring_name = "zenodo_keyring") {
  if (!requireNamespace("keyring", quietly = TRUE)) {
    stop("Package 'keyring' required. Install with: install.packages('keyring')")
  }

  service <- if (sandbox) "zenodo_sandbox" else "zenodo"

  # Use file backend explicitly

  kb <- keyring::backend_file$new()

  # Check if keyring exists
  available_keyrings <- kb$keyring_list()$keyring
  if (!(keyring_name %in% available_keyrings)) {
    stop(
      "Keyring '", keyring_name, "' not found.\n",
      "Create it with:\n",
      "  kb <- keyring::backend_file$new()\n",
      "  kb$keyring_create('", keyring_name, "')\n",
      "Then store your token with:\n",
      "  kb$set('", service, "', username = 'token', keyring = '", keyring_name, "')"
    )
  }

  # Get the token (will prompt for keyring password if locked)
  tryCatch(
    kb$get(service, username = "token", keyring = keyring_name),
    error = function(e) {
      stop(
        "Token not found in keyring for '", service, "'.\n",
        "Store it with:\n",
        "  kb <- keyring::backend_file$new()\n",
        "  kb$set('", service, "', username = 'token', keyring = '", keyring_name, "')\n",
        "Get your token from: https://",
        if (sandbox) "sandbox." else "", "zenodo.org/account/settings/applications/tokens/new/"
      )
    }
  )
}

# --- Inputs configuration ---
INPUT_BASE <- "./inputs"

# Directories to include (NULL = all)
# Your available directories:
#   aux_data, disease_burden, exposure_distributions, health_econ,
#   mortality, pop_estimates_lsoa, pop_projections, RR, synthpop
INPUT_DIRECTORIES <- c(
  "aux_data",
  "disease_burden",
  "exposure_distributions",
  "health_econ",
  "mortality",
  "pop_estimates_lsoa",
  "pop_projections"
  # Add more as needed
)

# Directories/subfolders to exclude (regex patterns)
INPUT_EXCLUDE_PATTERNS <- c(
  "^unprocessed$", # Exclude 'unprocessed' directories
  "_backup$", # Exclude anything ending in _backup
  "_old$", # Exclude anything ending in _old
  "scripts$", # Exclude anything ending in scripts
  "validation$" # Exclude anything ending in validation
)

# File patterns to exclude from archives (regex)
INPUT_EXCLUDE_FILE_PATTERNS <- c(
  "\\.R$",   # Exclude R scripts
  "\\.Rmd$", # Exclude R Markdown files
  "\\.yaml$", # Exclude YAML files
  "\\.html$" # Exclude HTML files
)

# Group subfolders by prefix within each directory?
# e.g., disease_burden/cancer_2020, disease_burden/cancer_2021 -> disease_burden_cancer.zip
INPUT_GROUP_BY_PREFIX <- TRUE

# --- Simulation configuration ---
SIMULATION_BASE <- "./simulation"

# Directories to include from simulation
# Available: parf (PAFs per disease), rr (relative risks)
SIMULATION_DIRECTORIES <- c(
  "parf",
  "rr"
)

# File patterns to exclude from simulation archives (regex)
SIMULATION_EXCLUDE_FILE_PATTERNS <- c(
  "\\.R$",   # Exclude R scripts
  "\\.csv$"  # Exclude CSV metadata files
)

# Group parf subfolders by disease prefix?
# e.g., PARF_af_<hash>, PARF_chd_<hash> -> parf_PARF_af.zip, parf_PARF_chd.zip
SIMULATION_GROUP_BY_PREFIX <- TRUE

# Archive storage location (NULL = tempdir)
ARCHIVE_DIR <- NULL  # Uses tempdir(); archives are deleted after successful upload

# ============================================================================
# CREATE MANAGER AND CONNECT
# ============================================================================

manager <- ZenodoAssetManager$new(
  hash_file = "./simulation/zenodo_manifest.csv",
  archive_dir = ARCHIVE_DIR,
  logs = TRUE,
  sandbox = USE_SANDBOX
)

# Enable progress callbacks
manager$set_progress_callback(upload = TRUE, download = TRUE)

# Connect to Zenodo using token from file-based keyring
ZENODO_TOKEN <- get_zenodo_token(sandbox = USE_SANDBOX, keyring_name = KEYRING_NAME)
manager$connect(token = ZENODO_TOKEN)

# Print manager status
print(manager)

# ============================================================================
# STEP 1: CREATE ARCHIVES (without uploading)
# ============================================================================

message("=== Creating input archives ===")

input_archives <- manager$create_input_archives(
  input_base = INPUT_BASE,
  directories = INPUT_DIRECTORIES,
  exclude_patterns = INPUT_EXCLUDE_PATTERNS,
  exclude_file_patterns = INPUT_EXCLUDE_FILE_PATTERNS,
  group_by_prefix = INPUT_GROUP_BY_PREFIX,
  compression_level = 6L,
  multicore = TRUE,
  n_cores = 20L, # TODO: use sim_prm$clusternumber
  update_gitignore = TRUE
)
print(input_archives)

message("\n=== Creating simulation archives ===")

simulation_archives <- manager$create_input_archives(
  input_base = SIMULATION_BASE,
  directories = SIMULATION_DIRECTORIES,
  exclude_file_patterns = SIMULATION_EXCLUDE_FILE_PATTERNS,
  group_by_prefix = SIMULATION_GROUP_BY_PREFIX,
  compression_level = 6L,
  multicore = TRUE,
  n_cores = 20L,
  update_gitignore = TRUE
)
print(simulation_archives)

# Combine all archives
archives <- rbind(input_archives, simulation_archives)

# Check archive sizes
message("\nInput archives size: ",
        round(sum(input_archives$size_bytes) / 1024^3, 2), " GB")
message("Simulation archives size: ",
        round(sum(simulation_archives$size_bytes) / 1024^3, 2), " GB")
message("Total archive size: ",
        round(sum(archives$size_bytes) / 1024^3, 2), " GB")

# View where archives are stored
message("Archives stored in: ", manager$archive_dir)

# ============================================================================
# STEP 2: CREATE ZENODO RECORD (first time only)
# ============================================================================

# Uncomment this section to create a new record

manager$create_new_record(
  title = "IMPACTncd England Model Input Data",
  description = paste(
    "Input and simulation data files for the IMPACTncd England microsimulation model.",
    "This dataset includes mortality rates, disease burden estimates,",
    "population projections, other required inputs,",
    "population attributable fractions (PAFs), and relative risks (RRs)."
  ),
  creators = list(
    list(
      firstname = "Chris",
      lastname = "Kypridemos",
      orcid = "0000-0002-0746-9229"  # Optional
    )
    # Add more creators as needed
  ),
  version = "1.0.0",
  keywords = c("microsimulation", "health", "IMPACTncd", "England", "NCD"),
  license = "cc-by-sa-4.0",
  publisher = "Zenodo"
)

# Save the concept DOI for future use!
message("\nIMPORTANT: Save this concept DOI for future uploads:")
message("  Concept DOI: ", manager$concept_doi)
# Concept DOI: 10.5072/zenodo.442996

# ============================================================================
# STEP 3: UPLOAD ARCHIVES
# ============================================================================

# Uncomment to upload (requires a record to be created or loaded first)
for (i in seq_len(nrow(archives))) {
  manager$upload_archive(archives$archive_path[i])
}
#
# List uploaded files
manager$list_remote_files()

# ============================================================================
# STEP 4: PUBLISH (IRREVERSIBLE!)
# ============================================================================

# Only uncomment when ready - this mints a permanent DOI!
# manager$publish_record()

# ============================================================================
# FOR SUBSEQUENT UPLOADS (new versions)
# ============================================================================

# If you already have a published record, use this workflow:
#
# # Set the concept DOI from your previous upload
# manager$set_concept_doi("10.5072/zenodo.442996")
#
# # Create a new version
# manager$create_new_version(version = "1.0.1", delete_previous_files = TRUE)
#
# # Upload new archives
# for (i in seq_len(nrow(archives))) {
#   manager$upload_archive(archives$archive_path[i])
# }
#
# # Publish the new version
# manager$publish_record()
