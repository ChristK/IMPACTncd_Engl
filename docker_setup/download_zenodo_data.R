#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# download_zenodo_data.R
#
# Downloads the IMPACTncd England model data (input data + pre-computed PARFs
# and compiled RR tables) from Zenodo into ./inputs and ./simulation.
#
# This runs in the Docker DATA image (Dockerfile.data.IMPACTncdENGL), which sits
# BENEATH the model-package layer, so it must NOT depend on the installed
# IMPACTncdEngland package. It therefore sources the ZenodoAssetManager class
# directly (its only needs — R6, data.table, zen4R, httr2, zip, digest — are in
# the prerequisite image). The published data is PUBLIC, so no Zenodo account or
# token is required.
#
# Run from the project root. Environment variables:
#   DOWNLOAD_DATA        "false"/"0"/"no" to skip (build a data-less image).
#   ZENODO_CONCEPT_DOI   Override the Zenodo concept DOI (default: the published
#                        IMPACTncd England input-data record).
# -----------------------------------------------------------------------------

if (tolower(Sys.getenv("DOWNLOAD_DATA", "true")) %in% c("false", "0", "no", "off")) {
  message("DOWNLOAD_DATA disabled - skipping Zenodo data download.")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages(library(data.table))
# Source the asset manager directly (no package install needed at this layer).
source("Rpackage/IMPACTncd_England_model_pkg/R/ZenodoManager_class.R")

concept_doi <- Sys.getenv("ZENODO_CONCEPT_DOI", unset = "10.5281/zenodo.20812409")
message("Downloading IMPACTncd England data from Zenodo (", concept_doi, ") ...")

zam <- ZenodoAssetManager$new(logs = TRUE)
zam$connect(token = NULL)              # anonymous; published public data
zam$set_concept_doi(concept_doi)

dir.create("inputs", showWarnings = FALSE)
dir.create("simulation", showWarnings = FALSE)

# 1) Input archives (everything except simulation parf/rr) -> ./inputs
zam$sync_inputs(
  input_base = "inputs", action = "download", overwrite = FALSE,
  exclude_archive_patterns = "^(parf|rr)"
)

# 2) Simulation archives (parf, rr) -> ./simulation (mirrors
#    Simulation$zenodo_download_PARFs_RRs(), which we cannot call here).
remote <- zam$list_remote_files()
sim <- remote[grepl("^(parf|rr)", filename)]
if (nrow(sim) > 0L) {
  ddir <- file.path(tempdir(), "sim_dl")
  dir.create(ddir, showWarnings = FALSE)
  for (i in seq_len(nrow(sim))) {
    f <- sim$filename[i]
    zam$download_file(f, ddir, overwrite = TRUE, checksum = sim$checksum[i])
    zam$extract_archive(file.path(ddir, f), "simulation", overwrite = FALSE)
  }
  unlink(ddir, recursive = TRUE)
}

message("Zenodo data download complete.")
