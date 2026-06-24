#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# download_zenodo_data.R
#
# Downloads the IMPACTncd England model data (input data + pre-computed PARFs
# and compiled RR tables) from Zenodo. This is run during the Docker build
# (see Dockerfile.IMPACTncdENGL) so the image does not need ~13GB of data
# bundled into the build context — the canonical, versioned data comes from the
# published Zenodo record instead.
#
# The published data is PUBLIC, so no Zenodo account or token is required.
#
# Run from the project root. Environment variables:
#   DOWNLOAD_DATA        "false"/"0"/"no" to skip the download (build a
#                        code-only image and download later at runtime).
#   ZENODO_CONCEPT_DOI   Override the Zenodo concept DOI (default: the published
#                        IMPACTncd England input-data record).
# -----------------------------------------------------------------------------

if (tolower(Sys.getenv("DOWNLOAD_DATA", "true")) %in% c("false", "0", "no", "off")) {
  message("DOWNLOAD_DATA is disabled - skipping Zenodo data download.")
  quit(save = "no", status = 0)
}

suppressPackageStartupMessages(library(IMPACTncdEngland))

concept_doi <- Sys.getenv("ZENODO_CONCEPT_DOI", unset = "10.5281/zenodo.20812409")
message("Downloading IMPACTncd England data from Zenodo (", concept_doi, ") ...")

# Build a Simulation object (skeleton mode — data is not present yet), connect
# anonymously to the public record, and download everything needed to run the
# model. zenodo_download_all() reconstructs ./inputs and ./simulation/{parf,rr}.
IMPACTncd <- Simulation$new("inputs/sim_design.yaml")
IMPACTncd$zenodo_connect(concept_doi = concept_doi)  # anonymous; no token needed
IMPACTncd$zenodo_download_all()

message("Zenodo data download complete.")
