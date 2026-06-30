#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# print_zenodo_version.R
#
# Prints "ZENODO_VERSION=<version>" for the published Zenodo record of a concept
# DOI, so build_push_data.sh can tag the data image with the data version.
#
# Run from the project root (it sources the ZenodoAssetManager class directly,
# like download_zenodo_data.R, so it does not need the installed package).
# Published data is PUBLIC -> anonymous, no token needed.
#
# Environment:
#   ZENODO_CONCEPT_DOI   Override the concept DOI (default: published record).
#
# The ZENODO_VERSION= marker lets the caller extract the version (e.g.
# `grep -oE 'ZENODO_VERSION=[^[:space:]]+' | cut -d= -f2`) even amid zen4R chatter.
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
source("Rpackage/IMPACTncd_England_model_pkg/R/ZenodoManager_class.R")

concept_doi <- Sys.getenv("ZENODO_CONCEPT_DOI", unset = "10.5281/zenodo.20812409")

zam <- ZenodoAssetManager$new(logs = FALSE)
# Suppress zen4R's record tree, which it prints to stdout, so the marker line is
# the only thing on stdout.
invisible(utils::capture.output({
  zam$connect(token = NULL)          # anonymous; published public data
  zam$get_record(concept_doi)
}))

ver <- zam$record$metadata$version
if (is.null(ver) || !nzchar(ver)) {
  stop("Could not determine the Zenodo record version for ", concept_doi)
}

cat(sprintf("ZENODO_VERSION=%s\n", ver))
