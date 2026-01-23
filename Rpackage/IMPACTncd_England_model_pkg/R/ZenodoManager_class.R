## IMPACTncdEngland is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
##
## Copyright (C) 2018-2026 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngland is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

#' R6 Class for managing Zenodo data archives
#'
#' @description
#' The `ZenodoAssetManager` class provides a complete pipeline for managing
#' input data assets via Zenodo. It supports uploading, downloading, hashing,
#' and versioning of parquet and other data files.
#'
#' @details
#' This class wraps the zen4R package to provide:
#' \itemize{
#'   \item **Archive Creation**: Groups files into zip/tar archives by directory
#'   \item **Hash Verification**: Uses digest to detect file changes
#'   \item **Version Control**: Manages DOI-based versioning on Zenodo
#'   \item **Cross-Platform**: Works on Linux, Windows, and macOS
#' }
#'
#' @section Zenodo Authentication:
#' Requires a Zenodo personal access token. Set via:
#' \itemize{
#'   \item Environment variable: `ZENODO_TOKEN`
#'   \item Or pass directly to `connect()` method
#' }
#'
#' Get your token at: https://zenodo.org/account/settings/applications/tokens/new/
#'
#' @export
ZenodoAssetManager <- R6::R6Class(
  classname = "ZenodoAssetManager",
  lock_objects = TRUE,
  lock_class = TRUE,

  # Public methods ----
  public = list(
    #' @field zenodo A ZenodoManager object from zen4R package.
    zenodo = NULL,

    #' @field concept_doi The concept DOI for the record (shared across versions).
    concept_doi = NULL,

    #' @field record The current Zenodo record object.
    record = NULL,

    #' @field hash_file Path to the hash manifest file.
    hash_file = NULL,

    #' @field archive_dir Directory for storing temporary archives.
    archive_dir = NULL,

    #' @field logs Enable verbose logging.
    logs = TRUE,

    #' @field sandbox Use Zenodo sandbox for testing.
    sandbox = FALSE,

    # initialize ----
    #' @description Create a new ZenodoAssetManager object.
    #' @param hash_file Path to store/read hash manifest (default: "./simulation/zenodo_manifest.csv").
    #' @param archive_dir Temporary directory for archives (default: tempdir()).
    #' @param logs Enable verbose logging.
    #' @param sandbox Use Zenodo sandbox (for testing).
    #' @return A new `ZenodoAssetManager` object.
    initialize = function(
      hash_file = "./simulation/zenodo_manifest.csv",
      archive_dir = NULL,
      logs = TRUE,
      sandbox = FALSE
    ) {
      self$hash_file <- hash_file
      self$archive_dir <- archive_dir %||% file.path(tempdir(), "zenodo_archives")
      self$logs <- logs
      self$sandbox <- sandbox

      # Ensure archive directory exists
      if (!dir.exists(self$archive_dir)) {
        dir.create(self$archive_dir, recursive = TRUE)
      }

      invisible(self)
    },

    # connect ----
    #' @description Connect to Zenodo API.
    #' @param token Zenodo personal access token. If NULL, reads from ZENODO_TOKEN env var.
    #' @param sandbox Use Zenodo sandbox for testing.
    #' @return The invisible self for chaining.
    connect = function(token = NULL, sandbox = NULL) {
      if (!requireNamespace("zen4R", quietly = TRUE)) {
        stop(
          "Package 'zen4R' is required for Zenodo integration. ",
          "Install it with: install.packages('zen4R')"
        )
      }

      if (is.null(token)) {
        token <- Sys.getenv("ZENODO_TOKEN", unset = NA)
        if (is.na(token) || token == "") {
          stop(
            "Zenodo token not provided. Set ZENODO_TOKEN environment variable ",
            "or pass token to connect(). Get your token at: ",
            "https://zenodo.org/account/settings/applications/tokens/new/"
          )
        }
      }

      if (!is.null(sandbox)) {
        self$sandbox <- sandbox
      }

      logger_level <- if (self$logs) "INFO" else NULL

      self$zenodo <- zen4R::ZenodoManager$new(
        token = token,
        sandbox = self$sandbox,
        logger = logger_level
      )

      if (self$logs) {
        message(
          "Connected to Zenodo ",
          if (self$sandbox) "(sandbox)" else "(production)"
        )
      }

      invisible(self)
    },

    # set_concept_doi ----
    #' @description Set the concept DOI for the record.
    #' @param concept_doi The concept DOI (shared across all versions of a record).
    #' @return The invisible self for chaining.
    set_concept_doi = function(concept_doi) {
      if (!is.null(concept_doi) && !grepl("^10\\.", concept_doi)) {
        stop("Invalid DOI format. DOI should start with '10.'")
      }
      self$concept_doi <- concept_doi
      invisible(self)
    },

    # get_record ----
    #' @description Get the Zenodo record by concept DOI.
    #' @param concept_doi Optional concept DOI. Uses stored value if not provided.
    #' @return The Zenodo record object.
    get_record = function(concept_doi = NULL) {
      private$check_connection()

      doi <- concept_doi %||% self$concept_doi
      if (is.null(doi)) {
        stop("No concept DOI set. Use set_concept_doi() first or provide concept_doi.")
      }

      self$record <- self$zenodo$getRecordByConceptDOI(doi)
      if (is.null(self$record)) {
        stop("Record not found for DOI: ", doi)
      }

      if (self$logs) {
        message("Retrieved record: ", self$record$metadata$title)
        message("  DOI: ", self$record$pids$doi$identifier)
        message("  Version: ", self$record$metadata$version %||% "N/A")
      }

      self$record
    },

    # list_remote_files ----
    #' @description List files in the current Zenodo record.
    #' @param concept_doi Optional concept DOI.
    #' @return A data.table with file information (name, size, checksum).
    list_remote_files = function(concept_doi = NULL) {
      if (is.null(self$record)) {
        self$get_record(concept_doi)
      }

      files <- self$record$listFiles(pretty = TRUE)
      if (is.null(files) || length(files) == 0) {
        return(data.table::data.table(
          filename = character(),
          size = numeric(),
          checksum = character()
        ))
      }

      data.table::as.data.table(files)
    },

    # create_archive ----
    #' @description Create a compressed zip archive from a directory.
    #' @param input_dir Path to the directory to archive.
    #' @param archive_name Name of the archive file (without extension).
    #' @param include_patterns File patterns to include (default: all files).
    #' @return Path to the created archive.
    create_archive = function(
      input_dir,
      archive_name = NULL,
      include_patterns = NULL
    ) {

      if (!dir.exists(input_dir)) {
        stop("Directory does not exist: ", input_dir)
      }

      # Default archive name from directory
      if (is.null(archive_name)) {
        archive_name <- basename(input_dir)
      }

      # Build archive path (always use zip for cross-platform compatibility)
      archive_path <- file.path(self$archive_dir, paste0(archive_name, ".zip"))

      # Get files to archive
      files <- list.files(
        input_dir,
        recursive = TRUE,
        full.names = TRUE,
        all.files = FALSE
      )

      if (!is.null(include_patterns)) {
        pattern <- paste(include_patterns, collapse = "|")
        files <- files[grepl(pattern, files, ignore.case = TRUE)]
      }

      if (length(files) == 0) {
        stop("No files found in directory: ", input_dir)
      }

      if (self$logs) {
        message("Creating zip archive: ", archive_path)
        message("  Source: ", input_dir)
        message("  Files: ", length(files))
      }

      # Create zip archive (cross-platform)
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(dirname(input_dir))

      # Remove existing archive if exists
      if (file.exists(archive_path)) {
        file.remove(archive_path)
      }

      zip::zip(
        zipfile = archive_path,
        files = basename(input_dir),
        recurse = TRUE,
        compression_level = 0  # No compression for speed; just merge files
      )

      if (!file.exists(archive_path)) {
        stop("Failed to create archive: ", archive_path)
      }

      if (self$logs) {
        size_mb <- round(file.size(archive_path) / 1024^2, 2)
        message("  Archive size: ", size_mb, " MB")
      }

      archive_path
    },

    # create_input_archives ----
    #' @description Create zip archives for all input directories.
    #' @param input_base Base directory containing inputs (default: "./inputs").
    #' @param directories Specific subdirectories to archive. If NULL, archives all.
    #' @return A data.table with archive information.
    create_input_archives = function(
      input_base = "./inputs",
      directories = NULL
    ) {

      if (!dir.exists(input_base)) {
        stop("Input directory does not exist: ", input_base)
      }

      # Get subdirectories
      if (is.null(directories)) {
        directories <- list.dirs(input_base, full.names = FALSE, recursive = FALSE)
        directories <- directories[directories != ""]
      }

      if (length(directories) == 0) {
        stop("No directories found in: ", input_base)
      }

      if (self$logs) {
        message("Creating archives for ", length(directories), " directories:")
        message("  ", paste(directories, collapse = ", "))
      }

      # Create archives
      results <- lapply(directories, function(dir_name) {
        dir_path <- file.path(input_base, dir_name)
        if (!dir.exists(dir_path)) {
          warning("Directory not found, skipping: ", dir_path)
          return(NULL)
        }

        archive_path <- tryCatch(
          self$create_archive(dir_path),
          error = function(e) {
            warning("Failed to archive ", dir_name, ": ", e$message)
            NULL
          }
        )

        if (!is.null(archive_path)) {
          data.table::data.table(
            directory = dir_name,
            archive_path = archive_path,
            size_bytes = file.size(archive_path),
            hash = private$compute_file_hash(archive_path)
          )
        }
      })

      results <- data.table::rbindlist(results[!sapply(results, is.null)])
      results
    },

    # compute_manifest ----
    #' @description Compute hash manifest for a directory.
    #' @param input_dir Directory to scan.
    #' @param recursive Scan recursively.
    #' @param file_patterns File patterns to include.
    #' @return A data.table with file hashes.
    compute_manifest = function(
      input_dir,
      recursive = TRUE,
      file_patterns = c("\\.parquet$", "\\.fst$", "\\.csv$", "\\.csvy$", "\\.yaml$")
    ) {
      if (!dir.exists(input_dir)) {
        stop("Directory does not exist: ", input_dir)
      }

      # Get all files
      all_files <- list.files(
        input_dir,
        recursive = recursive,
        full.names = TRUE,
        all.files = FALSE
      )

      # Filter by patterns
      if (!is.null(file_patterns)) {
        pattern <- paste(file_patterns, collapse = "|")
        all_files <- all_files[grepl(pattern, all_files, ignore.case = TRUE)]
      }

      if (length(all_files) == 0) {
        return(data.table::data.table(
          path = character(),
          relative_path = character(),
          hash = character(),
          size_bytes = numeric(),
          mtime = as.POSIXct(character())
        ))
      }

      if (self$logs) {
        message("Computing hashes for ", length(all_files), " files...")
      }

      # Compute hashes (with progress for large numbers)
      results <- lapply(seq_along(all_files), function(i) {
        f <- all_files[i]
        if (self$logs && i %% 100 == 0) {
          message("  Progress: ", i, "/", length(all_files))
        }

        data.table::data.table(
          path = f,
          relative_path = gsub(paste0("^", normalizePath(input_dir), "/?"), "", f),
          hash = private$compute_file_hash(f),
          size_bytes = file.size(f),
          mtime = file.mtime(f)
        )
      })

      data.table::rbindlist(results)
    },

    # save_manifest ----
    #' @description Save manifest to file.
    #' @param manifest A data.table manifest.
    #' @param path Path to save (default: self$hash_file).
    #' @return The invisible self for chaining.
    save_manifest = function(manifest, path = NULL) {
      path <- path %||% self$hash_file

      # Ensure directory exists
      dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

      data.table::fwrite(manifest, path)

      if (self$logs) {
        message("Manifest saved to: ", path)
      }

      invisible(self)
    },

    # load_manifest ----
    #' @description Load manifest from file.
    #' @param path Path to load (default: self$hash_file).
    #' @return A data.table manifest or NULL if file doesn't exist.
    load_manifest = function(path = NULL) {
      path <- path %||% self$hash_file

      if (!file.exists(path)) {
        if (self$logs) {
          message("Manifest file not found: ", path)
        }
        return(NULL)
      }

      manifest <- data.table::fread(path)

      if (self$logs) {
        message("Loaded manifest with ", nrow(manifest), " entries from: ", path)
      }

      manifest
    },

    # compare_manifests ----
    #' @description Compare two manifests to find differences.
    #' @param local_manifest Local file manifest.
    #' @param remote_manifest Remote file manifest.
    #' @return A list with added, removed, and modified files.
    compare_manifests = function(local_manifest, remote_manifest) {
      if (is.null(local_manifest) || nrow(local_manifest) == 0) {
        return(list(
          added = character(),
          removed = if (!is.null(remote_manifest)) remote_manifest$relative_path else character(),
          modified = character(),
          unchanged = character()
        ))
      }

      if (is.null(remote_manifest) || nrow(remote_manifest) == 0) {
        return(list(
          added = local_manifest$relative_path,
          removed = character(),
          modified = character(),
          unchanged = character()
        ))
      }

      # Merge manifests
      merged <- merge(
        local_manifest[, .(relative_path, local_hash = hash)],
        remote_manifest[, .(relative_path, remote_hash = hash)],
        by = "relative_path",
        all = TRUE
      )

      list(
        added = merged[is.na(remote_hash), relative_path],
        removed = merged[is.na(local_hash), relative_path],
        modified = merged[
          !is.na(local_hash) & !is.na(remote_hash) & local_hash != remote_hash,
          relative_path
        ],
        unchanged = merged[
          !is.na(local_hash) & !is.na(remote_hash) & local_hash == remote_hash,
          relative_path
        ]
      )
    },

    # upload_archive ----
    #' @description Upload an archive to Zenodo.
    #' @param archive_path Path to the archive file.
    #' @param record Zenodo record object. If NULL, uses self$record.
    #' @return The invisible self for chaining.
    upload_archive = function(archive_path, record = NULL) {
      private$check_connection()

      if (!file.exists(archive_path)) {
        stop("Archive file not found: ", archive_path)
      }

      rec <- record %||% self$record
      if (is.null(rec)) {
        stop("No record available. Call get_record() or create_new_version() first.")
      }

      if (self$logs) {
        message("Uploading: ", basename(archive_path))
        size_mb <- round(file.size(archive_path) / 1024^2, 2)
        message("  Size: ", size_mb, " MB")
      }

      self$zenodo$uploadFile(archive_path, record = rec)

      if (self$logs) {
        message("  Upload complete.")
      }

      invisible(self)
    },

    # download_file ----
    #' @description Download a file from Zenodo record.
    #' @param filename Name of the file to download.
    #' @param dest_dir Destination directory.
    #' @param overwrite Overwrite existing files.
    #' @return Path to downloaded file.
    download_file = function(filename, dest_dir, overwrite = FALSE) {
      private$check_connection()

      if (is.null(self$record)) {
        stop("No record loaded. Call get_record() first.")
      }

      dest_path <- file.path(dest_dir, filename)

      if (file.exists(dest_path) && !overwrite) {
        if (self$logs) {
          message("File exists, skipping: ", dest_path)
        }
        return(dest_path)
      }

      # Create destination directory
      dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

      if (self$logs) {
        message("Downloading: ", filename)
      }

      self$record$downloadFiles(path = dest_dir)

      if (!file.exists(dest_path)) {
        stop("Download failed for: ", filename)
      }

      if (self$logs) {
        message("  Downloaded to: ", dest_path)
      }

      dest_path
    },

    # extract_archive ----
    #' @description Extract a zip archive to a directory.
    #' @param archive_path Path to the zip archive.
    #' @param dest_dir Destination directory.
    #' @param overwrite Overwrite existing files.
    #' @return Path to extracted directory.
    extract_archive = function(archive_path, dest_dir, overwrite = FALSE) {
      if (!file.exists(archive_path)) {
        stop("Archive not found: ", archive_path)
      }

      if (!grepl("\\.zip$", archive_path, ignore.case = TRUE)) {
        stop("Only zip archives are supported: ", archive_path)
      }

      dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

      if (self$logs) {
        message("Extracting: ", basename(archive_path))
        message("  To: ", dest_dir)
      }

      zip::unzip(archive_path, exdir = dest_dir, overwrite = overwrite)

      if (self$logs) {
        message("  Extraction complete.")
      }

      dest_dir
    },

    # create_new_record ----
    #' @description Create a new Zenodo record (deposit).
    #' @param title Record title.
    #' @param description Record description.
    #' @param creators List of creators (each with firstname, lastname, optional orcid).
    #' @param version Version string.
    #' @param keywords Vector of keywords.
    #' @param license License identifier (default: "cc-by-4.0").
    #' @return The created record object.
    create_new_record = function(
      title,
      description,
      creators,
      version = "1.0.0",
      keywords = c("microsimulation", "health", "IMPACTncd"),
      license = "cc-by-4.0"
    ) {
      private$check_connection()

      if (self$logs) {
        message("Creating new Zenodo record: ", title)
      }

      # Create empty record
      self$record <- zen4R::ZenodoRecord$new()

      # Set metadata
      self$record$setTitle(title)
      self$record$setDescription(description)
      self$record$setResourceType("dataset")
      self$record$setPublicationDate(Sys.Date())
      self$record$setVersion(version)
      self$record$setLicense(license, sandbox = self$sandbox)
      self$record$setKeywords(keywords)

      # Add creators
      for (creator in creators) {
        self$record$addCreator(
          firstname = creator$firstname,
          lastname = creator$lastname,
          orcid = creator$orcid %||% NULL
        )
      }

      # Deposit the record
      self$record <- self$zenodo$depositRecord(self$record, publish = FALSE)

      if (self$logs) {
        message("Record created with ID: ", self$record$id)
        message("  Draft DOI: ", self$record$pids$doi$identifier %||% "pending")
      }

      # Store concept DOI if available
      if (!is.null(self$record$conceptdoi)) {
        self$concept_doi <- self$record$conceptdoi
      }

      self$record
    },

    # create_new_version ----
    #' @description Create a new version of an existing record.
    #' @param version New version string.
    #' @param delete_previous_files Delete files from previous version.
    #' @return The new version record object.
    create_new_version = function(version, delete_previous_files = TRUE) {
      private$check_connection()

      if (is.null(self$record) && is.null(self$concept_doi)) {
        stop("No record loaded. Call get_record() or provide concept_doi.")
      }

      # Get latest record if not loaded
      if (is.null(self$record)) {
        self$get_record()
      }

      if (self$logs) {
        message("Creating new version: ", version)
        message("  From record ID: ", self$record$id)
      }

      # Update version metadata
      self$record$setVersion(version)

      # Create new version
      self$record <- self$zenodo$depositRecordVersion(
        self$record,
        delete_latest_files = delete_previous_files,
        publish = FALSE
      )

      if (self$logs) {
        message("New version created with ID: ", self$record$id)
      }

      self$record
    },

    # publish_record ----
    #' @description Publish the current record.
    #' @return The published record object.
    publish_record = function() {
      private$check_connection()

      if (is.null(self$record)) {
        stop("No record to publish. Create or load a record first.")
      }

      if (self$logs) {
        message("Publishing record ID: ", self$record$id)
        message("  WARNING: This action is irreversible!")
      }

      self$record <- self$zenodo$publishRecord(self$record$id)

      if (self$logs) {
        message("Record published!")
        message("  DOI: ", self$record$pids$doi$identifier)
      }

      # Update concept DOI
      if (!is.null(self$record$conceptdoi)) {
        self$concept_doi <- self$record$conceptdoi
      }

      self$record
    },

    # get_versions ----
    #' @description Get all versions of the current record.
    #' @return A data.frame with version information.
    get_versions = function() {
      if (is.null(self$record) && is.null(self$concept_doi)) {
        stop("No record loaded. Call get_record() first.")
      }

      if (is.null(self$record)) {
        self$get_record()
      }

      self$record$getVersions()
    },

    # sync_inputs ----
    #' @description Synchronize local inputs with Zenodo.
    #' @param input_base Base input directory.
    #' @param directories Directories to sync (NULL for all).
    #' @param action Action to perform: "check", "download", "upload".
    #' @param overwrite Overwrite local files when downloading.
    #' @param version Version for new uploads.
    #' @return Sync status information.
    sync_inputs = function(
      input_base = "./inputs",
      directories = NULL,
      action = c("check", "download", "upload"),
      overwrite = FALSE,
      version = NULL
    ) {
      action <- match.arg(action)

      if (action == "upload") {
        return(private$sync_upload(input_base, directories, version))
      }

      # For check and download, we need a record
      if (is.null(self$record) && !is.null(self$concept_doi)) {
        self$get_record()
      }

      if (is.null(self$record)) {
        stop("No record available. Set concept_doi and call get_record() first.")
      }

      # Get remote file list
      remote_files <- self$list_remote_files()

      if (nrow(remote_files) == 0) {
        message("No files found in remote record.")
        return(invisible(NULL))
      }

      if (action == "check") {
        return(private$sync_check(remote_files, input_base, directories))
      }

      if (action == "download") {
        return(private$sync_download(remote_files, input_base, directories, overwrite))
      }
    },

    # print ----
    #' @description Print summary of the manager.
    print = function() {
      cat("ZenodoAssetManager\n")
      cat("  Sandbox: ", self$sandbox, "\n")
      cat("  Connected: ", !is.null(self$zenodo), "\n")
      cat("  Concept DOI: ", self$concept_doi %||% "not set", "\n")
      cat("  Record loaded: ", !is.null(self$record), "\n")
      if (!is.null(self$record)) {
        cat("    Record ID: ", self$record$id, "\n")
        cat("    Title: ", self$record$metadata$title %||% "N/A", "\n")
      }
      cat("  Hash file: ", self$hash_file, "\n")
      cat("  Archive dir: ", self$archive_dir, "\n")
      invisible(self)
    }
  ),

  # Private methods ----
  private = list(
    # check_connection ----
    check_connection = function() {
      if (is.null(self$zenodo)) {
        stop("Not connected to Zenodo. Call connect() first.")
      }
    },

    # compute_file_hash ----
    compute_file_hash = function(file_path) {
      if (!file.exists(file_path)) {
        return(NA_character_)
      }
      digest::digest(file = file_path, algo = "xxhash64")
    },

    # sync_check ----
    sync_check = function(remote_files, input_base, directories) {
      if (self$logs) {
        message("\n=== Checking sync status ===")
      }

      # Identify which archives correspond to which directories
      results <- lapply(seq_len(nrow(remote_files)), function(i) {
        fname <- remote_files$filename[i]
        archive_name <- gsub("\\.(zip|tar\\.gz|tgz)$", "", fname, ignore.case = TRUE)

        # Check if corresponding local directory exists
        local_dir <- file.path(input_base, archive_name)
        local_exists <- dir.exists(local_dir)

        data.table::data.table(
          archive = fname,
          directory = archive_name,
          local_exists = local_exists,
          remote_checksum = remote_files$checksum[i] %||% NA_character_
        )
      })

      status <- data.table::rbindlist(results)

      if (self$logs) {
        message("Remote archives: ", nrow(status))
        message("  Local directories found: ", sum(status$local_exists))
        message("  Missing locally: ", sum(!status$local_exists))

        if (sum(!status$local_exists) > 0) {
          message("\nMissing directories:")
          for (dir in status[local_exists == FALSE, directory]) {
            message("  - ", dir)
          }
        }
      }

      status
    },

    # sync_download ----
    sync_download = function(remote_files, input_base, directories, overwrite) {
      if (self$logs) {
        message("\n=== Downloading from Zenodo ===")
      }

      # Create temp directory for downloads
      download_dir <- file.path(self$archive_dir, "downloads")
      dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

      # Download all files from record
      self$record$downloadFiles(path = download_dir)

      # Extract each archive
      for (i in seq_len(nrow(remote_files))) {
        fname <- remote_files$filename[i]
        archive_name <- gsub("\\.(zip|tar\\.gz|tgz)$", "", fname, ignore.case = TRUE)

        # Check if this directory should be processed
        if (!is.null(directories) && !(archive_name %in% directories)) {
          if (self$logs) message("Skipping (not in list): ", archive_name)
          next
        }

        archive_path <- file.path(download_dir, fname)
        if (!file.exists(archive_path)) {
          warning("Archive not found after download: ", fname)
          next
        }

        local_dir <- file.path(input_base, archive_name)

        if (dir.exists(local_dir) && !overwrite) {
          if (self$logs) {
            message("Directory exists, skipping (use overwrite=TRUE to replace): ", archive_name)
          }
          next
        }

        # Extract
        self$extract_archive(archive_path, input_base, overwrite = overwrite)
      }

      # Cleanup downloads
      unlink(download_dir, recursive = TRUE)

      if (self$logs) {
        message("Download complete.")
      }

      invisible(self)
    },

    # sync_upload ----
    sync_upload = function(input_base, directories, version) {
      if (self$logs) {
        message("\n=== Uploading to Zenodo ===")
      }

      # Create zip archives
      archives <- self$create_input_archives(
        input_base = input_base,
        directories = directories
      )

      if (nrow(archives) == 0) {
        message("No archives to upload.")
        return(invisible(NULL))
      }

      # If we have an existing record, create new version
      if (!is.null(self$concept_doi) && !is.null(version)) {
        if (is.null(self$record)) {
          self$get_record()
        }
        self$create_new_version(version, delete_previous_files = TRUE)
      }

      if (is.null(self$record)) {
        stop(
          "No record available for upload. ",
          "Either create_new_record() or set concept_doi and specify version."
        )
      }

      # Upload each archive
      for (i in seq_len(nrow(archives))) {
        self$upload_archive(archives$archive_path[i])
      }

      # Save manifest
      manifest <- data.table::data.table(
        directory = archives$directory,
        hash = archives$hash,
        size_bytes = archives$size_bytes,
        upload_time = Sys.time()
      )
      self$save_manifest(manifest)

      if (self$logs) {
        message("\nUpload complete. ", nrow(archives), " archives uploaded.")
        message("NOTE: Record is still in DRAFT state. Use publish_record() to publish.")
      }

      invisible(archives)
    }
  )
)
