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
#'   \item **Archive Creation**: Groups files into zip archives by directory
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
# Standalone worker function for parallel archive creation.
# Defined outside the R6 class so it can be serialized to PSOCK workers
# without dragging along non-serializable objects (e.g., HTTP connections).
.zenodo_archive_worker <- function(task) {
  archive_dir <- task$archive_dir
  logs <- task$logs

  # Helper: compute combined hash of source files
  compute_source_hash <- function(file_paths) {
    if (length(file_paths) == 0) return(NA_character_)
    file_paths <- sort(file_paths)
    hashes <- vapply(file_paths, function(f) {
      if (file.exists(f)) digest::digest(file = f, algo = "xxhash64") else ""
    }, character(1), USE.NAMES = FALSE)
    digest::digest(paste(hashes, collapse = "|"), algo = "xxhash64", serialize = FALSE)
  }

  # Helper: create a single zip archive from one directory or grouped directories
  create_zip <- function(input_paths, archive_name, exclude_file_patterns,
                         compression_level) {
    # Collect files
    all_files <- character(0)
    for (p in input_paths) {
      all_files <- c(all_files, list.files(p, recursive = TRUE, full.names = TRUE,
                                           all.files = FALSE))
    }

    # Exclude files matching patterns
    if (!is.null(exclude_file_patterns) && length(exclude_file_patterns) > 0) {
      exclude_pattern <- paste(exclude_file_patterns, collapse = "|")
      all_files <- all_files[!grepl(exclude_pattern, all_files, ignore.case = TRUE)]
    }

    if (length(all_files) == 0) {
      stop("No files found in: ", paste(input_paths, collapse = ", "))
    }

    source_hash <- compute_source_hash(all_files)

    # Build archive path
    if (is.null(archive_name)) archive_name <- basename(input_paths[1])
    archive_path <- normalizePath(
      file.path(archive_dir, paste0(archive_name, ".zip")),
      mustWork = FALSE
    )

    if (logs) {
      message("Creating zip archive: ", archive_path)
      message("  Files: ", length(all_files))
    }

    # Remove existing archive
    if (file.exists(archive_path)) file.remove(archive_path)

    # Create zip (setwd to common parent)
    parent_dir <- dirname(input_paths[1])
    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(parent_dir)

    zip::zip(
      zipfile = archive_path,
      files = basename(input_paths),
      recurse = TRUE,
      compression_level = compression_level
    )

    if (!file.exists(archive_path)) {
      stop("Failed to create archive: ", archive_path)
    }

    if (logs) {
      size_mb <- round(file.size(archive_path) / 1024^2, 2)
      message("  Archive size: ", size_mb, " MB")
    }

    list(archive_path = archive_path, source_hash = source_hash,
         file_count = length(all_files))
  }

  # Dispatch based on task type
  archive_result <- tryCatch({
    input_paths <- if (task$type == "single") task$dir_path else task$dir_paths
    create_zip(input_paths, task$archive_name, task$exclude_file_patterns,
               task$compression_level)
  }, error = function(e) {
    warning("Failed to archive ", task$dir_name,
            if (!is.na(task$group_prefix)) paste0("/", task$group_prefix) else "",
            ": ", e$message)
    NULL
  })

  if (is.null(archive_result)) return(NULL)

  data.table::data.table(
    directory = task$dir_name,
    group_prefix = task$group_prefix,
    archive_path = archive_result$archive_path,
    source_hash = archive_result$source_hash,
    file_count = archive_result$file_count,
    size_bytes = file.size(archive_result$archive_path)
  )
}

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

    #' @field upload_progress Callback function for upload progress.
    #' Signature: function(bytes_transferred, total_bytes, filename)
    upload_progress = NULL,

    #' @field download_progress Callback function for download progress.
    #' Signature: function(bytes_transferred, total_bytes, filename)
    download_progress = NULL,

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
      self$logs <- logs
      self$sandbox <- sandbox

      # Set archive directory (use tempdir if not provided)
      archive_dir <- archive_dir %||% file.path(tempdir(), "zenodo_archives")

      # Ensure archive directory exists
      if (!dir.exists(archive_dir)) {
        dir.create(archive_dir, recursive = TRUE)
      }

      # Store as absolute path to avoid issues when changing working directory
      self$archive_dir <- normalizePath(archive_dir, mustWork = TRUE)

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

    # set_progress_callback ----
    #' @description Set progress callback functions for uploads and/or downloads.
    #' @param upload Callback for upload progress. Signature: function(bytes, total, filename).
    #'   Set to NULL to disable, or TRUE for a default console progress bar.
    #' @param download Callback for download progress. Signature: function(bytes, total, filename).
    #'   Set to NULL to disable, or TRUE for a default console progress bar.
    #' @return The invisible self for chaining.
    #' @examples
    #' \dontrun{
    #' # Use default console progress bars

    #' manager$set_progress_callback(upload = TRUE, download = TRUE)
    #'
    #' # Use custom callback
    #' manager$set_progress_callback(
    #'   upload = function(bytes, total, filename) {
    #'     pct <- round(100 * bytes / total, 1)
    #'     cat(sprintf("\r%s: %.1f%%", filename, pct))
    #'   }
    #' )
    #'
    #' # Disable progress callbacks
    #' manager$set_progress_callback(upload = NULL, download = NULL)
    #' }
    set_progress_callback = function(upload = NULL, download = NULL) {
      # Handle upload callback
      if (isTRUE(upload)) {
        self$upload_progress <- private$default_progress_callback
      } else if (is.function(upload) || is.null(upload)) {
        self$upload_progress <- upload
      } else {
        stop("upload must be TRUE, NULL, or a function")
      }

      # Handle download callback
      if (isTRUE(download)) {
        self$download_progress <- private$default_progress_callback
      } else if (is.function(download) || is.null(download)) {
        self$download_progress <- download
      } else {
        stop("download must be TRUE, NULL, or a function")
      }

      invisible(self)
    },

    # get_record ----
    #' @description Get the Zenodo record by concept DOI.
    #' Searches published records first, then falls back to draft records.
    #' @param concept_doi Optional concept DOI. Uses stored value if not provided.
    #' @return The Zenodo record object.
    get_record = function(concept_doi = NULL) {
      private$check_connection()

      doi <- concept_doi %||% self$concept_doi
      if (is.null(doi)) {
        stop("No concept DOI set. Use set_concept_doi() first or provide concept_doi.")
      }

      # Try published record first
      self$record <- tryCatch(
        self$zenodo$getRecordByConceptDOI(doi),
        error = function(e) NULL
      )

      if (!is.null(self$record)) {
        if (self$logs) {
          message("Retrieved published record: ", self$record$metadata$title)
          message("  DOI: ", self$record$pids$doi$identifier)
          message("  Version: ", self$record$metadata$version %||% "N/A")
        }
        return(self$record)
      }

      # Fallback: search user's depositions for draft records
      concept_id <- sub("^.*/zenodo\\.", "", doi)
      if (self$logs) {
        message("No published record found. Searching drafts for concept ID: ",
                concept_id)
      }

      deposits <- tryCatch(
        self$zenodo$getDepositions(size = 100L),
        error = function(e) NULL
      )

      if (!is.null(deposits) && length(deposits) > 0L) {
        # Collect ALL deposits matching this concept ID
        matching <- list()
        for (dep in deposits) {
          dep_concept_id <- as.character(
            dep$parent$id %||% dep$conceptrecid %||% ""
          )
          if (dep_concept_id == concept_id) {
            matching[[length(matching) + 1L]] <- dep
          }
        }

        if (length(matching) > 0L) {
          if (self$logs && length(matching) > 1L) {
            ids <- vapply(matching, function(d) as.character(d$id), character(1))
            message("Found ", length(matching),
                    " draft records for concept ID ", concept_id,
                    ": ", paste(ids, collapse = ", "))
          }

          # Prefer the record that has files; fall back to the first match
          self$record <- matching[[1L]]
          for (dep in matching) {
            dep_files <- tryCatch(
              dep$listFiles(pretty = TRUE),
              error = function(e) NULL
            )
            if (!is.null(dep_files) && length(dep_files) > 0L) {
              self$record <- dep
              break
            }
          }

          if (self$logs) {
            message("Retrieved draft record ID: ", self$record$id)
            message("  Title: ", self$record$metadata$title %||% "N/A")
          }
          return(self$record)
        }
      }

      stop(
        "Record not found for DOI: ", doi,
        " (checked both published and draft records)"
      )
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
    #' @description Create a zip archive from a directory.
    #' @param input_dir Path to the directory to archive.
    #' @param archive_name Name of the archive file (without extension).
    #' @param include_patterns File patterns to include (default: all files).
    #' @param exclude_file_patterns Regex patterns to exclude files (e.g., "\\.R$").
    #' @param compression_level Integer 0-9. 0 = no compression (fastest),
    #'   9 = maximum compression (smallest files). Default: 6.
    #' @return A list with archive_path, source_hash, and file_count.
    create_archive = function(
      input_dir,
      archive_name = NULL,
      include_patterns = NULL,
      exclude_file_patterns = NULL,
      compression_level = 6L
    ) {

      if (!dir.exists(input_dir)) {
        stop("Directory does not exist: ", input_dir)
      }

      # Ensure archive directory exists
      if (!dir.exists(self$archive_dir)) {
        dir.create(self$archive_dir, recursive = TRUE, showWarnings = FALSE)
      }

      # Default archive name from directory
      if (is.null(archive_name)) {
        archive_name <- basename(input_dir)
      }

      # Build archive path (always use zip for cross-platform compatibility)
      # Use absolute path to avoid issues when changing working directory
      archive_path <- normalizePath(
        file.path(self$archive_dir, paste0(archive_name, ".zip")),
        mustWork = FALSE
      )

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

      # Exclude files matching patterns
      if (!is.null(exclude_file_patterns) && length(exclude_file_patterns) > 0) {
        exclude_pattern <- paste(exclude_file_patterns, collapse = "|")
        excluded_count <- sum(grepl(exclude_pattern, files, ignore.case = TRUE))
        files <- files[!grepl(exclude_pattern, files, ignore.case = TRUE)]
        if (self$logs && excluded_count > 0) {
          message("  Excluded ", excluded_count, " files matching: ", exclude_pattern)
        }
      }

      if (length(files) == 0) {
        stop("No files found in directory: ", input_dir)
      }

      # Compute source hash BEFORE archiving (for change detection)
      source_hash <- private$compute_source_hash(files)

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
        compression_level = compression_level
      )

      if (!file.exists(archive_path)) {
        stop("Failed to create archive: ", archive_path)
      }

      if (self$logs) {
        size_mb <- round(file.size(archive_path) / 1024^2, 2)
        message("  Archive size: ", size_mb, " MB")
      }

      # Return both archive path and source hash for manifest
      list(
        archive_path = archive_path,
        source_hash = source_hash,
        file_count = length(files)
      )
    },

    # create_grouped_archive ----
    #' @description Create a zip archive from multiple directories.
    #' @param input_dirs Character vector of full paths to directories to archive.
    #' @param archive_name Name of the archive file (without extension).
    #' @param include_patterns File patterns to include (default: all files).
    #' @param exclude_file_patterns Regex patterns to exclude files (e.g., "\\.R$").
    #' @param compression_level Integer 0-9. 0 = no compression (fastest),
    #'   9 = maximum compression (smallest files). Default: 6.
    #' @return A list with archive_path, source_hash, and file_count.
    create_grouped_archive = function(
      input_dirs,
      archive_name,
      include_patterns = NULL,
      exclude_file_patterns = NULL,
      compression_level = 6L
    ) {
      # Validate all directories exist
      missing <- input_dirs[!dir.exists(input_dirs)]
      if (length(missing) > 0) {
        stop("Directories do not exist: ", paste(missing, collapse = ", "))
      }

      # Ensure archive directory exists
      if (!dir.exists(self$archive_dir)) {
        dir.create(self$archive_dir, recursive = TRUE, showWarnings = FALSE)
      }

      # Use absolute path to avoid issues when changing working directory
      archive_path <- normalizePath(
        file.path(self$archive_dir, paste0(archive_name, ".zip")),
        mustWork = FALSE
      )

      # Collect all files from all directories
      all_files <- character(0)
      for (input_dir in input_dirs) {
        files <- list.files(
          input_dir,
          recursive = TRUE,
          full.names = TRUE,
          all.files = FALSE
        )
        all_files <- c(all_files, files)
      }

      if (!is.null(include_patterns)) {
        pattern <- paste(include_patterns, collapse = "|")
        all_files <- all_files[grepl(pattern, all_files, ignore.case = TRUE)]
      }

      # Exclude files matching patterns
      if (!is.null(exclude_file_patterns) && length(exclude_file_patterns) > 0) {
        exclude_pattern <- paste(exclude_file_patterns, collapse = "|")
        excluded_count <- sum(grepl(exclude_pattern, all_files, ignore.case = TRUE))
        all_files <- all_files[!grepl(exclude_pattern, all_files, ignore.case = TRUE)]
        if (self$logs && excluded_count > 0) {
          message("  Excluded ", excluded_count, " files matching: ", exclude_pattern)
        }
      }

      if (length(all_files) == 0) {
        stop("No files found in directories: ", paste(input_dirs, collapse = ", "))
      }

      # Compute source hash BEFORE archiving (for change detection)
      source_hash <- private$compute_source_hash(all_files)

      if (self$logs) {
        message("Creating grouped zip archive: ", archive_path)
        message("  Sources: ", paste(basename(input_dirs), collapse = ", "))
        message("  Files: ", length(all_files))
      }

      # Remove existing archive if exists
      if (file.exists(archive_path)) {
        file.remove(archive_path)
      }

      # Get common parent directory
      parent_dir <- dirname(input_dirs[1])

      # Create zip archive with all directories
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(parent_dir)

      zip::zip(
        zipfile = archive_path,
        files = basename(input_dirs),
        recurse = TRUE,
        compression_level = compression_level
      )

      if (!file.exists(archive_path)) {
        stop("Failed to create archive: ", archive_path)
      }

      if (self$logs) {
        size_mb <- round(file.size(archive_path) / 1024^2, 2)
        message("  Archive size: ", size_mb, " MB")
      }

      # Return both archive path and source hash for manifest
      list(
        archive_path = archive_path,
        source_hash = source_hash,
        file_count = length(all_files)
      )
    },

    # create_input_archives ----
    #' @description Create zip archives for input directories.
    #' @param input_base Base directory containing inputs (default: "./inputs").
    #' @param directories Specific subdirectories to archive. If NULL, archives all.
    #' @param exclude_patterns Regex patterns to exclude directories/subfolders.
    #' @param exclude_file_patterns Regex patterns to exclude files (default: "\\.R$").
    #' @param group_by_prefix If TRUE, groups subfolders within each directory by prefix
    #'   (everything before the last underscore) into separate archives.
    #' @param compression_level Integer 0-9 passed to \code{zip::zip()}.
    #'   0 = no compression (fastest), 9 = maximum compression (smallest files).
    #'   Default: 6.
    #' @param multicore If TRUE, create archives in parallel across directories
    #'   using PSOCK clusters (parallelly/parallel) on all platforms.
    #' @param n_cores Number of cores for parallel execution
    #'   (default: \code{parallel::detectCores(logical = FALSE)}).
    #'   Capped to the number of archive tasks.
    #' @param update_gitignore If TRUE, automatically adds archived source files
    #'   to .gitignore (at the git repository root) unless they are already
    #'   covered by an existing pattern. Uses efficient directory/extension
    #'   patterns where possible (e.g., \code{inputs/mortality/*.fst}).
    #'   Requires git to be installed and the project to be in a git repository.
    #' @return A data.table with archive information.
    #' @examples
    #' \dontrun{
    #' # Archive specific directories
    #' manager$create_input_archives(directories = c("mortality", "pop"))
    #'
    #' # Exclude certain patterns
    #' manager$create_input_archives(exclude_patterns = c("^old_", "_backup$"))
    #'
    #' # Group subfolders within each directory by prefix
    #' # e.g., disease_burden/cancer_2020, disease_burden/cancer_2021 -> disease_burden_cancer.zip
    #' manager$create_input_archives(group_by_prefix = TRUE)
    #'
    #' # Exclude R files (default) and other patterns
    #' manager$create_input_archives(exclude_file_patterns = c("\\.R$", "\\.Rmd$"))
    #'
    #' # Parallel archive creation
    #' manager$create_input_archives(multicore = TRUE, n_cores = 4L)
    #'
    #' # Auto-update .gitignore with archived files
    #' manager$create_input_archives(update_gitignore = TRUE)
    #' }
    create_input_archives = function(
      input_base = "./inputs",
      directories = NULL,
      exclude_patterns = NULL,
      exclude_file_patterns = "\\.R$",
      group_by_prefix = FALSE,
      compression_level = 6L,
      multicore = FALSE,
      n_cores = parallel::detectCores(logical = FALSE),
      update_gitignore = FALSE
    ) {

      if (!dir.exists(input_base)) {
        stop("Input directory does not exist: ", input_base)
      }

      # Ensure archive directory exists
      if (!dir.exists(self$archive_dir)) {
        dir.create(self$archive_dir, recursive = TRUE, showWarnings = FALSE)
      }

      # Get subdirectories
      if (is.null(directories)) {
        directories <- list.dirs(input_base, full.names = FALSE, recursive = FALSE)
        directories <- directories[directories != ""]
      }

      if (length(directories) == 0) {
        stop("No directories found in: ", input_base)
      }

      # Apply exclusion patterns to top-level directories
      if (!is.null(exclude_patterns)) {
        exclude_regex <- paste(exclude_patterns, collapse = "|")
        excluded <- grepl(exclude_regex, directories)
        if (self$logs && any(excluded)) {
          message("Excluding directories: ", paste(directories[excluded], collapse = ", "))
        }
        directories <- directories[!excluded]

        if (length(directories) == 0) {
          stop("No directories remaining after applying exclude_patterns")
        }
      }

      # Phase 1: Collect archive tasks (fast, sequential)
      archive_tasks <- list()

      for (dir_name in directories) {
        dir_path <- file.path(input_base, dir_name)
        if (!dir.exists(dir_path)) {
          warning("Directory not found, skipping: ", dir_path)
          next
        }

        if (group_by_prefix) {
          # Get subfolders within this directory
          subfolders <- list.dirs(dir_path, full.names = FALSE, recursive = FALSE)
          subfolders <- subfolders[subfolders != ""]

          # Apply exclusion patterns to subfolders
          if (!is.null(exclude_patterns) && length(subfolders) > 0) {
            excluded <- grepl(exclude_regex, subfolders)
            if (self$logs && any(excluded)) {
              message("  Excluding subfolders in ", dir_name, ": ",
                      paste(subfolders[excluded], collapse = ", "))
            }
            subfolders <- subfolders[!excluded]
          }

          if (length(subfolders) == 0) {
            # No subfolders, archive the directory itself
            if (self$logs) {
              message("No subfolders in ", dir_name, ", archiving entire directory")
            }
            archive_tasks[[length(archive_tasks) + 1L]] <- list(
              type = "single",
              dir_name = dir_name,
              group_prefix = NA_character_,
              dir_path = dir_path,
              archive_name = NULL,
              exclude_file_patterns = exclude_file_patterns,
              compression_level = compression_level,
              archive_dir = self$archive_dir,
              logs = self$logs
            )
          } else {
            # Group subfolders by prefix (everything before last underscore)
            prefixes <- vapply(subfolders, function(sf) {
              if (grepl("_", sf)) {
                sub("_[^_]+$", "", sf)
              } else {
                sf  # No underscore, use full name as prefix
              }
            }, character(1))

            groups <- split(subfolders, prefixes)

            if (self$logs) {
              message("Grouping ", length(subfolders), " subfolders in ", dir_name,
                      " into ", length(groups), " archives:")
              for (prefix in names(groups)) {
                message("  ", prefix, ": ", paste(groups[[prefix]], collapse = ", "))
              }
            }

            # Create a task for each group
            for (prefix in names(groups)) {
              group_subfolders <- groups[[prefix]]
              subfolder_paths <- file.path(dir_path, group_subfolders)
              archive_name <- paste0(dir_name, "_", prefix)

              if (length(group_subfolders) == 1L) {
                archive_tasks[[length(archive_tasks) + 1L]] <- list(
                  type = "single",
                  dir_name = dir_name,
                  group_prefix = prefix,
                  dir_path = subfolder_paths,
                  archive_name = archive_name,
                  exclude_file_patterns = exclude_file_patterns,
                  compression_level = compression_level,
                  archive_dir = self$archive_dir,
                  logs = self$logs
                )
              } else {
                archive_tasks[[length(archive_tasks) + 1L]] <- list(
                  type = "grouped",
                  dir_name = dir_name,
                  group_prefix = prefix,
                  dir_paths = subfolder_paths,
                  archive_name = archive_name,
                  exclude_file_patterns = exclude_file_patterns,
                  compression_level = compression_level,
                  archive_dir = self$archive_dir,
                  logs = self$logs
                )
              }
            }
          }
        } else {
          # Original behavior: one archive per directory
          if (self$logs) {
            message("Creating archive for: ", dir_name)
          }

          archive_tasks[[length(archive_tasks) + 1L]] <- list(
            type = "single",
            dir_name = dir_name,
            group_prefix = NA_character_,
            dir_path = dir_path,
            archive_name = NULL,
            exclude_file_patterns = exclude_file_patterns,
            compression_level = compression_level,
            archive_dir = self$archive_dir,
            logs = self$logs
          )
        }
      }

      if (length(archive_tasks) == 0L) {
        return(data.table::data.table(
          directory = character(),
          group_prefix = character(),
          archive_path = character(),
          source_hash = character(),
          file_count = integer(),
          size_bytes = numeric()
        ))
      }

      # Phase 2: Execute archive creation
      # Check if parallel packages are available when multicore requested
      # Note: Always use PSOCK clusters (not forking) because zip::zip() C
      # internals are not fork-safe and segfault under mclapply.
      has_parallel_pkgs <- multicore &&
        requireNamespace("parallelly", quietly = TRUE) &&
        requireNamespace("parallel", quietly = TRUE)

      if (multicore && !has_parallel_pkgs) {
        warning(
          "Parallel packages (parallelly, parallel) not available. ",
          "Falling back to sequential execution."
        )
      }

      use_parallel <- has_parallel_pkgs && length(archive_tasks) > 1L
      actual_cores <- if (use_parallel) min(n_cores, length(archive_tasks)) else 1L

      if (use_parallel) {
        if (self$logs) {
          message("Creating ", length(archive_tasks), " archives in parallel (",
                  actual_cores, " cores)")
        }

        # PSOCK clusters on all platforms (separate R processes, not forks)
        cl <- parallelly::makeClusterPSOCK(
          actual_cores,
          dryrun = FALSE,
          quiet = !self$logs
        )
        on.exit(parallel::stopCluster(cl), add = TRUE)

        results <- parallel::parLapplyLB(
          cl = cl,
          X = archive_tasks,
          fun = .zenodo_archive_worker
        )
      } else {
        # Sequential execution
        if (self$logs && multicore && length(archive_tasks) <= 1L) {
          message("Only 1 archive task, running sequentially")
        }
        results <- lapply(archive_tasks, .zenodo_archive_worker)
      }

      # Filter out NULLs (failed tasks)
      results <- results[!vapply(results, is.null, logical(1))]

      if (length(results) == 0L) {
        return(data.table::data.table(
          directory = character(),
          group_prefix = character(),
          archive_path = character(),
          source_hash = character(),
          file_count = integer(),
          size_bytes = numeric()
        ))
      }

      results_dt <- data.table::rbindlist(results)

      # Update .gitignore with archived source files
      if (update_gitignore && nrow(results_dt) > 0L) {
        private$update_gitignore(archive_tasks)
      }

      results_dt
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

    # check_source_changes ----
    #' @description Check if local source files have changed since last upload.
    #' @param input_base Base directory containing inputs.
    #' @param exclude_file_patterns Regex patterns to exclude files (should match what was used during upload).
    #' @return A data.table showing which archives have changed.
    #' @examples
    #' \dontrun{
    #' changes <- manager$check_source_changes("./inputs")
    #' # Only re-archive and upload changed directories
    #' changed_dirs <- changes[changed == TRUE, directory]
    #' }
    check_source_changes = function(
      input_base = "./inputs",
      exclude_file_patterns = "\\.R$"
    ) {
      # Load existing manifest
      manifest <- self$load_manifest()
      if (is.null(manifest) || nrow(manifest) == 0) {
        if (self$logs) {
          message("No manifest found. All directories will be considered new.")
        }
        return(NULL)
      }

      if (self$logs) {
        message("Checking for source file changes...")
      }

      results <- lapply(seq_len(nrow(manifest)), function(i) {
        row <- manifest[i]
        dir_name <- row$directory
        group_prefix <- row$group_prefix
        stored_hash <- row$source_hash

        # Reconstruct the source path(s)
        if (is.na(group_prefix)) {
          # Single directory archive
          source_path <- file.path(input_base, dir_name)
          if (!dir.exists(source_path)) {
            return(data.table::data.table(
              directory = dir_name,
              group_prefix = group_prefix,
              archive_name = row$archive_name,
              stored_hash = stored_hash,
              current_hash = NA_character_,
              changed = NA,
              status = "missing"
            ))
          }

          # Get files (matching the archive creation logic)
          files <- list.files(source_path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
        } else {
          # Grouped archive - need to find subfolders matching the prefix
          dir_path <- file.path(input_base, dir_name)
          if (!dir.exists(dir_path)) {
            return(data.table::data.table(
              directory = dir_name,
              group_prefix = group_prefix,
              archive_name = row$archive_name,
              stored_hash = stored_hash,
              current_hash = NA_character_,
              changed = NA,
              status = "missing"
            ))
          }

          # Find subfolders that match this prefix
          subfolders <- list.dirs(dir_path, full.names = FALSE, recursive = FALSE)
          matching_subfolders <- subfolders[grepl(paste0("^", group_prefix, "_"), subfolders) |
                                              subfolders == group_prefix]

          if (length(matching_subfolders) == 0) {
            return(data.table::data.table(
              directory = dir_name,
              group_prefix = group_prefix,
              archive_name = row$archive_name,
              stored_hash = stored_hash,
              current_hash = NA_character_,
              changed = NA,
              status = "missing"
            ))
          }

          # Get all files from matching subfolders
          files <- character(0)
          for (sf in matching_subfolders) {
            sf_path <- file.path(dir_path, sf)
            sf_files <- list.files(sf_path, recursive = TRUE, full.names = TRUE, all.files = FALSE)
            files <- c(files, sf_files)
          }
        }

        # Apply file exclusion patterns
        if (!is.null(exclude_file_patterns) && length(exclude_file_patterns) > 0) {
          exclude_pattern <- paste(exclude_file_patterns, collapse = "|")
          files <- files[!grepl(exclude_pattern, files, ignore.case = TRUE)]
        }

        if (length(files) == 0) {
          return(data.table::data.table(
            directory = dir_name,
            group_prefix = group_prefix,
            archive_name = row$archive_name,
            stored_hash = stored_hash,
            current_hash = NA_character_,
            changed = NA,
            status = "empty"
          ))
        }

        # Compute current source hash
        current_hash <- private$compute_source_hash(files)
        changed <- current_hash != stored_hash

        data.table::data.table(
          directory = dir_name,
          group_prefix = group_prefix,
          archive_name = row$archive_name,
          stored_hash = stored_hash,
          current_hash = current_hash,
          changed = changed,
          status = if (changed) "modified" else "unchanged"
        )
      })

      result <- data.table::rbindlist(results)

      if (self$logs) {
        n_changed <- sum(result$changed, na.rm = TRUE)
        n_unchanged <- sum(!result$changed, na.rm = TRUE)
        n_missing <- sum(is.na(result$changed))
        message("  Changed: ", n_changed)
        message("  Unchanged: ", n_unchanged)
        if (n_missing > 0) message("  Missing/Empty: ", n_missing)
      }

      result
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

      # Use progress-enabled upload if callback is set
      if (!is.null(self$upload_progress)) {
        private$upload_with_progress(archive_path, rec)
      } else {
        if (self$logs) {
          message("Uploading: ", basename(archive_path))
          size_mb <- round(file.size(archive_path) / 1024^2, 2)
          message("  Size: ", size_mb, " MB")
        }

        self$zenodo$uploadFile(archive_path, record = rec)

        if (self$logs) {
          message("  Upload complete.")
        }
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

      # Use progress-enabled download if callback is set
      if (!is.null(self$download_progress)) {
        download_url <- private$get_file_download_url(filename)
        private$download_with_progress(download_url, dest_path, filename)
      } else {
        if (self$logs) {
          message("Downloading: ", filename)
        }

        self$record$downloadFiles(path = dest_dir)

        if (self$logs) {
          message("  Downloaded to: ", dest_path)
        }
      }

      if (!file.exists(dest_path)) {
        stop("Download failed for: ", filename)
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

      # Extract concept DOI from parent ID (zen4R stores it in parent$id)
      # The concept DOI is the DOI of the parent record (shared across versions)
      concept_id <- self$record$parent$id
      if (!is.null(concept_id)) {
        doi_prefix <- if (self$sandbox) "10.5072/zenodo." else "10.5281/zenodo."
        self$concept_doi <- paste0(doi_prefix, concept_id)
      } else if (!is.null(self$record$conceptdoi)) {
        # Fallback for older zen4R versions
        self$concept_doi <- self$record$conceptdoi
      }

      if (self$logs) {
        message("Record created with ID: ", self$record$id)
        message("  Draft DOI: ", self$record$pids$doi$identifier %||% "pending")
        message("  Concept DOI: ", self$concept_doi %||% "not available")
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

      # Update concept DOI from parent ID
      concept_id <- self$record$parent$id
      if (!is.null(concept_id)) {
        doi_prefix <- if (self$sandbox) "10.5072/zenodo." else "10.5281/zenodo."
        self$concept_doi <- paste0(doi_prefix, concept_id)
      } else if (!is.null(self$record$conceptdoi)) {
        # Fallback for older zen4R versions
        self$concept_doi <- self$record$conceptdoi
      }

      if (self$logs) {
        message("Record published!")
        message("  DOI: ", self$record$pids$doi$identifier)
        message("  Concept DOI: ", self$concept_doi %||% "not available")
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
      cat("  Progress callbacks:\n")
      cat("    Upload: ", if (!is.null(self$upload_progress)) "enabled" else "disabled", "\n")
      cat("    Download: ", if (!is.null(self$download_progress)) "enabled" else "disabled", "\n")
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

    # compute_source_hash ----
    # Compute a combined hash of multiple source files (for change detection)
    # This hash represents the content of source files BEFORE archiving
    compute_source_hash = function(file_paths) {
      if (length(file_paths) == 0) {
        return(NA_character_)
      }

      # Sort paths for consistent ordering
      file_paths <- sort(file_paths)

      # Compute individual file hashes
      hashes <- vapply(file_paths, function(f) {
        if (file.exists(f)) {
          digest::digest(file = f, algo = "xxhash64")
        } else {
          ""
        }
      }, character(1), USE.NAMES = FALSE)

      # Combine hashes into a single hash (hash of hashes)
      digest::digest(paste(hashes, collapse = "|"), algo = "xxhash64", serialize = FALSE)
    },

    # update_gitignore ----
    # Add archived source files to .gitignore unless already covered.
    # Uses `git check-ignore` for accurate pattern matching, then generates
    # efficient dir/ext patterns for any uncovered files.
    update_gitignore = function(archive_tasks) {
      # Find git root
      git_root <- tryCatch(
        trimws(system2("git", c("rev-parse", "--show-toplevel"),
                        stdout = TRUE, stderr = FALSE)),
        error = function(e) NULL,
        warning = function(w) NULL
      )
      if (is.null(git_root) || length(git_root) == 0L || !nzchar(git_root)) {
        if (self$logs) message("Not in a git repository; skipping .gitignore update")
        return(invisible(NULL))
      }

      # Reconstruct the list of all archived source files from task descriptors
      all_files <- character(0)
      for (task in archive_tasks) {
        input_paths <- if (task$type == "single") task$dir_path else task$dir_paths
        for (p in input_paths) {
          files <- list.files(p, recursive = TRUE, full.names = TRUE,
                              all.files = FALSE)
          if (!is.null(task$exclude_file_patterns) &&
              length(task$exclude_file_patterns) > 0) {
            exclude_pattern <- paste(task$exclude_file_patterns, collapse = "|")
            files <- files[!grepl(exclude_pattern, files, ignore.case = TRUE)]
          }
          all_files <- c(all_files, files)
        }
      }

      if (length(all_files) == 0L) return(invisible(NULL))

      # Convert to paths relative to git root
      git_root_norm <- normalizePath(git_root, winslash = "/")
      rel_files <- sub(
        paste0("^", gsub("([.\\\\|(){}^$*+?\\[\\]])", "\\\\\\1", git_root_norm), "/?"),
        "", normalizePath(all_files, winslash = "/")
      )

      # Use git check-ignore to find which files are already covered
      already_ignored <- tryCatch({
        system2("git", c("-C", git_root, "check-ignore", "--stdin"),
                input = rel_files, stdout = TRUE, stderr = FALSE)
      }, error = function(e) character(0),
         warning = function(w) character(0))

      not_ignored <- setdiff(rel_files, already_ignored)

      if (length(not_ignored) == 0L) {
        if (self$logs) {
          message("All ", length(rel_files),
                  " archived files already covered by .gitignore")
        }
        return(invisible(NULL))
      }

      # Generate efficient patterns grouped by directory + extension
      dirs <- dirname(not_ignored)
      exts <- tools::file_ext(not_ignored)

      # Build dir/*.ext patterns where 2+ files share dir+ext;
      # individual paths otherwise
      dt <- data.table::data.table(path = not_ignored, dir = dirs, ext = exts)
      patterns <- dt[, {
        if (.N >= 2L && nzchar(ext[1L])) {
          list(pattern = paste0(dir[1L], "/*.", ext[1L]))
        } else {
          list(pattern = path)
        }
      }, by = .(dir, ext)]$pattern
      patterns <- unique(patterns)

      # Read existing .gitignore and filter out patterns already present
      gitignore_path <- file.path(git_root, ".gitignore")
      existing <- if (file.exists(gitignore_path)) {
        trimws(readLines(gitignore_path, warn = FALSE))
      } else {
        character(0)
      }
      new_patterns <- patterns[!patterns %in% existing]

      if (length(new_patterns) == 0L) {
        if (self$logs) {
          message("All gitignore patterns already present")
        }
        return(invisible(NULL))
      }

      # Append new patterns under a header
      header <- "# Zenodo-archived input data (auto-generated by ZenodoAssetManager)"
      has_header <- any(grepl(
        "^# Zenodo-archived input data", existing, fixed = FALSE
      ))

      if (has_header) {
        # Append after the existing section (find last Zenodo-archived line)
        lines <- readLines(gitignore_path, warn = FALSE)
        header_idx <- max(grep("^# Zenodo-archived input data", lines))
        # Find the end of the section (next blank line or end of file)
        end_idx <- header_idx
        while (end_idx < length(lines) && nzchar(trimws(lines[end_idx + 1L]))) {
          end_idx <- end_idx + 1L
        }
        lines <- c(lines[1:end_idx], new_patterns,
                    if (end_idx < length(lines)) lines[(end_idx + 1L):length(lines)])
        writeLines(lines, gitignore_path)
      } else {
        # Append at end
        cat(c("", header, new_patterns),
            file = gitignore_path, sep = "\n", append = TRUE)
      }

      if (self$logs) {
        message("Added ", length(new_patterns),
                " patterns to .gitignore for archived files")
      }
      invisible(new_patterns)
    },

    # default_progress_callback ----
    # Default progress bar for console output
    default_progress_callback = function(bytes, total, filename) {
      if (total > 0) {
        pct <- round(100 * bytes / total, 1)
        bar_width <- 30
        filled <- round(bar_width * bytes / total)
        bar <- paste0(
          "[",
          paste(rep("=", filled), collapse = ""),
          paste(rep(" ", bar_width - filled), collapse = ""),
          "]"
        )
        size_mb <- round(bytes / 1024^2, 1)
        total_mb <- round(total / 1024^2, 1)
        cat(sprintf(
          "\r  %s %s %.1f/%.1f MB (%.1f%%)    ",
          basename(filename), bar, size_mb, total_mb, pct
        ))
        if (bytes >= total) cat("\n")
      }
    },

    # get_api_base_url ----
    # Get the base URL for Zenodo API
    get_api_base_url = function() {
      if (self$sandbox) {
        "https://sandbox.zenodo.org/api"
      } else {
        "https://zenodo.org/api"
      }
    },

    # get_token ----
    # Get the stored token from the ZenodoManager
    get_token = function() {
      if (is.null(self$zenodo)) {
        stop("Not connected. Call connect() first.")
      }
      # Access token from zen4R's ZenodoManager
      self$zenodo$.__enclos_env__$private$token
    },

    # upload_with_progress ----
    # Upload a file via zen4R (which handles the correct
    # Zenodo API protocol internally). Logs progress.
    upload_with_progress = function(file_path, record) {
      filename <- basename(file_path)
      file_size <- file.size(file_path)

      if (self$logs) {
        message(
          "Uploading: ", filename,
          " (", round(file_size / 1024^2, 2), " MB)"
        )
      }

      # Delegate to zen4R which handles API details
      self$zenodo$uploadFile(file_path, record = record)

      if (self$logs) {
        message("  Upload complete.")
      }

      invisible(TRUE)
    },

    # download_with_progress ----
    # Download a file with progress callback using httr2
    download_with_progress = function(url, dest_path, filename = NULL, total_size = NULL) {
      if (!requireNamespace("httr2", quietly = TRUE)) {
        stop("Package 'httr2' required for progress callbacks. Install with: install.packages('httr2')")
      }

      if (is.null(filename)) {
        filename <- basename(dest_path)
      }

      token <- private$get_token()

      # Try to get file size from record metadata if not provided
      if (is.null(total_size) && !is.null(self$record)) {
        files_info <- self$record$listFiles(pretty = TRUE)
        if (!is.null(files_info) && nrow(files_info) > 0) {
          file_row <- files_info[files_info$filename == filename, ]
          if (nrow(file_row) > 0 && "size" %in% names(file_row)) {
            total_size <- as.numeric(file_row$size[1])
          }
        }
      }
      total_size <- total_size %||% 0

      if (self$logs) {
        message("Downloading: ", filename)
        if (total_size > 0) {
          message("  Size: ", round(total_size / 1024^2, 2), " MB")
        }
      }

      # Build request
      req <- httr2::request(url) |>
        httr2::req_headers("Authorization" = paste("Bearer", token))

      # Perform request and save to file
      resp <- httr2::req_perform(req, path = dest_path)

      if (httr2::resp_status(resp) >= 400) {
        stop("Download failed with status ", httr2::resp_status(resp))
      }

      if (self$logs) {
        size_mb <- round(file.size(dest_path) / 1024^2, 2)
        message("  Downloaded: ", size_mb, " MB")
      }

      invisible(dest_path)
    },

    # get_file_download_url ----
    # Get download URL for a specific file in a record
    get_file_download_url = function(filename) {
      if (is.null(self$record)) {
        stop("No record loaded.")
      }

      files_info <- self$record$listFiles(pretty = TRUE)
      if (is.null(files_info) || nrow(files_info) == 0) {
        stop("No files in record.")
      }

      # Find matching file
      file_row <- files_info[files_info$filename == filename, ]
      if (nrow(file_row) == 0) {
        stop("File not found in record: ", filename)
      }

      # Get download link
      file_row$download_link[1]
    },

    # sync_check ----
    sync_check = function(remote_files, input_base, directories) {
      if (self$logs) {
        message("\n=== Checking sync status ===")
      }

      # Identify which archives correspond to which directories
      results <- lapply(seq_len(nrow(remote_files)), function(i) {
        fname <- remote_files$filename[i]
        archive_name <- gsub("\\.zip$", "", fname, ignore.case = TRUE)

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

      # Download files - use progress if callback is set
      if (!is.null(self$download_progress)) {
        # Download each file individually with progress
        for (i in seq_len(nrow(remote_files))) {
          fname <- remote_files$filename[i]
          archive_name <- gsub("\\.zip$", "", fname, ignore.case = TRUE)

          # Check if this directory should be processed
          if (!is.null(directories) && !(archive_name %in% directories)) {
            if (self$logs) message("Skipping (not in list): ", archive_name)
            next
          }

          # Get download URL and file size
          download_url <- private$get_file_download_url(fname)
          file_size <- if ("size" %in% names(remote_files)) {
            as.numeric(remote_files$size[i])
          } else {
            NULL
          }

          dest_path <- file.path(download_dir, fname)
          private$download_with_progress(download_url, dest_path, fname, file_size)
        }
      } else {
        # Download all files from record at once (no progress)
        self$record$downloadFiles(path = download_dir)
      }

      # Extract each archive
      for (i in seq_len(nrow(remote_files))) {
        fname <- remote_files$filename[i]
        archive_name <- gsub("\\.zip$", "", fname, ignore.case = TRUE)

        # Check if this directory should be processed
        if (!is.null(directories) && !(archive_name %in% directories)) {
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

      # Save manifest with source hashes (for change detection)
      manifest <- data.table::data.table(
        directory = archives$directory,
        group_prefix = archives$group_prefix,
        archive_name = basename(archives$archive_path),
        source_hash = archives$source_hash,
        file_count = archives$file_count,
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
