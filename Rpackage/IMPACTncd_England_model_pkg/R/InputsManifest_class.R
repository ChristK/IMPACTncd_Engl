## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by
## Chris Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2026 University of Liverpool, Chris Kypridemos
##
## IMPACTncd_Engl is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by the Free
## Software Foundation; either version 3 of the License, or (at your option) any
## later version. This program is distributed in the hope that it will be
## useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
## Public License for more details. You should have received a copy of the GNU
## General Public License along with this program; if not, see
## <http://www.gnu.org/licenses/> or write to the Free Software Foundation,
## Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


#' Hash a single file for the inputs manifest
#'
#' @description
#' Standalone worker function for parallel hashing of input files. Defined
#' outside the R6 class so it can be serialized to PSOCK cluster workers
#' without dragging along non-serializable R6 environments or HTTP connections.
#'
#' @param file_path Character. Absolute path to the file to hash.
#' @param algo Character. Hash algorithm passed to [digest::digest()]
#'   (default: `"xxhash64"`).
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{file_path}{The input path (unchanged).}
#'     \item{hash}{The computed hash string.}
#'     \item{size_bytes}{File size in bytes.}
#'     \item{mtime}{File modification time as ISO 8601 string.}
#'   }
#'
#' @details
#' This function follows the project convention for PSOCK-safe worker
#' functions: it receives only plain serializable data (character scalars) and
#' returns a plain list. See also `.zenodo_archive_worker()` in
#' `ZenodoManager_class.R` for the same pattern.
#'
#' @seealso `InputsManifest` which calls this function via
#'   [parallel::parLapplyLB()].
#'
#' @keywords internal
.manifest_hash_worker <- function(file_path, algo = "xxhash64") {
  list(
    file_path = file_path,
    hash = digest::digest(file = file_path, algo = algo),
    size_bytes = file.size(file_path),
    mtime = format(file.mtime(file_path), "%Y-%m-%dT%H:%M:%S")
  )
}


#' R6 Class for unified input file tracking and staleness detection
#'
#' @description
#' The `InputsManifest` class provides a unified system for tracking all input
#' data files used by the IMPACTncd simulation. It detects file changes and
#' determines which downstream artifacts (synthpops, compiled RR `.fst` files,
#' PARF files) need regeneration based on the model's causal structure.
#'
#' @details
#' This class replaces the previous fragmented tracking mechanisms
#' (`fileversion.csv`, synthpop meta YAML, `zenodo_manifest.csv`) with a
#' single manifest CSV that records relative paths and `xxhash64` checksums
#' for all data files in the `inputs/` directory.
#'
#' ## Manifest CSV format
#'
#' The manifest is stored at `simulation/inputs_manifest.csv` with columns:
#'
#' | Column | Type | Description |
#' |--------|------|-------------|
#' | `relative_path` | character | Path relative to `inputs/` (portable) |
#' | `hash` | character | xxhash64 checksum of file contents |
#' | `size_bytes` | numeric | File size in bytes |
#' | `mtime` | character | ISO 8601 modification time |
#' | `category` | character | Top-level subdirectory name |
#' | `zenodo_managed` | logical | Whether the file is uploaded to Zenodo |
#'
#' ## Hashing strategy
#'
#' All hashes use `xxhash64` via [digest::digest()]. This is ~10x faster than
#' MD5, which is critical when scanning 10,000+ parquet files in
#' `exposure_distributions/`. The `compare()` method applies an **mtime
#' pre-filter**: files whose modification time has not changed since the last
#' manifest generation are assumed unchanged and not re-hashed, making
#' incremental checks very fast.
#'
#' ## Causal invalidation
#'
#' The `detect_changes()` method maps file changes to downstream artifacts
#' using the model's causal structure:
#'
#' \itemize{
#'   \item Changes in `pop_estimates_lsoa/`, `exposure_distributions/`,
#'     `mortality/`, or `pop_projections/` invalidate **synthpops**
#'   \item Changes in `RR/*.csvy` files invalidate the corresponding
#'     **compiled RR `.fst`** files and downstream **PARF** files
#'   \item Changes in `disease_burden/` data files (not `.yaml`) invalidate
#'     that disease's **PARF** files
#'   \item A BFS traversal through the causal graph (built from
#'     `design$RR` names) finds **transitively affected** diseases
#' }
#'
#' ## Output artifact checksums
#'
#' Synthpop and PARF files are too large to re-hash at every simulation
#' start. Instead, their checksums are embedded in their filenames (e.g.
#' `synthpop_<hash>_<mc>.fst`, `PARF_<disease>_<hash>.qs`). When upstream
#' inputs change, the checksum changes, the filename changes, and the old
#' file is simply not found — triggering automatic regeneration. The
#' `verify_artifact()` method provides an opt-in diagnostic that actually
#' hashes the file on disk and compares against the filename-embedded hash.
#'
#' ## Zenodo integration
#'
#' Not all files in `inputs/` are uploaded to Zenodo — some are tracked by
#' git instead. The `zenodo_managed` column distinguishes between the two.
#' Use `set_zenodo_dirs()` to mark which directories are Zenodo-managed and
#' `get_zenodo_format()` to extract the Zenodo-managed subset for comparison
#' with the upload snapshot.
#'
#' ## Parallel hashing
#'
#' When `parallel = TRUE` is passed to `generate()`, the class uses
#' [parallelly::makeClusterPSOCK()] + [parallel::parLapplyLB()] with the
#' standalone `.manifest_hash_worker()` function. This follows the project
#' convention of never forking (which causes segfaults with `zip::zip()`
#' and cannot serialize R6 HTTP connections).
#'
#' @section Typical workflow:
#' The manifest is managed automatically by [Simulation$new()][Simulation]:
#' \enumerate{
#'   \item On first run, `generate()` scans `inputs/` and creates the manifest
#'   \item On subsequent runs, `detect_changes()` compares the filesystem
#'     against the stored manifest
#'   \item If changes are found, stale artifacts are deleted and the manifest
#'     is regenerated
#'   \item Existing generation logic (in SynthPop, ExposureEffect, Disease)
#'     automatically recreates deleted artifacts
#' }
#'
#' @examples
#' \dontrun{
#' # Create and generate a manifest
#' manifest <- InputsManifest$new(
#'   inputs_dir = "./inputs",
#'   manifest_path = "./simulation/inputs_manifest.csv"
#' )
#' manifest$generate(parallel = TRUE, n_workers = 4L)
#' print(manifest)
#'
#' # Check for changes since last generation
#' diff <- manifest$compare()
#' cat("Added:", length(diff$added), "\n")
#' cat("Modified:", length(diff$modified), "\n")
#'
#' # Detect which artifacts are affected (requires Design object)
#' changes <- manifest$detect_changes(design = my_design)
#' if (changes$affected_synthpops) {
#'   message("Synthpops need regeneration")
#' }
#'
#' # Get aggregate hash for synthpop-relevant inputs
#' hash <- manifest$get_synthpop_input_hash()
#'
#' # Mark Zenodo-managed directories
#' manifest$set_zenodo_dirs(c("RR", "exposure_distributions"))
#' zenodo_subset <- manifest$get_zenodo_format()
#'
#' # Verify integrity of a large artifact file (opt-in diagnostic)
#' result <- manifest$verify_artifact("simulation/parf/PARF_chd_8300ca5a.qs")
#' cat("Valid:", result$valid, "\n")
#' }
#'
#' @seealso
#' \itemize{
#'   \item `Simulation` — creates and manages the manifest during initialisation
#'   \item `ZenodoAssetManager` — uses the manifest for upload comparison via
#'     `compare_with_inputs_manifest()`
#'   \item `.manifest_hash_worker()` — standalone worker function for PSOCK
#'     parallel hashing
#' }
#'
#' @family data-management
#'
#' @export
InputsManifest <- R6::R6Class(
  classname = "InputsManifest",

  public = list(

    #' @field manifest data.table with columns: relative_path, hash,
    #'   size_bytes, mtime, category, zenodo_managed
    manifest = NULL,

    # initialize ----
    #' @description Create a new InputsManifest object.
    #'
    #' Initialises the manifest paths but does **not** scan files or load an
    #' existing manifest. Call `generate()` to scan or `load()` to read a
    #' previously saved manifest.
    #'
    #' @param inputs_dir Character. Path to inputs directory (default:
    #'   `"./inputs"`). Will be normalised internally.
    #' @param manifest_path Character. Path to the manifest CSV file
    #'   (default: `"./simulation/inputs_manifest.csv"`).
    #' @return A new `InputsManifest` object (invisibly).
    initialize = function(
      inputs_dir = "./inputs",
      manifest_path = "./simulation/inputs_manifest.csv"
    ) {
      private$inputs_dir <- normalizePath(inputs_dir, mustWork = FALSE)
      private$manifest_path <- manifest_path
      invisible(self)
    },

    # generate ----
    #' @description Scan the inputs directory and compute hashes for all data
    #'   files, then save the manifest to disk.
    #'
    #' This is the primary method for creating or refreshing the manifest.
    #' It lists all files matching `file_patterns` under the inputs directory,
    #' hashes each file with `xxhash64`, categorises them by top-level
    #' subdirectory, and writes the result to `manifest_path`.
    #'
    #' When `parallel = TRUE` and more than 100 files are found, hashing is
    #' distributed across PSOCK workers using `.manifest_hash_worker()`.
    #'
    #' @param file_patterns Character vector of regex patterns for file
    #'   extensions to include (default: parquet, fst, csv, csvy, qs, yaml).
    #' @param parallel Logical. If `TRUE`, use a PSOCK cluster for hashing
    #'   (recommended for large input sets).
    #' @param n_workers Integer. Number of PSOCK workers (default: `4L`).
    #' @param zenodo_dirs Character vector of directory names (relative to
    #'   `inputs_dir`) that are managed by Zenodo. Files in these directories
    #'   will have `zenodo_managed = TRUE`. Pass `NULL` to leave all as `FALSE`.
    #' @return Invisible self (for method chaining).
    generate = function(
      file_patterns = c(
        "\\.parquet$", "\\.fst$", "\\.csv$", "\\.csvy$",
        "\\.qs$", "\\.yaml$"
      ),
      parallel = FALSE,
      n_workers = 4L,
      zenodo_dirs = NULL
    ) {
      all_files <- private$scan_files(file_patterns)

      if (length(all_files) == 0L) {
        self$manifest <- data.table::data.table(
          relative_path = character(0),
          hash = character(0),
          size_bytes = numeric(0),
          mtime = character(0),
          category = character(0),
          zenodo_managed = logical(0)
        )
        self$save()
        return(invisible(self))
      }

      if (parallel && length(all_files) > 100L) {
        cl <- parallelly::makeClusterPSOCK(
          n_workers,
          rscript_startup = quote(library(digest))
        )
        on.exit(parallel::stopCluster(cl), add = TRUE)

        results <- parallel::parLapplyLB(
          cl, all_files,
          function(f) .manifest_hash_worker(f, algo = "xxhash64")
        )
      } else {
        results <- lapply(all_files, function(f) {
          .manifest_hash_worker(f, algo = "xxhash64")
        })
      }

      inputs_dir_norm <- private$inputs_dir
      self$manifest <- data.table::rbindlist(lapply(results, function(r) {
        rel_path <- private$to_relative(r$file_path)
        data.table::data.table(
          relative_path = rel_path,
          hash = r$hash,
          size_bytes = r$size_bytes,
          mtime = r$mtime,
          category = private$categorize_file(rel_path),
          zenodo_managed = FALSE
        )
      }))

      # Set zenodo_managed for configured directories
      if (!is.null(zenodo_dirs) && length(zenodo_dirs) > 0L) {
        self$set_zenodo_dirs(zenodo_dirs)
      }

      data.table::setkey(self$manifest, "relative_path")
      self$save()
      invisible(self)
    },

    # load ----
    #' @description Load an existing manifest from disk.
    #' @return Invisible self.
    load = function() {
      if (!file.exists(private$manifest_path)) {
        stop(
          "Manifest file not found: ", private$manifest_path,
          ". Run generate() first."
        )
      }
      self$manifest <- data.table::fread(
        private$manifest_path,
        colClasses = c(
          relative_path = "character",
          hash = "character",
          size_bytes = "numeric",
          mtime = "character",
          category = "character",
          zenodo_managed = "logical"
        )
      )
      # Backward compatibility: add zenodo_managed if missing
      if (!"zenodo_managed" %in% names(self$manifest)) {
        self$manifest[, zenodo_managed := FALSE]
      }
      data.table::setkey(self$manifest, "relative_path")
      invisible(self)
    },

    # save ----
    #' @description Write current manifest to disk.
    #' @return Invisible self.
    save = function() {
      if (is.null(self$manifest)) {
        stop("No manifest to save. Call generate() first.")
      }
      dir.create(
        dirname(private$manifest_path),
        showWarnings = FALSE,
        recursive = TRUE
      )
      data.table::fwrite(self$manifest, private$manifest_path)
      invisible(self)
    },

    # compare ----
    #' @description Compare the current filesystem state against the stored
    #'   manifest to detect added, removed, and modified files.
    #'
    #' Uses an **mtime pre-filter** for performance: files whose modification
    #' time has not changed since the manifest was generated are assumed
    #' unchanged and not re-hashed. Only files with a different mtime are
    #' re-hashed to confirm whether their content actually changed.
    #'
    #' Requires a manifest to be loaded (via `load()` or `generate()`).
    #'
    #' @param file_patterns Character vector of regex patterns for file
    #'   extensions (same format as `generate()`).
    #' @return A named list with four character vectors:
    #'   \describe{
    #'     \item{added}{Files present on disk but not in the manifest.}
    #'     \item{removed}{Files in the manifest but no longer on disk.}
    #'     \item{modified}{Files whose hash has changed since the manifest.}
    #'     \item{unchanged}{Files whose hash matches the manifest.}
    #'   }
    compare = function(
      file_patterns = c(
        "\\.parquet$", "\\.fst$", "\\.csv$", "\\.csvy$",
        "\\.qs$", "\\.yaml$"
      )
    ) {
      if (is.null(self$manifest) || nrow(self$manifest) == 0L) {
        stop("No manifest loaded. Call generate() or load() first.")
      }

      current_files <- private$scan_files(file_patterns)
      current_rel <- vapply(
        current_files, private$to_relative, character(1),
        USE.NAMES = FALSE
      )

      stored_paths <- self$manifest$relative_path

      added <- setdiff(current_rel, stored_paths)
      removed <- setdiff(stored_paths, current_rel)
      common <- intersect(current_rel, stored_paths)

      modified <- character(0)
      unchanged <- character(0)

      if (length(common) > 0L) {
        # Use mtime pre-filter: only re-hash files whose mtime changed
        common_dt <- data.table::data.table(relative_path = common)
        common_dt <- self$manifest[common_dt, on = "relative_path"]

        for (i in seq_len(nrow(common_dt))) {
          f <- common_dt$relative_path[i]
          full_path <- file.path(private$inputs_dir, f)
          current_mtime <- format(
            file.mtime(full_path), "%Y-%m-%dT%H:%M:%S"
          )

          if (identical(current_mtime, common_dt$mtime[i])) {
            # mtime unchanged -> assume content unchanged (fast path)
            unchanged <- c(unchanged, f)
          } else {
            # mtime changed -> re-hash to confirm
            current_hash <- digest::digest(
              file = full_path, algo = "xxhash64"
            )
            if (identical(current_hash, common_dt$hash[i])) {
              unchanged <- c(unchanged, f)
            } else {
              modified <- c(modified, f)
            }
          }
        }
      }

      list(
        added = added,
        removed = removed,
        modified = modified,
        unchanged = unchanged
      )
    },

    # get_category_hash ----
    #' @description Get aggregate hash for all files in a top-level category
    #'   (subdirectory of `inputs/`).
    #'
    #' The aggregate hash is computed by sorting the individual file hashes
    #' and digesting them together, so it changes when any file in the
    #' category is added, removed, or modified.
    #'
    #' @param category Character. Category name matching a top-level
    #'   subdirectory (e.g., `"RR"`, `"exposure_distributions"`,
    #'   `"mortality"`).
    #' @return Character hash string, or `NA_character_` if the category
    #'   contains no files.
    get_category_hash = function(category) {
      if (is.null(self$manifest)) return(NA_character_)
      cat_ <- category  # avoid data.table column name collision
      rows <- self$manifest[category == cat_]
      if (nrow(rows) == 0L) return(NA_character_)
      private$compute_aggregate_hash(sort(rows$hash))
    },

    # get_directory_hash ----
    #' @description Get aggregate hash for a directory path.
    #' @param dir_path Character. Directory path relative to inputs_dir.
    #' @return Character hash, or NA_character_ if no files found.
    get_directory_hash = function(dir_path) {
      if (is.null(self$manifest)) return(NA_character_)
      # Ensure trailing slash for prefix matching
      prefix <- sub("/$", "", dir_path)
      rows <- self$manifest[grepl(paste0("^", prefix, "/"), relative_path)]
      if (nrow(rows) == 0L) return(NA_character_)
      private$compute_aggregate_hash(sort(rows$hash))
    },

    # get_file_hash ----
    #' @description Get hash for a single file.
    #' @param file_path Character. Relative path from inputs_dir.
    #' @return Character hash, or NA_character_ if not found.
    get_file_hash = function(file_path) {
      if (is.null(self$manifest)) return(NA_character_)
      row <- self$manifest[relative_path == file_path]
      if (nrow(row) == 0L) return(NA_character_)
      row$hash[1]
    },

    # get_synthpop_input_hash ----
    #' @description Compute a single aggregate hash representing all input
    #'   files that synthetic populations depend on.
    #'
    #' This hash is included in the synthpop checksum (see
    #' `SynthPop$gen_checksum()`), so that when any upstream data file
    #' changes, the synthpop filename changes and regeneration is triggered
    #' automatically.
    #'
    #' The categories included are: `pop_estimates_lsoa`,
    #' `exposure_distributions`, `mortality`, and `pop_projections`.
    #'
    #' @return Character hash string, or `""` if none of the categories
    #'   exist in the manifest.
    get_synthpop_input_hash = function() {
      categories <- c(
        "pop_estimates_lsoa", "exposure_distributions",
        "mortality", "pop_projections"
      )
      hashes <- vapply(
        categories, self$get_category_hash, character(1)
      )
      # Remove NAs for categories that don't exist
      hashes <- hashes[!is.na(hashes)]
      if (length(hashes) == 0L) return("")
      private$compute_aggregate_hash(hashes)
    },

    # detect_changes ----
    #' @description Detect file changes and determine which downstream
    #'   simulation artifacts are affected, using the model's causal
    #'   structure.
    #'
    #' This is the primary method called by `Simulation$new()` to decide
    #' what needs regeneration. It:
    #' \enumerate{
    #'   \item Calls `compare()` to find added/removed/modified files
    #'   \item Maps changed files to categories (RR, exposure_distributions,
    #'     disease_burden, etc.)
    #'   \item For RR changes, parses `exposure~outcome` from filenames
    #'   \item Uses BFS through the causal graph (from `design$RR` names)
    #'     to find transitively affected diseases
    #'   \item `disease_burden/*.yaml` files are tracked but do **not**
    #'     trigger invalidation
    #' }
    #'
    #' @param design A `Design` object containing the `RR` list whose names
    #'   define the causal structure (e.g., `"bmi~chd"`, `"chd~hf"`). Pass
    #'   `NULL` to skip causal-structure-based invalidation (only category
    #'   mapping will be used).
    #' @return A named list:
    #'   \describe{
    #'     \item{added}{Character vector of newly added files.}
    #'     \item{removed}{Character vector of deleted files.}
    #'     \item{modified}{Character vector of modified files.}
    #'     \item{all_changed}{Union of added, removed, and modified.}
    #'     \item{affected_synthpops}{Logical: do synthpops need regeneration?}
    #'     \item{affected_rr}{Character vector of affected `exposure~outcome`
    #'       RR relationships.}
    #'     \item{affected_diseases}{Character vector of disease names whose
    #'       PARFs need regeneration (includes transitively affected).}
    #'     \item{changes_dt}{A `data.table` with columns `relative_path`,
    #'       `status`, and `category` for detailed inspection.}
    #'   }
    detect_changes = function(design = NULL) {
      diff <- self$compare()
      all_changed <- c(diff$added, diff$removed, diff$modified)

      if (length(all_changed) == 0L) {
        return(list(
          added = character(0),
          removed = character(0),
          modified = character(0),
          all_changed = character(0),
          affected_synthpops = FALSE,
          affected_rr = character(0),
          affected_diseases = character(0),
          changes_dt = data.table::data.table(
            relative_path = character(0),
            status = character(0),
            category = character(0)
          )
        ))
      }

      changes_dt <- data.table::data.table(
        relative_path = all_changed,
        status = c(
          rep("added", length(diff$added)),
          rep("removed", length(diff$removed)),
          rep("modified", length(diff$modified))
        )
      )
      changes_dt[, category := vapply(
        relative_path, private$categorize_file, character(1)
      )]

      # Determine affected synthpops
      synthpop_categories <- c(
        "pop_estimates_lsoa", "exposure_distributions",
        "mortality", "pop_projections"
      )
      affected_synthpops <- any(changes_dt$category %in% synthpop_categories)

      # Determine affected RR relationships
      rr_changes <- changes_dt[category == "RR"]
      affected_rr <- character(0)
      if (nrow(rr_changes) > 0L) {
        # Extract exposure~outcome from filename (e.g., "RR/bmi~chd.csvy")
        affected_rr <- gsub(
          "\\.csvy$", "",
          basename(rr_changes$relative_path)
        )
        # Filter out files from unused/ subdirectory
        unused_mask <- grepl("^RR/unused/", rr_changes$relative_path)
        affected_rr <- affected_rr[!unused_mask]
      }

      # Determine affected diseases using causal structure
      affected_diseases <- character(0)

      # From disease_burden changes
      db_changes <- changes_dt[
        category == "disease_burden" &
          !grepl("\\.yaml$", relative_path)
      ]
      if (nrow(db_changes) > 0L) {
        # Extract disease name: disease_burden/{disease_name}_dur/... or
        # disease_burden/{disease_name}_incd/...
        dirs <- vapply(
          strsplit(db_changes$relative_path, "/"),
          function(x) if (length(x) >= 2L) x[2] else "",
          character(1)
        )
        # Strip _dur, _incd suffixes to get disease name
        disease_names <- unique(sub("_(dur|incd).*$", "", dirs))
        affected_diseases <- disease_names[nzchar(disease_names)]
      }

      # From RR changes: trace through causal graph
      if (length(affected_rr) > 0L && !is.null(design)) {
        # The outcome (disease) is the part after ~
        disease_from_rr <- unique(sub("^.*~", "", affected_rr))

        # Use design$RR names to build causal edges for transitive lookup
        if (!is.null(design$RR) && length(design$RR) > 0L) {
          rr_names <- names(design$RR)
          # Build adjacency: exposure -> outcome
          edges <- strsplit(rr_names, "~")
          # Find all diseases transitively affected
          # (e.g., if bmi~chd changes and chd~hf exists, hf is also affected)
          all_outcomes <- unique(vapply(
            edges, function(e) e[2], character(1)
          ))
          # Simple BFS through causal chain
          frontier <- disease_from_rr
          visited <- character(0)
          while (length(frontier) > 0L) {
            visited <- unique(c(visited, frontier))
            next_frontier <- character(0)
            for (d in frontier) {
              # Find edges where d is the exposure
              downstream <- vapply(
                edges,
                function(e) if (e[1] == d) e[2] else NA_character_,
                character(1)
              )
              downstream <- downstream[!is.na(downstream)]
              next_frontier <- c(
                next_frontier,
                downstream[!downstream %in% visited]
              )
            }
            frontier <- unique(next_frontier)
          }
          affected_diseases <- unique(c(affected_diseases, visited))
        } else {
          affected_diseases <- unique(c(affected_diseases, disease_from_rr))
        }
      }

      list(
        added = diff$added,
        removed = diff$removed,
        modified = diff$modified,
        all_changed = all_changed,
        affected_synthpops = affected_synthpops,
        affected_rr = affected_rr,
        affected_diseases = affected_diseases,
        changes_dt = changes_dt
      )
    },

    # set_zenodo_dirs ----
    #' @description Mark directories as Zenodo-managed.
    #'
    #' Sets `zenodo_managed = TRUE` for all manifest entries whose
    #' `relative_path` starts with one of the given directory names. Files
    #' not in these directories are set to `FALSE`. This controls which files
    #' are included in the Zenodo upload comparison (see `get_zenodo_format()`
    #' and `ZenodoAssetManager$compare_with_inputs_manifest()`).
    #'
    #' @param dirs Character vector of directory names relative to
    #'   `inputs_dir` (e.g., `c("RR", "exposure_distributions",
    #'   "disease_burden")`).
    #' @return Invisible self (for method chaining).
    set_zenodo_dirs = function(dirs) {
      if (is.null(self$manifest)) {
        stop("No manifest loaded. Call generate() or load() first.")
      }
      self$manifest[, zenodo_managed := FALSE]
      for (d in dirs) {
        prefix <- sub("/$", "", d)
        self$manifest[
          grepl(paste0("^", prefix, "(/|$)"), relative_path),
          zenodo_managed := TRUE
        ]
      }
      invisible(self)
    },

    # get_zenodo_format ----
    #' @description Return the Zenodo-managed subset of the manifest.
    #'
    #' Filters to files with `zenodo_managed == TRUE` and returns only the
    #' columns needed for comparison with the Zenodo upload snapshot. Used
    #' by `ZenodoAssetManager$compare_with_inputs_manifest()`.
    #'
    #' @return A `data.table` with columns `relative_path` and `hash`, or
    #'   `NULL` if no manifest is loaded.
    get_zenodo_format = function() {
      if (is.null(self$manifest)) return(NULL)
      self$manifest[
        zenodo_managed == TRUE,
        .(relative_path, hash)
      ]
    },

    # verify_artifact ----
    #' @description Verify integrity of a large output artifact by comparing
    #'   the hash embedded in its filename against the actual file hash.
    #'
    #' This is an **opt-in diagnostic** — it is never called during normal
    #' simulation startup because synthpop and PARF files can be several
    #' gigabytes. Use it before a Zenodo upload, after a storage migration,
    #' or when debugging suspected data corruption.
    #'
    #' Supported filename patterns:
    #' \itemize{
    #'   \item `synthpop_<hash>_<mc>.fst`
    #'   \item `PARF_<disease>_<hash>.qs`
    #'   \item `PARF_<disease>_<hash>/` (parquet directory)
    #' }
    #'
    #' @param path Character. Path to the artifact file to verify.
    #' @param algo Character. Hash algorithm (default: `"xxhash64"`). Must
    #'   match the algorithm used when the artifact was created.
    #' @return A named list:
    #'   \describe{
    #'     \item{expected_hash}{The hash extracted from the filename.}
    #'     \item{actual_hash}{The hash computed from the file contents.}
    #'     \item{valid}{Logical: do the hashes match? `NA` if the filename
    #'       pattern was not recognised.}
    #'   }
    verify_artifact = function(path, algo = "xxhash64") {
      if (!file.exists(path)) {
        stop("File not found: ", path)
      }

      # Extract expected hash from filename
      # Patterns: synthpop_<hash>_<mc>.fst, PARF_<disease>_<hash>.qs,
      #           PARF_<disease>_<hash>/ (parquet dir)
      bname <- basename(path)
      expected_hash <- NA_character_

      if (grepl("^synthpop_", bname)) {
        # synthpop_<hash>_<mc>.fst
        parts <- strsplit(sub("\\.fst$", "", bname), "_")[[1]]
        if (length(parts) >= 3L) {
          expected_hash <- parts[2]
        }
      } else if (grepl("^PARF_", bname)) {
        # PARF_<disease>_<hash>.qs or PARF_<disease>_<hash>
        parts <- strsplit(sub("\\.(qs|parquet)$", "", bname), "_")[[1]]
        if (length(parts) >= 3L) {
          expected_hash <- parts[length(parts)]
        }
      }

      actual_hash <- digest::digest(file = path, algo = algo)

      list(
        expected_hash = expected_hash,
        actual_hash = actual_hash,
        valid = if (is.na(expected_hash)) {
          NA
        } else {
          identical(expected_hash, actual_hash)
        }
      )
    },

    # get_manifest_path ----
    #' @description Get the path to the manifest CSV file.
    #' @return Character path.
    get_manifest_path = function() {
      private$manifest_path
    },

    # get_inputs_dir ----
    #' @description Get the path to the inputs directory.
    #' @return Character path.
    get_inputs_dir = function() {
      private$inputs_dir
    },

    # print ----
    #' @description Print a summary of the manifest.
    print = function() {
      cat("InputsManifest\n")
      cat("  Inputs dir: ", private$inputs_dir, "\n")
      cat("  Manifest:   ", private$manifest_path, "\n")
      if (is.null(self$manifest)) {
        cat("  Status:     Not loaded\n")
      } else {
        cat("  Files:      ", nrow(self$manifest), "\n")
        cat("  Categories:\n")
        counts <- self$manifest[, .N, by = category]
        data.table::setorder(counts, -N)
        for (i in seq_len(nrow(counts))) {
          cat(sprintf(
            "    %-30s %d files\n",
            counts$category[i], counts$N[i]
          ))
        }
        zenodo_n <- sum(self$manifest$zenodo_managed)
        if (zenodo_n > 0L) {
          cat("  Zenodo-managed: ", zenodo_n, " files\n")
        }
      }
      invisible(self)
    }
  ),

  private = list(

    inputs_dir = NULL,
    manifest_path = NULL,

    # to_relative ----
    # Convert absolute path to relative path from inputs_dir.
    # Uses fixed string matching (no regex) to avoid escaping issues
    # with special characters in directory paths.
    to_relative = function(abs_path) {
      prefix <- paste0(private$inputs_dir, "/")
      if (startsWith(abs_path, prefix)) {
        return(substring(abs_path, nchar(prefix) + 1L))
      }
      if (startsWith(abs_path, private$inputs_dir)) {
        out <- substring(abs_path, nchar(private$inputs_dir) + 1L)
        return(sub("^/", "", out))
      }
      abs_path
    },

    # scan_files ----
    # List all data files in inputs_dir matching patterns
    scan_files = function(
      file_patterns = c(
        "\\.parquet$", "\\.fst$", "\\.csv$", "\\.csvy$",
        "\\.qs$", "\\.yaml$"
      )
    ) {
      if (!dir.exists(private$inputs_dir)) {
        warning("Inputs directory not found: ", private$inputs_dir)
        return(character(0))
      }

      all_files <- list.files(
        private$inputs_dir,
        recursive = TRUE,
        full.names = TRUE,
        all.files = FALSE
      )

      # Filter to matching patterns
      pattern <- paste(file_patterns, collapse = "|")
      all_files <- all_files[grepl(pattern, all_files, ignore.case = TRUE)]

      # Exclude scripts and non-data directories
      all_files <- all_files[!grepl("/0-scripts/", all_files)]
      all_files <- all_files[!grepl("/1-validation/", all_files)]
      all_files <- all_files[!grepl("\\.Rhistory$", all_files)]

      all_files
    },

    # categorize_file ----
    # Determine category from relative path (first path component)
    categorize_file = function(relative_path) {
      parts <- strsplit(relative_path, "/")[[1]]
      if (length(parts) >= 2L) {
        parts[1]
      } else {
        "root"
      }
    },

    # compute_aggregate_hash ----
    # Combine multiple hashes into a single aggregate hash
    compute_aggregate_hash = function(hashes) {
      digest::digest(paste(hashes, collapse = ";"), serialize = FALSE)
    },

    # compute_file_hash ----
    compute_file_hash = function(file_path, algo = "xxhash64") {
      digest::digest(file = file_path, algo = algo)
    }
  )
)
