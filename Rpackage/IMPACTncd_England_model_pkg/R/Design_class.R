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

#' R6 Class representing a simulation design
#'
#' @description
#' The `Design` class is responsible for managing the configuration and parameters of the simulation.
#' It handles the loading of parameters from a YAML file or a list, validates them, and performs
#' necessary preprocessing such as topological sorting of diseases based on their dependencies.
#'
#' @details
#' The `Design` class ensures that the simulation is set up correctly before execution.
#' It performs the following key tasks:
#' \itemize{
#'   \item **Parameter Loading:** Reads simulation parameters from a YAML file or a list.
#'   \item **Validation:** Checks for the existence of required parameters and validates their values (e.g., non-negative numbers).
#'   \item **Environment Configuration:** Detects if the simulation is running inside a Docker container and adjusts output directories accordingly.
#'   \item **Dependency Management:** Reorders the list of diseases based on their dependencies (topological sort) to ensure correct initialization order.
#' }
#' @export
Design <-
  R6::R6Class(
    classname = "Design",

    # public ------------------------------------------------------------------
    public = list(
      #' @field sim_prm A list containing all the simulation parameters.
      sim_prm = NA,

      #' @field exposures A named list of Exposure objects for generating synthetic population exposures.
      exposures = list(),

      #' @field RR A named list of ExposureEffect objects for relative risk effects.
      RR = list(),

      #' @field diseases A named list of Disease objects.
      diseases = list(),

      #' @field data_loaded Logical indicating whether input data has been loaded.
      data_loaded = FALSE,

      #' @description
      #' Create a new `Design` object.
      #'
      #' @param sim_prm Either a path to a YAML configuration file or a list containing the simulation parameters.
      #'   The parameters must include the following keys:
      #'   \itemize{
      #'     \item `clusternumber`: Integer. Number of CPU cores/clusters to use for parallel processing.
      #'     \item `clusternumber_export`: Integer. Number of CPU cores/clusters to use for exporting summaries.
      #'     \item `logs`: Logical. If TRUE, enables verbose logging of simulation progress.
      #'     \item `export_xps`: Logical. If TRUE, exports cross-sectional population states (xps) to disk.
      #'     \item `export_PARF`: Logical. If TRUE, exports Population Attributable Risk Fraction (PARF) files.
      #'     \item `n`: Integer. Size of the synthetic population to generate or simulate.
      #'     \item `num_chunks`: Integer. Number of chunks to split the simulation into for memory management.
      #'     \item `init_year_long`: Integer. The starting year of the simulation (e.g., 2015).
      #'     \item `sim_horizon_max`: Integer. The final year of the simulation horizon.
      #'     \item `ageL`: Integer. Lower bound of the age range for the simulation.
      #'     \item `ageH`: Integer. Upper bound of the age range for the simulation.
      #'     \item `apply_RR_to_mrtl2`: Logical. If TRUE, disease mortality is influenced by exposures like incidence.
      #'     \item `calibrate_to_incd_trends`: Logical. If TRUE, use incidence calibration multipliers.
      #'     \item `calibrate_to_ftlt_trends`: Logical. If TRUE, use fatality calibration multipliers.
      #'     \item `init_year_incd_calibration`: Logical. Indicates whether parf*p0 is calibrated to incd for the initial year.
      #'     \item `incd_uncertainty_distr`: String. Distribution of incidence uncertainty ("beta" or "uniform").
      #'     \item `prvl_uncertainty_distr`: String. Distribution of prevalence uncertainty ("beta" or "uniform").
      #'     \item `ftlt_uncertainty_distr`: String. Distribution of case fatality uncertainty ("beta" or "uniform").
      #'     \item `output_dir`: String. Path to the directory where simulation outputs will be saved.
      #'     \item `synthpop_dir`: String. Path to the directory containing or saving synthetic population files.
      #'     \item `diseases`: List. Definitions of diseases, including incidence, mortality, and disability weights.
      #'     \item `maxlag`: Integer. Maximum lag period (in years) for exposure effects.
      #'     \item `jumpiness`: Numeric. Parameter controlling the volatility or "jumpiness" of time-variant trends.
      #'     \item `keep_simulants_rn`: Logical. If TRUE, keep random numbers used for exposure generation.
      #'     \item `load_simulants_rn`: Logical. If TRUE, load random numbers used for exposure generation.
      #'     \item `decision_aid`: Logical. If TRUE, enables features related to decision aid outputs.
      #'     \item `stochastic`: Logical. If TRUE, enables stochastic uncertainty in the model parameters.
      #'     \item `kismet`: Logical. If TRUE, enables "kismet" (fate/randomness) in individual life courses.
      #'     \item `max_prvl_for_outputs`: Numeric. Threshold for maximum prevalence to report in outputs (to filter rare conditions).
      #'     \item `iteration_n_max`: Integer. Maximum number of iterations allowed (safety limit).
      #'     \item `cols_for_output`: Character vector. Names of columns to include in the output files.
      #'     \item `strata_for_output`: Character vector. Variables to stratify the output by (e.g., "sex", "age_group").
      #'     \item `exposures_for_output`: List. Definitions of risk factor exposures included in the model.
      #'   }
      #'
      #' @return A new `Design` object.
      #'
      #' @examples
      #' design <- Design$new("./validation/design_for_trends_validation.yaml")
      # initialize ----
      initialize = function(sim_prm) {
        data_type <- typeof(sim_prm)
        if (data_type == "character") {
          sim_prm <- read_yaml(base::normalizePath(sim_prm, mustWork = TRUE))
        } else if (data_type != "list") {
          stop(
            "You can initialise the object only with an R object of
                     type `list` or a path to a YAML configuration file"
          )
        }

        # Validate simulation parameters
        required_params <- c(
          "locality",
          "clusternumber",
          "logs",
          "cols_for_output",
          "strata_for_output",
          "exposures_for_output",
          "exposure_definitions",
          "n",
          "init_year_long",
          "sim_horizon_max",
          "ageL",
          "ageH",
          "diseases",
          "maxlag",
          "smoking_relapse_limit",
          "stochastic",
          "kismet",
          "jumpiness",
          "statin_adherence",
          "bpmed_adherence",
          "export_xps",
          "simsmok_calibration",
          "output_dir",
          "synthpop_dir",
          "validation",
          "iteration_n_max",
          "num_chunks"
        )

        stopifnot(
          required_params %in% names(sim_prm),
          vapply(sim_prm, function(x) {
            if (is.numeric(x)) {
              all(x >= 0)
            } else {
              TRUE
            }
          }, FUN.VALUE = logical(1))
        )

        # strata_for_output columns must be present in the lifecourse parquet,
        # otherwise DuckDB summary queries fail with an opaque Binder Error.
        # The lifecourse table also receives: mc/mc_chunk (added in
        # simulate_synthpop()), scenario (Hive partition key), and any
        # disease-derived columns matching the patterns below.
        implicit_cols <- c("mc", "mc_chunk", "scenario")
        disease_col_pattern <- "^cms_|_prvl$|_dgns$|_mrtl$"
        missing_strata <- setdiff(
          sim_prm$strata_for_output,
          c(sim_prm$cols_for_output, implicit_cols)
        )
        missing_strata <- grep(
          disease_col_pattern,
          missing_strata,
          invert = TRUE,
          value = TRUE
        )
        if (length(missing_strata) > 0L) {
          stop(sprintf(
            "strata_for_output references column(s) not in cols_for_output: %s. Add them to cols_for_output (and rerun the lifecourse step), or remove them from strata_for_output.",
            paste(missing_strata, collapse = ", ")
          ))
        }

        sim_prm$sim_horizon_max <- as.integer(sim_prm$sim_horizon_max -
          sim_prm$init_year_long)
        sim_prm$init_year <- as.integer(sim_prm$init_year_long - 2000L)
        # place holders to be updated from self$update_fromGUI(parameters)

        # Default values
        sim_prm$national_qimd <- TRUE
        sim_prm$init_year_fromGUI <- sim_prm$init_year
        sim_prm$sim_horizon_fromGUI <- sim_prm$sim_horizon_max

        # Force LAD population projection if smaller localities are used
        if (!"England" %in% sim_prm$locality) {
          sim_prm$calibrate_to_pop_projections_by_LAD <- TRUE
        }

        # change output_dir & synthpop_dir if inside a docker container created by create_env.sh (or ps1)
        # But NOT if running in GitHub Actions (where workspace is mounted differently)
        if (private$is_in_docker() && Sys.getenv("GITHUB_ACTIONS") == "") {
          # if in docker
          if (sim_prm$logs) {
            message(
              "R runs within docker.\nSetting output_dir and synthpop_dir set to /outputs and /synthpop."
            )
          }
          # set the output_dir and synthpop_dir to the docker container paths
          sim_prm$output_dir <- "/outputs"
          sim_prm$synthpop_dir <- "/synthpop"
        } else {
          # if not in docker
          sim_prm$output_dir <- normalizePath(
            sim_prm$output_dir,
            mustWork = FALSE
          )
          sim_prm$synthpop_dir <- normalizePath(
            sim_prm$synthpop_dir,
            mustWork = FALSE
          )
        }

        # Reorder the diseases based on dependencies
        sim_prm$diseases <- private$reorder_diseases(sim_prm)

        # Store the simulation parameters
        self$sim_prm = sim_prm

        # Load data if input files are present; otherwise create
        # a skeleton object that supports zenodo_download_all()
        if (self$has_input_data()) {
          self$load_data()
        } else {
          message(
            "Input data not found. Object created in skeleton mode.\n",
            "Download data with: $zenodo_connect() then $zenodo_download_all()"
          )
        }

        invisible(self)
      },

      #' @description
      #' Load exposure effects from RR files
      #' @return The invisible self for chaining.
      # load_exposure_effects ----
      load_exposure_effects = function() {
        message("Loading exposure effects.")
        # Create a named list of ExposureEffect objects for the files in ./inputs/RR
        fl <- list.files(
          path = "./inputs/RR",
          pattern = ".csvy$",
          full.names = TRUE
        )
        self$RR <- lapply(fl, ExposureEffect$new, design = self)
        names(self$RR) <- vapply(self$RR, function(x) x$get_name(), FUN.VALUE = character(1))

        # Generate stochastic effects sequentially to avoid registry race conditions
        # This is a one-time initialization cost, and the files are cached afterward
        invisible(lapply(self$RR, function(x) {
          x$gen_stochastic_effect(self, overwrite = FALSE, smooth = FALSE)
        }))
        # NOTE smooth cannot be exported to Design for now, because the first
        # time this parameter changes we need logic to overwrite unsmoothed file

        # Clean up orphaned RR .fst files not belonging to
        # any current ExposureEffect object
        valid_files <- unlist(lapply(
          self$RR, function(rr_obj) {
            suffix <- paste0(
              rr_obj$name, "~", rr_obj$outcome
            )
            c(
              file.path(
                "./simulation/rr",
                paste0("rr_", suffix, "_l.fst")
              ),
              file.path(
                "./simulation/rr",
                paste0("rr_", suffix, "_indx.fst")
              )
            )
          }
        ))
        if (is.null(valid_files)) valid_files <- character(0)
        valid_files <- normalizePath(
          valid_files, mustWork = FALSE
        )

        rr_dir <- "./simulation/rr"
        if (dir.exists(rr_dir)) {
          all_physical_files <- normalizePath(
            list.files(
              rr_dir, pattern = "\\.fst$",
              full.names = TRUE
            ),
            mustWork = FALSE
          )
          orphaned_files <- setdiff(
            all_physical_files, valid_files
          )

          if (length(orphaned_files) > 0) {
            if (self$sim_prm$logs) {
              message(
                "Removing ",
                length(orphaned_files),
                " orphaned RR files."
              )
            }
            file.remove(orphaned_files)
          }
        }

        # Clean up generation record to remove entries
        # for files that no longer exist
        rec_path <- file.path(
          getwd(), "simulation",
          "rr_generation_record.csv"
        )
        if (file.exists(rec_path)) {
          rec <- fread(rec_path)
          if (
            nrow(rec) > 0 &&
              "file_path" %in% names(rec)
          ) {
            rec <- rec[
              normalizePath(
                file_path, mustWork = FALSE
              ) %in% valid_files
            ]
            # Ensure paths are relative for portability
            # (handles prefixes from any user/machine)
            rec[, file_path := sub(
              "^.*/(?=simulation/)", "", file_path,
              perl = TRUE
            )]
            if ("source_file" %in% names(rec)) {
              rec[, source_file := sub(
                "^.*/(?=inputs/)", "", source_file,
                perl = TRUE
              )]
            }
            fwrite(rec, rec_path)
          }
        }

        # Remove legacy fileversion.csv if it exists
        # (replaced by inputs_manifest.csv and
        # *_generation_record.csv)
        legacy_registry <- file.path(
          getwd(), "simulation", "fileversion.csv"
        )
        if (file.exists(legacy_registry)) {
          if (self$sim_prm$logs) {
            message(
              "Removing legacy fileversion.csv ",
              "(replaced by inputs manifest system)."
            )
          }
          file.remove(legacy_registry)
        }

        invisible(self)
      },

      #' @description
      #' Create Exposure objects from YAML configuration
      #' @return The invisible self for chaining.
      # load_exposures ----
      load_exposures = function() {
        # Check if exposure data files exist (not just the directory, which
        # may contain only scripts). The first exposure definition's data
        # directory serves as a proxy for all data availability.
        first_exp <- self$sim_prm$exposure_definitions[[1]]
        first_path <- file.path("./inputs/exposure_distributions", first_exp$file_name)
        if (!dir.exists(first_path) && !file.exists(first_path)) {
          warning(
            "Exposure data not found (", first_path, ").\n",
            "Skipping exposure loading. Download inputs first via zenodo_download_all().",
            call. = FALSE
          )
          return(invisible(self))
        }
        message("Creating exposure generators from configuration.")
        self$exposures <- lapply(
          self$sim_prm$exposure_definitions,
          function(exp_def) {
            # Convert YAML configuration to Exposure$new() arguments
            args <- list(
              name = exp_def$name,
              file_name = exp_def$file_name,
              var_name = exp_def$var_name,
              rank_var = exp_def$rank_var,
              distribution = exp_def$distribution
            )

            # Add optional parameters if present
            if (!is.null(exp_def$qfun)) {
              args$qfun <- exp_def$qfun
            }
            if (!is.null(exp_def$pfun)) {
              args$pfun <- exp_def$pfun
            }
            if (!is.null(exp_def$min_value)) {
              args$min_value <- exp_def$min_value
            }
            if (!is.null(exp_def$max_value)) {
              args$max_value <- exp_def$max_value
            }
            if (!is.null(exp_def$qparams)) {
              args$qparams <- exp_def$qparams
            }
            if (!is.null(exp_def$qargs)) {
              args$qargs <- exp_def$qargs
            }
            if (!is.null(exp_def$thresholds)) {
              args$thresholds <- exp_def$thresholds
            }
            if (!is.null(exp_def$invert)) {
              args$invert <- exp_def$invert
            }
            if (!is.null(exp_def$lower_tail)) {
              args$lower_tail <- exp_def$lower_tail
            }
            if (!is.null(exp_def$transform_fn)) {
              args$transform_fn <- exp_def$transform_fn
            }
            if (!is.null(exp_def$offset)) {
              args$offset <- exp_def$offset
            }
            if (!is.null(exp_def$additional_rank_vars)) {
              args$additional_rank_vars <- exp_def$additional_rank_vars
            }
            if (!is.null(exp_def$join_fn)) {
              args$join_fn <- exp_def$join_fn
            }

            # Create the Exposure object
            exp_obj <- do.call(Exposure$new, args)

            # Store factorisation parameters if present (to be applied after generation)
            if (!is.null(exp_def$factorise)) {
              exp_obj$factorise_params <- exp_def$factorise
            }

            # Store additional metadata
            if (!is.null(exp_def$cast_to_int)) {
              exp_obj$cast_to_int <- exp_def$cast_to_int
            }
            if (!is.null(exp_def$invert_ratio)) {
              exp_obj$invert_ratio <- exp_def$invert_ratio
            }

            exp_obj
          }
        )
        names(self$exposures) <- vapply(self$exposures, function(x) x$name, FUN.VALUE = character(1))

        # Reorder exposures based on dependencies (topological sort)
        self$exposures <- private$reorder_exposures(self$exposures)

        invisible(self)
      },

      #' @description
      #' Create Disease objects from YAML configuration
      #' @return The invisible self for chaining.
      # load_diseases ----
      load_diseases = function() {
        message("Loading diseases.")
        self$diseases <- lapply(self$sim_prm$diseases, function(x) {
          x[["design_"]] <- self
          x[["RR"]] <- self$RR
          do.call(Disease$new, x)
        })
        names(self$diseases) <- vapply(self$sim_prm$diseases, `[[`, "name", FUN.VALUE = character(1))

        invisible(self)
      },

      #' @description
      #' Check whether input data files are present.
      #' @return Logical; TRUE if disease burden directories exist in
      #'   `./inputs/`.
      # has_input_data ----
      has_input_data = function() {
        dirs <- list.dirs("./inputs", recursive = FALSE, full.names = FALSE)
        any(grepl("^disease_burden", dirs))
      },

      #' @description
      #' Load all input data (exposure effects, exposures, diseases).
      #'
      #' Called automatically by the constructor when data is present.
      #' Also called by `Simulation$zenodo_download_all()` after
      #' downloading data, to complete initialization without needing
      #' to create a new object.
      #'
      #' @param force Logical; if TRUE, reload even if already loaded.
      #' @return The invisible self for chaining.
      # load_data ----
      load_data = function(force = FALSE) {
        if (self$data_loaded && !force) {
          if (self$sim_prm$logs) {
            message("Data already loaded. Use force = TRUE to reload.")
          }
          return(invisible(self))
        }
        if (!self$has_input_data()) {
          stop(
            "Input data not found in ./inputs/.\n",
            "Download inputs first via zenodo_download_all().",
            call. = FALSE
          )
        }
        self$load_exposure_effects()
        self$load_exposures()
        self$load_diseases()
        self$data_loaded <- TRUE
        invisible(self)
      },

      #' Save Simulation Design to Disk
      #' @description
      #' Serializes the current simulation parameters (`sim_prm`) and writes them to a YAML file.
      #' This can be useful for saving configurations that can later be reloaded or shared.
      #'
      #' @param path A character string specifying the full file path (including filename and extension)
      #' where the simulation design should be saved.
      #'
      #' @return The `Design` object, invisibly.
      #'
      #' @examples
      #' design <- Design$new("inputs/sim_design.yaml")
      #' design$save_to_disk("outputs/saved_sim_design.yaml")
      #'
      #' @importFrom yaml write_yaml
      # save_to_disk ----
      save_to_disk = function(path) {
        write_yaml(self$sim_prm, base::normalizePath(path, mustWork = FALSE))

        invisible(self)
      },

      #' @description
      #' Print the simulation parameters.
      #' @return The `Design` object, invisibly.
      # print ----
      print = function() {
        print(self$sim_prm)

        invisible(self)
      }
    ),

    # private ------------------------------------------------------------------
    private = list(
      mc_aggr = NA,

      # Reorder Exposures by Dependency
      # @description
      # Reorders the list of exposures using topological sorting based on column dependencies.
      # An exposure A depends on exposure B if B's var_name appears in A's parquet column names.
      # This ensures exposures are generated in the correct order.
      #
      # @param exposures A named list of Exposure objects.
      #
      # @return A reordered list of exposures in topological order.
      #
      # @keywords internal
      # reorder_exposures ----
      reorder_exposures = function(exposures) {
        # Get all exposure variable names
        exp_var_names <- names(exposures)

        # Build dependency graph edges
        # An edge from B -> A means A depends on B (B must come before A)
        edge_list <- lapply(exp_var_names, function(exp_name) {
          exp_obj <- exposures[[exp_name]]
          col_names <- exp_obj$get_col_names()

          # Check which other exposure var_names appear in this exposure's columns
          dependencies <- setdiff(intersect(col_names, exp_var_names), exp_name)

          if (length(dependencies) > 0) {
            # Return list with edges and messages
            list(
              edges = as.vector(rbind(dependencies, exp_name)),
              messages = paste0(exp_name, " depends on ", dependencies)
            )
          } else {
            list(edges = character(0), messages = character(0))
          }
        })
        edges <- do.call(c, lapply(edge_list, `[[`, "edges"))
        dep_messages <- do.call(c, lapply(edge_list, `[[`, "messages"))

        # If no dependencies, return as-is
        if (length(edges) == 0) {
          return(exposures)
        }

        # Log detected dependencies
        if (length(dep_messages) > 0) {
          message(
            "Exposure dependencies detected: ",
            paste(dep_messages, collapse = ", ")
          )
        }

        # Create directed graph
        g <- make_graph(edges, directed = TRUE)

        # Check for cycles and report them
        Cycles <- NULL
        for (v1 in V(g)) {
          for (v2 in neighbors(g, v1, mode = "out")) {
            Cycles <- c(
              Cycles,
              lapply(
                all_simple_paths(g, v2, v1, mode = "out"),
                function(p) c(v1, p)
              )
            )
          }
        }
        # Remove duplicates
        if (length(Cycles) > 0) {
          cycle_mins <- vapply(Cycles, min, FUN.VALUE = integer(1))
          cycle_firsts <- vapply(Cycles, `[`, 1, FUN.VALUE = integer(1))
          Cycles <- Cycles[cycle_mins == cycle_firsts]
          cycle_lengths <- vapply(Cycles, length, FUN.VALUE = integer(1))
          cycles_of_interest <- Cycles[which(cycle_lengths >= 3)]
          if (length(cycles_of_interest) > 0) {
            cycle_strs <- vapply(cycles_of_interest, function(c) {
              paste(names(c), collapse = " -> ")
            }, FUN.VALUE = character(1))
            warning(
              "Circular dependencies found in exposures: ",
              paste(cycle_strs, collapse = "; ")
            )
          }
        }

        # Perform topological sort
        o <- topo_sort(g)
        sorted_names <- names(o)

        # Add any exposures that weren't in the graph (no dependencies)
        missing_names <- setdiff(exp_var_names, sorted_names)
        final_order <- c(missing_names, sorted_names)

        # Log the reordering
        message("Exposures reordered: ", paste(final_order, collapse = " -> "))

        # Reorder exposures
        exposures[final_order]
      },

      # Reorder Diseases by Dependency
      # @description
      # Reorders the list of diseases in the simulation parameters using topological sorting.
      # This ensures that diseases are processed in the correct order based on their causal dependencies,
      # such that diseases influencing others appear earlier in the list.
      #
      # This is critical for simulations that depend on proper initialization and evaluation
      # of disease incidence affected by the presence of other diseases.
      #
      # @param sim_prm A list of simulation parameters. Must contain a named list of diseases
      # where each disease has a 'name' field and may have dependencies listed under
      # meta$incidence$influenced_by_disease_name.
      #
      # @return A reordered list of diseases in topological order.
      #
      # @importFrom igraph make_graph topo_sort
      # @keywords internal
      # reorder_diseases ----
      reorder_diseases = function(sim_prm) {
        # Reorder the diseases so dependencies are always calculated first
        # (topological ordering). This is crucial for init_prevalence
        # first name the list and
        sim_prm$diseases <-
          setNames(
            sim_prm$diseases,
            vapply(sim_prm$diseases, function(x) x$name, FUN.VALUE = character(1))
          )

        # Build graph structure using lapply to avoid vector growth
        ds <- names(sim_prm$diseases)
        edge_list <- lapply(seq_along(ds), function(i) {
          dep <- sim_prm[["diseases"]][[i]][["meta"]][["incidence"]][[
            "influenced_by_disease_name"
          ]]
          if (length(dep) > 0L) {
            # Create edge pairs: dep -> ds[i] for each dependency
            as.vector(rbind(dep, ds[i]))
          } else {
            character(0)
          }
        })
        out <- do.call(c, edge_list)

        # Handle case where no dependencies exist
        if (length(out) == 0) {
          return(sim_prm$diseases)
        }

        g <- make_graph(out, directed = TRUE)
        # stopifnot(is_dag(g)) # Removed for diseases depend on diseases
        # get all cycles in the graph
        Cycles <- NULL
        for (v1 in V(g)) {
          for (v2 in neighbors(g, v1, mode = "out")) {
            Cycles <- c(
              Cycles,
              lapply(
                all_simple_paths(g, v2, v1, mode = "out"),
                function(p) c(v1, p)
              )
            )
          }
        }
        # remove duplicates
        if (length(Cycles) > 0) {
          cycle_mins <- vapply(Cycles, min, FUN.VALUE = integer(1))
          cycle_firsts <- vapply(Cycles, `[`, 1, FUN.VALUE = integer(1))
          Cycles <- Cycles[cycle_mins == cycle_firsts]
          # find cycles of length i.e. 3 (i.e. chd -> t2dm -> chd)
          cycle_lengths <- vapply(Cycles, length, FUN.VALUE = integer(1))
          Cycles[which(cycle_lengths >= 3)]
        }

        if (sim_prm$logs && length(Cycles) > 0) {
          message("Cycles found: ", Cycles)
        }

        o <- topo_sort(g)

        # then reorder based on the topological ordering
        sim_prm$diseases[order(match(names(sim_prm$diseases), names(o)))]
      },

      # @description
      # Check whether the R session is running inside a Docker container.
      #
      # @details
      # This method detects whether the current R session is running in a Docker
      # container by checking for the presence of the special file `/.dockerenv` and
      # examining system-specific files for identifiers associated with
      # Docker or Kubernetes. The function is cross-platform compatible.
      #
      # @return A logical value: `TRUE` if inside a Docker container, otherwise `FALSE`.
      # is_in_docker ----
      is_in_docker = function() {
        # Check for the standard Docker environment file (works on all platforms)
        if (file.exists("/.dockerenv")) {
          return(TRUE)
        }

        # Platform-specific checks
        os_type <- Sys.info()[["sysname"]]

        if (os_type == "Linux") {
          # Linux: Check /proc/1/cgroup for docker/kubepods
          cgroup_file <- "/proc/1/cgroup"
          if (file.exists(cgroup_file)) {
            tryCatch(
              {
                cgroup_content <- readLines(cgroup_file, warn = FALSE)
                return(any(grepl("docker|kubepods", cgroup_content)))
              },
              error = function(e) {
                return(FALSE)
              }
            )
          }
        } else if (os_type == "Windows") {
          # Windows: Check for Docker-specific environment variables
          docker_vars <- c(
            "DOCKER_CONTAINER",
            "DOCKER_HOST",
            "DOCKER_MACHINE_NAME"
          )
          has_docker_var <- vapply(
            docker_vars,
            function(x) Sys.getenv(x) != "",
            FUN.VALUE = logical(1)
          )
          if (any(has_docker_var)) {
            return(TRUE)
          }

          # Check for Windows container indicators
          tryCatch(
            {
              # Check if running in Windows container by looking for container-specific registry
              system_output <- suppressWarnings(system(
                "reg query HKLM\\SYSTEM\\CurrentControlSet\\Control\\ContainerManager",
                intern = TRUE,
                ignore.stderr = TRUE
              ))
              return(
                length(system_output) > 0 && !any(grepl("ERROR", system_output))
              )
            },
            error = function(e) {
              return(FALSE)
            }
          )
        } else if (os_type == "Darwin") {
          # macOS: Check for Docker-specific environment variables and processes
          docker_vars <- c(
            "DOCKER_CONTAINER",
            "DOCKER_HOST",
            "DOCKER_MACHINE_NAME"
          )
          has_docker_var <- vapply(
            docker_vars,
            function(x) Sys.getenv(x) != "",
            FUN.VALUE = logical(1)
          )
          if (any(has_docker_var)) {
            return(TRUE)
          }

          # Check for macOS-specific container indicators
          tryCatch(
            {
              # Check if we're in a container by examining process hierarchy
              ps_output <- system(
                "ps -p 1 -o comm=",
                intern = TRUE,
                ignore.stderr = TRUE
              )
              return(
                length(ps_output) > 0 &&
                  any(grepl("docker|container", ps_output, ignore.case = TRUE))
              )
            },
            error = function(e) {
              return(FALSE)
            }
          )
        }

        # Fallback: Check common environment variables that might indicate Docker
        docker_env_vars <- c(
          "DOCKER_CONTAINER",
          "CONTAINER",
          "KUBERNETES_SERVICE_HOST"
        )
        has_env_var <- vapply(
          docker_env_vars,
          function(x) Sys.getenv(x) != "",
          FUN.VALUE = logical(1)
        )
        return(any(has_env_var))
      }
    ) # end of private
  ) # end of R6Class
