## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by
## Chris Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2025 University of Liverpool, Chris Kypridemos
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

#' @title Design Class
#' @name Design
#'
#' @description
#' An R6 class representing a simulation design for the IMPACTncd framework.
#' This class encapsulates simulation configuration parameters (`sim_prm`),
#' enforces validation, supports GUI updates, and provides tools for detecting and resolving
#' dependency structures among diseases.
#'
#' @details
#' The `Design` class is used to initialize and manage the structure of a simulation run.
#' It can load simulation parameters from a YAML configuration file or directly from an R list.
#' The class performs validation of required fields, sets default values, and computes derived fields
#' like initial year and simulation horizon. It also handles the topological ordering of diseases
#' based on their dependencies and can detect cyclic relationships using `igraph`.
#'
#' @section Initialization:
#' The class is initialized with either a file path (YAML) or a list:
#'
#' ```r
#' # From YAML file
#' design <- Design$new("inputs/sim_design.yaml")
#' 
#' # From list
#' config_list <- list(
#'   simulation_files_overwrite = FALSE,
#'   sTag = "my_simulation",
#'   # ... other required fields
#' )
#' design <- Design$new(config_list)
#' ```
#'
#' @section Public Fields:
#' \describe{
#'   \item{\code{sim_prm}}{A named list of validated and processed simulation parameters.}
#' }
#'
#' @section Public Methods:
#' \describe{
#'   \item{\code{initialize(sim_prm)}}{Initializes the object. Accepts a list or path to YAML file.}
#'   \item{\code{save_to_disk(path)}}{Saves the current configuration to disk as a YAML file.}
#'   \item{\code{print()}}{Prints the current design parameters.}
#' }
#'
#' @section Private Methods:
#' \describe{
#'   \item{\code{reorder_diseases(sim_prm)}}{Sorts diseases in topological order based on their dependencies.}
#'   \item{\code{is_in_docker()}}{Checks if R session is running inside a Docker container.}
#'   \item{\code{update_fromGUI(GUI_prm)}}{Updates the design based on GUI-provided inputs.}
#' }
#'
#' @section Required Configuration Fields:
#' The following fields are required in the simulation parameters:
#' \describe{
#'   \item{\code{simulation_files_overwrite}}{Whether to overwrite existing simulation files.}
#'   \item{\code{sTag}}{Simulation tag identifier.}
#'   \item{\code{bOverwriteFilesOnDeploy}}{Whether to overwrite files on deployment.}
#'   \item{\code{RootDirPath}}{Root directory path for the simulation.}
#'   \item{\code{sToken}}{Security token for the simulation.}
#'   \item{\code{locality}}{Geographic locality for the simulation.}
#'   \item{\code{clusternumber}}{Cluster number for parallel processing.}
#'   \item{\code{logs}}{Whether to enable logging.}
#'   \item{\code{scenarios}}{List of scenarios to simulate.}
#'   \item{\code{cols_for_output}}{Columns to include in output.}
#'   \item{\code{strata_for_output}}{Stratification variables for output.}
#'   \item{\code{exposures}}{List of exposure variables.}
#'   \item{\code{n}}{Sample size for simulation.}
#'   \item{\code{init_year_long}}{Initial year for simulation (full year format).}
#'   \item{\code{sim_horizon_max}}{Maximum simulation horizon year.}
#'   \item{\code{ageL}}{Lower age limit.}
#'   \item{\code{ageH}}{Upper age limit.}
#'   \item{\code{diseases}}{List of diseases to model.}
#'   \item{\code{maxlag}}{Maximum lag for disease models.}
#'   \item{\code{smoking_relapse_limit}}{Limit for smoking relapse modeling.}
#'   \item{\code{stochastic}}{Whether to use stochastic modeling.}
#'   \item{\code{kismet}}{Random seed or fate parameter.}
#'   \item{\code{jumpiness}}{Parameter controlling model jumpiness.}
#'   \item{\code{statin_adherence}}{Statin adherence parameter.}
#'   \item{\code{bpmed_adherence}}{Blood pressure medication adherence.}
#'   \item{\code{decision_aid}}{Decision aid parameter.}
#'   \item{\code{export_xps}}{Whether to export exposures.}
#'   \item{\code{simsmok_calibration}}{Smoking simulation calibration parameter.}
#'   \item{\code{output_dir}}{Output directory path.}
#'   \item{\code{synthpop_dir}}{Synthetic population directory path.}
#'   \item{\code{validation}}{Validation parameters.}
#'   \item{\code{iteration_n_max}}{Maximum number of iterations.}
#'   \item{\code{n_synthpop_aggregation}}{Number for synthetic population aggregation.}
#' }
#'
#' @section Docker Support:
#' The class automatically detects if running inside a Docker container and adjusts
#' output and synthpop directory paths accordingly. When in Docker, paths are set to:
#' \describe{
#'   \item{\code{output_dir}}{"/IMPACTncd_Engl/output"}
#'   \item{\code{synthpop_dir}}{"/IMPACTncd_Engl/synthpop"}
#' }
#'
#' @section Disease Dependencies:
#' The class handles disease dependencies by:
#' \itemize{
#'   \item Detecting cycles in disease dependency graphs
#'   \item Performing topological sorting to ensure proper disease ordering
#'   \item Logging any cycles found for debugging purposes
#' }
#'
#' @section Automatic Adjustments:
#' The class automatically performs several adjustments during initialization:
#' \itemize{
#'   \item Converts `sim_horizon_max` to relative horizon by subtracting `init_year_long`
#'   \item Converts `init_year_long` to short format by subtracting 2000
#'   \item Sets default values for GUI-related parameters
#'   \item Forces LAD population projection for non-England localities
#'   \item Normalizes directory paths (except when in Docker)
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom yaml read_yaml write_yaml
#' @importFrom igraph make_graph V neighbors all_simple_paths topo_sort
#'
#' @seealso \code{\link[yaml]{read_yaml}}, \code{\link[igraph]{make_graph}}
#'
#' @examples
#' \dontrun{
#' # Initialize from YAML file
#' design <- Design$new("inputs/sim_design.yaml")
#' 
#' # Initialize from list
#' config_list <- list(
#'   simulation_files_overwrite = FALSE,
#'   sTag = "my_simulation",
#'   # ... other required fields
#' )
#' design <- Design$new(config_list)
#' 
#' # Save configuration
#' design$save_to_disk("outputs/my_config.yaml")
#' 
#' # Print current configuration
#' design$print()
#' }
#'
#' @export
Design <-
  R6::R6Class(
    classname = "Design",

    # public ------------------------------------------------------------------
    public = list(
      #' @field sim_prm The simulation parameters.
      sim_prm = NA,

      #' @description Create a new design object.
      #' @param sim_prm Either a path to a YAML file or a list with
      #'   appropriate format containing simulation parameters.
      #' @return A new `Design` object.
      #' @examples
      #' \dontrun{
      #' design <- Design$new("inputs/sim_design.yaml")
      #' }
      # initialize ----
      initialize = function(sim_prm) {
        # Determine the type of the input parameter
        data_type <- typeof(sim_prm)

        # Read YAML file if the input is a character string
        if (data_type == "character") {
          sim_prm <- read_yaml(base::normalizePath(sim_prm, mustWork = TRUE))
        } else if (data_type != "list") {
          stop(
            "You can initialise the object only with an R object of
                     type `list` or a path to a YAML configuration file"
          )
        }

        # Validate simulation parameters
        required_fields <- c(
          "simulation_files_overwrite",
          "sTag",
          "bOverwriteFilesOnDeploy",
          "RootDirPath",
          "sToken",
          "locality",
          "clusternumber",
          "logs",
          "scenarios",
          "cols_for_output",
          "strata_for_output",
          "exposures",
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
          "decision_aid",
          "export_xps",
          "simsmok_calibration",
          "output_dir",
          "synthpop_dir",
          "validation",
          "iteration_n_max",
          "n_synthpop_aggregation"
        )

        stopifnot(
          required_fields %in% names(sim_prm),
          sapply(sim_prm, function(x) if (is.numeric(x)) x >= 0 else TRUE)
        )

        # Setup simulation parameters with default values and adjustments
        sim_prm$sim_horizon_max <- sim_prm$sim_horizon_max -
          sim_prm$init_year_long
        sim_prm$init_year <- sim_prm$init_year_long - 2000L

        # Default values
        sim_prm$national_qimd <- TRUE
        sim_prm$init_year_fromGUI <- sim_prm$init_year
        sim_prm$sim_horizon_fromGUI <- sim_prm$sim_horizon_max

        # Force LAD population projection if smaller localities are used
        if (!"England" %in% sim_prm$locality) {
          sim_prm$calibrate_to_pop_projections_by_LAD <- TRUE
        }

        # Normalize output and directory path
        # change output_dir & synthpop_dir if inside a docker container created by create_env.sh (or ps1)
        if (private$is_in_docker()) {
          # if in docker
          if (sim_prm$logs) {
            message(
              "R runs within docker./nSetting output_dir and synthpop_dir to /IMPACTncd_Engl/output and /IMPACTncd_Engl/synthpop."
            )
          }
          # set the output_dir and synthpop_dir to the docker container paths
          sim_prm$output_dir <- "/IMPACTncd_Engl/output"
          sim_prm$synthpop_dir <- "/IMPACTncd_Engl/synthpop"
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
            sapply(sim_prm$diseases, function(x) {
              x$name
            })
          )

        out <- vector() # will hold graph structure
        ds <- names(sim_prm$diseases)
        for (i in seq_along(ds)) {
          ds_ <- ds[i]
          dep <-
            sim_prm[["diseases"]][[i]][["meta"]][["incidence"]][[
              "influenced_by_disease_name"
            ]]
          if (length(dep) > 0L) {
            # dep <- gsub("_prvl", "", dep)
            for (j in seq_along(dep)) {
              out <- c(out, dep[[j]], ds_)
            }
          }
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
        Cycles <- Cycles[sapply(Cycles, min) == sapply(Cycles, `[`, 1)]
        # find cycles of length i.e. 3 (i.e. chd -> t2dm -> chd)
        Cycles[which(sapply(Cycles, length) >= 3)]

        if (sim_prm$logs && length(Cycles) > 0)
          message("Cycles found: ", Cycles)

        o <- topo_sort(g)

        # then reorder based on the topological ordering
        sim_prm$diseases[order(match(names(sim_prm$diseases), names(o)))]
      },

      # Check if running in Docker
      # @description
      # Check whether the R session is running inside a Docker container.
      #
      # @details
      # This method detects whether the current R session is running in a Docker
      # container by checking for the presence of the special file `/.dockerenv` and
      # examining the `/proc/1/cgroup` system file for identifiers associated with
      # Docker or Kubernetes.
      #
      # @return A logical value: `TRUE` if inside a Docker container, otherwise `FALSE`.
      #
      # @keywords internal
      is_in_docker = function() {
        file.exists("/.dockerenv") ||
          any(grepl(
            "docker|kubepods",
            readLines("/proc/1/cgroup", warn = FALSE)
          ))
      },

      # Update Design from GUI
      # @description Updates the design object from GUI.
      # @param GUI_prm A GUI parameter object.
      # @return The `Design` object, invisibly.
      # update_fromGUI ----
      update_fromGUI = function(GUI_prm) {
        self$sim_prm$national_qimd <- GUI_prm$national_qimd_checkbox
        # T = use national qimd, F = use local qimd
        self$sim_prm$init_year_fromGUI <-
          fromGUI_timeframe(GUI_prm)["init year"] - 2000L
        self$sim_prm$sim_horizon_fromGUI <-
          fromGUI_timeframe(GUI_prm)["horizon"]
        self$sim_prm$locality <- GUI_prm$locality_select

        if (!GUI_prm$national_qimd_checkbox) {
          self$sim_prm$cols_for_output <-
            c(setdiff(self$sim_prm$cols_for_output, "lqimd"), "nqimd")
        }

        self$sim_prm$iteration_n <- GUI_prm$iteration_n_gui
        self$sim_prm$iteration_n_final <- GUI_prm$iteration_n_final_gui
        self$sim_prm$n <- GUI_prm$n_gui
        self$sim_prm$n_synthpop_aggregation <- GUI_prm$n_synthpop_aggregation_gui
        self$sim_prm$n_primers <- GUI_prm$n_primers_gui
        self$sim_prm$cancer_cure <- GUI_prm$cancer_cure_gui
        self$sim_prm$jumpiness <- GUI_prm$jumpiness_gui
        self$sim_prm$statin_adherence <- GUI_prm$statin_adherence_gui
        self$sim_prm$bpmed_adherence <- GUI_prm$bpmed_adherence_gui
        self$sim_prm$decision_aid <- GUI_prm$decision_aid_gui
        self$sim_prm$logs <- GUI_prm$logs_gui

        invisible(self)
      }
    ) # end of private
  ) # end of R6Class
