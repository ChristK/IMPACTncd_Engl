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



# Check whether the current session can open an interactive graphics device.
# Returns TRUE if a device is already open, or if the platform can open one
# (RStudio plot pane, Windows device, or X11 with a valid DISPLAY).
.has_plot_device <- function() {
  # A file-based device (png, pdf, etc.) already open counts
  if (dev.cur() > 1L) return(TRUE)
  # RStudio has its own plot pane
  if (identical(.Platform$GUI, "RStudio")) return(TRUE)
  # Windows can always open a device
  if (.Platform$OS.type == "windows") return(TRUE)
  # On Unix, X11 needs a DISPLAY
  if (capabilities("X11") && nzchar(Sys.getenv("DISPLAY"))) return(TRUE)
  FALSE
}

# From
# https://stackoverflow.com/questions/33424233/how-do-i-tell-an-r6-class-what-to-do-with-square-brackets
# Allows data.table syntax to the R6class object directly. Assumes it has a
# field 'output' that is a data.table

#' @export
`[.Simulation` <- function(x, ...) x$output[...]

#' R6 Class representing a simulation environment
#' @description A simulation environment for running microsimulation models
#'   of chronic disease epidemiology and health economics in England.
#'
#' @details The Simulation class orchestrates the microsimulation, managing
#'   the synthetic population, disease modules, exposure distributions, and
#'   policy scenarios. It provides methods for running simulations, calibrating
#'   disease parameters, and exporting summary statistics including QALYs and costs.
#'
#' @section Implementation note:
#' For maintainability the public methods `calibrate_incd_ftlt()`,
#' `export_summaries()`, and `export_tables()` (and their private helpers)
#' are defined in satellite files (`Simulation_class_calibration.R`,
#' `Simulation_class_summaries.R`, `Simulation_class_tables.R`) and added
#' to the class via `Simulation$set()`. They appear below in the "Methods"
#' section like any directly-defined method.
#'
#' @section QALY Calculation:
#' QALYs are calculated using EQ-5D-5L utility values based on:
#' \itemize{
#'   \item Janssen & Szende (2014) population utility norms by single year of age
#'     \url{https://link.springer.com/book/10.1007/978-94-007-7596-1} (Table 3.7)
#'   \item Sullivan et al. (2011) utility decrements for chronic conditions
#'     \url{https://doi.org/10.1177/0272989X11401031} (Supplementary Tables 4 & 5)
#' }
#'
#' The EQ-5D-5L utility is computed as:
#' \deqn{EQ5D = PopNorm(age) + Income + Education + NCC + Sex + Ethnicity - \sum Diseases}
#'
#' Components:
#' \describe{
#'   \item{Population norms}{Age-specific baseline utilities (1.0 for ages 0-15,
#'     decreasing to 0.706 for ages 75+)}
#'   \item{Income}{Utility adjustment by income quintile (0 to +0.041)}
#'   \item{Education}{Utility adjustment by education level (0 to +0.007)}
#'   \item{NCC}{Number of chronic conditions adjustment based on cms_count (clamped 0-10)}
#'   \item{Sex}{+0.001 for men}
#'   \item{Ethnicity}{Decrements for non-white ethnicities (-0.00009 to -0.0015)}
#'   \item{Disease decrements}{27 condition-specific utility decrements including:
#'     T2DM (-0.071), CHD (-0.067), stroke (-0.106), dementia (-0.217),
#'     pain (-0.341), and others}
#' }
#'
#' For individuals who die during the year, the QALY is halved. Final values
#' are clamped to the range `[0, 1]`.
#'
#' Required columns in lifecourse data: \code{age}, \code{sex}, \code{income},
#' \code{education}, \code{ethnicity}, \code{cms_count}, \code{all_cause_mrtl},
#' and disease prevalence columns (e.g., \code{chd_prvl}, \code{stroke_prvl}).
#'
#' @section Cost Calculation:
#' Costs and economic output are calculated using the following components:
#' \describe{
#'   \item{healthcare_cost}{Direct healthcare costs based on disease prevalence,
#'     inflation-adjusted from 2019 values. Higher values = more healthcare spending.}
#'   \item{socialcare_cost}{Formal social care costs by age and disease status.
#'     Higher values = more social care spending.}
#'   \item{informalcare_cost}{Costs of unpaid care by family/friends, estimated
#'     using a two-stage regression model based on health utility and comorbidities.
#'     Higher values = greater informal care burden.}
#'   \item{economic_output}{Monetary value of economic production (paid and unpaid work),
#'     adjusted for health-related productivity based on EQ-5D utility.
#'     Higher values = more economic production (this is a benefit, not a cost).}
#' }
#'
#' Aggregated metrics (from a societal perspective):
#' \itemize{
#'   \item indirect_cost = socialcare + informalcare - economic_output
#'   \item total_cost = healthcare + socialcare + informalcare - economic_output
#' }
#'
#' @section Interpreting Net Costs:
#' Net costs are calculated as: intervention scenario minus baseline scenario.
#' \describe{
#'   \item{healthcare_cost, socialcare_cost, informalcare_cost}{
#'     Negative net value = cost savings (intervention reduces spending).
#'     Positive net value = additional costs (intervention increases spending).}
#'   \item{economic_output}{
#'     Positive net value = productivity gains (intervention increases economic production).
#'     Negative net value = productivity losses (intervention reduces economic production).}
#'   \item{indirect_cost, total_cost}{
#'     Negative net value = net societal benefit (savings exceed any lost production,
#'     or savings plus productivity gains).
#'     Positive net value = net societal cost (costs exceed benefits).}
#' }
#'
#' For a beneficial public health intervention, you would typically expect:
#' \itemize{
#'   \item Negative net healthcare/socialcare/informalcare costs (savings)
#'   \item Positive net economic_output (productivity gains from healthier population)
#'   \item Negative net total_cost (overall societal benefit)
#' }
#'
#' @references
#' Janssen MF, Szende A, Ramos-Goñi JM (2014). "Data and Methods." In
#'   \emph{Self-Reported Population Health: An International Perspective based
#'   on EQ-5D}, Chapter 3. Springer. \doi{10.1007/978-94-007-7596-1}
#'
#' Sullivan PW, Ghushchyan V (2006). "Preference-Based EQ-5D Index Scores for
#'   Chronic Conditions in the United States." \emph{Medical Decision Making},
#'   26(4), 410-420. \doi{10.1177/0272989X06290495}
#'
#' Sullivan PW, Slejko JF, Sculpher MJ, Ghushchyan V (2011). "Catalogue of
#'   EQ-5D Scores for the United Kingdom." \emph{Medical Decision Making},
#'   31(6), 800-804. \doi{10.1177/0272989X11401031}
#'
#' @export
Simulation <-
  R6::R6Class(
    classname = "Simulation",
    lock_objects = TRUE, # allows primary prevention scenario to be updated
    lock_class = FALSE, # allows adding methods via $set() from other files
    # public ------------------------------------------------------------------
    public = list(
      #' @field design A Design object.
      design = NA,

      #' @field diseases A list of Disease objects.
      diseases = NA,

      #' @field RR A list of RR for the simulated exposures.
      RR = NA,

      #' @field exposures A named list of Exposure objects for generating synthetic population exposures.
      exposures = NA,

      #' @field scenarios A list of scenario objects.
      scenarios = NA,

      # initialise ----
      #' @description Create a new simulation object.
      #' @param design Either a path to a yaml file or a Design object.
      #' @return A new `Simulation` object.
      initialize = function(design) {
        if (is.character(design)) {
          self$design <- Design$new(design)
        } else if (inherits(design, "Design")) {
          self$design <- design$clone(deep = TRUE)
        } else {
          stop(
            "design need to be a path to an appropriate yaml file or a Design object"
          )
        }

        # FIX: rocky uses OPENBLAS-OPENMP that creates issues with forking if
        # more than 1 thread is used before a foreach loop. The following code
        # fixes the issue by setting the number of threads to 1 from the
        # beginning of the simulation.
        blas_info <- sessionInfo()$BLAS
        is_openmp_blas <- grepl(
          "OPENMP|MKL",
          blas_info,
          ignore.case = TRUE
        )
        if (is_openmp_blas) {
          if (self$design$sim_prm$logs) {
            message(
              "Detected OPENMP/MKL BLAS library (",
              blas_info,
              "). Setting number of threads to 1 to avoid issues with forking."
            )
          }
          data.table::setDTthreads(
            threads = 1L,
            restore_after_fork = FALSE
          )
          fst::threads_fst(
            nr_of_threads = 1L,
            reset_after_fork = FALSE
          )
          arrow::set_cpu_count(1L) # limit Arrow's internal threading
        } else {
          if (self$design$sim_prm$logs) {
            message(
              "BLAS library detected: ",
              blas_info,
              ". No need to adjust number of threads."
            )
          }
          data.table::setDTthreads(
            threads = self$design$sim_prm$clusternumber,
            restore_after_fork = NULL
          )
          fst::threads_fst(
            nr_of_threads = self$design$sim_prm$clusternumber,
            reset_after_fork = NULL
          )
          arrow::set_cpu_count(self$design$sim_prm$clusternumber) # limit Arrow's internal threading
        }

        # Create folders if don't exist
        private$create_output_folder_structure()

        # European standardised population 2013 (esp) weights
        tt <- data.table(
          agegrp = agegrp_name(0, 99),
          wt_esp = c(
            1000,
            4000,
            5500,
            5500,
            5500,
            6000,
            6000,
            6500,
            7000,
            7000,
            7000,
            7000,
            6500,
            6000,
            5500,
            5000,
            4000,
            2500,
            1500,
            800,
            200
          )
        )
        esp <- CJ(
          agegrp = agegrp_name(0, 99),
          sex = c("men", "women"),
          dimd = c("1 most deprived", as.character(2:9), "10 least deprived")
        )

        private$esp_weights <- copy(absorb_dt(esp, tt))

        private$primary_prevention_scn <- function(synthpop) NULL # default for baseline scenario
        private$secondary_prevention_scn <- function(synthpop) NULL # default for baseline scenario

        # Complete data-dependent initialization if data is loaded
        if (self$design$data_loaded) {
          private$complete_data_init()
        }

        invisible(self)
      },

      # update_primary_prevention_scn ----
      #' @description Updates the primary prevention policy scenario
      #' @param method a function with synthpop as an argument that models the primary prevention policy.
      #' @return The invisible self for chaining.
      update_primary_prevention_scn = function(method) {
        private$primary_prevention_scn <- method
        environment(private$primary_prevention_scn) <- environment(
          private$update_primary_prevention_scn
        )
      },

      # get_primary_prevention_scn ----
      #' @description Get the primary prevention policy scenario
      #' @return The primary prevention policy scenario.
      get_primary_prevention_scn = function() {
        private$primary_prevention_scn
      },

      # update_secondary_prevention_scn ----
      #' @description Updates the secondary prevention policy scenario
      #' @param method a function with synthpop as an argument that models the secondary prevention policy.
      #' @return The invisible self for chaining.
      update_secondary_prevention_scn = function(method) {
        private$secondary_prevention_scn <- method
        environment(private$secondary_prevention_scn) <- environment(
          private$update_secondary_prevention_scn
        )
      },

      # get_secondary_prevention_scn ----
      #' @description Get the secondary prevention policy scenario
      #' @return The secondary prevention policy scenario.
      get_secondary_prevention_scn = function() {
        private$secondary_prevention_scn
      },

      # run ----
      #' @description Runs a simulation
      #' @param mc A positive sequential integer vector with the Monte Carlo
      #'   iterations of synthetic population to simulate, or a scalar.
      #' @param multicore If TRUE run the simulation in parallel.
      #' @param scenario_nam A string for the scenario name (e.g. "sc0",
      #'   "sc1"). Use only letters, digits, and underscores (e.g. "sc1",
      #'   "bmi_reduction_10pc"). Avoid hyphens, spaces, or other special
      #'   characters as these can cause issues in internal SQL queries.
      #'   The baseline scenario must always be named "sc0".
      #' @return The invisible self for chaining.
      run = function(mc, multicore = TRUE, scenario_nam) {
        if (!private$data_initialized) {
          stop(
            "Simulation not fully initialized. Input data files are missing.\n",
            "Download data first:\n",
            "  self$zenodo_connect(token, concept_doi, sandbox)\n",
            "  self$zenodo_download_all()",
            call. = FALSE
          )
        }
        if (!is.integer(mc)) {
          stop("mc need to be an integer")
        }
        if (any(mc <= 0)) {
          stop("mc need to be positive integer")
        }

        # Validate scenario_nam
        if (!missing(scenario_nam) && nzchar(scenario_nam) &&
            grepl("[^A-Za-z0-9_]", scenario_nam)) {
          warning(
            "scenario_nam '", scenario_nam, "' contains characters other ",
            "than letters, digits, and underscores. This is not recommended ",
            "and may cause issues. Please use only [A-Za-z0-9_] characters ",
            "(e.g. 'sc1', 'bmi_reduction_10pc').",
            call. = FALSE
          )
        }

        # recombine the chunks of large files
        # TODO logic to delete these files
        self$reconstruct_large_files()

        # check if sequential vector. Necessary if
        # design$sim_prm$num_chunks > 1
        if (
          anyNA(mc) ||
            any(is.infinite(mc)) ||
            length(mc) < 1L ||
            (length(mc) > 1L && diff(mc[1:2]) == 0) ||
            (length(mc) > 1L &&
              diff(range(diff(mc))) > sqrt(.Machine$double.eps))
        ) {
          stop("mc need to be a sequential integer vector, or a scalar")
        }
        # NOTE mc is in fact mc_aggr. mc_ is the mc of the synthpop
        mc_sp <-
          (min(mc) *
            self$design$sim_prm$num_chunks -
            self$design$sim_prm$num_chunks +
            1L):(max(mc) * self$design$sim_prm$num_chunks)

        # Create folders if don't exist (necessary for when output_dir in the
        # design.yaml is changed between scenarios i.e.)
        private$create_output_folder_structure()

        # Check for leftover results from a previous run of THIS scenario.
        # Each run writes lifecourse output to
        #   <output_dir>/lifecourse/mc=<n>/scenario=<nam>/<m>_lifecourse.parquet
        # so we only look at partitions belonging to scenario_nam — output
        # from other scenarios is legitimate and must not trigger the warning.
        lifecourse_dir <- file.path(
          self$design$sim_prm$output_dir,
          "lifecourse"
        )
        scenario_nam_eff <- if (missing(scenario_nam) || !nzchar(scenario_nam)) {
          "sc0"
        } else {
          scenario_nam
        }
        if (dir.exists(lifecourse_dir)) {
          scenario_files <- Sys.glob(file.path(
            lifecourse_dir,
            "mc=*",
            paste0("scenario=", scenario_nam_eff),
            "*"
          ))
          if (length(scenario_files) > 0L) {
            message(
              "Output from a previous run of scenario '", scenario_nam_eff,
              "' already exists. Matching files will be overwritten. ",
              "Remove them manually if this was unintentional."
            )
          }
        }

        # Generate PARF files if they don't exist. Note that generation is
        # multicore
        lapply(self$design$diseases, function(x) {
          x$gen_parf_files(self$design, self$design$diseases)
        })

        if (multicore) {
          if (self$design$sim_prm$logs) {
            private$time_mark("Start of parallelisation")
          }

          arrow::set_cpu_count(1L)
          data.table::setDTthreads(
            threads = 1L,
            restore_after_fork = NULL
          )
          fst::threads_fst(
            nr_of_threads = 1L,
            reset_after_fork = NULL
          )

          if (.Platform$OS.type == "windows") {
            cl <-
              makeClusterPSOCK(
                self$design$sim_prm$clusternumber,
                dryrun = FALSE,
                quiet = FALSE,
                rscript_startup = quote(local({
                  library(CKutils)
                  library(IMPACTncdEngland)
                  library(digest)
                  library(fst)
                  library(qs2)
                  library(wrswoR)
                  library(gamlss.dist)
                  library(dqrng)
                  library(data.table)
                })),
                rscript_args = c(
                  "--no-init-file",
                  "--no-site-file",
                  "--no-environ"
                ),
                setup_strategy = "parallel"
              ) # used for clustering. Windows compatible

            on.exit(if (exists("cl")) stopCluster(cl))

            xps_dt <- parLapplyLB(
              cl = cl,
              X = mc_sp,
              fun = function(x) private$run_sim(mc_ = x, scenario_nam)
            )
          } else {
            # used for forking. Only Linux/OSX compatible
            registerDoParallel(self$design$sim_prm$clusternumber)

            xps_dt <- foreach(
              mc_iter = mc_sp,
              .inorder = FALSE,
              .options.multicore = list(preschedule = FALSE),
              .verbose = self$design$sim_prm$logs,
              .packages = c(
                "R6",
                "digest",
                "qs2",
                "wrswoR",
                "gamlss.dist",
                "dqrng",
                "CKutils",
                "IMPACTncdEngland",
                "fst",
                "data.table"
              ),
              .export = ls(envir = globalenv()),
              .noexport = NULL # c("time_mark")
            ) %dopar%
              {
                private$run_sim(mc_ = mc_iter, scenario_nam)
              }
          }
          if (self$design$sim_prm$logs) {
            private$time_mark("End of parallelisation")
          }
        } else {
          # if multicore = FALSE
          if (self$design$sim_prm$logs) {
            private$time_mark("Start of single-core run")
          }
          # all implicit parallelisation
          arrow::set_cpu_count(self$design$sim_prm$clusternumber)
          data.table::setDTthreads(
            threads = self$design$sim_prm$clusternumber,
            restore_after_fork = NULL
          )
          fst::threads_fst(
            nr_of_threads = self$design$sim_prm$clusternumber,
            reset_after_fork = NULL
          )
          lapply(mc_sp, private$run_sim, scenario_nam)

          if (self$design$sim_prm$logs) {
            private$time_mark("End of single-core run")
          }
        }

        # Close any lingering sinks (safety cleanup)
        # Note: sink.number(type = "message") returns 0 or 2 (not a count),
        # so we use a single conditional call, not a while loop
        if (sink.number(type = "message") > 0L) {
          sink(type = "message")
        }
        while (sink.number(type = "output") > 0L) {
          sink(type = "output")
        }

        invisible(self)
      },

      # get_causal_structure ----

      #' @description Returns the causality matrix and optionally plots the
      #'   causality structure.
      #' @param processed If `TRUE` generates the causality matrix from the
      #'   graph.
      #' @param print_plot If `TRUE` prints the causal structure graph.
      #' @param focus If missing the whole causal structure is returned.
      #'  Otherwise, if a named node only the subgraph of the 1st order
      #'  neighbours that point to the given vertrice is returned.
      #' @param mode Character. When focus is specified, determines direction:
      #'   "in" for neighbors pointing to the node, "out" for neighbors
      #'   stemming from the node, "all" for both directions.
      #' @param order Integer. When focus is specified, determines how many
      #'   steps away to include neighbors. Use Inf for all possible orders.
      #' @return The processed causality matrix if `processed = TRUE` or the
      #'   graph otherwise.
      get_causal_structure = function(
        processed = TRUE,
        print_plot = FALSE,
        focus = FALSE,
        mode = "in",
        order = 1
      ) {
        if (!inherits(private$causality_structure, "igraph")) {
          stop(
            "Causality structure not initialised. ",
            "Please recreate the Simulation object."
          )
        }
        if (missing(focus)) {
          graph <- private$causality_structure
        } else {
          if (length(focus) > 1L) {
            stop("focus need to be scalar string.")
          }
          if (!focus %in% self$get_node_names()) {
            stop(
              "focus need to be a node name. Use get_node_names() to get the list of eligible values."
            )
          }
          if (!mode %in% c("in", "out", "all")) {
            stop("mode must be one of 'in', 'out', or 'all'.")
          }
          if (!is.numeric(order) || order < 1) {
            stop("order must be a positive integer or Inf.")
          }

          # Handle Inf order by using the graph diameter (maximum possible distance)
          if (is.infinite(order)) {
            # Use diameter of the graph as maximum order, or a large number if disconnected
            graph_diameter <- diameter(
              private$causality_structure,
              directed = TRUE
            )
            if (is.infinite(graph_diameter)) {
              # If graph is disconnected, use number of vertices as upper bound
              actual_order <- vcount(private$causality_structure)
            } else {
              actual_order <- graph_diameter
            }
          } else {
            actual_order <- order
          }

          graph <- make_ego_graph(
            private$causality_structure,
            order = actual_order,
            nodes = focus,
            mode = mode
          )[[1]]
        }
        if (print_plot) {
          if (vcount(graph) == 0L) {
            message("No causal structure to plot (no RR data loaded).")
          } else if (!.has_plot_device()) {
            message(
              "No interactive graphics device available. ",
              "To save the plot to a file, use:\n",
              "  png(\"causal_structure.png\", width = 1200, height = 800)\n",
              "  get_causal_structure(print_plot = TRUE)\n",
              "  dev.off()"
            )
          } else {
            lo <- layout_components(graph)
            plot.igraph(
              graph,
              vertex.shape = "none",
              edge.arrow.size = .3,
              vertex.label.font = 2,
              vertex.label.color = "gray40",
              edge.arrow.width = .5,
              vertex.label.cex = .7,
              edge.color = "gray85",
              layout = lo
            )
          }
        }

        if (processed) {
          graph <- as.matrix(as_adjacency_matrix(graph))
          n <- vapply(
            self$design$diseases,
            `[[`,
            "name",
            FUN.VALUE = character(1)
          )
          graph <- graph[rowSums(graph) > 0, colnames(graph) %in% n]
        }

        return(graph)
      },

      # get_node_names ----

      #' @description Returns the names of all exposures and diseases.
      #' @return A string vector.
      get_node_names = function() {
        return(V(private$causality_structure)$name)
      },

      # get_causal_path ----

      #' @description Returns the causal paths between an exposure and an outcome (disease).
      #' @param from the beginning of the path (an exposure) as a string. Use `get_node_names` for available nodes.
      #' @param to the end of the path (a disease) as a string. Use `get_node_names` for available nodes.
      #' @param shortest_paths Boolean. If true, only returns the paths with the smallest number of nodes. Else, all possible paths (excluding multiple and loop edges) are returned.
      #' @return A list with all the possible paths between exposure and disease.
      get_causal_path = function(from, to, shortest_paths = FALSE) {
        nm <- V(private$causality_structure)$name
        from <- which(nm == from)
        to <- which(nm == to)
        if (shortest_paths) {
          out <- get.all.shortest.paths(
            private$causality_structure,
            from,
            to,
            mode = "out"
          )
        } else {
          out <- all_simple_paths(
            private$causality_structure,
            from,
            to,
            mode = "out"
          )
        }
        return(out)
      },

      # get_inputs_manifest ----

      #' @description Returns the InputsManifest object for inspecting
      #'   tracked input files and their hashes.
      #' @return An InputsManifest object, or NULL if not initialised.
      get_inputs_manifest = function() {
        private$inputs_manifest
      },

      # update_design ----

      #' @description Updates the Design object that is stored in the Simulation
      #'   object.
      #' @param new_design A design object with the simulation parameters.
      #' @return The invisible self for chaining.
      update_design = function(new_design) {
        if (!inherits(new_design, "Design")) {
          stop("Argument new_design needs to be a Design object.")
        }

        self$design <- new_design

        invisible(self)
      },

      # del_outputs ----

      #' @description Delete all output files and folders below the first level.
      #' @return The invisible self for chaining.
      del_outputs = function() {
        if (dir.exists(self$design$sim_prm$output_dir)) {
          # Get all files in output_dir (including nested files)
          fl <- list.files(
            self$design$sim_prm$output_dir,
            full.names = TRUE,
            recursive = TRUE
          )

          # Remove all files
          file.remove(fl)

          # Get first-level directories (e.g., summaries/, logs/, tables/)
          first_level_dirs <- list.dirs(
            self$design$sim_prm$output_dir,
            full.names = TRUE,
            recursive = FALSE
          )
          # Remove the output_dir itself from the list
          first_level_dirs <- first_level_dirs[
            first_level_dirs != self$design$sim_prm$output_dir
          ]

          # For each first-level directory, get and remove all subdirectories
          subdirs_removed <- 0
          for (dir in first_level_dirs) {
            subdirs <- list.dirs(
              dir,
              full.names = TRUE,
              recursive = TRUE
            )
            # Remove the parent directory itself from the list (keep only subdirectories)
            subdirs <- subdirs[subdirs != dir]

            if (length(subdirs) > 0) {
              unlink(subdirs, recursive = TRUE)
              subdirs_removed <- subdirs_removed + length(subdirs)
            }
          }

          if (
            (length(fl) > 0 || subdirs_removed > 0) && self$design$sim_prm$logs
          ) {
            message(paste(
              "Output files deleted:",
              length(fl),
              "| Subdirectories removed:",
              subdirs_removed,
              "| First-level directories preserved:",
              length(first_level_dirs)
            ))
          }
        } else {
          message("Output folder doesn't exist.")
        }

        invisible(self)
      },

      # del_summaries ----
      #' @description Delete all output summary files and subdirectories while preserving first-level directory structure.
      #' @return The invisible self for chaining.
      del_summaries = function() {
        pth <- file.path(self$design$sim_prm$output_dir, "summaries")
        if (dir.exists(pth)) {
          # Get all items in the summaries directory
          all_items <- list.files(pth, full.names = TRUE, include.dirs = TRUE)

          # Separate files and directories
          files <- all_items[!dir.exists(all_items)]
          dirs <- all_items[dir.exists(all_items)]

          files_deleted <- 0
          subdirs_removed <- 0

          # Remove all files in the main summaries directory
          if (length(files) > 0) {
            file.remove(files)
            files_deleted <- length(files)
          }

          # For each subdirectory, remove it entirely (including all contents)
          if (length(dirs) > 0) {
            for (dir_path in dirs) {
              unlink(dir_path, recursive = TRUE)
              subdirs_removed <- subdirs_removed + 1
            }
          }

          if (
            self$design$sim_prm$logs &&
              (files_deleted > 0 || subdirs_removed > 0)
          ) {
            msg_parts <- character()
            if (files_deleted > 0) {
              msg_parts <- c(msg_parts, paste(files_deleted, "files deleted"))
            }
            if (subdirs_removed > 0) {
              msg_parts <- c(
                msg_parts,
                paste(subdirs_removed, "subdirectories removed")
              )
            }
            msg_parts <- c(msg_parts, "summaries directory preserved")
            message(
              "Output summary cleanup: ",
              paste(msg_parts, collapse = ", "),
              "."
            )
          }
        } else {
          message("Output summaries folder doesn't exist.")
        }

        invisible(self)
      },

      # del_logs ----
      #' @description Delete log files.
      #' @return The invisible self for chaining.
      del_logs = function() {
        fl <- list.files(private$output_dir("logs/"), full.names = TRUE)

        file.remove(fl)

        if (length(fl) > 0 && self$design$sim_prm$logs) {
          message("Log files deleted.")
        }

        invisible(self)
      },

      # del_parfs ----
      #' @description Delete all files in the ./simulation/parf folder.
      #' @return The invisible self for chaining.
      del_parfs = function() {
        fl <- list.files("./simulation/parf", full.names = TRUE)

        unlink(fl, recursive = TRUE)

        if (length(fl) > 0 && self$design$sim_prm$logs) {
          message("Parf files deleted.")
        }

        invisible(self)
      },

      # del_RR_cache ----
      #' @description Delete all files in the ./simulation/rr folder.
      #' @return The invisible self for chaining.
      del_RR_cache = function() {
        fl <- list.files("./simulation/rr", full.names = TRUE)

        file.remove(fl)

        if (length(fl) > 0 && self$design$sim_prm$logs) {
          message("RR cache files deleted.")
        }

        invisible(self)
      },

      # del_synthpops ----
      #' @description Delete all files in the synthpop folder.
      #' @return The invisible self for chaining.
      del_synthpops = function() {
        fl <- list.files(self$design$sim_prm$synthpop_dir, full.names = TRUE)

        unlink(fl, recursive = TRUE)

        if (length(fl) > 0 && self$design$sim_prm$logs) {
          message("Synthpop files deleted.")
        }

        invisible(self)
      },

      # get_esp ----

      #' @description Get the European Standardised Population 2013 by sex and
      #'   dimd.
      #' @return A data.table with the European Standardised Population 2013.
      get_esp = function() {
        private$esp_weights
      },
      # get_mm_weights ----

      #' @description Get the disease multimorbidity weights (i.e. Cambridge
      #'   Morbidity Score weights).
      #' @return A named vector with disease weights.
      get_mm_weights = function() {
        unlist(sapply(self$design$diseases, function(x) x$meta$diagnosis$mm_wt))
      },

      #' @description Internal validation of the disease burden.
      #' @return The invisible self for chaining.
      # validate ----
      validate = function() {
        HEIGHT <- 5
        WIDTH <- 10

        # TODO: Update path for England-specific population data
        data_pop <- read_fst(
          "./inputs/pop_projections/combined_population_england.fst",
          columns = c("year", "age", "sex", "pops"),
          as.data.table = TRUE
        )
        data_pop_agegrp <- copy(data_pop)
        to_agegrp(data_pop_agegrp, 5, 99)
        data_pop_agegrp <- data_pop_agegrp[,
          .(pops = sum(pops)),
          keyby = .(year, agegrp, sex)
        ]

        # MRTL
        mdd <- open_dataset(file.path(
          self$design$sim_prm$output_dir,
          "summaries",
          "dis_mrtl_scaled_up"
        )) %>%
          filter(scenario == "sc0") %>%
          collect()
        setDT(mdd)
        mdd[, `:=`(
          nonmodelled_mrtl_rate = nonmodelled_deaths / popsize,
          chd_mrtl_rate = chd_deaths / popsize,
          stroke_mrtl_rate = stroke_deaths / popsize
        )]
        mdd <- mdd[,
          .(
            nonmodelled_mrtl_rate = quantile(nonmodelled_mrtl_rate, p = 0.500),
            nonmodelled_mrtl_rate_low = quantile(
              nonmodelled_mrtl_rate,
              p = 0.025
            ),
            nonmodelled_mrtl_rate_upp = quantile(
              nonmodelled_mrtl_rate,
              p = 0.975
            ),
            chd_mrtl_rate = quantile(chd_mrtl_rate, p = 0.500),
            chd_mrtl_rate_low = quantile(chd_mrtl_rate, p = 0.025),
            chd_mrtl_rate_upp = quantile(chd_mrtl_rate, p = 0.957),
            stroke_mrtl_rate = quantile(stroke_mrtl_rate, p = 0.500),
            stroke_mrtl_rate_low = quantile(stroke_mrtl_rate, p = 0.025),
            stroke_mrtl_rate_upp = quantile(stroke_mrtl_rate, p = 0.975),
            type = "modelled"
          ),
          keyby = .(year, agegrp, sex)
        ]

        obs <- read_fst(
          paste0("./inputs/disease_burden/", "chd_ftlt.fst"),
          columns = c("age", "year", "sex", "mu2", "mu_lower", "mu_upper"),
          as.data.table = TRUE
        )
        setnames(
          obs,
          c("mu2", "mu_lower", "mu_upper"),
          c("chd_mrtl_rate", "chd_mrtl_rate_low", "chd_mrtl_rate_upp")
        )
        tt <- read_fst(
          paste0("./inputs/disease_burden/", "stroke_ftlt.fst"),
          columns = c("age", "year", "sex", "mu2", "mu_lower", "mu_upper"),
          as.data.table = TRUE
        )
        setnames(
          tt,
          c("mu2", "mu_lower", "mu_upper"),
          c("stroke_mrtl_rate", "stroke_mrtl_rate_low", "stroke_mrtl_rate_upp")
        )
        absorb_dt(obs, tt)
        tt <- read_fst(
          paste0("./inputs/disease_burden/", "nonmodelled_ftlt.fst"),
          columns = c("age", "year", "sex", "mu2", "mu_lower", "mu_upper"),
          as.data.table = TRUE
        )
        setnames(
          tt,
          c("mu2", "mu_lower", "mu_upper"),
          c(
            "nonmodelled_mrtl_rate",
            "nonmodelled_mrtl_rate_low",
            "nonmodelled_mrtl_rate_upp"
          )
        )
        absorb_dt(obs, tt)
        absorb_dt(obs, data_pop)
        to_agegrp(obs, 5, 99)
        obs <- obs[,
          lapply(.SD, weighted.mean, w = pops),
          .SDcols = -c("pops", "age"),
          keyby = .(agegrp, year, sex)
        ]
        obs[, type := "observed"]
        dt <- rbindlist(list(obs, mdd), use.names = TRUE)

        # Generate mortality plots
        for (sex_val in c("men", "women")) {
          for (metric in c("chd", "stroke", "nonmodelled")) {
            p <- ggplot() +
              geom_line(
                data = dt[sex == sex_val],
                aes(
                  x = year,
                  y = get(paste0(metric, "_mrtl_rate")),
                  color = type
                )
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(
                  x = year,
                  y = get(paste0(metric, "_mrtl_rate_low")),
                  color = type
                ),
                linetype = "dashed"
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(
                  x = year,
                  y = get(paste0(metric, "_mrtl_rate_upp")),
                  color = type
                ),
                linetype = "dashed"
              ) +
              facet_wrap(. ~ factor(agegrp), scales = "free") +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
              ggtitle(
                paste0(toupper(metric), " mrtl rate"),
                tools::toTitleCase(sex_val)
              )
            ggsave(
              file.path(
                self$design$sim_prm$output_dir,
                "plots",
                paste0(toupper(metric), "_as_", sex_val, "_mrtl.jpg")
              ),
              p,
              height = HEIGHT,
              width = WIDTH
            )
          }
        }

        # MRTL by sex alone
        absorb_dt(dt, data_pop_agegrp)
        dt <- dt[,
          lapply(.SD, weighted.mean, w = pops),
          .SDcols = -c("pops", "agegrp"),
          keyby = .(type, year, sex)
        ]

        for (metric in c("chd", "stroke", "nonmodelled")) {
          p <- ggplot() +
            geom_line(
              data = dt,
              aes(x = year, y = get(paste0(metric, "_mrtl_rate")), color = type)
            ) +
            geom_line(
              data = dt,
              aes(
                x = year,
                y = get(paste0(metric, "_mrtl_rate_low")),
                color = type
              ),
              linetype = "dashed"
            ) +
            geom_line(
              data = dt,
              aes(
                x = year,
                y = get(paste0(metric, "_mrtl_rate_upp")),
                color = type
              ),
              linetype = "dashed"
            ) +
            facet_wrap(. ~ factor(sex), scales = "free") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
            ggtitle(paste0(toupper(metric), " mrtl rate"), "By sex")
          ggsave(
            file.path(
              self$design$sim_prm$output_dir,
              "plots",
              paste0(toupper(metric), "_s_mrtl.jpg")
            ),
            p,
            height = HEIGHT,
            width = WIDTH
          )
        }

        # INCD
        mdd <- open_dataset(file.path(
          self$design$sim_prm$output_dir,
          "summaries",
          "incd_scaled_up"
        )) %>%
          filter(scenario == "sc0") %>%
          select(
            "mc",
            "scenario",
            "year",
            "agegrp",
            "sex",
            "popsize",
            "chd_incd",
            "stroke_incd"
          ) %>%
          collect()
        setDT(mdd)
        mdd[, `:=`(
          chd_incd_rate = chd_incd / popsize,
          stroke_incd_rate = stroke_incd / popsize
        )]
        mdd <- mdd[,
          .(
            chd_incd_rate = quantile(chd_incd_rate, p = 0.500),
            chd_incd_rate_low = quantile(chd_incd_rate, p = 0.025),
            chd_incd_rate_upp = quantile(chd_incd_rate, p = 0.957),
            stroke_incd_rate = quantile(stroke_incd_rate, p = 0.500),
            stroke_incd_rate_low = quantile(stroke_incd_rate, p = 0.025),
            stroke_incd_rate_upp = quantile(stroke_incd_rate, p = 0.975),
            type = "modelled"
          ),
          keyby = .(year, agegrp, sex)
        ]

        obs <- read_fst(
          paste0("./inputs/disease_burden/", "chd_incd.fst"),
          columns = c("age", "year", "sex", "mu", "mu_lower", "mu_upper"),
          as.data.table = TRUE
        )
        setnames(
          obs,
          c("mu", "mu_lower", "mu_upper"),
          c("chd_incd_rate", "chd_incd_rate_low", "chd_incd_rate_upp")
        )
        tt <- read_fst(
          paste0("./inputs/disease_burden/", "stroke_incd.fst"),
          columns = c("age", "year", "sex", "mu", "mu_lower", "mu_upper"),
          as.data.table = TRUE
        )
        setnames(
          tt,
          c("mu", "mu_lower", "mu_upper"),
          c("stroke_incd_rate", "stroke_incd_rate_low", "stroke_incd_rate_upp")
        )
        absorb_dt(obs, tt)
        absorb_dt(obs, data_pop)
        to_agegrp(obs, 5, 99)
        obs <- obs[,
          lapply(.SD, weighted.mean, w = pops),
          .SDcols = -c("pops", "age"),
          keyby = .(agegrp, year, sex)
        ]
        obs[, type := "observed"]
        dt <- rbindlist(list(obs, mdd), use.names = TRUE)

        # Generate incidence plots
        for (sex_val in c("men", "women")) {
          for (metric in c("chd", "stroke")) {
            p <- ggplot() +
              geom_line(
                data = dt[sex == sex_val],
                aes(
                  x = year,
                  y = get(paste0(metric, "_incd_rate")),
                  color = type
                )
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(
                  x = year,
                  y = get(paste0(metric, "_incd_rate_low")),
                  color = type
                ),
                linetype = "dashed"
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(
                  x = year,
                  y = get(paste0(metric, "_incd_rate_upp")),
                  color = type
                ),
                linetype = "dashed"
              ) +
              facet_wrap(. ~ factor(agegrp), scales = "free") +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
              ggtitle(
                paste0(toupper(metric), " incd rate"),
                tools::toTitleCase(sex_val)
              )
            ggsave(
              file.path(
                self$design$sim_prm$output_dir,
                "plots",
                paste0(toupper(metric), "_as_", sex_val, "_incd.jpg")
              ),
              p,
              height = HEIGHT,
              width = WIDTH
            )
          }
        }

        # INCD by sex alone
        absorb_dt(dt, data_pop_agegrp)
        dt <- dt[,
          lapply(.SD, weighted.mean, w = pops),
          .SDcols = -c("pops", "agegrp"),
          keyby = .(type, year, sex)
        ]

        for (metric in c("chd", "stroke")) {
          p <- ggplot() +
            geom_line(
              data = dt,
              aes(x = year, y = get(paste0(metric, "_incd_rate")), color = type)
            ) +
            geom_line(
              data = dt,
              aes(
                x = year,
                y = get(paste0(metric, "_incd_rate_low")),
                color = type
              ),
              linetype = "dashed"
            ) +
            geom_line(
              data = dt,
              aes(
                x = year,
                y = get(paste0(metric, "_incd_rate_upp")),
                color = type
              ),
              linetype = "dashed"
            ) +
            facet_wrap(. ~ factor(sex), scales = "free") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
            ggtitle(paste0(toupper(metric), " incd rate"), "By sex")
          ggsave(
            file.path(
              self$design$sim_prm$output_dir,
              "plots",
              paste0(toupper(metric), "_s_incd.jpg")
            ),
            p,
            height = HEIGHT,
            width = WIDTH
          )
        }

        # PRVL
        mdd <- open_dataset(file.path(
          self$design$sim_prm$output_dir,
          "summaries",
          "prvl_scaled_up"
        )) %>%
          filter(scenario == "sc0") %>%
          select(
            "mc",
            "scenario",
            "year",
            "agegrp",
            "sex",
            "popsize",
            "chd_prvl",
            "stroke_prvl"
          ) %>%
          collect()
        setDT(mdd)

        mdd[, `:=`(
          chd_prvl_rate = chd_prvl / popsize,
          stroke_prvl_rate = stroke_prvl / popsize
        )]
        mdd <- mdd[,
          .(
            chd_prvl_rate = quantile(chd_prvl_rate, p = 0.500),
            chd_prvl_rate_low = quantile(chd_prvl_rate, p = 0.025),
            chd_prvl_rate_upp = quantile(chd_prvl_rate, p = 0.957),
            stroke_prvl_rate = quantile(stroke_prvl_rate, p = 0.500),
            stroke_prvl_rate_low = quantile(stroke_prvl_rate, p = 0.025),
            stroke_prvl_rate_upp = quantile(stroke_prvl_rate, p = 0.975),
            type = "modelled"
          ),
          keyby = .(year, agegrp, sex)
        ]

        obs <- read_fst(
          paste0("./inputs/disease_burden/", "chd_prvl.fst"),
          columns = c(
            "age",
            "year",
            "sex",
            "mu",
            "mu_lower",
            "mu_upper",
            "prvl_mltp"
          ),
          as.data.table = TRUE
        )
        setnames(
          obs,
          c("mu", "mu_lower", "mu_upper"),
          c("chd_prvl_rate", "chd_prvl_rate_low", "chd_prvl_rate_upp")
        )
        obs[, `:=`(
          chd_prvl_rate = chd_prvl_rate * prvl_mltp,
          chd_prvl_rate_low = chd_prvl_rate_low * prvl_mltp,
          chd_prvl_rate_upp = chd_prvl_rate_upp * prvl_mltp,
          prvl_mltp = NULL
        )]
        tt <- read_fst(
          paste0("./inputs/disease_burden/", "stroke_prvl.fst"),
          columns = c(
            "age",
            "year",
            "sex",
            "mu",
            "mu_lower",
            "mu_upper",
            "prvl_mltp"
          ),
          as.data.table = TRUE
        )
        setnames(
          tt,
          c("mu", "mu_lower", "mu_upper"),
          c("stroke_prvl_rate", "stroke_prvl_rate_low", "stroke_prvl_rate_upp")
        )
        tt[, `:=`(
          stroke_prvl_rate = stroke_prvl_rate * prvl_mltp,
          stroke_prvl_rate_low = stroke_prvl_rate_low * prvl_mltp,
          stroke_prvl_rate_upp = stroke_prvl_rate_upp * prvl_mltp,
          prvl_mltp = NULL
        )]
        absorb_dt(obs, tt)
        absorb_dt(obs, data_pop)
        to_agegrp(obs, 5, 99)
        obs <- obs[,
          lapply(.SD, weighted.mean, w = pops),
          .SDcols = -c("pops", "age"),
          keyby = .(agegrp, year, sex)
        ]
        obs[, type := "observed"]
        dt <- rbindlist(list(obs, mdd), use.names = TRUE)

        # Generate prevalence plots
        for (sex_val in c("men", "women")) {
          for (metric in c("chd", "stroke")) {
            p <- ggplot() +
              geom_line(
                data = dt[sex == sex_val],
                aes(
                  x = year,
                  y = get(paste0(metric, "_prvl_rate")),
                  color = type
                )
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(
                  x = year,
                  y = get(paste0(metric, "_prvl_rate_low")),
                  color = type
                ),
                linetype = "dashed"
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(
                  x = year,
                  y = get(paste0(metric, "_prvl_rate_upp")),
                  color = type
                ),
                linetype = "dashed"
              ) +
              facet_wrap(. ~ factor(agegrp), scales = "free") +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
              ggtitle(
                paste0(toupper(metric), " prvl rate"),
                tools::toTitleCase(sex_val)
              )
            ggsave(
              file.path(
                self$design$sim_prm$output_dir,
                "plots",
                paste0(toupper(metric), "_as_", sex_val, "_prvl.jpg")
              ),
              p,
              height = HEIGHT,
              width = WIDTH
            )
          }
        }

        # prvl by sex alone
        absorb_dt(dt, data_pop_agegrp)
        dt <- dt[,
          lapply(.SD, weighted.mean, w = pops),
          .SDcols = -c("pops", "agegrp"),
          keyby = .(type, year, sex)
        ]

        for (metric in c("chd", "stroke")) {
          p <- ggplot() +
            geom_line(
              data = dt,
              aes(x = year, y = get(paste0(metric, "_prvl_rate")), color = type)
            ) +
            geom_line(
              data = dt,
              aes(
                x = year,
                y = get(paste0(metric, "_prvl_rate_low")),
                color = type
              ),
              linetype = "dashed"
            ) +
            geom_line(
              data = dt,
              aes(
                x = year,
                y = get(paste0(metric, "_prvl_rate_upp")),
                color = type
              ),
              linetype = "dashed"
            ) +
            facet_wrap(. ~ factor(sex), scales = "free") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
            ggtitle(paste0(toupper(metric), " prvl rate"), "By sex")
          ggsave(
            file.path(
              self$design$sim_prm$output_dir,
              "plots",
              paste0(toupper(metric), "_s_prvl.jpg")
            ),
            p,
            height = HEIGHT,
            width = WIDTH
          )
        }
        invisible(self)
      },

      # split_large_files ----
      #' @description Splits files larger than 50Mb into chunks of 49Mb.
      #' @details The function splits files larger than 50Mb into chunks of 49Mb
      #' so they can be tracked by GitHub. The large files are deleted and an
      #' index is created at "./simulation/large_files_indx.csv" so they can be
      #' reconstructed. The function also adds the large files to `.gitignore`.
      #' It works on Linux and Windows. Untested on Mac.
      #' @return The invisible `Simulation` object.
      split_large_files = function() {
        # identify large files
        fl <- list.files(".", full.names = TRUE, recursive = TRUE)
        fl <- sort(fl[file.size(fl) / (1024^2) >= 50])
        fl <- grep(
          "/synthpop/|/outputs/",
          fl,
          ignore.case = TRUE,
          value = TRUE,
          invert = TRUE
        )
        if (length(fl) == 0) {
          # no large files. Early escape.
          return(invisible(self))
        }

        # Merge with existing large files from previous uploads. This ensures
        # that if any previously large files are now smaller than 50MB, they are
        # still processed as if they are still more than 50Mb. This is perhaps
        # inefficient but has less side effects.  Otherwise code for special
        # case of files above the threshold that shrinked for whatever reason
        # below the threshold is needed (i.e. remove them from .gitignore).
        if (file.exists("./simulation/large_files_indx.csv")) {
          fl <- sort(unique(c(
            fread("./simulation/large_files_indx.csv")$pths,
            fl
          )))
        }
        fwrite(list(pths = fl), "./simulation/large_files_indx.csv")

        # add large files to .gitignore
        excl <- readLines("./.gitignore")
        for (i in seq_along(fl)) {
          file <- gsub("^./", "", fl[i])
          if (file %in% excl) {
            next
          }
          write(file, file = "./.gitignore", append = TRUE)
        }

        # split the files into 50MB chunks
        for (i in seq_along(fl)) {
          file <- fl[i]
          if (!file.exists(file)) {
            next
          }

          # split the file into 49MB chunks
          if (.Platform$OS.type == "unix") {
            system(paste0("split -b 49m ", file, " ", file, ".chunk"))
          } else if (.Platform$OS.type == "windows") {
            # For windows split and cat are from https://unxutils.sourceforge.net/
            shell(paste0("split -b 49m ", file, " ", file, ".chunk"))
          } else {
            stop("Operating system is not supported.")
          }
          # remove the original file
          unlink(file, recursive = TRUE)
        }

        invisible(self)
      },

      # reconstruct_large_files ----
      #' @description Reconstructs large files from chunks.
      #' @details The function reconstructs large files from chunks. The path of
      #' the files are stored in "./simulation/large_files_indx.csv". It works
      #' on Linux and Windows. Untested on Mac.
      #' @return The invisible `Simulation` object.
      reconstruct_large_files = function() {
        if (file.exists("./simulation/large_files_indx.csv")) {
          fl <- fread("./simulation/large_files_indx.csv")$pths
          for (i in seq_along(fl)) {
            if (file.exists(fl[i])) {
              next
            }
            file <- fl[i]
            # recombine the chunks
            if (.Platform$OS.type == "unix") {
              system(paste0("cat ", file, ".chunk?? > ", file, ""))
            } else if (.Platform$OS.type == "windows") {
              # For windows split and cat are from https://unxutils.sourceforge.net/
              shell(paste0("cat ", file, ".chunk?? > ", file, ""))
            } else {
              stop("Operating system is not supported.")
            }
          }
        }
        invisible(self)
      },

      # del_large_files ----
      #' @description Deletes large files.
      #' @details The function deletes large files that are stored in chunks.
      #' The path of the files are stored in
      #' "./simulation/large_files_indx.csv".
      #' @return The invisible `Simulation` object.
      del_large_files = function() {
        if (file.exists("./simulation/large_files_indx.csv")) {
          fl <- fread("./simulation/large_files_indx.csv")$pths
          unlink(fl, recursive = TRUE)
        }
        invisible(self)
      },

      # zenodo_connect ----
      #' @description Connect to Zenodo for uploading and downloading input
      #'   data assets.
      #'
      #' @details
      #' This is the first step in any Zenodo workflow. It creates an internal
      #' \code{\link{ZenodoAssetManager}} object, authenticates with the Zenodo
      #' API, and (optionally) sets the concept DOI so that subsequent methods
      #' know which record to work with.
      #'
      #' \strong{What is a concept DOI?}
      #' Every Zenodo record can have multiple versions. Each version receives
      #' its own DOI (e.g. \code{10.5281/zenodo.12346}), but they all share a
      #' single \emph{concept DOI} (e.g. \code{10.5281/zenodo.12345}) that
      #' always resolves to the latest published version. The concept DOI is the
      #' one you should store in your code and share with collaborators.
      #'
      #' \strong{Token management:}
      #' You need a personal access token from Zenodo. You can either:
      #' \enumerate{
      #'   \item Pass it directly: \code{token = "your_token_here"}
      #'   \item Set the environment variable \code{ZENODO_TOKEN} (recommended
      #'     for security — add it to \code{.Renviron})
      #'   \item Use a file-based keyring (see the
      #'     \code{vignette("zenodo_data_management")} for setup instructions)
      #' }
      #'
      #' \strong{Sandbox vs production:}
      #' Set \code{sandbox = TRUE} for testing. The sandbox uses a separate
      #' token from \url{https://sandbox.zenodo.org/account/settings/applications/tokens/new/}.
      #' Published sandbox records are not permanent and do not count as real
      #' publications.
      #'
      #' @param token Character. Zenodo personal access token. If \code{NULL}
      #'   (default), reads from the \code{ZENODO_TOKEN} environment variable;
      #'   if that is also unset, connects in \strong{anonymous read-only mode}.
      #'   A token is only needed to upload or publish — downloading the
      #'   published public data works without one.
      #' @param concept_doi Character. The concept DOI for your data record.
      #'   Defaults to the published IMPACTncd England input-data record:
      #'   \code{"10.5281/zenodo.20812409"} for production, or
      #'   \code{"10.5072/zenodo.442996"} when \code{sandbox = TRUE}. The
      #'   concept DOI always resolves to the latest published version. Set to
      #'   \code{NULL} to skip (e.g. for a brand-new first upload).
      #' @param sandbox Logical. If \code{TRUE}, connect to the Zenodo sandbox
      #'   for testing. Default: \code{FALSE}.
      #' @param archive_dir Character. Directory for storing zip archives
      #'   during upload. Archives are automatically deleted after
      #'   successful upload. If \code{NULL} (default), uses
      #'   \code{tempdir()}.
      #' @param progress Logical. If \code{TRUE} (default), enable console
      #'   progress bars for uploads and downloads.
      #' @return The invisible \code{Simulation} object (for method chaining).
      #'
      #' @examples
      #' # Connect to sandbox for testing
      #' IMPACTncd$zenodo_connect(
      #'   token   = Sys.getenv("ZENODO_SANDBOX_TOKEN"),
      #'   sandbox = TRUE
      #' )
      #'
      #' # Connect to production (defaults to the published input-data record;
      #' # no token needed to download published data)
      #' IMPACTncd$zenodo_connect()
      #'
      #' @seealso \code{\link{zenodo_upload_inputs}},
      #'   \code{\link{zenodo_download_inputs}},
      #'   \code{\link{zenodo_list_files}}
      zenodo_connect = function(
        token = NULL,
        concept_doi = if (sandbox) "10.5072/zenodo.442996" else "10.5281/zenodo.20812409",
        sandbox = FALSE,
        archive_dir = NULL,
        progress = TRUE
      ) {
        if (is.null(private$zenodo_manager)) {
          private$zenodo_manager <- ZenodoAssetManager$new(
            hash_file = "./simulation/zenodo_manifest.csv",
            archive_dir = archive_dir,
            logs = self$design$sim_prm$logs,
            sandbox = sandbox
          )
        }

        private$zenodo_manager$connect(token = token, sandbox = sandbox)

        if (isTRUE(progress)) {
          private$zenodo_manager$set_progress_callback(
            upload = TRUE,
            download = TRUE
          )
        }

        if (!is.null(concept_doi)) {
          private$zenodo_manager$set_concept_doi(concept_doi)
        }

        invisible(self)
      },

      # zenodo_set_doi ----
      #' @description Set the Zenodo concept DOI for data assets.
      #'
      #' @details
      #' Use this to switch to a different Zenodo record after connecting.
      #' The concept DOI is the DOI shared across all versions of a record
      #' (prefix \code{10.5281} for production, \code{10.5072} for sandbox).
      #'
      #' @param concept_doi Character. The concept DOI string, e.g.
      #'   \code{"10.5281/zenodo.12345"}.
      #' @return The invisible \code{Simulation} object (for method chaining).
      #' @seealso \code{\link{zenodo_connect}}
      zenodo_set_doi = function(concept_doi) {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }
        private$zenodo_manager$set_concept_doi(concept_doi)
        invisible(self)
      },

      # zenodo_check_inputs ----
      #' @description Check whether local input directories are in sync with
      #'   the Zenodo record.
      #'
      #' @details
      #' Compares the local \code{input_base} directory against the file list
      #' on the Zenodo record. Returns a \code{data.table} showing, for each
      #' remote archive, whether the corresponding local directory exists.
      #' This is a read-only operation — nothing is downloaded or uploaded.
      #'
      #' Requires a prior call to \code{zenodo_connect()} with a concept DOI.
      #'
      #' @param input_base Character. Path to the base directory containing
      #'   input subdirectories. Default: \code{"./inputs"}.
      #' @param directories Character vector. Specific subdirectories to check.
      #'   If \code{NULL} (default), checks all remote archives.
      #' @return A \code{data.table} with columns: \code{archive},
      #'   \code{directory}, \code{local_exists}, and \code{remote_checksum}.
      #'
      #' @examples
      #' IMPACTncd$zenodo_connect(sandbox = TRUE)
      #' IMPACTncd$zenodo_check_inputs()
      #'
      #' @seealso \code{\link{zenodo_download_inputs}},
      #'   \code{\link{zenodo_list_files}}
      zenodo_check_inputs = function(
        input_base = "./inputs",
        directories = NULL
      ) {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }
        private$zenodo_manager$sync_inputs(
          input_base = input_base,
          directories = directories,
          action = "check"
        )
      },

      # zenodo_download_inputs ----
      #' @description Download input data archives from Zenodo and extract
      #'   them into local directories.
      #'
      #' @details
      #' Downloads zip archives from the Zenodo record and extracts them into
      #' \code{input_base}. Extraction is per-file: by default missing data
      #' files are added and existing files are kept (so a fresh git clone,
      #' whose data directories already exist holding git-tracked \code{.yaml}/
      #' \code{.R} files, still receives all the data) — pass
      #' \code{overwrite = TRUE} to replace existing files instead.
      #'
      #' Use \code{zenodo_check_inputs()} first (a read-only check) to see
      #' exactly which directories would be downloaded.
      #'
      #' Progress bars are shown if \code{progress = TRUE} was set in
      #' \code{zenodo_connect()} (the default).
      #'
      #' @param input_base Character. Destination directory for extracted
      #'   archives. Default: \code{"./inputs"}.
      #' @param directories Character vector. Specific subdirectories to
      #'   download. If \code{NULL} (default), downloads all remote archives.
      #' @param overwrite Logical. If \code{TRUE}, overwrite existing local
      #'   directories. Default: \code{FALSE}.
      #' @return The invisible \code{Simulation} object (for method chaining).
      #'
      #' @examples
      #' # Download all missing inputs
      #' IMPACTncd$zenodo_connect(sandbox = TRUE)
      #' IMPACTncd$zenodo_download_inputs()
      #'
      #' # Force re-download of specific directories
      #' IMPACTncd$zenodo_download_inputs(
      #'   directories = c("mortality", "population"),
      #'   overwrite = TRUE
      #' )
      #'
      #' @seealso \code{\link{zenodo_check_inputs}},
      #'   \code{\link{zenodo_connect}}
      zenodo_download_inputs = function(
        input_base = "./inputs",
        directories = NULL,
        overwrite = FALSE
      ) {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }

        if (self$design$sim_prm$logs) {
          message("Downloading inputs from Zenodo...")
        }

        # Filter out simulation archives (parf/rr) — those are handled
        # by zenodo_download_PARFs_RRs() which extracts to ./simulation/
        sim_pattern <- "^(parf|rr)"

        # Check current status
        status <- private$zenodo_manager$sync_inputs(
          input_base = input_base,
          directories = directories,
          action = "check",
          exclude_archive_patterns = sim_pattern
        )

        # Note existing directories. Extraction is per-file: missing data
        # files are added and existing files are kept (use overwrite = TRUE to
        # replace existing files). This matters on a fresh git clone, where
        # data directories already exist holding git-tracked .yaml/.R files.
        if (nrow(status) > 0L && any(status$local_exists) && !overwrite) {
          existing_dirs <- status[local_exists == TRUE, directory]
          message(
            "These directories already exist; missing data files will be ",
            "added and existing files kept (use overwrite=TRUE to replace):\n  ",
            paste(existing_dirs, collapse = "\n  ")
          )
        }

        # Download missing/all
        private$zenodo_manager$sync_inputs(
          input_base = input_base,
          directories = directories,
          action = "download",
          overwrite = overwrite,
          exclude_archive_patterns = sim_pattern
        )

        if (self$design$sim_prm$logs) {
          message("Input download complete.")
        }

        invisible(self)
      },

      # zenodo_download_PARFs_RRs ----
      #' @description Download simulation output archives (PARF, RR files)
      #'   from the Zenodo record and extract them into local directories.
      #'
      #' @details
      #' This method downloads zip archives of simulation-generated files
      #' (Population Attributable Risk Fractions and compiled Relative Risk
      #' tables) from Zenodo and extracts them to \code{simulation_base}.
      #'
      #' The method identifies simulation archives by matching archive names
      #' against the \code{directories} parameter (e.g. archives whose name
      #' starts with \code{"parf"} or equals \code{"rr"}).
      #'
      #' Extraction is per-file: by default missing data files are added and
      #' existing files are kept — pass \code{overwrite = TRUE} to replace
      #' existing files instead.
      #'
      #' \strong{When to use:}
      #' \itemize{
      #'   \item When you need pre-computed PARFs and RR tables from a
      #'     colleague's upload
      #'   \item When setting up a new machine and want to skip the
      #'     time-consuming PARF/RR generation step
      #' }
      #'
      #' @param simulation_base Character. Destination directory for
      #'   extracted simulation archives. Default: \code{"./simulation"}.
      #' @param directories Character vector. Subdirectory prefixes to
      #'   download. Archives whose name starts with any of these prefixes
      #'   are included. Default: \code{c("parf", "rr")}.
      #' @param overwrite Logical. If \code{TRUE}, overwrite existing local
      #'   directories. Default: \code{FALSE}.
      #' @return The invisible \code{Simulation} object (for method chaining).
      #'
      #' @examples
      #' IMPACTncd$zenodo_connect(sandbox = TRUE)
      #'
      #' # Download missing simulation outputs
      #' IMPACTncd$zenodo_download_PARFs_RRs()
      #'
      #' # Force re-download of PARF files only
      #' IMPACTncd$zenodo_download_PARFs_RRs(
      #'   directories = "parf",
      #'   overwrite   = TRUE
      #' )
      #'
      #' @seealso \code{\link{zenodo_download_inputs}},
      #'   \code{\link{zenodo_download_all}},
      #'   \code{\link{zenodo_upload_simulation}}
      zenodo_download_PARFs_RRs = function(
        simulation_base = "./simulation",
        directories = c("parf", "rr"),
        overwrite = FALSE
      ) {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }

        if (self$design$sim_prm$logs) {
          message("Downloading simulation data from Zenodo...")
        }

        # Get remote file list
        remote_files <- private$zenodo_manager$list_remote_files()

        if (nrow(remote_files) == 0L) {
          message("No files found in remote record.")
          return(invisible(self))
        }

        # Filter to simulation archives: match names starting with any
        # of the requested directory prefixes
        dir_pattern <- paste0("^(", paste(directories, collapse = "|"), ")")
        sim_mask <- grepl(dir_pattern, remote_files$filename)

        sim_files <- remote_files[sim_mask, ]

        if (nrow(sim_files) == 0L) {
          message(
            "No simulation archives found matching: ",
            paste(directories, collapse = ", ")
          )
          return(invisible(self))
        }

        if (self$design$sim_prm$logs) {
          message("Found ", nrow(sim_files), " simulation archives to download.")
        }

        # Download and extract each matching archive individually
        download_dir <- file.path(tempdir(), "zenodo_sim_downloads")
        dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)

        for (i in seq_len(nrow(sim_files))) {
          fname <- sim_files$filename[i]

          # download_file() uses per-file URL, never bulk downloadFiles()
          private$zenodo_manager$download_file(
            filename  = fname,
            dest_dir  = download_dir,
            overwrite = TRUE
          )

          archive_path <- file.path(download_dir, fname)
          if (!file.exists(archive_path)) {
            warning("Archive not found after download: ", fname)
            next
          }

          if (self$design$sim_prm$logs) {
            message("Extracting: ", fname, " -> ", simulation_base)
          }

          private$zenodo_manager$extract_archive(
            archive_path, simulation_base, overwrite = overwrite
          )
        }

        # Cleanup
        unlink(download_dir, recursive = TRUE)

        if (self$design$sim_prm$logs) {
          message("Simulation download complete.")
        }

        invisible(self)
      },

      # zenodo_download_all ----
      #' @description Download both input data and simulation outputs from
      #'   Zenodo in a single call.
      #'
      #' @details
      #' This is a convenience method that combines
      #' \code{zenodo_download_inputs()} and
      #' \code{zenodo_download_PARFs_RRs()} into a single operation.
      #'
      #' \strong{Typical use case:} A new team member setting up their
      #' environment for the first time. One call downloads everything
      #' needed to run the simulation.
      #'
      #' @param input_base Character. Destination for input data. Default:
      #'   \code{"./inputs"}.
      #' @param input_directories Character vector. Input subdirectories to
      #'   download. \code{NULL} = all.
      #' @param simulation_base Character. Destination for simulation data.
      #'   Default: \code{"./simulation"}.
      #' @param simulation_directories Character vector. Simulation archive
      #'   prefixes to download. Default: \code{c("parf", "rr")}.
      #' @param overwrite Logical. Overwrite existing directories? Default:
      #'   \code{FALSE}.
      #' @return The invisible \code{Simulation} object (for method chaining).
      #'
      #' @examples
      #' # New team member setup — download everything
      #' IMPACTncd$zenodo_connect(
      #'   concept_doi = "10.5281/zenodo.20812409"
      #' )
      #' IMPACTncd$zenodo_download_all()
      #'
      #' # Force re-download everything
      #' IMPACTncd$zenodo_download_all(overwrite = TRUE)
      #'
      #' @seealso \code{\link{zenodo_download_inputs}},
      #'   \code{\link{zenodo_download_PARFs_RRs}}
      zenodo_download_all = function(
        input_base = "./inputs",
        input_directories = NULL,
        simulation_base = "./simulation",
        simulation_directories = c("parf", "rr"),
        overwrite = FALSE
      ) {
        self$zenodo_download_inputs(
          input_base  = input_base,
          directories = input_directories,
          overwrite   = overwrite
        )

        self$zenodo_download_PARFs_RRs(
          simulation_base = simulation_base,
          directories     = simulation_directories,
          overwrite       = overwrite
        )

        # Auto-complete initialization after download
        if (!self$design$data_loaded) {
          self$design$load_data()
        }
        if (!private$data_initialized) {
          private$complete_data_init()
        }

        invisible(self)
      },

      # zenodo_upload_inputs ----
      #' @description Create zip archives of input directories and upload
      #'   them to Zenodo.
      #'
      #' @details
      #' This method handles the complete upload pipeline:
      #' \enumerate{
      #'   \item Creates zip archives from your input directories
      #'   \item Creates a new Zenodo record (first time) or a new version
      #'     (subsequent times)
      #'   \item Uploads all archives to the record
      #'   \item Optionally publishes the record
      #' }
      #'
      #' \strong{First-time upload (new record):}
      #' When no concept DOI has been set (i.e. you did not pass one to
      #' \code{zenodo_connect()}), you \emph{must} provide \code{title},
      #' \code{description}, and \code{creators}. A new Zenodo deposit is
      #' created in DRAFT state. Review it on the Zenodo website before
      #' publishing with \code{zenodo_publish()}.
      #'
      #' \strong{Subsequent uploads (new version):}
      #' If a concept DOI is already set (from \code{zenodo_connect()}),
      #' provide a \code{version} string (e.g. \code{"1.1.0"}) and a new
      #' version of the existing record is created automatically.
      #'
      #' \strong{Archive creation parameters:}
      #' \itemize{
      #'   \item \code{exclude_patterns}: Regex patterns to skip entire
      #'     \emph{directories}. For example, \code{"^unprocessed$"} skips a
      #'     directory named exactly "unprocessed".
      #'   \item \code{exclude_file_patterns}: Regex patterns to skip
      #'     \emph{individual files}. For example, \code{"\\.R$"} skips all
      #'     R scripts.
      #'   \item \code{group_by_prefix}: When \code{TRUE} (default),
      #'     subdirectories sharing a common prefix are combined into a single
      #'     archive. For example, \code{disease_burden/cancer_2020/} and
      #'     \code{disease_burden/cancer_2021/} become
      #'     \code{disease_burden_cancer.zip}.
      #'   \item \code{multicore}: When \code{TRUE} (default), archives are
      #'     created in parallel using PSOCK clusters (safe on all platforms).
      #' }
      #'
      #' @param input_base Character. Path to the base directory containing
      #'   input subdirectories. Default: \code{"./inputs"}.
      #' @param directories Character vector. Specific subdirectories to
      #'   upload. If \code{NULL} (default), uploads all subdirectories in
      #'   \code{input_base}.
      #' @param version Character. Version string for the upload (e.g.
      #'   \code{"1.0.0"}). Required for creating new versions of existing
      #'   records. For first-time uploads, defaults to \code{"1.0.0"} if
      #'   not specified.
      #' @param title Character. Record title (required for new records).
      #' @param description Character. Record description (required for new
      #'   records).
      #' @param creators List of lists. Each inner list must have
      #'   \code{firstname} and \code{lastname}, and optionally \code{orcid}.
      #'   Required for new records.
      #' @param keywords Character vector. Metadata keywords for
      #'   discoverability. Default: \code{c("microsimulation", "health",
      #'   "IMPACTncd", "England", "NCD")}.
      #' @param license Character. SPDX license identifier. Default:
      #'   \code{"cc-by-sa-4.0"} (Creative Commons Attribution-ShareAlike).
      #' @param publisher Character. Publisher name (required by Zenodo for
      #'   DOI registration). Default: \code{"Zenodo"}.
      #' @param exclude_patterns Character vector of regex patterns for
      #'   directories to skip. Default: \code{c("^unprocessed$", "_backup$",
      #'   "_old$", "scripts$", "validation$", "^RR$")}.
      #' @param exclude_file_patterns Character vector of regex patterns for
      #'   files to skip. Default: \code{c("\\.R$", "\\.Rmd$")}.
      #' @param group_by_prefix Logical. Group subdirectories by shared
      #'   prefix into combined archives. Default: \code{TRUE}.
      #' @param compression_level Integer 0-9. Zip compression level (0 =
      #'   store only, 9 = maximum compression). Default: \code{6}.
      #' @param multicore Logical. Create archives in parallel using PSOCK
      #'   clusters. Default: \code{TRUE}.
      #' @param n_cores Integer. Number of CPU cores for parallel archive
      #'   creation. Default: \code{clusternumber} from the simulation design
      #'   YAML.
      #' @param update_gitignore Logical. Automatically add archived data files
      #'   to \code{.gitignore} and remove any already-tracked data files from
      #'   the git index (\code{git rm --cached}). Files are kept on disk —
      #'   only git tracking is removed. Default: \code{TRUE}.
      #' @param publish Logical. If \code{TRUE}, publish the record
      #'   immediately after upload. \strong{WARNING: publishing is
      #'   IRREVERSIBLE} — the DOI becomes permanent and files cannot be
      #'   modified. Default: \code{FALSE}.
      #' @return The invisible \code{Simulation} object (for method chaining).
      #'
      #' @examples
      #' # --- First-time upload (creates new record) ---
      #' IMPACTncd$zenodo_connect(
      #'   token   = Sys.getenv("ZENODO_SANDBOX_TOKEN"),
      #'   sandbox = TRUE,
      #'   concept_doi = NULL  # no existing record
      #' )
      #' IMPACTncd$zenodo_upload_inputs(
      #'   title       = "IMPACTncd England Input Data",
      #'   description = "Input data for IMPACTncd England microsimulation",
      #'   creators    = list(
      #'     list(firstname = "Chris", lastname = "Kypridemos",
      #'          orcid = "0000-0002-0746-9229")
      #'   ),
      #'   version     = "1.0.0"
      #' )
      #' # Review the draft on Zenodo, then:
      #' IMPACTncd$zenodo_publish()
      #'
      #' # --- Subsequent upload (new version) ---
      #' IMPACTncd$zenodo_connect(sandbox = TRUE)
      #' IMPACTncd$zenodo_upload_inputs(version = "1.1.0")
      #'
      #' @seealso \code{\link{zenodo_connect}},
      #'   \code{\link{zenodo_upload_simulation}},
      #'   \code{\link{zenodo_upload_all}},
      #'   \code{\link{zenodo_publish}},
      #'   \code{\link{zenodo_create_archives}}
      zenodo_upload_inputs = function(
        input_base = "./inputs",
        directories = NULL,
        version = NULL,
        title = NULL,
        description = NULL,
        creators = NULL,
        keywords = c("microsimulation", "health", "IMPACTncd",
                      "England", "NCD"),
        license = "cc-by-sa-4.0",
        publisher = "Zenodo",
        exclude_patterns = c("^unprocessed$", "_backup$", "_old$",
                              "scripts$", "validation$", "^RR$"),
        exclude_file_patterns = c("\\.R$", "\\.Rmd$"),
        group_by_prefix = TRUE,
        compression_level = 6L,
        multicore = TRUE,
        n_cores = self$design$sim_prm$clusternumber,
        update_gitignore = TRUE,
        publish = FALSE
      ) {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }

        # Check if we need to create a new record or a new version
        has_doi <- !is.null(private$zenodo_manager$concept_doi)

        if (!has_doi) {
          # Creating a brand-new record — metadata is required
          if (is.null(title) || is.null(description) || is.null(creators)) {
            stop(
              "For new records, you must provide: title, description, ",
              "and creators.\n",
              "Example:\n",
              "  creators = list(\n",
              "    list(firstname = 'Chris', lastname = 'Kypridemos',\n",
              "         orcid = '0000-0002-0746-9229')\n",
              "  )"
            )
          }

          if (self$design$sim_prm$logs) {
            message("Creating new Zenodo record...")
          }

          private$zenodo_manager$create_new_record(
            title = title,
            description = description,
            creators = creators,
            version = version %||% "1.0.0",
            keywords = keywords,
            license = license,
            publisher = publisher
          )
        } else if (!is.null(version)) {
          # Creating new version of existing record
          if (self$design$sim_prm$logs) {
            message("Creating new version of existing record...")
          }
          private$zenodo_manager$get_record()
          private$zenodo_manager$create_new_version(
            version = version,
            delete_previous_files = TRUE
          )
        }

        # Upload the archives (pass all archive creation params through)
        if (self$design$sim_prm$logs) {
          message("Uploading input archives to Zenodo...")
        }

        archives <- private$zenodo_manager$sync_inputs(
          input_base = input_base,
          directories = directories,
          action = "upload",
          version = version,
          exclude_patterns = exclude_patterns,
          exclude_file_patterns = exclude_file_patterns,
          group_by_prefix = group_by_prefix,
          compression_level = compression_level,
          multicore = multicore,
          n_cores = n_cores,
          update_gitignore = update_gitignore
        )

        # Snapshot the current inputs manifest so that
        # compare_with_inputs_manifest() can detect changes since this upload
        private$zenodo_manager$snapshot_manifest_for_upload()

        if (publish) {
          if (self$design$sim_prm$logs) {
            message("Publishing record (this is IRREVERSIBLE)...")
          }
          private$zenodo_manager$publish_record()
        } else {
          if (self$design$sim_prm$logs) {
            message(
              "\nRecord is in DRAFT state.",
              "\nReview at Zenodo, then call:",
              "\n  IMPACTncd$zenodo_publish()",
              "\nto make it publicly available."
            )
          }
        }

        invisible(self)
      },

      # zenodo_upload_simulation ----
      #' @description Upload simulation output archives (PARF, RR files) to
      #'   the current Zenodo record.
      #'
      #' @details
      #' This method archives and uploads simulation-generated files
      #' (Population Attributable Risk Fractions and compiled Relative Risk
      #' tables) to the \emph{same} Zenodo record used for input data. Call
      #' this \strong{after} \code{zenodo_upload_inputs()} so that a record
      #' already exists.
      #'
      #' These files live under \code{./simulation/} rather than
      #' \code{./inputs/}, and use different default exclusion patterns
      #' (e.g. CSV manifest files are excluded).
      #'
      #' @param simulation_base Character. Path to the simulation output
      #'   directory. Default: \code{"./simulation"}.
      #' @param directories Character vector. Subdirectories to archive and
      #'   upload. Default: \code{c("parf", "rr")}.
      #' @param exclude_file_patterns Character vector of regex patterns for
      #'   files to skip. Default: \code{c("\\.R$", "\\.csv$")}.
      #' @param group_by_prefix Logical. Group subdirectories by prefix.
      #'   Default: \code{TRUE}.
      #' @param compression_level Integer 0-9. Zip compression level.
      #'   Default: \code{6}.
      #' @param multicore Logical. Create archives in parallel. Default:
      #'   \code{TRUE}.
      #' @param n_cores Integer. Number of CPU cores. Default:
      #'   \code{clusternumber} from the simulation design YAML.
      #' @param update_gitignore Logical. Auto-update \code{.gitignore} and
      #'   remove tracked data files from git index. Default: \code{TRUE}.
      #' @return The invisible \code{Simulation} object (for method chaining).
      #'
      #' @examples
      #' # Upload inputs first, then simulation outputs
      #' IMPACTncd$zenodo_connect(sandbox = TRUE)
      #' IMPACTncd$zenodo_upload_inputs(version = "1.0.0", title = "...",
      #'   description = "...", creators = list(...))
      #' IMPACTncd$zenodo_upload_simulation()
      #' IMPACTncd$zenodo_publish()
      #'
      #' @seealso \code{\link{zenodo_upload_inputs}},
      #'   \code{\link{zenodo_upload_all}}
      zenodo_upload_simulation = function(
        simulation_base = "./simulation",
        directories = c("parf", "rr"),
        exclude_file_patterns = c("\\.R$", "\\.csv$"),
        group_by_prefix = TRUE,
        compression_level = 6L,
        multicore = TRUE,
        n_cores = self$design$sim_prm$clusternumber,
        update_gitignore = TRUE
      ) {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }
        if (is.null(private$zenodo_manager$record)) {
          stop(
            "No Zenodo record available. Call zenodo_upload_inputs() first ",
            "to create a record, or use zenodo_connect() with a concept DOI."
          )
        }

        if (self$design$sim_prm$logs) {
          message("Creating simulation archives...")
        }

        archives <- private$zenodo_manager$create_input_archives(
          input_base = simulation_base,
          directories = directories,
          exclude_file_patterns = exclude_file_patterns,
          group_by_prefix = group_by_prefix,
          compression_level = compression_level,
          multicore = multicore,
          n_cores = n_cores,
          update_gitignore = update_gitignore
        )

        if (nrow(archives) == 0L) {
          message("No simulation archives to upload.")
          return(invisible(self))
        }

        if (self$design$sim_prm$logs) {
          message("Uploading ", nrow(archives), " simulation archives...")
        }

        private$zenodo_manager$upload_archives(archives$archive_path)

        # Clean up temporary archive files after successful upload
        for (ap in archives$archive_path) {
          if (file.exists(ap)) file.remove(ap)
        }

        if (self$design$sim_prm$logs) {
          total_mb <- round(sum(archives$size_bytes) / 1024^2, 1)
          message("Simulation upload complete. Total: ", total_mb, " MB")
        }

        invisible(self)
      },

      # zenodo_upload_all ----
      #' @description Upload both input data and simulation outputs to
      #'   Zenodo in a single call.
      #'
      #' @details
      #' This is a convenience method that combines
      #' \code{zenodo_upload_inputs()} and \code{zenodo_upload_simulation()}
      #' into a single operation. It:
      #' \enumerate{
      #'   \item Creates and uploads input archives (from \code{./inputs/})
      #'   \item Creates and uploads simulation archives (from
      #'     \code{./simulation/parf/} and \code{./simulation/rr/})
      #'   \item Optionally publishes the record
      #' }
      #'
      #' This mirrors the complete workflow demonstrated in the test script.
      #'
      #' @param input_base Character. Path to input directory. Default:
      #'   \code{"./inputs"}.
      #' @param input_directories Character vector. Input subdirectories.
      #'   \code{NULL} = all.
      #' @param simulation_base Character. Path to simulation output
      #'   directory. Default: \code{"./simulation"}.
      #' @param simulation_directories Character vector. Simulation
      #'   subdirectories. Default: \code{c("parf", "rr")}.
      #' @param version Character. Version string (e.g. \code{"1.0.0"}).
      #' @param title Character. Record title (required for new records).
      #' @param description Character. Record description (required for new
      #'   records).
      #' @param creators List of lists with creator metadata (required for
      #'   new records).
      #' @param keywords Character vector. Default:
      #'   \code{c("microsimulation", "health", "IMPACTncd", "England",
      #'   "NCD")}.
      #' @param license Character. SPDX license ID. Default:
      #'   \code{"cc-by-sa-4.0"}.
      #' @param publisher Character. Default: \code{"Zenodo"}.
      #' @param exclude_patterns Character vector. Directory exclusion
      #'   patterns for inputs. Default:
      #'   \code{c("^unprocessed$", "_backup$", "_old$", "scripts$",
      #'   "validation$", "^RR$")}. \code{^RR$} excludes
      #'   \code{inputs/RR/} because the source CSVY files are tracked in
      #'   git rather than uploaded to Zenodo.
      #' @param input_exclude_file_patterns Character vector. File exclusion
      #'   patterns for inputs. Default: \code{c("\\.R$", "\\.Rmd$")}.
      #' @param simulation_exclude_file_patterns Character vector. File
      #'   exclusion patterns for simulation archives. Default:
      #'   \code{c("\\.R$", "\\.csv$")}.
      #' @param group_by_prefix Logical. Default: \code{TRUE}.
      #' @param compression_level Integer 0-9. Default: \code{6}.
      #' @param multicore Logical. Default: \code{TRUE}.
      #' @param n_cores Integer. Default: \code{clusternumber} from YAML.
      #' @param update_gitignore Logical. Auto-update \code{.gitignore} and
      #'   remove tracked data files from git index. Default: \code{TRUE}.
      #' @param publish Logical. Publish after upload? Default: \code{FALSE}.
      #' @return The invisible \code{Simulation} object (for method chaining).
      #'
      #' @examples
      #' # Complete first-time upload of everything
      #' IMPACTncd$zenodo_connect(
      #'   token   = Sys.getenv("ZENODO_SANDBOX_TOKEN"),
      #'   sandbox = TRUE,
      #'   concept_doi = NULL
      #' )
      #' IMPACTncd$zenodo_upload_all(
      #'   title       = "IMPACTncd England Input Data",
      #'   description = "All input and simulation data for IMPACTncd England",
      #'   creators    = list(
      #'     list(firstname = "Chris", lastname = "Kypridemos",
      #'          orcid = "0000-0002-0746-9229")
      #'   ),
      #'   version = "1.0.0"
      #' )
      #'
      #' # Verify what was uploaded
      #' IMPACTncd$zenodo_list_files()
      #'
      #' # Publish when satisfied
      #' IMPACTncd$zenodo_publish()
      #'
      #' @seealso \code{\link{zenodo_upload_inputs}},
      #'   \code{\link{zenodo_upload_simulation}},
      #'   \code{\link{zenodo_publish}}
      zenodo_upload_all = function(
        input_base = "./inputs",
        input_directories = NULL,
        simulation_base = "./simulation",
        simulation_directories = c("parf", "rr"),
        version = NULL,
        title = NULL,
        description = NULL,
        creators = NULL,
        keywords = c("microsimulation", "health", "IMPACTncd",
                      "England", "NCD"),
        license = "cc-by-sa-4.0",
        publisher = "Zenodo",
        exclude_patterns = c("^unprocessed$", "_backup$", "_old$",
                              "scripts$", "validation$", "^RR$"),
        input_exclude_file_patterns = c("\\.R$", "\\.Rmd$"),
        simulation_exclude_file_patterns = c("\\.R$", "\\.csv$"),
        group_by_prefix = TRUE,
        compression_level = 6L,
        multicore = TRUE,
        n_cores = self$design$sim_prm$clusternumber,
        update_gitignore = TRUE,
        publish = FALSE
      ) {
        # Step 1: Upload inputs (handles record creation / versioning)
        self$zenodo_upload_inputs(
          input_base = input_base,
          directories = input_directories,
          version = version,
          title = title,
          description = description,
          creators = creators,
          keywords = keywords,
          license = license,
          publisher = publisher,
          exclude_patterns = exclude_patterns,
          exclude_file_patterns = input_exclude_file_patterns,
          group_by_prefix = group_by_prefix,
          compression_level = compression_level,
          multicore = multicore,
          n_cores = n_cores,
          update_gitignore = update_gitignore,
          publish = FALSE  # defer — publish after simulation upload
        )

        # Step 2: Upload simulation archives to the same record
        self$zenodo_upload_simulation(
          simulation_base = simulation_base,
          directories = simulation_directories,
          exclude_file_patterns = simulation_exclude_file_patterns,
          group_by_prefix = group_by_prefix,
          compression_level = compression_level,
          multicore = multicore,
          n_cores = n_cores,
          update_gitignore = update_gitignore
        )

        # Step 3: Publish if requested
        if (publish) {
          self$zenodo_publish()
        }

        invisible(self)
      },

      # zenodo_publish ----
      #' @description Publish the current Zenodo record draft.
      #'
      #' @details
      #' \strong{CAUTION: Publishing is IRREVERSIBLE.} Once published:
      #' \itemize{
      #'   \item The DOI becomes permanent and publicly resolvable
      #'   \item Files in this version \strong{cannot} be modified
      #'   \item The record appears in search results and is citable
      #' }
      #'
      #' You can still create \emph{new versions} of a published record by
      #' calling \code{zenodo_upload_inputs(version = "2.0.0")}. Each version
      #' gets its own DOI, but they all share the concept DOI.
      #'
      #' Always review the draft on the Zenodo website before publishing.
      #'
      #' @return The invisible \code{Simulation} object (for method chaining).
      #'
      #' @examples
      #' IMPACTncd$zenodo_connect(sandbox = TRUE)
      #' IMPACTncd$zenodo_upload_inputs(...)
      #' # Review the draft at Zenodo, then:
      #' IMPACTncd$zenodo_publish()
      #'
      #' @seealso \code{\link{zenodo_upload_inputs}},
      #'   \code{\link{zenodo_connect}}
      zenodo_publish = function() {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }

        if (self$design$sim_prm$logs) {
          message("WARNING: Publishing is IRREVERSIBLE!")
          message("The DOI will be permanent and files cannot be modified.")
        }

        private$zenodo_manager$publish_record()

        if (self$design$sim_prm$logs) {
          message("Record published successfully.")
          message("DOI: ", private$zenodo_manager$record$pids$doi$identifier)
        }

        invisible(self)
      },

      # zenodo_list_files ----
      #' @description List all files in the current Zenodo record.
      #'
      #' @details
      #' Returns a \code{data.table} with one row per file, including
      #' filename, size, and checksum. Useful for verifying that uploads
      #' completed successfully.
      #'
      #' Works with both published records and drafts (uses a direct API
      #' fallback for drafts where zen4R's built-in method may return
      #' empty results).
      #'
      #' @param concept_doi Character. Optional concept DOI. If \code{NULL}
      #'   (default), uses the DOI set via \code{zenodo_connect()}.
      #' @return A \code{data.table} with file metadata (key, size, checksum,
      #'   etc.).
      #'
      #' @examples
      #' IMPACTncd$zenodo_connect(sandbox = TRUE)
      #' IMPACTncd$zenodo_list_files()
      #'
      #' @seealso \code{\link{zenodo_check_inputs}},
      #'   \code{\link{zenodo_connect}}
      zenodo_list_files = function(concept_doi = NULL) {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }
        private$zenodo_manager$list_remote_files(concept_doi)
      },

      # zenodo_get_versions ----
      #' @description Get the version history of the Zenodo record.
      #'
      #' @details
      #' Returns a data.frame listing all published versions of the record,
      #' including each version's DOI, version label, and publication date.
      #' Useful for auditing which data versions have been published.
      #'
      #' @return A \code{data.frame} with version information.
      #'
      #' @examples
      #' IMPACTncd$zenodo_connect(sandbox = TRUE)
      #' IMPACTncd$zenodo_get_versions()
      #'
      #' @seealso \code{\link{zenodo_connect}}
      zenodo_get_versions = function() {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }
        private$zenodo_manager$get_versions()
      },

      # zenodo_create_archives ----
      #' @description Create zip archives from input directories without
      #'   uploading them.
      #'
      #' @details
      #' Creates archives exactly as \code{zenodo_upload_inputs()} would,
      #' but stops before uploading. This lets you inspect archive sizes
      #' and contents before committing to an upload. Archives are saved
      #' to \code{archive_dir} (set in \code{zenodo_connect()}, or
      #' \code{tempdir()} by default).
      #'
      #' Does not require a Zenodo connection — only needs the simulation
      #' design to be loaded (for \code{n_cores} default).
      #'
      #' @param input_base Character. Path to the base directory containing
      #'   input subdirectories. Default: \code{"./inputs"}.
      #' @param directories Character vector. Specific subdirectories.
      #'   If \code{NULL} (default), processes all.
      #' @param exclude_patterns Character vector. Directory exclusion
      #'   patterns. Default: \code{c("^unprocessed$", "_backup$", "_old$",
      #'   "scripts$", "validation$", "^RR$")}.
      #' @param exclude_file_patterns Character vector. File exclusion
      #'   patterns. Default: \code{c("\\.R$", "\\.Rmd$")}.
      #' @param group_by_prefix Logical. Default: \code{TRUE}.
      #' @param compression_level Integer 0-9. Default: \code{6}.
      #' @param multicore Logical. Default: \code{TRUE}.
      #' @param n_cores Integer. Default: \code{clusternumber} from YAML.
      #' @param update_gitignore Logical. Default: \code{FALSE} (unlike
      #'   \code{zenodo_upload_inputs()}, since this is a preview step).
      #' @return A \code{data.table} with archive metadata: path, size,
      #'   file count, source hash.
      #'
      #' @examples
      #' # Preview what would be created
      #' archives <- IMPACTncd$zenodo_create_archives()
      #' print(archives[, .(archive = basename(archive_path),
      #'                     size_mb = round(size_bytes / 1024^2, 1),
      #'                     files = file_count)])
      #'
      #' @seealso \code{\link{zenodo_upload_inputs}}
      zenodo_create_archives = function(
        input_base = "./inputs",
        directories = NULL,
        exclude_patterns = c("^unprocessed$", "_backup$", "_old$",
                              "scripts$", "validation$", "^RR$"),
        exclude_file_patterns = c("\\.R$", "\\.Rmd$"),
        group_by_prefix = TRUE,
        compression_level = 6L,
        multicore = TRUE,
        n_cores = self$design$sim_prm$clusternumber,
        update_gitignore = FALSE
      ) {
        private$ensure_zenodo_manager()
        private$zenodo_manager$create_input_archives(
          input_base = input_base,
          directories = directories,
          exclude_patterns = exclude_patterns,
          exclude_file_patterns = exclude_file_patterns,
          group_by_prefix = group_by_prefix,
          compression_level = compression_level,
          multicore = multicore,
          n_cores = n_cores,
          update_gitignore = update_gitignore
        )
      },

      # zenodo_create_manifest ----
      #' @description Create a hash manifest of local input files.
      #'
      #' @details
      #' Scans \code{input_base} recursively and computes an xxhash64
      #' checksum for every file. The resulting manifest can be saved to
      #' disk and later used to detect which files have changed.
      #'
      #' This is a local operation — it does not require a Zenodo
      #' connection.
      #'
      #' @param input_base Character. Path to the input directory. Default:
      #'   \code{"./inputs"}.
      #' @param save Logical. If \code{TRUE} (default), save the manifest
      #'   to \code{./simulation/zenodo_manifest.csv}.
      #' @return A \code{data.table} with columns: \code{relative_path},
      #'   \code{hash}, \code{size_bytes}, \code{mtime}.
      #'
      #' @examples
      #' manifest <- IMPACTncd$zenodo_create_manifest()
      #' print(manifest)
      #'
      #' @seealso \code{\link{zenodo_upload_inputs}}
      zenodo_create_manifest = function(input_base = "./inputs", save = TRUE) {
        private$ensure_zenodo_manager()

        manifest <- private$zenodo_manager$compute_manifest(input_base)

        if (save) {
          private$zenodo_manager$save_manifest(manifest)
        }

        manifest
      },

      # print ----
      #' @description Prints the simulation object metadata.
      #' @return The invisible `Simulation` object.
      print = function() {
        print(c(
          "TODO..."
        ))
        invisible(self)
      },
      # allow_universal_output_folder_access ----

      #' @description Make output folder available to all users (Linux specific).
      #' @return The invisible self for chaining.
      allow_universal_output_folder_access = function() {
        if (Sys.info()["sysname"] == "Linux") {
          system2("chmod", paste0("ugo+rwx ", self$design$sim_prm$output_dir))
        } else {
          message("This function is only available in Linux.")
        }
        invisible(self)
      },
      # allow_universal_synthpop_folder_access ----

      #' @description Make synthpop folder available to all users (Linux specific).
      #' @return The invisible self for chaining.
      allow_universal_synthpop_folder_access = function() {
        if (Sys.info()["sysname"] == "Linux") {
          system2("chmod", paste0("ugo+rwx ", self$design$sim_prm$synthpop_dir))
        } else {
          message("This function is only available in Linux.")
        }
        invisible(self)
      },
      # update_output_path ----

      #' @description Updates the output path.
      #' @param new_path A string with the new output path (absolute).
      #' @param carry_over_lifecourse_files_from A string with a previous output
      #' path (absolute) from which the lifecourse files will be copied to the
      #' new output folder defined in new_path argument. Overwritting is not
      #' allowed. If missing, no copy occurs.
      #' @return The invisible self for chaining.
      update_output_path = function(
        new_path,
        carry_over_lifecourse_files_from
      ) {
        if (!is.character(new_path)) {
          stop("new_path needs to be a string.")
        }
        if (
          !missing(carry_over_lifecourse_files_from) &&
            !is.character(carry_over_lifecourse_files_from)
        ) {
          stop("carry_over_lifecourse_files_from needs to be a string.")
        }
        if (
          !missing(carry_over_lifecourse_files_from) &&
            !dir.exists(carry_over_lifecourse_files_from)
        ) {
          stop(
            "Folder defined with carry_over_lifecourse_files_from does not exist."
          )
        }

        if (!missing(carry_over_lifecourse_files_from)) {
          fl <- list.files(
            file.path(carry_over_lifecourse_files_from, "lifecourse"),
            pattern = "lifecourse",
            full.names = TRUE
          )
        }
        self$design$sim_prm$output_dir <- new_path
        private$create_output_folder_structure()
        if (!missing(carry_over_lifecourse_files_from)) {
          file.copy(fl, file.path(self$design$sim_prm$output_dir, "lifecourse"))
        }
        invisible(self)
      },
      # update_synthpop_path ----

      #' @description Updates the synthpop path.
      #' @param new_path A string with the new synthpop path (absolute).
      #' @return The invisible self for chaining.
      update_synthpop_path = function(new_path) {
        if (!is.character(new_path)) {
          stop("new_path needs to be a string.")
        }
        self$design$sim_prm$synthpop_dir <- new_path
        invisible(self)
      }
    ),
    # private -----------------------------------------------------------------
    private = list(
      synthpop_dir = NA,
      causality_structure = NA,
      death_codes = NA,
      # diseasenam_hlp = NA,
      esp_weights = data.table(),
      # Models a primary prevention policy scenario
      primary_prevention_scn = NULL,
      # Models a secondary prevention policy scenario
      secondary_prevention_scn = NULL,
      # Zenodo asset manager for input data
      zenodo_manager = NULL,
      # Inputs manifest for tracking file versions
      inputs_manifest = NULL,
      # Whether data-dependent initialization is complete
      data_initialized = FALSE,

      # complete_data_init ----
      # Performs data-dependent Phase 2 initialization:
      # builds the causality graph, extracts death codes, and
      # runs the inputs manifest check. Called automatically by
      # the constructor when data is present, or by
      # zenodo_download_all() after downloading data.
      complete_data_init = function() {
        if (private$data_initialized) return(invisible(NULL))

        message("Generating microsimulation structure.")
        # Generate the graph with the causality structure
        if (length(self$design$RR) > 0) {
          ds <- unlist(strsplit(names(self$design$RR), "~"))
          ds[grep("^smok_", ds)] <- "smoking"
          ds <- gsub("_prvl$", "", ds)

          ds1 <- ds[as.logical(seq_along(ds) %% 2)]
          ds2 <- ds[!as.logical(seq_along(ds) %% 2)]
          ds <- unique(data.table(ds1, ds2))

          private$causality_structure <- make_graph(
            unlist(transpose(ds)),
            directed = TRUE
          )
        } else {
          private$causality_structure <- make_empty_graph(directed = TRUE)
        }

        private$death_codes <- unlist(lapply(
          self$design$diseases, function(x) {
            x$meta$mortality$code
          }
        ))
        private$death_codes[["alive"]] <- 0L

        # --- INPUTS MANIFEST CHECK ---
        private$inputs_manifest <- InputsManifest$new(
          inputs_dir = "./inputs",
          manifest_path = "./simulation/inputs_manifest.csv"
        )

        manifest_path <- "./simulation/inputs_manifest.csv"
        if (file.exists(manifest_path)) {
          private$inputs_manifest$load()
          changes <- private$inputs_manifest$detect_changes(
            self$design
          )

          if (length(changes$all_changed) > 0L) {
            message("=== INPUT FILE CHANGES DETECTED ===")
            message(sprintf(
              "  %d added, %d removed, %d modified",
              length(changes$added),
              length(changes$removed),
              length(changes$modified)
            ))

            if (changes$affected_synthpops) {
              message(
                "  -> Synthpops will be regenerated",
                " (upstream data changed)"
              )
            }
            if (length(changes$affected_rr) > 0L) {
              message(
                "  -> RR files to regenerate: ",
                paste(changes$affected_rr, collapse = ", ")
              )
              # Delete stale compiled RR .fst files so
              # ExposureEffect regenerates them
              for (rr_name in changes$affected_rr) {
                rr_files <- list.files(
                  "./simulation/rr",
                  pattern = paste0(
                    "^rr_",
                    rr_name,
                    "_(l|indx)\\.fst$"
                  ),
                  full.names = TRUE
                )
                if (length(rr_files) > 0L) {
                  file.remove(rr_files)
                }
              }
            }
            if (length(changes$affected_diseases) > 0L) {
              message(
                "  -> PARF files to regenerate: ",
                paste(
                  changes$affected_diseases,
                  collapse = ", "
                )
              )
              # Delete stale PARF files so Disease
              # regenerates them
              parf_dir <- "./simulation/parf"
              if (dir.exists(parf_dir)) {
                for (d in changes$affected_diseases) {
                  parf_files <- list.files(
                    parf_dir,
                    pattern = paste0("^PARF_", d, "_"),
                    full.names = TRUE
                  )
                  if (length(parf_files) > 0L) {
                    file.remove(parf_files)
                  }
                }
              }
            }

            # Update manifest to reflect current state
            private$inputs_manifest$generate(
              parallel = self$design$sim_prm$clusternumber > 1L,
              n_workers = self$design$sim_prm$clusternumber
            )
            message(
              "  Manifest updated.",
              " Stale artifacts will be regenerated."
            )
            message("===================================")
          }
        } else {
          message(
            "No inputs manifest found.",
            " Generating initial manifest..."
          )
          private$inputs_manifest$generate(
            parallel = self$design$sim_prm$clusternumber > 1L,
            n_workers = self$design$sim_prm$clusternumber
          )
          message(sprintf(
            "Inputs manifest created: %d files tracked.",
            nrow(private$inputs_manifest$manifest)
          ))
        }

        private$data_initialized <- TRUE
        invisible(NULL)
      },

      # ensure_zenodo_manager ----
      # Lazily create a ZenodoAssetManager without connecting to Zenodo.
      # Used by local-only methods like zenodo_create_archives() and
      # zenodo_create_manifest() that don't need an API connection.
      ensure_zenodo_manager = function() {
        if (is.null(private$zenodo_manager)) {
          private$zenodo_manager <- ZenodoAssetManager$new(
            hash_file = "./simulation/zenodo_manifest.csv",
            logs = self$design$sim_prm$logs
          )
        }
      },

      # Helper function to execute database query and write to disk with retry logic
      # execute_db_diskwrite_with_retry ----
      execute_db_diskwrite_with_retry = function(
        duckdb_con,
        query,
        output_path,
        max_retries = 5
      ) {
        retry_count <- 0
        success <- FALSE

        # Normalize path for cross-platform compatibility
        output_path <- normalizePath(output_path, mustWork = FALSE)

        # Ensure parent directory exists and is accessible
        parent_dir <- dirname(output_path)
        if (!dir.exists(parent_dir)) {
          dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
        }

        # Add a small delay to ensure directory is fully created (Windows issue)
        Sys.sleep(0.05)

        # Check if we're running in Docker on Windows (SMB/Plan9 mount issue)
        # Docker Desktop on Windows uses SMB/Plan9 mounts which don't support atomic operations
        is_docker_env <- file.exists("/.dockerenv")
        is_docker_windows <- (is_docker_env &&
          (.Platform$OS.type == "unix" ||
            Sys.getenv("DOCKER_DESKTOP") != "" ||
            # Check for typical Windows Docker mount patterns
            any(grepl("/host_mnt/", getwd(), fixed = TRUE))))

        # Use safer method for Docker environments (especially Windows Docker Desktop)
        # or any Windows environment where atomic operations might not be reliable
        if (
          is_docker_windows || is_docker_env || .Platform$OS.type == "windows"
        ) {
          env_type <- if (is_docker_windows) {
            "Docker on Windows"
          } else if (is_docker_env) {
            "Docker"
          } else {
            "Windows"
          }
          if (self$design$sim_prm$logs) {
            cat(sprintf(
              "Using safe write method for %s (%s environment detected)\n",
              basename(output_path),
              env_type
            ))
          }

          # Use dbGetQuery + arrow::write_parquet instead of DuckDB COPY
          # This avoids the atomic rename issue with SMB/Plan9 mounts
          while (!success && retry_count < max_retries) {
            retry_count <- retry_count + 1
            if (self$design$sim_prm$logs) {
              cat(sprintf(
                "Safe write attempt %d for %s\n",
                retry_count,
                basename(output_path)
              ))
            }

            tryCatch(
              {
                # Get the data directly into R
                result_data <- private$query_sql(duckdb_con, query)
                if (self$design$sim_prm$logs) {
                  cat(sprintf(
                    "Query returned %d rows for %s\n",
                    nrow(result_data),
                    basename(output_path)
                  ))
                }

                if (nrow(result_data) > 0) {
                  # Write directly with arrow, avoiding DuckDB's temporary file mechanism
                  arrow::write_parquet(result_data, output_path)

                  # Verify the write with multiple checks
                  Sys.sleep(0.3) # Allow SMB sync time

                  if (file.exists(output_path)) {
                    file_size <- file.size(output_path)
                    if (file_size > 0) {
                      # Additional verification: try to read back a sample
                      tryCatch(
                        {
                          test_read <- arrow::read_parquet(
                            output_path,
                            as_data_frame = FALSE
                          )
                          if (test_read$num_rows > 0) {
                            success <- TRUE
                            if (self$design$sim_prm$logs) {
                              cat(sprintf(
                                "Successfully wrote %s (%d bytes, %d rows)\n",
                                basename(output_path),
                                file_size,
                                nrow(result_data)
                              ))
                            }
                          }
                        },
                        error = function(e) {
                          cat(sprintf(
                            "File %s exists but failed verification read: %s\n",
                            basename(output_path),
                            conditionMessage(e)
                          ))
                        }
                      )
                    } else {
                      cat(sprintf(
                        "File %s exists but has 0 bytes\n",
                        basename(output_path)
                      ))
                    }
                  }
                } else {
                  # Handle empty results properly
                  if (self$design$sim_prm$logs) {
                    cat(sprintf(
                      "Query returned 0 rows for %s - creating empty file\n",
                      basename(output_path)
                    ))
                  }

                  # Get column structure by limiting to 0 rows
                  empty_structure <- private$query_sql(
                    duckdb_con,
                    sprintf("SELECT * FROM (%s) subq LIMIT 0", query)
                  )
                  arrow::write_parquet(empty_structure, output_path)

                  if (file.exists(output_path)) {
                    success <- TRUE
                    if (self$design$sim_prm$logs) {
                      cat(sprintf(
                        "Created empty parquet file %s\n",
                        basename(output_path)
                      ))
                    }
                  }
                }
              },
              error = function(e) {
                cat(sprintf(
                  "Safe write attempt %d failed for %s: %s\n",
                  retry_count,
                  basename(output_path),
                  conditionMessage(e)
                ))

                # Clean up failed file
                if (file.exists(output_path)) {
                  try(unlink(output_path, recursive = TRUE), silent = TRUE)
                }

                # Progressive backoff
                Sys.sleep(0.5 * retry_count)
              }
            )
          }
        } else {
          # Use original DuckDB COPY method for native Linux environments
          cat(sprintf(
            "Using DuckDB COPY method for %s (native Linux environment)\n",
            basename(output_path)
          ))

          while (!success && retry_count < max_retries) {
            retry_count <- retry_count + 1
            if (self$design$sim_prm$logs) {
              cat(sprintf(
                "DuckDB COPY attempt %d for %s\n",
                retry_count,
                basename(output_path)
              ))
            }

            tryCatch(
              {
                # Use forward slashes for DuckDB COPY command
                db_path <- gsub("\\\\", "/", output_path)

                copy_command <- sprintf(
                  "COPY (%s) TO '%s' (FORMAT PARQUET);",
                  query,
                  db_path
                )
                result <- private$execute_sql(duckdb_con, copy_command)

                Sys.sleep(0.2) # Allow file system sync

                if (file.exists(output_path) && file.size(output_path) > 0) {
                  success <- TRUE
                  if (self$design$sim_prm$logs) {
                    cat(sprintf(
                      "DuckDB COPY succeeded for %s\n",
                      basename(output_path)
                    ))
                  }
                }
              },
              error = function(e) {
                cat(sprintf(
                  "DuckDB COPY attempt %d failed for %s: %s\n",
                  retry_count,
                  basename(output_path),
                  conditionMessage(e)
                ))

                if (file.exists(output_path) && file.size(output_path) == 0) {
                  try(unlink(output_path, recursive = TRUE), silent = TRUE)
                }

                Sys.sleep(0.2 * retry_count)
              }
            )
          }
        }

        # Final validation
        if (!success) {
          final_check_exists <- file.exists(output_path)
          final_check_size <- if (final_check_exists) {
            file.size(output_path)
          } else {
            0
          }

          error_msg <- sprintf(
            "Failed to write output to %s after %d attempts. Final state: exists=%s, size=%d",
            output_path,
            max_retries,
            final_check_exists,
            final_check_size
          )

          if (self$design$sim_prm$logs) {
            cat(paste("ERROR:", error_msg, "\n"))
          }
          stop(error_msg)
        }

        return(invisible(TRUE))
      },

      # Helper function to execute SQL, with error reporting
      # execute_sql ----
      execute_sql = function(con, sql, context = "") {
        tryCatch(
          {
            dbExecute(con, sql)
          },
          error = function(e) {
            stop(paste(
              "Error executing SQL for",
              context,
              ":",
              e$message,
              "\nSQL:\n",
              sql
            ))
          }
        )
      },

      # Helper function to query SQL, with error reporting
      # query_sql ----
      query_sql = function(con, sql, context = "") {
        tryCatch(
          {
            dbGetQuery(con, sql)
          },
          error = function(e) {
            stop(paste(
              "Error querying SQL for",
              context,
              ":",
              e$message,
              "\nSQL:\n",
              sql
            ))
          }
        )
      },

      # profile_sql_view ----
      profile_sql_view = function(con, view_name) {
        private$query_sql(
          con = con,
          sql = paste0("EXPLAIN ANALYZE SELECT * FROM ", view_name),
          context = ""
        )
      },

      # profile_sql_query ----
      profile_sql_query = function(con, sql_query) {
        private$query_sql(
          con = con,
          sql = paste0("EXPLAIN ANALYZE ", sql_query),
          context = ""
        )
      },

      # run_sim ----
      # Runs the simulation in one core. mc is scalar
      run_sim = function(mc_, scenario_nam = "") {
        if (!nzchar(scenario_nam)) {
          scenario_nam <- "sc0"
        }

        # Staggered worker startup to reduce memory pressure during parallel runs.
        # Only apply to the first batch of workers (iterations 1 to clusternumber).
        # After the first batch, workers finish at different times and naturally
        # pick up new iterations in a staggered manner.
        n_cores <- self$design$sim_prm$clusternumber
        if (n_cores > 1L && mc_ <= n_cores) {
          stagger_delay_sec <- 2L # seconds between worker startups
          worker_position <- mc_ - 1L
          if (worker_position > 0L) {
            Sys.sleep(worker_position * stagger_delay_sec)
          }
        }

        if (self$design$sim_prm$logs) {
          private$time_mark(paste0("Start mc iteration ", mc_))
          log_file <- private$output_dir(paste0("logs/log", mc_, ".txt"))
          log_con <- file(log_file, open = "at") # append text
          sink(log_con, type = "output", split = FALSE)
          sink(log_con, type = "message")
        }

        sp <- SynthPop$new(mc_, self$design)
        e <- read_fst("./inputs/mortality/mrtl_clb.fst", as.data.table = TRUE) # mortality calibration
        lookup_dt(
          sp$pop,
          e,
          check_lookup_tbl_validity = self$design$sim_prm$logs
        )
        setnafill(sp$pop, "const", 1, cols = "mrtl_clbr")
        rm(e)

        # make pop weight available to scenarios
        # message("Updating weights")
        if (scenario_nam != "sc0") {
          sp$update_pop_weights(scenario_nam)
        }
        # message("Updating weights finished")

        # Isolate tha rng state for the user defines scenarios
        rs <- .Random.seed
        dqrs <- dqrng_get_state()
        sdn <- digest2int(paste0("primary", scenario_nam), sp$mc_aggr) # sp$mc_aggr ensures same seed for sp batches that stem from the same sp.
        # consequently if the user needs different seeds for different batches
        # they have to explicitly use a new seed, generate rn and then restore
        # the seed.
        set.seed(sdn) # set seed based on scenario name
        dqset.seed(sdn)
        # Note that above does nor guarantee that different scenario_name/sp$mc combination
        # always generate different seed. But the probability of collision is
        # very low. export_summaries() checks if collision happened and warns user.

        private$primary_prevention_scn(sp) # apply primary prevention scenario
        # message("scenario finished")
        set.seed(rs)
        dqrng_set_state(dqrs)

        # gen_parf + set_init_prvl run AFTER the primary-prevention scenario
        # so that initial disease prevalence reflects the scenario world.
        # set_init_prvl internally buffers year == init_year columns with
        # year == init_year - 1 values for the duration of its prevalence
        # calculation; see ?set_init_prvl for the contract.
        lapply(self$design$diseases, function(x) {
          if (self$design$sim_prm$logs) {
            print(x$name)
          }
          x$gen_parf(sp, self$design, self$design$diseases)$set_init_prvl(
            sp = sp,
            design_ = self$design
          )
        })

        lapply(self$design$diseases, function(x) {
          x$set_rr(sp, self$design)$set_incd_prb(sp, self$design)$set_dgns_prb(
            sp,
            self$design
          )$set_mrtl_prb(sp, self$design)
        })
        # message("incd finished")

        # Isolate the rng state for the user defines scenarios
        rs <- .Random.seed
        dqrs <- dqrng_get_state()
        sdn <- digest2int(paste0("secondary", scenario_nam), sp$mc_aggr) # sp$mc_aggr ensures same seed for sp batches that stem from the same sp.
        # consequently if the user needs different seeds for different batches
        # they have to explicitly use a new seed, generate rn and then restore
        # the seed.
        set.seed(sdn) # set seed based on scenario name
        dqset.seed(sdn) # Note that above does nor guarantee that different scenario_name/sp$mc combination
        # always generate different seed. But the probability of collision is
        # very low. export_summaries() checks if collision happened and warns user.
        private$secondary_prevention_scn(sp) # apply secondary prevention scenario
        # message("2nd scenario finished")
        # message("scenario finished")
        set.seed(rs)
        dqrng_set_state(dqrs)

        # ds <- copy(self$diseases) # Necessary for parallelisation
        # lapply(self$diseases, function(x) {
        #   if (self$design$sim_prm$logs) print(x$name)
        #   x$gen_parf(sp, self$design, self$diseases)$
        #     set_init_prvl(sp, self$design)$
        #     set_rr(sp, self$design)$
        #     set_incd_prb(sp, self$design)$
        #     set_dgns_prb(sp, self$design)$
        #     set_mrtl_prb(sp, self$design)
        # })

        l <- private$mk_scenario_init(sp, scenario_nam)
        if (!identical(key(sp$pop), c("pid", "year"))) {
          stop("synthpop key is not as expected")
        }
        simcpp(sp$pop, l, sp$mc)
        # message("cpp finished")
        # it doesn't matter if mc or mc_aggr is used in the above, because it is
        # only used for the RNG stream and the pid are different in each mc_aggr
        # pop

        sp$update_pop_weights(scenario_nam)

        # Prune pop (NOTE that assignment in the function env makes this
        # data.table local)
        sp$pop <- sp$pop[
          all_cause_mrtl >= 0L &
            year >= self$design$sim_prm$init_year &
            between(age, self$design$sim_prm$ageL, self$design$sim_prm$ageH),
        ]
        setkey(sp$pop, pid, year)
        sp$pop[, pid_mrk := mk_new_simulant_markers(pid)]

        # apply ESP weights
        to_agegrp(sp$pop, 5, 99)
        absorb_dt(sp$pop, private$esp_weights, exclude_col = "wt_esp")
        sp$pop[,
          wt_esp := wt_esp * unique(wt_esp) / sum(wt_esp),
          by = .(year, agegrp, sex, dimd)
        ] # NOTE keyby changes the key

        if (self$design$sim_prm$export_xps) {
          if (self$design$sim_prm$logs) {
            message("Exporting exposures...")
          }
          private$export_xps(sp, scenario_nam)
        }

        # Columns kept in the lifecourse output are those listed in
        # `cols_for_output` (YAML) OR whose name matches one of the supported
        # suffix conventions below. Users can create custom columns inside a
        # scenario function (`synthpop$pop[, my_col := ...]`) and have them flow
        # into summaries automatically by naming them with the right suffix:
        #   *_prvl  -> duration counter: prevalence is SUM(wt) where col > 0
        #              (export_prvl_summaries) AND incidence is SUM(wt) where
        #              col = 1 (export_incd_summaries). One column drives both,
        #              mirroring the disease *_prvl columns.
        #   *_contd -> continuous:  population-weighted mean (export_contd_summaries)
        #   *_costs -> economic cost: SUM(col * wt)          (export_costs_summaries)
        #   *_dgns  -> diagnosis flag; cms_* / *_mrtl -> kept for internal use.
        # See vignette "custom-scenario-columns" for the full mechanism.
        nam <- c(
          self$design$sim_prm$cols_for_output,
          grep("^cms_|_prvl$|_contd$|_costs$|_dgns$|_mrtl$", names(sp$pop), value = TRUE)
        )
        nam <- grep("^prb_", nam, value = TRUE, invert = TRUE) # exclude prb_ ... _dgns
        sp$pop[, setdiff(names(sp$pop), nam) := NULL]
        sp$pop[, `:=`(mc = sp$mc_aggr, mc_chunk = sp$mc)]

        # TODO add logic for the years of having MM. Currently 1 is not the real
        # incidence. It is still prevalence
        sp$pop[, `:=`(
          cms1st_cont_prvl = carry_forward_incr(
            as.integer(cms_count == 1),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          ),
          cmsmm0_prvl = carry_forward_incr(
            as.integer(cms_score > 0),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          ),
          cmsmm1_prvl = carry_forward_incr(
            as.integer(cms_score > 1),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          ),
          cmsmm1.5_prvl = carry_forward_incr(
            as.integer(cms_score > 1.5),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          ),
          cmsmm2_prvl = carry_forward_incr(
            as.integer(cms_score > 2),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          ),
          cmscs1_prvl = carry_forward_incr(
            as.integer(cms_score < 0.09),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          ),
          cmscs2_prvl = carry_forward_incr(
            as.integer(cms_score >= 0.09 & cms_score < 0.69),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          ),
          cmscs3_prvl = carry_forward_incr(
            as.integer(cms_score >= 0.69 & cms_score < 1.59),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          ),
          cmscs4_prvl = carry_forward_incr(
            as.integer(cms_score >= 1.59 & cms_score < 2.95),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          ),
          cmscs5_prvl = carry_forward_incr(
            as.integer(cms_score >= 2.95),
            pid_mrk,
            TRUE,
            1L,
            byref = TRUE
          )
        )]

        # sp$pop[, scenario := scenario_nam]

        setkeyv(sp$pop, c("pid", "year"))

        # Write lifecourse
        if (self$design$sim_prm$logs) {
          message("Exporting lifecourse...")
        }

        fileformat <- "parquet"
        fnam <- private$output_dir(file.path(
          "lifecourse",
          paste0("mc=", sp$mc_aggr),
          paste0("scenario=", scenario_nam),
          paste0(sp$mc, "_lifecourse.", fileformat)
        ))
        # NOTE parquet format about 30 times smaller but about 50% slower in writting to disk
        write_dataset(dataset = sp$pop, path = fnam, format = fileformat)

        if (self$design$sim_prm$logs) {
          private$time_mark(paste0("End mc iteration ", mc_))
          sink(type = "message") # close message sink first
          sink(type = "output") # close output sink
          if (
            exists("log_con") &&
              inherits(log_con, "connection") &&
              isOpen(log_con)
          ) {
            close(log_con)
          }
        }

        NULL
      },

      # creates the list that is used in c++ side sp is needed for sp$mc_aggr in
      # to_cpp()

      # mk_scenario_init ----
      mk_scenario_init = function(sp, scenario_name) {
        # scenario_suffix_for_pop <- paste0("_", scenario_name)

        # TODO the next line counteracts the commented line above. This is
        # intentional until we finalise the scenario mechanism
        scenario_suffix_for_pop <- ""

        list(
          "exposures" = self$design$sim_prm$exposures_for_output,
          "scenarios" = self$design$sim_prm$scenarios, # to be generated programmatically
          "scenario" = scenario_name,
          "kismet" = self$design$sim_prm$kismet, # If TRUE random numbers are the same for each scenario.
          "init_year" = self$design$sim_prm$init_year,
          "pids" = "pid",
          "years" = "year",
          "ages" = "age",
          "sexs" = "sex",
          "dimds" = "dimd",
          "ageL" = self$design$sim_prm$ageL,
          "all_cause_mrtl" = paste0("all_cause_mrtl", scenario_suffix_for_pop),
          "cms_score" = paste0("cms_score", scenario_suffix_for_pop),
          "cms_count" = paste0("cms_count", scenario_suffix_for_pop),
          # "strata_for_outputs" = c("pid", "year", "age", "sex", "dimd"),
          "diseases" = lapply(self$design$diseases, function(x) {
            x$to_cpp(sp, self$design, scenario_name, scenario_suffix_for_pop)
          })
        )
      },
      # export_xps ----
      export_xps = function(sp, scenario_nam) {
        # NOTE no need to check validity of inputs here as it is only used
        # internally

        logs_enabled <- self$design$sim_prm$logs

        tryCatch(
          {
            if (logs_enabled) {
              message("  export_xps: Starting age grouping...")
            }

            to_agegrp(
              sp$pop,
              grp_width = 20L,
              max_age = self$design$sim_prm$ageH,
              min_age = self$design$sim_prm$ageL,
              age_colname = "age",
              agegrp_colname = "agegrp20",
              to_factor = TRUE
            )

            sp$pop[,
              smok_never_curr_xps := fifelse(
                smok_status_curr_xps == "1",
                1L,
                0L
              )
            ]
            sp$pop[,
              smok_active_curr_xps := fifelse(
                smok_status_curr_xps == "4",
                1L,
                0L
              )
            ]

            xps <- grep("_curr_xps$", names(sp$pop), value = TRUE)
            xps <- grep("_prvl_curr_xps$", xps, value = TRUE, invert = TRUE)
            # NOTE: Using setdiff() instead of -which() to avoid R gotcha where
            # x[-integer(0)] returns empty vector instead of x
            xps <- setdiff(
              xps,
              c("smok_status_curr_xps", "met_curr_xps", "bpmed_curr_xps")
            )

            # Add any extra columns from exposures_for_output that exist in sp$pop
            extra_xps <- intersect(
              self$design$sim_prm$exposures_for_output,
              names(sp$pop)
            )
            # Only keep numeric columns (weighted.mean requires numeric input)
            extra_xps <- extra_xps[
              vapply(extra_xps, function(col) is.numeric(sp$pop[[col]]), logical(1))
            ]
            xps <- unique(c(xps, extra_xps))

            # Defensive check: ensure xps is not empty
            if (length(xps) == 0L) {
              warning(
                "export_xps: No exposure columns found matching '_curr_xps$' pattern. ",
                "Available columns: ",
                paste(head(names(sp$pop), 20), collapse = ", "),
                if (length(names(sp$pop)) > 20) "..." else ""
              )
              return(NULL)
            }
            if (logs_enabled) {
              message(
                "  export_xps: Found ",
                length(xps),
                " exposure columns: ",
                paste(head(xps, 5), collapse = ", "),
                if (length(xps) > 5) "..." else ""
              )
            }

            sp$pop[
              smok_status_curr_xps == "1",
              `:=`(
                smok_packyrs_curr_xps = NA,
                smok_quit_yrs_curr_xps = NA,
                smok_dur_curr_xps = NA,
                smok_cig_curr_xps = NA
              )
            ]
            sp$pop[
              smok_status_curr_xps == "4",
              `:=`(
                smok_quit_yrs_curr_xps = NA
              )
            ]

            # Build xps stratification variables dynamically from strata_for_output
            # Mapping: agegrp -> agegrp20 (20-year bands for xps), dimd -> qimd
            xps_strata_raw <- setdiff(
              self$design$sim_prm$strata_for_output,
              c("scenario", "year")
            )
            xps20_var_map <- c("agegrp" = "agegrp20", "dimd" = "qimd")
            xps20_strata <- vapply(
              xps_strata_raw,
              function(v) {
                if (v %in% names(xps20_var_map)) xps20_var_map[[v]] else v
              },
              character(1),
              USE.NAMES = FALSE
            )
            # Generate all subset combinations for groupingsets
            xps20_subsets <- unlist(
              lapply(seq_along(xps20_strata), function(k) {
                combn(xps20_strata, k, simplify = FALSE)
              }),
              recursive = FALSE
            )
            xps20_sets <- c(
              list("year"),
              lapply(xps20_subsets, function(s) c("year", s))
            )
            xps20_by <- c("year", xps20_strata)

            if (logs_enabled) {
              message("  export_xps: Computing groupingsets for xps20...")
            }
            out_xps20 <- groupingsets(
              sp$pop[
                all_cause_mrtl >= 0L &
                  year >= self$design$sim_prm$init_year &
                  age >= self$design$sim_prm$ageL,
              ],
              j = lapply(.SD, weighted.mean, wt, na.rm = TRUE),
              by = xps20_by,
              .SDcols = xps,
              sets = xps20_sets
            )[, `:=`(year = year + 2000L, mc = sp$mc, scenario = scenario_nam)]
            # TODO above mc could also be mc_aggr. Getting the uncertainty right here is tricky

            for (j in seq_len(ncol(out_xps20))) {
              set(out_xps20, which(is.na(out_xps20[[j]])), j, "All")
            }
            setkey(out_xps20, year)

            fileformat <- "parquet"
            fnam <- private$output_dir(file.path(
              "xps",
              "xps20",
              paste0("mc=", sp$mc_aggr),
              paste0("scenario=", scenario_nam),
              paste0(sp$mc, "_xps20.", fileformat)
            ))

            # Ensure parent directory exists
            dir.create(dirname(fnam), recursive = TRUE, showWarnings = FALSE)

            if (logs_enabled) {
              message("  export_xps: Writing xps20 to ", fnam)
            }
            # NOTE parquet format about 30 times smaller but about 50% slower in writting to disk
            write_dataset(dataset = out_xps20, path = fnam, format = fileformat)

            # xps5 (ESP): same strata but without age groups (standardised by age)
            xps5_strata <- setdiff(xps20_strata, "agegrp20")
            xps5_subsets <- unlist(
              lapply(seq_along(xps5_strata), function(k) {
                combn(xps5_strata, k, simplify = FALSE)
              }),
              recursive = FALSE
            )
            xps5_sets <- c(
              list("year"),
              lapply(xps5_subsets, function(s) c("year", s))
            )
            xps5_by <- c("year", xps5_strata)

            if (logs_enabled) {
              message("  export_xps: Computing groupingsets for xps5...")
            }
            out_xps5 <- groupingsets(
              sp$pop[
                all_cause_mrtl >= 0L &
                  year >= self$design$sim_prm$init_year &
                  age >= self$design$sim_prm$ageL,
              ],
              j = lapply(.SD, weighted.mean, wt_esp, na.rm = TRUE),
              by = xps5_by,
              .SDcols = xps,
              sets = xps5_sets
            )[, `:=`(year = year + 2000L, mc = sp$mc, scenario = scenario_nam)]
            for (j in seq_len(ncol(out_xps5))) {
              set(out_xps5, which(is.na(out_xps5[[j]])), j, "All")
            }
            setkey(out_xps5, year)

            fnam <- private$output_dir(file.path(
              "xps",
              "xps5",
              paste0("mc=", sp$mc_aggr),
              paste0("scenario=", scenario_nam),
              paste0(sp$mc, "_xps_esp.", fileformat)
            ))

            # Ensure parent directory exists
            dir.create(dirname(fnam), recursive = TRUE, showWarnings = FALSE)

            if (logs_enabled) {
              message("  export_xps: Writing xps5 to ", fnam)
            }
            # NOTE parquet format about 30 times smaller but about 50% slower in writting to disk
            write_dataset(dataset = out_xps5, path = fnam, format = fileformat)

            if (logs_enabled) {
              message("  export_xps: Tidying up temporary columns...")
            }
            # Tidy up
            sp$pop[,
              c(
                "agegrp20",
                "smok_never_curr_xps",
                "smok_active_curr_xps"
              ) := NULL
            ]
            sp$pop[
              smok_status_curr_xps == "1",
              `:=`(
                smok_packyrs_curr_xps = 0,
                smok_quit_yrs_curr_xps = 0,
                smok_dur_curr_xps = 0,
                smok_cig_curr_xps = 0
              )
            ]
            sp$pop[
              smok_status_curr_xps == "4",
              `:=`(
                smok_quit_yrs_curr_xps = 0
              )
            ]

            if (logs_enabled) {
              message("  export_xps: Completed successfully")
            }
            NULL
          },
          error = function(e) {
            # Log the error with context for debugging
            err_msg <- paste0(
              "export_xps FAILED for mc=",
              sp$mc,
              ", scenario=",
              scenario_nam,
              "\n  Error: ",
              conditionMessage(e),
              "\n  Call: ",
              deparse(conditionCall(e))
            )
            warning(err_msg, immediate. = TRUE)
            # Re-throw the error so it propagates
            stop(e)
          }
        )
      },

      # Function for timing log
      time_mark = function(text_id) {
        sink(
          file = private$output_dir("logs/times.txt"),
          append = TRUE,
          type = "output",
          split = FALSE
        )
        cat(paste0(text_id, " at: ", Sys.time(), "\n"))
        sink()
      },
      output_dir = function(x = "") {
        file.path(self$design$sim_prm$output_dir, x)
      },

      # deep_clone ----
      # Special deep copy for data.table. Use POP$clone(deep = TRUE) to
      # dispatch. Otherwise a reference is created
      deep_clone = function(name, value) {
        if ("data.table" %in% class(value)) {
          data.table::copy(value)
        } else if ("R6" %in% class(value)) {
          value$clone()
        } else {
          # For everything else, just return it. This results in a shallow copy
          # of s3.
          value
        }
      },

      # create_new_folder ----
      # @description Create folder if doesn't exist. Stops on failure.
      # @param sDirPathName String folder path and name.
      # @param bReport Bool report folder creation. Default is design$sim_prm$logs.
      create_new_folder = function(
        sDirPathName,
        bReport = self$design$sim_prm$logs
      ) {
        # Normalize path for cross-platform compatibility
        sDirPathName <- normalizePath(sDirPathName, mustWork = FALSE)

        if (!dir.exists(sDirPathName)) {
          # Try creating directory with retries for Windows compatibility
          max_retries <- 3
          retry_count <- 0
          bSuccess <- FALSE

          while (!bSuccess && retry_count < max_retries) {
            retry_count <- retry_count + 1

            tryCatch(
              {
                bSuccess <- dir.create(sDirPathName, recursive = TRUE)

                # If creation returned FALSE, check if it exists (race condition handling)
                if (!bSuccess && dir.exists(sDirPathName)) {
                  bSuccess <- TRUE
                }

                # Add small delay for Windows file system sync
                if (.Platform$OS.type == "windows") {
                  Sys.sleep(0.05)
                }

                # Verify directory was actually created
                if (!dir.exists(sDirPathName)) {
                  bSuccess <- FALSE
                }
              },
              error = function(e) {
                cat(sprintf(
                  "Directory creation attempt %d failed: %s\n",
                  retry_count,
                  conditionMessage(e)
                ))
                bSuccess <- FALSE
              }
            )

            if (!bSuccess && retry_count < max_retries) {
              Sys.sleep(0.1 * retry_count) # Progressive backoff
            }
          }

          if (!bSuccess) {
            stop(paste0(
              "Failed creating directory ",
              sDirPathName,
              " after ",
              max_retries,
              " attempts. Check permissions."
            ))
          } else {
            if (bReport) {
              message(paste0("Folder ", sDirPathName, " was created"))
            }
          }
        }

        # Check if directory is writable
        if (file.access(sDirPathName, mode = 2) != 0) {
          stop(paste0(
            "Directory ",
            sDirPathName,
            " exists but is not writable. Check permissions."
          ))
        }

        invisible(TRUE)
      },

      # create_output_folder_structure ----
      # Create output folder structure
      create_output_folder_structure = function() {
        if (self$design$sim_prm$logs) {
          message("Creating output subfolders.")
        }
        if (
          file.exists(self$design$sim_prm$output_dir) &&
            file.access(self$design$sim_prm$output_dir, mode = 2) == -1L
        ) {
          stop(
            "You don't have write access to the output folder. Please change the permissions or the path. If you are using Linux you can use i.e. IMPACTncd$allow_universal_output_folder_access() to allow write access to the output folder."
          )
        }
        if (
          file.exists(self$design$sim_prm$synthpop_dir) &&
            file.access(self$design$sim_prm$synthpop_dir, mode = 2) == -1L
        ) {
          stop(
            "You don't have write access to the synthpop folder. Please change the permissions or the path.  If you are using Linux you can use i.e. IMPACTncd$allow_universal_synthpop_folder_access() to allow write access to the synthpop folder."
          )
        }
        private$create_new_folder(
          self$design$sim_prm$output_dir,
          self$design$sim_prm$logs
        )
        private$create_new_folder(
          private$output_dir("summaries/"),
          self$design$sim_prm$logs
        )
        private$create_new_folder(
          private$output_dir("tables/"),
          self$design$sim_prm$logs
        )
        private$create_new_folder(
          private$output_dir("plots/"),
          self$design$sim_prm$logs
        )
        private$create_new_folder(
          private$output_dir("lifecourse/"),
          self$design$sim_prm$logs
        )
        if (self$design$sim_prm$export_PARF) {
          private$create_new_folder(
            private$output_dir("parf/"),
            self$design$sim_prm$logs
          )
        }
        if (self$design$sim_prm$export_xps) {
          private$create_new_folder(
            private$output_dir("xps/"),
            self$design$sim_prm$logs
          )
        }
        if (self$design$sim_prm$logs) {
          private$create_new_folder(
            private$output_dir("logs/"),
            self$design$sim_prm$logs
          )
        }
        # NOTE code below is duplicated in Synthpop class. This is intentional
        private$create_new_folder(
          self$design$sim_prm$synthpop_dir,
          self$design$sim_prm$logs
        )
      }
    ) # End of private methods
  ) # End of class
