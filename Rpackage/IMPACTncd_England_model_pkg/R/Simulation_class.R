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



# From
# https://stackoverflow.com/questions/33424233/how-do-i-tell-an-r6-class-what-to-do-with-square-brackets
# Allows data.table syntax to the R6class object directly. Assumes it has a
# field 'output' that is a data.table

#' @export
`[.Simulation` <- function(x, ...) x$output[...]

#' R6 Class representing a simulation environment
#' @description A simulation environment.
#' @details To be completed...
#'
#' @section Methods defined in satellite files:
#' The following public methods are implemented in separate files for maintainability
#' and added to this class via \code{$set()}:
#'
#' \describe{
#'   \item{\code{calibrate_incd_ftlt(mc, replace = FALSE)}}{
#'     Calibrates incidence and case fatality rates. See Simulation_class_calibration.R.
#'   }
#'   \item{\code{export_summaries(multicore, type, single_year_of_age)}}{
#'     Exports simulation summaries. See Simulation_class_summaries.R.
#'   }
#'   \item{\code{export_tables(baseline_year_for_change_outputs, prbl, comparator_scenario, two_agegrps)}}{
#'     Exports summary tables. See Simulation_class_tables.R.
#'   }
#' }
#' @export
Simulation <-
  R6::R6Class(
    classname = "Simulation",
    lock_objects = TRUE, # allows primary prevention scenario to be updated
    lock_class = FALSE,  # allows adding methods via $set() from other files
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
        
        # private$create_empty_calibration_prms_file(replace = FALSE)
        message("Generating microsimulation structure.")
        # Generate the graph with the causality structure
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

        private$death_codes <- unlist(lapply(self$design$diseases, function(x) {
          x$meta$mortality$code
        }))
        private$death_codes[["alive"]] <- 0L

        private$primary_prevention_scn <- function(synthpop) NULL # default for baseline scenario
        private$secondary_prevention_scn <- function(synthpop) NULL # default for baseline scenario

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
      #' @param scenario_nam A string for the scenario name (i.e. sc1)
      #' @return The invisible self for chaining.
      run = function(mc, multicore = TRUE, scenario_nam) {
        if (!is.integer(mc)) {
          stop("mc need to be an integer")
        }
        if (any(mc <= 0)) {
          stop("mc need to be positive integer")
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

        # TODO better logic as this is always true for the non baseline scenario
        if (
          any(file.exists(
            # TODO fix when lifecourse is not saved
            file.path(
              self$design$sim_prm$output_dir,
              "lifecourse"
            ),
            recursive = TRUE
          ))
        ) {
          # stop("Results from a previous simulation exists in the output
          #      folder. Please remove them before run a new one.")
          if (self$design$sim_prm$logs) {
          message(
              "Results from a previous simulation exists in the output folder. Please remove them if this was unintentional."
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
          print(
            plot(
              graph,
              vertex.shape = "none",
              edge.arrow.size = .3,
              vertex.label.font = 2,
              vertex.label.color = "gray40",
              edge.arrow.width = .5,
              vertex.label.cex = .7,
              edge.color = "gray85",
              layout = layout_components
            )
          )
        }

        if (processed) {
          graph <- as.matrix(as_adjacency_matrix(graph))
          n <- vapply(self$design$diseases, `[[`, "name", FUN.VALUE = character(1))
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
                aes(x = year, y = get(paste0(metric, "_mrtl_rate")), color = type)
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(x = year, y = get(paste0(metric, "_mrtl_rate_low")), color = type),
                linetype = "dashed"
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(x = year, y = get(paste0(metric, "_mrtl_rate_upp")), color = type),
                linetype = "dashed"
              ) +
              facet_wrap(. ~ factor(agegrp), scales = "free") +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
              ggtitle(paste0(toupper(metric), " mrtl rate"), tools::toTitleCase(sex_val))
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
            geom_line(data = dt, aes(x = year, y = get(paste0(metric, "_mrtl_rate")), color = type)) +
            geom_line(
              data = dt,
              aes(x = year, y = get(paste0(metric, "_mrtl_rate_low")), color = type),
              linetype = "dashed"
            ) +
            geom_line(
              data = dt,
              aes(x = year, y = get(paste0(metric, "_mrtl_rate_upp")), color = type),
              linetype = "dashed"
            ) +
            facet_wrap(. ~ factor(sex), scales = "free") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
            ggtitle(paste0(toupper(metric), " mrtl rate"), "By sex")
          ggsave(
            file.path(self$design$sim_prm$output_dir, "plots", paste0(toupper(metric), "_s_mrtl.jpg")),
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
                aes(x = year, y = get(paste0(metric, "_incd_rate")), color = type)
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(x = year, y = get(paste0(metric, "_incd_rate_low")), color = type),
                linetype = "dashed"
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(x = year, y = get(paste0(metric, "_incd_rate_upp")), color = type),
                linetype = "dashed"
              ) +
              facet_wrap(. ~ factor(agegrp), scales = "free") +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
              ggtitle(paste0(toupper(metric), " incd rate"), tools::toTitleCase(sex_val))
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
            geom_line(data = dt, aes(x = year, y = get(paste0(metric, "_incd_rate")), color = type)) +
            geom_line(
              data = dt,
              aes(x = year, y = get(paste0(metric, "_incd_rate_low")), color = type),
              linetype = "dashed"
            ) +
            geom_line(
              data = dt,
              aes(x = year, y = get(paste0(metric, "_incd_rate_upp")), color = type),
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
                aes(x = year, y = get(paste0(metric, "_prvl_rate")), color = type)
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(x = year, y = get(paste0(metric, "_prvl_rate_low")), color = type),
                linetype = "dashed"
              ) +
              geom_line(
                data = dt[sex == sex_val],
                aes(x = year, y = get(paste0(metric, "_prvl_rate_upp")), color = type),
                linetype = "dashed"
              ) +
              facet_wrap(. ~ factor(agegrp), scales = "free") +
              theme_bw() +
              theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) +
              ggtitle(paste0(toupper(metric), " prvl rate"), tools::toTitleCase(sex_val))
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
            geom_line(data = dt, aes(x = year, y = get(paste0(metric, "_prvl_rate")), color = type)) +
            geom_line(
              data = dt,
              aes(x = year, y = get(paste0(metric, "_prvl_rate_low")), color = type),
              linetype = "dashed"
            ) +
            geom_line(
              data = dt,
              aes(x = year, y = get(paste0(metric, "_prvl_rate_upp")), color = type),
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
      #' @description Connect to Zenodo for asset management.
      #' @param token Zenodo personal access token. If NULL, reads from ZENODO_TOKEN env var.
      #' @param concept_doi Optional concept DOI for the data record.
      #' @param sandbox Use Zenodo sandbox for testing.
      #' @return The invisible self for chaining.
      #' @details This method initializes the Zenodo asset manager for uploading/downloading
      #'   input data files. Get your token at: https://zenodo.org/account/settings/applications/tokens/new/
      zenodo_connect = function(token = NULL, concept_doi = NULL, sandbox = FALSE) {
        if (is.null(private$zenodo_manager)) {
          private$zenodo_manager <- ZenodoAssetManager$new(
            hash_file = "./simulation/zenodo_manifest.csv",
            logs = self$design$sim_prm$logs,
            sandbox = sandbox
          )
        }

        private$zenodo_manager$connect(token = token, sandbox = sandbox)

        if (!is.null(concept_doi)) {
          private$zenodo_manager$set_concept_doi(concept_doi)
        }

        invisible(self)
      },

      # zenodo_set_doi ----
      #' @description Set the Zenodo concept DOI for data assets.
      #' @param concept_doi The concept DOI (shared across all versions of a record).
      #' @return The invisible self for chaining.
      zenodo_set_doi = function(concept_doi) {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }
        private$zenodo_manager$set_concept_doi(concept_doi)
        invisible(self)
      },

      # zenodo_check_inputs ----
      #' @description Check sync status of local inputs with Zenodo.
      #' @param input_base Base directory containing inputs (default: "./inputs").
      #' @param directories Specific subdirectories to check. If NULL, checks all.
      #' @return A data.table with sync status information.
      #' @details Compares local input directories against the Zenodo record to identify
      #'   missing or potentially modified files. Requires prior call to zenodo_connect()

      #'   and zenodo_set_doi().
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
      #' @description Download input data from Zenodo.
      #' @param input_base Base directory for inputs (default: "./inputs").
      #' @param directories Specific subdirectories to download. If NULL, downloads all.
      #' @param overwrite If TRUE, overwrite existing local directories.
      #' @return The invisible self for chaining.
      #' @details Downloads and extracts input data archives from the Zenodo record.
      #'   Only downloads missing directories unless overwrite=TRUE.
      #'   Use zenodo_check_inputs() first to see what would be downloaded.
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

        # Check current status
        status <- private$zenodo_manager$sync_inputs(
          input_base = input_base,
          directories = directories,
          action = "check"
        )

        # Warn about existing directories
        if (any(status$local_exists) && !overwrite) {
          existing_dirs <- status[local_exists == TRUE, directory]
          message(
            "The following directories already exist and will be SKIPPED:\n  ",
            paste(existing_dirs, collapse = "\n  "),
            "\nUse overwrite=TRUE to replace them."
          )
        }

        # Download missing/all
        private$zenodo_manager$sync_inputs(
          input_base = input_base,
          directories = directories,
          action = "download",
          overwrite = overwrite
        )

        if (self$design$sim_prm$logs) {
          message("Input download complete.")
        }

        invisible(self)
      },

      # zenodo_upload_inputs ----
      #' @description Upload input data to Zenodo.
      #' @param input_base Base directory containing inputs (default: "./inputs").
      #' @param directories Specific subdirectories to upload. If NULL, uploads all.
      #' @param version Version string for the new upload (required for new versions).
      #' @param title Record title (required for new records).
      #' @param description Record description (required for new records).
      #' @param creators List of creators for new records.
      #' @param publish If TRUE, publish the record after upload. Use with caution!
      #' @return The invisible self for chaining.
      #' @details Creates zip archives of input directories and uploads them to Zenodo.
      #'   For new records, provide title, description, and creators.
      #'   For new versions of existing records, just provide version.
      #'   The record remains in DRAFT state unless publish=TRUE.
      zenodo_upload_inputs = function(
        input_base = "./inputs",
        directories = NULL,
        version = NULL,
        title = NULL,
        description = NULL,
        creators = NULL,
        publish = FALSE
      ) {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }

        # Check if we need to create a new record or a new version
        has_doi <- !is.null(private$zenodo_manager$concept_doi)

        if (!has_doi) {
          # Creating a new record
          if (is.null(title) || is.null(description) || is.null(creators)) {
            stop(
              "For new records, you must provide: title, description, and creators.\n",
              "Example:\n",
              "  creators = list(\n",
              "    list(firstname = 'Chris', lastname = 'Kypridemos', orcid = '0000-0002-0746-9229')\n",
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
            version = version %||% "1.0.0"
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

        # Upload the archives
        if (self$design$sim_prm$logs) {
          message("Uploading input archives to Zenodo...")
        }

        archives <- private$zenodo_manager$sync_inputs(
          input_base = input_base,
          directories = directories,
          action = "upload",
          version = version
        )

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

      # zenodo_publish ----
      #' @description Publish the current Zenodo record.
      #' @return The invisible self for chaining.
      #' @details CAUTION: Publishing is irreversible. Once published, the DOI is
      #'   permanent and files cannot be modified (only new versions can be created).
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

      # zenodo_get_versions ----
      #' @description Get all versions of the Zenodo record.
      #' @return A data.frame with version information (date, version, DOI).
      zenodo_get_versions = function() {
        if (is.null(private$zenodo_manager)) {
          stop("Not connected to Zenodo. Call zenodo_connect() first.")
        }
        private$zenodo_manager$get_versions()
      },

      # zenodo_create_manifest ----
      #' @description Create a hash manifest for input files.
      #' @param input_base Base directory containing inputs.
      #' @param save If TRUE, save the manifest to disk.
      #' @return A data.table with file hashes.
      #' @details Creates a manifest file containing hashes of all input files.
      #'   This is used to detect changes between local and remote files.
      zenodo_create_manifest = function(input_base = "./inputs", save = TRUE) {
        if (is.null(private$zenodo_manager)) {
          # Create manager without connection for local operations
          private$zenodo_manager <- ZenodoAssetManager$new(
            hash_file = "./simulation/zenodo_manifest.csv",
            logs = self$design$sim_prm$logs
          )
        }

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
                          warning(sprintf(
                            "File %s exists but failed verification read: %s",
                            basename(output_path),
                            e$message
                          ))
                        }
                      )
                    } else {
                      warning(sprintf(
                        "File %s exists but has 0 bytes",
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
                warning(sprintf(
                  "Safe write attempt %d failed for %s: %s",
                  retry_count,
                  basename(output_path),
                  e$message
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
                warning(sprintf(
                  "DuckDB COPY attempt %d failed for %s: %s",
                  retry_count,
                  basename(output_path),
                  e$message
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
          stagger_delay_sec <- 2L  # seconds between worker startups
          worker_position <- mc_ - 1L
          if (worker_position > 0L) {
            Sys.sleep(worker_position * stagger_delay_sec)
          }
        }

        if (self$design$sim_prm$logs) {
          private$time_mark(paste0("Start mc iteration ", mc_))
          log_file <- private$output_dir(paste0("logs/log", mc_, ".txt"))
          log_con <- file(log_file, open = "at")  # append text
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

        # From Karl, somehow scenario_fn() makes init_prvl different. The
        # following code solves the problem.
        # TODO: investigate the root cause

        lapply(self$design$diseases, function(x) {
          if (self$design$sim_prm$logs) {
            print(x$name)
          }
          x$gen_parf(sp, self$design, self$design$diseases)$set_init_prvl(
            sp = sp,
            design_ = self$design
          )
        })

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

        lapply(self$design$diseases, function(x) {
          x$set_rr(sp, self$design)$set_incd_prb(sp, self$design)$set_dgns_prb(
            sp,
            self$design
          )$set_mrtl_prb(sp, self$design)
        })
        # message("incd finished")

        # Isolate tha rng state for the user defines scenarios
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
        private$secondary_prevention_scn(sp) # apply secondary pevention scenario
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

        nam <- c(
          self$design$sim_prm$cols_for_output,
          grep("^cms_|_prvl$|_dgns$|_mrtl$", names(sp$pop), value = TRUE)
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
          sink(type = "message")  # close message sink first
          sink(type = "output")   # close output sink
          if (exists("log_con") && inherits(log_con, "connection") && isOpen(log_con)) {
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

        tryCatch({
          if (logs_enabled) message("  export_xps: Starting age grouping...")

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
          smok_never_curr_xps := fifelse(smok_status_curr_xps == "1", 1L, 0L)
        ]
        sp$pop[,
          smok_active_curr_xps := fifelse(smok_status_curr_xps == "4", 1L, 0L)
        ]

        xps <- grep("_curr_xps$", names(sp$pop), value = TRUE)
        xps <- grep("_prvl_curr_xps$", xps, value = TRUE, invert = TRUE)
        # NOTE: Using setdiff() instead of -which() to avoid R gotcha where
        # x[-integer(0)] returns empty vector instead of x
        xps <- setdiff(xps, c("smok_status_curr_xps", "met_curr_xps", "bpmed_curr_xps"))

        # Defensive check: ensure xps is not empty
        if (length(xps) == 0L) {
          warning("export_xps: No exposure columns found matching '_curr_xps$' pattern. ",
                  "Available columns: ", paste(head(names(sp$pop), 20), collapse = ", "),
                  if (length(names(sp$pop)) > 20) "..." else "")
          return(NULL)
        }
        if (logs_enabled) {
          message("  export_xps: Found ", length(xps), " exposure columns: ",
                  paste(head(xps, 5), collapse = ", "),
                  if (length(xps) > 5) "..." else "")
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

        if (logs_enabled) message("  export_xps: Computing groupingsets for xps20...")
        out_xps20 <- groupingsets(
          sp$pop[
            all_cause_mrtl >= 0L &
              year >= self$design$sim_prm$init_year &
              age >= self$design$sim_prm$ageL,
          ],
          j = lapply(.SD, weighted.mean, wt, na.rm = TRUE),
          by = c("year", "sex", "agegrp20", "qimd"), # "ethnicity", "sha"
          .SDcols = xps,
          sets = list(
            "year",
            c("year", "agegrp20"),
            c("year", "sex"),
            c("year", "qimd"),
            c("year", "agegrp20", "sex"),
            c("year", "sex", "agegrp20", "qimd")

            # c("year", "ethnicity"),
            # c("year", "sha")
          )
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

        if (logs_enabled) message("  export_xps: Writing xps20 to ", fnam)
        # NOTE parquet format about 30 times smaller but about 50% slower in writting to disk
        write_dataset(dataset = out_xps20, path = fnam, format = fileformat)

        # TODO link strata in the outputs to the design.yaml
        if (logs_enabled) message("  export_xps: Computing groupingsets for xps5...")
        out_xps5 <- groupingsets(
          sp$pop[
            all_cause_mrtl >= 0L &
              year >= self$design$sim_prm$init_year &
              age >= self$design$sim_prm$ageL,
          ],
          j = lapply(.SD, weighted.mean, wt_esp, na.rm = TRUE),
          by = c("year", "sex", "qimd"), # "ethnicity", "sha"
          .SDcols = xps,
          sets = list(
            "year",
            c("year", "sex"),
            c("year", "qimd"),
            c("year", "sex", "qimd")
            # c("year", "ethnicity"),
            # c("year", "sha")
          )
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

        if (logs_enabled) message("  export_xps: Writing xps5 to ", fnam)
        # NOTE parquet format about 30 times smaller but about 50% slower in writting to disk
        write_dataset(dataset = out_xps5, path = fnam, format = fileformat)




        if (logs_enabled) message("  export_xps: Tidying up temporary columns...")
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

        if (logs_enabled) message("  export_xps: Completed successfully")
        NULL

        }, error = function(e) {
          # Log the error with context for debugging
          err_msg <- paste0(
            "export_xps FAILED for mc=", sp$mc, ", scenario=", scenario_nam,
            "\n  Error: ", conditionMessage(e),
            "\n  Call: ", deparse(conditionCall(e))
          )
          warning(err_msg, immediate. = TRUE)
          # Re-throw the error so it propagates
          stop(e)
        })
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
                warning(sprintf(
                  "Directory creation attempt %d failed: %s",
                  retry_count,
                  e$message
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
        if (self$design$sim_prm$logs) message("Creating output subfolders.")
        if (file.exists(self$design$sim_prm$output_dir) && file.access(self$design$sim_prm$output_dir, mode = 2) == -1L) {
          stop("You don't have write access to the output folder. Please change the permissions or the path. If you are using Linux you can use i.e. IMPACTncd$allow_universal_output_folder_access() to allow write access to the output folder.")
        }
        if (file.exists(self$design$sim_prm$synthpop_dir) && file.access(self$design$sim_prm$synthpop_dir, mode = 2) == -1L) {
          stop("You don't have write access to the synthpop folder. Please change the permissions or the path.  If you are using Linux you can use i.e. IMPACTncd$allow_universal_synthpop_folder_access() to allow write access to the synthpop folder.")
        }
        private$create_new_folder(self$design$sim_prm$output_dir, self$design$sim_prm$logs)
        private$create_new_folder(private$output_dir("summaries/"), self$design$sim_prm$logs)
        private$create_new_folder(private$output_dir("tables/"), self$design$sim_prm$logs)
        private$create_new_folder(private$output_dir("plots/"), self$design$sim_prm$logs)
        private$create_new_folder(private$output_dir("lifecourse/"), self$design$sim_prm$logs)
        if (self$design$sim_prm$export_PARF) {
           private$create_new_folder(private$output_dir("parf/"), self$design$sim_prm$logs)
        }
        if (self$design$sim_prm$export_xps) {
           private$create_new_folder(private$output_dir("xps/"), self$design$sim_prm$logs)
        }
        if (self$design$sim_prm$logs) {
           private$create_new_folder(private$output_dir("logs/"), self$design$sim_prm$logs)
        }
        # NOTE code below is duplicated in Synthpop class. This is intentional
        private$create_new_folder(self$design$sim_prm$synthpop_dir, self$design$sim_prm$logs)
      }
    ) # End of private methods
  ) # End of class
