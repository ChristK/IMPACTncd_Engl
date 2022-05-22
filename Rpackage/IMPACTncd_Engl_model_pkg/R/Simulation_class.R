## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngl is free software; you can redistribute it and/or modify it under
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
# field 'pop' that is a data.table
#' @export
`[.Simulation` <- function(x, ...) x$output[...]

#' R6 Class representing a simulation environment
#'
#' @description
#' A simulation environment.
#'
#' @details
#' To be completed...
#'
#' @export
Simulation <-
  R6::R6Class(
    classname = "Simulation",

    # public ------------------------------------------------------------------
    public = list(
      #' @field design A Design object with .
      design = NA,

      #' @field diseases A list of Disease objects.
      diseases = NA,

      #' @field scenarios A list of scenario objects.
      scenarios = NA,


      #' @description Create a new simulation object.
      #' @param sim_prm Either a path to a yaml file or a Design object.
      #' @return A new `Simulation` object.
      initialize = function(sim_prm) {
        if (is.character(sim_prm))
          self$design <- Design$new(sim_prm)
        else if (inherits(sim_prm, "Design"))
          self$design <- sim_prm$clone(deep = TRUE)
        else
          stop("sim_prm need to be a path to an appropriate yaml file or a Design object")

        data.table::setDTthreads(threads = self$design$sim_prm$clusternumber, restore_after_fork = NULL)
        fst::threads_fst(nr_of_threads = self$design$sim_prm$clusternumber, reset_after_fork = NULL)


        # Create folders if don't exist
        if (!dir.exists(self$design$sim_prm$output_dir)) {
          dir.create(self$design$sim_prm$output_dir, recursive = TRUE)
          if (self$design$sim_prm$logs)
            message(paste0("Folder ", self$design$sim_prm$output_dir,
                           " was created"))
        }

        pth <- private$output_dir("lifecourse/")
        if (!dir.exists(pth)) {
          dir.create(pth)
          if (self$design$sim_prm$logs)
            message(paste0("Folder ", pth, " was created"))
        }

        if (self$design$sim_prm$export_PARF) {
          pth <- private$output_dir("parf/")
          if (!dir.exists(pth)) {
            dir.create(pth)
            if (self$design$sim_prm$logs)
              message(paste0("Folder ", pth, " was created"))
          }
        }

        if (self$design$sim_prm$export_xps) {
          pth <- private$output_dir("xps/")
          if (!dir.exists(pth)) {
            dir.create(pth)
            if (self$design$sim_prm$logs)
              message(paste0("Folder ", pth, " was created"))
          }
        }

        if (self$design$sim_prm$logs) {
          pth <- private$output_dir("logs/")
          if (!dir.exists(pth)) {
            dir.create(pth)
            message(paste0("Folder ", pth, " was created"))
          }
        }

        # NOTE code below is duplicated in Synthpop class. This is intentional
        if (!dir.exists(self$design$sim_prm$synthpop_dir)) {
          dir.create(self$design$sim_prm$synthpop_dir, recursive = TRUE)
          if (self$design$sim_prm$logs)
            message(paste0("Folder ", self$design$sim_prm$synthpop_dir,
                         " was created"))
        }

        # RR ----
        # Create a named list of Exposure objects for the files in ./inputs/RR
        fl <- list.files(path = "./inputs/RR", pattern = ".csvy$", full.names = TRUE)
        # RR <- future_lapply(fl, Exposure$new, future.seed = 950480304L)
        RR <- lapply(fl, Exposure$new, design = self$design)
        names(RR) <- sapply(RR, function(x) x$get_name())
        # invisible(future_lapply(RR, function(x) {
        #   x$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)
        # },
        # future.seed = 627524136L))
        invisible(lapply(RR, function(x) {
          x$gen_stochastic_effect(self$design, overwrite = FALSE, smooth = FALSE)
        }))
        # NOTE smooth cannot be exported to Design for now, because the first time
        # this parameter changes we need logic to overwrite unsmoothed files
        rm(fl)

        # Generate diseases ----
        self$diseases <- lapply(self$design$sim_prm$diseases, function(x) {
          x[["design_"]] <- self$design
          x[["RR"]] <- RR
          do.call(Disease$new, x)
        })
        names(self$diseases) <- sapply(self$design$sim_prm$diseases, `[[`, "name")

        # Generate the graph with the causality structure
        ds <- unlist(strsplit(names(RR), "~"))
        ds[grep("^smok_", ds)] <- "smoking"
        ds <- gsub("_prvl$", "", ds)

        ds1 <- ds[as.logical(seq_along(ds) %% 2)]
        ds2 <- ds[!as.logical(seq_along(ds) %% 2)]
        ds <- unique(data.table(ds1, ds2))

        private$causality_structure <- make_graph(unlist(transpose(ds)),
                                                  directed = TRUE)

        invisible(self)
      },

      #' @description
      #' Runs a simulation
      #' @param mc_ An integer vector with the Monte Carlo iterations of
      #'   synthetic population to simulate.
      #' @param multicore If TRUE run the simulation in parallel.
      #' @return The invisible self for chaining.
      run = function(mc_, multicore = TRUE) {


        # Check if results for these mc_ exist from previous simulation

        mc_aggr <-
          unique(as.integer(ceiling(
            mc_ / self$design$sim_prm$n_synthpop_aggregation
          )))

        if (any(file.exists(
          file.path(
            self$design$sim_prm$output_dir,
            "lifecourse",
            paste0(mc_aggr, "_lifecourse.csv")
          )
        ))) {
          stop("Results from a previous simulation exists in the output folder.
               Please remove them before run a new one.")
        }

        # self$del_logs()

        # Generate PARF files if they don't exist. Note that generation is multicore
        lapply(self$diseases, function(x) {
          x$gen_parf_files(self$design, self$diseases)
          })

        if (multicore) {

          if (Sys.info()["sysname"] == "Windows") {
            cl <-
              makeCluster(self$design$sim_prm$clusternumber) # used for clustering. Windows compatible
            registerDoParallel(cl)
          } else {
            registerDoParallel(self$design$sim_prm$clusternumber) # used for forking. Only Linux/OSX compatible
          }

          if (self$design$sim_prm$logs)
            private$time_mark("Start of parallelisation")

          xps_dt <- foreach(
            mc_iter = mc_,
            .inorder = FALSE,
            .verbose = self$design$sim_prm$logs,
            .packages = c(
              "R6",
              "gamlss.dist",
              "dqrng",
              "CKutils",
              "IMPACTncdEngl",
              "fst",
              "data.table"
            ),
            .export = NULL,
            .noexport = NULL # c("time_mark")
          ) %dopar% {

            private$run_sim(mc_ = mc_iter)

          }

          if (exists("cl")) stopCluster(cl)

          if (self$design$sim_prm$logs) private$time_mark("End of parallelisation")


        } else {
          if (self$design$sim_prm$logs)
            private$time_mark("Start of single-core run")

          lapply(mc_, private$run_sim)

          if (self$design$sim_prm$logs)
            private$time_mark("End of single-core run")

        }

        while (sink.number() > 0L) sink()


        invisible(self)
        },


      #' @description Returns the causality matrix and optionally plots the
      #' causality structure.
      #' @param processed If `TRUE` generates the causality matrix from the graph.
      #' @param print_plot If `TRUE` prints the causal structure graph.
      #' @return The processed causality matrix if `processed = TRUE` or the graph
      #'   otherwise.
      get_causal_structure = function(processed = TRUE, print_plot = FALSE) {
        if (print_plot) {
          print(
            plot(
              private$causality_structure,
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
          g <- as.matrix(as_adjacency_matrix(private$causality_structure))
          n <- sapply(self$diseases, `[[`, "name")
          g <- g[!rownames(g) %in% n, colnames(g) %in% n]
        } else {
          g <- private$causality_structure
        }
        return(g)
      },

      #' @description
      #' Updates the Design object that is stored in the Simulation object.
      #' @param new_design A design object with the simulation parameters.
      #' @return The invisible self for chaining.
      update_design = function(new_design) {
        if (!inherits(new_design, "Design"))
          stop("Argument new_design needs to be a Design object.")

        self$design <- new_design

        invisible(self)
      },


      #' @description
      #' Delete all output files.
      #' @return The invisible self for chaining.
      del_outputs = function() {

        fl <- list.files(self$design$sim_prm$output_dir, full.names = TRUE,
                         recursive = TRUE)

        file.remove(fl)

        if (length(fl) > 0 && self$design$sim_prm$logs)
          message("Output files deleted.")

        invisible(self)
      },

      #' @description
      #' Delete log files.
      #' @return The invisible self for chaining.
      del_logs = function() {

        fl <- list.files(private$output_dir("logs/"), full.names = TRUE)

        file.remove(fl)

        if (length(fl) > 0 && self$design$sim_prm$logs)
          message("Log files deleted.")

        invisible(self)
      },


      #' @description
      #' Prints the simulation object metadata.
      #' @return The invisible `SynthPop` object.
      print = function() {
        print(c(
          "TODO..."
        ))
        invisible(self)
      }
    ),



    # private -----------------------------------------------------------------
    private = list(
      synthpop_dir = NA,
      causality_structure = NA,

      # Runs the simulation in one core
      run_sim = function(mc_) {

        if (self$design$sim_prm$logs) {
          private$time_mark(paste0("Start mc iteration ", mc_))
          sink(
            file = private$output_dir(paste0("logs/log", mc_, ".txt")),
            append = TRUE,
            type = "output",
            split = FALSE
          )
        }

        sp <- SynthPop$new(mc_, self$design)
        # ds <- copy(self$diseases) # Necessary for parallelisation
        lapply(self$diseases, function(x) {
          print(x$name)
          x$gen_parf(sp, self$design, self$diseases)$
            set_init_prvl(sp, self$design)$
            set_rr(sp, self$design)$
            set_incd_prb(sp, self$design)$
            set_dgns_prb(sp, self$design)$
            set_mrtl_prb(sp, self$design)
        })

        l <- private$mk_scenario_init(sp, "") # TODO update with scenarios
        simcpp(sp$pop, l, sp$mc)

        sp$update_pop_weights()

        if (self$design$sim_prm$export_xps) {
          if (self$design$sim_prm$logs) message("Exporting exposures...")
          private$export_xps(sp)
        }

        nam <- c(self$design$sim_prm$cols_for_output,
                 grep("_prvl$|_mrtl$", names(sp$pop), value = TRUE))

        sp$pop[, mc := sp$mc_aggr]

        if (self$design$sim_prm$logs) message("Exporting lifecourse...")
        fwrite_safe(sp$pop[all_cause_mrtl >= 0L, ..nam],
                    private$output_dir(paste0("/lifecourse/", sp$mc_aggr, "_lifecourse.csv")))

        if (self$design$sim_prm$logs) {

          private$time_mark(paste0("End mc iteration ", mc_))

          sink()

          }

        NULL
      },

      # creates the list that is used in c++ side
      # sp is needed for sp$mc_aggr in to_cpp()
      mk_scenario_init = function(sp, scenario_name) {
        # scenario_suffix_for_pop <- paste0("_", scenario_name) # TODO get suffix from design
        scenario_suffix_for_pop <- scenario_name
        list(
          "exposures"          = self$design$sim_prm$exposures,
          "scenarios"          = self$design$sim_prm$scenarios, # to be generated programmatically
          "scenario"           = self$scenarios, # TODO update when implement scenarios
          "kismet"             = self$design$sim_prm$kismet, # If TRUE random numbers are the same for each scenario.
          "init_year"          = self$design$sim_prm$init_year,
          "pids"               = "pid",
          "years"              = "year",
          "ages"               = "age",
          "ageL"               = self$design$sim_prm$ageL,
          "all_cause_mrtl"     = paste0("all_cause_mrtl", scenario_suffix_for_pop),
          "strata_for_outputs" = c("pid", "year", "age", "sex", "dimd"),
          "diseases"           = lapply(self$diseases, function(x) x$to_cpp(sp, self$design))
        )
      },

      # Function to export xps
      export_xps = function(sp) {
        # NOTE no need to check validity of inputs here as it is only used internally

        to_agegrp(sp$pop, grp_width = 20L, max_age = self$design$sim_prm$ageH,
                  min_age = self$design$sim_prm$ageL, age_colname = "age",
                  agegrp_colname = "agegrp20", to_factor = TRUE)

        sp$pop[, smok_never_curr_xps := fifelse(smok_status_curr_xps == "1", 1L, 0L)]
        sp$pop[, smok_active_curr_xps := fifelse(smok_status_curr_xps == "4", 1L, 0L)]

        xps <- grep("_curr_xps$", names(sp$pop), value = TRUE)
        xps <- xps[-which(xps %in% c("smok_status_curr_xps", "met_curr_xps",
                                     "bpmed_curr_xps", "t2dm_prvl_curr_xps",
                                     "af_prvl_curr_xps"))]
        sp$pop[smok_status_curr_xps == "1", `:=` (
          smok_packyrs_curr_xps = NA,
          smok_quit_yrs_curr_xps = NA,
          smok_dur_curr_xps = NA,
          smok_cig_curr_xps = NA
        )]
        out_xps <- groupingsets(
          sp$pop[all_cause_mrtl >= 0L &
                   year >= self$design$sim_prm$init_year &
                   age >= self$design$sim_prm$ageL, ],
          j = lapply(.SD, weighted.mean, wt, na.rm = TRUE),
          by = c("year", "sex", "agegrp20", "qimd", "ethnicity", "sha"),
          .SDcols = xps,
          sets = list(
            c("year", "sex", "agegrp20", "qimd"),
            c("year", "sex"),
            c("year", "agegrp20"),
            c("year", "qimd"),
            c("year", "ethnicity"),
            c("year", "sha")
          )
        )[, `:=` (year = year + 2000L, mc = sp$mc)] # TODO this could also be mc_aggr. Getting the uncertainty right here is tricky
        for (j in seq_len(ncol(out_xps)))
          set(out_xps, which(is.na(out_xps[[j]])), j, "All")
        sp$pop[, c(
          "agegrp20",
          "smok_never_curr_xps",
          "smok_active_curr_xps"
        ) := NULL]
        sp$pop[smok_status_curr_xps == "1", `:=` (
          smok_packyrs_curr_xps = 0,
          smok_quit_yrs_curr_xps = 0,
          smok_dur_curr_xps = 0,
          smok_cig_curr_xps = 0
        )]
        setkey(out_xps, year)

        fwrite_safe(out_xps, private$output_dir("xps/xps.csv"))

        NULL
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

      # Special deep copy for data.table. Use POP$clone(deep = TRUE) to
      # dispatch. Otherwise a reference is created
      deep_clone = function(name, value) {
        if ("data.table" %in% class(value)) {
          data.table::copy(value)
        } else if ("R6" %in% class(value)) {
          value$clone()
        } else {
          # For everything else, just return it. This results in a shallow
          # copy of s3.
          value
        }
      }


    )
  )
