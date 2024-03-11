## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by
## Chris Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngl is free software; you can redistribute it and/or modify it
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
#' A design has a sim_prm list that holds the simulation parameters.
#' This R6 class represents a simulation design with associated parameters and methods.
#'
#' @details
#' To be completed...
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
      #' @param sim_prm Either a path to a yaml file or a list with
      #'   appropriate format.
      #' @return A new `Design` object.
      #' @examples
      #' design <- Design$new("./validation/design_for_trends_validation.yaml")
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

        # Validation
        stopifnot(
          c(
            "simulation_files_overwrite",
            "sTag"                   ,
            "bOverwriteFilesOnDeploy",
            "RootDirPath"  ,
            "sToken"                ,
            "locality"              ,
            "iteration_n"           ,
            "iteration_n_final"     ,
            "clusternumber"         ,
            "logs"                  ,
            "scenarios"             ,
            "cols_for_output"       ,
            "strata_for_output"     ,
            "exposures"             ,
            "n"                     ,
            "init_year_long"        ,
            "sim_horizon_max"       ,
            "ageL"                  ,
            "ageH"                  ,
            "diseases"              ,
            "maxlag"                ,
            "smoking_relapse_limit" ,
            "stochastic"            ,
            "kismet"                ,
            "jumpiness"             ,
            "statin_adherence"      ,
            "bpmed_adherence"       ,
            "decision_aid"          ,
            "export_xps"            ,
            "simsmok_calibration"   ,
            "output_dir"            ,
            "synthpop_dir"          ,
            "validation"            ,
            "max_prvl_for_outputs"  ,
            "iteration_n_max"       ,
            "n_synthpop_aggregation"
          ) %in% names(sim_prm),

          sapply(sim_prm, function(x)
            if (is.numeric(x))
              x >= 0
            else
              TRUE)
        )


        sim_prm$sim_horizon_max <- sim_prm$sim_horizon_max - sim_prm$init_year_long
        sim_prm$init_year <- sim_prm$init_year_long - 2000L
        # place holders to be updated from self$update_fromGUI(parameters)

        sim_prm$national_qimd       <- TRUE
        sim_prm$init_year_fromGUI   <- sim_prm$init_year
        sim_prm$sim_horizon_fromGUI <- sim_prm$sim_horizon_max

        # Force LAD population proj if smaller localities
        if (!"England" %in% sim_prm$locality) sim_prm$calibrate_to_pop_projections_by_LAD <- TRUE

        # Create synthpop_dir_ if it doesn't exists
        sim_prm$output_dir <-
          normalizePath(sim_prm$output_dir, mustWork = FALSE)


        # Reorder the diseases so dependencies are always calculated first
        # (topological ordering). This is crucial for init_prevalence
        # first name the list and
        sim_prm$diseases <-
          setNames(sim_prm$diseases, sapply(sim_prm$diseases, function(x)
            x$name))

        out <- vector() # will hold graph structure
        ds <- names(sim_prm$diseases)
        for (i in seq_along(ds)) {
          ds_ <- ds[i]
          dep <-
            sim_prm[["diseases"]][[i]][["meta"]][["incidence"]][["influenced_by_disease_name"]]
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
        Cycles = NULL
        for(v1 in V(g)) {
          for(v2 in neighbors(g, v1, mode="out")) {
            Cycles = c(Cycles,
                       lapply(all_simple_paths(g, v2,v1, mode = "out"), function(p) c(v1,p)))
          }
        }
        # remove duplicates
        Cycles <- Cycles[sapply(Cycles, min) == sapply(Cycles, `[`, 1)]
        # find cycles of length i.e. 3 (i.e. chd -> t2dm -> chd)
        Cycles[which(sapply(Cycles, length) >= 3)]

        if (sim_prm$logs && length(Cycles) > 0) message("Cycles found: ", Cycles)

        o <- topo_sort(g)

        # then reorder based on the topological ordering
        sim_prm$diseases <- sim_prm$diseases[order(match(names(sim_prm$diseases), names(o)))]


        self$sim_prm = sim_prm

        invisible(self)
      },

      #' @description Create a new design object.
      #' @param path Path including file name and extension to save a yaml
      #'   file with the simulation parameters.
      #' @return The `Design` object.
      save_to_disk = function(path) {
        write_yaml(self$sim_prm, base::normalizePath(path, mustWork = FALSE))

        invisible(self)
      },

      #' @description Updates the design object from GUI.
      #' @param GUI_prm A GUI parameter object.
      #' @return The `Design` object.
      update_fromGUI = function(GUI_prm) {
        self$sim_prm$national_qimd       <- GUI_prm$national_qimd_checkbox
        # T = use national qimd, F = use local qimd
        self$sim_prm$init_year_fromGUI   <-
          fromGUI_timeframe(GUI_prm)["init year"] - 2000L
        self$sim_prm$sim_horizon_fromGUI <-
          fromGUI_timeframe(GUI_prm)["horizon"]
        self$sim_prm$locality <- GUI_prm$locality_select
        if (!GUI_prm$national_qimd_checkbox) {
          self$sim_prm$cols_for_output <-
            c(setdiff(self$sim_prm$cols_for_output, "lqimd"), "nqimd")
        }
        self$sim_prm$iteration_n            <- GUI_prm$iteration_n_gui
        self$sim_prm$iteration_n_final      <- GUI_prm$iteration_n_final_gui
        self$sim_prm$n                      <- GUI_prm$n_gui
        self$sim_prm$n_synthpop_aggregation <- GUI_prm$n_synthpop_aggregation_gui
        self$sim_prm$n_primers              <- GUI_prm$n_primers_gui
        self$sim_prm$cancer_cure            <- GUI_prm$cancer_cure_gui
        self$sim_prm$jumpiness              <- GUI_prm$jumpiness_gui
        self$sim_prm$statin_adherence       <- GUI_prm$statin_adherence_gui
        self$sim_prm$bpmed_adherence        <- GUI_prm$bpmed_adherence_gui
        self$sim_prm$decision_aid           <- GUI_prm$decision_aid_gui
        self$sim_prm$logs                   <- GUI_prm$logs_gui

        invisible(self)
      },

      #' @description
      #' Print the simulation parameters.
      #' @return The `Design` object.
      print = function() {
        print(self$sim_prm)

        invisible(self)
      }
    ),

    # private ------------------------------------------------------------------
     private = list(
      mc_aggr = NA
    )
  )
