## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
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


#' R6 Class representing a simulation design
#'
#' @description
#' A design has a sim_prm list that holds the simulation parameters.
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
      #' @field lags_mc The disease lag times.
      lags_mc = NA,
      #' @field max_lag_mc The longest disease lag time.
      max_lag_mc = NA,

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
            "iteration_n"           ,
            "iteration_n_final"     ,
            "clusternumber"         ,
            "n_cpus"                ,
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
            "cancer_cure"           ,
            "validation"            ,
            "max_prvl_for_outputs"  ,
            "iteration_n_max"       ,
            # "n_primers"             ,
            "n_synthpop_aggregation"
          ) %in% names(sim_prm),

          any(sim_prm$clusternumber == 1L, sim_prm$n_cpus == 1L),
          sapply(sim_prm, function(x)
            if (is.numeric(x))
              x >= 0
            else
              TRUE)
        )

        # tt <- sim_prm[grepl("_lag$", names(sim_prm))]
        # l <- lapply(names(tt), function(x) {
        #   disease_enum <-
        #     (cumsum(unlist(tt) == tt[[x]]) * duplicated(tt))[[x]]
        #   if (disease_enum == 0L)
        #     disease_enum <- 1L
        #   names(disease_enum) <- NULL
        #   return(disease_enum)
        # })
        # names(l) <- paste0(names(tt), "_enum")

        sim_prm$sim_horizon_max <- sim_prm$sim_horizon_max - sim_prm$init_year_long
        sim_prm$init_year <- sim_prm$init_year_long - 2000L
        # place holders to be updated from self$update_fromGUI(parameters)

        sim_prm$national_qimd       <- TRUE
        sim_prm$init_year_fromGUI   <- sim_prm$init_year
        sim_prm$sim_horizon_fromGUI <- sim_prm$sim_horizon_max
        sim_prm$locality            <- "England"

        # Create synthpop_dir_ if it doesn't exists
        sim_prm$output_dir <-
          base::normalizePath(sim_prm$output_dir, mustWork = FALSE)
        if (!dir.exists(sim_prm$output_dir)) {
          dir.create(sim_prm$output_dir, recursive = TRUE)
          message(paste0("Directory ", sim_prm$output_dir, " was created"))
        }


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
        self$sim_prm$n_cpus                 <- GUI_prm$n_cpus_gui
        self$sim_prm$n                      <- GUI_prm$n_gui
        self$sim_prm$n_synthpop_aggregation <- GUI_prm$n_synthpop_aggregation_gui
        self$sim_prm$n_primers              <- GUI_prm$n_primers_gui
        # self$sim_prm$cvd_lag                <- GUI_prm$cvd_lag_gui
        # self$sim_prm$copd_lag               <- GUI_prm$copd_lag_gui
        # self$sim_prm$cancer_lag             <- GUI_prm$cancer_lag_gui
        # self$sim_prm$nonmodelled_lag        <- GUI_prm$nonmodelled_lag_gui
        self$sim_prm$cancer_cure            <- GUI_prm$cancer_cure_gui
        self$sim_prm$jumpiness              <- GUI_prm$jumpiness_gui
        self$sim_prm$statin_adherence       <- GUI_prm$statin_adherence_gui
        self$sim_prm$bpmed_adherence        <- GUI_prm$bpmed_adherence_gui
        self$sim_prm$decision_aid           <- GUI_prm$decision_aid_gui
        self$sim_prm$logs                   <- GUI_prm$logs_gui

        # Gen new lags enum after deleting existing ones
        self$sim_prm[grepl("_enum$", names(self$sim_prm))] <- NULL

        tt <- self$sim_prm[grepl("_lag$", names(self$sim_prm))]
        l <- lapply(names(tt), function(x) {
          disease_enum <-
            (cumsum(unlist(tt) == tt[[x]]) * duplicated(tt))[[x]]
          if (disease_enum == 0L)
            disease_enum <- 1L
          names(disease_enum) <- NULL
          return(disease_enum)
        })
        names(l) <- paste0(names(tt), "_enum")
        self$sim_prm <- c(self$sim_prm, l)

        invisible(self)
      },

      #' @description Generates the lag per disease for a given Monte
      #'   Carlo iteration (`mc`).
      #' @param mc_ An integer for the Monte Carlo iteration.
      #' @return The `Design` object invisibly. Updates the fields
      #'   `lags_mc` and `max_lag_mc`
      get_lags = function(mc_) {
        if ((!is.na(private$mc_aggr) && mc_ != private$mc_aggr) ||
            is.na(self$lags_mc) || is.na(self$max_lag_mc)) {
          # for n_synthpop_aggregation > 1 we need RR and disease epi to change
          # every n_synthpop_aggregation, not in every mc_
          # n_synthpop_aggregation is handled at runtime. The line below is not
          # necessary
          # mc <-
          #   ceiling(mc_ / self$sim_prm$n_synthpop_aggregation)
          mc <- mc_
          # argument checks in the private function get_lag_mc_hlp
          nam <-
            grep("_lag$", names(self$sim_prm), value = TRUE)

          if (self$sim_prm$stochastic) {
            out <- vector("list", length(nam))
            names(out) <- nam
            self$lags_mc <-
              invisible(lapply(nam, function(x) {
                private$get_lag_mc_hlp(mc, self$sim_prm[[paste0(x, "_enum")]], self$sim_prm[[x]])
              }))
          } else {
            self$lags_mc <- self$sim_prm[nam]
          }
          names(self$lags_mc) <- nam
          self$lags_mc$plco_lag <-
            6L # for lung ca plco formula
          self$lags_mc$statin_t2dm_lag <-
            4L # for lung ca plco formula
          self$max_lag_mc <- max(unlist(self$lags_mc))
          private$mc_aggr <- mc_
        }
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
      mc_aggr = NA,

      get_lag_mc_hlp =
        function(mc,
                 disease_enum, # 1:10
                 lag) {
          # 2:9
          stopifnot(between(disease_enum, 1, 10),
                    between(lag, 2, 9),
                    between(mc, 1, 1e3))
          colnam <- paste0("lag_", lag)
          filenam <-
            "./disease_epidemiology/disease_lags_l.fst"
          filenam_indx <-
            "./disease_epidemiology/disease_lags_indx.fst"
          indx <-
            read_fst(
              filenam_indx,
              from = mc,
              to = mc,
              as.data.table = TRUE
            )
          out  <-
            read_fst(
              filenam,
              colnam,
              from = indx$from,
              to = indx$to,
              as.data.table = TRUE
            )[disease_enum, get(colnam)]
          out
        }
    )
  )
