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
#' @export
Simulation <-
  R6::R6Class(
    classname = "Simulation",
    lock_objects = TRUE, # allows primary prevention scenario to be updated
    lock_class = TRUE,
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
        
        private$create_empty_calibration_prms_file(replace = FALSE)
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

        private$primary_prevention_scn = function(synthpop) NULL # default for baseline scenario
        private$secondary_prevention_scn = function(synthpop) NULL # default for baseline scenario

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

      # calibrate_incd_ftlt ----
      #' @description Calibrate incidence and case fatality parameters
      #' @param mc A positive sequential integer vector with the Monte Carlo iterations
      #' @param replace If TRUE, replace the existing calibration file. If FALSE, resume from where left off.
      #' @return The invisible self for chaining.
      #' @details NOTE: This method requires England-specific disease burden files in ./inputs/disease_burden/
      #' including chd_incd.fst, stroke_incd.fst, chd_ftlt.fst, stroke_ftlt.fst, nonmodelled_ftlt.fst
      calibrate_incd_ftlt = function(mc, replace = FALSE) {
        # recombine the chunks of large files
        self$reconstruct_large_files()

        export_xps <- self$design$sim_prm$export_xps # save the original value to be restored later
        self$design$sim_prm$export_xps <- FALSE # turn off export_xps to speed up the calibration
        private$create_empty_calibration_prms_file(replace = replace)
        clbr <- fread(
          "./simulation/calibration_prms.csv",
          colClasses = list(
            numeric = c(
              "chd_incd_clbr_fctr",
              "stroke_incd_clbr_fctr",
              "chd_ftlt_clbr_fctr",
              "stroke_ftlt_clbr_fctr",
              "nonmodelled_ftlt_clbr_fctr"
            )
          )
        )

        memedian <- function(x) {
          out <- median(x)
          if (out == 0L) {
            # For rare events, consider alternative estimators
            nonzero_vals <- x[x > 0]
            if (length(nonzero_vals) > 0) {
              prob_occurrence <- length(nonzero_vals) / length(x)
              median_nonzero <- median(nonzero_vals)
              out <- prob_occurrence * median_nonzero
            } else {
              out <- 0
            }
          }
          out
        }

        if (replace) {
          age_start <- self$design$sim_prm$ageL
        } else {
          # if replace == FALSE
          # if all ages exist skip calibration
          if (
            dim(clbr[
              chd_incd_clbr_fctr == 1 |
                stroke_incd_clbr_fctr == 1 |
                chd_ftlt_clbr_fctr == 1 |
                stroke_ftlt_clbr_fctr == 1 |
                nonmodelled_ftlt_clbr_fctr == 1
            ])[1] ==
              0
          ) {
            if (self$design$sim_prm$logs) {
              message("All ages have been calibrated. Skipping calibration.")
            }
            return(invisible(self))
          }
          age_start <- clbr[
            chd_incd_clbr_fctr == 1 |
              stroke_incd_clbr_fctr == 1 |
              chd_ftlt_clbr_fctr == 1 |
              stroke_ftlt_clbr_fctr == 1 |
              nonmodelled_ftlt_clbr_fctr == 1,
            min(age)
          ]
          if (self$design$sim_prm$logs) {
            message("Starting calibration from age ", age_start, ".")
          }
        }

        # Run the simulation from min to max age
        for (age_ in age_start:self$design$sim_prm$ageH) {
          # Run the simulation and export summaries
          self$del_logs()$del_outputs()$run(
            mc,
            multicore = TRUE,
            "sc0"
          )$export_summaries(
            multicore = TRUE,
            type = c("incd", "prvl", "dis_mrtl"),
            single_year_of_age = TRUE
          )

          # Incidence calibration
          # load the uncalibrated results
          unclbr <- open_dataset(file.path(
            self$design$sim_prm$output_dir,
            "summaries",
            "incd_scaled_up"
          )) %>%
            filter(age == age_) %>%
            select(
              "year",
              "age",
              "sex",
              "mc",
              "popsize",
              "chd_incd",
              "stroke_incd"
            ) %>%
            collect()
          setDT(unclbr)

          unclbr <- unclbr[,
            .(
              age,
              sex,
              year,
              mc,
              chd_incd = chd_incd / popsize,
              stroke_incd = stroke_incd / popsize
            )
          ][,
            .(
              chd_incd = memedian(chd_incd),
              stroke_incd = memedian(stroke_incd)
            ),
            keyby = .(age, sex, year)
          ]

          # for CHD
          # fit a log-log linear model to the uncalibrated results and store the coefficients
          tt <- unclbr[
            chd_incd > 0,
            as.list(coef(lm(
              log(chd_incd) ~ log(year)
            ))),
            by = sex
          ]

          unclbr[
            tt,
            on = "sex",
            c("intercept_unclbr", "trend_unclbr") := .(
              `(Intercept)`,
              `log(year)`
            )
          ]
          rm(tt)

          # load benchmark
          benchmark <- read_fst(
            file.path("./inputs/disease_burden", "chd_incd.fst"),
            columns = c("age", "sex", "year", "mu"),
            as.data.table = TRUE
          )[age == age_, ]
          # fit a log-log linear model to the benchmark incidence and store the coefficients
          benchmark[
            year >= self$design$sim_prm$init_year_long,
            c("intercept_bnchmrk", "trend_bnchmrk") := as.list(coef(lm(
              log(mu) ~ log(year)
            ))),
            by = sex
          ]

          # calculate the calibration factors
          unclbr[
            benchmark[year == max(year)],
            chd_incd_clbr_fctr := exp(
              intercept_bnchmrk + trend_bnchmrk * log(year)
            ) /
              exp(intercept_unclbr + trend_unclbr * log(year)),
            on = c("age", "sex")
          ]
          unclbr[, c("intercept_unclbr", "trend_unclbr") := NULL]

          # Repeat for stroke
          tt <- unclbr[
            stroke_incd > 0,
            as.list(coef(lm(
              log(stroke_incd) ~ log(year)
            ))),
            by = sex
          ]

          unclbr[
            tt,
            on = "sex",
            c("intercept_unclbr", "trend_unclbr") := .(
              `(Intercept)`,
              `log(year)`
            )
          ]
          rm(tt)

          benchmark <- read_fst(
            file.path("./inputs/disease_burden", "stroke_incd.fst"),
            columns = c("age", "sex", "year", "mu"),
            as.data.table = TRUE
          )[age == age_, ]
          benchmark[
            year >= self$design$sim_prm$init_year_long,
            c(
              "intercept_bnchmrk",
              "trend_bnchmrk"
            ) := as.list(coef(lm(
              log(mu) ~ log(year)
            ))),
            by = sex
          ]

          unclbr[
            benchmark[year == max(year)],
            stroke_incd_clbr_fctr := exp(
              intercept_bnchmrk + trend_bnchmrk * log(year)
            ) /
              exp(intercept_unclbr + trend_unclbr * log(year)),
            on = c("age", "sex")
          ]

          unclbr[, `:=`(
            chd_prvl_correction = chd_incd * (chd_incd_clbr_fctr - 1),
            stroke_prvl_correction = stroke_incd * (stroke_incd_clbr_fctr - 1),
            chd_incd = NULL,
            stroke_incd = NULL,
            intercept_unclbr = NULL,
            trend_unclbr = NULL
          )]
          clbr[
            unclbr,
            on = c("year", "age", "sex"),
            `:=`(
              chd_incd_clbr_fctr = i.chd_incd_clbr_fctr,
              stroke_incd_clbr_fctr = i.stroke_incd_clbr_fctr
            )
          ]

          # Case fatality calibration
          prvl <- open_dataset(file.path(
            self$design$sim_prm$output_dir,
            "summaries",
            "prvl_scaled_up"
          )) %>%
            filter(age == age_) %>%
            select(
              "year",
              "age",
              "sex",
              "mc",
              "popsize",
              "chd_prvl",
              "stroke_prvl"
            ) %>%
            collect()
          setDT(prvl)

          prvl <- prvl[, .(
            chd_prvl = chd_prvl / popsize,
            stroke_prvl = stroke_prvl / popsize,
            popsize,
            age,
            sex,
            year,
            mc
          )][,
            .(
              chd_prvl = memedian(chd_prvl),
              stroke_prvl = memedian(stroke_prvl),
              popsize = memedian(popsize)
            ),
            keyby = .(age, sex, year)
          ]
          prvl[
            unclbr,
            on = c("year", "age", "sex"),
            `:=`(
              chd_prvl_correction = i.chd_prvl_correction,
              stroke_prvl_correction = i.stroke_prvl_correction
            )
          ]
          benchmark <- read_fst(
            file.path("./inputs/disease_burden", "chd_ftlt.fst"),
            columns = c("age", "sex", "year", "mu2"),
            as.data.table = TRUE
          )[age == age_, ]
          prvl[benchmark, on = c("age", "sex", "year"), chd_mrtl := mu2]
          benchmark <- read_fst(
            file.path("./inputs/disease_burden", "stroke_ftlt.fst"),
            columns = c("age", "sex", "year", "mu2"),
            as.data.table = TRUE
          )[age == age_, ]
          prvl[benchmark, on = c("age", "sex", "year"), stroke_mrtl := mu2]
          benchmark <- read_fst(
            file.path("./inputs/disease_burden", "nonmodelled_ftlt.fst"),
            columns = c("age", "sex", "year", "mu2"),
            as.data.table = TRUE
          )[age == age_, ]
          prvl[benchmark, on = c("age", "sex", "year"), nonmodelled_mrtl := mu2]

          prvl[, `:=`(
            chd_ftlt_clbr_fctr = 1 / (chd_prvl + chd_prvl_correction),
            stroke_ftlt_clbr_fctr = 1 / (stroke_prvl + stroke_prvl_correction),
            nonmodelled_ftlt_clbr_fctr = 1 / (1 - chd_mrtl - stroke_mrtl)
          )]

          # Fix the calibration factors for the ages that have been calibrated
          if (age_ > age_start) {
            mrtl <- open_dataset(file.path(
              self$design$sim_prm$output_dir,
              "summaries",
              "dis_mrtl_scaled_up"
            )) %>%
              filter(age == age_ - 1L) %>%
              select(
                "year",
                "age",
                "sex",
                "mc",
                "popsize",
                "chd_deaths",
                "stroke_deaths",
                "nonmodelled_deaths"
              ) %>%
              collect()
            setDT(mrtl)

            mrtl <- mrtl[, .(
              chd_mrtl = chd_deaths / popsize,
              stroke_mrtl = stroke_deaths / popsize,
              nonmodelled_mrtl = nonmodelled_deaths / popsize,
              popsize,
              age,
              sex,
              year,
              mc
            )][,
              .(
                chd_mrtl = memedian(chd_mrtl),
                stroke_mrtl = memedian(stroke_mrtl),
                nonmodelled_mrtl = memedian(nonmodelled_mrtl),
                popsize = memedian(popsize)
              ),
              keyby = .(age, sex, year)
            ]
            benchmark <- read_fst(
              file.path("./inputs/disease_burden", "chd_ftlt.fst"),
              columns = c("age", "sex", "year", "mu2"),
              as.data.table = TRUE
            )[age == age_ - 1L, ]
            mrtl[
              benchmark,
              on = c("age", "sex", "year"),
              chd_ftlt_clbr_fctr := mu2 / chd_mrtl
            ]
            benchmark <- read_fst(
              file.path("./inputs/disease_burden", "stroke_ftlt.fst"),
              columns = c("age", "sex", "year", "mu2"),
              as.data.table = TRUE
            )[age == age_ - 1L, ]
            mrtl[
              benchmark,
              on = c("age", "sex", "year"),
              stroke_ftlt_clbr_fctr := mu2 / stroke_mrtl
            ]
            benchmark <- read_fst(
              file.path("./inputs/disease_burden", "nonmodelled_ftlt.fst"),
              columns = c("age", "sex", "year", "mu2"),
              as.data.table = TRUE
            )[age == age_ - 1L, ]
            mrtl[
              benchmark,
              on = c("age", "sex", "year"),
              nonmodelled_ftlt_clbr_fctr := mu2 / nonmodelled_mrtl
            ]
            mrtl[chd_ftlt_clbr_fctr == Inf, chd_ftlt_clbr_fctr := 1]
            mrtl[stroke_ftlt_clbr_fctr == Inf, stroke_ftlt_clbr_fctr := 1]
            mrtl[
              nonmodelled_ftlt_clbr_fctr == Inf,
              nonmodelled_ftlt_clbr_fctr := 1
            ]

            clbr[
              mrtl,
              on = c("year", "age", "sex"),
              `:=`(
                chd_ftlt_clbr_fctr = i.chd_ftlt_clbr_fctr * chd_ftlt_clbr_fctr,
                stroke_ftlt_clbr_fctr = i.stroke_ftlt_clbr_fctr *
                  stroke_ftlt_clbr_fctr,
                nonmodelled_ftlt_clbr_fctr = i.nonmodelled_ftlt_clbr_fctr *
                  nonmodelled_ftlt_clbr_fctr
              )
            ]
          }

          if (age_ == self$design$sim_prm$ageH) {
            # shortcut for age == 99 hopefully with tiny bias
            mrtl[, age := age + 1L]
            prvl[
              mrtl,
              on = c("year", "age", "sex"),
              `:=`(
                chd_ftlt_clbr_fctr = i.chd_ftlt_clbr_fctr * chd_ftlt_clbr_fctr,
                stroke_ftlt_clbr_fctr = i.stroke_ftlt_clbr_fctr *
                  stroke_ftlt_clbr_fctr,
                nonmodelled_ftlt_clbr_fctr = i.nonmodelled_ftlt_clbr_fctr *
                  nonmodelled_ftlt_clbr_fctr
              )
            ]
          }

          clbr[
            prvl,
            on = c("year", "age", "sex"),
            `:=`(
              chd_ftlt_clbr_fctr = i.chd_ftlt_clbr_fctr,
              stroke_ftlt_clbr_fctr = i.stroke_ftlt_clbr_fctr,
              nonmodelled_ftlt_clbr_fctr = i.nonmodelled_ftlt_clbr_fctr
            )
          ]

          fwrite(clbr, "./simulation/calibration_prms.csv")
        } # end loop over ages

        self$design$sim_prm$export_xps <- export_xps # restore the original value
        invisible(self)
      },

      # validate ----
      #' @description Validate the simulation results against observed data
      #' @return The invisible self for chaining.
      #' @details NOTE: This method requires England-specific population files.
      #' Update the file path below to point to your England population data.
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

        # check if sequential vector. Necessary if
        # design$sim_prm$n_synthpop_aggregation > 1
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
            self$design$sim_prm$n_synthpop_aggregation -
            self$design$sim_prm$n_synthpop_aggregation +
            1L):(max(mc) * self$design$sim_prm$n_synthpop_aggregation)

        # Create folders if don't exist (necessary for when output_dir in the
        # design.yaml is changed between scenarios i.e.)
        private$create_output_folder_structure()

        # TODO better logic as this is always true for the non baseline scenario
        if (
          any(file.exists(
            # TODO fix when lifecourse is not saved
            file.path(
              self$design$sim_prm$output_dir,
              "lifecourse",
              paste0(mc, "_lifecourse.csv.gz")
            )
          ))
        ) {
          # stop("Results from a previous simulation exists in the output
          #      folder. Please remove them before run a new one.")
          message(
            "Results from a previous simulation exists in the output folder. This is usually results from a previous scenario. Please remove them if this was unintentional."
          )
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

          lapply(mc_sp, private$run_sim, scenario_nam)

          if (self$design$sim_prm$logs) {
            private$time_mark("End of single-core run")
          }
        }

        if (self$design$sim_prm$avoid_appending_csv) {
          message(
            "Collecting the fragmented lifecourse files. This may take some time. Please be patient..."
          )
          private$collect_files(
            "lifecourse",
            "_lifecourse.csv$",
            to_mc_aggr = TRUE
          )

          if (self$design$sim_prm$export_xps) {
            private$collect_files("xps", "_xps20.csv$", to_mc_aggr = FALSE)
            private$collect_files("xps", "_xps_esp.csv$", to_mc_aggr = FALSE)
          }

          if (self$design$sim_prm$logs) {
            private$time_mark("End of collecting mc lifecourse files")
          }
        }

        while (sink.number() > 0L) {
          sink()
        }

        invisible(self)
      },

      # export_summaries ----

      #' @description Process the lifecourse files
      #' @param multicore If TRUE run the simulation in parallel.
      #' @param type The type of summary to extract.
      #' @param single_year_of_age Export summaries by single year of age. Useful for the calibration proccess.
      #' @return The invisible self for chaining.
      export_summaries = function(
        multicore = TRUE,
        type = c(
          "le",
          "hle",
          "dis_char",
          "prvl",
          "incd",
          "dis_mrtl",
          "mrtl",
          "allcause_mrtl_by_dis",
          "cms"
        ),
        single_year_of_age = FALSE
      ) {
        fl <- list.files(private$output_dir("lifecourse"), full.names = TRUE)

        # logic to avoid inappropriate dual processing of already processed mcs
        # TODO take into account scenarios
        if ("le" %in% type) {
          file_pth <- private$output_dir("summaries/le_scaled_up.csv.gz")
        } else if ("hle" %in% type) {
          file_pth <- private$output_dir(
            "summaries/hle_1st_cond_scaled_up.csv.gz"
          )
        } else if ("cms" %in% type) {
          file_pth <- private$output_dir("summaries/cms_count_scaled_up.csv.gz")
        } else if ("mrtl" %in% type) {
          file_pth <- private$output_dir("summaries/mrtl_scaled_up.csv.gz")
        } else if ("dis_mrtl" %in% type) {
          file_pth <- private$output_dir("summaries/dis_mrtl_scaled_up.csv.gz")
        } else if ("dis_char" %in% type) {
          file_pth <- private$output_dir(
            "summaries/dis_characteristics_scaled_up.csv.gz"
          )
        } else if ("incd" %in% type) {
          file_pth <- private$output_dir("summaries/incd_scaled_up.csv.gz")
        } else if ("prvl" %in% type) {
          file_pth <- private$output_dir("summaries/prvl_scaled_up.csv.gz")
        } else if ("allcause_mrtl_by_dis" %in% type) {
          file_pth <- private$output_dir(
            "summaries/all_cause_mrtl_by_dis_scaled_up.csv.gz"
          )
        }

        if (file.exists(file_pth)) {
          tt <- unique(fread(file_pth, select = "mc")$mc)
          for (i in seq_along(tt)) {
            fl <- grep(
              paste0("/", tt[[i]], "_lifecourse.csv.gz$"),
              fl,
              value = TRUE,
              invert = TRUE
            )
          }
        }
        # end of logic

        if (multicore) {
          if (self$design$sim_prm$logs) {
            private$time_mark("Start exporting summaries")
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
                self$design$sim_prm$clusternumber_export,
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

            parLapplyLB(
              cl = cl,
              X = seq_along(fl),
              fun = function(i) {
                lc <- fread(
                  fl[i],
                  stringsAsFactors = TRUE,
                  key = c("scenario", "pid", "year")
                )
                private$export_summaries_hlpr(
                  lc,
                  type = type,
                  single_year_of_age = single_year_of_age
                )
                NULL
              }
            )
          } else {
            registerDoParallel(self$design$sim_prm$clusternumber_export) # used for forking. Only Linux/OSX compatible
            xps_dt <- foreach(
              i = seq_along(fl),
              .inorder = TRUE,
              .options.multicore = list(preschedule = FALSE),
              .verbose = self$design$sim_prm$logs,
              .packages = c(
                "R6",
                "CKutils",
                "IMPACTncdEngland",
                "data.table"
              ),
              .export = NULL,
              .noexport = NULL # c("time_mark")
            ) %dopar%
              {
                lc <- fread(
                  fl[i],
                  stringsAsFactors = TRUE,
                  key = c("scenario", "pid", "year")
                )
                private$export_summaries_hlpr(
                  lc,
                  type = type,
                  single_year_of_age = single_year_of_age
                )
                NULL
              }
          }

          if (self$design$sim_prm$logs) {
            private$time_mark("End of exporting summuries")
          }
        } else {
          if (self$design$sim_prm$logs) {
            private$time_mark("Start of single-core run")
          }

          lapply(seq_along(fl), function(i) {
            lc <- fread(fl[i], stringsAsFactors = TRUE, key = c("pid", "year"))
            private$export_summaries_hlpr(
              lc,
              type = type,
              single_year_of_age = single_year_of_age
            )
            NULL
          })

          if (self$design$sim_prm$logs) {
            private$time_mark("End of single-core run")
          }
        }

        if (self$design$sim_prm$avoid_appending_csv) {
          # collect the summary fragmentrd file
          if ("le" %in% type) {
            private$collect_files(
              "summaries",
              "_le_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_le_esp.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_le60_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_le60_esp.csv$",
              to_mc_aggr = FALSE
            )
          }
          if ("hle" %in% type) {
            private$collect_files(
              "summaries",
              "_hle_1st_cond_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_hle_1st_cond_esp.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_hle_cmsmm1.5_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_hle_cmsmm1.5_esp.csv$",
              to_mc_aggr = FALSE
            )
          }
          if ("cms" %in% type) {
            private$collect_files(
              "summaries",
              "_cms_score_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_cms_score_esp.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_cms_score_by_age_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_cms_score_by_age_esp.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_cms_count_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_cms_count_esp.csv$",
              to_mc_aggr = FALSE
            )
          }
          if ("mrtl" %in% type) {
            private$collect_files(
              "summaries",
              "_mrtl_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_mrtl_esp.csv$",
              to_mc_aggr = FALSE
            )
          }
          if ("dis_mrtl" %in% type) {
            private$collect_files(
              "summaries",
              "_dis_mrtl_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_dis_mrtl_esp.csv$",
              to_mc_aggr = FALSE
            )
          }
          if ("dis_char" %in% type) {
            private$collect_files(
              "summaries",
              "_dis_characteristics_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_dis_characteristics_esp.csv$",
              to_mc_aggr = FALSE
            )
          }
          if ("incd" %in% type) {
            private$collect_files(
              "summaries",
              "_incd_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_incd_esp.csv$",
              to_mc_aggr = FALSE
            )
          }
          if ("prvl" %in% type) {
            private$collect_files(
              "summaries",
              "_prvl_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_prvl_esp.csv$",
              to_mc_aggr = FALSE
            )
          }
          if ("allcause_mrtl_by_dis" %in% type) {
            private$collect_files(
              "summaries",
              "_all_cause_mrtl_by_dis_scaled_up.csv$",
              to_mc_aggr = FALSE
            )
            private$collect_files(
              "summaries",
              "_all_cause_mrtl_by_dis_esp.csv$",
              to_mc_aggr = FALSE
            )
          }

          if (self$design$sim_prm$logs) {
            private$time_mark("End of collecting mc_aggr summary files")
          }
        }

        while (sink.number() > 0L) {
          sink()
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
      #' @return The processed causality matrix if `processed = TRUE` or the
      #'   graph otherwise.
      get_causal_structure = function(
        processed = TRUE,
        print_plot = FALSE,
        focus = FALSE
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
          graph <- make_ego_graph(
            private$causality_structure,
            order = 1,
            nodes = focus,
            mode = "in"
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
          n <- sapply(self$design$diseases, `[[`, "name")
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

      #' @description Delete all output files.
      #' @return The invisible self for chaining.
      del_outputs = function() {
        if (dir.exists(self$design$sim_prm$output_dir)) {
          # Check for safety that folders /lifecourse, /tables, /plots, and /summaries exist to avoid accidental deletes of other folders
          if (
            dir.exists(file.path(
              self$design$sim_prm$output_dir,
              "lifecourse"
            )) &&
              dir.exists(file.path(
                self$design$sim_prm$output_dir,
                "summaries"
              )) &&
              dir.exists(file.path(self$design$sim_prm$output_dir, "tables")) &&
              dir.exists(file.path(self$design$sim_prm$output_dir, "plots"))
          ) {
            fl <- list.files(
              self$design$sim_prm$output_dir,
              full.names = TRUE,
              recursive = TRUE
            )
            file.remove(fl)

            if (length(fl) > 0 && self$design$sim_prm$logs) {
              message("Output files deleted.")
            }
          } else {
            message(
              "Output folder doesn't contain the expected subfolders. Please check the output folder path."
            )
          }
        } else {
          # If output folder doesn't exist
          message("Output folder doesn't exist.")
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

        file.remove(fl)

        if (length(fl) > 0 && self$design$sim_prm$logs) {
          message("Parf files deleted.")
        }

        invisible(self)
      },

      # del_synthpops ----
      #' @description Delete all files in the synthpop folder.
      #' @return The invisible self for chaining.
      del_synthpops = function() {
        fl <- list.files(self$design$sim_prm$synthpop_dir, full.names = TRUE)

        file.remove(fl)

        if (length(fl) > 0 && self$design$sim_prm$logs) {
          message("Synthpop files deleted.")
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
        for (i in 1:length(fl)) {
          file <- gsub("^./", "", fl[i])
          if (file %in% excl) {
            next
          }
          write(file, file = "./.gitignore", append = TRUE)
        }

        # split the files into 50MB chunks
        for (i in 1:length(fl)) {
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
          file.remove(file)
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
          for (i in 1:length(fl)) {
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
          file.remove(fl)
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

      # export_tables ----
      #' @description
      #' Export summary tables for the simulation results.
      #'
      #' This method generates and exports summary tables for the main simulation outputs,
      #' including prevalence, incidence, mortality, disease characteristics, and exposures.
      #' It calls modular helper methods for each type of summary, ensuring output directories
      #' are created as needed and that all tables are written to the appropriate locations.
      #'
      #' @param baseline_year_for_change_outputs Integer. The baseline year to use for change outputs (default: 2019L).
      #' @param prbl Numeric vector. The quantiles to use for summary statistics (default: c(0.5, 0.025, 0.975, 0.1, 0.9)).
      #'
      #' @details
      #' This method is a high-level wrapper that orchestrates the export of all main summary tables.
      #' It delegates the actual export logic to the following private helper methods:
      #' - \code{private$export_main_tables}
      #' - \code{private$export_all_cause_mrtl_tables}
      #' - \code{private$export_disease_characteristics_tables}
      #' - \code{private$export_xps_tables}
      #'
      #' Each helper method is responsible for a specific set of outputs and ensures that
      #' the results are saved in the correct format and location.
      #'
      #' @return The invisible self for chaining.
      #'
      #' @examples
      #' IMPACTncd$export_tables()
      export_tables = function(
        baseline_year_for_change_outputs = 2019L,
        prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)
      ) {
        private$export_main_tables(
          prbl,
          baseline_year_for_change_outputs,
          private$output_dir()
        )
        private$export_all_cause_mrtl_tables(
          prbl,
          private$output_dir("summaries"),
          private$output_dir("tables")
        )
        private$export_disease_characteristics_tables(
          prbl,
          private$output_dir("summaries"),
          private$output_dir("tables")
        )
        private$export_xps_tables(
          prbl,
          private$output_dir(),
          private$output_dir("tables")
        )

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
      },

      # print ----

      #' @description Prints the simulation object metadata.
      #' @return The invisible `Simulation` object.
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
      death_codes = NA,
      # diseasenam_hlp = NA,
      esp_weights = data.table(),
      #Models a primary prevention policy scenario
      primary_prevention_scn = NULL,
      #Models a secondary prevention policy scenario
      secondary_prevention_scn = NULL,

      # create_empty_calibration_prms_file ----
      # Helper function to create an empty calibration parameters file
      create_empty_calibration_prms_file = function(replace = FALSE) {
        if (replace || !file.exists("./simulation/calibration_prms.csv")) {
          clbr <- CJ(
            year = self$design$sim_prm$init_year_long:self$design$sim_prm$sim_horizon_max,
            age = self$design$sim_prm$ageL:self$design$sim_prm$ageH,
            sex = c("men", "women")
          )
          clbr[, `:=`(
            chd_incd_clbr_fctr = 1,
            stroke_incd_clbr_fctr = 1,
            chd_ftlt_clbr_fctr = 1,
            stroke_ftlt_clbr_fctr = 1,
            nonmodelled_ftlt_clbr_fctr = 1
          )]
          fwrite(clbr, "./simulation/calibration_prms.csv")
        }
      },

      # run_sim ----
      # Runs the simulation in one core. mc is scalar
      run_sim = function(mc_, scenario_nam = "") {
        if (!nzchar(scenario_nam)) {
          scenario_nam <- "sc0"
        }

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

        private$primary_prevention_scn(sp) # apply primary pevention scenario

        lapply(self$design$diseases, function(x) {
          x$set_rr(sp, self$design)$set_incd_prb(sp, self$design)$set_dgns_prb(
            sp,
            self$design
          )$set_mrtl_prb(sp, self$design)
        })

        private$secondary_prevention_scn(sp) # apply secondary pevention scenario

        l <- private$mk_scenario_init(sp, scenario_nam)
        if (!identical(key(sp$pop), c("pid", "year"))) {
          stop("synthpop key is not as expected")
        }
        simcpp(sp$pop, l, sp$mc)
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
        absorb_dt(sp$pop, private$esp_weights)
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

        sp$pop[, scenario := scenario_nam]

        setkeyv(sp$pop, c("pid", "year"))

        # Write lifecourse
        if (self$design$sim_prm$logs) {
          message("Exporting lifecourse...")
        }

        if (self$design$sim_prm$avoid_appending_csv) {
          fnam <- private$output_dir(paste0(
            "lifecourse/",
            sp$mc_aggr,
            "_",
            sp$mc,
            "_lifecourse.csv"
          ))
        } else {
          fnam <- private$output_dir(paste0(
            "lifecourse/",
            sp$mc_aggr,
            "_lifecourse.csv.gz"
          ))
        }
        fwrite_safe(sp$pop, fnam)

        if (self$design$sim_prm$logs) {
          private$time_mark(paste0("End mc iteration ", mc_))
          sink()
        }

        NULL
      },

      # creates the list that is used in c++ side sp is needed for sp$mc_aggr in
      # to_cpp()

      # mk_scenario_init ----
      # functions - ⁠carry_forward_inc;
      # Create initialization parameters for a simulation scenario.
      # @param sp A SynthPop object.
      # @param scenario_name A character string specifying the name of the scenario.
      #
      # @return A list of initialization parameters for the scenario.
      #
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
      # ? ARU check details
      # Export Exposure Data
      # @param sp A SynthPop object.
      # @param scenario_nam A character string specifying the scenario name.
      # @return NULL
      #
      # @details This function exports exposure data from a SynthPop object to CSV files.
      #
      export_xps = function(sp, scenario_nam) {
        # NOTE no need to check validity of inputs here as it is only used
        # internally

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
        xps <- xps[
          -which(
            xps %in% c("smok_status_curr_xps", "met_curr_xps", "bpmed_curr_xps")
          )
        ]
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
        if (self$design$sim_prm$avoid_appending_csv) {
          fwrite_safe(
            out_xps20,
            private$output_dir(paste0("xps/", sp$mc, "_xps20.csv"))
          )
        } else {
          fwrite_safe(out_xps20, private$output_dir("xps/xps20.csv.gz"))
        }

        # TODO link strata in the outputs to the design.yaml
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
        if (self$design$sim_prm$avoid_appending_csv) {
          fwrite_safe(
            out_xps5,
            private$output_dir(paste0("xps/", sp$mc, "_xps_esp.csv"))
          )
        } else {
          fwrite_safe(out_xps5, private$output_dir("xps/xps_esp.csv.gz"))
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

        NULL
      },

      # time_mark ----
      # Function for timing log
      # @param text_id A character string representing the text identifier.
      # @return NULL
      #
      # @details
      # The function logs a timestamp, along with the provided text identifier, to a file.
      # The file is specified as "logs/times.txt" within the output directory. The logging
      # is done in append mode, meaning new entries are added to the end of the existing file.
      # The logged information includes the text identifier, timestamp, and a newline character.
      #
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

      # output_dir ----
      output_dir = function(x = "") {
        file.path(self$design$sim_prm$output_dir, x)
      },

      # export_summaries_hlpr ----
      # function to export summaries from lifecourse files.
      # lc is a lifecourse file
      # ? ARU check short forms
      # Export Summaries Helper Function
      # @param lc Data table containing lifecourse information.
      # @param type Character vector specifying the types of summaries to export. Possible values include "le" (life expectancy), "hle" (healthy life expectancy), "dis_char" (disease characteristics), "prvl" (prevalence), "incd" (incidence), "mrtl" (mortality), "dis_mrtl" (disease-specific mortality), "allcause_mrtl_by_dis" (all-cause mortality by disease), and "cms" (CMS - Chronic Multimorbidity Score).
      # @param design An object representing the design of the simulation. It should contain simulation parameters.
      #
      # @details
      # The function exports various summaries based on the specified types. For example, it exports life expectancy, healthy life expectancy, disease characteristics, prevalence, incidence, mortality, disease-specific mortality, all-cause mortality by disease, and Chronic Multimorbidity Score (CMS) information.
      #
      # @return
      # Returns an invisible updated simulation object.
      #
      export_summaries_hlpr = function(
        lc,
        type = c(
          "le",
          "hle",
          "dis_char",
          "prvl",
          "incd",
          "mrtl",
          "dis_mrtl",
          "allcause_mrtl_by_dis",
          "cms"
        ),
        single_year_of_age = FALSE
      ) {
        if (self$design$sim_prm$logs) {
          message("Exporting summaries...")
        }

        strata <- c("mc", self$design$sim_prm$strata_for_output)
        strata_noagegrp <- c(
          "mc",
          setdiff(self$design$sim_prm$strata_for_output, c("agegrp"))
        )
        strata_age <- c(strata_noagegrp, "age")

        if (single_year_of_age) {
          strata <- strata_age
        } # used for calibrate_incd_ftlt

        setkeyv(lc, c("scenario", "pid", "year")) # necessary for age_onset

        mcaggr <- ifelse(
          self$design$sim_prm$avoid_appending_csv,
          paste0(lc$mc[1], "_"),
          ""
        )
        ext <- ifelse(
          self$design$sim_prm$avoid_appending_csv,
          ".csv",
          ".csv.gz"
        )

        # Life expectancy ----
        # NOTE for scaled_up LE weights need to apply from the very beginning.
        # Also note that currently this ignores the deaths for people younger
        # than min_age so not a true LE at birth
        if ("le" %in% type) {
          # fwrite_safe(lc[all_cause_mrtl > 0, .("popsize" = (.N), LE = mean(age)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0("summaries/", "le_out.csv.gz"
          #             )))
          fwrite_safe(
            lc[
              all_cause_mrtl > 0,
              .("popsize" = sum(wt), LE = weighted.mean(age, wt)),
              keyby = strata_noagegrp
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "le_scaled_up",
              ext
            ))
          )
          fwrite_safe(
            lc[
              all_cause_mrtl > 0,
              .("popsize" = sum(wt_esp), LE = weighted.mean(age, wt_esp)),
              keyby = strata_noagegrp
            ],
            private$output_dir(paste0("summaries/", mcaggr, "le_esp", ext))
          )
          # Life expectancy at 60 ----

          if (
            self$design$sim_prm$ageL < 60L && self$design$sim_prm$ageH > 60L
          ) {
            # fwrite_safe(lc[all_cause_mrtl > 0 & age > 60, .("popsize" = (.N), LE60 = mean(age)),  keyby = strata_noagegrp],
            #             private$output_dir(paste0("summaries/", "le60_out", ext
            #             )))
            fwrite_safe(
              lc[
                all_cause_mrtl > 0 & age > 60,
                .("popsize" = sum(wt), LE60 = weighted.mean(age, wt)),
                keyby = strata_noagegrp
              ],
              private$output_dir(paste0(
                "summaries/",
                mcaggr,
                "le60_scaled_up",
                ext
              ))
            )
            fwrite_safe(
              lc[
                all_cause_mrtl > 0 & age > 60,
                .("popsize" = sum(wt_esp), LE60 = weighted.mean(age, wt_esp)),
                keyby = strata_noagegrp
              ],
              private$output_dir(paste0("summaries/", mcaggr, "le60_esp", ext))
            )
          }
          # Note: for less aggregation use wtd.mean with popsize i.e le_out[,
          # weighted.mean(LE, popsize), keyby = year]
        }

        # Healthy life expectancy ----
        if ("hle" %in% type) {
          # TODO currently some individuals are counted more than once because
          # disease counter and score can be reduced.
          # Ideally only the first reach to the threshold should be counted
          # fwrite_safe(lc[cms_count == 1L, .("popsize" = (.N), HLE = mean(age)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0("summaries/", "hle_1st_cond_out", ext)))
          # fwrite_safe(lc[cms_count == 1L | (cms_count == 0L & all_cause_mrtl > 0),
          #                .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0(
          #               "summaries/", mcaggr, "hle_old_1st_cond_scaled_up", ext
          #             )))
          # fwrite_safe(lc[cms_count == 1L | (cms_count == 0L & all_cause_mrtl > 0),
          #                .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0("summaries/", mcaggr, "hle_old_1st_cond_esp", ext
          #             )))

          # # fwrite_safe(lc[cmsmm1.5_prvl == 1L, .("popsize" = (.N), HLE = mean(age)),
          # #                keyby = strata_noagegrp],
          # #             private$output_dir(paste0("summaries/", "hle_cmsmm1.5_out", ext)))
          # fwrite_safe(lc[cmsmm1.5_prvl == 1L | (cmsmm1.5_prvl == 0L & all_cause_mrtl > 0),
          #                .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0(
          #               "summaries/", mcaggr, "hle_old_cmsmm1.5_scaled_up", ext
          #             )))
          # fwrite_safe(lc[cmsmm1.5_prvl == 1L | (cmsmm1.5_prvl == 0L & all_cause_mrtl > 0),
          #                .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0("summaries/", mcaggr, "hle_old_cmsmm1.5_esp", ext
          #             )))

          fwrite_safe(
            lc[
              cms_count == 0L,
              .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
              keyby = strata_noagegrp
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "hle_1st_cond_scaled_up",
              ext
            ))
          )
          fwrite_safe(
            lc[
              cms_count == 0L,
              .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
              keyby = strata_noagegrp
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "hle_1st_cond_esp",
              ext
            ))
          )

          # fwrite_safe(lc[cmsmm1.5_prvl == 1L, .("popsize" = (.N), HLE = mean(age)),
          #                keyby = strata_noagegrp],
          #             private$output_dir(paste0("summaries/", "hle_cmsmm1.5_out", ext)))
          fwrite_safe(
            lc[
              cmsmm1.5_prvl == 0L,
              .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
              keyby = strata_noagegrp
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "hle_cmsmm1.5_scaled_up",
              ext
            ))
          )
          fwrite_safe(
            lc[
              cmsmm1.5_prvl == 0L,
              .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
              keyby = strata_noagegrp
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "hle_cmsmm1.5_esp",
              ext
            ))
          )
        }

        # Disease characteristics----
        if ("dis_char" %in% type) {
          nm <- grep("_prvl$", names(lc), value = TRUE)

          # tt <- rbindlist(lapply(nm, function(x) {
          #   # sr are the rows the 1st episode occurs per pid
          #   # Need to be sorted on year
          #   sr <- lc[get(x) > 0L, .I[match(1L, get(x))], by = .(pid, scenario)]$V1
          #   sr <- sr[!is.na(sr)]
          #   lc[sr, age_onset := age] # age at 1st ever event
          #
          #   lc[get(x) > 0L, .("disease" = gsub("_prvl$", "", x),
          #                     "cases" = (.N),
          #                     "mean_age_incd" = mean(age[get(x) == 1L]),
          #                     "mean_age_1st_onset" = mean(age_onset, na.rm = TRUE),
          #                     "mean_age_prvl" = mean(age),
          #                     "mean_duration" = mean(get(x)), # Note get(x) very slow here. Implementation with .SDcols also slow because of cases
          #                     "mean_cms_score" = mean(cms_score),
          #                     "mean_cms_count" = mean(cms_count)),
          #      keyby = strata_noagegrp]
          #   lc[, age_onset := NULL]
          # }))
          # tt <- rbindlist(lapply(nm, function(x) {
          #   lc[get(x) > 0L, lapply(.SD, mean), .SDcols = c(x, "age", "cms_score", "cms_count"), keyby = strata_noagegrp]
          # })) # Fast but without cases
          # tt <-
          #   dcast(tt, as.formula(paste0(paste(strata_noagegrp, collapse = "+"), "~disease")),
          #         fill = 0L, value.var = c("cases", "mean_duration", "mean_age_incd",
          #                                  "mean_age_prvl", "mean_cms_score",
          #                                  "mean_cms_count"))
          # fwrite_safe(tt,
          #             private$output_dir(paste0("summaries/", "dis_characteristics_out", ext
          #             )))

          tt <- rbindlist(lapply(nm, function(x) {
            # sr are the rows the 1st episode occurs per pid
            # Need to be sorted on year
            sr <- lc[
              get(x) > 0L,
              .I[match(1L, get(x))],
              by = .(pid, scenario)
            ]$V1
            sr <- sr[!is.na(sr)]
            lc[, wt1st := 0]
            lc[sr, `:=`(age_onset = age, wt1st = wt)] # age at 1st ever event

            ans <- lc[
              get(x) > 0L,
              .(
                "disease" = gsub("_prvl$", "", x),
                "cases" = sum(wt),
                "mean_age_incd" = weighted.mean(
                  age[get(x) == 1L],
                  wt[get(x) == 1L]
                ),
                "mean_age_1st_onset" = weighted.mean(
                  age_onset,
                  wt1st,
                  na.rm = TRUE
                ),

                "mean_age_prvl" = weighted.mean(age, wt),
                "mean_duration" = weighted.mean(get(x), wt), # Note get(x) very slow here. Implementation with .SDcols also slow because of cases
                "mean_cms_score" = weighted.mean(cms_score, wt),
                "mean_cms_count" = weighted.mean(cms_count, wt)
              ),
              keyby = strata_noagegrp
            ]
            lc[, c("age_onset", "wt1st") := NULL]
            ans
          }))
          tt <-
            dcast(
              tt,
              as.formula(paste0(
                paste(strata_noagegrp, collapse = "+"),
                "~disease"
              )),
              fill = 0L,
              value.var = c(
                "cases",
                "mean_duration",
                "mean_age_incd",
                "mean_age_1st_onset",
                "mean_age_prvl",
                "mean_cms_score",
                "mean_cms_count"
              )
            )
          fwrite_safe(
            tt,
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "dis_characteristics_scaled_up",
              ext
            ))
          )

          tt <- rbindlist(lapply(nm, function(x) {
            # sr are the rows the 1st episode occurs per pid
            # Need to be sorted on year
            sr <- lc[
              get(x) > 0L,
              .I[match(1L, get(x))],
              by = .(pid, scenario)
            ]$V1
            sr <- sr[!is.na(sr)]
            lc[, wt1st := 0]
            lc[sr, `:=`(age_onset = age, wt1st = wt_esp)] # age at 1st ever event

            ans <- lc[
              get(x) > 0L,
              .(
                "disease" = gsub("_prvl$", "", x),
                "cases" = sum(wt_esp),
                "mean_age_incd" = weighted.mean(
                  age[get(x) == 1L],
                  wt_esp[get(x) == 1L]
                ),
                "mean_age_1st_onset" = weighted.mean(
                  age_onset,
                  wt1st,
                  na.rm = TRUE
                ),
                "mean_age_prvl" = weighted.mean(age, wt_esp),
                "mean_duration" = weighted.mean(get(x), wt_esp), # Note get(x) very slow here. Implementation with .SDcols also slow because of cases
                "mean_cms_score" = weighted.mean(cms_score, wt_esp),
                "mean_cms_count" = weighted.mean(cms_count, wt_esp)
              ),
              keyby = strata_noagegrp
            ]
            lc[, c("age_onset", "wt1st") := NULL]
            ans
          }))
          tt <-
            dcast(
              tt,
              as.formula(paste0(
                paste(strata_noagegrp, collapse = "+"),
                "~disease"
              )),
              fill = 0L,
              value.var = c(
                "cases",
                "mean_duration",
                "mean_age_incd",
                "mean_age_1st_onset",
                "mean_age_prvl",
                "mean_cms_score",
                "mean_cms_count"
              )
            )
          fwrite_safe(
            tt,
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "dis_characteristics_esp",
              ext
            ))
          )
          rm(tt)
        }

        # prvl ----
        if ("prvl" %in% type) {
          # Note for mortality this exports qx directly. mx is defined as the
          # number of deaths during the year divided by the average number alive
          # during the year, i.e. This differs slightly from qx , which is the
          # number of deaths during the year divided by the number alive at the
          # beginning of the year.
          # fwrite_safe(lc[, c("popsize" = (.N),
          #                    lapply(.SD, function(x) sum(x > 0))),
          #                .SDcols = patterns("_prvl$"), keyby = strata],
          #             private$output_dir(paste0("summaries/", "prvl_out", ext
          #             )))
          fwrite_safe(
            lc[,
              c(
                "popsize" = sum(wt),
                lapply(.SD, function(x, wt) sum((x > 0) * wt), wt)
              ),
              .SDcols = patterns("_prvl$"),
              keyby = strata
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "prvl_scaled_up",
              ext
            ))
          )
          fwrite_safe(
            lc[,
              c(
                "popsize" = sum(wt_esp),
                lapply(.SD, function(x, wt) sum((x > 0) * wt), wt_esp)
              ),
              .SDcols = patterns("_prvl$"),
              keyby = strata
            ],
            private$output_dir(paste0("summaries/", mcaggr, "prvl_esp", ext))
          )
        }

        # incd ----
        if ("incd" %in% type) {
          # NOTE incd includes prevalent cases in denominator
          # fwrite_safe(lc[, c("popsize" = (.N),
          #                    lapply(.SD, function(x) sum(x == 1))),
          #                .SDcols = patterns("_prvl$"), keyby = strata],
          #             private$output_dir(paste0("summaries/", "incd_out", ext
          #             )))
          incdtbl <- lc[,
            c(
              "popsize" = sum(wt),
              lapply(.SD, function(x, wt) sum((x == 1) * wt), wt)
            ),
            .SDcols = patterns("_prvl$"),
            keyby = strata
          ]
          nm <- grep("_prvl$", names(incdtbl), value = TRUE)
          setnames(incdtbl, nm, gsub("_prvl$", "_incd", nm))
          fwrite_safe(
            incdtbl,
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "incd_scaled_up",
              ext
            ))
          )

          incdtbl <- lc[,
            c(
              "popsize" = sum(wt_esp),
              lapply(.SD, function(x, wt) sum((x == 1) * wt), wt_esp)
            ),
            .SDcols = patterns("_prvl$"),
            keyby = strata
          ]
          nm <- grep("_prvl$", names(incdtbl), value = TRUE)
          setnames(incdtbl, nm, gsub("_prvl$", "_incd", nm))
          fwrite_safe(
            incdtbl,
            private$output_dir(paste0("summaries/", mcaggr, "incd_esp", ext))
          )

          rm(incdtbl, nm)
        }

        # mrtl ----
        if ("mrtl" %in% type) {
          # fwrite_safe(lc[, .("popsize" = (.N),
          #                    "all_cause_mrtl" = sum(all_cause_mrtl > 0)),
          #                keyby = strata],
          #             private$output_dir(paste0("summaries/", "mrtl_out", ext
          #             )))
          fwrite_safe(
            lc[,
              .(
                "popsize" = sum(wt),
                "all_cause_mrtl" = sum((all_cause_mrtl > 0) * wt)
              ),
              keyby = strata
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "mrtl_scaled_up",
              ext
            ))
          )
          fwrite_safe(
            lc[,
              .(
                "popsize" = sum(wt_esp),
                "all_cause_mrtl" = sum((all_cause_mrtl > 0) * wt_esp)
              ),
              keyby = strata
            ],
            private$output_dir(paste0("summaries/", mcaggr, "mrtl_esp", ext))
          )
        }

        # disease specific mortality ----
        if ("dis_mrtl" %in% type) {
          # dis_mrtl_out <-
          #   dcast(
          #     lc[, .("deaths" = (.N)),
          #        keyby = c(strata, "all_cause_mrtl")],
          #     formula = as.formula(paste0(
          #       paste(strata, collapse = "+"), "~all_cause_mrtl"
          #     )),
          #     fill = 0L,
          #     value.var = "deaths"
          #   )
          #
          # setnames(dis_mrtl_out, as.character(private$death_codes),
          #          names(private$death_codes), skip_absent = TRUE)
          # dis_mrtl_out[, `:=` (
          #   popsize = Reduce(`+`, .SD),
          #   alive = NULL
          # ), .SDcols = !strata]
          # fwrite_safe(dis_mrtl_out,
          #             private$output_dir(paste0("summaries/", "dis_mrtl_out", ext
          #             )))

          dis_mrtl_out <- # scale up
            dcast(
              lc[, .("deaths" = sum(wt)), keyby = c(strata, "all_cause_mrtl")],
              formula = as.formula(paste0(
                paste(strata, collapse = "+"),
                "~all_cause_mrtl"
              )),
              fill = 0L,
              value.var = "deaths"
            )

          setnames(
            dis_mrtl_out,
            as.character(private$death_codes),
            paste0(names(private$death_codes), "_deaths"),
            skip_absent = TRUE
          )
          dis_mrtl_out[,
            `:=`(
              popsize = Reduce(`+`, .SD), # it includes alive so it is the pop at the start of the year
              alive_deaths = NULL
            ),
            .SDcols = !strata
          ]
          fwrite_safe(
            dis_mrtl_out,
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "dis_mrtl_scaled_up",
              ext
            ))
          )

          dis_mrtl_out <- # scale up esp
            dcast(
              lc[,
                .("deaths" = sum(wt_esp)),
                keyby = c(strata, "all_cause_mrtl")
              ],
              formula = as.formula(paste0(
                paste(strata, collapse = "+"),
                "~all_cause_mrtl"
              )),
              fill = 0L,
              value.var = "deaths"
            )

          setnames(
            dis_mrtl_out,
            as.character(private$death_codes),
            paste0(names(private$death_codes), "_deaths"),
            skip_absent = TRUE
          )
          dis_mrtl_out[,
            `:=`(
              popsize = Reduce(`+`, .SD),
              alive_deaths = NULL
            ),
            .SDcols = !strata
          ]
          fwrite_safe(
            dis_mrtl_out,
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "dis_mrtl_esp",
              ext
            ))
          )
          rm(dis_mrtl_out)
        }

        # All-cause mrtl by disease ----
        if ("allcause_mrtl_by_dis" %in% type) {
          nm <- grep("_prvl$", names(lc), value = TRUE)

          # tt <- lapply(nm, function(x) {
          #   lc[get(x) > 0L, .("disease" = gsub("_prvl$", "", x), "cases" = (.N), "deaths" = sum(all_cause_mrtl > 0)), keyby = strata]
          # })
          # tt <-
          # dcast(rbindlist(tt), as.formula(paste0(paste(strata, collapse = "+"), "~disease")),
          #       fill = 0L, value.var = c("deaths", "cases"))
          # fwrite_safe(tt,
          #             private$output_dir(paste0("summaries/", "all_cause_mrtl_by_dis_out", ext
          #             )))

          tt <- lapply(nm, function(x) {
            lc[
              get(x) > 0L,
              .(
                "disease" = gsub("_prvl$", "", x),
                "cases" = sum(wt),
                "deaths" = sum(wt * (all_cause_mrtl > 0))
              ),
              keyby = strata
            ]
          })
          tt <-
            dcast(
              rbindlist(tt),
              as.formula(paste0(paste(strata, collapse = "+"), "~disease")),
              fill = 0L,
              value.var = c("deaths", "cases")
            )
          fwrite_safe(
            tt,
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "all_cause_mrtl_by_dis_scaled_up",
              ext
            ))
          )

          tt <- lapply(nm, function(x) {
            lc[
              get(x) > 0L,
              .(
                "disease" = gsub("_prvl$", "", x),
                "cases" = sum(wt_esp),
                "deaths" = sum(wt_esp * (all_cause_mrtl > 0))
              ),
              keyby = strata
            ]
          })
          tt <-
            dcast(
              rbindlist(tt),
              as.formula(paste0(paste(strata, collapse = "+"), "~disease")),
              fill = 0L,
              value.var = c("deaths", "cases")
            )
          fwrite_safe(
            tt,
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "all_cause_mrtl_by_dis_esp",
              ext
            ))
          )
          rm(tt)
        }

        # CMS mean ----
        if ("cms" %in% type) {
          # fwrite_safe(lc[, .("popsize" = (.N), cms_score = mean(cms_score)), keyby = strata],
          #             private$output_dir(paste0("summaries/", "cms_score_out", ext
          #             )))
          fwrite_safe(
            lc[,
              .("popsize" = sum(wt), cms_score = weighted.mean(cms_score, wt)),
              keyby = strata
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "cms_score_scaled_up",
              ext
            ))
          )
          fwrite_safe(
            lc[,
              .("popsize" = sum(wt), cms_score = weighted.mean(cms_score, wt)),
              keyby = strata_age
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "cms_score_by_age_scaled_up",
              ext
            ))
          )

          fwrite_safe(
            lc[,
              .(
                "popsize" = sum(wt_esp),
                cms_score = weighted.mean(cms_score, wt_esp)
              ),
              keyby = strata
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "cms_score_esp",
              ext
            ))
          )

          fwrite_safe(
            lc[,
              .(
                "popsize" = sum(wt_esp),
                cms_score = weighted.mean(cms_score, wt_esp)
              ),
              keyby = strata_age
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "cms_score_by_age_esp",
              ext
            ))
          )

          # CMS count ----
          # fwrite_safe(lc[, .("popsize" = (.N), cms_count = mean(cms_count)), keyby = strata],
          #             private$output_dir(paste0("summaries/", "cms_count_out", ext
          #             )))
          fwrite_safe(
            lc[,
              .("popsize" = sum(wt), cms_count = weighted.mean(cms_count, wt)),
              keyby = strata
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "cms_count_scaled_up",
              ext
            ))
          )
          fwrite_safe(
            lc[,
              .(
                "popsize" = sum(wt_esp),
                cms_count = weighted.mean(cms_count, wt_esp)
              ),
              keyby = strata
            ],
            private$output_dir(paste0(
              "summaries/",
              mcaggr,
              "cms_count_esp",
              ext
            ))
          )
        }

        if (!self$design$sim_prm$keep_lifecourse) {
          file.remove(pth)
        }

        return(invisible(self))
      },

      # Special deep copy for data.table. Use POP$clone(deep = TRUE) to
      # dispatch. Otherwise a reference is created
      # @param name A character string representing the name or identifier of the object being cloned.
      # @param value The object to be cloned.
      #
      # @details
      # The function checks the class of the input object and performs the appropriate cloning method. For data tables, it uses data.table::copy, and for R6 objects, it calls the clone() method. For other object types, a shallow copy is returned.
      #
      # @return
      # A deep clone of the input object.
      #
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

      # calc_QALYs ----
      # Memory-optimised version of calc_QALYs using DuckDB SQL.
      # Creates a temporary view named 'output_view_name' in DuckDB,
      # which is the 'input_table_name' (filtered by mcaggr) augmented with EQ5D5L and HUI3 columns.
      calc_QALYs = function(
        duckdb_con,
        mcaggr,
        input_table_name,
        output_view_name,
        include_non_significant = FALSE
      ) {
        eq5d5l_expr <- "
          0.989
          + CASE agegrp
              WHEN '20-24' THEN -0.018 WHEN '25-29' THEN -0.018
              WHEN '30-34' THEN -0.019 WHEN '35-39' THEN -0.019
              WHEN '40-44' THEN -0.018 WHEN '45-49' THEN -0.018
              WHEN '50-54' THEN -0.028 WHEN '55-59' THEN -0.028
              WHEN '60-64' THEN -0.021 WHEN '65-69' THEN -0.021
              WHEN '70-74' THEN -0.057 WHEN '75-79' THEN -0.057
              WHEN '80-84' THEN -0.129 WHEN '85-89' THEN -0.129
              WHEN '90-94' THEN -0.129 WHEN '95-99' THEN -0.129
              ELSE 0.0
            END
          + CASE WHEN sex = 'women' THEN -0.011 ELSE 0.0 END
          + CASE WHEN chd_prvl = 0 THEN 0.0 ELSE -0.073 END
          + CASE WHEN stroke_prvl = 0 THEN 0.0 ELSE -0.265 END
          + CASE WHEN t2dm_prvl = 0 THEN 0.0 ELSE -0.046 END
        "

        hui3_expr <- "
          0.897
          + CASE agegrp
              WHEN '20-24' THEN -0.023 WHEN '25-29' THEN -0.023
              WHEN '30-34' THEN -0.018 WHEN '35-39' THEN -0.018
              WHEN '40-44' THEN -0.004 WHEN '45-49' THEN -0.004
              WHEN '50-54' THEN -0.021 WHEN '55-59' THEN -0.021
              WHEN '60-64' THEN -0.013 WHEN '65-69' THEN -0.013
              WHEN '70-74' THEN -0.042 WHEN '75-79' THEN -0.042
              WHEN '80-84' THEN -0.145 WHEN '85-89' THEN -0.145
              WHEN '90-94' THEN -0.145 WHEN '95-99' THEN -0.145
              ELSE 0.0
            END
          + CASE WHEN sex = 'women' THEN 0.011 ELSE 0.0 END
          + CASE WHEN chd_prvl = 0 THEN 0.0 ELSE -0.081 END
          + CASE WHEN stroke_prvl = 0 THEN 0.0 ELSE -0.293 END
          + CASE WHEN t2dm_prvl = 0 THEN 0.0 ELSE -0.055 END
        "

        if (!include_non_significant) {
          # Original logic: if (!include_non_significant == TRUE)
          eq5d5l_expr <- paste0(
            eq5d5l_expr,
            " + CASE WHEN htn_prvl = 0 THEN 0.0 ELSE -0.005 END",
            " + CASE WHEN obesity_prvl = 0 THEN 0.0 ELSE -0.034 END"
          )
          hui3_expr <- paste0(
            hui3_expr,
            " + CASE WHEN htn_prvl = 0 THEN 0.0 ELSE -0.006 END",
            " + CASE WHEN obesity_prvl = 0 THEN 0.0 ELSE 0.019 END"
          )
        }

        create_view_sql <- sprintf(
          "
          CREATE OR REPLACE TEMP VIEW %s AS
          SELECT
            *,
            (%s) AS EQ5D5L,
            (%s) AS HUI3
          FROM %s
          WHERE mc = %d;
        ",
          output_view_name,
          eq5d5l_expr,
          hui3_expr,
          input_table_name,
          mcaggr
        )

        private$execute_sql(duckdb_con, create_view_sql)

        NULL
      },

      # calc_costs ----
      # Memory-optimised version of calc_costs using DuckDB SQL.
      # Creates a temporary view named 'output_view_name' in DuckDB,
      # which is the 'input_table_name' (filtered by mcaggr) augmented with calculated cost columns.
      # NOTE: This is adapted from the Japan model. You may need to update the
      # input file paths and cost parameters to match England-specific data.
      calc_costs = function(
        duckdb_con,
        mcaggr,
        input_table_name,
        output_view_name
      ) {
        # get scenario names
        scnams <- gsub(
          "^scenario=",
          "",
          list.dirs(
            private$output_dir(file.path("lifecourse", paste0("mc=", mcaggr))),
            full.names = FALSE,
            recursive = FALSE
          )
        )

        # --- Inflation Factors ---
        prod_informal_inflation_factor <- 1.025
        direct_costs_inflation_factor <- 99.6 / 99.7

        # --- Step 1: Create baseline aggregation views using SQL only ---
        base_agg_sql <- "
          CREATE OR REPLACE TEMP VIEW %s AS
          SELECT agegrp, sex, ROUND(SUM(CASE WHEN %s THEN wt ELSE 0 END)) AS V1
          FROM %s WHERE year = %d AND scenario = 'sc0' AND mc = %d GROUP BY agegrp, sex
          "

        # Create all baseline aggregation views
        aggregation_configs <- list(
          list("chd_prvl_2016_agg_view", "chd_dgns > 0", 2016),
          list("chd_prvl_2019_agg_view", "chd_dgns > 0", 2019),
          list("stroke_prvl_2016_agg_view", "stroke_dgns > 0", 2016),
          list("stroke_prvl_2019_agg_view", "stroke_dgns > 0", 2019),
          list("chd_mrtl_2016_initial_view", "all_cause_mrtl = 2", 2016),
          list("stroke_mrtl_2016_initial_view", "all_cause_mrtl = 3", 2016)
        )

        for (config in aggregation_configs) {
          private$execute_sql(
            duckdb_con,
            sprintf(
              base_agg_sql,
              config[[1]],
              config[[2]],
              input_table_name,
              config[[3]],
              mcaggr
            ),
            config[[1]]
          )
        }

        # --- Step 2: Memory-efficient mortality data handling ---
        # TODO: Update file paths for England-specific data
        # Load only essential columns and filter immediately

        # Load observed population (minimal columns)
        obs_pop_2016 <- read_fst(
          "inputs/pop_estimates/observed_population_england.fst",
          columns = c("year", "age", "sex", "pops"),
          as.data.table = TRUE
        )[year == 2016L]

        # Process CHD mortality efficiently
        chd_ftlt_2016 <- read_fst(
          "inputs/disease_burden/chd_ftlt.fst",
          columns = c("year", "age", "sex", "mu2"),
          as.data.table = TRUE
        )[year == 2016]

        chd_joined <- chd_ftlt_2016[
          obs_pop_2016,
          on = c("age", "sex"),
          nomatch = 0L
        ][, `:=`(
          deaths_calc = mu2 * pops,
          agegrp = fcase(
            age %between% c(30, 34),
            "30-34",
            age %between% c(35, 39),
            "35-39",
            age %between% c(40, 44),
            "40-44",
            age %between% c(45, 49),
            "45-49",
            age %between% c(50, 54),
            "50-54",
            age %between% c(55, 59),
            "55-59",
            age %between% c(60, 64),
            "60-64",
            age %between% c(65, 69),
            "65-69",
            age %between% c(70, 74),
            "70-74",
            age %between% c(75, 79),
            "75-79",
            age %between% c(80, 84),
            "80-84",
            age %between% c(85, 89),
            "85-89",
            age %between% c(90, 94),
            "90-94",
            age >= 95,
            "95-99",
            default = NA_character_
          )
        )][
          !is.na(agegrp),
          .(calculated_deaths = round(sum(deaths_calc))),
          keyby = .(agegrp, sex)
        ]

        # Register minimal table
        dbWriteTable(
          duckdb_con,
          "chd_ftlt_ext_2016_table",
          chd_joined,
          overwrite = TRUE
        )
        rm(chd_joined) # Immediate cleanup

        # Process stroke mortality efficiently
        stroke_ftlt_2016 <- read_fst(
          "inputs/disease_burden/stroke_ftlt.fst",
          columns = c("year", "age", "sex", "mu2"),
          as.data.table = TRUE
        )[year == 2016]

        stroke_joined <- stroke_ftlt_2016[
          obs_pop_2016,
          on = c("age", "sex"),
          nomatch = 0L
        ][, `:=`(
          deaths_calc = mu2 * pops,
          agegrp = fcase(
            age %between% c(30, 34),
            "30-34",
            age %between% c(35, 39),
            "35-39",
            age %between% c(40, 44),
            "40-44",
            age %between% c(45, 49),
            "45-49",
            age %between% c(50, 54),
            "50-54",
            age %between% c(55, 59),
            "55-59",
            age %between% c(60, 64),
            "60-64",
            age %between% c(65, 69),
            "65-69",
            age %between% c(70, 74),
            "70-74",
            age %between% c(75, 79),
            "75-79",
            age %between% c(80, 84),
            "80-84",
            age %between% c(85, 89),
            "85-89",
            age %between% c(90, 94),
            "90-94",
            age >= 95,
            "95-99",
            default = NA_character_
          )
        )][
          !is.na(agegrp),
          .(calculated_deaths = round(sum(deaths_calc))),
          keyby = .(agegrp, sex)
        ]

        dbWriteTable(
          duckdb_con,
          "stroke_ftlt_ext_2016_table",
          as.data.frame(stroke_joined),
          overwrite = TRUE
        )
        rm(stroke_joined) # Immediate cleanup

        # Cleanup large intermediate objects immediately
        rm(obs_pop_2016, chd_ftlt_2016, stroke_ftlt_2016)

        # Update mortality views
        private$execute_sql(
          duckdb_con,
          "
          CREATE OR REPLACE TEMP VIEW chd_mrtl_2016_agg_view AS
          SELECT i.agegrp, i.sex, 
                 CASE WHEN i.V1 = 0 THEN COALESCE(f.calculated_deaths, i.V1) ELSE i.V1 END AS V1
          FROM chd_mrtl_2016_initial_view i 
          LEFT JOIN chd_ftlt_ext_2016_table f ON i.agegrp = f.agegrp AND i.sex = f.sex
        ",
          "chd_mrtl_2016_agg_view"
        )

        private$execute_sql(
          duckdb_con,
          "
          CREATE OR REPLACE TEMP VIEW stroke_mrtl_2016_agg_view AS
          SELECT i.agegrp, i.sex, 
                 CASE WHEN i.V1 = 0 THEN COALESCE(f.calculated_deaths, i.V1) ELSE i.V1 END AS V1
          FROM stroke_mrtl_2016_initial_view i 
          LEFT JOIN stroke_ftlt_ext_2016_table f ON i.agegrp = f.agegrp AND i.sex = f.sex  
        ",
          "stroke_mrtl_2016_agg_view"
        )

        # --- Step 3: Memory-efficient cost parameter calculation ---
        # Create parameter tables directly in SQL to avoid R object creation

        # Employee parameters - create as SQL view to avoid R data.table
        # TODO: Update employee counts for England
        private$execute_sql(
          duckdb_con,
          "
          CREATE OR REPLACE TEMP VIEW employee_params_view AS
          SELECT agegrp, sex, CAST(employees AS DOUBLE) AS employees FROM (VALUES
            ('30-34', 'men', 1683780), ('35-39', 'men', 1829610), ('40-44', 'men', 2174550), ('45-49', 'men', 2057710),
            ('50-54', 'men', 1702470), ('55-59', 'men', 1425510), ('60-64', 'men', 963430), ('65-69', 'men', 369640),
            ('70-74', 'men', 106850), ('75-79', 'men', 0), ('80-84', 'men', 0), ('85-89', 'men', 0),
            ('90-94', 'men', 0), ('95-99', 'men', 0),
            ('30-34', 'women', 919700), ('35-39', 'women', 894770), ('40-44', 'women', 1049490), ('45-49', 'women', 1037140),
            ('50-54', 'women', 854970), ('55-59', 'women', 685040), ('60-64', 'women', 376370), ('65-69', 'women', 132470),
            ('70-74', 'women', 44050), ('75-79', 'women', 0), ('80-84', 'women', 0), ('85-89', 'women', 0),
            ('90-94', 'women', 0), ('95-99', 'women', 0)
          ) AS t(agegrp, sex, employees)
        ",
          "employee_params_view"
        )

        # CHD informal care parameters
        # TODO: Update informal care hours for England
        private$execute_sql(
          duckdb_con,
          "
          CREATE OR REPLACE TEMP VIEW chd_infm_care_view AS
          SELECT agegrp, sex, CAST(infm_care_hrs AS DOUBLE) AS infm_care_hrs FROM (VALUES
            ('30-34', 'men', 0.030), ('35-39', 'men', 0.030), ('40-44', 'men', 0.030), ('45-49', 'men', 0.030),
            ('50-54', 'men', 0.030), ('55-59', 'men', 0.030), ('60-64', 'men', 0.030), ('65-69', 'men', 0.200),
            ('70-74', 'men', 0.200), ('75-79', 'men', 0.200), ('80-84', 'men', 0), ('85-89', 'men', 0),
            ('90-94', 'men', 0), ('95-99', 'men', 0),
            ('30-34', 'women', 0.030), ('35-39', 'women', 0.030), ('40-44', 'women', 0.030), ('45-49', 'women', 0.030),
            ('50-54', 'women', 0.030), ('55-59', 'women', 0.030), ('60-64', 'women', 0.030), ('65-69', 'women', 0.200),
            ('70-74', 'women', 0.200), ('75-79', 'women', 0.200), ('80-84', 'women', 0), ('85-89', 'women', 0),
            ('90-94', 'women', 0), ('95-99', 'women', 0)
          ) AS t(agegrp, sex, infm_care_hrs)
        ",
          "chd_infm_care_view"
        )

        # Stroke informal care parameters
        # TODO: Update informal care hours for England
        private$execute_sql(
          duckdb_con,
          "
          CREATE OR REPLACE TEMP VIEW stroke_infm_care_view AS
          SELECT agegrp, sex, CAST(infm_care_hrs AS DOUBLE) AS infm_care_hrs FROM (VALUES
            ('30-34', 'men', 5.20), ('35-39', 'men', 5.20), ('40-44', 'men', 5.20), ('45-49', 'men', 5.20),
            ('50-54', 'men', 5.20), ('55-59', 'men', 5.20), ('60-64', 'men', 5.20), ('65-69', 'men', 5.03),
            ('70-74', 'men', 5.03), ('75-79', 'men', 5.03), ('80-84', 'men', 9.23), ('85-89', 'men', 9.23),
            ('90-94', 'men', 9.23), ('95-99', 'men', 9.23),
            ('30-34', 'women', 5.20), ('35-39', 'women', 5.20), ('40-44', 'women', 5.20), ('45-49', 'women', 5.20),
            ('50-54', 'women', 5.20), ('55-59', 'women', 5.20), ('60-64', 'women', 5.20), ('65-69', 'women', 5.03),
            ('70-74', 'women', 5.03), ('75-79', 'women', 5.03), ('80-84', 'women', 9.23), ('85-89', 'women', 9.23),
            ('90-94', 'women', 9.23), ('95-99', 'women', 9.23)
          ) AS t(agegrp, sex, infm_care_hrs)
        ",
          "stroke_infm_care_view"
        )

        # Direct cost parameters
        # TODO: Update direct costs for England
        private$execute_sql(
          duckdb_con,
          "
          CREATE OR REPLACE TEMP VIEW chd_direct_tcost_view AS
          SELECT agegrp2, sex, CAST(tcost_val AS DOUBLE) AS tcost_val FROM (VALUES
            ('30-44', 'men', 10300000000.0), ('45-64', 'men', 121000000000.0), ('65-69', 'men', 72300000000.0),
            ('70-74', 'men', 90100000000.0), ('75-99', 'men', 197000000000.0),
            ('30-44', 'women', 2500000000.0), ('45-64', 'women', 22500000000.0), ('65-69', 'women', 18900000000.0),
            ('70-74', 'women', 30300000000.0), ('75-99', 'women', 132800000000.0)
          ) AS t(agegrp2, sex, tcost_val)
        ",
          "chd_direct_tcost_view"
        )

        private$execute_sql(
          duckdb_con,
          "
          CREATE OR REPLACE TEMP VIEW stroke_direct_tcost_view AS
          SELECT agegrp2, sex, CAST(tcost_val AS DOUBLE) AS tcost_val FROM (VALUES
            ('30-44', 'men', 24900000000.0), ('45-64', 'men', 186400000000.0), ('65-69', 'men', 109000000000.0),
            ('70-74', 'men', 144100000000.0), ('75-99', 'men', 465600000000.0),
            ('30-44', 'women', 18000000000.0), ('45-64', 'women', 106800000000.0), ('65-69', 'women', 60900000000.0),
            ('70-74', 'women', 95700000000.0), ('75-99', 'women', 606800000000.0)
          ) AS t(agegrp2, sex, tcost_val)
        ",
          "stroke_direct_tcost_view"
        )

        # Optimised cost parameter calculation - single SQL statement approach
        cost_param_sql <- "
          CREATE OR REPLACE TEMP VIEW %s AS
          WITH joined_data AS (
            SELECT p.agegrp, p.sex, p.%s AS factor_col, agg.V1 
            FROM %s p 
            JOIN %s agg ON p.agegrp = agg.agegrp AND p.sex = agg.sex
          ), 
          weighted_data AS (
            SELECT agegrp, sex, (factor_col * V1) AS weighted_factor 
            FROM joined_data
          ), 
          total_weighted_sum AS (
            SELECT SUM(weighted_factor) AS total_wt_sum FROM weighted_data
          )
          SELECT wd.agegrp, wd.sex, 
                 (%.2f * wd.weighted_factor / NULLIF(tws.total_wt_sum, 0)) * %.6f / NULLIF(jd.V1, 0) AS cost_param
          FROM weighted_data wd
          CROSS JOIN total_weighted_sum tws
          JOIN joined_data jd ON wd.agegrp = jd.agegrp AND wd.sex = jd.sex
        "

        # Create all cost parameter views efficiently
        # TODO: Update total cost values for England
        cost_configs <- list(
          list(
            "chd_prvl_prdv_cost_param_view",
            "employees",
            "employee_params_view",
            "chd_prvl_2016_agg_view",
            141000000000.00
          ),
          list(
            "stroke_prvl_prdv_cost_param_view",
            "employees",
            "employee_params_view",
            "stroke_prvl_2016_agg_view",
            322000000000.00
          ),
          list(
            "chd_mrtl_prdv_cost_param_view",
            "employees",
            "employee_params_view",
            "chd_mrtl_2016_agg_view",
            2257000000000.00
          ),
          list(
            "stroke_mrtl_prdv_cost_param_view",
            "employees",
            "employee_params_view",
            "stroke_mrtl_2016_agg_view",
            1352000000000.00
          ),
          list(
            "chd_informal_cost_param_view",
            "infm_care_hrs",
            "chd_infm_care_view",
            "chd_prvl_2016_agg_view",
            291000000000.00
          ),
          list(
            "stroke_informal_cost_param_view",
            "infm_care_hrs",
            "stroke_infm_care_view",
            "stroke_prvl_2016_agg_view",
            1651000000000.00
          )
        )

        for (config in cost_configs) {
          private$execute_sql(
            duckdb_con,
            sprintf(
              cost_param_sql,
              config[[1]],
              config[[2]],
              config[[3]],
              config[[4]],
              config[[5]],
              prod_informal_inflation_factor
            ),
            config[[1]]
          )
        }

        # Direct cost parameters with optimised SQL
        direct_cost_sql <- "
          CREATE OR REPLACE TEMP VIEW %s AS
          WITH lc_with_agegrp2 AS (
            SELECT agegrp, sex, V1,
              CASE 
                WHEN agegrp IN ('30-34', '35-39', '40-44') THEN '30-44'
                WHEN agegrp IN ('45-49', '50-54', '55-59', '60-64') THEN '45-64'
                WHEN agegrp = '65-69' THEN '65-69'
                WHEN agegrp = '70-74' THEN '70-74'
                ELSE '75-99'
              END AS agegrp2
            FROM %s
          ), 
          agg_by_agegrp2 AS (
            SELECT agegrp2, sex, SUM(V1) AS V1_sum 
            FROM lc_with_agegrp2 
            GROUP BY agegrp2, sex
          )
          SELECT orig.agegrp, orig.sex, 
                 (tc.tcost_val * %.6f / NULLIF(agg.V1_sum, 0)) AS cost_param
          FROM %s orig
          JOIN lc_with_agegrp2 lwa ON orig.agegrp = lwa.agegrp AND orig.sex = lwa.sex
          JOIN agg_by_agegrp2 agg ON lwa.agegrp2 = agg.agegrp2 AND lwa.sex = agg.sex
          JOIN %s tc ON agg.agegrp2 = tc.agegrp2 AND agg.sex = tc.sex
        "

        private$execute_sql(
          duckdb_con,
          sprintf(
            direct_cost_sql,
            "chd_direct_cost_param_view",
            "chd_prvl_2019_agg_view",
            direct_costs_inflation_factor,
            "chd_prvl_2019_agg_view",
            "chd_direct_tcost_view"
          ),
          "chd_direct_cost_param_view"
        )

        private$execute_sql(
          duckdb_con,
          sprintf(
            direct_cost_sql,
            "stroke_direct_cost_param_view",
            "stroke_prvl_2019_agg_view",
            direct_costs_inflation_factor,
            "stroke_prvl_2019_agg_view",
            "stroke_direct_tcost_view"
          ),
          "stroke_direct_cost_param_view"
        )

        # --- Step 4: Create Final Output View with All Cost Columns ---
        # One view per scenario named paste0(output_view_name, "_", scnams, "_view")
        final_view_creation_sql <- sprintf(
          "
          CREATE OR REPLACE TEMP VIEW %s AS
          WITH base_filtered AS (
            SELECT mc, scenario, year, agegrp, sex, chd_dgns, all_cause_mrtl, stroke_dgns, wt, wt_esp
            FROM %s 
            WHERE mc = %d AND scenario = %s
            ),
            chd_costs AS (
              SELECT agegrp, sex,
                COALESCE(cppc.cost_param, 0) AS chd_prvl_prdv,
                COALESCE(cpmc.cost_param, 0) AS chd_mrtl_prdv,
                COALESCE(cic.cost_param, 0) AS chd_informal,
                COALESCE(cdc.cost_param, 0) AS chd_direct
              FROM chd_prvl_prdv_cost_param_view cppc
              LEFT JOIN chd_mrtl_prdv_cost_param_view cpmc USING(agegrp, sex)
              LEFT JOIN chd_informal_cost_param_view cic USING(agegrp, sex)
              LEFT JOIN chd_direct_cost_param_view cdc USING(agegrp, sex)
            ),
            stroke_costs AS (
              SELECT agegrp, sex,
                COALESCE(sppc.cost_param, 0) AS stroke_prvl_prdv,
                COALESCE(spmc.cost_param, 0) AS stroke_mrtl_prdv,
                COALESCE(sic.cost_param, 0) AS stroke_informal,
                COALESCE(sdc.cost_param, 0) AS stroke_direct
              FROM stroke_prvl_prdv_cost_param_view sppc
              LEFT JOIN stroke_mrtl_prdv_cost_param_view spmc USING(agegrp, sex)
              LEFT JOIN stroke_informal_cost_param_view sic USING(agegrp, sex)
              LEFT JOIN stroke_direct_cost_param_view sdc USING(agegrp, sex)
            ),
            basic_costs AS (
              SELECT
                m.mc, m.scenario, m.year, m.agegrp, m.sex,
                m.wt, m.wt_esp,
              
                -- CHD basic cost components
                CASE WHEN m.chd_dgns > 0 THEN cc.chd_prvl_prdv ELSE 0 END AS chd_prvl_prdv_costs,
                CASE WHEN m.all_cause_mrtl = 2 THEN cc.chd_mrtl_prdv ELSE 0 END AS chd_mrtl_prdv_costs,
                CASE WHEN m.chd_dgns > 0 THEN cc.chd_informal ELSE 0 END AS chd_informal_costs,
                CASE WHEN m.chd_dgns > 0 THEN cc.chd_direct ELSE 0 END AS chd_direct_costs,
              
                -- Stroke basic cost components
                CASE WHEN m.stroke_dgns > 0 THEN sc.stroke_prvl_prdv ELSE 0 END AS stroke_prvl_prdv_costs,
                CASE WHEN m.all_cause_mrtl = 3 THEN sc.stroke_mrtl_prdv ELSE 0 END AS stroke_mrtl_prdv_costs,
                CASE WHEN m.stroke_dgns > 0 THEN sc.stroke_informal ELSE 0 END AS stroke_informal_costs,
                CASE WHEN m.stroke_dgns > 0 THEN sc.stroke_direct ELSE 0 END AS stroke_direct_costs
              
              FROM base_filtered m
              LEFT JOIN chd_costs cc ON m.agegrp = cc.agegrp AND m.sex = cc.sex
              LEFT JOIN stroke_costs sc ON m.agegrp = sc.agegrp AND m.sex = sc.sex
            )
            SELECT
              mc, scenario, year, agegrp, sex, wt, wt_esp,
            
              -- Basic cost components (already calculated)
              chd_prvl_prdv_costs,
              chd_mrtl_prdv_costs,
              chd_informal_costs,
              chd_direct_costs,
              stroke_prvl_prdv_costs,
              stroke_mrtl_prdv_costs,
              stroke_informal_costs,
              stroke_direct_costs,
            
              -- Aggregated productivity costs
              (chd_prvl_prdv_costs + chd_mrtl_prdv_costs) AS chd_productivity_costs,
              (stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs) AS stroke_productivity_costs,
              (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs) AS cvd_productivity_costs,
             
              -- Aggregated informal costs
              (chd_informal_costs + stroke_informal_costs) AS cvd_informal_costs,

              -- Aggregated direct costs
              (chd_direct_costs + stroke_direct_costs) AS cvd_direct_costs,

              -- Aggregated indirect costs (productivity + informal)
              (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + chd_informal_costs) AS chd_indirect_costs,
              (stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs + stroke_informal_costs) AS stroke_indirect_costs,
              (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + chd_informal_costs + stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs + stroke_informal_costs) AS cvd_indirect_costs,
            
              -- Total costs (indirect + direct)
              (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + chd_informal_costs + chd_direct_costs) AS chd_total_costs,
              (stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs + stroke_informal_costs + stroke_direct_costs) AS stroke_total_costs,
              (chd_prvl_prdv_costs + chd_mrtl_prdv_costs + chd_informal_costs + chd_direct_costs + stroke_prvl_prdv_costs + stroke_mrtl_prdv_costs + stroke_informal_costs + stroke_direct_costs) AS cvd_total_costs
            
            FROM basic_costs;
        ",
          paste0(output_view_name, "_", scnams, "_view"),
          input_table_name,
          mcaggr,
          paste0("'", scnams, "'")
        )

        sapply(final_view_creation_sql, function(sql) {
          private$execute_sql(
            duckdb_con,
            sql,
            paste0("Final cost views per scenario:", output_view_name)
          )
        })

        return(invisible(NULL))
      }, # end calc_costs

      # collect_files ----
      # Collect files written by mc_aggr or mc_aggr_mc in a folder and combine
      # them into one file
      # @param folder_name A character string specifying the name of the folder containing the files to be collected.
      # @param pattern A regular expression pattern to filter files for collection. Only files matching this pattern will be considered.
      # @param to_mc_aggr Logical, indicating whether to aggregate files for Monte Carlo iteration (to_mc_aggr = TRUE) or not.
      #
      # @details
      # The function reads each file from the specified folder, aggregates its content, and, if needed, compresses the aggregated result. The aggregation process involves rewriting the files by applying certain string patterns to the filenames.
      #
      # @return
      # This function does not return a value directly. It collects and aggregates files, and compresses them if specified, cleaning up the original files in the process.
      #
      collect_files = function(
        folder_name,
        pattern = NULL,
        to_mc_aggr = FALSE
      ) {
        if (self$design$sim_prm$logs) {
          message("Collecting mc files...")
        }
        if (to_mc_aggr) {
          string1 <- "_[0-9]+_"
          string2 <- "_"
        } else {
          string1 <- "[0-9]+_"
          string2 <- ""
        }
        sapply(
          list.files(
            path = private$output_dir(folder_name),
            pattern = pattern,
            full.names = TRUE
          ),
          function(fnam) {
            fwrite_safe(fread(fnam), file = sub(string1, string2, fnam))
            file.remove(fnam)
          }
        )
        # gzip the .csv files to .csv.gz (faster than using gzip() and same speed/compression as with fst 80. But fst reads faster)
        if (self$design$sim_prm$logs) {
          message("Compressing aggregated files...")
        }
        sapply(
          list.files(
            path = private$output_dir(folder_name),
            pattern = sub("^_", "", pattern),
            full.names = TRUE
          ),
          function(fnam) {
            fwrite_safe(fread(fnam), file = gsub(".csv$", ".csv.gz", fnam))
            file.remove(fnam)
          }
        )
        NULL
      },

      # create_new_folder ----
      # @details
      # The function checks if the specified folder exists. If it does not exist, it creates the directory using the \code{dir.create} function. The \code{recursive} argument is set to TRUE to create parent directories if needed.
      #
      # @return
      # This function does not return a value directly. It creates the specified folder and reports the status if specified.
      #
      # @description Create folder if doesn't exist. Stops on failure.
      # @param sDirPathName String folder path and name.
      # @param bReport Bool report folder creation.
      create_new_folder = function(sDirPathName, bReport) {
        if (!dir.exists(sDirPathName)) {
          bSuccess <- dir.create(sDirPathName, recursive = TRUE)
          if (!bSuccess) {
            stop(paste("Failed creating directory", sDirPathName))
          }
          if (bReport) message(paste0("Folder ", sDirPathName, " was created"))
        }
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
        is_docker_env <- file.exists("/.dockerenv")
        is_docker_windows <- (is_docker_env &&
          (.Platform$OS.type == "unix" ||
            Sys.getenv("DOCKER_DESKTOP") != "" ||
            any(grepl("/host_mnt/", getwd(), fixed = TRUE))))

        # Use safer method for Docker environments (especially Windows Docker Desktop)
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
                result_data <- private$query_sql(duckdb_con, query)
                if (self$design$sim_prm$logs) {
                  cat(sprintf(
                    "Query returned %d rows for %s\n",
                    nrow(result_data),
                    basename(output_path)
                  ))
                }

                if (nrow(result_data) > 0) {
                  arrow::write_parquet(result_data, output_path)
                  Sys.sleep(0.3)

                  if (file.exists(output_path)) {
                    file_size <- file.size(output_path)
                    if (file_size > 0) {
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
                    }
                  }
                } else {
                  if (self$design$sim_prm$logs) {
                    cat(sprintf(
                      "Query returned 0 rows for %s - creating empty file\n",
                      basename(output_path)
                    ))
                  }
                  empty_structure <- private$query_sql(
                    duckdb_con,
                    sprintf("SELECT * FROM (%s) subq LIMIT 0", query)
                  )
                  arrow::write_parquet(empty_structure, output_path)
                  if (file.exists(output_path)) {
                    success <- TRUE
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
                if (file.exists(output_path)) {
                  try(file.remove(output_path), silent = TRUE)
                }
                Sys.sleep(0.5 * retry_count)
              }
            )
          }
        } else {
          # Use original DuckDB COPY method for native Linux environments
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
                db_path <- gsub("\\\\", "/", output_path)
                copy_command <- sprintf(
                  "COPY (%s) TO '%s' (FORMAT PARQUET);",
                  query,
                  db_path
                )
                result <- private$execute_sql(duckdb_con, copy_command)
                Sys.sleep(0.2)

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
                  try(file.remove(output_path), silent = TRUE)
                }
                Sys.sleep(0.2 * retry_count)
              }
            )
          }
        }

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

      # Life Expectancy (LE) Export Section
      #
      # This block exports life expectancy (LE) at minAge and life expectancy at age 60 (LE60)
      # summaries for each Monte Carlo (mc) iteration. It uses a configuration list
      # and a loop to avoid code repetition. Each configuration specifies:
      #   - prefix: output file prefix ("le" or "le60")
      #   - age_filter: SQL WHERE clause for age (e.g., only ages > 60 for LE60)
      #   - weight_col: which weight column to use ("wt" or "wt_esp")
      #   - suffix: output file suffix (e.g., "_scaled_up" or "_esp")
      #
      # The code dynamically constructs the SQL query and output path for each
      # configuration, then executes the query and writes the result to a Parquet file.
      #
      # Directories for all output types are created if they do not exist.
      # The LE60 calculation is skipped if the simulation age range does not cover age 60.
      #
      # Args:
      #   type: character vector, must include "le" to trigger this block
      #   duckdb_con: DuckDB connection to the lifecourse data
      #   mc: Monte Carlo iteration number
      #   strata_noagegrp: grouping columns for summary
      #   ext: file extension (usually "parquet")
      #
      # Output:
      #   Writes LE and LE60 summary files to the appropriate summaries/ subfolders.
      #
      # Example output files:
      #   summaries/le_scaled_up/1_le_scaled_up.parquet
      #   summaries/le_esp/1_le_esp.parquet
      #   summaries/le60_scaled_up/1_le60_scaled_up.parquet
      #   summaries/le60_esp/1_le60_esp.parquet
      #
      # --- End documentation ---
      # export_le_summaries ----
      export_le_summaries = function(duckdb_con, mcaggr, strata_noagegrp, ext) {
        lapply(
          paste0(rep(c("le", "le60"), each = 2), "_", c("scaled_up", "esp")),
          function(subdir) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir
            )))
          }
        )

        group_by_cols <- paste(strata_noagegrp, collapse = ", ")

        # Define configurations for LE calculations
        le_configs <- list(
          list(
            prefix = "le",
            age_filter = "",
            weight_col = "wt",
            suffix = "_scaled_up"
          ),
          list(
            prefix = "le",
            age_filter = "",
            weight_col = "wt_esp",
            suffix = "_esp"
          ),
          list(
            prefix = "le60",
            age_filter = "AND age > 60",
            weight_col = "wt",
            suffix = "_scaled_up"
          ),
          list(
            prefix = "le60",
            age_filter = "AND age > 60",
            weight_col = "wt_esp",
            suffix = "_esp"
          )
        )

        for (config in le_configs) {
          # Skip LE60 calculation if age range doesn't cover 60
          if (
            startsWith(config$prefix, "le60") &&
              !(self$design$sim_prm$ageL < 60L &&
                self$design$sim_prm$ageH > 60L)
          ) {
            next
          }

          query <- sprintf(
            "SELECT %s, SUM(%s) AS popsize, SUM(age * %s) / NULLIF(SUM(%s), 0) AS LE
            FROM lc_table
            WHERE mc = %d AND all_cause_mrtl > 0 %s
            GROUP BY %s
            ORDER BY %s",
            group_by_cols,
            config$weight_col,
            config$weight_col,
            config$weight_col, # Use weight_col 3 times
            mcaggr,
            config$age_filter,
            group_by_cols,
            group_by_cols
          )

          output_path <- private$output_dir(
            paste0(
              "summaries/",
              config$prefix,
              config$suffix,
              "/",
              mcaggr,
              "_",
              config$prefix,
              config$suffix,
              ".",
              ext
            )
          )

          # Execute query and write result with retry logic
          private$execute_db_diskwrite_with_retry(
            duckdb_con,
            query,
            output_path
          )
        }
        NULL
      }, # end export_le_summaries

      # Healthy Life Expectancy (HLE) Export Section
      #
      # This section exports healthy life expectancy (HLE) summaries for each Monte Carlo (mc)
      # iteration, using different health state definitions and weighting schemes.
      # It uses a configuration list and a loop to avoid code repetition.
      #
      # For each configuration, the following parameters are specified:
      #   - prefix: output file prefix (e.g., "hle_1st_cond", "hle_cmsmm1.5")
      #   - condition: SQL WHERE clause defining the 'healthy' state (e.g., "cms_count = 1")
      #   - weight_col: which weight column to use ("wt" or "wt_esp")
      #   - suffix: output file suffix (e.g., "_scaled_up" or "_esp")
      #
      # The code dynamically constructs the SQL query and output path for each configuration,
      # executes the query using DuckDB, and writes the result to a Parquet file.
      # Output directories are created if they do not exist.
      #
      # Args:
      #   type: character vector, must include "hle" to trigger this block
      #   duckdb_con: DuckDB connection to the lifecourse data
      #   mc: Monte Carlo iteration number
      #   strata_noagegrp: grouping columns for summary
      #   ext: file extension (usually "parquet")
      #
      # Output:
      #   Writes HLE summary files to the appropriate summaries/ subfolders.
      #
      # Example output files:
      #   summaries/hle_1st_cond_scaled_up/1_hle_1st_cond_scaled_up.parquet
      #   summaries/hle_1st_cond_esp/1_hle_1st_cond_esp.parquet
      #   summaries/hle_cmsmm1.5_scaled_up/1_hle_cmsmm1.5_scaled_up.parquet
      #   summaries/hle_cmsmm1.5_esp/1_hle_cmsmm1.5_esp.parquet
      #
      # --- End documentation ---
      # export_hle_summaries ----
      export_hle_summaries = function(
        duckdb_con,
        mcaggr,
        strata_noagegrp,
        ext
      ) {
        # TODO currently some individuals are counted more than once because
        # disease counter and score can be reduced.
        # Ideally only the first reach to the threshold should be counted
        lapply(
          paste0(
            rep(c("hle_1st_cond", "hle_cmsmm1.5"), each = 2),
            "_",
            c("scaled_up", "esp")
          ),
          function(subdir) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir
            )))
          }
        )

        group_by_cols <- paste(strata_noagegrp, collapse = ", ")

        # Define configurations for HLE calculations
        hle_configs <- list(
          list(
            prefix = "hle_1st_cond",
            condition = "cms_count = 1",
            weight_col = "wt",
            suffix = "_scaled_up"
          ),
          list(
            prefix = "hle_1st_cond",
            condition = "cms_count = 1",
            weight_col = "wt_esp",
            suffix = "_esp"
          ),
          list(
            prefix = "hle_cmsmm1.5",
            condition = '"cmsmm1.5_prvl" = 1',
            weight_col = "wt",
            suffix = "_scaled_up"
          ),
          list(
            prefix = "hle_cmsmm1.5",
            condition = '"cmsmm1.5_prvl" = 1',
            weight_col = "wt_esp",
            suffix = "_esp"
          )
        )

        for (config in hle_configs) {
          query <- sprintf(
            "SELECT %s, SUM(%s) AS popsize, SUM(age * %s) / NULLIF(SUM(%s), 0) AS HLE
               FROM lc_table
               WHERE mc = %d AND %s
               GROUP BY %s
               ORDER BY %s",
            group_by_cols,
            config$weight_col,
            config$weight_col,
            config$weight_col,
            mcaggr,
            config$condition,
            group_by_cols,
            group_by_cols
          )

          output_path <- private$output_dir(
            paste0(
              "summaries/",
              config$prefix,
              config$suffix,
              "/",
              mcaggr,
              "_",
              config$prefix,
              config$suffix,
              ".",
              ext
            )
          )
          private$execute_db_diskwrite_with_retry(
            duckdb_con,
            query,
            output_path
          )

          NULL
        }
      }, # end export_hle_summaries

      # The code block is responsible for exporting disease
      # characteristics summaries for each Monte Carlo (mc) iteration. It
      # dynamically queries the DuckDB lifecourse table to calculate, for each
      # disease (identified by columns ending in _prvl), the number of cases,
      # mean age at incidence, mean age at first onset, mean age at
      # prevalence, mean duration, mean CMS score, and mean CMS count, grouped
      # by the specified strata. The results are written to Parquet files in
      # the summaries/dis_characteristics_scaled_up and
      # summaries/dis_characteristics_esp directories, using both standard and
      # ESP weights. The code uses SQL PIVOT to reshape the output so that
      # each disease's metrics become columns, and the output is suitable for
      # further analysis or reporting. This process is repeated for both
      # standard and ESP-weighted results.
      # export_dis_char_summaries ----
      export_dis_char_summaries = function(
        duckdb_con,
        mcaggr,
        strata_noagegrp,
        ext
      ) {
        # Create output directories
        lapply(
          paste0(
            rep(c("dis_characteristics"), each = 2),
            "_",
            c("scaled_up", "esp")
          ),
          function(subdir) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir
            )))
          }
        )

        # Get disease prevalence columns from DuckDB schema
        lc_table_name <- "lc_table"
        all_cols <- dbListFields(duckdb_con, lc_table_name)
        nm <- grep("_prvl$", all_cols, value = TRUE)

        # Construct the SQL query dynamically
        # Ensure agegrp is excluded from grouping to aggregate over it
        select_cols_no_agegrp <- paste(
          setdiff(strata_noagegrp, "agegrp"),
          collapse = ", "
        )

        if (length(nm) > 0) {
          # Memory-efficient approach: Process diseases one at a time
          disease_results_scaled_up <- list()
          disease_results_esp <- list()

          # Get disease names without _prvl suffix
          disease_names <- gsub("_prvl$", "", nm)

          # Process each disease individually to avoid memory issues
          for (i in seq_along(nm)) {
            disease_col <- nm[i]
            disease_name <- disease_names[i]
            quoted_disease_col <- paste0('"', disease_col, '"')

            # Query for current disease (scaled_up version)
            disease_query_scaled_up <- sprintf(
              "SELECT %s, '%s' AS disease,
               SUM(wt) AS cases,
               SUM(CASE WHEN %s = 1 THEN age * wt ELSE 0 END) / NULLIF(SUM(CASE WHEN %s = 1 THEN wt ELSE 0 END), 0) AS mean_age_incd,
               SUM(age * wt) / NULLIF(SUM(wt), 0) AS mean_age_prvl,
               SUM(%s * wt) / NULLIF(SUM(wt), 0) AS mean_duration,
               SUM(cms_score * wt) / NULLIF(SUM(wt), 0) AS mean_cms_score,
               SUM(cms_count * wt) / NULLIF(SUM(wt), 0) AS mean_cms_count
               FROM %s
               WHERE mc = %d AND %s > 0
               GROUP BY %s",
              select_cols_no_agegrp,
              disease_name,
              quoted_disease_col,
              quoted_disease_col,
              quoted_disease_col,
              lc_table_name,
              mcaggr,
              quoted_disease_col,
              select_cols_no_agegrp
            )

            # Execute and store result for scaled_up
            result_scaled_up <- private$query_sql(
              duckdb_con,
              disease_query_scaled_up,
              paste("Disease characteristics for", disease_name, "(scaled_up)")
            )
            if (nrow(result_scaled_up) > 0) {
              disease_results_scaled_up[[i]] <- result_scaled_up
            }

            # Query for current disease (ESP version)
            disease_query_esp <- gsub("wt", "wt_esp", disease_query_scaled_up)

            # Execute and store result for ESP
            result_esp <- private$query_sql(
              duckdb_con,
              disease_query_esp,
              paste("Disease characteristics for", disease_name, "(ESP)")
            )
            if (nrow(result_esp) > 0) {
              disease_results_esp[[i]] <- result_esp
            }
          }

          # Combine and process scaled_up results
          if (length(disease_results_scaled_up) > 0) {
            # Remove NULL entries
            disease_results_scaled_up <- disease_results_scaled_up[
              !sapply(disease_results_scaled_up, is.null)
            ]

            if (length(disease_results_scaled_up) > 0) {
              combined_scaled_up <- rbindlist(disease_results_scaled_up)
              setDT(combined_scaled_up)

              # Get strata columns without agegrp for formula
              strata_no_agegrp <- setdiff(strata_noagegrp, "agegrp")

              # Create pivot transformation using data.table dcast
              pivot_result_scaled_up <- dcast(
                combined_scaled_up,
                formula = as.formula(paste(
                  paste(strata_no_agegrp, collapse = " + "),
                  "~ disease"
                )),
                value.var = c(
                  "cases",
                  "mean_age_incd",
                  "mean_age_prvl",
                  "mean_duration",
                  "mean_cms_score",
                  "mean_cms_count"
                ),
                fill = 0
              )

              # Order by strata columns (excluding agegrp)
              setkeyv(pivot_result_scaled_up, strata_no_agegrp)

              # Write scaled_up version
              output_path <- private$output_dir(
                paste0(
                  "summaries/dis_characteristics_scaled_up/",
                  mcaggr,
                  "_dis_characteristics_scaled_up.",
                  ext
                )
              )
              arrow::write_parquet(pivot_result_scaled_up, output_path)
            }
          }

          # Combine and process ESP results
          if (length(disease_results_esp) > 0) {
            # Remove NULL entries
            disease_results_esp <- disease_results_esp[
              !sapply(disease_results_esp, is.null)
            ]

            if (length(disease_results_esp) > 0) {
              combined_esp <- rbindlist(disease_results_esp)
              setDT(combined_esp)

              # Get strata columns without agegrp for formula
              strata_no_agegrp <- setdiff(strata_noagegrp, "agegrp")

              # Create pivot transformation using data.table dcast
              pivot_result_esp <- dcast(
                combined_esp,
                formula = as.formula(paste(
                  paste(strata_no_agegrp, collapse = " + "),
                  "~ disease"
                )),
                value.var = c(
                  "cases",
                  "mean_age_incd",
                  "mean_age_prvl",
                  "mean_duration",
                  "mean_cms_score",
                  "mean_cms_count"
                ),
                fill = 0
              )

              # Order by strata columns (excluding agegrp)
              setkeyv(pivot_result_esp, strata_no_agegrp)

              # Write ESP version
              output_path_esp <- private$output_dir(
                paste0(
                  "summaries/dis_characteristics_esp/",
                  mcaggr,
                  "_dis_characteristics_esp.",
                  ext
                )
              )
              arrow::write_parquet(pivot_result_esp, output_path_esp)
            }
          }
        }

        NULL
      }, # end export_dis_char_summaries

      # export_prvl_summaries ----
      export_prvl_summaries = function(duckdb_con, mcaggr, strata, ext) {
        # TODO currently some individuals are counted more than once because
        # disease counter and score can be reduced.
        # Ideally only the first reach to the threshold should be counted
        lapply(
          paste0(rep(c("prvl"), each = 2), "_", c("scaled_up", "esp")),
          function(subdir) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir
            )))
          }
        )

        # Get disease prevalence columns from DuckDB schema
        lc_table_name <- "lc_table" # Assuming the view/table name in DuckDB is lc_table
        all_cols <- dbListFields(duckdb_con, lc_table_name)
        nm_prvl <- grep("_prvl$", all_cols, value = TRUE)

        # Construct the SQL query dynamically
        select_cols <- paste(strata, collapse = ", ")
        sum_cases_cols <- paste(
          sprintf(
            'SUM(CASE WHEN "%s" > 0 THEN wt ELSE 0 END) AS "%s"',
            nm_prvl,
            nm_prvl
          ),
          collapse = ", "
        )

        sql_query <- sprintf(
          "SELECT %s, SUM(wt) AS popsize, %s
                  FROM %s
                  WHERE mc = %d
                  GROUP BY %s
                  ORDER BY %s",
          select_cols,
          sum_cases_cols,
          lc_table_name,
          mcaggr,
          select_cols,
          select_cols
        )

        output_path <- private$output_dir(
          paste0("summaries/prvl_scaled_up/", mcaggr, "_prvl_scale_up.", ext)
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query,
          output_path
        )

        # esp
        sql_query_esp <- gsub("wt", "wt_esp", sql_query)
        output_path_esp <- private$output_dir(
          paste0("summaries/prvl_esp/", mcaggr, "_prvl_esp.", ext)
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query_esp,
          output_path_esp
        )

        NULL
      }, # end of export_prvl_summaries

      # export_incd_summaries ----
      export_incd_summaries = function(duckdb_con, mcaggr, strata, ext) {
        # TODO currently some individuals are counted more than once because
        # disease counter and score can be reduced.
        # Ideally only the first reach to the threshold should be counted
        lapply(
          paste0(rep(c("incd"), each = 2), "_", c("scaled_up", "esp")),
          function(subdir) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir
            )))
          }
        )

        # Get disease prevalence columns from DuckDB schema
        lc_table_name <- "lc_table" # Assuming the view/table name in DuckDB is lc_table
        all_cols <- dbListFields(duckdb_con, lc_table_name)
        nm_prvl <- grep("_prvl$", all_cols, value = TRUE)

        # Construct the SQL query dynamically
        select_cols <- paste(strata, collapse = ", ")
        sum_cases_cols <- paste(
          sprintf(
            'SUM(CASE WHEN "%s" = 1 THEN wt ELSE 0 END) AS "%s"',
            nm_prvl,
            gsub("_prvl$", "_incd", nm_prvl)
          ),
          collapse = ", "
        )

        sql_query <- sprintf(
          "SELECT %s, SUM(wt) AS popsize, %s
                  FROM %s
                  WHERE mc = %d
                  GROUP BY %s
                  ORDER BY %s",
          select_cols,
          sum_cases_cols,
          lc_table_name,
          mcaggr,
          select_cols,
          select_cols
        )

        output_path <- private$output_dir(
          paste0("summaries/incd_scaled_up/", mcaggr, "_incd_scale_up.", ext)
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query,
          output_path
        )

        # esp
        sql_query_esp <- gsub("wt", "wt_esp", sql_query)
        output_path_esp <- private$output_dir(
          paste0("summaries/incd_esp/", mcaggr, "_incd_esp.", ext)
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query_esp,
          output_path_esp
        )

        NULL
      }, # end of export_incd_summaries

      # export_mrtl_summaries ----
      export_mrtl_summaries = function(duckdb_con, mcaggr, strata, ext) {
        lapply(
          paste0(rep(c("mrtl"), each = 2), "_", c("scaled_up", "esp")),
          function(subdir) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir
            )))
          }
        )

        # Get disease prevalence columns from DuckDB schema
        lc_table_name <- "lc_table" # Assuming the view/table name in DuckDB is lc_table
        # Construct the SQL query dynamically for mortality summary
        select_cols_mrtl <- paste(strata, collapse = ", ")
        sql_query <- sprintf(
          "SELECT %s,
                    SUM(wt) AS popsize,
                    SUM(CASE WHEN all_cause_mrtl > 0 THEN wt ELSE 0 END) AS all_cause_mrtl
             FROM %s
             WHERE mc = %d
             GROUP BY %s
             ORDER BY %s", # Add ORDER BY to match data.table's keyby behavior
          select_cols_mrtl,
          lc_table_name,
          mcaggr,
          select_cols_mrtl, # Quoted strata cols for GROUP BY
          select_cols_mrtl # Order by the grouping columns
        )
        output_path <- private$output_dir(
          paste0("summaries/mrtl_scaled_up/", mcaggr, "_mrtl_scale_up.", ext)
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query,
          output_path
        )

        # esp
        sql_query_esp <- gsub("wt", "wt_esp", sql_query)
        output_path_esp <- private$output_dir(
          paste0("summaries/mrtl_esp/", mcaggr, "_mrtl_esp.", ext)
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query_esp,
          output_path_esp
        )

        NULL
      }, # end of export_mrtl_summaries

      # export_dis_mrtl_summaries ----
      # this exports cause specific mortality summaries for each Monte Carlo (mc)
      # iteration. It dynamically queries the DuckDB lifecourse table.
      export_dis_mrtl_summaries = function(duckdb_con, mcaggr, strata, ext) {
        lapply(
          paste0(rep(c("dis_mrtl"), each = 2), "_", c("scaled_up", "esp")),
          function(subdir) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir
            )))
          }
        )

        # Get disease prevalence columns from DuckDB schema
        lc_table_name <- "lc_table" # Assuming the view/table name in DuckDB is lc_table
        # Construct the SQL query dynamically for mortality summary
        # Define strata and death codes for the query
        quoted_strata <- paste0('"', strata, '"')
        strata_cols_sql <- paste(quoted_strata, collapse = ", ")

        # All pivoted column names, including potentially 'alive_deaths'
        pivoted_death_col_names_all <- paste0(
          names(private$death_codes),
          "_deaths"
        )
        quoted_pivoted_death_col_names_all <- paste0(
          '"',
          pivoted_death_col_names_all,
          '"'
        )

        # Correctly format the IN clause with quoted aliases for PIVOT (includes 'alive_deaths')
        death_codes_pivot_sql <- paste0(
          "'",
          private$death_codes,
          "' AS ", # Use single quotes for values, double for aliases
          quoted_pivoted_death_col_names_all, # Use quoted aliases
          collapse = ", "
        )

        # Pivoted column names EXCLUDING 'alive_deaths' for the final SELECT list
        pivoted_death_col_names_final <- setdiff(
          pivoted_death_col_names_all,
          "alive_deaths"
        )
        quoted_pivoted_death_col_names_final <- paste0(
          '"',
          pivoted_death_col_names_final,
          '"'
        )

        # Create COALESCE expressions for the final SELECT list (excludes 'alive_deaths')
        # Use quoted column names and aliases
        coalesce_select_sql_final <- paste0(
          "COALESCE(t.",
          quoted_pivoted_death_col_names_final,
          ", 0) AS ",
          quoted_pivoted_death_col_names_final, # Use quoted aliases
          collapse = ", "
        )

        # Sum of ALL coalesced columns for popsize (includes 'alive_deaths')
        # Note: This uses quoted_pivoted_death_col_names_all
        death_cols_sum_sql <- paste0(
          "COALESCE(t.",
          quoted_pivoted_death_col_names_all,
          ", 0)",
          collapse = " + "
        )

        # Select strata columns prefixed with t. and quoted
        select_strata_sql <- paste0("t.", quoted_strata, collapse = ", ")

        # Construct the PIVOT SQL query (inner query) - uses all death codes
        sql_query <- sprintf(
          "WITH AggregatedDeaths AS (
           SELECT
             %s,
             all_cause_mrtl,
             SUM(wt) AS deaths
           FROM %s
           WHERE mc = %d -- Filter by mc
           GROUP BY %s, all_cause_mrtl
           )
           PIVOT AggregatedDeaths
           ON all_cause_mrtl IN (%s) -- Specify codes and quoted aliases to pivot
           USING SUM(deaths) -- Aggregation function
           GROUP BY %s -- Columns to keep (quoted)
           ", # Removed trailing semicolon
          strata_cols_sql, # Quoted strata cols for SELECT
          lc_table_name,
          mcaggr,
          strata_cols_sql, # Quoted strata cols for GROUP BY
          death_codes_pivot_sql, # Use pivot definition including quoted aliases
          strata_cols_sql # Quoted strata cols for GROUP BY
        )

        # Construct the final SQL query with COALESCE (excluding alive_deaths) and popsize calculation (including alive_deaths)
        sql_query_final <- sprintf(
          "SELECT
           %s, -- Select quoted strata columns
           %s, -- Select coalesced death columns (excluding alive_deaths, quoted)
           %s AS popsize -- Calculate popsize from ALL coalesced columns (including alive_deaths, quoted)
           FROM (%s) AS t
           ORDER BY %s", # Order by quoted strata columns, no semicolon
          select_strata_sql,
          coalesce_select_sql_final,
          death_cols_sum_sql,
          sql_query,
          strata_cols_sql
        ) # Use quoted strata_cols_sql for ORDER BY

        output_path <- private$output_dir(
          paste0(
            "summaries/dis_mrtl_scaled_up/",
            mcaggr,
            "_dis_mrtl_scale_up.",
            ext
          )
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query_final,
          output_path
        )

        # esp
        sql_query_esp <- gsub("wt", "wt_esp", sql_query_final)
        output_path_esp <- private$output_dir(
          paste0("summaries/dis_mrtl_esp/", mcaggr, "_dis_mrtl_esp.", ext)
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query_esp,
          output_path_esp
        )

        NULL
      }, # end of export_dis_mrtl_summaries

      # export_all_cause_mrtl_by_dis_summaries ----
      export_all_cause_mrtl_by_dis_summaries = function(
        duckdb_con,
        mcaggr,
        strata,
        ext
      ) {
        lapply(
          paste0(
            rep(c("all_cause_mrtl_by_dis"), each = 2),
            "_",
            c("scaled_up", "esp")
          ),
          function(subdir) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir
            )))
          }
        )

        # Get disease prevalence columns from DuckDB schema
        lc_table_name <- "lc_table"
        schema <- private$query_sql(
          duckdb_con,
          sprintf("DESCRIBE %s;", lc_table_name)
        )
        prvl_cols <- schema$column_name[endsWith(schema$column_name, "_prvl")]
        disease_names <- gsub("_prvl$", "", prvl_cols)

        # Quote strata columns for safety
        quoted_strata <- paste0('"', strata, '"')
        strata_cols_sql <- paste(quoted_strata, collapse = ", ")

        # Construct CASE WHEN statements for each disease's cases and deaths
        case_statements <- sapply(disease_names, function(dis) {
          prvl_col <- paste0('"', dis, '_prvl"')
          sprintf(
            'SUM(CASE WHEN %s > 0 THEN wt ELSE 0 END) AS "cases_%s"',
            prvl_col,
            dis
          )
        })

        death_statements <- sapply(disease_names, function(dis) {
          prvl_col <- paste0('"', dis, '_prvl"')
          sprintf(
            'SUM(CASE WHEN %s > 0 AND all_cause_mrtl > 0 THEN wt ELSE 0 END) AS "deaths_%s"',
            prvl_col,
            dis
          )
        })

        # Combine all select parts
        select_parts_sql <- paste(
          c(strata_cols_sql, case_statements, death_statements),
          collapse = ",\n  "
        )

        # Construct the full SQL query
        sql_query_dis_char <- sprintf(
          "
        SELECT
          %s
        FROM %s
        WHERE mc = %d -- Filter for the specific mc iteration if needed
        GROUP BY %s
        ORDER BY %s
        ",
          select_parts_sql,
          lc_table_name,
          mcaggr,
          strata_cols_sql,
          strata_cols_sql
        )

        output_path <- private$output_dir(
          paste0(
            "summaries/all_cause_mrtl_by_dis_scaled_up/",
            mcaggr,
            "_all_cause_mrtl_by_dis_scale_up.",
            ext
          )
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query_dis_char,
          output_path
        )

        # esp
        sql_query_esp <- gsub("wt", "wt_esp", sql_query_dis_char)
        output_path_esp <- private$output_dir(
          paste0(
            "summaries/all_cause_mrtl_by_dis_esp/",
            mcaggr,
            "_all_cause_mrtl_by_dis_esp.",
            ext
          )
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          sql_query_esp,
          output_path_esp
        )

        NULL
      }, # end of export_all_cause_mrtl_by_dis_summaries

      # export_cms_summaries ----
      export_cms_summaries = function(
        duckdb_con,
        mcaggr,
        strata,
        strata_age,
        ext
      ) {
        lapply(
          paste0(
            rep(c("cms_score", "cms_score_by_age", "cms_count"), each = 2),
            "_",
            c("scaled_up", "esp")
          ),
          function(subdir) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir
            )))
          }
        )

        # Define configurations for CMS calculations
        # strata_sql is strata; strata_age_sql is strata_age
        strata_sql <- paste(strata, collapse = ", ")
        strata_age_sql <- paste(strata_age, collapse = ", ")

        cms_configs <- list(
          list(
            metric = "cms_score",
            group_cols_sql = strata_sql,
            group_cols_r = strata,
            weight_col = "wt",
            suffix = "_scaled_up",
            file_group_suffix = ""
          ),
          list(
            metric = "cms_score",
            group_cols_sql = strata_age_sql,
            group_cols_r = strata_age,
            weight_col = "wt",
            suffix = "_by_age_scaled_up",
            file_group_suffix = "_by_age"
          ),
          list(
            metric = "cms_score",
            group_cols_sql = strata_sql,
            group_cols_r = strata,
            weight_col = "wt_esp",
            suffix = "_esp",
            file_group_suffix = ""
          ),
          list(
            metric = "cms_score",
            group_cols_sql = strata_age_sql,
            group_cols_r = strata_age,
            weight_col = "wt_esp",
            suffix = "_by_age_esp",
            file_group_suffix = "_by_age"
          ),
          list(
            metric = "cms_count",
            group_cols_sql = strata_sql,
            group_cols_r = strata,
            weight_col = "wt",
            suffix = "_scaled_up",
            file_group_suffix = ""
          ),
          list(
            metric = "cms_count",
            group_cols_sql = strata_sql,
            group_cols_r = strata,
            weight_col = "wt_esp",
            suffix = "_esp",
            file_group_suffix = ""
          )
        )

        lc_table_name <- "lc_table"

        for (config in cms_configs) {
          query <- sprintf(
            "SELECT %s, SUM(%s) AS popsize, SUM(%s * %s) / SUM(%s) AS %s
             FROM %s
             WHERE mc = %d
             GROUP BY %s
             ORDER BY %s",
            config$group_cols_sql,
            config$weight_col,
            config$metric,
            config$weight_col,
            config$weight_col,
            config$metric,
            lc_table_name,
            mcaggr,
            config$group_cols_sql,
            config$group_cols_sql
          )

          output_path <- private$output_dir(
            paste0(
              "summaries/",
              config$metric,
              config$file_group_suffix,
              ifelse(config$weight_col == "wt_esp", "_esp", "_scaled_up"),
              "/",
              mcaggr,
              "_",
              config$metric,
              config$file_group_suffix,
              ifelse(config$weight_col == "wt_esp", "_esp", "_scaled_up"),
              ".",
              ext
            )
          )

          # Ensure output directory for the specific file exists (handles cases like cms_count_scaled_up)
          private$create_new_folder(dirname(output_path))

          private$execute_db_diskwrite_with_retry(
            duckdb_con,
            query,
            output_path
          )
        }

        NULL
      }, # end of export_cms_summaries

      # export_qalys_summaries ----
      export_qalys_summaries = function(duckdb_con, mcaggr, strata, ext) {
        lapply(
          paste0(rep(c("qalys"), each = 2), "_", c("scaled_up", "esp")),
          function(subdir_suffix) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir_suffix
            )))
          }
        )

        lc_table_name <- "lc_table"

        # Define the name for the temporary view that calc_QALYs will create
        qaly_view_name <- "lc_with_qalys_view"

        # Call calc_QALYs to create/replace the temporary view with EQ5D5L and HUI3 columns.
        # This view will be based on lc_table and filtered for the current mcaggr.
        private$calc_QALYs(
          duckdb_con = duckdb_con,
          mcaggr = mcaggr,
          input_table_name = lc_table_name,
          output_view_name = qaly_view_name,
          include_non_significant = FALSE
        )

        # Prepare strata columns for SQL query (quoted)
        # 'strata' is defined earlier in export_summaries_hlpr
        quoted_strata_cols_sql <- paste(
          sprintf('"%s"', strata),
          collapse = ", "
        )

        # Define QALY metrics for SELECT statement
        qaly_metrics_select_wt <- 'SUM("EQ5D5L" * wt) AS "EQ5D5L", SUM("HUI3" * wt) AS "HUI3"'
        qaly_metrics_select_wt_esp <- 'SUM("EQ5D5L" * wt_esp) AS "EQ5D5L", SUM("HUI3" * wt_esp) AS "HUI3"'

        # --- Scaled-up QALYs ---
        query_scaled_up <- sprintf(
          "SELECT %s, SUM(wt) AS popsize, %s
             FROM %s -- This view is already filtered by mcaggr by calc_QALYs
             GROUP BY %s
             ORDER BY %s",
          quoted_strata_cols_sql,
          qaly_metrics_select_wt,
          qaly_view_name,
          quoted_strata_cols_sql,
          quoted_strata_cols_sql
        )
        output_path_scaled_up <- private$output_dir(
          paste0("summaries/qalys_scaled_up/", mcaggr, "_qalys_scaled_up.", ext)
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          query_scaled_up,
          output_path_scaled_up
        )

        # --- ESP QALYs ---
        query_esp <- sprintf(
          "SELECT %s, SUM(wt_esp) AS popsize, %s
             FROM %s -- This view is already filtered by mcaggr by calc_QALYs
             GROUP BY %s
             ORDER BY %s",
          quoted_strata_cols_sql,
          qaly_metrics_select_wt_esp,
          qaly_view_name,
          quoted_strata_cols_sql,
          quoted_strata_cols_sql
        )
        output_path_esp <- private$output_dir(
          paste0("summaries/qalys_esp/", mcaggr, "_qalys_esp.", ext)
        )

        private$execute_db_diskwrite_with_retry(
          duckdb_con,
          query_esp,
          output_path_esp
        )

        # drop the temporary view if it's no longer needed for this mcaggr
        private$execute_sql(
          duckdb_con,
          sprintf("DROP VIEW IF EXISTS %s;", qaly_view_name)
        )

        NULL
      }, # end of export_qalys_summaries

      # export_costs_summaries ----
      export_costs_summaries = function(duckdb_con, mcaggr, strata, ext) {
        # Create output directories for scaled-up and ESP-weighted summaries
        lapply(
          paste0(rep(c("costs"), each = 2), "_", c("scaled_up", "esp")),
          function(subdir_suffix) {
            private$create_new_folder(private$output_dir(paste0(
              "summaries/",
              subdir_suffix
            )))
          }
        )

        # get scenario names
        scnams <- gsub(
          "^scenario=",
          "",
          list.dirs(
            private$output_dir(file.path("lifecourse", paste0("mc=", mcaggr))),
            full.names = FALSE,
            recursive = FALSE
          )
        )

        lc_table_name <- "lc_table"

        # Define the name for the temporary view that calc_costs will create
        costs_view_name <- "lc_with_costs"

        costs_scn_views <- paste0(costs_view_name, "_", scnams, "_view")

        # Call calc_costs to create/replace the temporary view with cost columns.
        # This calculates cost parameters using sc0 baseline scenario data,
        # but applies them to all scenarios for the current mcaggr.
        private$calc_costs(
          duckdb_con = duckdb_con,
          mcaggr = mcaggr,
          input_table_name = lc_table_name,
          output_view_name = costs_view_name
        )

        # Prepare strata columns for SQL query (quoted)
        # The strata includes scenario as the first element
        quoted_strata_cols_sql <- paste(
          sprintf('"%s"', strata),
          collapse = ", "
        )

        # Define cost metrics for SELECT statement
        cost_metrics_select_wt_esp <- paste(
          'SUM(chd_direct_costs * wt_esp) AS chd_direct_costs',
          'SUM(stroke_direct_costs * wt_esp) AS stroke_direct_costs',
          'SUM(cvd_direct_costs * wt_esp) AS cvd_direct_costs',
          'SUM(chd_productivity_costs * wt_esp) AS chd_productivity_costs',
          'SUM(stroke_productivity_costs * wt_esp) AS stroke_productivity_costs',
          'SUM(cvd_productivity_costs * wt_esp) AS cvd_productivity_costs',
          'SUM(chd_informal_costs * wt_esp) AS chd_informal_costs',
          'SUM(stroke_informal_costs * wt_esp) AS stroke_informal_costs',
          'SUM(cvd_informal_costs * wt_esp) AS "cvd_informal_costs"',
          'SUM(chd_indirect_costs * wt_esp) AS chd_indirect_costs',
          'SUM(stroke_indirect_costs * wt_esp) AS stroke_indirect_costs',
          'SUM(cvd_indirect_costs * wt_esp) AS cvd_indirect_costs',
          'SUM(chd_total_costs * wt_esp) AS chd_total_costs',
          'SUM(stroke_total_costs * wt_esp) AS stroke_total_costs',
          'SUM(cvd_total_costs * wt_esp) AS cvd_total_costs',
          sep = ", "
        )

        # Memory-efficient approach: Process scenarios one at a time
        # Initialize result collectors for both ESP and scaled-up results
        esp_results <- list()
        scaled_up_results <- list()

        # Process each scenario individually to avoid loading all scenarios in RAM
        for (i in seq_along(scnams)) {
          scnam <- scnams[i]
          view_name <- costs_scn_views[i]

          # ESP query for current scenario
          query_esp_scenario <- sprintf(
            "SELECT %s, SUM(wt_esp) AS popsize, %s
             FROM %s
             GROUP BY %s",
            quoted_strata_cols_sql,
            cost_metrics_select_wt_esp,
            view_name,
            quoted_strata_cols_sql
          )

          # Execute and store result for ESP
          esp_results[[i]] <- as.data.table(private$query_sql(
            duckdb_con,
            query_esp_scenario,
            paste("ESP costs for scenario", scnam)
          ))

          # Scaled-up query for current scenario (replace wt_esp with wt)
          cost_metrics_select_wt <- gsub(
            "wt_esp",
            "wt",
            cost_metrics_select_wt_esp
          )
          query_scaled_up_scenario <- sprintf(
            "SELECT %s, SUM(wt) AS popsize, %s
             FROM %s
             GROUP BY %s",
            quoted_strata_cols_sql,
            cost_metrics_select_wt,
            view_name,
            quoted_strata_cols_sql
          )

          # Execute and store result for scaled-up
          scaled_up_results[[i]] <- as.data.table(private$query_sql(
            duckdb_con,
            query_scaled_up_scenario,
            paste("Scaled-up costs for scenario", scnam)
          ))
        }

        # Combine all ESP results and write to disk
        if (length(esp_results) > 0) {
          combined_esp <- rbindlist(esp_results, use.names = TRUE, fill = TRUE)
          rm(esp_results)
          # Order by strata columns
          setkeyv(combined_esp, strata)

          output_path_esp <- private$output_dir(
            paste0("summaries/costs_esp/", mcaggr, "_costs_esp.", ext)
          )

          arrow::write_parquet(combined_esp, output_path_esp)
        }

        # Combine all scaled-up results and write to disk
        if (length(scaled_up_results) > 0) {
          combined_scaled_up <- rbindlist(
            scaled_up_results,
            use.names = TRUE,
            fill = TRUE
          )
          rm(scaled_up_results)
          # Order by strata columns
          setkeyv(combined_scaled_up, strata)

          output_path_scaled_up <- private$output_dir(
            paste0(
              "summaries/costs_scaled_up/",
              mcaggr,
              "_costs_scaled_up.",
              ext
            )
          )

          arrow::write_parquet(combined_scaled_up, output_path_scaled_up)
        }

        # Drop the temporary views if they're no longer needed for this mcaggr
        sapply(costs_scn_views, function(view_name) {
          private$execute_sql(
            duckdb_con,
            sprintf("DROP VIEW IF EXISTS %s;", view_name)
          )
        })

        NULL
      }, # end of export_costs_summaries

      # Stub methods for export_tables - to be implemented
      # export_main_tables ----
      export_main_tables = function(prbl, baseline_year, output_dir) {
        # TODO: Implement main tables export
        if (self$design$sim_prm$logs) {
          message("export_main_tables: Not yet implemented")
        }
        NULL
      },

      # export_all_cause_mrtl_tables ----
      export_all_cause_mrtl_tables = function(prbl, summaries_dir, tables_dir) {
        # TODO: Implement all-cause mortality tables export
        if (self$design$sim_prm$logs) {
          message("export_all_cause_mrtl_tables: Not yet implemented")
        }
        NULL
      },

      # export_disease_characteristics_tables ----
      export_disease_characteristics_tables = function(prbl, summaries_dir, tables_dir) {
        # TODO: Implement disease characteristics tables export
        if (self$design$sim_prm$logs) {
          message("export_disease_characteristics_tables: Not yet implemented")
        }
        NULL
      },

      # export_xps_tables ----
      export_xps_tables = function(prbl, output_dir, tables_dir) {
        # TODO: Implement exposures tables export
        if (self$design$sim_prm$logs) {
          message("export_xps_tables: Not yet implemented")
        }
        NULL
      }
    )
  )

