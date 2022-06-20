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

        pth <- private$output_dir("summaries/")
        if (!dir.exists(pth)) {
          dir.create(pth)
          if (self$design$sim_prm$logs)
            message(paste0("Folder ", pth, " was created"))
        }

        pth <- private$output_dir("tables/")
        if (!dir.exists(pth)) {
          dir.create(pth)
          if (self$design$sim_prm$logs)
            message(paste0("Folder ", pth, " was created"))
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

        # European standardised population 2013 (esp) weights
        tt <- data.table(agegrp = agegrp_name(0, 99),
                         wt_esp  = c(1000, 4000, 5500, 5500, 5500, 6000, 6000, 6500,
                                     7000, 7000, 7000, 7000, 6500, 6000, 5500, 5000,
                                     4000, 2500, 1500, 800, 200))
        esp <- CJ(agegrp = agegrp_name(0, 99),
                  sex = c("men", "women"),
                  dimd = c("1 most deprived", as.character(2:9), "10 least deprived")
        )

        private$esp_weights <- copy(absorb_dt(esp, tt))

        # disnam <- paste0(names(self$diseases), "_prvl")
        # disnam <- disnam[disnam != "nonmodelled_prvl"]
        # disnam <-
        #   c(
        #     grep(
        #       "_ca_prvl$|^ctd_prvl$|^ra_prvl$|^t1dm_prvl$|^t2dm_prvl$|^obesity_prvl$",
        #       disnam,
        #       value = TRUE,
        #       invert = TRUE
        #     ),
        #     "cancer_prvl", "ctdra_prvl", "dm_prvl"
        #   )
        # private$diseasenam_hlp <- disnam

        # Fill cms_weights with 0 for conditions that are not in CMS
        # for (i in setdiff(private$diseasenam_hlp, names(private$cms_weights)))
        #   private$cms_weights[[i]] <- 0

        private$death_codes <- unlist(lapply(self$diseases, function(x)
          x$meta$mortality$code))
        private$death_codes[["alive"]] <- 0L

        invisible(self)
      },

      #' @description
      #' Runs a simulation
      #' @param mc A positive sequential integer vector with the Monte Carlo
      #'   iterations of synthetic population to simulate.
      #' @param multicore If TRUE run the simulation in parallel.
      #' @return The invisible self for chaining.
      run = function(mc, multicore = TRUE) {

        if (!is.integer(mc)) stop("mc need to be an integer vector")
        if (any(mc <= 0)) stop("mc need to be positive")

        # check if sequential vector. Necessary if design$sim_prm$n_synthpop_aggregation > 1
        if (anyNA(mc) || any(is.infinite(mc)) || length(mc) <= 1 ||
            diff(mc[1:2]) == 0 || diff(range(diff(mc))) > sqrt(.Machine$double.eps))
              stop("mc need to be a sequential integer vector")
        # NOTE mc is in fact mc_aggr. mc_ is the mc of the synthpop
        mc_sp <-
          (
            min(mc) * self$design$sim_prm$n_synthpop_aggregation -
              self$design$sim_prm$n_synthpop_aggregation + 1L
          ):(max(mc) * self$design$sim_prm$n_synthpop_aggregation)



        if (any(file.exists( # TODO fix when lifecourse is not saved
          file.path(
            self$design$sim_prm$output_dir,
            "lifecourse",
            paste0(mc, "_lifecourse.csv")
          )
        ))) {
          stop("Results from a previous simulation exists in the output folder.
               Please remove them before run a new one.")
        }



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
            mc_iter = mc_sp,
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

          lapply(mc_sp, private$run_sim)

          if (self$design$sim_prm$logs)
            private$time_mark("End of single-core run")

        }

        # process summarised outputs

        if (multicore) {

          if (Sys.info()["sysname"] == "Windows") {
            cl <-
              makeCluster(self$design$sim_prm$clusternumber) # used for clustering. Windows compatible
            registerDoParallel(cl)
          } else {
            registerDoParallel(self$design$sim_prm$clusternumber) # used for forking. Only Linux/OSX compatible
          }

          if (self$design$sim_prm$logs)
            private$time_mark("Start exporting summaries")

          xps_dt <- foreach(
             i = mc,
            .inorder = TRUE,
            .verbose = self$design$sim_prm$logs,
            .packages = c(
              "R6",
              "CKutils",
              "IMPACTncdEngl",
              "data.table"
            ),
            .export = NULL,
            .noexport = NULL # c("time_mark")
          ) %dopar% {

            pth <- private$output_dir(paste0("lifecourse/", i, "_lifecourse.csv.gz"))
            lc <-   fread(pth, stringsAsFactors = TRUE, key = c("pid", "year"))
            private$export_summaries(lc)
            NULL
          }

          if (exists("cl")) stopCluster(cl)

          if (self$design$sim_prm$logs) private$time_mark("End of exporting summuries")


        } else {
          if (self$design$sim_prm$logs)
            private$time_mark("Start of single-core run")

          lapply(mc, function(i) {
            pth <- private$output_dir(paste0("lifecourse/", i, "_lifecourse.csv.gz"))
            lc <-   fread(pth, stringsAsFactors = TRUE, key = c("pid", "year"))
            private$export_summaries(lc)
            NULL
          })

          if (self$design$sim_prm$logs)
            private$time_mark("End of single-core run")

        }



        for (i in mc) {



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
      #' Get the European Standardised Population 2013 by sex and dimd.
      #' @return A data.table with the European Standardised Population 2013.
      get_esp = function() {
        private$esp_weights
      },

      #' @description
      #' Get the disease multimorbidity weifgts (i.e. Cambridge Morbidity Score weights).
      #' @return A named vector with disease weights.
      get_mm_weights = function() {
        unlist(sapply(self$diseases, function(x) x$meta$diagnosis$mm_wt))
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
      death_codes = NA,
      # diseasenam_hlp = NA,
      esp_weights = data.table(),
      # cms_weights = c(
      #   "htn_prvl"      = 0.08,
      #   "andep_prvl"    = 0.5,
      #   "pain_prvl"     = 0.92,
      #   "helo_prvl"     = 0.09, # hearing loss
      #   "ibs_prvl"      = 0.21,
      #   "asthma_prvl"   = 0.19,
      #   "dm_prvl"       = 0.75, # t2dm + t1dm
      #   "chd_prvl"      = 0.49,
      #   "ckd_prvl"      = 0.53,
      #   "ckd4_prvl"     = 0.53,
      #   "ckd45_prvl"    = 0.53,
      #   "af_prvl"       = 1.34,
      #   "constip_prvl"  = 1.12,
      #   "stroke_prvl"   = 0.8,
      #   "copd_prvl"     = 1.46,
      #   "ctdra_prvl"    = 0.43, # connective tissue disorders + rheumatoid arthritis
      #   "cancer_prvl"   = 1.53,
      #   "alcpr_prvl"    = 0.65,
      #   "hf_prvl"       = 1.18,
      #   "dementia_prvl" = 2.50,
      #   "psychos_prvl"  = 0.64,
      #   "epilepsy_prvl" = 0.92
      # ),

      # helper function to calculate cms score for each row
      # disprvl the vector with the disease colnames private$diseasenam_hlp.
      # dt the sp$pop
      # cmswt the private$cms_weights
      # cms_hlpfn = function(disprvl, dt, cmswt = cmswt)
      #   clamp(dt[[disprvl]]) * cmswt[[disprvl]],



      # Runs the simulation in one core. mc is scalar
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
        # it doesn't matter if mc or mc_aggr is used in the above, because it is
        # only used for the RNG stream and the pid are different in each mc_aggr
        # pop

        sp$update_pop_weights()

        if (self$design$sim_prm$export_xps) {
          if (self$design$sim_prm$logs) message("Exporting exposures...")
          private$export_xps(sp)
        }

        nam <- c(self$design$sim_prm$cols_for_output,
                 grep("^cms_|_prvl$|_dgns$|_mrtl$", names(sp$pop), value = TRUE))
        sp$pop[, mc := sp$mc_aggr]

        # Prune pop (NOTE that assignment in the function env makes this data.table local)
        sp$pop <- sp$pop[all_cause_mrtl >= 0L &
                 year >= self$design$sim_prm$init_year &
                 between(age, self$design$sim_prm$ageL, self$design$sim_prm$ageH), ..nam]
        setkey(sp$pop, pid, year)
        sp$pop[, pid_mrk := mk_new_simulant_markers(pid)]

        # apply ESP weights
        to_agegrp(sp$pop, 5, 99)
        absorb_dt(sp$pop, private$esp_weights)
        sp$pop[, wt_esp := wt_esp * unique(wt_esp) / sum(wt_esp),
           keyby = .(year, agegrp, sex, dimd)]


        # combine all cancers (moved to C++)
        # sp$pop[, cancer_prvl := clamp(Reduce(`+`, .SD)), .SDcols = patterns("_ca_prvl$")]
        # sp$pop[, cancer_prvl := carry_forward_incr(cancer_prvl, pid_mrk, TRUE, 1L)]

        # combine ctd & ra (moved to C++)
        # sp$pop[, ctdra_prvl := clamp(Reduce(`+`, .SD)),
        #        .SDcols = patterns("^ctd_prvl$|^ra_prvl$")]
        # sp$pop[, ctdra_prvl := carry_forward_incr(ctdra_prvl, pid_mrk, TRUE, 1L)]

        # combine t1dm & t2dm (moved to C++)
        # sp$pop[, dm_prvl := clamp(Reduce(`+`, .SD)),
        #    .SDcols = patterns("^t1dm_prvl$|^t2dm_prvl$")]
        # sp$pop[, dm_prvl := carry_forward_incr(dm_prvl, pid_mrk, TRUE, 1L)]

        # sp$pop[, cms_score := Reduce(`+`, lapply(private$diseasenam_hlp,
        #                                    private$cms_hlpfn, sp$pop,
        #                                    private$cms_weights))]
        sp$pop[, `:=` (
          cms1st_cont_prvl   = carry_forward_incr(as.integer(cms_count == 1),
                                             pid_mrk, TRUE, 1L),
          cmsmm1_prvl   = carry_forward_incr(as.integer(cms_score >= 1),
                                             pid_mrk, TRUE, 1L),
          cmsmm1.5_prvl = carry_forward_incr(as.integer(cms_score >= 1.5),
                                             pid_mrk, TRUE, 1L),
          cmsmm2_prvl   = carry_forward_incr(as.integer(cms_score >= 2),
                                             pid_mrk, TRUE, 1L)
        )]


        # Write lifecourse
          if (self$design$sim_prm$logs) message("Exporting lifecourse...")
          fwrite_safe(sp$pop,
                      private$output_dir(paste0(
                        "lifecourse/", sp$mc_aggr, "_lifecourse.csv.gz"
                      )))

        # Export summaries
        # if (self$design$sim_prm$logs) message("Exporting summaries...")
        # strata <- setdiff(self$design$sim_prm$cols_for_output, c("age", "pid", "wt"))
        #
        # # Life expectancy NOTE for scale_up LE weights need to apply from the very beginning
        # fwrite_safe(sp$pop[all_cause_mrtl > 0, .("popsize" = (.N), LE = mean(age)),  keyby = strata],
        #             private$output_dir(paste0("summaries/", "le_out.csv"
        #             )))
        # fwrite_safe(sp$pop[all_cause_mrtl > 0, .("popsize" = sum(wt), LE = weighted.mean(age, wt)),  keyby = strata],
        #             private$output_dir(paste0("summaries/", "le_scale_up.csv"
        #             )))
        # fwrite_safe(sp$pop[all_cause_mrtl > 0, .("popsize" = sum(wt_esp), LE = weighted.mean(age, wt_esp)),  keyby = strata],
        #             private$output_dir(paste0("summaries/", "le_esp.csv"
        #             )))
        # # Life expectancy at 60
        #
        # fwrite_safe(sp$pop[all_cause_mrtl > 0 & age > 60, .("popsize" = (.N), LE60 = mean(age)),  keyby = strata],
        #             private$output_dir(paste0("summaries/", "le60_out.csv"
        #             )))
        # fwrite_safe(sp$pop[all_cause_mrtl > 0 & age > 60, .("popsize" = sum(wt), LE60 = weighted.mean(age, wt)),  keyby = strata],
        #             private$output_dir(paste0("summaries/", "le60_scale_up.csv"
        #             )))
        # fwrite_safe(sp$pop[all_cause_mrtl > 0 & age > 60, .("popsize" = sum(wt_esp), LE60 = weighted.mean(age, wt_esp)),  keyby = strata],
        #             private$output_dir(paste0("summaries/", "le60_esp.csv"
        #             )))
        # # Note: for less aggregation use wtd.mean with popsize i.e le_out[, weighted.mean(LE, popsize), keyby = year]
        #
        # # Healthy life expectancy
        # # TODO currently some individuals are counted more than once because
        # # disease counter and score can be reduced. Ideally only the first reach
        # # to the threshold should be counted
        # fwrite_safe(sp$pop[cms_count == 1L, .("popsize" = (.N), HLE = mean(age)),
        #                    keyby = strata],
        #             private$output_dir(paste0("summaries/", "hle_1st_cond_out.csv")))
        # fwrite_safe(sp$pop[cms_count == 1L,
        #                    .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
        #                    keyby = strata],
        #             private$output_dir(paste0(
        #               "summaries/", "hle_1st_cond_scale_up.csv"
        #             )))
        # fwrite_safe(sp$pop[cms_count == 1L,
        #                    .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
        #                    keyby = strata],
        #             private$output_dir(paste0("summaries/", "hle_1st_cond_esp.csv"
        #             )))
        #
        # fwrite_safe(sp$pop[cmsmm1.5_prvl == 1L, .("popsize" = (.N), HLE = mean(age)),
        #                    keyby = strata],
        #             private$output_dir(paste0("summaries/", "hle_cmsmm1.5_out.csv")))
        # fwrite_safe(sp$pop[cmsmm1.5_prvl == 1L,
        #                    .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
        #                    keyby = strata],
        #             private$output_dir(paste0(
        #               "summaries/", "hle_cmsmm1.5_scale_up.csv"
        #             )))
        # fwrite_safe(sp$pop[cmsmm1.5_prvl == 1L,
        #                    .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
        #                    keyby = strata],
        #             private$output_dir(paste0("summaries/", "hle_cmsmm1.5_esp.csv"
        #             )))
        #
        # strata <- c("agegrp", strata) # Need to be after LE
        #
        # fwrite_safe(sp$pop[, c("popsize" = (.N),
        #                        lapply(.SD, function(x) sum(x > 0))),
        #                    .SDcols = patterns("_prvl$"), keyby = strata],
        #             private$output_dir(paste0("summaries/", "prvl_out.csv"
        #             )))
        # fwrite_safe(sp$pop[, c("popsize" = sum(wt),
        #                        lapply(.SD, function(x, wt) sum((x > 0) * wt), wt)),
        #                    .SDcols = patterns("_prvl$"), keyby = strata],
        #             private$output_dir(paste0("summaries/", "prvl_scale_up.csv"
        #             )))
        # fwrite_safe(sp$pop[, c("popsize" = sum(wt_esp),
        #                        lapply(.SD, function(x, wt) sum((x > 0) * wt), wt_esp)),
        #                    .SDcols = patterns("_prvl$"), keyby = strata],
        #             private$output_dir(paste0("summaries/", "prvl_esp.csv"
        #             )))
        # fwrite_safe(sp$pop[, c("popsize" = (.N),
        #                        lapply(.SD, function(x) sum(x == 1))),
        #                    .SDcols = patterns("_prvl$"), keyby = strata],
        #             private$output_dir(paste0("summaries/", "incd_out.csv"
        #             )))
        # fwrite_safe(sp$pop[, c("popsize" = sum(wt),
        #                        lapply(.SD, function(x, wt) sum((x == 1) * wt), wt)),
        #                    .SDcols = patterns("_prvl$"), keyby = strata],
        #             private$output_dir(paste0("summaries/", "incd_scale_up.csv"
        #             )))
        # fwrite_safe(sp$pop[, c("popsize" = sum(wt_esp),
        #                        lapply(.SD, function(x, wt) sum((x == 1) * wt), wt_esp)),
        #                    .SDcols = patterns("_prvl$"), keyby = strata],
        #             private$output_dir(paste0("summaries/", "incd_esp.csv"
        #             )))
        # fwrite_safe(sp$pop[, c("popsize" = (.N),
        #                        "all_cause_mrtl" = sum(all_cause_mrtl > 0)),
        #                    keyby = strata],
        #             private$output_dir(paste0("summaries/", "mrtl_out.csv"
        #             )))
        # fwrite_safe(sp$pop[, c("popsize" = sum(wt),
        #                        "all_cause_mrtl" = sum((all_cause_mrtl > 0) * wt)),
        #                    keyby = strata],
        #             private$output_dir(paste0("summaries/", "mrtl_scale_up.csv"
        #             )))
        # fwrite_safe(sp$pop[, c("popsize" = sum(wt_esp),
        #                        "all_cause_mrtl" = sum((all_cause_mrtl > 0) * wt_esp)),
        #                    keyby = strata],
        #             private$output_dir(paste0("summaries/", "mrtl_esp.csv"
        #             )))
        #
        #
        #
        # # disease specific mortality
        # dis_mrtl_out <-
        #   dcast(
        #     sp$pop[, .("deaths" = (.N)),
        #            keyby = c(strata, "all_cause_mrtl")],
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
        #             private$output_dir(paste0("summaries/", "dis_mrtl_out.csv"
        #             )))
        #
        # dis_mrtl_out <- # scale up
        #   dcast(
        #     sp$pop[, .("deaths" = sum(wt)),
        #            keyby = c(strata, "all_cause_mrtl")],
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
        #             private$output_dir(paste0("summaries/", "dis_mrtl_scale_up.csv"
        #             )))
        #
        # dis_mrtl_out <- # scale up esp
        #   dcast(
        #     sp$pop[, .("deaths" = sum(wt_esp)),
        #            keyby = c(strata, "all_cause_mrtl")],
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
        #             private$output_dir(paste0("summaries/", "dis_mrtl_esp.csv"
        #             )))


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
          "cms_score"          = paste0("cms_score", scenario_suffix_for_pop),
          "cms_count"          = paste0("cms_count", scenario_suffix_for_pop),
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

      # function to export summaries from lifecourse files
      # lc is a lifecourse file
      export_summaries = function(lc) {
        if (self$design$sim_prm$logs) message("Exporting summaries...")
        strata <- setdiff(self$design$sim_prm$cols_for_output, c("age", "pid", "wt"))

        # Life expectancy NOTE for scale_up LE weights need to apply from the very beginning
        fwrite_safe(lc[all_cause_mrtl > 0, .("popsize" = (.N), LE = mean(age)),  keyby = strata],
                    private$output_dir(paste0("summaries/", "le_out.csv.gz"
                    )))
        fwrite_safe(lc[all_cause_mrtl > 0, .("popsize" = sum(wt), LE = weighted.mean(age, wt)),  keyby = strata],
                    private$output_dir(paste0("summaries/", "le_scale_up.csv.gz"
                    )))
        fwrite_safe(lc[all_cause_mrtl > 0, .("popsize" = sum(wt_esp), LE = weighted.mean(age, wt_esp)),  keyby = strata],
                    private$output_dir(paste0("summaries/", "le_esp.csv.gz"
                    )))
        # Life expectancy at 60

        fwrite_safe(lc[all_cause_mrtl > 0 & age > 60, .("popsize" = (.N), LE60 = mean(age)),  keyby = strata],
                    private$output_dir(paste0("summaries/", "le60_out.csv.gz"
                    )))
        fwrite_safe(lc[all_cause_mrtl > 0 & age > 60, .("popsize" = sum(wt), LE60 = weighted.mean(age, wt)),  keyby = strata],
                    private$output_dir(paste0("summaries/", "le60_scale_up.csv.gz"
                    )))
        fwrite_safe(lc[all_cause_mrtl > 0 & age > 60, .("popsize" = sum(wt_esp), LE60 = weighted.mean(age, wt_esp)),  keyby = strata],
                    private$output_dir(paste0("summaries/", "le60_esp.csv.gz"
                    )))
        # Note: for less aggregation use wtd.mean with popsize i.e le_out[, weighted.mean(LE, popsize), keyby = year]

        # Healthy life expectancy
        # TODO currently some individuals are counted more than once because
        # disease counter and score can be reduced. Ideally only the first reach
        # to the threshold should be counted
        fwrite_safe(lc[cms_count == 1L, .("popsize" = (.N), HLE = mean(age)),
                       keyby = strata],
                    private$output_dir(paste0("summaries/", "hle_1st_cond_out.csv.gz")))
        fwrite_safe(lc[cms_count == 1L,
                       .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
                       keyby = strata],
                    private$output_dir(paste0(
                      "summaries/", "hle_1st_cond_scale_up.csv.gz"
                    )))
        fwrite_safe(lc[cms_count == 1L,
                       .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
                       keyby = strata],
                    private$output_dir(paste0("summaries/", "hle_1st_cond_esp.csv.gz"
                    )))

        fwrite_safe(lc[cmsmm1.5_prvl == 1L, .("popsize" = (.N), HLE = mean(age)),
                       keyby = strata],
                    private$output_dir(paste0("summaries/", "hle_cmsmm1.5_out.csv.gz")))
        fwrite_safe(lc[cmsmm1.5_prvl == 1L,
                       .("popsize" = sum(wt), HLE = weighted.mean(age, wt)),
                       keyby = strata],
                    private$output_dir(paste0(
                      "summaries/", "hle_cmsmm1.5_scale_up.csv.gz"
                    )))
        fwrite_safe(lc[cmsmm1.5_prvl == 1L,
                       .("popsize" = sum(wt_esp), HLE = weighted.mean(age, wt_esp)),
                       keyby = strata],
                    private$output_dir(paste0("summaries/", "hle_cmsmm1.5_esp.csv.gz"
                    )))

        strata <- c("agegrp", strata) # Need to be after LE

        fwrite_safe(lc[, c("popsize" = (.N),
                           lapply(.SD, function(x) sum(x > 0))),
                       .SDcols = patterns("_prvl$"), keyby = strata],
                    private$output_dir(paste0("summaries/", "prvl_out.csv.gz"
                    )))
        fwrite_safe(lc[, c("popsize" = sum(wt),
                           lapply(.SD, function(x, wt) sum((x > 0) * wt), wt)),
                       .SDcols = patterns("_prvl$"), keyby = strata],
                    private$output_dir(paste0("summaries/", "prvl_scale_up.csv.gz"
                    )))
        fwrite_safe(lc[, c("popsize" = sum(wt_esp),
                           lapply(.SD, function(x, wt) sum((x > 0) * wt), wt_esp)),
                       .SDcols = patterns("_prvl$"), keyby = strata],
                    private$output_dir(paste0("summaries/", "prvl_esp.csv.gz"
                    )))
        fwrite_safe(lc[, c("popsize" = (.N),
                           lapply(.SD, function(x) sum(x == 1))),
                       .SDcols = patterns("_prvl$"), keyby = strata],
                    private$output_dir(paste0("summaries/", "incd_out.csv.gz"
                    )))
        fwrite_safe(lc[, c("popsize" = sum(wt),
                           lapply(.SD, function(x, wt) sum((x == 1) * wt), wt)),
                       .SDcols = patterns("_prvl$"), keyby = strata],
                    private$output_dir(paste0("summaries/", "incd_scale_up.csv.gz"
                    )))
        fwrite_safe(lc[, c("popsize" = sum(wt_esp),
                           lapply(.SD, function(x, wt) sum((x == 1) * wt), wt_esp)),
                       .SDcols = patterns("_prvl$"), keyby = strata],
                    private$output_dir(paste0("summaries/", "incd_esp.csv.gz"
                    )))
        fwrite_safe(lc[, c("popsize" = (.N),
                           "all_cause_mrtl" = sum(all_cause_mrtl > 0)),
                       keyby = strata],
                    private$output_dir(paste0("summaries/", "mrtl_out.csv.gz"
                    )))
        fwrite_safe(lc[, c("popsize" = sum(wt),
                           "all_cause_mrtl" = sum((all_cause_mrtl > 0) * wt)),
                       keyby = strata],
                    private$output_dir(paste0("summaries/", "mrtl_scale_up.csv.gz"
                    )))
        fwrite_safe(lc[, c("popsize" = sum(wt_esp),
                           "all_cause_mrtl" = sum((all_cause_mrtl > 0) * wt_esp)),
                       keyby = strata],
                    private$output_dir(paste0("summaries/", "mrtl_esp.csv.gz"
                    )))



        # disease specific mortality
        dis_mrtl_out <-
          dcast(
            lc[, .("deaths" = (.N)),
               keyby = c(strata, "all_cause_mrtl")],
            formula = as.formula(paste0(
              paste(strata, collapse = "+"), "~all_cause_mrtl"
            )),
            fill = 0L,
            value.var = "deaths"
          )

        setnames(dis_mrtl_out, as.character(private$death_codes),
                 names(private$death_codes), skip_absent = TRUE)
        dis_mrtl_out[, `:=` (
          popsize = Reduce(`+`, .SD),
          alive = NULL
        ), .SDcols = !strata]
        fwrite_safe(dis_mrtl_out,
                    private$output_dir(paste0("summaries/", "dis_mrtl_out.csv.gz"
                    )))

        dis_mrtl_out <- # scale up
          dcast(
            lc[, .("deaths" = sum(wt)),
               keyby = c(strata, "all_cause_mrtl")],
            formula = as.formula(paste0(
              paste(strata, collapse = "+"), "~all_cause_mrtl"
            )),
            fill = 0L,
            value.var = "deaths"
          )

        setnames(dis_mrtl_out, as.character(private$death_codes),
                 names(private$death_codes), skip_absent = TRUE)
        dis_mrtl_out[, `:=` (
          popsize = Reduce(`+`, .SD),
          alive = NULL
        ), .SDcols = !strata]
        fwrite_safe(dis_mrtl_out,
                    private$output_dir(paste0("summaries/", "dis_mrtl_scale_up.csv.gz"
                    )))

        dis_mrtl_out <- # scale up esp
          dcast(
            lc[, .("deaths" = sum(wt_esp)),
               keyby = c(strata, "all_cause_mrtl")],
            formula = as.formula(paste0(
              paste(strata, collapse = "+"), "~all_cause_mrtl"
            )),
            fill = 0L,
            value.var = "deaths"
          )

        setnames(dis_mrtl_out, as.character(private$death_codes),
                 names(private$death_codes), skip_absent = TRUE)
        dis_mrtl_out[, `:=` (
          popsize = Reduce(`+`, .SD),
          alive = NULL
        ), .SDcols = !strata]
        fwrite_safe(dis_mrtl_out,
                    private$output_dir(paste0("summaries/", "dis_mrtl_esp.csv.gz"
                    )))


        if (!self$design$sim_prm$keep_lifecourse) file.remove(pth)

        return(invisible(self))
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
