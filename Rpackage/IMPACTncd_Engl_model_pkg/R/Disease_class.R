## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by
## Chris Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2022 University of Liverpool, Chris Kypridemos
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


#' R6 Class representing an disease
#'
#' @description
#' A disease has a sim_prm list that holds the simulation parameters.
#'
#' @details
#' To be completed...
#'
#' @export
Disease <-
  R6::R6Class(
    classname = "Disease",
    # cloneable = FALSE, # cloneable is necessary for multi threading
    # public ------------------------------------------------------------------
    public = list(
      #' @field name The name of the disease.
      name = NA_character_,
      #' @field friendly_name A friendly name for the disease.
      friendly_name = NA_character_,
      #' @field meta Disease metadata including type.
      meta = NA_integer_,
      #' @field notes Any notes regarding the disease.
      notes = NA_character_,

      #' @description Create a new disease object.
      #' @param name A string with disease name.
      #' @param friendly_name A string with disease friendly name.
      #' @param RR A list of exposure objects.
      #' @param meta A list with the disease type and other information for
      #'   incidence, diagnosis, and mortality.
      #' @param notes A string with any notes.
      #' @param design_ A design object with the simulation parameters.
      #' @return A new `Disease` object.

      initialize = function(name, friendly_name, meta, notes = NA_character_,
                            design_ = design, RR) {
        if (!inherits(design_, "Design")) {
          stop("Argument design needs to be a Design object.")
        }
        if (!(is.character(name) && is.character(friendly_name))) {
          stop("Both arguments need to be strings.")
        }
        if (!all(sapply(RR, inherits, "Exposure"))) {
          stop("Argument RR needs to be a list of exposure object.")
        }

        self$name          <- name
        self$friendly_name <- friendly_name
        self$meta          <- meta
        self$notes         <- notes

        # Generate unique name using the relevant RR and lags
        rr <- RR[sapply(RR, `[[`, "outcome") == self$name]
        # above only contains exposures for this disease object
        # Reorder risk so smok_status & smok_cig is calculated before endsmok
        private$rr <-
          rr[order(match(sapply(rr, `[[`, "name"),
                         c("smok_status", "smok_cig")))]

        dqRNGkind("pcg64")
        private$seed <- abs(digest2int(name, seed = 230565490L))


        if (is.numeric(meta$incidence$type) && meta$incidence$type > 1L) {
          private$filenams$incd <- file.path(
            getwd(), "inputs", "disease_burden",
            paste0(self$name, "_incd", ".fst")
          )
          private$filenams$prvl <- file.path(
            getwd(), "inputs", "disease_burden",
            paste0(self$name, "_prvl", ".fst")
          )
          private$filenams$dur <- file.path(
            getwd(), "inputs", "disease_burden",
            paste0(self$name, "_dur", ".fst")
          )
        }
        if (is.numeric(meta$mortality$type)) {
          private$filenams$ftlt <- file.path(
            getwd(), "inputs", "disease_burden",
            paste0(self$name, "_ftlt", ".fst")
          )
        }


        if (is.numeric(meta$incidence$type) && meta$incidence$type > 1L) {
          keys <- lapply(private$filenams[names(private$filenams) != "dur"],
                         function(x) metadata_fst(x)$keys[[1]])

          if (!all(sapply(keys, identical, "year")))
            stop("1st key need to be year")

          private$incd_indx <-
            read_fst(private$filenams$incd, as.data.table = TRUE
                     )[, .(from = min(.I), to = max(.I)), keyby = "year"]
          private$prvl_indx <-
            read_fst(private$filenams$prvl, as.data.table = TRUE
                     )[, .(from = min(.I), to = max(.I)), keyby = "year"]
        }
        if (is.numeric(meta$mortality$type)) {
          private$ftlt_indx <-
            read_fst(private$filenams$ftlt, as.data.table = TRUE
                     )[, .(from = min(.I), to = max(.I)), keyby = "year"]
        }

        # TODO add check for stop('For type 1 incidence aggregation of RF need
        # to be "any" or "all".')

        nam <- paste0("prb_", self$name, "_incd")

        if (self$meta$incidence$type == 3) {
          nam <- paste0(
            nam, "_no",
            paste(self$meta$incidence$influenced_by_disease_name,
              collapse = ""
            )
          )
        }
        private$incd_colnam <- nam
        private$dgns_colnam <- paste0("prb_", self$name, "_dgns")

        invisible(self)
      },

      #' @description Get disease incident probability.
      #' @param design_ A design object with the simulation parameters.
      #' @param diseases_ A list of Disease objects
      #' @param popsize The population size for each stratum
      #' @param check Check for NAs in parf_dt.
      #' @return The invisible self for chaining.

      gen_parf = function(design_ = design, diseases_ = diseases,
                          popsize = 100, check = design_$sim_prm$logs) {

        if (is.numeric(self$meta$incidence$type) && self$meta$incidence$type == 1L) {
          # Early break for type 1 incidence
          return(invisible(self))
        }


        if (!inherits(design_, "Design")) {
          stop("Argument design_ needs to be a Design object.")
        }
        if (!all(sapply(diseases_, inherits, "Disease"))) {
          stop("Argument diseases_ needs to be a list of disease object.")
        }

        checksum <-
          digest(list(
            design_$sim_prm[c("init_year", "ageL", "ageH")],
            lapply(private$rr, function(x) x$get_input_rr()),
            lapply(private$rr, `[[`, "lag"),
            lapply(private$rr, `[[`, "distribution")
          ),
          seed = 305834490L
          )

        dir_path <- file.path(getwd(), "simulation", "parameters")

        parf_filenam <- file.path(
          dir_path,
          paste0("PARF_", self$name, "_", checksum, ".fst")
        )

        if (!file.exists(parf_filenam)) {
          # start if file not exist
          # NOTE multicore crashes
          if (sum(dim(private$incd_indx)) > 0) {
            ttt <- self$get_incd(design_$sim_prm$init_year)
          } else {
            ttt <- self$get_ftlt(design_$sim_prm$init_year)
          }

          ttt <- CJ(
            age = seq(design_$sim_prm$ageL, design_$sim_prm$ageH),
            sex = levels(ttt$sex),
            dimd = levels(ttt$dimd),
            ethnicity = levels(ttt$ethnicity),
            sha = levels(ttt$sha),
            year = design_$sim_prm$init_year
          )

          ttt[, qimd := fcase(
            dimd == "1 most deprived",  "1 most deprived",
            dimd == "2",  "1 most deprived",
            dimd == "3",  "2",
            dimd == "4",  "2",
            dimd == "5",  "3",
            dimd == "6",  "3",
            dimd == "7",  "4",
            dimd == "8",  "4",
            dimd == "9",  "5 least deprived",
            dimd == "10 least deprived", "5 least deprived"
          )]
          ttt[, qimd := factor(qimd,
            levels = c(
              "1 most deprived",
              as.character(2:4),
              "5 least deprived"
            )
          )]

          ttt <- clone_dt(ttt, 10)

          # NOTE future and mclapply do not work here for some reason
          if (Sys.info()["sysname"] == "Windows") {
            cl <-
              makeCluster(design_$sim_prm$clusternumber) # used for clustering. Windows compatible
            registerDoParallel(cl)
          } else {
            registerDoParallel(design_$sim_prm$clusternumber) # used for forking. Only Linux/OSX compatible
          }
          xps_dt <- foreach(
            mc_iter = seq(1, (popsize/10L)),
            .inorder = FALSE,
            .verbose = design_$sim_prm$logs,
            .packages = c(
              "R6",
              "gamlss.dist",
              "dqrng",
              "CKutils",
              "IMPACTncdEngl",
              "data.table"
            ),
            .export = NULL,
            .noexport = NULL # c("time_mark")
          ) %dopar% {
            private$gen_sp_forPARF(mc_iter, ff = ttt, design_ = design_)

          }
          if (exists("cl")) stopCluster(cl)

          if (design_$sim_prm$logs) message("End of parallelisation for PARF.")

          ans <- list()
          setattr(ans, "class", "SynthPop") # to dispatch
          ans$pop <- rbindlist(xps_dt)
          ans$pop [, `:=` (pid = .I, pid_mrk = TRUE)] # TODO add check to avoid intmax limit
          ans$mc <- 0L
          ans$mc_aggr <- 0L

          # NOTE xps_dt does not contain disease init prevalence. I simulate
          # here as set_init_prvl expects a synthpop and not a data.table as
          # input. All relevant risk factors for the diseases need to be
          # available in xps_dt. Currently this requires manual setup of the
          # logic if RF for relevant diseases altered.
          # TODO implement automation of the logic above based on topological
          # ordering.
          # Generate diseases that act as exposures
          xps <- "t2dm_prvl"
          if (xps %in% sapply(private$rr, `[[`, "name")) {
            lag <- private$rr[[paste0(xps, "~", self$name)]]$lag
            ans$pop[, year := year - lag]
            design_$sim_prm$init_year <- design_$sim_prm$init_year - lag
            diseases_$t2dm$set_init_prvl(ans, design_)
            design_$sim_prm$init_year <- design_$sim_prm$init_year + lag
            ans$pop[, year := year + lag]
          }
          xps <- "af_prvl"
          if (xps %in% sapply(private$rr, `[[`, "name")) {
            lag <- private$rr[[paste0(xps, "~", self$name)]]$lag
            ans$pop[, year := year - lag]
            design_$sim_prm$init_year <- design_$sim_prm$init_year - lag
            diseases_$af$set_init_prvl(ans, design_)
            design_$sim_prm$init_year <- design_$sim_prm$init_year + lag
            ans$pop[, year := year + lag]
          }

          self$set_rr(ans, design_, forPARF = TRUE)

          nam <- grep("_rr$", names(ans$pop), value = TRUE)

          parf_dt <-
            ans$pop[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH),
              .(parf = 1 - 1 / (sum(Reduce(`*`, mget(nam))) / .N)),
              keyby = .(age, sex, dimd, ethnicity, sha)
            ]

          if (sum(dim(private$incd_indx)) > 0) {
            tt <- self$get_incd(design_$sim_prm$init_year)
            nam <- "p0"
          } else {
            tt <- self$get_ftlt(design_$sim_prm$init_year)
            setnames(tt, "mu2", "mu")
            nam <- "m0"
          }

          absorb_dt(parf_dt, tt)
          parf_dt[, (nam) := mu * (1 - parf)]
          if (check && anyNA(parf_dt)) {
            warning("NAs in parf. Replacing any missing p0 with incidence.")
            parf_dt[is.na(get(nam)), (nam) := mu]
          }
          parf_dt[, c("year", "mu") := NULL]

          if (!dir.exists(dir_path)) {
            dir.create(dir_path, showWarnings = FALSE, recursive = TRUE)
          }
          write_fst(parf_dt, parf_filenam)
        } # end if file not exist

        private$parf <- read_fst(parf_filenam, as.data.table = TRUE)

        invisible(self)
      },

      #' @description Deletes the PARF file from disk. If both arguments
      #'   missing, deletes ALL parf files for this disease.
      #' @param design_ A design object with the simulation parameters.
      #' @param invert deletes all other disease relevant PARF file except those
      #'   that are associated to the curent settings.
      #' @return The invisible self for chaining.

      del_parf_file = function(design_ = design, invert = FALSE) {
        if (!missing(design_)) {
          if (!inherits(design_, "Design")) {
            stop("Argument design_ needs to be a Design object.")
          }

          checksum <-
            digest(list(
              design_$sim_prm[c("init_year", "ageL", "ageH")],
              lapply(private$rr, function(x) x$get_input_rr()),
              lapply(private$rr, `[[`, "lag"),
              lapply(private$rr, `[[`, "distribution")
            ),
            seed = 305834490L
            )
          dir_path <- file.path(getwd(), "simulation", "parameters")

          parf_filenam <- file.path(
            dir_path,
            paste0("PARF_", self$name, "_", checksum, ".fst")
          )

          if (invert) {
            parf_filenam2 <- list.files(dir_path,
                                        pattern = paste0(
                                          "^PARF_",
                                          self$name, ".*\\.fst$"
                                        ),
                                        full.names = TRUE
            )

            parf_filenam <- setdiff(parf_filenam2, parf_filenam)
          }

          file.remove(parf_filenam)
        } else {
          dir_path <- file.path(getwd(), "simulation", "parameters")
          parf_filenams <- list.files(dir_path,
            pattern = paste0(
              "^PARF_",
              self$name, ".*\\.fst$"
            ),
            full.names = TRUE
          )
          file.remove(parf_filenams)
        }
        invisible(self)
      },


      #' @description Set disease prevalence & diagnosisin a new col in sp$pop.
      #' @param sp A synthetic population.
      #' @param design_ A design object with the simulation parameters.
      #' @return The invisible self for chaining.
      set_init_prvl = function(sp, design_ = design) {
        # TODO correlate with other diseases prevalence
        if (is.numeric(self$meta$incidence$type)) {
          if (!inherits(sp, "SynthPop")) {
            stop("Argument sp needs to be a SynthPop object.")
          }
          if (!inherits(design_, "Design")) {
            stop("Argument design_ needs to be a Design object.")
          }

          namprvl <- paste0(self$name, "_prvl")
          if (namprvl %in% names(sp$pop)) {
            stop(
              "A column named ",
              namprvl,
              " already exists in sp$pop.
              Please delete it and run set_init_prvl() afterwards."
            )
          }

          if (self$meta$incidence$type == 1L) {
            self$set_rr(sp, design_, forPARF = FALSE)
            riskcolnam <- grep("_rr$",
                               names(private$risks),
                               value = TRUE,
                               perl = TRUE)
            if (length(riskcolnam) == 1L) {
              thresh <- as.integer(private$risks[[riskcolnam]])
            }
            if (length(riskcolnam) > 1L && self$meta$incidence$aggregation == "any") {
              thresh <- as.integer(private$risks[, do.call(pmax, .SD), .SDcols = riskcolnam])
            }
            if (length(riskcolnam) > 1L && self$meta$incidence$aggregation == "all") {
              thresh <- as.integer(private$risks[, Reduce(`*`, .SD), .SDcols = riskcolnam])
            }

            set(sp$pop, NULL, namprvl, thresh)

            sp$pop[year > design_$sim_prm$init_year & age > design_$sim_prm$ageL,
                   (namprvl) := 0L]

            sp$pop[, carry_forward_incr(get((namprvl)), pid_mrk,
                                        recur = self$meta$incidence$can_recur,
                                        y = 1L, byref = TRUE)]
            invisible(self)

          } else { # if incidence type not 1

            dqRNGkind("pcg64")
            dqset.seed(private$seed, stream = sp$mc * 10 + 1L) # not mc_aggr
            set.seed(private$seed + sp$mc * 10 + 1L) # for sample_int_expj
            # First find out how many prevalent cases by pop subgroup
            tbl <- self$get_prvl(design_$sim_prm$init_year
            )[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH)]
            strata <- setdiff(names(tbl), c("mu")) # age removed to avoid very small groups. To be replaced by agegroups below

            sp$pop[tbl, on = .NATURAL, mu := i.mu]
            sp$pop[
              year == design_$sim_prm$init_year,
              (namprvl) := as.integer(dqrunif(.N) < mu)
            ]
            sp$pop[, mu := NULL]

            # for new entries we assume no time trends on initial prevalence
            tbl <- tbl[age == design_$sim_prm$ageL][, year := NULL]
            sp$pop[tbl, on = .NATURAL, mu := i.mu]
            sp$pop[
              year > design_$sim_prm$init_year &
                age == design_$sim_prm$ageL,
              (namprvl) := as.integer(dqrunif(.N) < mu)
            ]
            setnafill(sp$pop, "c", 0L, cols = namprvl)
            sp$pop[, mu := NULL]

            # Then select individuals with higher risk to have higher probability
            # to be a prevalent case. I estimate weights based on relevant RF for
            # each disease. Note that RR and sampling weights are equivalent here.
            # The RR are for incident. One would expect prevalent cases to be
            # healthier because of survival of the fittest and because some may
            # have better RF control post diagnosis. For that I will arbitrarily
            # assume that the risk for prevalence is half of that for incidence.


            # ncases is the number of prevalent cases expected in each stratum
            sp$pop[year >= design_$sim_prm$init_year,
                   ncases := sum(get(namprvl)), by = eval(strata)]
            # Generate unique name using the relevant RR and lags

            # TODO below should apply only to strata with ncases > 0 for efficiency
            self$set_rr(sp, design_, forPARF = TRUE)

            nam <- grep("_rr$", names(sp$pop), value = TRUE) # necessary because above forPARF = TRUE

            sp$pop[, disease_wt := (Reduce(`*`, .SD)), .SDcols = nam]
            sp$pop[, (nam) := NULL]
            # adjust for prevalent risk half of incident risk
            sp$pop[, disease_wt := ((disease_wt - 1) * 0.5) + 1]
            ss <- sp$pop[ncases > 0,][, .("pid" = pid[sample_int_expj(unique(.N), unique(ncases), disease_wt)]), by = eval(strata)]
            sp$pop[, (namprvl) := 0L]
            sp$pop[ss, on = c("pid", "year"), (namprvl) := 1L]
            sp$pop[, c("ncases", "disease_wt") := NULL]

            # set duration
            dqset.seed(private$seed, stream = sp$mc * 10 + 2L) # not mc_aggr
            tbl <- read_fst(private$filenams$dur, as.data.table = TRUE)
            col_nam <- setdiff(names(tbl), intersect(names(sp$pop), names(tbl)))
            tbl[, (namprvl) := 1L]
            absorb_dt(sp$pop, tbl)
            fn <- paste0("q", self$meta$diagnosis$duration_distr)
            sp$pop[get(namprvl) == 1L,
                   (namprvl) := 1L + do.call(fn, c(p = list(dqrunif(.N)), .SD)),
                   .SDcols = col_nam]
            sp$pop[, (col_nam) := NULL]

            sp$pop[, (namprvl) := carry_backward_decr(get(namprvl), pid_mrk)] # necessary for c++

          } # End if incidence type not 1


          # TODO this only makes sense when probability of diagnosis is 1
          namdgns <- paste0(self$name, "_dgns")
          set(sp$pop, NULL, namdgns, 0L)
          sp$pop[
            get(namprvl) > 0 & year == design_$sim_prm$init_year & dqrunif(.N) < self$meta$diagnosis$probability,
            (namdgns) := get(namprvl)
          ]
          sp$pop[, (namdgns) := carry_backward_decr(get(namdgns), pid_mrk)]
        }

        invisible(self)
      },



      #' @description Get disease incident probability.
      #' @param year_ A vector of years to return. All if missing.
      #' @return A data.table with disease incident probabilities unless
      #'   incidence type: Universal when it returns data.table(NULL).
      get_incd = function(year_) {
        if (sum(dim(private$incd_indx)) > 0) {
          if (missing(year_)) {
            ro <- list(from = 1, to = NULL)
          } else {
            ro <- private$incd_indx[
              year %in% sort(year_),
              .("from" = min(from), "to" = max(to))
            ]
          }
          out <- read_fst(
            private$filenams$incd,
            as.data.table = TRUE,
            from = ro$from,
            to = ro$to
          )
        } else {
          message("Incidence type: ", self$meta$incidence$type)
          out <- data.table(NULL)
        }
        return(out)
      },


      #' @description Set disease incidence probability in a new col in sp$pop.
      #' @param sp A synthetic population.
      #' @param design_ A design object with the simulation parameters.
      #' @param checkNAs If `TRUE`, prints the table of NAs before they get
      #'   overwritten with 1. Note that for some exposures, NAs are expected
      #'   for certain levels of exposure (i.e. for active days).
      #' @param forPARF Set TRUE when applied on the specialised forPARF
      #'   SynthPop
      #' @return The invisible self for chaining.

      set_rr = function(sp, design_ = design,
                        checkNAs = design_$sim_prm$logs, forPARF = FALSE) {
        # For incd type 1 forPARF = TRUE is meaningless but gen_parf() skips
        # this type so we are good here.


        if (!inherits(sp, "SynthPop")) {
          stop("Argument sp needs to be a SynthPop object.")
        }
        if (!inherits(design_, "Design")) {
          stop("Argument design_ needs to be a Design object.")
        }

        lapply(private$rr, function(x)
          x$xps_to_rr(sp, design_, checkNAs = checkNAs, forPARF = forPARF))

        if (!forPARF) {
          nam <- grep("_rr$", names(sp$pop), value = TRUE)
          private$risks <- sp$pop[, .SD, .SDcols = c("pid", "year", nam)]
          sp$pop[, (nam) := NULL]
        }
        return(invisible(self))
      },


      #' @description Set disease incident probability in a new col in sp$pop.
      #' @param sp A synthetic population.
      #' @param design_ A design object with the simulation parameters.
      #' @return The invisible self for chaining.
      set_incd_prb = function(sp, design_ = design) {

        if (!inherits(design_, "Design")) {
          stop("Argument design_ needs to be a Design object.")
        }
        if (!inherits(sp, "SynthPop")) {
          stop("Argument sp needs to be a SynthPop object.")
        }

        if (is.numeric(self$meta$incidence$type)) {
          # check that the stored risks are for the sp (i.e. no rows deleted)
          if (!identical(sp$pop$pid, private$risks$pid)) {
            stop("Stored risks are for a different synthetic population")
          }


          if (private$incd_colnam %in% names(sp$pop)) {
            stop(
              "A column named ", nam,
              " already exists in sp$pop. ",
              "Please delete it and run set_init_prvl() afterwards."
            )
          }

          # Get colnames in risk that end with _rr but exclude the influence by
          # diseases
          if (self$meta$incidence$type == 3) {
            riskcolnam <- grep(
              paste0(
                "^((?!",
                paste(self$meta$incidence$influenced_by_disease_name,
                  collapse = "|"
                ),
                ").)*_rr$"
              ),
              names(private$risks),
              value = TRUE,
              perl = TRUE
            )
          } else {
            riskcolnam <- grep("_rr$",
              names(private$risks),
              value = TRUE,
              perl = TRUE
            )
          }

          if (self$meta$incidence$type == 1L) {
            if (length(riskcolnam) == 1L) {
              thresh <- as.numeric(private$risks[[riskcolnam]])
            }
            if (length(riskcolnam) > 1L && self$meta$incidence$aggregation == "any") {
              thresh <- as.numeric(private$risks[, do.call(pmax, .SD), .SDcols = riskcolnam])
            }
            if (length(riskcolnam) > 1L && self$meta$incidence$aggregation == "all") {
              thresh <- as.numeric(private$risks[, Reduce(`*`, .SD), .SDcols = riskcolnam])
            }

            set(sp$pop, NULL, private$incd_colnam, thresh)

            } else {
            # if incident$type not 1
            risk_product <-
              private$risks[, Reduce(`*`, .SD), .SDcols = riskcolnam]
            absorb_dt(sp$pop, private$parf)
            set(sp$pop,
                NULL,
                private$incd_colnam,
                clamp(sp$pop$p0 * risk_product))
            # NOTE product above not expected to be = to incidence because p0
            # estimated using mean lags and RR, while each mc run samples from
            # their distribution.

            if (design_$sim_prm$export_PARF) {
              path <- file.path(design_$sim_prm$output_dir, "parf")
              filenam <- file.path(path, "parf.csv")
              if (!dir.exists(path)) {
                dir.create(path, showWarnings = FALSE, recursive = TRUE)
              }
              parf_dt <-
                sp$pop[, .(parf = unique(parf), pop_size = sum(wt)),
                       keyby = .(age, sex, dimd, ethnicity, sha)]
              parf_dt[, `:=`(disease = self$name, mc = sp$mc)] # not sp$mc_aggr
              # EXAMPLE parf[, weighted.mean(parf, pop_size), keyby = sex]
              fwrite_safe(parf_dt, filenam)
            }

            sp$pop[, c("p0", "parf") := NULL]
          } # End of incident$type not 1
        } # End if incident$type numeric
        invisible(self)
      },

      #' @description Get disease duration distribution parameters.
      #' @return A data.table with duration distribution parameters. unless
      #'   incidence type: Universal when it returns data.table(NULL).
      get_dur = function() {
        if (!is.null(private$filenams$dur)) {
          out <- read_fst(
            private$filenams$dur,
            as.data.table = TRUE
          )
        } else {
          message("Incidence type: ", self$meta$incidence$type)
          out <- data.table(NULL)
        }
        return(out)
      },

      #' @description Get disease prevalent probability.
      #' @param year_ A vector of years to return. All if missing.
      #' @return A data.table with disease prevalent probabilities unless
      #'   incidence type: Universal when it returns data.table(NULL).
      get_prvl = function(year_) {
        if (sum(dim(private$prvl_indx)) > 0) {
          if (missing(year_)) {
            ro <- list(from = 1, to = NULL)
          } else {
            ro <- private$prvl_indx[
              year %in% sort(year_),
              .("from" = min(from), "to" = max(to))
            ]
          }
          out <- read_fst(
            private$filenams$prvl,
            as.data.table = TRUE,
            from = ro$from,
            to = ro$to
          )
        } else {
          message("Incidence type: ", self$meta$incidence$type)
          out <- data.table(NULL)
        }
        return(out)
      },



      #' @description Get disease case fatality probability.
      #' @param year_ A vector of years to return. All if missing.
      #' @return A data.table with disease case fatality probabilities unless
      #'   mortality type: Non-fatal when it returns data.table(NULL).
      get_ftlt = function(year_) {
        if (sum(dim(private$ftlt_indx)) > 0) {
          if (missing(year_)) {
            ro <- list(from = 1, to = NULL)
          } else {
            ro <- private$ftlt_indx[
              year %in% sort(year_),
              .("from" = min(from), "to" = max(to))
            ]
          }
          out <- read_fst(
            private$filenams$ftlt,
            as.data.table = TRUE,
            from = ro$from,
            to = ro$to
          )
        } else {
          message("Mortality type: ", self$meta$mortality$type)
          out <- data.table(NULL)
        }
        return(out)
      },


      #' @description Set disease case fatality in a new col in sp$pop.
      #' @param sp A synthetic population.
      #' @param design_ A design object with the simulation parameters.
      #' @return The invisible self for chaining.

      set_mrtl_prb = function(sp, design_ = design) {
        if (is.numeric(self$meta$mortality$type)) {
          if (!inherits(sp, "SynthPop")) {
            stop("Argument sp needs to be a SynthPop object.")
          }
          if (!inherits(design_, "Design")) {
            stop("Argument design_ needs to be a Design object.")
          }

          tt <-
            self$get_ftlt(seq(design_$sim_prm$init_year, sp$pop[, max(year)]))
          hlp <- grep("^mu", names(tt), value = TRUE)
          if (length(hlp) == 2L) {
            nam <- paste0("prb_", self$name, "_mrtl", 1:2)
            setnames(tt, c("mu1", "mu2"), nam)
            private$mrtl2flag <- TRUE
          } else {
            nam <- paste0("prb_", self$name, "_mrtl2")
            setnames(tt, hlp, nam)
          }

          absorb_dt(sp$pop, tt)
          setnafill(sp$pop,
            type = "const",
            fill = 0,
            cols = nam
          )
        }
        invisible(self)
      },

      #' @description Set diagnosis probability in a new col in sp$pop.
      #' @param sp A synthetic population.
      #' @param design_ A design object with the simulation parameters.
      #' @return The invisible self for chaining.

      set_dgns_prb = function(sp, design_ = design) {
        if (is.numeric(self$meta$diagnosis$type)) {
          if (!inherits(sp, "SynthPop")) {
            stop("Argument sp needs to be a SynthPop object.")
          }
          if (!inherits(design_, "Design")) {
            stop("Argument design_ needs to be a Design object.")
          }

          set(sp$pop, NULL, private$dgns_colnam, self$meta$diagnosis$probability)
        }

        invisible(self)
      },


      #' @description Get seed for RNG.
      #' @return A seed for the RNG that is produced by the digest of disease
      #'   name and outcome.
      get_seed = function() {
        private$seed
      },

      #' @description Get the list of rr for all relevant exposures.
      #' @return A list of exposure objects.
      get_rr = function() {
        private$rr
      },


      #' @description Get the risks for all individuals in a synthetic
      #'   population.
      #' @return A data.table with columns for pid, year, and all associated
      #'   risks.
      get_risks = function() {
        private$risks
      },

      #' @description Get the PARF by age/sex/dimd/ethnicity/sha.
      #' @return A data.table with PARF.
      get_parf = function() {
        private$parf
      },

      #' @description Harmonises classes and levels between the synthetic
      #'   population and the incidence/prevalence/fatality tables. It saves the
      #'   harmonised table to disk, overwriting the existing one.
      #' @param sp A synthetic population.
      #' @return The invisible self for chaining.

      harmonise_epi_tables = function(sp) {
        if (!inherits(sp, "SynthPop")) {
          stop("Argument sp needs to be a SynthPop object.")
        }

        for (i in seq_along(private$filenams)) {
          tbl <- read_fst(private$filenams[[i]], as.data.table = TRUE)
          k <- key(tbl)
          for (j in names(tbl)) {
            if (j %in% names(sp$pop)) {
              if (inherits(sp$pop[[j]], "integer") && !inherits(tbl[[j]], "integer")) {
                tbl[, (j) := as.integer(get(j))]
              }
              if (inherits(sp$pop[[j]], "numeric") && !inherits(tbl[[j]], "numeric")) {
                tbl[, (j) := as.numeric(get(j))]
              }
              if (inherits(sp$pop[[j]], "character") && !inherits(tbl[[j]], "character")) {
                tbl[, (j) := as.character(get(j))]
              }
              if (inherits(sp$pop[[j]], "factor")) {
                # irrespective of class(j) to make sure that levels are the same
                # and in the right order.
                if (j == "dimd" && levels(tbl[[j]])[1] != "1 most deprived") {
                  tbl[, (j) := factor(get(j), levels = as.character(10:1), labels = levels(sp$pop[[j]]))]
                } else {
                  tbl[, (j) := factor(get(j), levels = levels(sp$pop[[j]]))]
                }
              }
            }
          }

          if (grepl("_ftlt.fst$", private$filenams[[i]]) &&
            "mu" %in% names(tbl)) {
            setnames(tbl, "mu", "mu2")
          }

          if (anyNA(tbl)) stop("NAs in ", private$filenams[[i]])
          setkeyv(tbl, k)
          write_fst(tbl, private$filenams[[i]])
        }
      },

      #' @description Returns a list to pass to the C++ side for Chris' parser.
      #' @param sp A synthetic population.
      #' @param design_ A design object with the simulation parameters.
      #' @param scenario_suffix the suffix to identify columns from different
      #'   scenarios
      #' @return A list.

      to_cpp = function(sp, design_ = design, scenario_suffix = "") {
        if (!inherits(sp, "SynthPop")) {
          stop("Argument sp needs to be a SynthPop object.")
        }
        if (!inherits(design_, "Design")) {
          stop("Argument design_ needs to be a Design object.")
        }

        out <- list()
        out <-
          list(
            "incidence" = NULL,
            "diagnosis" = NULL,
            "mortality" = NULL,
            "seed" = fifelse(
              design_$sim_prm$kismet,
              private$seed,
              abs(digest2int(
                paste0(self$name, scenario_suffix),
                seed = 230565490L
              ))
            )
          )


        out[["incidence"]] <- list(
          "type" = fifelse(
            is.numeric(self$meta$incidence$type),
            paste0("Type", self$meta$incidence$type),
            as.character(self$meta$incidence$type)
          ),
          "prevalence" = paste0(self$name, "_prvl", scenario_suffix),
          "probability" = paste0(private$incd_colnam, scenario_suffix),
          "can_recur" = self$meta$incidence$can_recur
        )
        if (is.null(out$incidence$can_recur)) out$incidence <- within(out$incidence, rm("can_recur"))
        if (out$incidence$type == "Universal") out$incidence <- within(out$incidence, rm("prevalence", "probability"))

        if (self$meta$incidence$type == 3) {
          influenced_by_incd <- list()

          for (i in self$meta$incidence$influenced_by_disease_name) {
            influenced_by_incd[[paste0(i, "_prvl")]] <-
              list(
                "multiplier" = paste0(self$name, "_incd_", i, "_prvl_mltp"),
                "lag" = private$rr[[paste0(i, "_prvl", "~", self$name)]]$
                  get_lag(fifelse(design_$sim_prm$stochastic, sp$mc_aggr, 0L)))
          } # end for loop over influenced_by_disease_name
          out[["incidence"]][["influenced_by"]] <- influenced_by_incd
        } # end if incidence type 3


        if (!is.null(self$meta$diagnosis$type)) {
          out[["diagnosis"]] <- list(
            "type" = fifelse(
              is.numeric(self$meta$diagnosis$type),
              paste0("Type", self$meta$diagnosis$type),
              as.character(self$meta$diagnosis$type)
            ),
            "diagnosed" = paste0(self$name, "_dgns", scenario_suffix),
            "probability" = paste0(
              private$dgns_colnam,
              scenario_suffix
            )
          )
        } else {
          out <- within(out, rm("diagnosis"))
        }

        if (!is.null(self$meta$mortality$type)) {
          out[["mortality"]] <- list(
            "type" = fifelse(
              is.numeric(self$meta$mortality$type),
              paste0("Type", self$meta$mortality$type),
              as.character(self$meta$mortality$type)
            ),
            "probability" = paste0("prb_", self$name, "_mrtl2", scenario_suffix),
            "code" = self$meta$mortality$code
          )

          if (private$mrtl2flag) {
            out[["mortality"]][["probability1styear"]] <-
              paste0("prb_", self$name, "_mrtl1", scenario_suffix)
          }


          if (self$meta$mortality$type == 3) {
            influenced_by_mrtl <- list()
            for (i in self$meta$mortality$influenced_by_disease_name) {
              influenced_by_mrtl[[paste0(i, "_prvl")]] <-
                list(
                  "multiplier" = paste0(self$name, "_mrtl_", i, "_prvl_mltp"),
                  "lag" = private$rr[[paste0(i, "_prvl", "~", self$name)]]$
                    get_lag(fifelse(design_$sim_prm$stochastic, sp$mc_aggr, 0L))
                )
            } # end for loop over influenced_by_disease_name
            out[["mortality"]][["influenced_by"]] <-
              influenced_by_mrtl
          } # end if mortality type 3
        } else { # end of not null mortality
          out <- within(out, rm("mortality"))
        } # end of null mortality
        out
      },

      #' @description Returns a list to pass to the C++ side for Peter's parser.
      #' @param sp A synthetic population.
      #' @param design_ A design object with the simulation parameters.
      #' @param scenario_suffix the suffix to identify columns from different
      #'   scenarios
      #' @return A list.

      to_cpp_Peter = function(sp, design_ = design, scenario_suffix = "") {
        if (!inherits(sp, "SynthPop")) {
          stop("Argument sp needs to be a SynthPop object.")
        }
        if (!inherits(design_, "Design")) {
          stop("Argument design_ needs to be a Design object.")
        }

        out <- list()
        out <- list("incidence" = NULL, "diagnosis" = NULL, "mortality" = NULL)


        out[["incidence"]] <- list(
          "type" = fifelse(
            is.numeric(self$meta$incidence$type),
            paste0("Type", self$meta$incidence$type),
            as.character(self$meta$incidence$type)
          ),
          "prevalence" = list("column_name" = paste0(self$name, "_prvl", scenario_suffix)),
          "incidence" = list("column_name" = paste0(private$incd_colnam, scenario_suffix)),
          "rn" = list("column_name" = paste0(
            "rn_", self$name, "_incd",
            fifelse(design_$sim_prm$kismet, "", scenario_suffix)
          )),
          "can_recur" = self$meta$incidence$can_recur
        )

        if (self$meta$incidence$type == 3) {
          influenced_by_incd <- list()

          for (i in self$meta$incidence$influenced_by_disease_name) {
            influenced_by_incd[[i]] <-
              list(
                "multiplier" = list(
                  "column_name" = paste0(self$name, "_incd_", i, "_prvl_mltp")),
                  "lag" = private$rr[[paste0(i, "_prvl", "~", self$name)]]$
                  get_lag(fifelse(design_$sim_prm$stochastic, sp$mc_aggr, 0L)))
          } # end for loop over influenced_by_disease_name
          out[["incidence"]][["influenced_by"]] <- influenced_by_incd
        } # end if incidence type 3


        if (!is.null(self$meta$diagnosis$type)) {
          out[["diagnosis"]] <- list(
            "type" = fifelse(
              is.numeric(self$meta$diagnosis$type),
              paste0("Type", self$meta$diagnosis$type),
              as.character(self$meta$diagnosis$type)
            ),
            "diagnosed" = list("column_name" = paste0(self$name, "_dgns", scenario_suffix)),
            "probability" = list("column_name" = paste0(private$dgns_colnam, scenario_suffix)),
            "rn" = list("column_name" = paste0(
              "rn_",
              self$name,
              "_dgn",
              fifelse(design_$sim_prm$kismet, "", scenario_suffix)
            ))
          )
        } else {
          out <- within(out, rm("diagnosis"))
        }

        out[["mortality"]] <- list(
          "type" = fifelse(
            is.numeric(self$meta$mortality$type),
            paste0("Type", self$meta$mortality$type),
            as.character(self$meta$mortality$type)
          ),
          "probability" = list("column_name" = paste0(
            "prb_", self$name, "_mrtl2", scenario_suffix
          )),
          "rn" = list("column_name" = paste0(
            "rn_",
            self$name,
            "_mrtl",
            fifelse(design_$sim_prm$kismet, "", scenario_suffix)
          )),
          "code" = self$meta$mortality$code
        )

        if (self$meta$mortality$type == 3) {
          influenced_by_mrtl <- list()
          for (i in self$meta$mortality$influenced_by_disease_name) {
            influenced_by_mrtl[[i]] <-
              list(
                "multiplier" = list("column_name" = paste0(self$name, "_mrtl_",
                                                           i, "_prvl_mltp")),
                "lag" = private$rr[[paste0(i, "_prvl", "~", self$name)]]$
                  get_lag(fifelse(design_$sim_prm$stochastic, sp$mc_aggr, 0L)))
          } # end for loop over influenced_by_disease_name
          out[["mortality"]][["influenced_by"]] <- influenced_by_mrtl
        } # end if mortality type 3
        out
      },

      #' @description Print the simulation parameters.
      #' @return The invisible self for chaining.
      print = function() {
        print(paste0("Disease name:       ", self$name))
        print(paste0("Meta incidence:     ", self$meta$incidence))
        print(paste0("Meta diagnosis:     ", self$meta$diagnosis))
        print(paste0("Meta mortality:     ", self$meta$mortality))
        print(paste0("Notes:              ", self$notes))

        invisible(self)
      }
    ), # end of public

    # private ------------------------------------------------------------------
    private = list(
      seed = NA_integer_,
      filenams = list(),
      mrtl2flag = FALSE, # TRUE if separate mortality for 1st year of disease is present
      incd_colnam = NA,
      dgns_colnam = NA,
      mrtl_colnam = NA,
      incd_indx = data.table(NULL),
      prvl_indx = data.table(NULL),
      ftlt_indx = data.table(NULL),
      parf = data.table(NULL),
      risks = data.table(NULL), # holds the risks for all individuals
      rr = list(), # holds the list of relevant RR
      gen_sp_forPARF =
        function(mc_, ff, design_) {

          dqRNGkind("pcg64")
          set.seed(private$seed + mc_)
          dqset.seed(private$seed, mc_)

          cm_mean <- as.matrix(
            read_fst(
              "./inputs/exposure_distributions/exposure_corr_mean.fst",
              as.data.table = TRUE
            ),
            rownames = "rn"
          )

          r <-
            which(
              rownames(cm_mean) %in% c(
                "education_r", "income_r", "bp_med_r", "af_r",
                "famcvd_r", "ckd_r", "dm_r", "dm_dgn_r"
              )
            )

          rank_mtx <- generate_corr_unifs(nrow(ff), cm_mean[-r, -r])

          # Restrict the range of some RNs to avoid unrealistic exposures
          # This scaling does not affect correlations
          # /0.999 because I multiplied all the columns below
          rank_mtx <- rank_mtx * 0.999
          rank_mtx[, "frtpor_r"] <- rank_mtx[, "frtpor_r"] * 0.99 / 0.999
          rank_mtx[, "vegpor_r"] <- rank_mtx[, "vegpor_r"] * 0.93 / 0.999
          rank_mtx[, "smok_cig_ex_r"] <-
            rank_mtx[, "smok_cig_ex_r"] * 0.99 / 0.999
          rank_mtx[, "totalwu_r"] <- rank_mtx[, "totalwu_r"] * 0.99 / 0.999
          rank_mtx[, "smok_quit_yrs_r"] <-
            rank_mtx[, "smok_quit_yrs_r"] * 0.99 / 0.999
          rank_mtx[, "smok_dur_ex_r"] <-
            rank_mtx[, "smok_dur_ex_r"] * 0.99 / 0.999
          rank_mtx[, "smok_dur_curr_r"] <-
            rank_mtx[, "smok_dur_curr_r"] * 0.88 / 0.999

          # sum((cor(rank_mtx) - cm_mean) ^ 2)

          rank_mtx <- data.table(rank_mtx)

          # NOTE rankstat_* is unaffected by the RW. Stay constant through the lifecourse
          ff[, c(
            # "rank_education",
            # "rank_income",
            "rank_pa",
            "rank_fruit",
            "rank_veg",
            "rankstat_smok",
            "rankstat_smok_quit_yrs",
            "rankstat_smok_dur_ex",
            "rankstat_smok_dur_curr",
            "rankstat_smok_cig_ex",
            "rankstat_smok_cig_curr",
            "rank_ets",
            "rank_alcohol",
            "rank_bmi",
            "rank_sbp",
            "rank_tchol",
            "rank_hdl",
            "rank_statin_px"
          ) := rank_mtx]

          rm(rank_mtx)


          setkeyv(ff, c("year", "age", "sex", "sha", "dimd", "ethnicity"))


          if (max(ff$age) > 90L) {
            ff[, age100 := age]
            ff[age > 90L, age := 90L]
          }

          # Generate active days ----
          xps <- c("active_days", "met", "t2dm_prvl") # t2dm require MET
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/active_days_table.fst",
                as.data.table = TRUE
              )
            nam <- intersect(names(ff), names(tbl))
            ff[tbl, active_days_curr_xps := (rank_pa > pa0) + (rank_pa > pa1) + (rank_pa > pa2) +
              (rank_pa > pa3) + (rank_pa > pa4) + (rank_pa > pa5) + (rank_pa > pa6),
            on = nam
            ]
            ff[, rank_pa := NULL]
            ff[, year := year + lag]
          }

          # Generate MET ----
          xps <- c("met", "t2dm_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]

            ff[, met_curr_xps := as.integer(floor(active_days_curr_xps * (3L + qbinom(dqrunif(.N), 8, 3 / 11)) *
              (30 + qexp(dqrunif(.N), 1 / 7)) / 100))]
            ff[, year := year + lag]
          }

          # Generate fruit consumption (ZISICHEL) ----
          xps <- c("fruit", "t2dm_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]

            tbl <-
              read_fst("./inputs/exposure_distributions/frtpor_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, fruit_curr_xps :=
              my_qZISICHEL(rank_fruit,
                mu, sigma, nu, tau,
                n_cpu = design_$sim_prm$n_cpu
              ) * 80L] # g/d
            ff[, (col_nam) := NULL]
            ff[, rank_fruit := NULL]
            ff[, year := year + lag]
          }
          # Generate veg consumption (DEL) ----
          xps <- "veg"
          if (xps %in% sapply(private$rr, `[[`, "name")) {
            lag <- private$rr[[paste0(xps, "~", self$name)]]$lag
            ff[, year := year - lag]

            tbl <-
              read_fst("./inputs/exposure_distributions/vegpor_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, veg_curr_xps :=
              my_qDEL(rank_veg, mu, sigma, nu, n_cpu = design_$sim_prm$n_cpu) * 80L] # g/d
            ff[, (col_nam) := NULL]
            ff[, rank_veg := NULL]
            ff[, year := year + lag]
          }

          # Smoking ----
          xps <- c("smok_status", "smok_cig", "smok_packyrs",
                   "ets", "t2dm_prvl", "af_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (all(xps %in% sapply(private$rr, `[[`, "name"))) {
              stop("smok_status & smok_cig cannot be both risk factors for a disease.")
            }

            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            }

            if (xps[[2]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[2]], "~", self$name)]]$lag
            }

            if (xps[[3]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[3]], "~", self$name)]]$lag
            }

            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_status_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, smok_status_curr_xps := my_qMN4(rankstat_smok, mu, sigma, nu)] # for calibration
            ff[, (col_nam) := NULL]

            # Assign smok_quit_yrs when pid_mrk == true (the first year an individual enters the simulation)
            # I could use these estimates for calibration but I need to calculate mortality first
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_quit_yrs_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            set(ff, NULL, "smok_quit_yrs", 0L)
            ff[
              smok_status_curr_xps %in% 2:3,
              smok_quit_yrs := my_qDPO(rankstat_smok_quit_yrs, mu, sigma)
            ]
            ff[, rankstat_smok_quit_yrs := NULL]
            ff[, (col_nam) := NULL]

            # Assign smok_dur_ex when pid_mrk == true (the first year an individual enters the simulation)
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_dur_ex_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            set(ff, NULL, "smok_dur", 0L)
            ff[smok_status_curr_xps %in% 2:3, smok_dur := my_qDPO(rankstat_smok_dur_ex, mu, sigma)]
            ff[, rankstat_smok_dur_ex := NULL]
            ff[, (col_nam) := NULL]

            # Assign smok_dur_curr when pid_mrk == true (the first year an individual enters the simulation)
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_dur_curr_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[smok_status_curr_xps == 4, smok_dur := as.integer(round(qNBI(rankstat_smok_dur_curr, mu, sigma)))]
            ff[, rankstat_smok_dur_curr := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }

          # Assign smok_cig_curr when pid_mrk == true (the first year an individual enters the simulation)
          xps <- c("smok_cig", "smok_packyrs",
                   "t2dm_prvl", "af_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            set(ff, NULL, "smok_cig_curr_xps", 0L)

            tbl <-
              read_fst("./inputs/exposure_distributions/smok_cig_curr_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[
              smok_status_curr_xps == 4L,
              smok_cig_curr_xps :=
                qZINBI(rankstat_smok_cig_curr, mu, sigma, nu)
            ]
            ff[, (col_nam) := NULL]

            # Assign smok_cig_ex when pid_mrk == true (the first year an individual enters the simulation)
            # ff[smok_status == 2, smok_cig := 1L]


            tbl <-
              read_fst("./inputs/exposure_distributions/smok_cig_ex_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[
              smok_status_curr_xps == 3L,
              smok_cig_curr_xps := my_qZABNB(rankstat_smok_cig_ex,
                mu,
                sigma,
                nu,
                tau,
                n_cpu = design_$sim_prm$n_cpu
              )
            ]
            ff[, (col_nam) := NULL]
            ff[, smok_status_curr_xps := factor(smok_status_curr_xps)]
            ff[, year := year + lag]
          }

          if ("smok_packyrs" %in% sapply(private$rr, `[[`, "name"))
            ff[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps *
                                          smok_dur_curr_xps / 20))]

          # Generate ETS (BI) ----

          # Note at the moment this is independent of smoking prevalence TODO
          # calculate how many each smoker pollutes by year, SHA (not qimd) to
          # be used in scenarios. Ideally correct for mortality
          xps <- c("ets", "t2dm_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/ets_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, ets_curr_xps := as.integer(rank_ets < mu)]
            ff[, rank_ets := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }
          # Generate alcohol (ZINBI) ----
          xps <- c("alcohol", "t2dm_prvl", "af_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/alcohol_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, alcohol_curr_xps := as.integer(
              qZINBI(rank_alcohol, mu, sigma, nu))]
            ff[, rank_alcohol := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }

          # Generate BMI (BCPEo) ----
          xps <- c("bmi", "t2dm_prvl", "af_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/bmi_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, bmi_curr_xps := my_qBCPEo(rank_bmi, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
            ff[, rank_bmi := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }
          # Generate SBP (BCPEo) ----
          xps <- c("sbp", "af_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/sbp_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, sbp_curr_xps := my_qBCPEo(rank_sbp, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
            ff[, rank_sbp := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }
          # Generate tchol (BCT) ----
          xps <- c("tchol", "statin_px", "t2dm_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/tchol_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, tchol_curr_xps := my_qBCT(rank_tchol, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
            ff[, rank_tchol := NULL]
            ff[, (col_nam) := NULL]

            # Generate HDL (to tchol ratio) (GB1) ----

            # NOTE this very highly correlated with hdl level (~0.76) and
            #  highly to tchol (~-0.47). The latter is captured by the correlated RNs
            tbl <-
              read_fst("./inputs/exposure_distributions/hdl_to_tchol_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, tchol_hdl_ratio := 1 / qGB1(rank_hdl, mu, sigma, nu, tau)]
            ff[, rank_hdl := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }
          # Generate statins medication (BI) -----
          xps <- c("statin_px", "t2dm_prvl")
          if (any(xps %in% sapply(private$rr, `[[`, "name"))) {
            if (xps[[1]] %in% sapply(private$rr, `[[`, "name")) {
              lag <- private$rr[[paste0(xps[[1]], "~", self$name)]]$lag
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]

            ff[, `:=`(tchol = round(clamp(tchol_curr_xps, 2, 12), 0))]
            tbl <-
              read_fst("./inputs/exposure_distributions/statin_px_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            absorb_dt(ff, tbl)
            ff[, statin_px_curr_xps := as.integer(rank_statin_px < mu)]
            ff[, rank_statin_px := NULL]
            ff[, (col_nam) := NULL]
            ff[, `:=`(tchol = NULL)]
            ff[, year := year + lag]
          }

          nam <- grep("rank", names(ff), value = T)
          if (length(nam) > 0) ff[, (nam) := NULL]

          if ("age100" %in% names(dt)) {
            dt[, age := NULL]
            setnames(dt, "age100", "age")
          }
          invisible(ff)
        }
    ) # end of private
  )
