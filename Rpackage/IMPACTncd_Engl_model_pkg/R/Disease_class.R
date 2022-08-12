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
                            design_, RR) {
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
        # Reorder risk so smok_status & smok_cig is calculated before quit_yrs
        private$rr <-
          rr[order(match(sapply(rr, `[[`, "name"),
                         c("smok_status", "smok_cig", "smok_packyrs")))]

        dqRNGkind("pcg64")
        private$seed <- abs(digest2int(name, seed = 230565490L))

        db <- file.path(getwd(), "inputs", "disease_burden")
        it <- vector("character")
        if (is.numeric(meta$incidence$type) && meta$incidence$type > 1L) {
          private$filenams$incd <- file.path(db,
            paste0(self$name, "_incd", ".fst")
          )
          private$filenams$incd_indx <- file.path(db,
            paste0(self$name, "_incd_indx", ".fst")
          )
          private$filenams$prvl <- file.path(db,
            paste0(self$name, "_prvl", ".fst")
          )
          private$filenams$prvl_indx <- file.path(db,
            paste0(self$name, "_prvl_indx", ".fst")
          )
          private$filenams$dur <- file.path(db,
            paste0(self$name, "_dur", ".fst")
          )
          it <- c(it, "incd", "prvl")
        }

        if (is.numeric(meta$mortality$type)) {
          private$filenams$ftlt <- file.path(db,
            paste0(self$name, "_ftlt", ".fst")
          )
          private$filenams$ftlt_indx <- file.path(db,
            paste0(self$name, "_ftlt_indx", ".fst")
          )
          it <- c(it, "ftlt")

        }


        # TODO add check for stop('For type 1 incidence aggregation of RF need
        # to be "any" or "all".')

        private$incd_colnam  <- paste0("prb_", self$name, "_incd")
        private$dgns_colnam  <- paste0("prb_", self$name, "_dgns")
        private$mrtl_colnam2 <- paste0("prb_", self$name, "_mrtl2") # Only for mrtl 2

        private$chksum <-
          digest(list(
            design_$sim_prm[c("init_year", "ageL", "ageH", "apply_RR_to_mrtl2",
                              "model_trends_in_redidual_incd")],
            lapply(private$rr, function(x)
              x$get_input_rr()),
            lapply(private$rr, `[[`, "lag"),
            lapply(private$rr, `[[`, "distribution")
          ))

        private$parf_dir <- file.path(getwd(), "simulation", "parf")
        if (!dir.exists(private$parf_dir)) dir.create(private$parf_dir)
        private$parf_filenam <- file.path(
          private$parf_dir,
          paste0("PARF_", self$name, "_", private$chksum, ".fst")
        )

        keys <- sapply(private$filenams[names(private$filenams) %in% it],
                       function(x) metadata_fst(x)$keys[[1]])

        if (length(keys) > 0 && !all(sapply(keys, identical, "year")))
          stop("1st key need to be year")

        # Logic to ensure _indx files are up to date.
        snfile <- file.path(db, paste0(".", self$name, "_file_snapshot.qs"))
        if (file.exists(snfile)) snapshot <- changedFiles(qread(snfile))

        if (!file.exists(snfile) ||
            any(nzchar(snapshot$added),
                nzchar(snapshot$deleted),
                nzchar(snapshot$changed))) {

          for (i in it) {
            finn <- paste0(i, "_indx")
            if (file.exists(private$filenams[[finn]]))
              file.remove(private$filenams[[finn]])

            private[[finn]] <-
              read_fst(private$filenams[[i]], as.data.table = TRUE, columns = "year"
              )[, .(from = min(.I), to = max(.I)), keyby = "year"]

            write_fst(private[[finn]], private$filenams[[finn]], 100L)
          }

          if (file.exists(snfile)) file.remove(snfile)
          qsave(
            fileSnapshot(
              private$parf_dir,
              timestamp = NULL,
              md5sum = TRUE,
              recursive = FALSE,
              pattern = self$name
            ), snfile
          )
        } else { # if indx up to date
          for (i in it) {
            private[[paste0(i, "_indx")]] <-
              read_fst(private$filenams[[paste0(i, "_indx")]], as.data.table = TRUE)
          }
        }

        invisible(self)
      },

      #' @description Generates PARF and stores it to disk if one doesn not
      #'   exists already.
      #' @param design_ A design object with the simulation parameters.
      #' @param diseases_ A list of Disease objects
      #' @param popsize The population size for each stratum
      #' @param check Check for NAs in parf_dt.
      #' @param keep_intermediate_file Whether to keep the intermediate synthpop file
      #' @return The PARF data.table if it was created, otherwise `NULL`.

      gen_parf_files = function(design_ = design, diseases_ = diseases,
                                popsize = 100, check = design_$sim_prm$logs,
                                keep_intermediate_file = TRUE) {

        if ((is.numeric(self$meta$incidence$type) &&
             self$meta$incidence$type < 2L) ||
            length(private$rr) == 0L) {
          # Early break for type 1 incidence and diseases with no RF
          return(NULL)
        }

        if (!inherits(design_, "Design")) {
          stop("Argument design_ needs to be a Design object.")
        }
        if (!all(sapply(diseases_, inherits, "Disease"))) {
          stop("Argument diseases_ needs to be a list of disease object.")
        }
        if (!file.exists(private$parf_filenam)) {
        tmpfile <- file.path(private$parf_dir,
                             paste0("PARF_", self$name, "_", digest(sort(
                               sapply(private$rr, `[[`, "name")
                             )), ".qs"))

        if (file.exists(tmpfile)) {
          if (design_$sim_prm$logs) message("Reading file from cache.")
          ans <- qread(tmpfile)
          setDT(ans$pop)
        } else {
          if (design_$sim_prm$logs) message("No available cached file.")

          self$del_parf_file(invert = TRUE) # Delete old versions

          # start if file not exist
          if (sum(dim(private$incd_indx)) > 0) {
            ttt <- self$get_incd(design_$sim_prm$init_year)
          } else {
            ttt <- self$get_ftlt(design_$sim_prm$init_year)
          }

          ttt <- CJ(
            age = seq(design_$sim_prm$ageL, design_$sim_prm$ageH),
            sex = ttt$sex,
            dimd = ttt$dimd,
            ethnicity = ttt$ethnicity,
            sha = ttt$sha,
            year = design_$sim_prm$init_year,
            unique = TRUE
          )

          lv <- c("1 most deprived", as.character(2:4), "5 least deprived")

          ttt[, qimd := fcase(
            dimd == "1 most deprived", factor(lv[[1]], levels = lv),
            dimd == "2", factor(lv[[1]], levels = lv),
            dimd == "3", factor(lv[[2]], levels = lv),
            dimd == "4", factor(lv[[2]], levels = lv),
            dimd == "5", factor(lv[[3]], levels = lv),
            dimd == "6", factor(lv[[3]], levels = lv),
            dimd == "7", factor(lv[[4]], levels = lv),
            dimd == "8", factor(lv[[4]], levels = lv),
            dimd == "9", factor(lv[[5]], levels = lv),
            dimd == "10 least deprived", factor(lv[[5]], levels = lv)
          )]

          ttt <- clone_dt(ttt, 10, idcol = NULL)

          # NOTE future and mclapply do not work here for some reason
          if (Sys.info()["sysname"] == "Windows") {
            cl <-
              makeCluster(design_$sim_prm$clusternumber) # used for clustering. Windows compatible
            registerDoParallel(cl)
          } else {
            registerDoParallel(design_$sim_prm$clusternumber) # used for forking. Only Linux/OSX compatible
          }
          xps_dt <- foreach(
            mc_iter = seq(1, (popsize / 10L)),
            .inorder = FALSE,
            .options.multicore = list(preschedule = FALSE),
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
            private$gen_sp_forPARF(mc_iter, ff = ttt, design_ = design_,
                                   diseases_ = diseases_)

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
          # available in xps_dt.

          # Generate diseases that act as exposures
          xps_dep <- private$get_xps_dependency_tree(x = self$name, dssl = diseases_)
          xps_dep <- xps_dep[grepl("_prvl$", xpscol)]
          setkey(xps_dep, xpscol)

          for (xps in paste0(names(diseases_), "_prvl")) {
            if (xps %in% xps_dep$xpscol) {
              lag <- xps_dep[xps, max(lag)] # See note below in line ~ 1681 about max
              ans$pop[, year := year - lag]
              design_$sim_prm$init_year <-
                design_$sim_prm$init_year - lag
              disnam <- gsub("_prvl$", "", xps)
              diseases_[[disnam]]$set_init_prvl(ans, design_)
              design_$sim_prm$init_year <-
                design_$sim_prm$init_year + lag
              ans$pop[, year := year + lag]
            }
          }
          if (design_$sim_prm$logs) message("Saving parf cache.")
          qsave(ans, tmpfile)
        } # end tmpfile bypass
        self$set_rr(ans, design_, forPARF = TRUE)

        nam <- grep("_rr$", names(ans$pop), value = TRUE)
        # risks <- ans$pop[, .SD, .SDcols = c("pid", "year", nam)]

        parf_dt <-
          ans$pop[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH),
                  .(parf = 1 - .N / sum(Reduce(`*`, .SD))),
                  keyby = .(age, sex, dimd, ethnicity, sha), .SDcols = nam
          ]

        if (sum(dim(private$incd_indx)) > 0) {
          if (design_$sim_prm$logs) message("Estimating p0.")

          if (design_$sim_prm$model_trends_in_redidual_incd) {
            # TODO calculate PARF per year. Now it is assumed static as in init year
            yrs <- seq(
              design_$sim_prm$init_year,
              design_$sim_prm$init_year + design_$sim_prm$sim_horizon_max
            )
          } else {
            yrs <- design_$sim_prm$init_year
          }

          tt <- self$get_incd(yrs)
          nam <- "p0"

          parf_dt <- clone_dt(parf_dt, length(yrs)) # works even if length(yrs) == 1
          parf_dt[, `:=` (
            year = .id - 1L + design_$sim_prm$init_year,
            .id = NULL
          )]
          lookup_dt(parf_dt, tt, check_lookup_tbl_validity = FALSE)
          setnafill(parf_dt, "c", fill = 0, cols = "mu") # fix for prostate and breast cancer
          parf_dt[, (nam) := mu * (1 - parf)]
          parf_dt[, "mu" := NULL]

        }

        if (sum(dim(private$ftlt_indx)) > 0) {
          if (design_$sim_prm$logs) message("Estimating m0.")

          # Re-estimate parf without exposures and only for diseases when
          # apply_rr_to_mrtl2 is FALSE
          if (!design_$sim_prm$apply_RR_to_mrtl2) {
            riskcolnam <- grep(
              paste0(
                "^((?!",
                paste(self$meta$mortality$influenced_by_disease_name,
                      collapse = "|"),
                ").)*_rr$"
              ),
              names(nam),
              value = TRUE,
              perl = TRUE
            )

            if (length(riskcolnam) > 0) {
              parf_dt_mrtl <-
                ans$pop[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH),
                        .(parf_mrtl = 1 - .N / sum(Reduce(`*`, .SD))),
                        keyby = .(age, sex, dimd, ethnicity, sha), .SDcols = riskcolnam
                ]
              absorb_dt(parf_dt, parf_dt_mrtl)
            }
          }


          yrs <- seq(
            design_$sim_prm$init_year,
            design_$sim_prm$init_year + design_$sim_prm$sim_horizon_max
          )
          tt <- self$get_ftlt(yrs)
          if ("mu" %in% names(tt)) stop("mu in ftlt file. Please rename to mu2.")

          setnames(tt, "mu2", "mu")
          if ("mu1" %in% names(tt)) tt[, mu1 := NULL]
          # nam <- "m0"

          if (!all(yrs %in% unique(parf_dt$years))) { # TODO safer logic here
            parf_dt <- clone_dt(parf_dt, length(yrs))
            parf_dt[, year := .id - 1L + design_$sim_prm$init_year]
            parf_dt[, .id := NULL]
          }
          lookup_dt(parf_dt, tt, check_lookup_tbl_validity = FALSE)
          setnafill(parf_dt, "c", fill = 0, cols = "mu") # fix for prostate and breast cancer

          if ("parf_mrtl" %in% names(parf_dt)) {
            parf_dt[, "m0" := mu * (1 - parf_mrtl)]
            parf_dt[, parf_mrtl := NULL]
          } else {
            parf_dt[, "m0" := mu * (1 - parf)]
          }
          parf_dt[, "mu" := NULL]
        }

        if (uniqueN(parf_dt$year) == 1L) parf_dt[, year := NULL]

        if (check && anyNA(parf_dt)) {
          warning("NAs in parf.")
          print(summary(parf_dt))

          # parf_dt[is.na(get(nam)), (nam) := mu]
        }
        if (design_$sim_prm$logs) message("Saving parf file ", private$parf_filenam)
        write_fst(parf_dt, private$parf_filenam, 100L)
        if (!keep_intermediate_file) file.remove(tmpfile)

        parf_dt
        } else NULL
      },

      #' @description Read PARF file from disk. If missing, generates PARF and
      #'   writes it to disk.
      #' @param sp A synthpop object
      #' @param design_ A design object with the simulation parameters.
      #' @param diseases_ A list of Disease objects
      #' @param popsize The population size for each stratum
      #' @param check Check for NAs in parf_dt.
      #' @param keep_intermediate_file Whether to keep the intermediate synthpop file
      #' @return The invisible self for chaining.

      gen_parf = function(sp = sp, design_ = design, diseases_ = diseases,
                          popsize = 100, check = design_$sim_prm$logs,
                          keep_intermediate_file = TRUE) {

        # TODO add logic to delete the intermediate synthpop file outside this
        # function

        if ((is.numeric(self$meta$incidence$type) &&
            self$meta$incidence$type == 1L) ||
            length(private$rr) == 0L) {
          # Early break for type 1 incidence and diseases with no RF
          return(invisible(self))
        }


        if (!inherits(sp, "SynthPop")) {
          stop("Argument sp needs to be a SynthPop object.")
        }
        if (!inherits(design_, "Design")) {
          stop("Argument design_ needs to be a Design object.")
        }
        if (!all(sapply(diseases_, inherits, "Disease"))) {
          stop("Argument diseases_ needs to be a list of disease object.")
        }

        # Logic to ensure parf files are regenerated when disease incidence
        # change.
        # TODO this is currently too strict. Changes in prvl and dur files
        # shouldn't trigger the deletion of the parf file
        snfile <- file.path(getwd(), "inputs", "disease_burden",
                            paste0(".", self$name, "_file_snapshot_parf.qs"))
        if (file.exists(snfile)) snapshot <- changedFiles(qread(snfile))

        if (!file.exists(snfile) ||
            any(nzchar(snapshot$added),
                nzchar(snapshot$deleted),
                nzchar(snapshot$changed))) {

            self$del_parf_file()

          if (file.exists(snfile)) file.remove(snfile)
          qsave(
            fileSnapshot(
              private$parf_dir,
              timestamp = NULL,
              md5sum = TRUE,
              recursive = FALSE,
              pattern = paste0("^", self$name, "_incd|^", self$name, "_ftlt")
            ), snfile
          )
          }

        if (file.exists(private$parf_filenam)) {
          if (design_$sim_prm$logs) message("Reading parf file from disk.")

          tt <- read_fst(private$parf_filenam, as.data.table = TRUE)
          colnam <-
            setdiff(names(tt), intersect(names(sp$pop), names(tt)))
          private$parf <- tt[sp$pop, on = .NATURAL, ..colnam]

        } else { # if file not exist
          if (design_$sim_prm$logs) message("Generating new parf file.")

          # shortcut for when parallel part is successful but function crashes
          # after it. This logic caches the synthpop for parf that is
          # independent of the RR but only depends on the causal structure.
          # Hence this is using another checksam that only depends on the names
          # of xps stored in rr

          parf_dt <- self$gen_parf_files(design_, diseases_,
                                     popsize, check,
                                     keep_intermediate_file)
          colnam <-
            setdiff(names(parf_dt), intersect(names(sp$pop), names(parf_dt)))
          private$parf <- parf_dt[sp$pop, on = .NATURAL, ..colnam]
        } # end if file not exist


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

          if (self$meta$incidence$type == 0L) {
            set(sp$pop, NULL, namprvl, 0L)
          } else if (self$meta$incidence$type == 1L) {
            self$set_rr(sp, design_, forPARF = FALSE)
            riskcolnam <- grep("_rr$",
                               names(sp$get_risks(self$name)),
                               value = TRUE,
                               perl = TRUE)
            if (length(riskcolnam) == 1L) {
              thresh <- as.integer(sp$get_risks(self$name)[[riskcolnam]])
            }
            if (length(riskcolnam) > 1L && self$meta$incidence$aggregation == "any") {
              thresh <- as.integer(sp$get_risks(self$name)[, do.call(pmax, .SD),
                                                 .SDcols = riskcolnam])
            }
            if (length(riskcolnam) > 1L && self$meta$incidence$aggregation == "all") {
              thresh <- as.integer(sp$get_risks(self$name)[, Reduce(`*`, .SD),
                                                 .SDcols = riskcolnam])
            }

            set(sp$pop, NULL, namprvl, thresh)

            sp$pop[year > design_$sim_prm$init_year & age > design_$sim_prm$ageL,
                   (namprvl) := 0L]

            sp$pop[, carry_forward_incr(get((namprvl)), pid_mrk,
                                        recur = self$meta$incidence$can_recur,
                                        y = 1L, byref = TRUE)]
            invisible(self)

          } else if (self$meta$incidence$type > 1L) { # if incidence type not 0 or 1

            dqRNGkind("pcg64")
            dqset.seed(private$seed, stream = sp$mc * 10 + 1L) # not mc_aggr
            set.seed(private$seed + sp$mc * 10 + 1L) # for sample_int_expj
            # First find out how many prevalent cases by pop subgroup
            tbl <- rbind(
              self$get_prvl(design_$sim_prm$init_year
            )[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH)],
            self$get_prvl(seq(design_$sim_prm$init_year + 1L,
                              design_$sim_prm$init_year +
                                design_$sim_prm$sim_horizon_max)
            )[age == design_$sim_prm$ageL])

            strata <- setdiff(names(tbl), c("mu"))
            absorb_dt(sp$pop, tbl) # no lookup_dt here as tbl not proper lu_tbl

            setnafill(sp$pop, "c", fill = 0, cols = "mu") # fix for prostate and breast cancer
            sp$pop[, (namprvl) := as.integer(dqrunif(.N) < mu)]
            sp$pop[, mu := NULL]


            # Then select individuals with higher risk to have higher
            # probability to be a prevalent case. I estimate weights based on
            # relevant RF for each disease. Note that RR and sampling weights
            # are equivalent here. The RR are for incident. One would expect
            # prevalent cases to be healthier because of survival of the fittest
            # and because some may have better RF control post diagnosis. For
            # that I will arbitrarily assume that the risk for prevalence is
            # half of that for incidence.

            if (length(private$rr) > 0L && any(sp$pop[[namprvl]] > 0L)) {
              # ncases is the number of prevalent cases expected in each stratum
              tt <- sp$pop[get(namprvl) > 0, .("ncases" = .N), by = eval(strata)]
              absorb_dt(sp$pop, tt, on = strata) # no lookup_dt as tt not a lu_tbl
              setnafill(sp$pop, "c", fill = 0, cols = "ncases")

              # Generate unique name using the relevant RR and lags

              # TODO below should apply only to strata with ncases > 0 for efficiency
              self$set_rr(sp, design_, forPARF = TRUE,
                          checkNAs = design_$sim_prm$logs)

              # Delete rr on self (i.e. for asthma)
              if (paste0(self$name, "_prvl_rr") %in% names(sp$pop))
                sp$pop[, (paste0(self$name, "_prvl_rr")) := NULL]

              nam <- grep("_rr$", names(sp$pop), value = TRUE) # necessary because above forPARF = TRUE

              sp$pop[, disease_wt := (Reduce(`*`, .SD)), .SDcols = nam]
              sp$pop[, (nam) := NULL]
              # adjust for prevalent risk half of incident risk to account for reverse causality and survival of the fittest

              sp$pop[, disease_wt := ((disease_wt - 1) * 0.5) + 1]

              ss <- sp$pop[ncases > 0,][, .("pid" = pid[sample_int_expj(unique(.N), unique(ncases), disease_wt)]), by = eval(strata)]
              sp$pop[, (namprvl) := 0L]
              sp$pop[ss, on = c("pid", "year"), (namprvl) := 1L]
              sp$pop[, c("ncases", "disease_wt") := NULL]
            }


            # set duration
            dqset.seed(private$seed, stream = sp$mc * 10 + 2L) # not mc_aggr
            tbl <- read_fst(private$filenams$dur, as.data.table = TRUE)
            col_nam <- setdiff(names(tbl), intersect(names(sp$pop), names(tbl)))
            tbl[, (namprvl) := 1L]
            lookup_dt(sp$pop, tbl, check_lookup_tbl_validity = FALSE)
            fn <- paste0("q", self$meta$diagnosis$duration_distr_backwards)
            sp$pop[get(namprvl) == 1L,
                   (namprvl) := 2L + do.call(fn, c(p = list(dqrunif(.N)), .SD)),
                   .SDcols = col_nam]
            sp$pop[, (col_nam) := NULL]

            if (!is.null(self$meta$mortality$cure) && self$meta$mortality$cure > 0L)
              sp$pop[get(namprvl) > self$meta$mortality$cure,
                     (namprvl) := self$meta$mortality$cure]

            sp$pop[, (namprvl) := carry_backward_decr(get(namprvl), pid_mrk)] # necessary for c++

          } # End if incidence type not 0 or 1


          # TODO this only makes sense when probability of diagnosis is 1
          namdgns <- paste0(self$name, "_dgns")
          set(sp$pop, NULL, namdgns, 0L)
          sp$pop[
            get(namprvl) > 0 & year >= design_$sim_prm$init_year & dqrunif(.N) < self$meta$diagnosis$probability,
            (namdgns) := get(namprvl)
          ]

          sp$pop[, (namdgns) := carry_backward_decr(get(namdgns), pid_mrk)]
        }

        invisible(self)
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

        if (!forPARF && length(private$rr) > 0) sp$store_risks(self$name)

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

        if (is.numeric(self$meta$incidence$type) && self$meta$incidence$type > 0L) {

          if (private$incd_colnam %in% names(sp$pop)) {
            stop(
              "A column named ", private$incd_colnam,
              " already exists in sp$pop. ",
              "Please delete it and run set_incd_prb() afterwards."
            )
          }

          # Get colnames in risk that end with _rr but exclude the influence by
          # diseases
          if (self$meta$incidence$type == 3) {
            # private$rr never NULL here but riskcolnam can be empty if disease
            # only influenced by other diseases but not exposures
            riskcolnam <- grep(
              paste0(
                "^((?!",
                paste(self$meta$incidence$influenced_by_disease_name,
                      collapse = "|"
                ),
                ").)*_rr$"
              ),
              names(sp$get_risks(self$name)),
              value = TRUE,
              perl = TRUE
            )
          } else { # private$rr may be NULL here but I will cover this case below
            riskcolnam <- grep("_rr$",
                               names(sp$get_risks(self$name)),
                               value = TRUE,
                               perl = TRUE
            )
          }

          if (self$meta$incidence$type == 1L) {
            if (length(riskcolnam) == 1L) {
              thresh <- as.numeric(sp$get_risks(self$name)[[riskcolnam]])
            }
            if (length(riskcolnam) > 1L && self$meta$incidence$aggregation == "any") {
              thresh <- as.numeric(sp$get_risks(self$name)[, do.call(pmax, .SD),
                                                           .SDcols = riskcolnam])
            }
            if (length(riskcolnam) > 1L && self$meta$incidence$aggregation == "all") {
              thresh <- as.numeric(sp$get_risks(self$name)[, Reduce(`*`, .SD),
                                                           .SDcols = riskcolnam])
            }

            set(sp$pop, NULL, private$incd_colnam, thresh)

          } else if (length(private$rr) > 0L) {
            # if incidence$type not 1 and at least 1 associated RF
            if (length(riskcolnam) > 0) {
              risk_product <-
                sp$get_risks(self$name)[, Reduce(`*`, .SD), .SDcols = riskcolnam]
            } else {
              risk_product <- 1
            }

            # Calibrate estimated incidence prbl to init year incidence
            tbl <- self$get_incd(design_$sim_prm$init_year
            )[between(age, design_$sim_prm$ageL,
                      design_$sim_prm$ageH)]
            lookup_dt(sp$pop, tbl)
            sp$pop[, rp := private$parf$p0 * sp$get_risks(self$name)[, Reduce(`*`, .SD),
                                                                     .SDcols = patterns("_rr$")]]
            setnafill(sp$pop, "c", 0, cols = c("rp", "mu"))
            # Above rp includes rr from diseases that risk_product doesn't have
            tbl <- sp$pop[year == design_$sim_prm$init_year &
                            get(paste0(self$name, "_prvl")) == 0L,
                          .(clbfctr = sum(mu)/sum(rp)),
                          keyby = .(sex)] # ,ethnicity, sha

            lookup_dt(sp$pop, tbl)
            setnafill(sp$pop, "c", 1, cols = "clbfctr")

            # End of calibration

            set(sp$pop, NULL, private$incd_colnam,
                clamp(private$parf$p0 * risk_product * sp$pop$clbfctr))
            # NOTE product above not expected to be equal to incidence because
            # p0 estimated using mean lags and RR, while each mc run samples
            # from their distribution.

            # setnames(sp$pop, "clbfctr", paste0(self$name, "_clbfctr"))
            # sp$pop[, (paste0(self$name, "_risk_product")) := risk_product]
            # sp$pop[, (paste0(self$name, "_p0")) := private$parf$p0]
            sp$pop[, c("mu", "rp", "clbfctr") := NULL]

            if (design_$sim_prm$export_PARF) {
              path <- file.path(design_$sim_prm$output_dir, "parf")
              filenam <- file.path(path, "parf.csv")
              if (!dir.exists(path)) {
                dir.create(path, showWarnings = FALSE, recursive = TRUE)
              }

              parf_dt <-
                cbind(sp$pop[, .(wt_immrtl, age, sex, dimd, ethnicity, sha, year)],
                      "parf" = private$parf$parf
                      )[year == design_$sim_prm$init_year]
              parf_dt <-
                parf_dt[!is.na(parf), .(parf = unique(parf), pop_size = sum(wt_immrtl)),
                       keyby = .(age, sex, dimd, ethnicity, sha)]
              parf_dt[, `:=`(disease = self$name, mc = sp$mc)] # not sp$mc_aggr
              # EXAMPLE parf[, weighted.mean(parf, pop_size), keyby = sex]
              fwrite_safe(parf_dt, filenam)
            } # End export PARF

          } else { # End of incident$type not 1 and no associated RF
            tbl <- self$get_incd(seq(design_$sim_prm$init_year,
                                     design_$sim_prm$init_year +
                                       design_$sim_prm$sim_horizon_max)
            )[between(age, design_$sim_prm$ageL, design_$sim_prm$ageH)]
            setnames(tbl, "mu", private$incd_colnam)
            lookup_dt(sp$pop, tbl, check_lookup_tbl_validity = FALSE)
          } # End if no associated RF

          setnafill(sp$pop, "c", 0, cols = private$incd_colnam)
        } # End if incident$type numeric
        invisible(self)
      },


      #' @description Set diagnosis probability in a new col in sp$pop.
      #' @param sp A synthetic population.
      #' @param design_ A design object with the simulation parameters.
      #' @return The invisible self for chaining.

      set_dgns_prb = function(sp, design_ = design) {
        if (is.numeric(self$meta$diagnosis$type) && self$meta$diagnosis$type > 0L) {
          if (!inherits(sp, "SynthPop")) {
            stop("Argument sp needs to be a SynthPop object.")
          }
          if (!inherits(design_, "Design")) {
            stop("Argument design_ needs to be a Design object.")
          }

          if (private$dgns_colnam %in% names(sp$pop)) {
            stop(
              "A column named ", private$dgns_colnam,
              " already exists in sp$pop. ",
              "Please delete it and run set_dgns_prb() afterwards."
            )
          }

          set(sp$pop, NULL, private$dgns_colnam, self$meta$diagnosis$probability)
        }

        invisible(self)
      },

      #' @description Set disease case fatality when relevant, in a new col in
      #'   sp$pop.
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
          if (private$mrtl_colnam2 %in% names(sp$pop))
            stop("Column ", private$mrtl_colnam2, " exists already in sp$pop.")

          ftlt <-
            self$get_ftlt(
              seq(
                design_$sim_prm$init_year,
                design_$sim_prm$init_year + design_$sim_prm$sim_horizon_max
              )
            )

          if (!"mu2" %in% names(ftlt))
            stop("mu2 need to be present in the ftlt file.")

          # Deal with mrtl1 if present as it is unaffected by the logic below
          if ("mu1" %in% names(ftlt)) {
            private$mrtl2flag <- TRUE
            nam <- paste0("prb_", self$name, "_mrtl", 1)
            if (nam %in% names(sp$pop))
              stop("Column ", nam, " exists already in sp$pop.")
            setnames(ftlt, "mu1", nam)
            lookup_dt(sp$pop, ftlt[, .SD, .SDcols = !"mu2"],
                      check_lookup_tbl_validity = FALSE)
            setnafill(sp$pop,
                      type = "const",
                      fill = 0,
                      cols = nam)
            ftlt[, (nam) := NULL]
          } else {
            private$mrtl2flag <- FALSE
          }

          # if (not apply_RR_to_mrtl2 and not depend on other diseases) or ----
          # length(private$rr) == 0L
          if ((!design_$sim_prm$apply_RR_to_mrtl2 &&
               !self$meta$mortality$type %in% 3:4) ||
              length(private$rr) == 0L) {
            setnames(ftlt, "mu2", private$mrtl_colnam2)
            lookup_dt(sp$pop, ftlt, check_lookup_tbl_validity = FALSE)
          } else {
            if (self$meta$mortality$type %in% 3:4) {
              # private$rr never NULL here but riskcolnam can be empty if disease
              # only influenced by other diseases but not exposures
              riskcolnam <- grep(
                paste0(
                  "^((?!",
                  paste(
                    self$meta$mortality$influenced_by_disease_name,
                    collapse = "|"
                  ),
                  ").)*_rr$"
                ),
                names(sp$get_risks(self$name)),
                value = TRUE,
                perl = TRUE
              )
            } else {
              # private$rr may be NULL here but I will cover this case
              # below
              riskcolnam <- grep("_rr$",
                                 names(sp$get_risks(self$name)),
                                 value = TRUE,
                                 perl = TRUE)
            }

            if (length(riskcolnam) > 0) {
              risk_product <-
                sp$get_risks(self$name)[, Reduce(`*`, .SD), .SDcols = riskcolnam]
            } else {
              risk_product <- 1
            }

            # Calibrate estimated fatality prbl to init year incidence
            tbl <- self$get_ftlt(design_$sim_prm$init_year
            )[between(age, design_$sim_prm$ageL,
                      design_$sim_prm$ageH)]
            if ("mu1" %in% names(tbl)) tbl[, mu1 := NULL]
            lookup_dt(sp$pop, tbl)
            sp$pop[, rp := private$parf$m0 * sp$get_risks(self$name)[, Reduce(`*`, .SD),
                                                                     .SDcols = patterns("_rr$")]]
            setnafill(sp$pop, "c", 0, cols = c("rp", "mu2"))
            # Above rp includes rr from diseases that risk_product doesn't have

            if (self$name == "nonmodelled") {
              tbl <- sp$pop[year == design_$sim_prm$init_year,
                            .(clbfctr = sum(mu2)/sum(rp)),
                            keyby = .(sex)]

            } else {
              tbl <- sp$pop[year == design_$sim_prm$init_year &
                              get(paste0(self$name, "_prvl")) > 0L,
                            .(clbfctr = sum(mu2)/sum(rp)),
                            keyby = .(sex)]

            }
            # NOTE the above excludes incident cases. Not appropriate when
            # private$mrtl2flag == FALSE
            lookup_dt(sp$pop, tbl)
            setnafill(sp$pop, "c", 1, cols = "clbfctr")

            # End of calibration

            set(sp$pop, NULL, private$mrtl_colnam2,
                clamp(private$parf$m0 * risk_product * sp$pop$clbfctr))

            # setnames(sp$pop, "clbfctr", paste0(self$name, "_clbfctr_mrtl"))
            # sp$pop[, (paste0(self$name, "_risk_product_mrtl")) := risk_product]
            # sp$pop[, (paste0(self$name, "_m0")) := private$parf$m0]
            sp$pop[, c("mu2", "rp", "clbfctr") := NULL]

          }

          setnafill(sp$pop, type = "const", fill = 0,
                    cols = private$mrtl_colnam2)

          } # End numeric mortality type

        invisible(self)
      },


      #' @description Deletes the PARF file from disk.
      #' @param invert deletes all other disease relevant PARF file except those
      #'   that are associated to the current settings.
      #' @return The invisible self for chaining.

      del_parf_file = function(invert = FALSE) {
        stopifnot(is.logical(invert))

        if (invert) {
          parf_filenam2 <- list.files(
            private$parf_dir,
            pattern = paste0("^PARF_", self$name, ".*\\.fst$"),
            full.names = TRUE
          )

          parf_filenam <-
            setdiff(parf_filenam2, private$parf_filenam)

          if (any(file.exists(parf_filenam))) file.remove(parf_filenam)

        } else { # if not invert

          if (file.exists(private$parf_filenam)) file.remove(private$parf_filenam)
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

      #' @description Deletes the stochastic effect files and indices from disk
      #'   for all relevant RR.
      #' @return The invisible self for chaining.
      del_stochastic_effect = function() {
        file.remove(private$filenam)
        file.remove(private$filenam_indx)
        lapply(private$rr, function(x) x$del_stochastic_effect)
        invisible(self)
      },



      #' @description Get the PARF by age/sex/dimd/ethnicity/sha.
      #' @param what Columns to return (p0, m0, or parf)
      #' @return A data.table with PARF.
      get_parf = function(what) {
        if (sum(dim(private$parf)) > 0L && !missing(what)) {
          private$parf[, ..what]
        } else {
          private$parf
        }
      },

      #' @description Get the PARF filename.
      #' @return A data.table with PARF.
      get_parf_filename = function() {
          private$parf_filenam
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

        # TODO add logic to track file changes

        for (i in seq_along(private$filenams)) {
          if (grepl("_indx.fst$", private$filenams[[i]])) next
          print(private$filenams[[i]])
          tbl <- read_fst(private$filenams[[i]], as.data.table = TRUE)
          val <- setdiff(names(tbl), names(sp$pop))
          com <- sort(intersect(names(tbl), names(sp$pop)))
          # Ensure year is always the first key if present
          com <- com[order(match(com, "year"))]


          # TODO add a property on yaml to recognise diseases apply to only one sex
          if (!"sex" %in% names(tbl) || uniqueN(tbl$sex) == 1L) {
            tbl2 <- copy(tbl)
            if (self$name == "breast_ca") {
              tbl[, sex := factor("women", levels = c("men", "women"))]
              tbl2[, sex := factor("men", levels = c("men", "women"))]
            }

            if (self$name == "prostate_ca") {
              tbl[, sex := factor("men", levels = c("men", "women"))]
              tbl2[, sex := factor("women", levels = c("men", "women"))]
            }

            for (j in val) set(tbl2, NULL, j, 0)

            if (!"sex" %in% com) com <- c(com, "sex")
            tbl <- rbind(tbl, tbl2)
          }

          for (j in com) {
            if (inherits(sp$pop[[j]], "integer") &&
                !inherits(tbl[[j]], "integer")) {
              tbl[, (j) := as.integer(get(j))]
            }
            if (inherits(sp$pop[[j]], "numeric") &&
                !inherits(tbl[[j]], "numeric")) {
              tbl[, (j) := as.numeric(get(j))]
            }
            if (inherits(sp$pop[[j]], "character") &&
                !inherits(tbl[[j]], "character")) {
              tbl[, (j) := as.character(get(j))]
            }
            if (inherits(sp$pop[[j]], "factor")) {
              # irrespective of class(j) to make sure that levels are the same
              # and in the right order.
              if (j == "dimd" &&
                  levels(tbl[[j]])[1] != "1 most deprived") {
                tbl[, (j) := factor(get(j),
                                    levels = as.character(10:1),
                                    labels = levels(sp$pop[[j]]))]
              } else {
                tbl[, (j) := factor(get(j), levels = levels(sp$pop[[j]]))]
              }
            }
          }

          if (grepl("_ftlt.fst$", private$filenams[[i]]) &&
            "mu" %in% names(tbl)) {
            setnames(tbl, "mu", "mu2")
          }
          if (self$name %in% c("breast_ca", "prostate_ca") && anyNA(tbl))
            setnafill(tbl, "c", fill = 0, cols = val)

          # if (!self$name %in% c("breast_ca", "prostate_ca") &&  anyNA(tbl))
          if (anyNA(tbl))
            stop("NAs in ", private$filenams[[i]])

          setkeyv(tbl, com)

          is_valid_lookup_tbl(tbl, com) # stops if not

          if (!identical(read_fst(private$filenams[[i]], as.data.table = TRUE), tbl))
            write_fst(tbl, private$filenams[[i]])
          } # End loop over relevant files
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
        if (is.null(out$incidence$can_recur))
          out$incidence <- within(out$incidence, rm("can_recur"))
        if (out$incidence$type == "Universal")
          out$incidence <- within(out$incidence, rm("prevalence", "probability"))
        if (out$incidence$type == "Type0")
          out$incidence <- within(out$incidence, rm("probability"))

        # TODO resolve influenced by disease automatically from
        # paste0(names(design_$sim_prm$diseases), "_prvl") %in% private$rr

        if (self$meta$incidence$type == 0L) {
          influenced_by_incd <- list()

          for (i in self$meta$incidence$influenced_by_disease_name) {
            influenced_by_incd[[paste0(i, "_prvl")]] <- list("lag" = 0L)
          } # end for loop over influenced_by_disease_name
          out[["incidence"]][["influenced_by"]] <- influenced_by_incd
        } # end if incidence type 0

        if (self$meta$incidence$type == 3L) {
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
            ),
            "mm_wt" = self$meta$diagnosis$mm_wt
          )

          if (!is.null(self$meta$diagnosis$duration_distr_forwards))
            out[["diagnosis"]][["duration_distr_forwards"]] <-
              read_yaml(self$meta$diagnosis$duration_distr_forwards)

          if (self$meta$diagnosis$type == 0L) {
            influenced_by_dgns <- list()

            for (i in self$meta$incidence$influenced_by_disease_name) {
              influenced_by_dgns[[paste0(i, "_dgns")]] <- list("lag" = 0L)
            } # end for loop over influenced_by_disease_name

            out[["diagnosis"]][["influenced_by"]] <- influenced_by_dgns

            out$diagnosis <- within(out$diagnosis, rm("probability"))
          }

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

          if (!is.null(self$meta$mortality$cure)) {
            out[["mortality"]][["cure"]] <- self$meta$mortality$cure
          }

          if (private$mrtl2flag) {
            out[["mortality"]][["probability1styear"]] <-
              paste0("prb_", self$name, "_mrtl1", scenario_suffix)
          }


          if (self$meta$mortality$type %in% 3:4) {
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
      mrtl_colnam2 = NA,
      incd_indx = data.table(NULL),
      prvl_indx = data.table(NULL),
      ftlt_indx = data.table(NULL),
      chksum = NA,
      parf_dir = NA,
      parf_filenam = NA,
      parf = data.table(NULL),
      rr = list(), # holds the list of relevant RR

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
      },

      # helper function to get the tree of dependencies to exposures
      # x is a disease name string i.e. x = "other_ca"
      # diseases_ is a list of disease objects
      # TODO test that if (!i %in% out$ds)) allows interdependency (i.e. chd
      # causes t2dm and t2dm causes chd). Current approach will lead to an
      # infinite loop
      get_xps_dependency_tree = function(x = self$name, dssl = diseases_) {
        tr <- sapply(dssl[[x]]$get_rr(), `[[`, "name")
        if (length(tr) == 0L) {
          out <- data.table(
            xpscol = character(),
            lag = numeric(),
            ds = character())
          return(out)
        }
        out <- data.table(
          xpscol = tr,
          lag = sapply(dssl[[x]]$get_rr(), `[[`, "lag"),
          ds = x)
        allds <-
          unique(
            c(
              dssl[[x]]$meta$incidence$influenced_by_disease_name,
              dssl[[x]]$meta$mortality$influenced_by_disease_name
            )
          )
        if (!is.null(allds)) {
          for (i in allds) {
            if (!i %in% out$ds)
              out <- unique(rbind(out, private$get_xps_dependency_tree(i, dssl)))
          }
        }
        return(out)
      },

      gen_sp_forPARF =
        function(mc_, ff, design_, diseases_) {

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


          xps_dep <- private$get_xps_dependency_tree(x = self$name, dssl = diseases_)
          xps_dep <- xps_dep[!grepl("_prvl$", xpscol)]
          setkey(xps_dep, xpscol)

          # NOTE an exposure may appear multiple times because it is a risk for
          # multiple conditions, with potentially different lags. However the
          # current approach allows only one lag time per exposure for
          # efficiency. This is a limitation of this approach which. To resolve
          # this for now, I will pick the highest based on the fact that this
          # is a disease chain. I.e BMI ->5y-> t2dm ->9y-> breast_ca ->1y-> other_ca
          # TODO implement a better solution for the above

          # Generate active days ----
          xps <- c("active_days", "met") # t2dm require MET
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
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
          xps <- "met"
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]

            ff[, met_curr_xps := as.integer(floor(active_days_curr_xps * (3L + qbinom(dqrunif(.N), 8, 3 / 11)) *
              (30 + qexp(dqrunif(.N), 1 / 7)) / 100))]
            ff[, year := year + lag]
          }

          # Generate fruit consumption (ZISICHEL) ----
          xps <- "fruit"
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]

            tbl <-
              read_fst("./inputs/exposure_distributions/frtpor_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
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
          if (xps %in% xps_dep$xpscol) {
            lag <- xps_dep[xps, max(lag)]
            ff[, year := year - lag]

            tbl <-
              read_fst("./inputs/exposure_distributions/vegpor_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            ff[, veg_curr_xps :=
              my_qDEL(rank_veg, mu, sigma, nu, n_cpu = design_$sim_prm$n_cpu) * 80L] # g/d
            ff[, (col_nam) := NULL]
            ff[, rank_veg := NULL]
            ff[, year := year + lag]
          }

          # Smoking ----
          xps <- c("smok_status", "smok_cig", "smok_packyrs", "ets")
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps_dep[xps[1:3], any(duplicated(ds))]) {
              stop("smok_status, smok_cig, & smok_packyrs cannot coexist as risk factors for a disease.")
            }

            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
            }

            if (xps[[2]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[2]], max(lag)]
            }

            if (xps[[3]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[3]], max(lag)]
            }

            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_status_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
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
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            set(ff, NULL, "smok_quit_yrs_curr_xps", 0L)
            ff[
              smok_status_curr_xps %in% 2:3,
              smok_quit_yrs_curr_xps := my_qDPO(rankstat_smok_quit_yrs, mu, sigma)
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
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            set(ff, NULL, "smok_dur_curr_xps", 0L)
            ff[smok_status_curr_xps %in% 2:3, smok_dur_curr_xps := my_qDPO(rankstat_smok_dur_ex, mu, sigma)]
            ff[, rankstat_smok_dur_ex := NULL]
            ff[, (col_nam) := NULL]

            # Assign smok_dur_curr when pid_mrk == true (the first year an individual enters the simulation)
            tbl <-
              read_fst("./inputs/exposure_distributions/smok_dur_curr_table.fst",
                as.data.table = TRUE
              )
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            ff[smok_status_curr_xps == 4, smok_dur_curr_xps := as.integer(round(qNBI(rankstat_smok_dur_curr, mu, sigma)))]
            ff[, rankstat_smok_dur_curr := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }

          # Assign smok_cig_curr when pid_mrk == true (the first year an individual enters the simulation)
          xps <- c("smok_cig", "smok_packyrs")
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
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
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
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
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
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

          if ("smok_packyrs" %in% xps_dep$xpscol)
            ff[, smok_packyrs_curr_xps := as.integer(round(smok_cig_curr_xps *
                                          smok_dur_curr_xps / 20))]

          # Generate ETS (BI) ----

          # Note at the moment this is independent of smoking prevalence TODO
          # calculate how many each smoker pollutes by year, SHA (not qimd) to
          # be used in scenarios. Ideally correct for mortality
          xps <- "ets"
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/ets_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            ff[, ets_curr_xps := as.integer(rank_ets < mu)]
            ff[, rank_ets := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }
          # Generate alcohol (ZINBI) ----
          xps <- "alcohol"
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/alcohol_table.fst",
                       as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            ff[, alcohol_curr_xps := as.integer(
              qZINBI(rank_alcohol, mu, sigma, nu))]
            ff[, rank_alcohol := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }

          # Generate BMI (BCPEo) ----
          xps <- "bmi"
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/bmi_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            ff[, bmi_curr_xps := my_qBCPEo(rank_bmi, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
            ff[, rank_bmi := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }
          # Generate SBP (BCPEo) ----
          xps <- "sbp"
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/sbp_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            ff[, sbp_curr_xps := my_qBCPEo(rank_sbp, mu, sigma, nu, tau, n_cpu = design_$sim_prm$n_cpu)]
            ff[, rank_sbp := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }
          # Generate tchol (BCT) ----
          xps <- c("tchol", "statin_px")
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
            } else {
              lag <- 0L
            }
            ff[, year := year - lag]
            tbl <-
              read_fst("./inputs/exposure_distributions/tchol_table.fst", as.data.table = TRUE)
            col_nam <-
              setdiff(names(tbl), intersect(names(ff), names(tbl)))
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
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
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            ff[, tchol_hdl_ratio := 1 / qGB1(rank_hdl, mu, sigma, nu, tau)]
            ff[, rank_hdl := NULL]
            ff[, (col_nam) := NULL]
            ff[, year := year + lag]
          }
          # Generate statins medication (BI) -----
          xps <- "statin_px"
          if (any(xps %in% xps_dep$xpscol)) {
            if (xps[[1]] %in% xps_dep$xpscol) {
              lag <- xps_dep[xps[[1]], max(lag)]
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
            lookup_dt(ff, tbl, check_lookup_tbl_validity = FALSE)
            ff[, statin_px_curr_xps := as.integer(rank_statin_px < mu)]
            ff[, rank_statin_px := NULL]
            ff[, (col_nam) := NULL]
            ff[, `:=`(tchol = NULL)]
            ff[, year := year + lag]
          }

          nam <- grep("rank", names(ff), value = TRUE)
          if (length(nam) > 0) ff[, (nam) := NULL]

          if ("age100" %in% names(ff)) {
            ff[, age := NULL]
            setnames(ff, "age100", "age")
          }

          invisible(ff)
        }


    ) # end of private
  )
