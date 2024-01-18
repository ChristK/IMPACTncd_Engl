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


#' R6 Class representing an exposure
#'
#' @description
#' An exposure object has the RR between an exposure and a disease.
#'
#' @details
#' To be completed...
#'
#' @export
Exposure <-
  R6::R6Class(
    classname = "Exposure",
    # cloneable = FALSE, # cloneable is necessary for multi-threading
    # public ------------------------------------------------------------------
    public = list(
      #' @field name The name of the exposure.
      name = NA_character_,
      #' @field outcome The name of the outcome the exposure influences.
      outcome = NA_character_,
      #' @field lag Disease lag.
      lag = NA_integer_,
      #' @field distribution The distribution to be used for Monte Carlo.
      #'   Currently lognormal and normal are supported.
      distribution = NA_character_,
      #' @field source Citation info for the effect size.
      source = NA_character_,
      #' @field notes Any notes regarding the exposure -> outcome relation.
      notes = NA_character_,


      # initialize ----
      #' @description Reads exposure parameter from file and creates a new exposure object..
      #' @param sRelativeRiskByPopulationSubsetForExposureFilePath string, path to .csvy file detailing relative risk (RR) by population subset (age, sex, maybe decile 'Index of Multiple Deprivation' DIMD) for a specific exposure. File header may contain other exposure parameters.
      #' @param design A design object with the simulation parameters.
      #' @return An `Exposure` object.
      #' @examples
      #' af_stroke$read_xps_prm("./inputs/RR/af_stroke.csvy", design)
      initialize = function(sRelativeRiskByPopulationSubsetForExposureFilePath, design) {

        sRelativeRiskByPopulationSubsetForExposureFilePath<- normalizePath(sRelativeRiskByPopulationSubsetForExposureFilePath, mustWork = TRUE)

        if (!inherits(design, "Design"))
          stop("Argument design needs to be a Design object.")

        if (file.exists(sRelativeRiskByPopulationSubsetForExposureFilePath)) {
          # TODO add some checks to ensure proper structure
          dtRelativeRiskByPopulationSubset<- fread(
            sRelativeRiskByPopulationSubsetForExposureFilePath,
            stringsAsFactors = TRUE, yaml = TRUE
            )

          # NOTE This needs manual update every time a new factor (other than
          # smok_status, ethnicity, dimd) appears as a new stratum in a RR file.
          # TODO add some checks and perhaps automate to some extent. Add
          # documentation to avoid silent errors
          if ("smok_status" %in% names(dtRelativeRiskByPopulationSubset))
            dtRelativeRiskByPopulationSubset[, smok_status :=
                               factor(smok_status, levels = 1:4)]
          if ("ethnicity" %in% names(dtRelativeRiskByPopulationSubset))
            dtRelativeRiskByPopulationSubset[, ethnicity :=
                     factor(
                       ethnicity,
                       levels = c(
                         "white", "indian", "pakistani", "bangladeshi",
                         "other asian", "black caribbean", "black african",
                         "chinese", "other"
                         )
                      )
                   ]

          if ("dimd" %in% names(dtRelativeRiskByPopulationSubset))
            dtRelativeRiskByPopulationSubset[, dimd :=
                     factor(
                       dimd,
                       levels = c(
                         "1 most deprived", 2:9, "10 least deprived"
                       )
                     )
            ]

          if ("qimd" %in% names(dtRelativeRiskByPopulationSubset))
            dtRelativeRiskByPopulationSubset[, qimd :=
                     factor(
                       qimd,
                       levels = c(
                         "1 most deprived", 2:4, "5 least deprived"
                       )
                     )
            ]


          dtRelativeRiskByPopulationSubset[, agegroup := relevel(agegroup, "<1")]

          nam <- setdiff(names(dtRelativeRiskByPopulationSubset), c("rr", "ci_rr"))
          # sex and agegroup always at the beginning
          nam <- nam[order(match(nam, c("sex", "agegroup")))]
          setkeyv(dtRelativeRiskByPopulationSubset, nam)

          metadata <- attr(dtRelativeRiskByPopulationSubset, "yaml_metadata")
          setattr(dtRelativeRiskByPopulationSubset, "yaml_metadata", NULL)
            # private$get_yaml_header(sRelativeRiskByPopulationSubsetForExposureFilePath, verbose = design$sim_prm$logs)
        } else {
          stop(
            "File does not exist for argument sRelativeRiskByPopulationSubsetForExposureFilePath path."
          )
        }

        # Used only for parf by xps
        if (length(design$sim_prm$ignore_xps) > 0L &&
            (identical(metadata$xps_name, design$sim_prm$ignore_xps) ||
            (metadata$xps_name == "met" && design$sim_prm$ignore_xps == "active_days") ||
            (grepl("^smok_", metadata$xps_name) && design$sim_prm$ignore_xps == "smoking")
            )) {
          set(dtRelativeRiskByPopulationSubset, NULL, "rr", 1)
          set(dtRelativeRiskByPopulationSubset, NULL, "ci_rr", 1)
        }

        # Validation
        stopifnot(
          c("xps_name", "outcome", "distribution", "lag") %in% names(metadata),
          metadata$lag >= 1 | (metadata$lag == 0L & metadata$incidence$type == 1L),
          metadata$lag <= design$sim_prm$maxlag,
          length(metadata$distribution) == 1L,
          metadata$distribution %in% c("lognormal", "normal")
        )

        self$name         <- metadata$xps_name
        self$outcome      <- metadata$outcome
        self$lag          <- metadata$lag
        self$distribution <- metadata$distribution
        self$source       <- metadata$source
        self$notes        <- metadata$notes

        private$nam_rr <- paste0(self$name, "_rr") # i.e. bmi_rr

        if ("apply_rr_extra_fn" %in% names(metadata) &&
            startsWith(metadata$apply_rr_extra_fn, "function(")) {
          # https://www.r-bloggers.com/2016/08/yaml-define-an-r-function/
          # NOTE check that string starts with function for an added layer
          # of security. In any case still very insecure.
          private$apply_rr_extra <- eval(str2lang(metadata$apply_rr_extra_fn))
        } else private$apply_rr_extra <- function(...) NULL

        dqRNGkind("pcg64")

        # NOTE that for the different dimensions of smoking per disease we need
        # the same lag, so quit-yrs applied with the correct lag. Hence self$name
        private$seed <- abs(digest2int(paste0(ifelse(
          grepl("^smok_", self$name), "smoking", self$name), # rename all smoking dimensions to smoking
          self$outcome), seed = 764529L))
        private$chksum <- digest(dtRelativeRiskByPopulationSubset) # NOTE not affected by metadata changes
        private$suffix <- paste0(self$name, "~", self$outcome)

        private$filedir <- file.path(getwd(), "simulation", "rr")
        if (!dir.exists(private$filedir)) dir.create(private$filedir,recursive=TRUE)
        private$filenam <- file.path(private$filedir,
                                     paste0("rr_", private$suffix, "_",
                                            private$chksum, "_l.fst"))
        private$filenam_indx <- file.path(private$filedir,
                                          paste0("rr_", private$suffix, "_",
                                                 private$chksum, "_indx.fst"))

        dqset.seed(private$seed, stream = 1) # stream = 1 for lags

        if (self$lag == 0L) { # Only allowed for incd type 1

           private$lag_mc <-
            rep(0L, times = design$sim_prm$iteration_n_max)

        } else {
          private$lag_mc <-
            1L + as.integer(qbinom(
              dqrunif(design$sim_prm$iteration_n_max),
              design$sim_prm$maxlag - 1L,
              (self$lag - 1) / (design$sim_prm$maxlag - 1L)
            )) # Note qbinom returns double not int
        }


        if ("ideal_xps_lvl_fn" %in% names(metadata) &&
            startsWith(metadata$ideal_xps_lvl_fn, "function(")) {
          # NOTE check that string starts with function for an added layer
          # of security. In any case still very insecure.
          foo <- eval(str2lang(metadata$ideal_xps_lvl_fn))
          private$ideal_xps_lvl_mc <- foo(design)
          rm(foo)
        }

        if (length(nam) == 2L) {
          tt <-
            CJ(age = (design$sim_prm$ageL - design$sim_prm$maxlag):design$sim_prm$ageH,
              sex = factor(c("men", "women")))
          setkeyv(tt, c("sex", "age"))
        } else if (length(nam) == 3L) {
          nam <- setdiff(nam, c("sex", "agegroup"))
          if (is.numeric(dtRelativeRiskByPopulationSubset[[nam]])) {
            t3 <- min((dtRelativeRiskByPopulationSubset[[nam]])):max((dtRelativeRiskByPopulationSubset[[nam]]))
            t4 <- sort(unique(dtRelativeRiskByPopulationSubset[[nam]]))
            # check if consecutive elements
            if (!isTRUE(all.equal(t3, t4))) interpolate <- TRUE

            tt <-
              CJ(age = (design$sim_prm$ageL - design$sim_prm$maxlag):design$sim_prm$ageH,
                sex = factor(c("men", "women")),
                V3 = t3)

          } else { # if not numeric
            tt <-
              CJ(age = (design$sim_prm$ageL - design$sim_prm$maxlag):design$sim_prm$ageH,
                sex = factor(c("men", "women")),
                V3 = unique(dtRelativeRiskByPopulationSubset[[nam]]))
          }

          setnames(tt, "V3", nam)
          setkeyv(tt, c("sex", "age", nam))
        } else { # length(nam) > 3
          stop("Only one additional column beyond sex and agegroup is currently supported.")
        }


        to_agegrp(tt, 5L, max(tt$age), "age", "agegroup")

        private$effect <- copy(dtRelativeRiskByPopulationSubset)
        private$metadata <- metadata
        private$xps_prm_file <- sRelativeRiskByPopulationSubsetForExposureFilePath

        private$input_rr <-
          setcolorder(dtRelativeRiskByPopulationSubset[tt, on = .NATURAL], c("age", "agegroup", "sex"))
        private$input_rr[, "ci_rr" := NULL] # NOTE agegroup necessary for gen_stochastic_RR
        if (exists("interpolate") && isTRUE(interpolate)) { # linear interpolation
          private$input_rr[, rr := approx(get(nam), rr, get(nam))$y,
                           by = .(age, sex)]
        }
        if ("smok_status" %in% names(private$input_rr))
          private$input_rr[, smok_status :=
                         factor(smok_status, levels = 1:4)]
        if ("ethnicity" %in% names(private$input_rr))
          private$input_rr[, ethnicity :=
                   factor(
                     ethnicity,
                     levels = c(
                       "white", "indian", "pakistani", "bangladeshi",
                       "other asian", "black caribbean", "black african",
                       "chinese", "other"
                     )
                   )
          ]
        invisible(self)
      },

      # gen_stochastic_effect ----
      #' @description Generates and write to disk the stochastic effect.
      #' @param design_ A design object with the simulation parameters.
      #' @param overwrite If TRUE overwrite the files. Else if files exist they
      #'   are not regenerated.
      #' @param smooth If true applies loess smoothing to the input relative risks
      #' @param ... Further arguments to be passed to `loess`, usually span =
      #'   0.7, degree = 1
      #' @return An `Exposure` object.
      #' @examples
      #' af_stroke <- Exposure$new("af", "stroke")
      #' af_stroke$read_xps_prm("./inputs/RR/af_stroke.csvy", design)
      #' af_stroke$gen_stochastic_effect(design, TRUE)
      gen_stochastic_effect = function(design_ = design,
        overwrite = FALSE, smooth, ...) {
        if (!inherits(design_, "Design"))
          stop("Argument design needs to be a Design object.")

        if (overwrite ||
            all(!overwrite, any(
              !file.exists(private$filenam),
              !file.exists(private$filenam_indx)
            )) || (max(read_fst(private$filenam, columns = "age")) <
                   design_$sim_prm$ageH)
            ) {
          if (design_$sim_prm$logs)
            message(paste0(self$name, "~", self$outcome, " simulate RR uncertainty."))

          self$del_stochastic_effect(invert = TRUE) # deletes previous versions

          stoch_effect <-
            private$generate_rr_l(private$effect,
              private$input_rr[, .SD, .SDcols = !c("rr")],
              design_$sim_prm$iteration_n_max,
              smooth,
              design_$sim_prm$logs, ...)
          if ("smok_status" %in% names(stoch_effect))
            stoch_effect[, smok_status :=
                           factor(smok_status, levels = 1:4)]
          if ("ethnicity" %in% names(stoch_effect))
            stoch_effect[, ethnicity :=
                               factor(
                                 ethnicity,
                                 levels = c(
                                   "white", "indian", "pakistani", "bangladeshi",
                                   "other asian", "black caribbean", "black african",
                                   "chinese", "other"
                                 )
                               )
            ]

          write_fst(stoch_effect, private$filenam, 100)
          # create a table with row numbers for each mc
          stoch_effect[, rn := .I]
          tt <-
            stoch_effect[, .(from = min(rn), to = max(rn)), keyby = mc]
          write_fst(tt, private$filenam_indx, 100L)

        } else {
          if (design_$sim_prm$logs)
            message(paste0(self$name, "~", self$outcome, " RR files exist already."))
        }

        invisible(self)
      },

      # del_stochastic_effect ----
      #' @description Deletes the stochastic effect file and index from disk.
      #' @param invert If TRUE keeps the  file with the current checksum and deletes
      #'   all other files for this exposure~outcome.
      #' @return The invisible self for chaining.
      del_stochastic_effect = function(invert = FALSE) {
        stopifnot(is.logical(invert))

        if (invert) {
          fl <- list.files(
            private$filedir,
            pattern = paste0("^rr_", private$suffix, ".*\\.fst$"),
            full.names = TRUE
          )

          fl <- setdiff(fl, c(private$filenam, private$filenam_indx))

          file.remove(fl)

        } else {

          file.remove(c(private$filenam, private$filenam_indx))

        }
        invisible(self)
      },

      # get_rr ----
      #' @description Get relative risks from disk.
      #' @param mc An integer that signifies the Monte Carlo iteration.
      #' @param design_ A design object with the simulation parameters.
      #' @param drop If `TRUE` returns a scalar numeric of the RR if the RR is
      #'   common for all age groups and both sexes. Otherwise, a data.table.
      #' @param plot_rr If TRUE, plots the relative risk
      #' @return A data.table with the stochastic relative risks, if stochastic
      #'   = TRUE; else, the deterministic relative risks.
      get_rr =
        function(mc,
          design_ = design,
          drop = FALSE,
          plot_rr = FALSE) {
          if (!inherits(design_, "Design"))
            stop("Argument design_ needs to be a Design object.")

          stopifnot(length(mc) == 1L, between(mc, 1L, design_$sim_prm$iteration_n_max))
          mc <- as.integer(ceiling(mc))
          if (identical(mc, private$cache_mc)) {
            out <- copy(private$cache) # copy for safety
          } else {
            if (design_$sim_prm$stochastic) {
              indx <-
                read_fst(
                  private$filenam_indx,
                  from = mc,
                  to = mc,
                  as.data.table = TRUE
                )
              out <-
                read_fst(
                  private$filenam,
                  from = indx$from,
                  to = indx$to,
                  as.data.table = TRUE
                )
              out[, mc := NULL]
            } else {
              # if deterministic
              out <- copy(private$input_rr)
              setnames(out, "rr", private$nam_rr)
            }

            private$cache <- copy(out) # copy for safety
            private$cache_mc <- mc
          }


          if (plot_rr) {
            out[, plot(
              age,
              get(private$nam_rr),
              col = sex,
              pch = ".",
              ylab = "RR",
              main = (private$nam_rr)
            )]
            points(private$input_rr$age,
              private$input_rr$rr,
              col = private$input_rr$sex)
          }


          if (drop && uniqueN(out[[private$nam_rr]]) == 1L)
            out <- unique(out[[private$nam_rr]])


          out[]
        },

      # clear_cache ----
      #' @description Clear the cache for get_rr.
      #' @return The `Exposure` object.
      clear_cache = function() {
        private$cache <- NA
        private$cache_mc <- NA
        message("Cache was cleared")

        invisible(self)
      },

      # xps_to_rr ----
      #' @description Apply the RR in a new column in sp$pop based on the
      #'   exposure level. I the case of smok_quit_yrs it modifies the risk of
      #'   smok_cig or packyrs and no new column is created.
      #' @param sp A synthetic population.
      #' @param design_ A design object with the simulation parameters.
      #' @param checkNAs If `TRUE`, prints the table of NAs before they get
      #'   overwritten with 1. Note that for some exposures, NAs are expected
      #'   for certain levels of exposure (i.e. for active days).
      #' @param forPARF Set TRUE when applied on the specialised forPARF
      #'   SynthPop
      #' @return The invisible self for chaining.

      xps_to_rr = function(sp, design_, checkNAs = design_$sim_prm$logs,
                            forPARF = FALSE) {
          if (!inherits(design_, "Design"))
            stop("Argument design_ needs to be a Design object.")
          if (!inherits(sp, "SynthPop"))
            stop("Argument sp needs to be a SynthPop object.")
          if (private$nam_rr %in% names(sp$pop))
            stop(private$nam_rr, " already present in the data.")

          xps_tolag <- paste0(self$name, "_curr_xps")

          if (self$name %in% names(sp$pop)) {
            # To prevent overwriting t2dm_prvl, af_prvl etc.
            if (!xps_tolag %in% names(sp$pop)) {
              set(sp$pop, NULL, xps_tolag, 0L) # Assume only missing for diseases
              sp$pop[get(self$name) > 0, (xps_tolag) := 1L]
            }
            setnames(sp$pop, self$name, paste0(self$name, "____"))
          }

          if (forPARF) {
            set(sp$pop, NULL, self$name, sp$pop[[xps_tolag]])
            lookup_dt(sp$pop, self$get_input_rr(),
                      check_lookup_tbl_validity = design_$sim_prm$logs)

          } else {
            if (inherits(sp$pop[[xps_tolag]], "numeric")) {
              rw <- 0
            } else if (inherits(sp$pop[[xps_tolag]], "integer")) {
              rw <- 0L
            } else if (inherits(sp$pop[[xps_tolag]], "factor")) {
              rw <- 1L # The first level
            } else {
              stop("Only numerics, integers, and factors are supported")
            }
            set(sp$pop, NULL, self$name, # column without _curr_xps is lagged
                shift_bypid(sp$pop[[xps_tolag]], self$get_lag(sp$mc_aggr), sp$pop$pid, rw))
            # setnafill(sp$pop, "nocb", cols = self$name)
            lookup_dt(sp$pop, self$get_rr(sp$mc_aggr, design_, drop = FALSE),
                      check_lookup_tbl_validity = design_$sim_prm$logs)
          }

          private$apply_rr_extra(sp)

          if (checkNAs && private$nam_rr %in% names(sp$pop)) {
            # NOTE in case of smok_quit_yrs the risk column is deleted from the
            # rr_extra_fn()
            print(self$name)
            print(sp$pop[!is.na(get(self$name)) & is.na(get(private$nam_rr)) &
                           age >= design_$sim_prm$ageL &
                           year >= design_$sim_prm$init_year,
                           table(get(self$name), useNA = "always")])
          }

          if (private$nam_rr %in% names(sp$pop) &&
              is.numeric(sp$pop[[private$nam_rr]]))
            setnafill(sp$pop, type = "const", 1, cols = private$nam_rr)

          sp$pop[, (self$name) := NULL]

          if (paste0(self$name, "____") %in% names(sp$pop)) {
            # To prevent overwriting t2dm_prvl
            sp$pop[, (xps_tolag) := NULL]
            setnames(sp$pop, paste0(self$name, "____"), self$name)
          }

        return(invisible(self))
      },

      # validate_rr ----
      #' @description Get relative risks from disk.
      #' @return A plot with the input and stochastic relative risks.
      validate_rr =
        function() {
          layout(matrix(1:2))
          out <-
            read_fst(private$filenam,
              as.data.table = TRUE)
          out[, mc := NULL]

          layout(matrix(1:2))
          out[sex == "men", plot(
            age,
            get(private$nam_rr),
            col = rgb(red = 0, green = 0, blue = 1, alpha = 0.05),
            pch = ".",
            ylab = "RR",
            main = (paste0(self$name, "_rr (men)"))
          )]
          private$input_rr[sex == "men", points(age, rr)]
          out[sex == "women", plot(
            age,
            get(private$nam_rr),
            col = rgb(red = 1, green = 0, blue = 0, alpha = 0.05),
            pch = ".",
            ylab = "RR",
            main = (paste0(self$name, "_rr (women)"))
          )]
          private$input_rr[sex == "women", points(age, rr)]
          layout(matrix(1:1))
          # points(private$input_rr$age,
          #   private$input_rr$rr,
          #   col = private$input_rr$sex)
      },

      # get_input_rr ----
      #' @description Get original input relative risks by age and sex.
      #' @return A copied data.table with the original relative risks.
      get_input_rr = function() {
        out <- copy(private$input_rr)
        out[, agegroup := NULL]
        setnames(out, "rr", private$nam_rr)
        private$cache <- copy(out)
        private$cache_mc <- "forPARF"
        out
        },

      # get_metadata ----
      #' @description Get original metadata.
      #' @return A list with the original metadata.
      get_metadata = function() {private$metadata},

      # get_seed ----
      #' @description Get seed for RNG.
      #' @return A seed for the RNG that is produced by the digest of exposure
      #'   name and outcome.
      get_seed = function() {private$seed},

      # write_xps_tmplte_file ----
      #' @description Writes a template exposure file to disk.
      #' @param file_path Path including file name and .csvy extension to write the
      #'   file with placeholder exposure parameters.
      #' @return The `Exposure` object.
      write_xps_tmplte_file =
        function(file_path = "./inputs/RR/template.csvy") {
          file_path <- normalizePath(file_path, mustWork = FALSE)
          metadata <- list(
            "xps_name"     = "Exposure name",
            "outcome"       = "Disease name",
            "distribution"  = "lognormal or normal",
            "source"        = "Citation for the effectsize",
            "notes"         = "Any additional notes"
          )
          effect <- CJ(
            agegroup = factor(agegrp_name(min_age = 0, max_age = 99, grp_width = 5)),
            sex = c("men", "women"),
            rr = 1,
            ci_rr = 1
          )
          setkey(effect, sex, agegroup)

          private$write_xps_prm_file(effect, metadata, file_path)
          invisible(self)
        },

      # convert_from_old_format ----
      #' @description Convert the old format .csv to the new format .csvy.
      #' @param old_file Path to the old format .csv file with the RR.
      #' @param metadata List with the metadata information.
      #' @param file_path Path including file name and .csvy extension to write the
      #'   new file.
      #' @param estimates Only used when old_file is missing. A vector of length
      #'   2 with the 1st element being the point estimate and the 2nd one of
      #'   the CI.
      #' @param second_estimate_is_p If `TRUE` and the `estimates` not `NULL` then the
      #'   second element of `estimates` is interpreted like a p-value for ratio
      #'   rather than a CI.
      #' @return The `Exposure` object.
      convert_from_old_format =
        function(old_file, metadata, file_path, estimates = NULL, second_estimate_is_p = FALSE) {
          file_path <- normalizePath(file_path, mustWork = FALSE)
          colby <- character(0)

          if (!missing(old_file)) {
            old_file <- normalizePath(old_file, mustWork = TRUE)

            effect <-
              fread(old_file,
                stringsAsFactors = TRUE)
            if ("agegroup" %in% names(effect)) {
              colby <- c(colby, "agegroup")
              setkeyv(effect, "agegroup")
            }
            if ("sex" %in% names(effect)) {
              effect[, sex := factor(sex, 1:2, c("men", "women"))]
              colby <- c(colby, "sex")
            }
          } else if (second_estimate_is_p) {
            effect <- CJ(
              agegroup = factor(agegrp_name(min_age = 30, max_age = 89, grp_width = 5)),
              sex = factor(c("men", "women")),
              mean.rr = estimates[[1]],
              ci.rr = private$ci_from_p_forratio(estimates[[1]], estimates[[2]])[[2]]
            )
            colby <- c(colby, "agegroup", "sex")
          } else {
            effect <- CJ(
              agegroup = factor(agegrp_name(min_age = 30, max_age = 89, grp_width = 5)),
              sex = factor(c("men", "women")),
              mean.rr = estimates[[1]],
              ci.rr = estimates[[2]]
            )
            colby <- c(colby, "agegroup", "sex")
          }


          nam <- names(effect)
          nam <- setdiff(nam, c(colby, "mean.rr", "ci.rr"))

          if (length(nam) == 0L) {
            full <- CJ(
              agegroup = factor(agegrp_name(min_age = 0, max_age = 99, grp_width = 5)),
              sex = factor(c("men", "women")),
              rr = 1,
              ci_rr = 1
            )
          } else if (length(nam) == 1L) {
            full <- CJ(
              agegroup = factor(agegrp_name(min_age = 0, max_age = 99, grp_width = 5)),
              sex = factor(c("men", "women")),
              unique(effect[[nam]]),
              rr = 1,
              ci_rr = 1
            )
            setnames(full, "V3", nam)
            colby <- c(colby, nam)
          } else {
            stop("The structure of the old file has more cols than expected")
          }

          full[effect, on = colby, `:=` (
            rr = i.mean.rr,
            ci_rr = i.ci.rr
          )]

          setkeyv(full, c(nam, "sex", "agegroup"))

          private$write_xps_prm_file(full, metadata, file_path)
          invisible(self)
        },

      # get_lag ----
      #' @description Get exposure lag.
      #' @param mc_ A vector of Monte Carlo iterations. If missing or 0 return median (= mean).
      #' @return An integer vector with exposure~disease lag.
      get_lag = function(mc_) {
        if (missing(mc_) || mc_ == 0) {
          self$lag
        } else {
          private$lag_mc[mc_]
        }
      },

      # get_ideal_xps_lvl ----
      #' @description Get ideal exposure level.
      #' @param mc_ A vector of Monte Carlo iterations. If missing or 0 return
      #'   the user input lag.
      #' @return An integer vector with exposure~disease lag.
      get_ideal_xps_lvl = function(mc_) {
        if (missing(mc_) || mc_ == 0) {
          mean(private$ideal_xps_lvl_mc)
        } else {
          private$ideal_xps_lvl_mc[mc_]
        }
      },

      # get_name ----
      #' @description Get name of the object.
      #' @return An string.
      get_name = function() {private$suffix},

      # print ----
      #' @description Print the simulation parameters.
      #' @return The `Exposure` object.
      print = function() {
        print(paste0("Exposure name:      ", self$name))
        print(paste0("Outcome influenced: ", self$outcome))
        print(paste0("Distribution:       ", self$distribution))
        print(paste0("Median lag:         ", self$lag))
        print(paste0("Source:             ", self$source))
        print(paste0("Notes:              ", self$notes))

        invisible(self)
      }
    ), # end of public

    # private ------------------------------------------------------------------
     private = list(
       xps_prm_file = NA_character_,
       seed = NA_integer_,
       lag_mc = NA_integer_,
       input_rr = NA,
       effect = NA,
       metadata = NA,
       apply_rr_extra = NA,
       # https://www.r-bloggers.com/2016/08/yaml-define-an-r-function/
       ideal_xps_lvl_mc = 0,
       suffix = NA_character_,
       filenam = NA_character_,
       filenam_indx = NA_character_,
       filedir = NA_character_,
       cache = NA,
       cache_mc = NA,
       nam_rr = NA,
       chksum = NA,

       # write_xps_prm_file ----
       # @param dt A data table containing synthetic population data.
       # @param metadata A list or data structure containing metadata information.
       # @param file_path The path to the file where the exposure parameter file will be written.
       #
       # @details
       # The function creates an exposure file by combining metadata and synthetic population data, writes it to the specified file path, and appends the data using the fwrite function.
       #
       write_xps_prm_file = function(dt, metadata, file_path) {
         y <- paste0("---\n", yaml::as.yaml(metadata), "---\n")
         con <- textConnection(y)
         on.exit(close(con))
         m <- readLines(con)
         y <- paste0("#", m[-length(m)], collapse = "\n")
         y <- c(y, "\n")
         cat(y, file = file_path)
         fwrite(x = dt, file = file_path, append = TRUE, col.names = TRUE)
       },

       # stochRRtabl ----
       # @param m Central estimate of the relative risk.
       # @param ci Confidence interval of the relative risk.
       # @param distribution The distribution type to use for stochastic calculation, either "lognormal" or "normal".
       #
       # @return A stochastic relative risk value.
       #
       # @details
       # The function calculates stochastic relative risks by sampling from a specified distribution (lognormal or normal) using the central estimate and confidence interval.
       #
       stochRRtabl = # need to run by id
         function(m, ci, distribution = c("lognormal", "normal")) {
           distribution <- match.arg(distribution)
           kk <- dqrunif(1)
           dnm <- qnorm(0.975) # ~1.96
           if (distribution == "lognormal") {
             rr <- exp(qnorm(kk, log(m), abs(log(m) - log(ci)) / dnm))
           } else if (distribution == "normal") {
             rr <- qnorm(kk, m, abs(m - ci) / dnm)
           } else {
             stop("Only lognormal or normal are allowed.")
           }
           # rr[!is.finite(rr)] <- 1 # fix for rtruncnorm above
           return(rr)
         },

       # generate_rr_l ----
       # @param dt A data table containing age-grouped information.
       # @param tt A data table containing age information in years.
       # @param mc_max Maximum number of Monte Carlo iteration
       # @param smooth Logical, indicating whether smoothing should be applied to the relative risks.
       # @param do_checks Logical, indicating whether to perform additional checks on the generated relative risks.
       # @param ... Additional arguments to be passed to the loess smoothing function.
       #
       # @return A data table with generated relative risks based on linear interpolation and optional smoothing.
       #
       # @details
       # The function combines linear interpolation for intermediate exposure levels with optional smoothing to generate relative risks. It also performs checks and constraints on the generated values.
       #
       generate_rr_l =
         function(dt, tt, mc_max, smooth, do_checks = FALSE, ...) {
           # dt -> agegrouped, tt -> age in years
           dt <- tt[dt, on = .NATURAL, nomatch = NULL, mult = "first"
             ][, age := NULL]
         # uses no smoothing
         if (dt[, all(rr >= 1)]) {
           constrain <- "above_1"
         } else if (dt[, all(rr <= 1)]) {
           constrain <- "below_1"
         } else constrain <- "mix"
         colnam <- private$nam_rr
         colby <- c("mc", "sex")

         dt <- clone_dt(dt, mc_max, "mc")
         dqRNGkind("pcg64")
         dqset.seed(private$seed, stream = NULL)
         dt[, (colnam) := private$stochRRtabl(rr, ci_rr, self$distribution), keyby = mc]
         dt[, c("rr", "ci_rr") := NULL]
         dt <- tt[dt, on = .NATURAL, allow.cartesian = TRUE]

         dt[, agegroup := NULL]

         nam <- setdiff(names(dt), colnam)
         # mc, sex, and agegroup always at the beginning
         nam <- nam[order(match(nam, c("mc", "sex", "age")))]
         setkeyv(dt, nam)
         setcolorder(dt, c("mc", "age", "sex"))

         # linear interpolation for intermediate exposure levels
         nam <- setdiff(nam, c("mc", "sex", "age"))
         colby <- c(colby, nam)

         if (length(nam) == 1L && is.numeric(dt[[nam]])) {

           t3 <- min((dt[[nam]])):max((dt[[nam]]))
           t4 <- sort(unique(dt[[nam]]))
           # check if consecutive elements and do interpolation if not
           if (!isTRUE(all.equal(t3, t4))) { # non consecutive
             tt <-
               CJ(mc = 1:mc_max,
                 age = min(tt$age):max(tt$age),
                 sex = factor(c("men", "women")),
                 V3 = t3)

             replace_from_table(
               tt,
               colname = "V3",
               from = t3,
               to = rep(t4, times = c(diff(t4), 1L)),
               newcolname = nam
             )

             dt <- tt[dt, on = .NATURAL]
             dt[V3 != get(nam), (colnam) := NA_real_]
             dt[, (colnam) := approx(get(nam), get(colnam), V3)$y, by = .(mc, age, sex)]
             dt[, (nam) := NULL]
            setnames(dt, "V3", nam)
           }
         }


         if (smooth) {
           dt[, (colnam) := predict(loess(as.formula(paste0(colnam, " ~ age")),
             .SD, ...)), by = eval(colby)]
         }
         if (constrain == "above_1") dt[get(colnam) < 1, (colnam) := 1]
         if (constrain == "below_1") dt[get(colnam) > 1, (colnam) := 1]


         if (do_checks) {
           if (dt[, anyNA(get(colnam))]) {
             warning("RR has NAs")
           } else if (dt[, all(get(colnam) >= 1)]) {
             message("All RR above 1")
           } else if (dt[, all(get(colnam) <= 1)]) {
             message("All RR below 1")
           } else
             message("RR crosses 1")
         }
         dt
       },

       # ci_from_p_forratio ----
       # calculate 95% CI of a ratio from a p value
       # @param mean_rr Mean relative risk.
       # @param p p-value.
       #
       # @return A list containing the lower and upper bounds of the confidence interval.
       #
       # @details
       # The function computes the confidence interval for a ratio using a mean relative risk and a p-value.
       #
       ci_from_p_forratio = function(mean_rr, p) {
         z <- -0.862 + sqrt(0.743 - 2.404 * log(p))
         est <- log(mean_rr)
         se <- abs(est/z)
         dnm <- qnorm(0.975, 0, 1) # ~1.96
         ci <- list()
         ci$lower_ci <- exp(est - dnm * se)
         ci$upper_ci <- exp(est + dnm * se)
         ci
       }

    ) # end of private
  )
