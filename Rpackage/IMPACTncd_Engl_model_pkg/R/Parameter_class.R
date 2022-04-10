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


#' R6 Class representing a parameter
#'
#' @description
#' A parameter object has a stochastic simulation parameter.
#'
#' @details
#' To be completed...
#'
#' @export
Parameter <- # TODO Needs extensive refinement
  R6::R6Class(
    classname = "Parameter",
    # cloneable = FALSE, # cloneable is necessary for multi-threading
    # public ------------------------------------------------------------------
    public = list(
      #' @field name The name of the parameter.
      name = NA_character_,
      #' @field outcome The name of the outcome the parameter influences.
      outcome = NA_character_,
      #' @field distribution The distribution to be used for Monte Carlo.
      #'   Currently lognormal and normal are supported.
      distribution = NA_character_,
      #' @field source Citation info for the effect size.
      source = NA_character_,
      #' @field notes Any notes regarding the parameter -> outcome relation.
      notes = NA_character_,

      #' @description Create a new parameter object.
      #' @param exps_name A string with parameter name.
      #' @param outcome A string with outcome name.
      #' @param distribution A string with distribution for Monte Carlo (lognormal or normal).
      #' @param source A string with source info.
      #' @param notes A string with any notes.
      #' @return A new `Parameter` object.

      initialize = function(exps_name, outcome, distribution = c("lognormal", "normal"),
        source = NA_character_, notes = NA_character_) {
        if (!(is.character(exps_name) && is.character(outcome)))
          stop("Both arguments need to be strings.")
        self$name <- exps_name
        self$outcome <- outcome
        self$distribution <- distribution[[1]]
        self$source <- source
        self$notes <- notes

        private$seed <- abs(digest2int(paste0(exps_name, outcome), seed = 764529L))
        private$suffix <- paste0(self$name, "~", self$outcome)
        private$filenam <- file.path(getwd(), "simulation", "parameters", paste0(private$suffix, "_prmtr_l.fst"))
        private$filenam_indx <- file.path(getwd(), "simulation", "parameters", paste0(private$suffix, "_prmtr_indx.fst"))

        invisible(self)
      },

      #' @description Reads parameter parameter from file.
      #' @param exps_prm A path to a csvy file with the parameter parameters.
      #' @param design_ A design object with the simulation parameters.
      #' @return A `Parameter` object.

      read_exps_prm = function(exps_prm, design_ = design) {

        exps_prm <- normalizePath(exps_prm, mustWork = TRUE)

        if (!inherits(design_, "Design"))
          stop("Argument design needs to be a Design object.")

        if (file.exists(exps_prm)) {
          # TODO add some checks to ensure proper structure
          effect <- fread(
            exps_prm,
            stringsAsFactors = TRUE
            )

          effect[, agegroup := relevel(agegroup, "<1")]

          nam <- setdiff(names(effect), c("rr", "ci_rr"))
          # sex and agegroup always at the beginning
          nam <- nam[order(match(nam, c("sex", "agegroup")))]
          setkeyv(effect, nam)

          metadata <-
            private$get_yaml_header(exps_prm, verbose = design_$sim_prm$logs)
        } else {
          stop(
            "File does not exist for argument exps_prm path."
          )
        }

        # Validation
        stopifnot(
          c(
            "exps_name"           ,
            "outcome"             ,
            "distribution"
          ) %in% names(metadata),
          all.equal(self$name, metadata$exps_name),
          all.equal(self$outcome, metadata$outcome)
          )

        if (length(nam) == 2L) {
          tt <-
            CJ(age = design_$sim_prm$ageL:design_$sim_prm$ageH,
              sex = factor(c("men", "women")))
          setkeyv(tt, c("sex", "age"))
        } else if (length(nam) == 3L) {
          nam <- setdiff(nam, c("sex", "agegroup"))
          if (is.numeric(effect[[nam]])) {
            t3 <- min((effect[[nam]])):max((effect[[nam]]))
            t4 <- sort(unique(effect[[nam]]))
            # check if consecutive elements
            if (!isTRUE(all.equal(t3, t4))) interpolate <- TRUE

            tt <-
              CJ(age = design_$sim_prm$ageL:design_$sim_prm$ageH,
                sex = factor(c("men", "women")),
                V3 = t3)

          } else { # if not numeric
            tt <-
              CJ(age = design_$sim_prm$ageL:design_$sim_prm$ageH,
                sex = factor(c("men", "women")),
                V3 = unique(effect[[nam]]))
          }

          setnames(tt, "V3", nam)
          setkeyv(tt, c("sex", "age", nam))
        } else { # length(nam) > 3
          stop("Only one additional column beyond sex and agegroup is currently supported.")
        }


        to_agegrp(tt, 5L, max(tt$age), "age", "agegroup")

        private$effect <- copy(effect)
        private$metadata <- metadata
        private$exps_prm_file <- exps_prm

        private$input_rr <-
          setcolorder(effect[tt, on = .NATURAL], c("age", "agegroup", "sex"))
        private$input_rr[, c("ci_rr") := NULL]
        if (exists("interpolate") && isTRUE(interpolate)) { # linear interpolation
          private$input_rr[, rr := approx(get(nam), rr, get(nam))$y, by = .(age, sex)]
        }

        self$distribution <- metadata$distribution
        self$source       <- metadata$source
        self$notes        <- metadata$notes

        invisible(self)
      },


      #' @description Generates and write to disk the stochastic effect.
      #' @param design_ A design object with the simulation parameters.
      #' @param overwrite If TRUE overwrite the files. Else if files exist they
      #'   are not regenerated.
      #' @param smooth If true applies loess smoothing to the input relative risks
      #' @param ... Further arguments to be passed to `loess`, usually span =
      #'   0.7, degree = 1
      #' @return A `Parameter` object.

      gen_stochastic_effect = function(design_ = design,
        overwrite = FALSE, smooth, ...) {
        if (!inherits(design_, "Design"))
          stop("Argument design needs to be a Design object.")

        if (overwrite ||
            all(!overwrite, any(
              !file.exists(private$filenam),
              !file.exists(private$filenam_indx)
            ))) {
          stoch_effect <-
            private$generate_rr_l(private$effect,
              private$input_rr[, .SD, .SDcols = !c("rr")],
              design_$sim_prm$iteration_n_max,
              smooth,
              design_$sim_prm$logs, ...)

          write_fst(stoch_effect, private$filenam, 100)
          # create a table with row numbers for each mc
          stoch_effect[, rn := .I]
          tt <-
            stoch_effect[, .(from = min(rn), to = max(rn)), keyby = mc]
          write_fst(tt, private$filenam_indx, 100L)

        } else {
          message(paste0(self$name, "~", self$outcome, " RR files exist already."))
        }
        invisible(self)
      },

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
          drop = TRUE,
          plot_rr = FALSE) {
          if (!inherits(design_, "Design"))
            stop("Argument design_ needs to be a Design object.")

          stopifnot(between(mc, 1, design_$sim_prm$iteration_n_max))
          mc <- floor(mc)
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
              private$input_rr[, c("agegroup") := NULL]
              setnames(out, "rr", paste0(self$name, "_rr"))
            }

            if ("smok_status" %in% names(out))
              out[, smok_status := factor(smok_status, levels = 1:4)]

            private$cache <- copy(out) # copy for safety
            private$cache_mc <- mc
          }


          if (plot_rr) {
            out[, plot(
              age,
              get(paste0(self$name, "_rr")),
              col = sex,
              pch = ".",
              ylab = "RR",
              main = (paste0(self$name, "_rr"))
            )]
            points(private$input_rr$age,
              private$input_rr$rr,
              col = private$input_rr$sex)
          }


          if (drop && uniqueN(out[[paste0(self$name, "_rr")]]) == 1L)
            out <- unique(out[[paste0(self$name, "_rr")]])


          out[]
        },

      #' @description Clear the cache for get_rr.
      #' @return The `Parameter` object.
      clear_cache = function() {
        private$cache <- NA
        private$cache_mc <- NA
        message("Cache was cleared")

        invisible(self)
      },


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
            get(paste0(self$name, "_rr")),
            col = rgb(red = 0, green = 0, blue = 1, alpha = 0.05),
            pch = ".",
            ylab = "RR",
            main = (paste0(self$name, "_rr (men)"))
          )]
          private$input_rr[sex == "men", points(age, rr)]
          out[sex == "women", plot(
            age,
            get(paste0(self$name, "_rr")),
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


      #' @description Get original input relative risks by age and sex.
      #' @return A copied data.table with the original relative risks.
      get_input_rr = function() {copy(private$input_rr)},


      #' @description Get original metadata.
      #' @return A list with the original metadata.
      get_metadata = function() {private$metadata},


      #' @description Get seed for RNG.
      #' @return A seed for the RNG that is produced by the digest of parameter
      #'   name and outcome.
      get_seed = function() {private$seed},


      #' @description Writes a template parameter file to disk.
      #' @param file_path Path including file name and .csvy extension to write the
      #'   file with placeholder parameter parameters.
      #' @return The `Parameter` object.
      write_xps_tmplte_file =
        function(file_path = "./inputs/RR/template.csvy") {
          file_path <- normalizePath(file_path, mustWork = FALSE)
          metadata <- list(
            "exps_name"     = "Parameter name",
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
      #' @return The `Parameter` object.
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



      #' @description Print the simulation parameters.
      #' @return The `Parameter` object.
      print = function() {
        print(paste0("Parameter name:      ", self$name))
        print(paste0("Outcome influenced: ", self$outcome))
        print(paste0("Distribution:       ", self$distribution))
        print(paste0("Source:             ", self$source))
        print(paste0("Notes:              ", self$notes))

        invisible(self)
      }
    ), # end of public

    # private ------------------------------------------------------------------
     private = list(
       exps_prm_file = NA_character_,
       seed = NA_integer_,
       input_rr = NA,
       effect = NA,
       metadata = NA,
       suffix = NA_character_,
       filenam = NA_character_,
       filenam_indx = NA_character_,
       cache = NA,
       cache_mc = NA,

       get_yaml_header =
         function (file, yaml_rxp = "^\\#*---[[:space:]]*$", verbose = TRUE) {
         con <- file(file, "r")
         on.exit(close(con))
         first_line <- readLines(con, n = 1L)
         if (!length(first_line) || !grepl(yaml_rxp, first_line)) {
           if (isTRUE(verbose)) {
             warning("No YAML header found")
           }
           return(NULL)
         }
         iline <- 1L
         closing_tag <- FALSE
         out <- character()
         while (!isTRUE(closing_tag)) {
           out[iline] <- readLines(con, n = 1L)
           if (grepl(yaml_rxp, out[iline])) {
             closing_tag <- TRUE
           }
           else {
             iline <- iline + 1L
           }
         }
         if (all(grepl("^#", out))) {
           out <- gsub("^#", "", out[-iline])
         }
         yaml::yaml.load(paste(out, collapse = "\n"))
       },

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

       stochRRtabl = # need to run by id
         function(m, ci, distribution = c("lognormal", "normal")) {
           distribution <- match.arg(distribution)
           kk <- dqrunif(1)

           if (distribution == "lognormal") {
             rr <- exp(qnorm(kk, log(m), abs(log(m) - log(ci)) / 1.96))
           } else if (distribution == "normal") {
             rr <- qnorm(kk, m, abs(m - ci) / 1.96)
           } else {
             stop("Only lognormal or normal are allowed.")
           }
           # rr[!is.finite(rr)] <- 1 # fix for rtruncnorm above
           return(rr)
         },

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
         colnam <- paste0(self$name, "_rr")
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

         # linear interpolation for intermediate parameter levels
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

       # calculate 95% CI of a ratio from a p value
       ci_from_p_forratio = function(mean_rr, p) {
         z <- -0.862 + sqrt(0.743 - 2.404 * log(p))
         est <- log(mean_rr)
         se <- abs(est/z)
         ci <- list()
         ci$lower_ci <- exp(est - 1.96 * se)
         ci$upper_ci <- exp(est + 1.96 * se)
         ci
       }



    ) # end of private
  )
