## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2025 University of Liverpool, Chris Kypridemos
##
## IMPACTncd_Engl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

#' @title Exposure Class
#' @name Exposure
#'
#' @description
#' An R6 class representing an exposure-disease relationship in the IMPACTncd framework.
#' This class encapsulates relative risk (RR) data between an exposure and a disease outcome,
#' supports stochastic modeling, and provides methods for generating, caching, and applying
#' relative risks in population health simulations.
#'
#' @details
#' The `Exposure` class manages exposure-disease relationships by:
#' \itemize{
#'   \item Loading relative risk data from CSVY files (CSV with YAML metadata)
#'   \item Generating stochastic relative risks using Monte Carlo methods
#'   \item Applying time lags between exposure and disease outcomes
#'   \item Supporting various distributions (lognormal, normal) for uncertainty modeling
#'   \item Caching computed relative risks for performance optimization
#'   \item Handling population stratification (age, sex, ethnicity, deprivation)
#' }
#'
#' The class supports both deterministic and stochastic modeling modes, with the ability
#' to apply exposure effects with specified time lags and uncertainty distributions.
#'
#' @section Initialization:
#' The class is initialized with a path to a CSVY file containing relative risk data
#' and a Design object with simulation parameters:
#'
#' ```r
#' # Initialize exposure object
#' bmi_chd <- Exposure$new("./inputs/RR/bmi_chd.csvy", design)
#' ```
#'
#' @section Public Fields:
#' \describe{
#'   \item{\code{name}}{Character. The name of the exposure (e.g., "bmi", "smoking").}
#'   \item{\code{outcome}}{Character. The name of the disease outcome the exposure influences.}
#'   \item{\code{lag}}{Integer. The time lag (in years) between exposure and disease onset.}
#'   \item{\code{distribution}}{Character. Distribution for Monte Carlo uncertainty ("lognormal" or "normal").}
#'   \item{\code{source}}{Character. Citation information for the effect size.}
#'   \item{\code{notes}}{Character. Additional notes regarding the exposure-outcome relationship.}
#' }
#'
#' @section Public Methods:
#' \describe{
#'   \item{\code{initialize(file_path, design)}}{Initializes the exposure object from a CSVY file.}
#'   \item{\code{gen_stochastic_effect(design, overwrite, smooth, ...)}}{Generates stochastic relative risks.}
#'   \item{\code{del_stochastic_effect(invert)}}{Deletes stochastic effect files from disk.}
#'   \item{\code{get_rr(mc, design, drop, plot_rr)}}{Retrieves relative risks for a Monte Carlo iteration.}
#'   \item{\code{clear_cache()}}{Clears the internal cache for relative risks.}
#'   \item{\code{xps_to_rr(sp, design, checkNAs, forPARF)}}{Applies relative risks to synthetic population.}
#'   \item{\code{validate_rr()}}{Validates and plots input vs stochastic relative risks.}
#'   \item{\code{get_input_rr()}}{Returns original input relative risks.}
#'   \item{\code{get_metadata()}}{Returns metadata from the CSVY file.}
#'   \item{\code{get_seed()}}{Returns the random seed used for this exposure.}
#'   \item{\code{write_xps_tmplte_file(file_path)}}{Writes a template exposure file.}
#'   \item{\code{convert_from_old_format(old_file, metadata, file_path, estimates, second_estimate_is_p)}}{Converts old CSV format to new CSVY format.}
#'   \item{\code{get_lag(mc)}}{Returns exposure lag for given Monte Carlo iteration.}
#'   \item{\code{get_ideal_xps_lvl(mc)}}{Returns ideal exposure level for given iteration.}
#'   \item{\code{get_name()}}{Returns the exposure-outcome combination name.}
#'   \item{\code{print()}}{Prints exposure information.}
#' }
#'
#' @section Private Methods:
#' \describe{
#'   \item{\code{write_xps_prm_file(dt, metadata, file_path)}}{Writes exposure parameter files.}
#'   \item{\code{stochRRtabl(m, ci, distribution)}}{Generates stochastic relative risks from distributions.}
#'   \item{\code{generate_rr_l(dt, tt, mc_max, smooth, do_checks, ...)}}{Generates age-specific relative risks with interpolation.}
#'   \item{\code{ci_from_p_forratio(mean_rr, p)}}{Calculates confidence intervals from p-values.}
#' }
#'
#' @section File Format:
#' The class expects CSVY files (CSV with YAML header) containing:
#' \describe{
#'   \item{YAML metadata}{exposure name, outcome, distribution, lag, source, notes}
#'   \item{CSV data}{relative risks by age group, sex, and optional stratification variables}
#' }
#'
#' @section Supported Stratification:
#' The class supports stratification by:
#' \itemize{
#'   \item Age groups (5-year bands from 0-99 years)
#'   \item Sex (men, women)
#'   \item Smoking status (4 levels: never, former, current light, current heavy)
#'   \item Ethnicity (9 categories: white, indian, pakistani, bangladeshi, other asian, black caribbean, black african, chinese, other)
#'   \item Deprivation (DIMD: 10 deciles, QIMD: 5 quintiles)
#' }
#'
#' @section Monte Carlo Features:
#' \itemize{
#'   \item Supports lognormal and normal distributions for uncertainty
#'   \item Generates deterministic seeds based on exposure-outcome combination
#'   \item Applies variable time lags using binomial distribution
#'   \item Caches results for performance optimization
#'   \item Supports smoothing with loess regression
#' }
#'
#' @section File Management:
#' \itemize{
#'   \item Automatically creates simulation/rr/ directory structure
#'   \item Uses content-based checksums for file versioning
#'   \item Stores separate files for relative risks and indexing
#'   \item Supports file cleanup and overwrite options
#' }
#'
#' @importFrom R6 R6Class
#' @importFrom data.table fread fwrite CJ setkeyv setnames setcolorder setnafill
#' @importFrom dqrng dqRNGkind dqset.seed dqrunif
#' @importFrom fst read_fst write_fst
#' @importFrom yaml as.yaml
#' @importFrom digest digest
#' @importFrom stats qnorm qbinom approx loess predict
#' @importFrom graphics layout plot points
#' @importFrom grDevices rgb
#'
#' @seealso \code{\link{Design}}, \code{\link[fst]{read_fst}}, \code{\link[yaml]{as.yaml}}
#'
#' @examples
#' \dontrun{
#' # Initialize exposure from CSVY file
#' bmi_chd <- Exposure$new("./inputs/RR/bmi_chd.csvy", design)
#'
#' # Generate stochastic effects
#' bmi_chd$gen_stochastic_effect(design, overwrite = TRUE)
#'
#' # Get relative risks for Monte Carlo iteration 1
#' rr_data <- bmi_chd$get_rr(mc = 1, design)
#'
#' # Apply to synthetic population
#' bmi_chd$xps_to_rr(synthpop, design)
#'
#' # Validate relative risks
#' bmi_chd$validate_rr()
#'
#' # Create template file
#' Exposure$new()$write_xps_tmplte_file("./template.csvy")
#' }
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
      #' @description
      #' Initialize a new Exposure object by reading exposure parameters from a CSVY file.
      #' This method loads relative risk data, validates the file structure, processes
      #' stratification variables, and sets up Monte Carlo parameters for stochastic modeling.
      #'
      #' @param sRelativeRiskByPopulationSubsetForExposureFilePath Character string.
      #'   Path to a CSVY file containing relative risk (RR) data by population subset
      #'   (age, sex, and optionally deprivation index, ethnicity, smoking status).
      #'   The file must have a YAML header with exposure metadata and CSV data with RR values.
      #' @param design A \code{Design} object containing simulation parameters including
      #'   maximum lag, iteration count, age limits, and logging preferences.
      #'
      #' @return A new \code{Exposure} object with populated fields and prepared data structures.
      #'
      #' @details
      #' The method performs several key operations:
      #' \itemize{
      #'   \item Validates file existence and Design object
      #'   \item Reads CSVY file using \code{fread} with YAML parsing
      #'   \item Processes factor levels for stratification variables (smoking, ethnicity, deprivation)
      #'   \item Validates metadata fields (name, outcome, distribution, lag)
      #'   \item Sets up Monte Carlo lag distribution using binomial sampling
      #'   \item Generates deterministic seed from exposure-outcome combination
      #'   \item Prepares file paths for stochastic effect storage
      #'   \item Expands relative risk data to full age range with interpolation
      #' }
      #'
      #' @examples
      #' \dontrun{
      #' # Initialize BMI-CHD exposure relationship
      #' bmi_chd <- Exposure$new("./inputs/RR/bmi_chd.csvy", design)
      #'
      #' # Initialize with smoking exposure
      #' smoking_stroke <- Exposure$new("./inputs/RR/smoking_stroke.csvy", design)
      #' }
      initialize = function(
        sRelativeRiskByPopulationSubsetForExposureFilePath,
        design
      ) {
        sRelativeRiskByPopulationSubsetForExposureFilePath <- normalizePath(
          sRelativeRiskByPopulationSubsetForExposureFilePath,
          mustWork = TRUE
        )

        if (!inherits(design, "Design"))
          stop("Argument design needs to be a Design object.")

        if (file.exists(sRelativeRiskByPopulationSubsetForExposureFilePath)) {
          # TODO add some checks to ensure proper structure
          dtRelativeRiskByPopulationSubset <- fread(
            sRelativeRiskByPopulationSubsetForExposureFilePath,
            stringsAsFactors = TRUE,
            yaml = TRUE
          )

          # NOTE This needs manual update every time a new factor (other than
          # smok_status, ethnicity, dimd) appears as a new stratum in a RR file.
          # TODO add some checks and perhaps automate to some extent. Add
          # documentation to avoid silent errors
          if ("smok_status" %in% names(dtRelativeRiskByPopulationSubset)) {
            # Validate that smok_status is numeric or factor
            if (!is.numeric(dtRelativeRiskByPopulationSubset$smok_status) && 
                !is.factor(dtRelativeRiskByPopulationSubset$smok_status)) {
              stop("Column 'smok_status' must be numeric or factor, but found: ", 
                   typeof(dtRelativeRiskByPopulationSubset$smok_status)[1])
            }
            dtRelativeRiskByPopulationSubset[,
              smok_status := factor(smok_status, levels = 1:(max(smok_status, na.rm = TRUE)))
            ]
          } # end of smok_status
          if ("ethnicity" %in% names(dtRelativeRiskByPopulationSubset)) {
            # Define expected ethnicity levels
            expected_ethnicity_levels <- c(
              "white",
              "indian", 
              "pakistani",
              "bangladeshi",
              "other asian",
              "black caribbean",
              "black african",
              "chinese",
              "other"
            )
            
            # Validate that ethnicity values match expected levels exactly
            unique_ethnicity <- unique(as.character(dtRelativeRiskByPopulationSubset$ethnicity))
            
            # Use symmetric difference to find any mismatched levels
            mismatched_levels <- setdiff(union(unique_ethnicity, expected_ethnicity_levels), 
                                       intersect(unique_ethnicity, expected_ethnicity_levels))
            
            if (length(mismatched_levels) > 0) {
              # Identify which are extra vs missing for detailed error message
              extra_levels <- setdiff(unique_ethnicity, expected_ethnicity_levels)
              missing_levels <- setdiff(expected_ethnicity_levels, unique_ethnicity)
              
              error_msg <- "Column 'ethnicity' validation failed:"
              if (length(extra_levels) > 0) {
                error_msg <- paste0(error_msg, "\n  - Invalid levels found: ", 
                                   paste(extra_levels, collapse = ", "))
              }
              if (length(missing_levels) > 0) {
                error_msg <- paste0(error_msg, "\n  - Missing expected levels: ", 
                                   paste(missing_levels, collapse = ", "))
              }
              error_msg <- paste0(error_msg, "\n  - Expected levels are: ", 
                                 paste(expected_ethnicity_levels, collapse = ", "))
              stop(error_msg)
            }
            
            dtRelativeRiskByPopulationSubset[,
              ethnicity := factor(
                ethnicity,
                levels = expected_ethnicity_levels
              )
            ]
          } # end of ethnicity

          if ("dimd" %in% names(dtRelativeRiskByPopulationSubset))
            dtRelativeRiskByPopulationSubset[,
              dimd := factor(
                dimd,
                levels = c(
                  "1 most deprived",
                  2:9,
                  "10 least deprived"
                )
              )
            ]

          if ("qimd" %in% names(dtRelativeRiskByPopulationSubset))
            dtRelativeRiskByPopulationSubset[,
              qimd := factor(
                qimd,
                levels = c(
                  "1 most deprived",
                  2:4,
                  "5 least deprived"
                )
              )
            ]

          dtRelativeRiskByPopulationSubset[,
            agegroup := relevel(agegroup, "<1")
          ]

          nam <- setdiff(
            names(dtRelativeRiskByPopulationSubset),
            c("rr", "ci_rr")
          )
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
        if (
          length(design$sim_prm$ignore_xps) > 0L &&
            (identical(metadata$xps_name, design$sim_prm$ignore_xps) ||
              (metadata$xps_name == "met" &&
                design$sim_prm$ignore_xps == "active_days") ||
              (grepl("^smok_", metadata$xps_name) &&
                design$sim_prm$ignore_xps == "smoking"))
        ) {
          set(dtRelativeRiskByPopulationSubset, NULL, "rr", 1)
          set(dtRelativeRiskByPopulationSubset, NULL, "ci_rr", 1)
        }

        # Validation
        stopifnot(
          c("xps_name", "outcome", "distribution", "lag") %in% names(metadata),
          metadata$lag >= 1 |
            (metadata$lag == 0L & metadata$incidence$type == 1L),
          metadata$lag <= design$sim_prm$maxlag,
          length(metadata$distribution) == 1L,
          metadata$distribution %in% c("lognormal", "normal")
        )

        self$name <- metadata$xps_name
        self$outcome <- metadata$outcome
        self$lag <- metadata$lag
        self$distribution <- metadata$distribution
        self$source <- metadata$source
        self$notes <- metadata$notes

        private$nam_rr <- paste0(self$name, "_rr") # i.e. bmi_rr

        if (
          "apply_rr_extra_fn" %in%
            names(metadata) &&
            startsWith(metadata$apply_rr_extra_fn, "function(")
        ) {
          # https://www.r-bloggers.com/2016/08/yaml-define-an-r-function/
          # NOTE check that string starts with function for an added layer
          # of security. In any case still very insecure.
          private$apply_rr_extra <- eval(str2lang(metadata$apply_rr_extra_fn))
        } else private$apply_rr_extra <- function(...) NULL

        dqRNGkind("pcg64")

        # NOTE that for the different dimensions of smoking per disease we need
        # the same lag, so quit-yrs applied with the correct lag. Hence self$name
        private$seed <- abs(digest2int(
          paste0(
            ifelse(
              grepl("^smok_", self$name),
              "smoking",
              self$name
            ), # rename all smoking dimensions to smoking
            self$outcome
          ),
          seed = 764529L
        ))
        private$chksum <- digest(dtRelativeRiskByPopulationSubset) # NOTE not affected by metadata changes
        private$suffix <- paste0(self$name, "~", self$outcome)

        private$filedir <- file.path(getwd(), "simulation", "rr")
        if (!dir.exists(private$filedir))
          dir.create(private$filedir, recursive = TRUE)
        private$filenam <- file.path(
          private$filedir,
          paste0("rr_", private$suffix, "_", private$chksum, "_l.fst")
        )
        private$filenam_indx <- file.path(
          private$filedir,
          paste0("rr_", private$suffix, "_", private$chksum, "_indx.fst")
        )

        dqset.seed(private$seed, stream = 1) # stream = 1 for lags

        if (self$lag == 0L) {
          # Only allowed for incd type 1

          private$lag_mc <-
            rep(0L, times = design$sim_prm$iteration_n_max)
        } else {
          private$lag_mc <-
            1L +
            as.integer(qbinom(
              dqrunif(design$sim_prm$iteration_n_max),
              design$sim_prm$maxlag - 1L,
              (self$lag - 1) / (design$sim_prm$maxlag - 1L)
            )) # Note qbinom returns double not int
        }

        if (
          "ideal_xps_lvl_fn" %in%
            names(metadata) &&
            startsWith(metadata$ideal_xps_lvl_fn, "function(")
        ) {
          # NOTE check that string starts with function for an added layer
          # of security. In any case still very insecure.
          foo <- eval(str2lang(metadata$ideal_xps_lvl_fn))
          private$ideal_xps_lvl_mc <- foo(design)
          rm(foo)
        }

        if (length(nam) == 2L) {
          tt <-
            CJ(
              age = (design$sim_prm$ageL -
                design$sim_prm$maxlag):design$sim_prm$ageH,
              sex = factor(c("men", "women"))
            )
          setkeyv(tt, c("sex", "age"))
        } else if (length(nam) == 3L) {
          nam <- setdiff(nam, c("sex", "agegroup"))
          if (is.numeric(dtRelativeRiskByPopulationSubset[[nam]])) {
            t3 <- min((dtRelativeRiskByPopulationSubset[[nam]])):max(
              (dtRelativeRiskByPopulationSubset[[nam]])
            )
            t4 <- sort(unique(dtRelativeRiskByPopulationSubset[[nam]]))
            # check if consecutive elements
            if (!isTRUE(all.equal(t3, t4))) interpolate <- TRUE

            tt <-
              CJ(
                age = (design$sim_prm$ageL -
                  design$sim_prm$maxlag):design$sim_prm$ageH,
                sex = factor(c("men", "women")),
                V3 = t3
              )
          } else {
            # if not numeric
            tt <-
              CJ(
                age = (design$sim_prm$ageL -
                  design$sim_prm$maxlag):design$sim_prm$ageH,
                sex = factor(c("men", "women")),
                V3 = unique(dtRelativeRiskByPopulationSubset[[nam]])
              )
          }

          setnames(tt, "V3", nam)
          setkeyv(tt, c("sex", "age", nam))
        } else {
          # length(nam) > 3
          stop(
            "Only one additional column beyond sex and agegroup is currently supported."
          )
        }

        to_agegrp(tt, 5L, max(tt$age), "age", "agegroup")

        private$effect <- copy(dtRelativeRiskByPopulationSubset)
        private$metadata <- metadata
        private$xps_prm_file <- sRelativeRiskByPopulationSubsetForExposureFilePath

        private$input_rr <-
          setcolorder(
            dtRelativeRiskByPopulationSubset[tt, on = .NATURAL],
            c("age", "agegroup", "sex")
          )
        private$input_rr[, "ci_rr" := NULL] # NOTE agegroup necessary for gen_stochastic_RR
        if (exists("interpolate") && isTRUE(interpolate)) {
          # linear interpolation
          private$input_rr[,
            rr := approx(get(nam), rr, get(nam))$y,
            by = .(age, sex)
          ]
        }
        if ("smok_status" %in% names(private$input_rr))
          private$input_rr[, smok_status := factor(smok_status, levels = 1:4)]
        if ("ethnicity" %in% names(private$input_rr))
          private$input_rr[,
            ethnicity := factor(
              ethnicity,
              levels = c(
                "white",
                "indian",
                "pakistani",
                "bangladeshi",
                "other asian",
                "black caribbean",
                "black african",
                "chinese",
                "other"
              )
            )
          ]
        invisible(self)
      },

      # gen_stochastic_effect ----
      #' @description
      #' Generate and write stochastic relative risk effects to disk for Monte Carlo simulations.
      #' This method creates uncertainty distributions around the point estimates of relative risks
      #' and stores them in optimized FST format files for fast access during simulations.
      #'
      #' @param design_ A \code{Design} object containing simulation parameters. Defaults to \code{design}.
      #' @param overwrite Logical. If \code{TRUE}, overwrites existing stochastic effect files.
      #'   If \code{FALSE}, skips generation if files already exist. Default is \code{FALSE}.
      #' @param smooth Logical. If \code{TRUE}, applies loess smoothing to relative risks across age.
      #'   This can help create smoother age-specific risk profiles.
      #' @param ... Additional arguments passed to \code{loess} function when \code{smooth = TRUE}.
      #'   Common arguments include \code{span} (smoothing parameter, default 0.7) and
      #'   \code{degree} (polynomial degree, default 1).
      #'
      #' @return The \code{Exposure} object, invisibly, for method chaining.
      #'
      #' @details
      #' The method generates stochastic relative risks by:
      #' \itemize{
      #'   \item Sampling from specified distributions (lognormal/normal) using confidence intervals
      #'   \item Expanding data to full age range with linear interpolation if needed
      #'   \item Applying optional loess smoothing across age groups
      #'   \item Enforcing constraints (RR >= 1 or RR <= 1) based on input data
      #'   \item Writing results to FST files with indexing for fast Monte Carlo access
      #'   \item Creating checksums to ensure data integrity and version control
      #' }
      #'
      #' Files are stored in \code{simulation/rr/} directory with names containing
      #' exposure-outcome combination and data checksum for versioning.
      #'
      #' @examples
      #' \dontrun{
      #' # Generate stochastic effects with default settings
      #' bmi_chd$gen_stochastic_effect(design)
      #'
      #' # Overwrite existing files with smoothing
      #' bmi_chd$gen_stochastic_effect(design, overwrite = TRUE, smooth = TRUE, span = 0.5)
      #'
      #' # Generate without smoothing, overwriting existing
      #' smoking_copd$gen_stochastic_effect(design, overwrite = TRUE, smooth = FALSE)
      #' }
      gen_stochastic_effect = function(
        design_ = design,
        overwrite = FALSE,
        smooth,
        ...
      ) {
        if (!inherits(design_, "Design"))
          stop("Argument design needs to be a Design object.")

        if (
          overwrite ||
            all(
              !overwrite,
              any(
                !file.exists(private$filenam),
                !file.exists(private$filenam_indx)
              )
            ) ||
            (max(read_fst(private$filenam, columns = "age")) <
              design_$sim_prm$ageH)
        ) {
          if (design_$sim_prm$logs)
            message(paste0(
              self$name,
              "~",
              self$outcome,
              " simulate RR uncertainty."
            ))

          self$del_stochastic_effect(invert = TRUE) # deletes previous versions

          stoch_effect <-
            private$generate_rr_l(
              private$effect,
              private$input_rr[, .SD, .SDcols = !c("rr")],
              design_$sim_prm$iteration_n_max,
              smooth,
              design_$sim_prm$logs,
              ...
            )
          if ("smok_status" %in% names(stoch_effect))
            stoch_effect[, smok_status := factor(smok_status, levels = 1:4)]
          if ("ethnicity" %in% names(stoch_effect))
            stoch_effect[,
              ethnicity := factor(
                ethnicity,
                levels = c(
                  "white",
                  "indian",
                  "pakistani",
                  "bangladeshi",
                  "other asian",
                  "black caribbean",
                  "black african",
                  "chinese",
                  "other"
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
            message(paste0(
              self$name,
              "~",
              self$outcome,
              " RR files exist already."
            ))
        }

        invisible(self)
      },

      # del_stochastic_effect ----
      #' @description
      #' Delete stochastic effect files from disk to free storage space or force regeneration.
      #' This method can either delete all files for this exposure-outcome combination or
      #' keep only the current version (based on data checksum) while removing outdated versions.
      #'
      #' @param invert Logical. If \code{FALSE} (default), deletes the current stochastic effect
      #'   files. If \code{TRUE}, keeps files with the current checksum and deletes all other
      #'   versions for this exposure-outcome combination (useful for cleanup).
      #'
      #' @return The \code{Exposure} object, invisibly, for method chaining.
      #'
      #' @details
      #' Files deleted include both the main relative risk data file and the indexing file
      #' used for fast Monte Carlo iteration access. When \code{invert = TRUE}, this serves
      #' as a cleanup function to remove outdated versions while preserving current data.
      #'
      #' @examples
      #' \dontrun{
      #' # Delete current stochastic effects
      #' bmi_chd$del_stochastic_effect()
      #'
      #' # Clean up old versions, keep current
      #' smoking_copd$del_stochastic_effect(invert = TRUE)
      #' }
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
      #' @description
      #' Retrieve relative risk data from disk for a specific Monte Carlo iteration.
      #' This method loads stochastic relative risks (if available) or returns deterministic
      #' values, with built-in caching for performance optimization.
      #'
      #' @param mc Integer. Monte Carlo iteration number between 1 and \code{design$iteration_n_max}.
      #'   Must be a single integer value.
      #' @param design_ A \code{Design} object containing simulation parameters. Defaults to \code{design}.
      #' @param drop Logical. If \code{TRUE} and all relative risks are identical across age/sex,
      #'   returns a scalar numeric value instead of a data.table. Default is \code{FALSE}.
      #' @param plot_rr Logical. If \code{TRUE}, creates a diagnostic plot showing relative risks
      #'   by age and sex, overlaying stochastic values with original input data. Default is \code{FALSE}.
      #'
      #' @return A data.table with relative risk values by age and sex (and other stratification
      #'   variables if present), or a scalar numeric if \code{drop = TRUE} and values are uniform.
      #'   Column names include age, sex, and the exposure-specific RR column (e.g., "bmi_rr").
      #'
      #' @details
      #' The method behavior depends on simulation mode:
      #' \itemize{
      #'   \item \strong{Stochastic mode}: Reads pre-generated stochastic effects from FST files
      #'   \item \strong{Deterministic mode}: Returns original input relative risks
      #'   \item \strong{Caching}: Stores last accessed data to avoid repeated file I/O
      #'   \item \strong{Validation}: Ensures Monte Carlo iteration is within valid range
      #' }
      #'
      #' The plot option provides visual validation of stochastic vs input relative risks,
      #' with points showing original data and scattered points showing stochastic variation.
      #'
      #' @examples
      #' \dontrun{
      #' # Get relative risks for Monte Carlo iteration 1
      #' rr_data <- bmi_chd$get_rr(mc = 1, design)
      #'
      #' # Get scalar value if uniform across age/sex
      #' uniform_rr <- simple_exposure$get_rr(mc = 1, design, drop = TRUE)
      #'
      #' # Create diagnostic plot
      #' rr_with_plot <- smoking_copd$get_rr(mc = 1, design, plot_rr = TRUE)
      #'
      #' # Access from cache (subsequent calls with same mc)
      #' cached_rr <- bmi_chd$get_rr(mc = 1, design)  # Fast access
      #' }
      get_rr = function(mc, design_ = design, drop = FALSE, plot_rr = FALSE) {
        if (!inherits(design_, "Design"))
          stop("Argument design_ needs to be a Design object.")

        stopifnot(
          length(mc) == 1L,
          between(mc, 1L, design_$sim_prm$iteration_n_max)
        )
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
          points(
            private$input_rr$age,
            private$input_rr$rr,
            col = private$input_rr$sex
          )
        }

        if (drop && uniqueN(out[[private$nam_rr]]) == 1L)
          out <- unique(out[[private$nam_rr]])

        out[]
      },

      # clear_cache ----
      #' @description
      #' Clear the internal cache used for storing recently accessed relative risk data.
      #' This method frees memory and forces fresh data loading on next access.
      #'
      #' @return The \code{Exposure} object, invisibly, for method chaining.
      #'
      #' @details
      #' The cache stores the most recently accessed Monte Carlo iteration data to avoid
      #' repeated file I/O operations. Clearing it may be useful when memory is limited
      #' or when ensuring fresh data access is critical.
      #'
      #' @examples
      #' \dontrun{
      #' # Clear cache to free memory
      #' bmi_chd$clear_cache()
      #'
      #' # Chain with other operations
      #' exposure$clear_cache()$get_rr(mc = 1, design)
      #' }
      clear_cache = function() {
        private$cache <- NA
        private$cache_mc <- NA
        message("Cache was cleared")

        invisible(self)
      },

      # xps_to_rr ----
      #' @description
      #' Apply relative risk values to a synthetic population based on exposure levels.
      #' This method creates a new column in the population data with relative risks
      #' corresponding to each individual's exposure level, age, sex, and other characteristics.
      #' It handles time lags and special cases like smoking cessation effects.
      #'
      #' @param sp A \code{SynthPop} object containing the synthetic population data.
      #'   Must have columns for the exposure variable, age, sex, and stratification variables.
      #' @param design_ A \code{Design} object containing simulation parameters including
      #'   age limits, logging preferences, and lag settings.
      #' @param checkNAs Logical. If \code{TRUE}, prints diagnostic tables showing NA patterns
      #'   before they are replaced with neutral values (RR = 1). Default follows \code{design$logs}.
      #'   Note that NAs may be expected for certain exposure levels.
      #' @param forPARF Logical. If \code{TRUE}, applies Population Attributable Risk Fraction
      #'   mode without time lags. Used for counterfactual analyses. Default is \code{FALSE}.
      #'
      #' @return The \code{Exposure} object, invisibly, for method chaining.
      #'   The \code{sp$pop} data.table is modified in place with a new relative risk column.
      #'
      #' @details
      #' The method performs several operations:
      #' \itemize{
      #'   \item Creates lagged exposure variables using \code{shift_bypid} (unless \code{forPARF = TRUE})
      #'   \item Matches population characteristics to relative risk lookup tables
      #'   \item Applies special transformation functions if defined in metadata (e.g., smoking quit years)
      #'   \item Handles missing values by setting them to neutral effect (RR = 1)
      #'   \item Preserves original exposure variables while creating temporary working copies
      #' }
      #'
      #' For smoking exposures, quit years may modify existing smoking relative risks rather
      #' than creating new columns. The method handles conflicts with disease prevalence
      #' variables by temporarily renaming them.
      #'
      #' @examples
      #' \dontrun{
      #' # Apply BMI relative risks to population
      #' bmi_chd$xps_to_rr(synthpop, design)
      #'
      #' # Apply with NA checking disabled
      #' smoking_copd$xps_to_rr(synthpop, design, checkNAs = FALSE)
      #'
      #' # Apply for PARF calculation (no lags)
      #' alcohol_stroke$xps_to_rr(parf_population, design, forPARF = TRUE)
      #'
      #' # Chain with other operations
      #' exposure$gen_stochastic_effect(design)$xps_to_rr(synthpop, design)
      #' }

      xps_to_rr = function(
        sp,
        design_,
        checkNAs = design_$sim_prm$logs,
        forPARF = FALSE
      ) {
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
          lookup_dt(
            sp$pop,
            self$get_input_rr(),
            check_lookup_tbl_validity = design_$sim_prm$logs
          )
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
          set(
            sp$pop,
            NULL,
            self$name, # column without _curr_xps is lagged
            shift_bypid(
              sp$pop[[xps_tolag]],
              self$get_lag(sp$mc_aggr),
              sp$pop$pid,
              rw
            )
          )
          # setnafill(sp$pop, "nocb", cols = self$name)
          lookup_dt(
            sp$pop,
            self$get_rr(sp$mc_aggr, design_, drop = FALSE),
            check_lookup_tbl_validity = design_$sim_prm$logs
          )
        }

        private$apply_rr_extra(sp)

        if (checkNAs && private$nam_rr %in% names(sp$pop)) {
          # NOTE in case of smok_quit_yrs the risk column is deleted from the
          # rr_extra_fn()
          print(self$name)
          print(sp$pop[
            !is.na(get(self$name)) &
              is.na(get(private$nam_rr)) &
              age >= design_$sim_prm$ageL &
              year >= design_$sim_prm$init_year,
            table(get(self$name), useNA = "always")
          ])
        }

        if (
          private$nam_rr %in%
            names(sp$pop) &&
            is.numeric(sp$pop[[private$nam_rr]])
        )
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
      #' @description
      #' Validate stochastic relative risks by creating diagnostic plots comparing
      #' input data with generated stochastic effects. This helps ensure the
      #' Monte Carlo generation process is working correctly.
      #'
      #' @return \code{NULL}. The method creates plots as a side effect.
      #'
      #' @details
      #' Creates a two-panel plot showing:
      #' \itemize{
      #'   \item Stochastic relative risks as semi-transparent points (showing uncertainty)
      #'   \item Original input relative risks as solid points (showing central estimates)
      #'   \item Separate panels for men and women
      #'   \item Age on x-axis, relative risk on y-axis
      #' }
      #'
      #' This visualization helps identify issues such as:
      #' \itemize{
      #'   \item Excessive or insufficient uncertainty
      #'   \item Bias in stochastic generation
      #'   \item Age-specific patterns in uncertainty
      #' }
      #'
      #' @examples
      #' \dontrun{
      #' # Validate after generating stochastic effects
      #' bmi_chd$gen_stochastic_effect(design)
      #' bmi_chd$validate_rr()
      #'
      #' # Quick validation workflow
      #' smoking_copd$gen_stochastic_effect(design, overwrite = TRUE)$validate_rr()
      #' }
      validate_rr = function() {
        layout(matrix(1:2))
        out <-
          read_fst(private$filenam, as.data.table = TRUE)
        out[, mc := NULL]

        layout(matrix(1:2))
        out[
          sex == "men",
          plot(
            age,
            get(private$nam_rr),
            col = rgb(red = 0, green = 0, blue = 1, alpha = 0.05),
            pch = ".",
            ylab = "RR",
            main = (paste0(self$name, "_rr (men)"))
          )
        ]
        private$input_rr[sex == "men", points(age, rr)]
        out[
          sex == "women",
          plot(
            age,
            get(private$nam_rr),
            col = rgb(red = 1, green = 0, blue = 0, alpha = 0.05),
            pch = ".",
            ylab = "RR",
            main = (paste0(self$name, "_rr (women)"))
          )
        ]
        private$input_rr[sex == "women", points(age, rr)]
        layout(matrix(1:1))
        # points(private$input_rr$age,
        #   private$input_rr$rr,
        #   col = private$input_rr$sex)
      },

      # get_input_rr ----
      #' @description
      #' Get the original input relative risk data by age and sex without Monte Carlo variation.
      #' This method returns the deterministic relative risks as loaded from the CSVY file,
      #' useful for comparisons with stochastic effects or PARF calculations.
      #'
      #' @return A data.table containing the original relative risks with columns for age,
      #'   sex, and the exposure-specific relative risk column (e.g., "bmi_rr").
      #'   The agegroup column is removed and data is cached for performance.
      #'
      #' @details
      #' This method:
      #' \itemize{
      #'   \item Returns deterministic (non-stochastic) relative risks
      #'   \item Removes the agegroup column (keeps only age)
      #'   \item Renames "rr" column to exposure-specific name
      #'   \item Caches result with special "forPARF" identifier
      #'   \item Provides baseline data for validation and PARF calculations
      #' }
      #'
      #' @examples
      #' \dontrun{
      #' # Get original input relative risks
      #' original_rr <- bmi_chd$get_input_rr()
      #'
      #' # Compare with stochastic version
      #' stoch_rr <- bmi_chd$get_rr(mc = 1, design)
      #' }
      get_input_rr = function() {
        out <- copy(private$input_rr)
        out[, agegroup := NULL]
        setnames(out, "rr", private$nam_rr)
        private$cache <- copy(out)
        private$cache_mc <- "forPARF"
        out
      },

      # get_metadata ----
      #' @description
      #' Get the original metadata associated with the exposure-disease relationship.
      #' Returns the complete metadata list as loaded from the CSVY file header,
      #' including exposure name, outcome, distribution, source, and notes.
      #'
      #' @return A list containing the original metadata with elements:
      #' \itemize{
      #'   \item \code{xps_name} - The exposure name
      #'   \item \code{outcome} - The disease outcome
      #'   \item \code{distribution} - Uncertainty distribution type
      #'   \item \code{source} - Citation or data source
      #'   \item \code{notes} - Additional notes
      #'   \item \code{lag} - Exposure lag in years (if specified)
      #' }
      #'
      #' @examples
      #' \dontrun{
      #' # Get metadata for inspection
      #' meta <- bmi_chd$get_metadata()
      #' print(meta$source)
      #'
      #' # Check distribution type
      #' if (smoking_copd$get_metadata()$distribution == "lognormal") {
      #'   message("Using lognormal distribution")
      #' }
      #' }
      get_metadata = function() {
        private$metadata
      },

      # get_seed ----
      #' @description
      #' Get the random number generator seed used for Monte Carlo simulations.
      #' The seed is generated as a hash digest of the exposure name and outcome,
      #' ensuring reproducible stochastic effects for the same exposure-outcome pair.
      #'
      #' @return An integer representing the RNG seed used for this exposure's
      #'   Monte Carlo simulations.
      #'
      #' @details
      #' The seed is automatically generated during initialization using:
      #' \itemize{
      #'   \item A digest hash of exposure name and outcome
      #'   \item Ensures reproducibility across simulation runs
      #'   \item Different exposure-outcome pairs get different seeds
      #'   \item Used by \code{dqset.seed()} for high-quality random number generation
      #' }
      #'
      #' @examples
      #' \dontrun{
      #' # Get the seed for reproducibility checks
      #' seed_val <- bmi_chd$get_seed()
      #'
      #' # Verify different exposures have different seeds
      #' bmi_seed <- bmi_chd$get_seed()
      #' smoking_seed <- smoking_copd$get_seed()
      #' identical(bmi_seed, smoking_seed)  # Should be FALSE
      #' }
      get_seed = function() {
        private$seed
      },

      # write_xps_tmplte_file ----
      #' @description
      #' Write a template exposure file to disk in CSVY format. This creates a
      #' skeleton file with placeholder metadata and a full age-sex grid of
      #' relative risks set to 1.0, which can be edited to create new exposures.
      #'
      #' @param file_path Character string. Path where the template file will be written,
      #'   including filename and .csvy extension. Defaults to "./inputs/RR/template.csvy".
      #'
      #' @return The \code{Exposure} object, invisibly, for method chaining.
      #'
      #' @details
      #' Creates a template file containing:
      #' \itemize{
      #'   \item YAML metadata header with placeholder values
      #'   \item Full age-sex grid from 0-99 years in 5-year groups
      #'   \item All relative risks initialized to 1.0
      #'   \item All confidence intervals initialized to 1.0
      #'   \item Proper factor levels for agegroup and sex
      #' }
      #'
      #' The template provides the correct structure for new exposure files and
      #' can be customized by replacing placeholder values with actual data.
      #'
      #' @examples
      #' \dontrun{
      #' # Create template file in default location
      #' exposure$write_xps_tmplte_file()
      #'
      #' # Create template with custom path
      #' exposure$write_xps_tmplte_file("./my_exposure_template.csvy")
      #'
      #' # Method chaining
      #' exposure$write_xps_tmplte_file()$print()
      #' }
      write_xps_tmplte_file = function(
        file_path = "./inputs/RR/template.csvy"
      ) {
        file_path <- normalizePath(file_path, mustWork = FALSE)
        metadata <- list(
          "xps_name" = "Exposure name",
          "outcome" = "Disease name",
          "distribution" = "lognormal or normal",
          "source" = "Citation for the effectsize",
          "notes" = "Any additional notes"
        )
        effect <- CJ(
          agegroup = factor(agegrp_name(
            min_age = 0,
            max_age = 99,
            grp_width = 5
          )),
          sex = c("men", "women"),
          rr = 1,
          ci_rr = 1
        )
        setkey(effect, sex, agegroup)

        private$write_xps_prm_file(effect, metadata, file_path)
        invisible(self)
      },

      # convert_from_old_format ----
      #' @description
      #' Convert legacy CSV format exposure files to the new CSVY format with metadata header.
      #' This method facilitates migration from older IMPACTncd versions by reading old
      #' CSV files and creating properly formatted CSVY files with YAML metadata.
      #'
      #' @param old_file Character string. Path to the legacy CSV file containing relative
      #'   risk data. If missing, uses the \code{estimates} parameter instead.
      #' @param metadata List containing metadata fields required for the new format:
      #'   xps_name, outcome, distribution, source, notes, and optionally lag.
      #' @param file_path Character string. Output path for the new CSVY file,
      #'   including filename and .csvy extension.
      #' @param estimates Numeric vector of length 2. Only used when \code{old_file} is missing.
      #'   First element is the point estimate, second is either CI or p-value.
      #' @param second_estimate_is_p Logical. If \code{TRUE} and \code{estimates} is provided,
      #'   interprets the second element as a p-value rather than confidence interval.
      #'
      #' @return The \code{Exposure} object, invisibly, for method chaining.
      #'
      #' @details
      #' The conversion process:
      #' \itemize{
      #'   \item Reads legacy CSV with age groups and sex stratification
      #'   \item Expands to full 0-99 age range with 5-year groups
      #'   \item Converts sex coding from numeric (1,2) to factor ("men","women")
      #'   \item Handles single estimate expansion across age/sex groups
      #'   \item Adds YAML metadata header for new format compatibility
      #'   \item Supports p-value to confidence interval conversion
      #' }
      #'
      #' @examples
      #' \dontrun{
      #' # Convert existing CSV file
      #' metadata <- list(
      #'   xps_name = "bmi", outcome = "chd", distribution = "lognormal",
      #'   source = "Meta-analysis 2020", notes = "Converted from legacy format"
      #' )
      #' exposure$convert_from_old_format(
      #'   old_file = "legacy_bmi_chd.csv",
      #'   metadata = metadata,
      #'   file_path = "bmi_chd.csvy"
      #' )
      #'
      #' # Create from estimates with p-value
      #' exposure$convert_from_old_format(
      #'   metadata = metadata,
      #'   file_path = "new_exposure.csvy",
      #'   estimates = c(1.2, 0.05),
      #'   second_estimate_is_p = TRUE
      #' )
      #' }
      convert_from_old_format = function(
        old_file,
        metadata,
        file_path,
        estimates = NULL,
        second_estimate_is_p = FALSE
      ) {
        file_path <- normalizePath(file_path, mustWork = FALSE)
        colby <- character(0)

        if (!missing(old_file)) {
          old_file <- normalizePath(old_file, mustWork = TRUE)

          effect <-
            fread(old_file, stringsAsFactors = TRUE)
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
            agegroup = factor(agegrp_name(
              min_age = 30,
              max_age = 89,
              grp_width = 5
            )),
            sex = factor(c("men", "women")),
            mean.rr = estimates[[1]],
            ci.rr = private$ci_from_p_forratio(estimates[[1]], estimates[[2]])[[
              2
            ]]
          )
          colby <- c(colby, "agegroup", "sex")
        } else {
          effect <- CJ(
            agegroup = factor(agegrp_name(
              min_age = 30,
              max_age = 89,
              grp_width = 5
            )),
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
            agegroup = factor(agegrp_name(
              min_age = 0,
              max_age = 99,
              grp_width = 5
            )),
            sex = factor(c("men", "women")),
            rr = 1,
            ci_rr = 1
          )
        } else if (length(nam) == 1L) {
          full <- CJ(
            agegroup = factor(agegrp_name(
              min_age = 0,
              max_age = 99,
              grp_width = 5
            )),
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

        full[
          effect,
          on = colby,
          `:=`(
            rr = i.mean.rr,
            ci_rr = i.ci.rr
          )
        ]

        setkeyv(full, c(nam, "sex", "agegroup"))

        private$write_xps_prm_file(full, metadata, file_path)
        invisible(self)
      },

      # get_lag ----
      #' @description
      #' Get the exposure lag for the exposure-disease relationship. Returns either
      #' the fixed lag value or Monte Carlo-specific lag values for stochastic simulations.
      #'
      #' @param mc_ Integer vector. Monte Carlo iteration numbers. If missing or 0,
      #'   returns the median (deterministic) lag value. Otherwise returns lag values
      #'   for the specified Monte Carlo iterations.
      #'
      #' @return Integer vector containing exposure lag(s) in years. Length matches
      #'   the length of \code{mc_} parameter, or single value if \code{mc_} is missing/0.
      #'
      #' @details
      #' Exposure lag represents the time delay between exposure and disease onset:
      #' \itemize{
      #'   \item For deterministic runs: returns the fixed lag value from metadata
      #'   \item For Monte Carlo runs: returns iteration-specific lag values
      #'   \item Lag affects when exposure changes influence disease incidence
      #'   \item Important for policy impact timing and lifecycle analysis
      #' }
      #'
      #' @examples
      #' \dontrun{
      #' # Get median lag
      #' median_lag <- bmi_chd$get_lag()
      #'
      #' # Get lag for specific Monte Carlo iterations
      #' mc_lags <- bmi_chd$get_lag(c(1, 5, 10))
      #'
      #' # Use in conditional logic
      #' if (smoking_copd$get_lag() > 10) {
      #'   message("Long lag exposure")
      #' }
      #' }
      get_lag = function(mc_) {
        if (missing(mc_) || mc_ == 0) {
          self$lag
        } else {
          private$lag_mc[mc_]
        }
      },

      # get_ideal_xps_lvl ----
      #' @description
      #' Get the ideal (counterfactual) exposure level used for population attributable
      #' fraction calculations. This represents the theoretical minimum risk exposure
      #' level for calculating preventable disease burden.
      #'
      #' @param mc_ Integer vector. Monte Carlo iteration numbers. If missing or 0,
      #'   returns the mean ideal exposure level. Otherwise returns iteration-specific
      #'   ideal exposure levels for stochastic calculations.
      #'
      #' @return Numeric vector containing ideal exposure level(s). Length matches
      #'   the length of \code{mc_} parameter, or single value if \code{mc_} is missing/0.
      #'
      #' @details
      #' The ideal exposure level is used for:
      #' \itemize{
      #'   \item Population Attributable Risk Fraction (PARF) calculations
      #'   \item Determining preventable disease burden
      #'   \item Setting counterfactual scenarios in policy simulations
      #'   \item Calculating theoretical minimum risk exposure distributions
      #' }
      #'
      #' Common ideal levels include:
      #' \itemize{
      #'   \item BMI: 21-23 kg/m² (optimal weight range)
      #'   \item Smoking: 0 (complete cessation)
      #'   \item Salt intake: WHO recommended levels
      #'   \item Physical activity: WHO recommended levels
      #' }
      #'
      #' @examples
      #' \dontrun{
      #' # Get mean ideal exposure level
      #' ideal_bmi <- bmi_chd$get_ideal_xps_lvl()
      #'
      #' # Get ideal levels for Monte Carlo iterations
      #' ideal_levels <- bmi_chd$get_ideal_xps_lvl(c(1, 5, 10))
      #'
      #' # Use in PARF calculations
      #' if (current_exposure > smoking_copd$get_ideal_xps_lvl()) {
      #'   # Calculate attributable risk
      #' }
      #' }
      get_ideal_xps_lvl = function(mc_) {
        if (missing(mc_) || mc_ == 0) {
          mean(private$ideal_xps_lvl_mc)
        } else {
          private$ideal_xps_lvl_mc[mc_]
        }
      },

      # get_name ----
      #' @description
      #' Get the internal suffix name used for file naming and identification.
      #' This is typically a combination of the exposure name and outcome,
      #' used for creating unique filenames and internal references.
      #'
      #' @return Character string containing the internal name/suffix used for
      #'   file identification and naming conventions.
      #'
      #' @details
      #' The internal name is used for:
      #' \itemize{
      #'   \item Creating unique filenames for stochastic effect files
      #'   \item Internal object identification and referencing
      #'   \item Avoiding naming conflicts between exposure-outcome pairs
      #'   \item File management and organization
      #' }
      #'
      #' Format typically follows: "exposure_outcome" pattern
      #' (e.g., "bmi_chd", "smoking_copd", "salt_stroke")
      #'
      #' @examples
      #' \dontrun{
      #' # Get internal name for file identification
      #' internal_name <- bmi_chd$get_name()
      #'
      #' # Use in file path construction
      #' file_path <- paste0("results/", smoking_copd$get_name(), "_output.csv")
      #'
      #' # Compare internal names
      #' if (exp1$get_name() == exp2$get_name()) {
      #'   warning("Duplicate exposure-outcome pairs detected")
      #' }
      #' }
      get_name = function() {
        private$suffix
      },

      # print ----
      #' @description
      #' Print a summary of the exposure object showing key parameters and metadata.
      #' Displays exposure name, outcome, distribution type, lag, source, and notes
      #' in a readable format.
      #'
      #' @return The \code{Exposure} object, invisibly, for method chaining.
      #'
      #' @details
      #' Prints formatted information including:
      #' \itemize{
      #'   \item Exposure name (e.g., "bmi", "smoking")
      #'   \item Disease outcome influenced
      #'   \item Uncertainty distribution type
      #'   \item Median lag in years
      #'   \item Source citation
      #'   \item Additional notes
      #' }
      #'
      #' @examples
      #' \dontrun{
      #' # Print exposure information
      #' bmi_chd$print()
      #'
      #' # Print in workflow
      #' smoking_copd$gen_stochastic_effect(design)$print()
      #' }
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

      # Write Exposure Parameter File
      # @description
      # Writes an exposure parameter file in CSVY format (CSV with YAML metadata header).
      # This method combines metadata and relative risk data into a single file that can
      # be read by the \code{initialize} method.
      #
      # @param dt A data.table containing relative risk data with columns for agegroup,
      #   sex, rr (relative risk), ci_rr (confidence interval), and optional stratification variables.
      # @param metadata A list containing metadata fields such as exposure name, outcome,
      #   distribution, lag, source, and notes.
      # @param file_path Character string. The complete file path where the CSVY file
      #   will be written, including filename and .csvy extension.
      #
      # @return \code{NULL}. The method writes to disk as a side effect.
      #
      # @details
      # The method creates a CSVY file by:
      # \itemize{
      #   \item Converting metadata to YAML format with comment prefixes
      #   \item Writing YAML header to the file
      #   \item Appending CSV data using \code{fwrite}
      # }
      #
      # @keywords internal
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

      # Generate Stochastic Relative Risk Table
      # @description
      # Generates a single stochastic relative risk value by sampling from a specified
      # probability distribution using the central estimate and confidence interval.
      #
      # @param m Numeric. Central estimate of the relative risk (point estimate).
      # @param ci Numeric. Confidence interval bound of the relative risk (upper or lower).
      # @param distribution Character. Distribution type for stochastic sampling:
      #   \code{"lognormal"} or \code{"normal"}.
      #
      # @return Numeric. A single stochastic relative risk value sampled from the distribution.
      #
      # @details
      # For lognormal distribution:
      # \itemize{
      #   \item Assumes log(RR) follows normal distribution
      #   \item Uses log transformation for sampling
      #   \item Appropriate for relative risks (always positive)
      # }
      #
      # For normal distribution:
      # \itemize{
      #   \item Samples directly from normal distribution
      #   \item May produce negative values (use with caution)
      # }
      #
      # Uses \code{dqRNG} for reproducible random number generation.
      #
      # @keywords internal
      # need to run by id
      stochRRtabl = function(m, ci, distribution = c("lognormal", "normal")) {
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
      # Generate Stochastic Relative Risk List
      # @description
      # Core method for generating age-specific stochastic relative risks through
      # Monte Carlo sampling, linear interpolation, and optional smoothing.
      # This method transforms age-grouped relative risks into single-year estimates
      # with uncertainty quantification.
      #
      # @param dt A data.table containing age-grouped relative risk information
      #   with columns for agegroup, sex, rr (relative risk), and ci_rr (confidence interval).
      # @param tt A data.table containing single-year age information to interpolate to.
      # @param mc_max Integer. Maximum number of Monte Carlo iterations to generate.
      # @param smooth Logical. Whether to apply loess smoothing to the interpolated
      #   relative risks across age groups.
      # @param do_checks Logical. Whether to perform validation checks on generated
      #   relative risks (default FALSE for performance).
      # @param ... Additional arguments passed to the loess smoothing function
      #   (e.g., span parameter for smoothing bandwidth).
      #
      # @return A data.table with generated relative risks containing columns for
      #   mc (Monte Carlo iteration), age (single years), sex, exposure-specific RR column,
      #   and any stratification variables from the input data.
      #
      # @details
      # The generation process involves:
      # \itemize{
      #   \item \strong{Constraint detection}: Determines if RRs should be ≥1, ≤1, or mixed
      #   \item \strong{Monte Carlo expansion}: Clones data for each MC iteration
      #   \item \strong{Stochastic sampling}: Uses \code{stochRRtabl()} for distribution sampling
      #   \item \strong{Age interpolation}: Linear interpolation from age groups to single years
      #   \item \strong{Exposure interpolation}: Handles non-consecutive exposure levels
      #   \item \strong{Optional smoothing}: Applies loess smoothing across age if requested
      #   \item \strong{Constraint enforcement}: Ensures RRs respect original constraints
      #   \item \strong{Validation}: Optional checks for NAs and constraint violations
      # }
      #
      # Interpolation handles:
      # \itemize{
      #   \item Age groups to single-year ages using linear interpolation
      #   \item Non-consecutive exposure levels (e.g., BMI categories 1,3,4 → 1,2,3,4)
      #   \item Multiple stratification variables (sex, smoking status, etc.)
      # }
      #
      # @keywords internal
      generate_rr_l = function(dt, tt, mc_max, smooth, do_checks = FALSE, ...) {
        # dt -> agegrouped, tt -> age in years
        dt <- tt[dt, on = .NATURAL, nomatch = NULL, mult = "first"][,
          age := NULL
        ]
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
        dt[,
          (colnam) := private$stochRRtabl(rr, ci_rr, self$distribution),
          keyby = mc
        ]
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
          if (!isTRUE(all.equal(t3, t4))) {
            # non consecutive
            tt <-
              CJ(
                mc = 1:mc_max,
                age = min(tt$age):max(tt$age),
                sex = factor(c("men", "women")),
                V3 = t3
              )

            replace_from_table(
              tt,
              colname = "V3",
              from = t3,
              to = rep(t4, times = c(diff(t4), 1L)),
              newcolname = nam
            )

            dt <- tt[dt, on = .NATURAL]
            dt[V3 != get(nam), (colnam) := NA_real_]
            dt[,
              (colnam) := approx(get(nam), get(colnam), V3)$y,
              by = .(mc, age, sex)
            ]
            dt[, (nam) := NULL]
            setnames(dt, "V3", nam)
          }
        }

        if (smooth) {
          dt[,
            (colnam) := predict(loess(
              as.formula(paste0(colnam, " ~ age")),
              .SD,
              ...
            )),
            by = eval(colby)
          ]
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
          } else message("RR crosses 1")
        }
        dt
      },

      # ci_from_p_forratio ----
      # Calculate 95% Confidence Interval from P-value for Ratios
      # @description
      # Converts a p-value and mean relative risk to approximate 95% confidence intervals
      # using the relationship between z-statistics and p-values. This is useful when
      # only p-values are available from published studies rather than explicit CIs.
      #
      # @param mean_rr Numeric. The mean (point estimate) relative risk.
      # @param p Numeric. The p-value associated with the relative risk estimate.
      #
      # @return A list containing:
      # \itemize{
      #   \item \code{lower_ci} - Lower bound of the 95% confidence interval
      #   \item \code{upper_ci} - Upper bound of the 95% confidence interval
      # }
      #
      # @details
      # The conversion process:
      # \itemize{
      #   \item Uses the approximation: z ≈ -0.862 + √(0.743 - 2.404 × ln(p))
      #   \item Calculates standard error from z-statistic and log(RR)
      #   \item Constructs symmetric confidence interval on log scale
      #   \item Exponentiates back to RR scale for final intervals
      # }
      #
      # This approximation assumes:
      # \itemize{
      #   \item Two-sided test with normal distribution
      #   \item Log-normal distribution for relative risks
      #   \item Symmetric confidence intervals on log scale
      # }
      #
      # @keywords internal
      ci_from_p_forratio = function(mean_rr, p) {
        z <- -0.862 + sqrt(0.743 - 2.404 * log(p))
        est <- log(mean_rr)
        se <- abs(est / z)
        dnm <- qnorm(0.975, 0, 1) # ~1.96
        ci <- list()
        ci$lower_ci <- exp(est - dnm * se)
        ci$upper_ci <- exp(est + dnm * se)
        ci
      }
    ) # end of private
  )
