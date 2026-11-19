## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2026 University of Liverpool, Chris Kypridemos
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

#' @title Exposure Variable Generator
#' @name Exposure
#'
#' @description
#' An R6 class for generating exposure variables in synthetic populations.
#' This class encapsulates the repetitive logic of reading distribution files,
#' joining with population data, generating exposure values from quantile functions,
#' and managing temporary column cleanup.
#'
#' @details
#' The `Exposure` class standardizes exposure variable generation by:
#' \itemize{
#'   \item Reading exposure distribution files from standardized location
#'   \item Managing join operations (lookup_dt, absorb_dt)
#'   \item Applying quantile functions with appropriate parameters
#'   \item Handling distribution-specific generation logic
#'   \item Handling conditional removal of rank variables
#'   \item Supporting custom post-generation transformations
#'   \item Providing consistent logging
#' }
#'
#' @section Distribution Types:
#' \describe{
#'   \item{ordinal}{Discrete ordered categories using threshold comparisons}
#'   \item{continuous}{Continuous distributions using quantile functions}
#'   \item{binary}{Binary (0/1) outcomes from probability thresholds}
#'   \item{custom}{User-defined generation function}
#' }
#'
#' @section Usage Pattern:
#' ```r
#' # Ordinal exposure with thresholds
#' # Note: thresholds auto-detected from parquet columns not in synthpop
#' Exposure$new("income", "income",
#'              var_name = "income", rank_var = "rank_income",
#'              distribution = "ordinal")$
#'   generate(pop, design)$
#'   factorise(levels = 1:5, labels = c("1 Highest", "2", "3", "4", "5 Lowest"))
#'
#' # Continuous exposure with quantile function
#' # Note: qparams auto-detected from parquet columns (mu, sigma, nu, tau)
#' Exposure$new("BMI", "bmi",
#'              var_name = "bmi", rank_var = "rank_bmi",
#'              distribution = "continuous",
#'              qfun = "my_qBCPEo",
#'              qargs = list(n_cpu = 1L))$
#'   generate(pop, design)
#'
#' # Binary exposure
#' Exposure$new("ETS", "ets",
#'              var_name = "ets", rank_var = "rank_ets",
#'              distribution = "binary",
#'              invert = TRUE)$  # rank_ets >= (1 - mu)
#'   generate(pop, design)
#' ```
#'
#' @export
Exposure <-
  R6::R6Class(
    classname = "Exposure",
    # public ----
    public = list(
      #' @field name The name of the exposure variable being generated (for logging)
      name = NULL,

      #' @field file_name The identifier of the partitioned parquet dataset
      #'   (bare name, YAML-relative path, or absolute path).
      file_name = NULL,

      #' @field file_path Full (resolved) path to the parquet dataset directory
      file_path = NULL,

      #' @field var_name The name of the variable to create in the data.table
      var_name = NULL,

      #' @field rank_var The name of the rank variable
      rank_var = NULL,

      #' @field distribution Type of distribution: "ordinal", "continuous", "binary", or "custom"
      distribution = NULL,

      #' @field qfun Quantile function name for continuous distributions
      qfun = NULL,

      #' @field pfun Probability function name (for generating qmin/qmax)
      pfun = NULL,

      #' @field min_value Lowest value a draw may take, INCLUSIVE: it becomes
      #'   `minq = P(Y < min_value)` (see `add_qmin_qmax()`), so a point mass
      #'   at `min_value` (e.g. the zeros of a count) is drawn with its
      #'   probability. Before 2026-11-18 it was `P(Y <= min_value)`, which
      #'   made the lower bound exclusive.
      min_value = NULL,

      #' @field max_value Highest value a draw may take, inclusive:
      #'   `maxq = P(Y <= max_value)`.
      max_value = NULL,

      #' @field family_type Support of the distribution, as gamlss.dist's
      #'   family `type`: "Discrete", "Continuous" or "Mixed". It decides how
      #'   `min_value` becomes `minq`. Inferred from `qfun` (e.g. `fqZINBI` ->
      #'   ZINBI -> "Discrete") unless declared (YAML `family_type:`, any
      #'   case). NULL when neither gives one.
      family_type = NULL,

      #' @field qparams Parameter names for quantile function
      qparams = NULL,

      #' @field qargs Additional arguments for quantile function
      qargs = NULL,

      #' @field thresholds Threshold variable names for ordinal distributions
      thresholds = NULL,

      #' @field invert For binary: if TRUE uses (1-mu), if FALSE uses mu as threshold
      invert = FALSE,

      #' @field lower_tail For binary: if TRUE uses rank < threshold, if FALSE uses rank > threshold
      lower_tail = FALSE,

      #' @field custom_fn Custom generation function
      custom_fn = NULL,

      #' @field transform_fn Transformation multiplier or function
      transform_fn = NULL,

      #' @field offset Numeric offset to add after transform_fn (e.g., -1 to shift 1:8 to 0:7)
      offset = NULL,

      #' @field additional_rank_vars Additional rank variables to remove
      additional_rank_vars = NULL,

      #' @field join_fn Join function: "lookup_dt" or "absorb_dt"
      join_fn = "lookup_dt",

      #' @field factorise_params Parameters for optional factorisation (levels, labels)
      factorise_params = NULL,

      #' @field cast_to_int Whether to cast result to integer
      cast_to_int = FALSE,

      #' @field invert_ratio Whether to invert ratio (1/x)
      invert_ratio = FALSE,

      #' @description
      #' Initialize a new Exposure generator
      #'
      #' @param name Character. Name of the exposure variable (for logging)
      #' @param file_name Character. Identifier of the partitioned parquet
      #'   dataset. Typically a bare directory name resolved under
      #'   inputs/exposure_distributions/, but may also be a YAML-relative or
      #'   absolute path (see `Design`'s "Exposure table resolution" section).
      #'   Retained here mainly for logging/lookup; the concrete location is
      #'   `file_path`.
      #' @param var_name Character. Variable name to create in data.table
      #' @param rank_var Character. Name of rank variable
      #' @param distribution Character. Distribution type: "ordinal", "continuous", "binary", or "custom"
      #' @param qfun Character. Quantile function name for continuous distributions
      #' @param pfun Character. Probability function name (derived from qfun if NULL, e.g., my_qBCPEo -> my_pBCPEo)
      #' @param min_value Numeric. Lowest value a draw may take, inclusive
      #'   (e.g. 14 for BMI; 0 for a count whose zeros are real, such as
      #'   alcohol). See the `min_value` field.
      #' @param max_value Numeric. Highest value a draw may take, inclusive
      #'   (e.g. 70 for BMI).
      #' @param family_type Character. "Discrete", "Continuous" or "Mixed"
      #'   (any case). NULL (default) infers it from `qfun` through
      #'   gamlss.dist; declare it only for a qfun that is not a gamlss.dist
      #'   family. Required, one way or the other, when both `min_value` and
      #'   `max_value` are set.
      #' @param qparams Character vector. Parameter names for quantile function (auto-detected from parquet if NULL)
      #' @param qargs List. Additional arguments for quantile function
      #' @param thresholds Character vector. Threshold names for ordinal distributions (auto-detected from parquet if NULL)
      #' @param invert Logical. For binary: use (1 - mu) instead of mu
      #' @param lower_tail Logical. For binary: if TRUE uses rank < threshold, if FALSE uses rank > threshold
      #' @param custom_fn Function. Custom generation function for complex cases
      #' @param transform_fn Numeric or function. Multiplier or transformation function
      #' @param offset Numeric. Value to add after transform_fn (e.g., -1 to shift 1:8 to 0:7)
      #' @param additional_rank_vars Character vector. Additional rank vars to remove
      #' @param join_fn Character. Join function: "lookup_dt" or "absorb_dt"
      #' @param base_path Character. Base path to exposure distributions
      #'   directory. Only used when `file_path` is NULL and `file_name` is not
      #'   an absolute path (i.e. the fall-back for direct callers).
      #' @param file_path Character. Pre-resolved full path to the parquet
      #'   dataset. When supplied (as `Design$load_exposures()` does) it is used
      #'   verbatim and takes precedence over `base_path`/`file_name`, keeping
      #'   path-resolution policy in one place (`Design`). Leave NULL to fall
      #'   back to combining `base_path` with `file_name`.
      #'
      #' @return A new Exposure object
      #'
      #' @examples
      #' # Ordinal
      #' income_exp <- Exposure$new("income", "income",
      #'                            var_name = "income", rank_var = "rank_income",
      #'                            distribution = "ordinal",
      #'                            thresholds = c("inc1", "inc2", "inc3", "inc4"))
      #' # Continuous
      #' bmi_exp <- Exposure$new("BMI", "bmi",
      #'                         var_name = "bmi", rank_var = "rank_bmi",
      #'                         distribution = "continuous",
      #'                         qfun = "my_qBCPEo",
      #'                         qparams = c("mu", "sigma", "nu", "tau"))
      # initialize ----
      initialize = function(
        name,
        file_name,
        var_name = NULL,
        rank_var = NULL,
        distribution = c("ordinal", "continuous", "binary", "custom"),
        qfun = NULL,
        pfun = NULL,
        min_value = NULL,
        max_value = NULL,
        qparams = NULL,
        qargs = NULL,
        thresholds = NULL,
        invert = FALSE,
        lower_tail = FALSE,
        custom_fn = NULL,
        transform_fn = NULL,
        offset = NULL,
        additional_rank_vars = NULL,
        join_fn = c("lookup_dt", "absorb_dt"),
        base_path = "./inputs/exposure_distributions/",
        file_path = NULL,
        family_type = NULL
      ) {
        self$name <- name
        self$file_name <- file_name
        # Path policy lives in Design$resolve_exposure_path(): when it hands us
        # a pre-resolved file_path, use it verbatim. Otherwise (direct callers,
        # examples) fall back to the historical rule: absolute file_name as-is,
        # else combined with base_path.
        if (!is.null(file_path)) {
          self$file_path <- file_path
        } else if (grepl("^/", file_name) || grepl("^[A-Za-z]:", file_name)) {
          self$file_path <- file_name
        } else {
          self$file_path <- file.path(base_path, file_name)
        }

        # Pre-read column names from parquet for efficiency
        # This avoids reading the full dataset just to get column names
        private$read_col_names()

        self$var_name <- var_name
        self$rank_var <- rank_var
        self$distribution <- match.arg(distribution)
        self$qfun <- qfun

        # Auto-derive pfun from qfun if not provided (qZINBI -> pZINBI,
        # my_qBCPEo -> my_pBCPEo, fqZINBI -> fpZINBI)
        self$pfun <- .exposure_pfun(qfun, pfun)

        self$min_value <- min_value
        self$max_value <- max_value
        # The support type decides how min_value becomes minq (inclusive
        # bounds, .exposure_cdf_bounds()). Without one, bounds cannot be
        # applied correctly, so refuse here rather than draw wrongly later.
        self$family_type <- .exposure_family_type(qfun, family_type)
        if (!is.null(min_value) && !is.null(max_value) &&
            is.null(self$family_type)) {
          stop(
            "Exposure '", name, "': cannot tell whether qfun '",
            format(qfun), "' is discrete, continuous or mixed, so its ",
            "min_value cannot be applied. Declare family_type ",
            "(discrete / continuous / mixed) in its exposure definition.",
            call. = FALSE
          )
        }
        self$qparams <- qparams
        self$qargs <- qargs
        self$thresholds <- thresholds
        self$invert <- invert
        self$lower_tail <- lower_tail
        self$custom_fn <- custom_fn
        self$transform_fn <- transform_fn
        self$offset <- offset
        self$additional_rank_vars <- additional_rank_vars
        self$join_fn <- match.arg(join_fn)
      },

      #' @description
      #' Read the exposure distribution table from disk.
      #'
      #' @param ... Additional arguments passed to read_parquet_dt()
      #'   (e.g., filter, columns, as_data_table).
      #'
      #' @return The object returned by read_parquet_dt().
      # get_table ----
      get_table = function(...) {
        read_parquet_dt(self$file_path, ...)
      },

      #' @description
      #' Generate exposure variable in population data.table
      #'
      #' @param pop data.table. The synthetic population data
      #' @param idx integer vector. Select rows in pop to generate exposure for (optional)
      #' @param design Design object. The simulation design
      #' @return Self (invisibly) for method chaining
      #'
      #' @examples
      #' # Ordinal exposure
      #' Exposure$new("income", "income",
      #'              var_name = "income", rank_var = "rank_income",
      #'              distribution = "ordinal",
      #'              thresholds = c("inc1", "inc2", "inc3", "inc4"))$
      #'   generate(pop, design)
      #'
      #' # Continuous exposure
      #' Exposure$new("BMI", "bmi",
      #'              var_name = "bmi", rank_var = "rank_bmi",
      #'              distribution = "continuous",
      #'              qfun = "my_qBCPEo",
      #'              qparams = c("mu", "sigma", "nu", "tau"))$
      #'   generate(pop, design)
      # generate ----
      generate = function(pop, design, idx) {
        # Log message
        if (design$sim_prm$logs) {
          message(paste0("Generate ", self$name, " exposure variable..."))
        }

        # Read distribution table from partitioned parquet
        tbl <- read_parquet_dt(self$file_path)

        # Add qmin/qmax columns if needed (for continuous distributions)
        if (self$distribution == "continuous") {
          self$add_qmin_qmax(tbl)
        }

        # Auto-detect qparams using cached column names if not provided (for continuous distributions)
        if (self$distribution == "continuous" && is.null(self$qparams)) {
          # Common parameter names in order of priority
          param_candidates <- c("mu", "sigma", "nu", "tau", "alpha", "beta")
          detected_params <- intersect(param_candidates, private$col_names)
          if (length(detected_params) > 0) {
            self$qparams <- detected_params
          }
        }

        # Auto-detect thresholds using cached column names if not provided (for ordinal distributions)
        if (self$distribution == "ordinal" && is.null(self$thresholds)) {
          # Common join keys and parameters to exclude
          exclude_cols <- c(
            "year",
            "age",
            "sex",
            "qimd",
            "dimd",
            "ethnicity",
            "sha",
            "education",
            "LAD17CD",
            "LAD17NM",
            "lsoa",
            "mc",
            "pid",
            "mu",
            "sigma",
            "nu",
            "tau",
            "alpha",
            "beta",
            # Quantile bounds are not thresholds. Tables written by
            # add_qmin_qmax(write_to_disk = TRUE) carry a constant minq = 0 /
            # maxq = 1 even for ordinal exposures; picked up as thresholds,
            # `rank > minq` is always TRUE and shifted every value by +1
            # (active_days was 1..8 instead of 0..7 until 2026-11-18).
            "minq",
            "maxq"
          )

          # Find columns in table that are not in population and not common join keys
          cols_in_tbl_only <- setdiff(
            private$col_names,
            c(names(pop), exclude_cols)
          )

          # These should be the threshold columns
          if (length(cols_in_tbl_only) > 0) {
            # Sort to ensure consistent ordering (e.g., inc1, inc2, inc3, inc4 or pa0, pa1, ...)
            self$thresholds <- sort(cols_in_tbl_only)
          }
        }

        kc <- sort(setdiff(
          names(tbl),
          c("mu", "sigma", "nu", "tau", "maxq", "minq", self$thresholds)
        ))
        kc <- kc[order(match(kc, "year"))]
        setcolorder(tbl, kc)
        setkeyv(tbl, kc)

        # Calculate columns to remove (columns in tbl that are not already in pop)
        # Use names(tbl) rather than private$col_names so runtime-added columns
        # (e.g. maxq/minq from add_qmin_qmax) are also cleaned up
        col_nam <- setdiff(
          names(tbl),
          intersect(names(pop), names(tbl))
        )

        # Harmonize column types between tbl and pop for join columns
        # This fixes issues where parquet stores factors but pop has integers
        join_cols <- intersect(names(tbl), names(pop))
        for (jc in join_cols) {
          tbl_type <- typeof(tbl[[jc]])
          pop_type <- typeof(pop[[jc]])
          tbl_is_factor <- is.factor(tbl[[jc]])
          pop_is_factor <- is.factor(pop[[jc]])

          if (tbl_is_factor && !pop_is_factor && pop_type == "integer") {
            # Convert factor to integer to match pop
            tbl[, (jc) := as.integer(as.character(get(jc)))]
          } else if (!tbl_is_factor && pop_is_factor && tbl_type == "integer") {
            # Convert integer to factor to match pop
            tbl[, (jc) := factor(get(jc), levels = levels(pop[[jc]]))]
          }
        }

        # Re-set key after type conversion
        setkeyv(tbl, kc)

        # Perform join
        if (self$join_fn == "lookup_dt") {
          lookup_dt(pop, tbl, check_lookup_tbl_validity = design$sim_prm$logs)
        } else {
          absorb_dt(pop, tbl)
        }

        # Generate based on distribution type
        if (self$distribution == "ordinal") {
          private$generate_ordinal(pop, tbl, idx)
        } else if (self$distribution == "continuous") {
          private$generate_continuous(pop, idx)
        } else if (self$distribution == "binary") {
          private$generate_binary(pop, idx)
        } else if (self$distribution == "custom") {
          private$generate_custom(pop, tbl, idx)
        }

        # Apply post-generation transformations automatically
        private$apply_transformations(pop)

        # Conditionally add rank variable to removal list
        if (!is.null(self$rank_var) && !design$sim_prm$keep_simulants_rn) {
          col_nam <- c(col_nam, self$rank_var)
        }

        # Add any additional rank variables
        if (!is.null(self$additional_rank_vars)) {
          col_nam <- c(col_nam, self$additional_rank_vars)
        }

        # Remove temporary columns
        pop[, (col_nam) := NULL]

        invisible(self)
      },

      #' @description
      #' Apply post-generation transformation (factorisation, scaling, etc.)
      #'
      #' @param pop data.table. The synthetic population data to transform
      #' @param transform_fn Function. Function that takes pop as argument and transforms it
      #'
      #' @return Self (invisibly) for method chaining
      #'
      #' @examples
      #' Exposure$new(...)$generate(pop, design)$transform(pop, function(p) {
      #'     p[, income := factor(income, levels = 1:5)]
      #'   })
      # transform ----
      transform = function(pop, transform_fn) {
        transform_fn(pop)
        invisible(self)
      },

      #' @description
      #' Factorise the generated variable (convenience method)
      #'
      #' @param pop data.table. The synthetic population data
      #' @param levels Levels for the factor
      #' @param labels Labels for the factor
      #'
      #' @return Self (invisibly) for method chaining
      # factorise ----
      factorise = function(pop, levels, labels) {
        pop[,
          (self$var_name) := factor(.var, levels = levels, labels = labels),
          env = list(.var = self$var_name)
        ]
        invisible(self)
      },

      #' @description
      #' Add qmin and qmax columns to distribution table if missing or incorrect.
      #' Validates existing values by sampling 0.01% of rows and comparing with expected values.
      #' Can optionally write the updated table back to disk as partitioned parquet.
      #'
      #' The bounds truncate draws to `[min_value, max_value]`, BOTH ends
      #' inclusive: `minq = P(Y < min_value)` and `maxq = P(Y <= max_value)`,
      #' by the exposure's `family_type` (see `.exposure_cdf_bounds()`). A
      #' table whose stored `minq` follows the rule used before 2026-11-18
      #' (`P(Y <= min_value)`, which never drew a point mass at `min_value`)
      #' fails the check and is recalculated in memory. The PARF generator
      #' (`Disease$gen_sp_forPARF()`) calls this method too, so both paths
      #' draw within the same bounds.
      #'
      #' @param tbl data.table. The distribution table
      #' @param force Logical. If TRUE, recalculate even if columns exist and are valid
      #' @param write_to_disk Logical. If TRUE, write updated table back to disk as partitioned parquet
      #' @param partition_cols Character vector. Column names to use for partitioning when writing
      #' @param tolerance Numeric. Tolerance for comparing existing vs expected values (default 1e-8)
      #' @param sample_idx Integer vector or NULL. Rows on which the stored
      #'   bounds are checked. NULL (default) draws 0.01% of the rows with
      #'   `sample.int()`, i.e. from the ambient base-R stream, as the
      #'   synthetic-population path always has. The PARF path passes evenly
      #'   spaced rows (`.even_sample_idx()`) so that it consumes no random
      #'   numbers.
      #' @param quiet Logical. If TRUE, stale bounds are recalculated without
      #'   the (once-per-object) message, and without marking it as given.
      #'
      #' @return Self (invisibly) for method chaining
      # add_qmin_qmax ----
      add_qmin_qmax = function(
        tbl,
        force = FALSE,
        write_to_disk = FALSE,
        partition_cols = "year",
        tolerance = 1e-8,
        sample_idx = NULL,
        quiet = FALSE
      ) {
        needs_recalc <- force

        # Auto-detect qparams from table if not already set
        if (is.null(self$qparams)) {
          param_candidates <- c("mu", "sigma", "nu", "tau", "alpha", "beta")
          detected_params <- intersect(param_candidates, names(tbl))
          if (length(detected_params) > 0) {
            self$qparams <- detected_params
          }
        }

        # Check if columns exist
        has_minq <- "minq" %in% names(tbl)
        has_maxq <- "maxq" %in% names(tbl)

        if (!has_minq || !has_maxq) {
          needs_recalc <- TRUE
        } else if (
          !force && !is.null(self$min_value) && !is.null(self$max_value)
        ) {
          # Validate existing values on 0.01% of rows. By default a random
          # sample from the ambient base-R stream (kept as it was, so the
          # synthetic population's stream is unchanged); the PARF path hands
          # in evenly spaced rows instead.
          n_rows <- nrow(tbl)
          if (is.null(sample_idx)) {
            sample_size <- max(1L, as.integer(ceiling(n_rows * 0.0001)))
            sample_idx <- sample.int(n_rows, sample_size)
          }
          sample_size <- length(sample_idx)

          # Build parameter list for pfun using sampled rows
          params_list_sample <- lapply(self$qparams, function(p) {
            tbl[[p]][sample_idx]
          })
          names(params_list_sample) <- self$qparams

          # Expected bounds for the sample, by the inclusive rule
          expected <- private$cdf_bounds(params_list_sample, sample_size)

          # Compare with existing values
          existing_maxq <- tbl[["maxq"]][sample_idx]
          existing_minq <- tbl[["minq"]][sample_idx]

          # Check if values match within tolerance
          maxq_mismatch <- any(
            abs(existing_maxq - expected$maxq) > tolerance,
            na.rm = TRUE
          )
          minq_mismatch <- any(
            abs(existing_minq - expected$minq) > tolerance,
            na.rm = TRUE
          )

          if (maxq_mismatch || minq_mismatch) {
            # Once per object: a stale table is met again by every synthpop
            # chunk and every PARF disease in the same process.
            if (!isTRUE(quiet) && !isTRUE(private$minq_msg_done)) {
              message(
                "Existing minq/maxq values do not match expected values for ",
                self$name,
                " (inclusive bounds: minq = P(Y < min_value), maxq = ",
                "P(Y <= max_value)). Recalculating in memory..."
              )
              private$minq_msg_done <- TRUE
            }
            needs_recalc <- TRUE
          }
        }

        if (needs_recalc) {
          # If min_value and max_value are specified, calculate from distribution
          if (!is.null(self$min_value) && !is.null(self$max_value)) {
            # Build parameter list for pfun
            params_list <- lapply(self$qparams, function(p) tbl[[p]])
            names(params_list) <- self$qparams

            bounds <- private$cdf_bounds(params_list, nrow(tbl))
            tbl[, maxq := bounds$maxq]
            tbl[, minq := bounds$minq]
          } else {
            # Default: use full quantile range [0, 1]
            tbl[, minq := 0]
            tbl[, maxq := 1]
          }

          # Write back to disk if requested
          if (write_to_disk) {
            message("Writing updated table with minq/maxq to: ", self$file_path)

            # Ensure partition columns exist in table
            valid_partition_cols <- intersect(partition_cols, names(tbl))
            if (length(valid_partition_cols) == 0) {
              warning(
                "No valid partition columns found. Writing without partitioning."
              )
              valid_partition_cols <- NULL
            }
            
            # Write using CKutils
            if (is.null(valid_partition_cols)) {
              # Write as single file when no partitioning
              CKutils::write_parquet_dt(
                tbl,
                file.path(self$file_path, "part-0.parquet")
              )
            } else {
              # Write as partitioned dataset
              CKutils::write_parquet_dt(
                tbl,
                self$file_path,
                partition = valid_partition_cols
              )
            }

            # Update cached column names to include minq/maxq
            private$read_col_names()
          }
        }

        invisible(self)
      },

      #' @description
      #' Get cached column names from the parquet dataset
      #'
      #' @return Character vector of column names
      # get_col_names ----
      get_col_names = function() {
        private$col_names
      }
    ),

    # private ----
    private = list(
      # Cached column names from parquet dataset (read once on initialization)
      col_names = NULL,

      # add_qmin_qmax() has already reported a stale stored minq / maxq
      minq_msg_done = FALSE,

      # cdf_bounds ----
      # minq and maxq for `n` rows of distribution parameters, both ends of
      # [min_value, max_value] inclusive: .exposure_cdf_bounds().
      cdf_bounds = function(params_list, n) {
        .exposure_cdf_bounds(params_list, n, self$pfun, self$family_type,
                             self$min_value, self$max_value, self$name)
      },

      # Read column names from parquet dataset efficiently (metadata only)
      # read_col_names ----
      read_col_names = function() {
        # Omit partitioning to let Arrow auto-detect Hive-style partitions
        if (!dir.exists(self$file_path) && !file.exists(self$file_path)) {
          stop(
            "Exposure distribution file not found: ", self$file_path,
            "\nPlease ensure the file exists or check the file_name parameter."
          )
        }
        ds <- arrow::open_dataset(self$file_path, format = "parquet")
        private$col_names <- names(ds)
        invisible(self)
      },

      # Generate ordinal variable using threshold comparisons
      # generate_ordinal ----
      generate_ordinal = function(pop, tbl, idx) {
        # Build expression using substitute2:
        # var := (rank > threshold1) + (rank > threshold2) + ... + 1L
        # Note: threshold columns are already in pop after lookup_dt/absorb_dt join
        rank_sym <- as.name(self$rank_var)
        threshold_syms <- lapply(self$thresholds, as.name)

        # Build the sum of comparisons
        comparison_calls <- lapply(threshold_syms, function(th) {
          substitute2((.rank > .th), env = list(.rank = rank_sym, .th = th))
        })

        # Reduce comparisons into sum expression
        if (length(comparison_calls) == 1L) {
          sum_expr <- comparison_calls[[1L]]
        } else {
          sum_expr <- Reduce(
            function(a, b) substitute2(.a + .b, list(.a = a, .b = b)),
            comparison_calls
          )
        }

        # Add 1L to get 1-based ordinal
        final_expr <- substitute2(.sum + 1L, list(.sum = sum_expr))

        # Execute - threshold columns are already in pop from prior join
        var_name <- self$var_name
        if (missing(idx)) {
          # Use bquote for expression splicing
          dt_call <- bquote(pop[, `:=`((.(var_name)), .(final_expr))])
          eval(dt_call)
        } else {
          dt_call <- bquote(pop[idx, `:=`((.(var_name)), .(final_expr))])
          eval(dt_call)
        }
      },

      # Generate continuous variable using quantile function
      # generate_continuous ----
      generate_continuous = function(pop, idx) {
        # Build quantile function call using substitute2
        rank_sym <- as.name(self$rank_var)
        param_syms <- lapply(self$qparams, as.name)

        # Build the quantile input: clamp(minq + rank * (maxq - minq), eps, 1-eps)
        # Clamping to (0, 1) exclusive avoids GAMLSS quantile errors at boundary
        # values, which occur when max_value truncation produces maxq ≈ 0
        q_input <- substitute2(
          pmin(pmax(minq + .rank * (maxq - minq), .eps), .one_m_eps),
          list(.rank = rank_sym, .eps = .Machine$double.eps, .one_m_eps = 1 - .Machine$double.eps)
        )

        # Build parameter list for the call
        call_args <- c(list(q_input), param_syms)

        # Add extra arguments if provided
        if (!is.null(self$qargs)) {
          for (arg_name in names(self$qargs)) {
            call_args[[arg_name]] <- self$qargs[[arg_name]]
          }
        }

        # Build the function call
        qfun_call <- as.call(c(as.name(self$qfun), call_args))

        # Execute using bquote for expression splicing
        var_name <- self$var_name
        if (missing(idx)) {
          dt_call <- bquote(pop[, `:=`((.(var_name)), .(qfun_call))])
          eval(dt_call)
        } else {
          dt_call <- bquote(pop[idx, `:=`((.(var_name)), .(qfun_call))])
          eval(dt_call)
        }
      },

      # Generate binary variable from probability threshold
      # generate_binary ----
      generate_binary = function(pop, idx) {
        rank_sym <- as.name(self$rank_var)

        # Build threshold: mu or (1 - mu)
        if (self$invert) {
          threshold <- quote(1 - mu)
        } else {
          threshold <- quote(mu)
        }

        # Build comparison: rank < threshold or rank > threshold
        if (self$lower_tail) {
          # rank < threshold - low rank gives 1 (e.g., statin)
          expr <- substitute2(
            as.integer(.rank < .threshold),
            list(.rank = rank_sym, .threshold = threshold)
          )
        } else {
          # rank > threshold - high rank gives 1 (e.g., ETS, BP med)
          expr <- substitute2(
            as.integer(.rank > .threshold),
            list(.rank = rank_sym, .threshold = threshold)
          )
        }
        # Execute using bquote for expression splicing
        var_name <- self$var_name
        if (missing(idx)) {
          dt_call <- bquote(pop[, `:=`((.(var_name)), .(expr))])
          eval(dt_call)
        } else {
          dt_call <- bquote(pop[idx, `:=`((.(var_name)), .(expr))])
          eval(dt_call)
        }
      },

      # Use custom function for complex generation
      # generate_custom ----
      generate_custom = function(pop, tbl, idx) {
        if (!is.null(self$custom_fn)) {
          self$custom_fn(pop, tbl, idx)
        } else if (!is.null(self$var_name) && "mu" %in% names(pop)) {
          # Default behavior: copy mu to var_name (for simple probability lookups)
          set(pop, j = self$var_name, value = pop[["mu"]])
        }
        # If custom_fn is NULL and no mu column, do nothing - just the join has been performed
      },

      # Apply post-generation transformations (invert_ratio, transform_fn, cast_to_int, factorise)
      # apply_transformations ----
      apply_transformations = function(pop) {
        var <- self$var_name
        var_sym <- as.name(var)

        # 1. Apply invert_ratio (1/x) if specified
        if (self$invert_ratio) {
          invert_expr <- substitute2(1 / .var, list(.var = var_sym))
          dt_call <- bquote(pop[, `:=`((.(var)), .(invert_expr))])
          eval(dt_call)
        }

        # 2. Apply transform_fn if specified
        if (!is.null(self$transform_fn)) {
          if (is.numeric(self$transform_fn)) {
            # Numeric multiplier
            mult <- self$transform_fn
            mult_expr <- substitute2(.var * .mult, list(.var = var_sym, .mult = mult))
            dt_call <- bquote(pop[, `:=`((.(var)), .(mult_expr))])
            eval(dt_call)
          } else if (is.function(self$transform_fn)) {
            # Function transformation - use set() for function objects
            set(pop, j = var, value = self$transform_fn(pop[[var]]))
          }
        }

        # 2b. Apply offset if specified (adds value, use negative to subtract)
        if (!is.null(self$offset)) {
          offs <- self$offset
          # Coerce offset to integer if the column is integer to preserve type
          if (is.integer(pop[[var]])) {
            offs <- as.integer(offs)
          }
          offs_expr <- substitute2(.var + .offs, list(.var = var_sym, .offs = offs))
          dt_call <- bquote(pop[, `:=`((.(var)), .(offs_expr))])
          eval(dt_call)
        }

        # 3. Apply cast_to_int if specified
        if (self$cast_to_int) {
          int_expr <- substitute2(as.integer(.var), list(.var = var_sym))
          dt_call <- bquote(pop[, `:=`((.(var)), .(int_expr))])
          eval(dt_call)
        }

        # 4. Apply factorise_params if specified
        if (!is.null(self$factorise_params)) {
          lvls <- self$factorise_params$levels
          lbls <- self$factorise_params$labels
          if (!is.null(lvls) && !is.null(lbls)) {
            # Use set() for factor conversion with labels
            set(pop, j = var, value = factor(pop[[var]], levels = lvls, labels = lbls))
          } else if (!is.null(lvls)) {
            set(pop, j = var, value = factor(pop[[var]], levels = lvls))
          }
        }

        invisible(NULL)
      }
    )
  )

# Exposure bounds: helpers shared by Exposure (synthetic population), the PARF
# generator (Disease$gen_sp_forPARF()) and the cache keys. ----

# .EXPOSURE_BOUNDS_RULE ----
# Version of the rule that turns min_value / max_value into minq / maxq. It
# enters the synthpop checksum (SynthPop gen_checksum()), the PARF directory
# key and the PARF intermediate .qs key (Disease), so a cache built under an
# earlier rule can never be matched. Change it whenever .exposure_cdf_bounds()
# changes what it returns.
#   (none)                 before 2026-11-18: minq = P(Y <= min_value)
#   "inclusive_2026-11-18" minq = P(Y < min_value), maxq = P(Y <= max_value)
.EXPOSURE_BOUNDS_RULE <- "inclusive_2026-11-18"

# .exposure_pfun ----
# The distribution function an exposure's minq / maxq come from: the declared
# `pfun`, else derived from `qfun` (qXXX -> pXXX, prefix_qXXX -> prefix_pXXX,
# fqXXX -> fpXXX, any other: the first q -> p). NULL without a qfun. Same
# derivation as Exposure$initialize() always used.
#' @keywords internal
#' @noRd
.exposure_pfun <- function(qfun, pfun = NULL) {
  if (is.null(qfun) || !is.null(pfun)) return(pfun)
  if (grepl("^q[A-Z]", qfun)) {
    sub("^q", "p", qfun)  # starts with qX (no prefix)
  } else if (grepl("_q[A-Z]", qfun)) {
    sub("_q", "_p", qfun) # prefix_qX
  } else {
    sub("q", "p", qfun)   # fallback: the first q (fqZINBI -> fpZINBI)
  }
}

# .exposure_cdf_bounds ----
# minq and maxq for `n` rows of distribution parameters (`params_list`, named
# by the qparams): the truncation to [min_value, max_value], BOTH ends
# inclusive,
#   minq = P(Y < min_value),   maxq = P(Y <= max_value).
#
# A draw is qfun(minq + rank * (maxq - minq)), and a quantile function returns
# the smallest y with F(y) >= p. Until 2026-11-18 minq was P(Y <= min_value),
# which made the lower end EXCLUSIVE: a point mass AT min_value was never
# drawn. Every count exposure declared `min_value: 0`, so the model had no
# non-drinkers and no one eating no fruit or veg; and because the whole rank
# range was rescaled onto (F(0), F(max)], EVERY simulant moved up the
# distribution, not only the would-be zeros. For a continuous family the two
# definitions are the same number. Ported from IMPACTncd_Engl2026, which fixed
# the same bug on 2026-10-11.
#   Continuous  F(min_value).
#   Discrete    (integer support 0, 1, ...) F(ceiling(min_value) - 1), and 0
#               below the support WITHOUT calling pfun: CKutils' fast pfuns
#               up to 0.1.35 stop on a negative q, and some truncate it toward
#               zero (fpDPO(-0.5) was P(Y <= 0)). maxq is F(floor(max_value)).
#   Mixed       (gamlss.dist's ZAGA, ZAIG, BEZI, BEINF*, BEOI: point masses at
#               0 and / or 1) 0 when min_value <= 0, else F(min_value). A
#               lower bound AT the mass at 1 is not handled; it would truncate
#               to a single value anyway.
#' @keywords internal
#' @noRd
.exposure_cdf_bounds <- function(params_list, n, pfun, family_type, min_value,
                                 max_value, name = "?") {
  typ <- family_type
  if (is.null(typ) || !typ %in% c("Discrete", "Continuous", "Mixed")) {
    stop("Exposure '", name, "' has no valid family_type (got '",
         format(typ), "'), so its min_value cannot be applied (see the ",
         "family_type field).", call. = FALSE)
  }
  pfun_obj <- if (is.function(pfun)) pfun else get(pfun)
  bounds_of <- function(pl, m) {
    pq <- function(q) do.call(pfun_obj, c(list(q = rep(q, m)), pl))
    a <- min_value
    b <- max_value
    list(
      minq = switch(typ,
        Continuous = pq(a),
        Discrete = if (ceiling(a) - 1 < 0) rep(0, m) else pq(ceiling(a) - 1),
        Mixed = if (a <= 0) rep(0, m) else pq(a)
      ),
      maxq = pq(if (identical(typ, "Discrete")) floor(b) else b)
    )
  }
  out <- bounds_of(params_list, n)
  # A valid CDF gives 0 <= minq < maxq. Rows that break it are first
  # re-evaluated ONE AT A TIME: some gamlss.dist pfuns are not vectorised over
  # a parameter whose sign picks the tail (6.1.11's pGG hands pGA a vector
  # `lower.tail`, so every row takes the first row's tail; the testing
  # design's hfss_food_weight had 5 of 72,240 rows with flipped bounds,
  # drawn rank-reversed until 2026-11-19).
  bad <- which(out$maxq < out$minq | out$maxq <= 0)
  if (length(bad)) {
    for (i in bad) {
      pl <- lapply(params_list, function(p) if (length(p) == 1L) p else p[i])
      one <- bounds_of(pl, 1L)
      out$minq[i] <- one$minq
      out$maxq[i] <- one$maxq
    }
    bad <- which(out$maxq < out$minq | out$maxq <= 0)
  }
  # Still invalid row by row: the pfun is not a CDF of the family (gamlss.dist
  # 6.1.11's pZABNB returns 0 for every q >= 1, which collapsed every
  # smok_cig_ex draw to 0 until 2026-11-18), or the bounds exclude the whole
  # distribution. Either way every draw would be wrong: refuse.
  if (length(bad)) {
    pname <- if (is.character(pfun)) pfun else "<function>"
    stop("Exposure '", name, "': its pfun '", pname,
         "' gives P(Y <= max_value) below P(Y < min_value), or 0, on ",
         length(bad), " of ", n, " rows (e.g. row ", bad[1L], ": minq ",
         format(out$minq[bad[1L]]), ", maxq ", format(out$maxq[bad[1L]]),
         "), even row by row. Either it is not a valid CDF for these ",
         "parameters or [min_value, max_value] lies outside the distribution.",
         if (grepl("ZABNB", pname)) {
           " gamlss.dist 6.1.11's pZABNB is broken: use qfun fqZABNB (CKutils)."
         } else "", call. = FALSE)
  }
  out
}

# .exposure_family_type ----
# The support of an exposure's distribution, as gamlss.dist's family `type`:
# "Discrete", "Continuous" or "Mixed" (Exposure's `family_type`, which decides
# how min_value becomes minq). A declared value wins, in any case. Otherwise
# the family is read off the qfun's name -- `qZABNB`, CKutils' fast `fqZINBI`,
# a prefixed `my_qBCPEo` -- and its constructor in gamlss.dist says. NULL when
# neither gives one.
#' @keywords internal
#' @noRd
.exposure_family_type <- function(qfun, declared = NULL) {
  types <- c("Discrete", "Continuous", "Mixed")
  if (!is.null(declared)) {
    out <- if (is.character(declared) && length(declared) == 1L) {
      types[match(tolower(declared), tolower(types))]
    } else NA_character_
    if (is.na(out)) {
      stop("family_type must be one of discrete, continuous, mixed; got '",
           paste(format(declared), collapse = ", "), "'.", call. = FALSE)
    }
    return(out)
  }
  if (!is.character(qfun) || length(qfun) != 1L || is.na(qfun)) return(NULL)
  fam <- regmatches(qfun, regexec("^(?:.*_)?f?q([A-Z][A-Za-z0-9]*)$", qfun,
                                  perl = TRUE))[[1L]]
  if (length(fam) != 2L) return(NULL)
  ctor <- get0(fam[2L], envir = asNamespace("gamlss.dist"), mode = "function",
               inherits = FALSE)
  if (is.null(ctor)) return(NULL)
  out <- tryCatch(suppressWarnings(ctor())$type, error = function(e) NULL)
  if (is.character(out) && length(out) == 1L && out %in% types) out else NULL
}

# .even_sample_idx ----
# 0.01% of `n` rows (at least one), evenly spaced: the rows add_qmin_qmax()
# checks stored bounds on when it must not consume random numbers.
#' @keywords internal
#' @noRd
.even_sample_idx <- function(n) {
  if (n < 1L) return(integer(0))
  size <- max(1L, as.integer(ceiling(n * 0.0001)))
  unique(as.integer(round(seq(1, n, length.out = size))))
}

# .refresh_xps_bounds ----
# Make a distribution table read by the PARF generator carry the SAME minq /
# maxq the synthetic population draws with: the design's Exposure object
# checks the stored bounds and recalculates them in memory when they do not
# follow its rule (add_qmin_qmax(), on evenly spaced rows, so no random
# numbers are consumed). Until 2026-11-18 the PARF path used the stored
# columns as they were, whatever min_value said. `tbl` is updated by
# reference. An exposure the design does not define keeps its stored bounds,
# with a warning (the synthetic population does not draw it at all then).
#' @keywords internal
#' @noRd
.refresh_xps_bounds <- function(tbl, design_, name) {
  e <- design_$exposures[[name]]
  if (is.null(e)) {
    warning("PARF: exposure '", name, "' is not in exposure_definitions, so ",
            "its table's stored minq / maxq are used as they are.",
            call. = FALSE)
    return(invisible(tbl))
  }
  if (identical(e$distribution, "continuous")) {
    # quiet: the stock alcohol / fruit / veg tables store the pre-2026-11-18
    # minq, so the recalculation is expected, and the PARF generator meets it
    # once per forked task. `quiet` (not suppressMessages()) also leaves the
    # object's report-once flag alone, so the synthetic-population path still
    # reports a stale table even after a PARF run in the same process.
    e$add_qmin_qmax(tbl, sample_idx = .even_sample_idx(nrow(tbl)),
                    quiet = TRUE)
  }
  invisible(tbl)
}

# .PARF_EXPOSURES ----
# Every exposure whose table Disease$gen_sp_forPARF() reads (through
# design_$xps_table_path()). The PARF key hashes only their definitions, so a
# design that differs only in exposures the PARF never draws (e.g. the testing
# design's hfss / energy / salt exposures) still shares the PARFs of the stock
# design instead of evicting them. A test pins this list against the reads in
# gen_sp_forPARF().
.PARF_EXPOSURES <- c(
  "active_days", "fruit", "veg", "smok_status", "smok_quit_yrs",
  "smok_dur_ex", "smok_dur_curr", "smok_cig_curr", "smok_cig_ex", "ets",
  "alcohol", "bmi", "sbp", "tchol", "hdl_to_tchol", "statin_px"
)

# .exposure_bounds_key ----
# What the PARF's drawn exposures depend on beyond the tables: the bounds rule
# and the definitions (qfun / pfun / family_type / min_value / max_value /
# file_name ...) of the exposures the PARF generator draws (.PARF_EXPOSURES).
# The PARF directory and intermediate .qs keys include it (Disease), so a
# change in either regenerates the PARFs instead of silently reusing draws
# made under the old definitions. (The synthpop checksum hashes ALL the
# definitions -- SynthPop get_unique_characteristics() -- because the synthpop
# carries every exposure, and adds the rule in gen_checksum().)
#' The definitions are canonicalised first (entries ordered by name, fields
#' by name, numbers as doubles), so a YAML that merely reorders keys or writes
#' `1` for `1.0` keeps the same key instead of evicting the other's PARFs.
#' @keywords internal
#' @noRd
.exposure_bounds_key <- function(design_) {
  canon <- function(x) {
    if (is.list(x)) {
      if (!is.null(names(x))) x <- x[order(names(x))]
      lapply(x, canon)
    } else if (is.numeric(x)) as.double(x) else x
  }
  defs <- design_$sim_prm$exposure_definitions
  nms <- vapply(defs, function(e) if (is.null(e$name)) "" else e$name, "")
  keep <- nms %in% .PARF_EXPOSURES
  defs <- lapply(defs[keep][order(nms[keep])], canon)
  list(
    exposure_bounds_rule = .EXPOSURE_BOUNDS_RULE,
    exposure_definitions = digest(defs, serialize = TRUE)
  )
}

# .SMOKING_COUNT_EXPOSURES ----
# The smoking count exposures whose values feed simsmok() and its scenario
# variants: quit years and durations must be >= 1 for an ex- / current smoker
# there (pr_relapse is read at column quit_yrs - 1; the scenario helpers take
# smok_dur - 1), and the synthpop recodes a cigarettes-per-day 0 of any ever-
# smoker to 1. Under the inclusive bounds (2026-11-18) `min_value: 0` draws
# zeros, so a design must say min_value + offset >= 1 for them
# (.check_legacy_exposure_bounds()).
.SMOKING_COUNT_EXPOSURES <- c(
  "smok_quit_yrs", "smok_dur_ex", "smok_dur_curr", "smok_cig_curr",
  "smok_cig_ex"
)

# .check_legacy_exposure_bounds ----
# Refuse, in ONE message, the two things a design written before 2026-11-18
# (an older copy of sim_design.yaml, a Docker user's scenario folder, another
# branch) gets wrong under the inclusive bounds:
#  * a smoking count with min_value + offset < 1: `min_value: 0` then excluded
#    the zeros (minq = P(Y <= 0)) and now draws them, silently changing every
#    initial smoking history (and quit years 0 would reach simsmok(), which now
#    stops on them). `min_value: 1` gives exactly the previous bounds.
#  * a qfun whose derived pfun is gamlss.dist's pZABNB (`qfun: qZABNB`): broken
#    in 6.1.11 (0 for every q >= 1), so its bounds cannot be computed and the
#    run would stop later, inside a forked worker. `qfun: fqZABNB` (CKutils).
#' @keywords internal
#' @noRd
.check_legacy_exposure_bounds <- function(exposure_definitions) {
  low <- character(0)
  zabnb <- character(0)
  for (e in exposure_definitions) {
    if (is.null(e$name)) next
    if (e$name %in% .SMOKING_COUNT_EXPOSURES && !is.null(e$min_value)) {
      off <- if (is.null(e$offset)) 0 else as.numeric(e$offset)
      if (as.numeric(e$min_value) + off < 1) low <- c(low, e$name)
    }
    if (!is.null(e$qfun) && !is.null(e$min_value) && !is.null(e$max_value) &&
        identical(.exposure_pfun(e$qfun, e$pfun), "pZABNB")) {
      zabnb <- c(zabnb, e$name)
    }
  }
  if (length(low) || length(zabnb)) {
    stop("exposure_definitions written for the bounds rule used before ",
         "2026-11-18 (min_value is INCLUSIVE since then). Edit the design:",
         if (length(low)) paste0(
           "\n  * `min_value: 1` for ", paste(low, collapse = ", "),
           " (min_value + offset is below 1: `min_value: 0` now draws zeros; ",
           "`min_value: 1` gives exactly the previous bounds)"),
         if (length(zabnb)) paste0(
           "\n  * `qfun: fqZABNB` for ", paste(zabnb, collapse = ", "),
           " (its pfun, gamlss.dist's pZABNB, is broken in 6.1.11)"),
         "\nThe stock designs (inputs/sim_design.yaml) show both.",
         call. = FALSE)
  }
  invisible(TRUE)
}
