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

      #' @field file_name The directory name of the partitioned parquet dataset
      file_name = NULL,

      #' @field file_path Full path to the parquet dataset directory
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

      #' @field min_value Minimum value for qmin calculation
      min_value = NULL,

      #' @field max_value Maximum value for qmax calculation
      max_value = NULL,

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
      #' @param file_name Character. Directory name of the partitioned parquet dataset in inputs/exposure_distributions/
      #' @param var_name Character. Variable name to create in data.table
      #' @param rank_var Character. Name of rank variable
      #' @param distribution Character. Distribution type: "ordinal", "continuous", "binary", or "custom"
      #' @param qfun Character. Quantile function name for continuous distributions
      #' @param pfun Character. Probability function name (derived from qfun if NULL, e.g., my_qBCPEo -> my_pBCPEo)
      #' @param min_value Numeric. Minimum value for qmin calculation (e.g., 75 for BMI)
      #' @param max_value Numeric. Maximum value for qmax calculation (e.g., 250 for BMI)
      #' @param qparams Character vector. Parameter names for quantile function (auto-detected from parquet if NULL)
      #' @param qargs List. Additional arguments for quantile function
      #' @param thresholds Character vector. Threshold names for ordinal distributions (auto-detected from parquet if NULL)
      #' @param invert Logical. For binary: use (1 - mu) instead of mu
      #' @param custom_fn Function. Custom generation function for complex cases
      #' @param transform_fn Numeric or function. Multiplier or transformation function
      #' @param offset Numeric. Value to add after transform_fn (e.g., -1 to shift 1:8 to 0:7)
      #' @param additional_rank_vars Character vector. Additional rank vars to remove
      #' @param join_fn Character. Join function: "lookup_dt" or "absorb_dt"
      #' @param base_path Character. Base path to exposure distributions directory
      #'
      #' @return A new Exposure object
      #'
      #' @examples
      #' \dontrun{
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
      #' }
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
        base_path = "./inputs/exposure_distributions/"
      ) {
        self$name <- name
        self$file_name <- file_name
        # If file_name is an absolute path (starts with /), use it as is
        # Otherwise, combine with base_path
        if (grepl("^/", file_name) || grepl("^[A-Za-z]:", file_name)) {
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

        # Auto-derive pfun from qfun if not provided
        if (!is.null(qfun) && is.null(pfun)) {
          # Convert qXXX to pXXX (e.g., my_qBCPEo -> my_pBCPEo, qZINBI -> pZINBI)
          if (grepl("^q[A-Z]", qfun)) {
            # Starts with qX (no prefix)
            self$pfun <- sub("^q", "p", qfun)
          } else if (grepl("_q[A-Z]", qfun)) {
            # Has prefix_qX format
            self$pfun <- sub("_q", "_p", qfun)
          } else {
            # Fallback: just replace first q with p
            self$pfun <- sub("q", "p", qfun)
          }
        } else {
          self$pfun <- pfun
        }

        self$min_value <- min_value
        self$max_value <- max_value
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
      #' \dontrun{
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
      #' }
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
            "beta"
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

        # Calculate columns to remove (use cached col_names for efficiency)
        col_nam <- setdiff(
          private$col_names,
          intersect(names(pop), private$col_names)
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
      #' @param transform_fn Function. Function that transforms the generated exposure
      #'
      #' @return Self (invisibly) for method chaining
      #'
      #' @examples
      #' \dontrun{
      #' Exposure$new(...)$generate(pop, design)$transform(function() {
      #'     pop[, income := factor(income, levels = 1:5)]
      #'   })
      #' }
      # transform ----
      transform = function(transform_fn) {
        transform_fn()
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
      #' @param tbl data.table. The distribution table
      #' @param force Logical. If TRUE, recalculate even if columns exist and are valid
      #' @param write_to_disk Logical. If TRUE, write updated table back to disk as partitioned parquet
      #' @param partition_cols Character vector. Column names to use for partitioning when writing
      #' @param tolerance Numeric. Tolerance for comparing existing vs expected values (default 1e-8)
      #'
      #' @return Self (invisibly) for method chaining
      # add_qmin_qmax ----
      add_qmin_qmax = function(
        tbl,
        force = FALSE,
        write_to_disk = FALSE,
        partition_cols = "year",
        tolerance = 1e-8
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
          # Validate existing values by sampling 1% of rows
          n_rows <- nrow(tbl)
          sample_size <- max(1L, as.integer(ceiling(n_rows * 0.0001)))
          sample_idx <- sample.int(n_rows, sample_size)

          # Build parameter list for pfun using sampled rows
          params_list_sample <- lapply(self$qparams, function(p) {
            tbl[[p]][sample_idx]
          })
          names(params_list_sample) <- self$qparams

          # Get the probability function
          pfun_obj <- get(self$pfun)

          # Calculate expected qmax for sample
          params_for_max <- c(
            list(q = rep(self$max_value, sample_size)),
            params_list_sample
          )
          expected_maxq <- do.call(pfun_obj, params_for_max)

          # Calculate expected qmin for sample
          params_for_min <- c(
            list(q = rep(self$min_value, sample_size)),
            params_list_sample
          )
          expected_minq <- do.call(pfun_obj, params_for_min)

          # Compare with existing values
          existing_maxq <- tbl[["maxq"]][sample_idx]
          existing_minq <- tbl[["minq"]][sample_idx]

          # Check if values match within tolerance
          maxq_mismatch <- any(
            abs(existing_maxq - expected_maxq) > tolerance,
            na.rm = TRUE
          )
          minq_mismatch <- any(
            abs(existing_minq - expected_minq) > tolerance,
            na.rm = TRUE
          )

          if (maxq_mismatch || minq_mismatch) {
            message(
              "Existing minq/maxq values do not match expected values for ",
              self$name,
              ". Recalculating..."
            )
            needs_recalc <- TRUE
          }
        }

        if (needs_recalc) {
          # If min_value and max_value are specified, calculate from distribution
          if (!is.null(self$min_value) && !is.null(self$max_value)) {
            # Build parameter list for pfun
            params_list <- lapply(self$qparams, function(p) tbl[[p]])
            names(params_list) <- self$qparams

            # Get the probability function
            pfun_obj <- get(self$pfun)

            # Calculate qmax using max_value
            params_for_max <- c(
              list(q = rep(self$max_value, nrow(tbl))),
              params_list
            )
            tbl[, maxq := do.call(pfun_obj, params_for_max)]

            # Calculate qmin using min_value
            params_for_min <- c(
              list(q = rep(self$min_value, nrow(tbl))),
              params_list
            )
            tbl[, minq := do.call(pfun_obj, params_for_min)]
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
              arrow::write_parquet(
                tbl,
                file.path(self$file_path, "part-0.parquet")
              )
            } else {
              # Write as partitioned dataset
              arrow::write_dataset(
                tbl,
                path = self$file_path,
                format = "parquet",
                partitioning = valid_partition_cols,
                existing_data_behavior = "overwrite"
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

      # Read column names from parquet dataset efficiently (metadata only)
      # read_col_names ----
      read_col_names = function() {
        # Omit partitioning to let Arrow auto-detect Hive-style partitions
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

        # Build the quantile input: minq + rank * (maxq - minq)
        q_input <- substitute2(
          minq + .rank * (maxq - minq),
          list(.rank = rank_sym)
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
