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

# -----------------------------------------------------------------------------
# Simulation class table export methods
# This file adds table export methods to the Simulation class using $set()
# -----------------------------------------------------------------------------


# safe_fquantile_byid ----
# Safe wrapper for CKutils::fquantile_byid that handles empty vectors.
# The underlying C++ function crashes with a segfault when passed empty vectors.
# This wrapper returns an empty data.table with the expected structure instead.
safe_fquantile_byid <- function(x, q, id, rounding = FALSE) {
  # Handle empty vector case - return empty data.table with expected structure
  if (length(x) == 0L) {
    result <- data.table::data.table(id = character(0))
    for (qi in q) {
      result[, (as.character(qi)) := numeric(0)]
    }
    return(result)
  }

  # Call the original function for non-empty vectors
  fquantile_byid(x, q, id, rounding)
}


# export_tables ----
# Exports summary tables from simulation summaries.
# See main class documentation in Simulation_class.R for details.
#
# @param strata A named list specifying stratification levels for different table types.
#   If NULL (default), uses built-in defaults matching process_out_for_NotinghamLA.R.
#   The list can contain:
#   - ons: List of character vectors for non-standardised main tables
#   - esp: List of character vectors for standardised (ESP) main tables
#   - mrtl_ons: List of character vectors for non-standardised all-cause mortality tables
#   - mrtl_esp: List of character vectors for standardised all-cause mortality tables
#   - disease_char: List of character vectors for disease characteristics tables
#   - xps_ons: List of character vectors for non-standardised exposure tables
#   - xps_esp: List of character vectors for standardised exposure tables
#
#   Valid stratification variables:
#   - year: Simulation year (always required)
#   - sex: Sex (men/women)
#   - agegrp: Age groups (5-year bands)
#   - dimd: IMD deciles (10 levels: "1 most deprived" to "10 least deprived")
#   - agegrp20: 20-year age groups (used in xps tables)
#   - qimd: IMD quintiles (5 levels, used in xps tables)
#
# @examples
# # Use default strata (includes dimd and qimd stratification)
# sim$export_tables()
#
# # Custom strata - minimal outputs
# sim$export_tables(strata = list(
#   ons = list("year", c("year", "sex")),
#   esp = list("year")
# ))
#
# # Full stratification with dimd
# sim$export_tables(strata = list(
#   ons = list("year", c("year", "sex"), c("year", "dimd"),
#              c("year", "agegrp"), c("year", "agegrp", "sex"),
#              c("year", "agegrp", "sex", "dimd")),
#   esp = list("year", c("year", "sex"), c("year", "dimd"),
#              c("year", "sex", "dimd"))
# ))
Simulation$set("public", "export_tables", function(
    baseline_year_for_change_outputs = 2019L,
    prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
    comparator_scenario = "sc0",
    two_agegrps = FALSE,
    strata = NULL,
    multicore = TRUE
) {
  # Ensure baseline year is in full format (e.g. 2019, not 19)
  # Data is converted to full year format in export_main_tables()
  if (baseline_year_for_change_outputs <= 100) {
    baseline_year_for_change_outputs <- baseline_year_for_change_outputs + 2000L
  }

  # Thread control for parallel execution
  if (multicore) {
    arrow::set_cpu_count(1L)
    data.table::setDTthreads(threads = 1L, restore_after_fork = NULL)
    fst::threads_fst(nr_of_threads = 1L, reset_after_fork = NULL)
  } else {
    arrow::set_cpu_count(self$design$sim_prm$clusternumber_export)
    data.table::setDTthreads(
      threads = self$design$sim_prm$clusternumber_export,
      restore_after_fork = NULL
    )
    fst::threads_fst(
      nr_of_threads = self$design$sim_prm$clusternumber_export,
      reset_after_fork = NULL
    )
  }

  # Build strata configuration (merge user-provided with defaults)
  strata_cfg <- private$build_strata_config(strata, two_agegrps)

  tables_subdir <- if (two_agegrps) "tables2agegrps" else "tables"
  tables_dir <- private$output_dir(tables_subdir)
  private$create_new_folder(tables_dir)

  # Build task list for parallel execution
  tasks <- list(
    list(
      id = 1L,
      type = "main",
      prbl = prbl,
      baseline_year = baseline_year_for_change_outputs,
      output_dir = private$output_dir(),
      tables_dir = tables_dir,
      comparator_scenario = comparator_scenario,
      two_agegrps = two_agegrps,
      strata_ons = strata_cfg$ons,
      strata_esp = strata_cfg$esp
    ),
    list(
      id = 2L,
      type = "all_cause_mrtl",
      prbl = prbl,
      summaries_dir = private$output_dir("summaries"),
      tables_dir = tables_dir,
      strata_ons = strata_cfg$mrtl_ons,
      strata_esp = strata_cfg$mrtl_esp
    ),
    list(
      id = 3L,
      type = "disease_char",
      prbl = prbl,
      summaries_dir = private$output_dir("summaries"),
      tables_dir = tables_dir,
      strata = strata_cfg$disease_char
    ),
    list(
      id = 4L,
      type = "xps",
      prbl = prbl,
      output_dir = private$output_dir(),
      tables_dir = tables_dir,
      strata_ons = strata_cfg$xps_ons,
      strata_esp = strata_cfg$xps_esp
    )
  )

  if (multicore) {
    if (self$design$sim_prm$logs) {
      private$time_mark("Start exporting tables (parallel)")
    }

    n_cores <- min(length(tasks), self$design$sim_prm$clusternumber_export)

    if (.Platform$OS.type == "windows") {
      cl <- parallelly::makeClusterPSOCK(
        n_cores,
        dryrun = FALSE,
        quiet = !self$design$sim_prm$logs,
        rscript_startup = quote(local({
          library(CKutils)
          library(IMPACTncdEngland)
          library(R6)
          library(data.table)
          library(scales)
        })),
        rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"),
        setup_strategy = "parallel"
      )
      on.exit(parallel::stopCluster(cl), add = TRUE)

      parallel::parLapplyLB(
        cl = cl,
        X = tasks,
        fun = function(task) {
          private$export_tables_hlpr(task, implicit_parallelism = FALSE)
          NULL
        }
      )
    } else {
      # Linux/macOS: forking
      doParallel::registerDoParallel(n_cores)
      foreach::foreach(
        task = tasks,
        .inorder = FALSE,
        .packages = c("R6", "CKutils", "IMPACTncdEngland", "data.table", "scales"),
        .verbose = self$design$sim_prm$logs
      ) %dopar% {
        private$export_tables_hlpr(task, implicit_parallelism = FALSE)
        NULL
      }
    }

    if (self$design$sim_prm$logs) {
      private$time_mark("End exporting tables (parallel)")
    }
  } else {
    # Sequential execution
    lapply(tasks, function(task) {
      private$export_tables_hlpr(task, implicit_parallelism = TRUE)
    })
  }

  invisible(self)
})


# export_tables_hlpr ----
# Helper function for parallel table export. Dispatches to the appropriate
# export function based on task type.
Simulation$set("private", "export_tables_hlpr", function(task, implicit_parallelism) {
  # Thread control
  if (implicit_parallelism) {
    arrow::set_cpu_count(self$design$sim_prm$clusternumber_export)
    data.table::setDTthreads(
      threads = self$design$sim_prm$clusternumber_export,
      restore_after_fork = NULL
    )
    fst::threads_fst(
      nr_of_threads = self$design$sim_prm$clusternumber_export,
      reset_after_fork = NULL
    )
  } else {
    arrow::set_cpu_count(1L)
    data.table::setDTthreads(threads = 1L, restore_after_fork = NULL)
    fst::threads_fst(nr_of_threads = 1L, reset_after_fork = NULL)
  }

  # Dispatch to appropriate export function
  switch(task$type,
    "main" = private$export_main_tables(
      prbl = task$prbl,
      baseline_year = task$baseline_year,
      output_dir = task$output_dir,
      tables_dir = task$tables_dir,
      comparator_scenario = task$comparator_scenario,
      two_agegrps = task$two_agegrps,
      strata_ons = task$strata_ons,
      strata_esp = task$strata_esp
    ),
    "all_cause_mrtl" = private$export_all_cause_mrtl_tables(
      prbl = task$prbl,
      summaries_dir = task$summaries_dir,
      tables_dir = task$tables_dir,
      strata_ons = task$strata_ons,
      strata_esp = task$strata_esp
    ),
    "disease_char" = private$export_disease_characteristics_tables(
      prbl = task$prbl,
      summaries_dir = task$summaries_dir,
      tables_dir = task$tables_dir,
      strata = task$strata
    ),
    "xps" = private$export_xps_tables(
      prbl = task$prbl,
      output_dir = task$output_dir,
      tables_dir = task$tables_dir,
      strata_ons = task$strata_ons,
      strata_esp = task$strata_esp
    )
  )

  gc(verbose = FALSE)
  invisible(NULL)
})


# build_strata_config ----
# Builds the strata configuration by merging user-provided strata with defaults.
# Defaults match the stratification from process_out_for_NotinghamLA.R
Simulation$set("private", "build_strata_config", function(user_strata, two_agegrps = FALSE) {
  # Default strata configurations (from process_out_for_NotinghamLA.R)
  if (two_agegrps) {
    # For two_agegrps mode, only include agegrp-based strata
    defaults <- list(
      ons = list(
        c("year", "agegrp"),
        c("year", "agegrp", "sex"),
        c("year", "agegrp", "sex", "dimd")
      ),
      esp = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      mrtl_ons = list(
        c("year", "agegrp"),
        c("year", "agegrp", "sex"),
        c("year", "agegrp", "sex", "dimd")
      ),
      mrtl_esp = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      disease_char = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      xps_ons = list(
        "year",
        c("year", "agegrp20"),
        c("year", "sex"),
        c("year", "qimd"),
        c("year", "agegrp20", "sex"),
        c("year", "agegrp20", "sex", "qimd")
      ),
      xps_esp = list(
        "year",
        c("year", "sex"),
        c("year", "qimd"),
        c("year", "sex", "qimd")
      )
    )
  } else {
    # Standard strata configurations
    defaults <- list(
      ons = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "agegrp"),
        c("year", "agegrp", "sex"),
        c("year", "agegrp", "sex", "dimd")
      ),
      esp = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      mrtl_ons = list(
        "year",
        c("year", "sex"),
        c("year", "agegrp"),
        c("year", "agegrp", "sex"),
        c("year", "agegrp", "sex", "dimd")
      ),
      mrtl_esp = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      disease_char = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      xps_ons = list(
        "year",
        c("year", "agegrp20"),
        c("year", "sex"),
        c("year", "qimd"),
        c("year", "agegrp20", "sex"),
        c("year", "agegrp20", "sex", "qimd")
      ),
      xps_esp = list(
        "year",
        c("year", "sex"),
        c("year", "qimd"),
        c("year", "sex", "qimd")
      )
    )
  }

  # If no user strata provided, return defaults
  if (is.null(user_strata)) {
    return(defaults)
  }

  # Merge user-provided strata with defaults (user overrides defaults)
  result <- defaults
  for (name in names(user_strata)) {
    if (name %in% names(defaults)) {
      result[[name]] <- user_strata[[name]]
    }
  }

  return(result)
})


# tbl_smmrs_core ----
# Core table summary logic adapted from process_out_Bradford.R tbl_smmrs()
# Handles aggregation, rate calculation, and quantile computation
Simulation$set("private", "tbl_smmrs_core", function(
    tt,                    # data.table with summary data
    what,                  # metric type
    population,            # "ons" or "esp"
    strata,                # strata list (e.g., list("year", c("year", "sex")))
    prbl,                  # quantile probabilities
    baseline_year,         # for _change calculations
    comparator_scenario,   # for comparison metrics
    comparison_starting_year,
    tables_dir,            # output directory
    two_agegrps = FALSE,
    qaly_discount_rate = 0,    # annual discount rate for QALYs (%)
    cost_discount_rate = 0,    # annual discount rate for costs (%)
    discount_from_year = NULL  # first year from which discounting starts
) {
  # String mappings for file paths and column patterns (from process_out_Bradford.R)
  str0 <- c(
    "prvl" = "prvl", "prvl_change_relative" = "prvl", "prvl_change_absolute" = "prvl",
    "incd" = "incd", "incd_change_relative" = "incd", "incd_change_absolute" = "incd",
    "ftlt" = "dis_mrtl", "ftlt_change_relative" = "dis_mrtl", "ftlt_change_absolute" = "dis_mrtl",
    "mrtl" = "mrtl", "mrtl_change_relative" = "mrtl", "mrtl_change_absolute" = "mrtl",
    "dis_mrtl" = "dis_mrtl", "dis_mrtl_change_relative" = "dis_mrtl", "dis_mrtl_change_absolute" = "dis_mrtl",
    "qalys" = "qalys", "net_qalys" = "qalys",
    "costs" = "costs", "net_costs" = "costs",
    "cypp" = "prvl", "cpp" = "incd", "dpp" = "mrtl",
    "pop" = "prvl"
  )

  # Column patterns for grep
  str2 <- c(
    "prvl" = "_prvl$|^popsize$", "prvl_change_relative" = "_prvl$|^popsize$", "prvl_change_absolute" = "_prvl$|^popsize$",
    "incd" = "_incd$|^popsize$", "incd_change_relative" = "_incd$|^popsize$", "incd_change_absolute" = "_incd$|^popsize$",
    "ftlt" = "_deaths$|_prvl$", "ftlt_change_relative" = "_deaths$|_prvl$", "ftlt_change_absolute" = "_deaths$|_prvl$",
    "mrtl" = "_mrtl$|^popsize$", "mrtl_change_relative" = "_mrtl$|^popsize$", "mrtl_change_absolute" = "_mrtl$|^popsize$",
    "dis_mrtl" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
    "dis_mrtl_change_relative" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
    "dis_mrtl_change_absolute" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
    "qalys" = "^EQ5D5L$", "net_qalys" = "^EQ5D5L$",
    "costs" = "_cost$|^economic_output$", "net_costs" = "_cost$|^economic_output$",
    "cypp" = "_prvl$", "cpp" = "_incd$", "dpp" = "_mrtl$",
    "pop" = "^popsize$"
  )

  # Output column name prefixes
  str3 <- c(
    "prvl" = "prvl_rate_", "prvl_change_relative" = "prct_change_relative_", "prvl_change_absolute" = "abs_change_",
    "incd" = "incd_rate_", "incd_change_relative" = "prct_change_relative_", "incd_change_absolute" = "abs_change_",
    "ftlt" = "ftlt_rate_", "ftlt_change_relative" = "ftlt_rate_", "ftlt_change_absolute" = "ftlt_abs_change_",
    "mrtl" = "mrtl_rate_", "mrtl_change_relative" = "mrtl_change_relative_", "mrtl_change_absolute" = "mrtl_abs_change_",
    "dis_mrtl" = "disease_mrtl_rate_", "dis_mrtl_change_relative" = "disease_mrtl_change_relative_", "dis_mrtl_change_absolute" = "disease_mrtl_abs_change_",
    "qalys" = "qalys_", "net_qalys" = "net_qalys_",
    "costs" = "costs_", "net_costs" = "net_costs_",
    "cypp" = "cypp_", "cpp" = "cpp_", "dpp" = "dpp_",
    "pop" = "pop_size_"
  )

  # Output file descriptions
  str4 <- c(
    "prvl" = "prevalence by ", "prvl_change_relative" = "prevalence relative change by ", "prvl_change_absolute" = "prevalence absolute change by ",
    "incd" = "incidence by ", "incd_change_relative" = "incidence relative change by ", "incd_change_absolute" = "incidence absolute change by ",
    "ftlt" = "case fatality by ", "ftlt_change_relative" = "case fatality relative change by ", "ftlt_change_absolute" = "case fatality absolute change by ",
    "mrtl" = "all-cause mortality by ", "mrtl_change_relative" = "all-cause mortality relative change by ", "mrtl_change_absolute" = "all-cause mortality absolute change by ",
    "dis_mrtl" = "disease-specific mortality by ",
    "dis_mrtl_change_relative" = "disease-specific mortality relative change by ",
    "dis_mrtl_change_absolute" = "disease-specific mortality absolute change by ",
    "qalys" = "QALYs by ", "net_qalys" = "net QALYs by ",
    "costs" = "costs by ", "net_costs" = "net costs by ",
    "cypp" = "case-years prevented or postponed by ",
    "cpp" = "cases prevented or postponed by ",
    "dpp" = "deaths prevented or postponed by ",
    "pop" = "pop size by "
  )

  # Add mc and scenario to strata
  strata <- lapply(strata, function(x) c("mc", "scenario", x))

  # Process each strata combination
  lapply(strata, function(x) {
    if (grepl("^qalys$", what)) {
      # QALYs processing (EQ5D5L only - HUI3 not implemented)
      d <- tt[, .("EQ5D5L" = sum(EQ5D5L)), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "scale", value.name = "QALYs")
      setkeyv(d, c(x[x != "year"], "scale", "year"))
      # Apply discounting: PV = FV / (1 + r)^(year - discount_from_year)
      if (qaly_discount_rate > 0 && !is.null(discount_from_year)) {
        d[, QALYs := QALYs / (1 + qaly_discount_rate / 100) ^ pmax(0, year - discount_from_year)]
      }
      d[, cumulative := cumsum(QALYs), keyby = c(setdiff(x, "year"), "scale")]
      d <- melt(d, id.vars = c(x, "scale"), variable.name = "type")
      d[, type := fifelse(type == "cumulative", "QALYs_cuml", "QALYs")]
      setkey(d, "type", "scale")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
             keyby = eval(setdiff(c(x, "scale"), "mc"))]
      setnames(d, c(setdiff(c(x, "scale"), "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(c(x, "scale"), "mc")))
      setcolorder(d, setdiff(c(x, "scale"), "mc"))

    } else if (grepl("^net_qalys$", what)) {
      # Net QALYs (intervention - baseline, EQ5D5L only - HUI3 not implemented)
      d <- tt[, .("EQ5D5L" = sum(EQ5D5L)), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "scale", value.name = "QALYs")
      # Apply discounting before calculating net QALYs
      if (qaly_discount_rate > 0 && !is.null(discount_from_year)) {
        d[, QALYs := QALYs / (1 + qaly_discount_rate / 100) ^ pmax(0, year - discount_from_year)]
      }
      d_sc0 <- d[scenario == comparator_scenario & year >= comparison_starting_year][, scenario := NULL]
      d <- d[scenario != comparator_scenario & year >= comparison_starting_year][
        d_sc0, on = c(setdiff(x, "scenario"), "scale"), net_QALYs := QALYs - i.QALYs]
      d[, QALYs := NULL]
      setkeyv(d, c(x[x != "year"], "scale", "year"))
      d[, cumulative := cumsum(net_QALYs), keyby = c(setdiff(x, "year"), "scale")]
      d <- melt(d, id.vars = c(x, "scale"), variable.name = "type")
      d[type == "cumulative", type := "net_QALYs_cuml"]
      setkey(d, "type", "scale")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
             keyby = eval(setdiff(c(x, "scale"), "mc"))]
      x <- c(x, "scale")
      setnames(d, c(setdiff(x, "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(x, "mc")))
      setcolorder(d, setdiff(x, "mc"))

    } else if (grepl("^costs$", what)) {
      # Costs processing
      d <- tt[, lapply(.SD, sum), .SDcols = patterns("_cost$|^economic_output$"), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "costs_type", value.name = "costs")
      # Apply discounting: PV = FV / (1 + r)^(year - discount_from_year)
      if (cost_discount_rate > 0 && !is.null(discount_from_year)) {
        d[, costs := costs / (1 + cost_discount_rate / 100) ^ pmax(0, year - discount_from_year)]
      }
      d[, cumulative := cumsum(costs), keyby = c(setdiff(x, "year"), "costs_type")]
      d <- melt(d, id.vars = c(x, "costs_type"), variable.name = "type")
      d[type == "cumulative", type := "costs_cuml"]
      setkey(d, "type", "costs_type")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = TRUE),
             keyby = eval(setdiff(c(x, "costs_type"), "mc"))]
      setnames(d, c(setdiff(c(x, "costs_type"), "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(c(x, "costs_type"), "mc")))
      setcolorder(d, setdiff(c(x, "costs_type"), "mc"))

    } else if (grepl("^net_costs$", what)) {
      # Net costs (intervention - baseline)
      d <- tt[, lapply(.SD, sum), .SDcols = patterns("_cost$|^economic_output$"), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "costs_type", value.name = "value")
      # Apply discounting before calculating net costs
      if (cost_discount_rate > 0 && !is.null(discount_from_year)) {
        d[, value := value / (1 + cost_discount_rate / 100) ^ pmax(0, year - discount_from_year)]
      }
      d_sc0 <- d[scenario == comparator_scenario & year >= comparison_starting_year][, scenario := NULL]
      d <- d[scenario != comparator_scenario & year >= comparison_starting_year][
        d_sc0, on = c(setdiff(x, "scenario"), "costs_type"), net_costs := value - i.value]
      d[, value := NULL]
      setkeyv(d, c(x[x != "year"], "costs_type", "year"))
      d[, cumulative := cumsum(net_costs), keyby = c(setdiff(x, "year"), "costs_type")]
      d <- melt(d, id.vars = c(x, "costs_type"), variable.name = "type")
      d[type == "cumulative", type := "net_costs_cuml"]
      setkey(d, "type", "costs_type")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = TRUE),
             keyby = eval(setdiff(c(x, "costs_type"), "mc"))]
      x <- c(x, "costs_type")
      setnames(d, c(setdiff(x, "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(x, "mc")))
      setcolorder(d, setdiff(x, "mc"))

    } else {
      # All other metrics (prvl, incd, mrtl, dis_mrtl, ftlt, pop, cypp, cpp, dpp)
      d <- tt[, lapply(.SD, sum), .SDcols = patterns(str2[[what]]), keyby = x]

      # Convert integer columns to numeric
      is_int <- sapply(d[, .SD, .SDcols = -x], is.integer)
      is_int <- names(is_int[is_int])
      if (length(is_int) > 0) {
        d[, (is_int) := lapply(.SD, as.numeric), .SDcols = is_int]
      }

      if (grepl("^ftlt", what)) {
        # Case fatality: deaths / prevalence
        nm <- names(d)
        nm <- grep("_deaths$", nm, value = TRUE)
        nm <- gsub("_deaths$", "", nm)
        nm <- setdiff(nm, "alive")
        for (i in nm) {
          set(d, NULL, paste0(i, "_ftlt"),
              d[[paste0(i, "_deaths")]] / d[[paste0(i, "_prvl")]])
        }
        nm <- names(d)
        nm <- grep("_deaths$|_prvl$", nm, value = TRUE)
        d[, (nm) := NULL]
        setnafill(d, "const", 0, cols = grep("_ftlt$", names(d), value = TRUE))
      } else if (!what %in% c("pop", "cypp", "cpp", "dpp")) {
        # Calculate rates (divide by popsize)
        d <- d[, lapply(.SD, function(y) y / popsize), keyby = x]
      }

      d <- melt(d, id.vars = x)

      if (grepl("_change_relative$", what)) {
        # Relative change from baseline year
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year"), "variable"), value := value / i.value]
      }

      if (grepl("_change_absolute$", what)) {
        # Absolute change from baseline year
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year"), "variable"), value := value - i.value]
      }

      if (grepl("^cypp$|^cpp$|^dpp$", what)) {
        # Comparison metrics: baseline - intervention
        d_sc0 <- d[scenario == comparator_scenario & year >= comparison_starting_year][, scenario := NULL]
        d <- d[scenario != comparator_scenario & year >= comparison_starting_year][
          d_sc0, on = c(setdiff(x, "scenario"), "variable"), value := i.value - value]
        d[, variable := gsub(paste0("_", str0[[what]]), "", variable)]
        setkeyv(d, c(x[x != "year"], "variable", "year"))
        d[, cumulative := cumsum(value), keyby = c(setdiff(x, "year"), "variable")]
        d <- melt(d, id.vars = c(x, "variable"), variable.name = "type")
        d[, type := fifelse(type == "cumulative", paste0(what, "_cuml"), what)]
        setkey(d, "type", "variable")
        d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable),
                                rounding = (what %in% c("pop", "cypp", "cpp", "dpp"))),
               keyby = eval(setdiff(c(x, "type"), "mc"))]
        x <- c(x, "type")
        setnames(d, c(setdiff(x, "mc"), "disease", scales::percent(prbl, prefix = str3[[what]])))
      } else {
        setkey(d, "variable")
        d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable),
                                rounding = what == "pop"),
               keyby = eval(setdiff(x, "mc"))]
        setnames(d, c(setdiff(x, "mc"), "disease", scales::percent(prbl, prefix = str3[[what]])))
      }

      if (what == "pop") {
        d[, disease := NULL]
      } else {
        if ("popsize" %in% d$disease) d <- d[disease != "popsize"]
      }
      setkeyv(d, setdiff(x, "mc"))
      setcolorder(d, setdiff(x, "mc"))
    }

    # Build output filename
    str5 <- c(
      "ons" = " (not standardised).csv",
      "esp" = paste0(" (", paste(setdiff(c("mc", "scenario", "year", "age", "sex"), x),
                                 collapse = "-"), " standardised).csv")
    )
    str6 <- paste0(
      str4[[what]],
      paste(setdiff(x, c("mc", "scenario", "type", "scale", "costs_type")), collapse = "-"),
      str5[[population]]
    )

    fwrite(d, file.path(tables_dir, str6))
  })

  invisible(NULL)
})


# export_main_tables ----
# Generate main summary tables (prevalence, incidence, mortality, etc.)
# Memory-optimized: processes by source dataset to minimize simultaneous memory
Simulation$set("private", "export_main_tables", function(
    prbl,
    baseline_year,
    output_dir,
    tables_dir,
    comparator_scenario = "sc0",
    two_agegrps = FALSE,
    strata_ons = NULL,
    strata_esp = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating main summary tables...")
  }

  # String mappings for source datasets
  str0 <- c(
    "prvl" = "prvl", "prvl_change_relative" = "prvl", "prvl_change_absolute" = "prvl",
    "incd" = "incd", "incd_change_relative" = "incd", "incd_change_absolute" = "incd",
    "ftlt" = "dis_mrtl", "ftlt_change_relative" = "dis_mrtl", "ftlt_change_absolute" = "dis_mrtl",
    "mrtl" = "mrtl", "mrtl_change_relative" = "mrtl", "mrtl_change_absolute" = "mrtl",
    "dis_mrtl" = "dis_mrtl", "dis_mrtl_change_relative" = "dis_mrtl", "dis_mrtl_change_absolute" = "dis_mrtl",
    "qalys" = "qalys", "net_qalys" = "qalys",
    "costs" = "costs", "net_costs" = "costs",
    "cypp" = "prvl", "cpp" = "incd", "dpp" = "mrtl",
    "pop" = "prvl"
  )
  str1 <- c("ons" = "scaled_up", "esp" = "esp")

  # Group metrics by source dataset for efficient memory usage
  source_to_metrics <- list(
    prvl = c("prvl", "prvl_change_relative", "prvl_change_absolute", "cypp", "pop"),
    incd = c("incd", "incd_change_relative", "incd_change_absolute", "cpp"),
    dis_mrtl = c("ftlt", "ftlt_change_relative", "ftlt_change_absolute", "dis_mrtl", "dis_mrtl_change_relative", "dis_mrtl_change_absolute"),
    mrtl = c("mrtl", "mrtl_change_relative", "mrtl_change_absolute", "dpp"),
    qalys = c("qalys", "net_qalys"),
    costs = c("costs", "net_costs")
  )

  # Process each source dataset group
  for (source_name in names(source_to_metrics)) {
    metrics_for_source <- source_to_metrics[[source_name]]

    # Process both populations for this source
    for (pop_name in c("ons", "esp")) {
      pop_key <- str1[[pop_name]]

      # Load dataset once for this source/population combination
      tt_base <- private$read_summary_dataset(source_name, pop_key)
      if (is.null(tt_base)) next
      # Convert year from short format (19) to full format (2019)
      tt_base[, year := year + 2000L]

      # Also load prvl for ftlt metrics (case fatality denominator)
      prvl_for_ftlt <- NULL
      if (source_name == "dis_mrtl") {
        prvl_for_ftlt <- private$read_summary_dataset("prvl", pop_key)
        if (!is.null(prvl_for_ftlt)) {
          prvl_for_ftlt[, year := year + 2000L]
        }
      }

      # Process each metric that uses this source
      for (what in metrics_for_source) {
        # Skip pop for esp
        if (what == "pop" && pop_name == "esp") next
        if (grepl("_age", what) && pop_name == "esp") next

        # Use configurable strata (already filtered by two_agegrps in build_strata_config)
        if (pop_name == "ons") {
          strata <- strata_ons
        } else {
          strata <- strata_esp
        }

        if (self$design$sim_prm$logs) {
          message(paste0("  ", what, "-", pop_name))
        }

        # Get a copy of the base dataset
        tt <- copy(tt_base)

        # Check if comparison metrics can be computed
        comparison_metrics <- c("cypp", "cpp", "dpp", "net_qalys", "net_costs")
        if (what %in% comparison_metrics) {
          available_scenarios <- unique(tt$scenario)
          non_comparator_scenarios <- setdiff(available_scenarios, comparator_scenario)
          if (length(non_comparator_scenarios) == 0) {
            if (self$design$sim_prm$logs) {
              message("    Skipping ", what, " - no intervention scenarios (only '",
                      comparator_scenario, "' found)")
            }
            rm(tt)
            next
          }
        }

        # Handle two_agegrps transformation
        if (two_agegrps && "agegrp" %in% names(tt)) {
          tt[agegrp %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64"),
             agegrp := "30-64"]
          tt[agegrp %in% c("65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99"),
             agegrp := "65-99"]
        }

        # For case fatality, add prevalence denominator
        if (grepl("^ftlt", what) && !is.null(prvl_for_ftlt)) {
          t1 <- copy(prvl_for_ftlt)
          setnames(t1, "popsize", "nonmodelled_prvl")
          if (two_agegrps && "agegrp" %in% names(t1)) {
            t1[agegrp %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64"),
               agegrp := "30-64"]
            t1[agegrp %in% c("65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99"),
               agegrp := "65-99"]
          }
          absorb_dt(tt, t1)
          tt <- tt[nonmodelled_prvl > 0]
          rm(t1)
        }

        # Get discount settings from design (with defaults for backwards compatibility)
        disc <- self$design$sim_prm$discounting
        qaly_rate <- if (!is.null(disc$qaly_discount_rate)) disc$qaly_discount_rate else 0
        cost_rate <- if (!is.null(disc$cost_discount_rate)) disc$cost_discount_rate else 0
        disc_year <- disc$discount_from_year

        # Generate tables
        private$tbl_smmrs_core(
          tt = tt,
          what = what,
          population = pop_name,
          strata = strata,
          prbl = prbl,
          baseline_year = baseline_year,
          comparator_scenario = comparator_scenario,
          comparison_starting_year = baseline_year,
          tables_dir = tables_dir,
          two_agegrps = two_agegrps,
          qaly_discount_rate = qaly_rate,
          cost_discount_rate = cost_rate,
          discount_from_year = disc_year
        )
        rm(tt)
      }

      # Cleanup after processing this source/population
      rm(tt_base)
      if (!is.null(prvl_for_ftlt)) rm(prvl_for_ftlt)
    }

    # Garbage collect after each source group
    gc(verbose = FALSE)
  }

  invisible(NULL)
})


# export_all_cause_mrtl_tables ----
# Generate all-cause mortality by disease tables
# Memory-optimized: loads datasets once and reuses them
Simulation$set("private", "export_all_cause_mrtl_tables", function(
    prbl,
    summaries_dir,
    tables_dir,
    strata_ons = NULL,
    strata_esp = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating all-cause mortality by disease tables...")
  }

  # Helper function to convert user strata to internal format with mc and scenario
  make_strata_configs <- function(strata_list, standardised = FALSE) {
    lapply(strata_list, function(s) {
      outstrata <- c("mc", s, "scenario")
      suffix <- paste(s, collapse = "-")
      # Map agegrp to agegroup in suffix for backwards compatibility
      suffix <- gsub("agegrp", "agegroup", suffix)
      if (standardised) {
        # Determine what was standardised by (variables NOT in strata)
        possible_vars <- c("age", "sex", "dimd")
        standardised_vars <- setdiff(possible_vars, s)
        std_suffix <- paste(standardised_vars, collapse = "-")
        list(strata = outstrata, suffix = suffix, std = std_suffix)
      } else {
        list(strata = outstrata, suffix = suffix)
      }
    })
  }

  # Load datasets once (avoid duplicate reads)
  tt_scaled <- private$read_summary_dataset("all_cause_mrtl_by_dis", "scaled_up")
  pp_scaled <- private$read_summary_dataset("prvl", "scaled_up")
  tt_esp <- private$read_summary_dataset("all_cause_mrtl_by_dis", "esp")

  # Convert year from short format (19) to full format (2019)
  if (!is.null(tt_scaled)) tt_scaled[, year := year + 2000L]
  if (!is.null(pp_scaled)) pp_scaled[, year := year + 2000L]
  if (!is.null(tt_esp)) tt_esp[, year := year + 2000L]

  # ---- Non-standardised with disease denominator ----
  if (!is.null(tt_scaled)) {
    strata_configs <- make_strata_configs(strata_ons, standardised = FALSE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- tt_scaled[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"), keyby = eval(outstrata)]
      d <- melt(d, id.vars = outstrata)
      cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
      d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
      d[cases, on = c(outstrata, "variable"), value := value / i.value]
      rm(cases)  # Cleanup intermediate
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "disease",
                    scales::percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("all-cause mortality given disease-", cfg$suffix, " (not standardised).csv")))
      rm(d)
    }
  }

  # ---- Non-standardised with population denominator ----
  if (!is.null(tt_scaled) && !is.null(pp_scaled)) {
    strata_configs <- make_strata_configs(strata_ons, standardised = FALSE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      cases <- pp_scaled[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
      d <- tt_scaled[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"), keyby = eval(outstrata)]
      d <- melt(d, id.vars = outstrata)
      d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
      d[cases, on = outstrata, value := value / popsize]
      rm(cases)  # Cleanup intermediate
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "disease",
                    scales::percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("all-cause mortality given disease-", cfg$suffix, " popdenom (not standardised).csv")))
      rm(d)
    }
  }

  # Cleanup scaled_up datasets before ESP processing
  rm(tt_scaled, pp_scaled)
  gc(verbose = FALSE)

  # ---- Standardised (ESP) ----
  if (!is.null(tt_esp)) {
    strata_configs <- make_strata_configs(strata_esp, standardised = TRUE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- tt_esp[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"), keyby = eval(outstrata)]
      d <- melt(d, id.vars = outstrata)
      cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
      d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
      d[cases, on = c(outstrata, "variable"), value := value / i.value]
      rm(cases)  # Cleanup intermediate
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "disease",
                    scales::percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("all-cause mortality given disease-", cfg$suffix,
                                 " (", cfg$std, " standardised).csv")))
      rm(d)
    }
  }

  rm(tt_esp)
  invisible(NULL)
})


# export_disease_characteristics_tables ----
# Generate disease characteristics tables (duration, age metrics, CMS)
Simulation$set("private", "export_disease_characteristics_tables", function(
    prbl,
    summaries_dir,
    tables_dir,
    strata = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating disease characteristics tables...")
  }

  tt <- private$read_summary_dataset("dis_characteristics", "scaled_up")
  if (is.null(tt)) return(invisible(NULL))

  # Convert year from short format (19) to full format (2019)
  tt[, year := year + 2000L]

  # Type conversion if needed
  if ("mean_cms_count_cms1st_cont" %in% names(tt)) {
    tt[, mean_cms_count_cms1st_cont := as.numeric(mean_cms_count_cms1st_cont)]
  }

  # Derive id variables dynamically from requested strata, keeping only
  # columns that actually exist in the data (e.g. agegrp is excluded from

  # disease characteristics summaries via strata_noagegrp)
  all_strata_vars <- unique(unlist(strata))
  id_vars <- intersect(unique(c("mc", "scenario", all_strata_vars)), names(tt))
  id_pattern <- paste(id_vars, collapse = "|")

  # Extract case counts for weighting
  d1 <- tt[, .SD, .SDcols = patterns(paste0(id_pattern, "|^cases_"))]
  d1 <- melt(d1, id.vars = id_vars)
  d1 <- unique(d1, by = c(id_vars, "variable"))
  d1[, disease := gsub("^cases_", "", variable)]
  d1[, variable := NULL]

  # Extract characteristics columns
  char_patterns <- paste0(id_pattern, "|^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_")
  tt <- tt[, .SD, .SDcols = patterns(char_patterns)]

  if ("mean_cms_count_cmsmm1" %in% names(tt)) {
    tt[, mean_cms_count_cmsmm1 := as.double(mean_cms_count_cmsmm1)]
  }

  tt <- melt(tt, id.vars = id_vars)
  tt[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
  tt[d1, on = c(id_vars, "disease"), cases := i.value]

  # Helper function to convert user strata to internal format
  make_strata_configs <- function(strata_list) {
    lapply(strata_list, function(s) {
      outstrata <- c("mc", s, "scenario")
      suffix <- paste(s, collapse = "-")
      list(strata = outstrata, suffix = suffix)
    })
  }

  # Process each strata
  strata_configs <- make_strata_configs(strata)

  for (cfg in strata_configs) {
    outstrata <- cfg$strata
    d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")]
    setkey(d, "variable")
    d <- d[, safe_fquantile_byid(V1, prbl, id = as.character(variable)),
           keyby = eval(setdiff(outstrata, "mc"))]
    setnames(d, c(setdiff(outstrata, "mc"), "variable", scales::percent(prbl, prefix = "value_")))

    # Parse variable name to extract disease and type
    d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
    d[grep("^mean_duration_", variable), type := "mean_duration"]
    d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
    d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
    d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
    d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
    d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
    d[, variable := NULL]
    setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
    setcolorder(d)

    fwrite(d, file.path(tables_dir,
                        paste0("disease characteristics by ", cfg$suffix, " (not standardised).csv")))
    rm(d)  # Cleanup after each output
  }

  # Cleanup
  rm(tt, d1)
  invisible(NULL)
})


# export_xps_tables ----
# Generate exposure summary tables
# Memory-optimized: adds cleanup between sections
Simulation$set("private", "export_xps_tables", function(
    prbl,
    output_dir,
    tables_dir,
    strata_ons = NULL,
    strata_esp = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating exposure tables...")
  }

  # Helper function to build filter expression from strata
  # Variables in strata -> filter on != "All"
  # Variables not in strata -> filter on == "All"
  make_xps_strata_configs <- function(strata_list, filterable_vars, standardised = FALSE) {
    lapply(strata_list, function(s) {
      outstrata <- c("mc", s, "scenario")
      suffix <- paste(setdiff(s, "year"), collapse = "-")
      if (suffix == "") suffix <- "year" else suffix <- paste0("year-", suffix)
      # Map agegrp20 to agegroup in suffix for backwards compatibility
      suffix <- gsub("agegrp20", "agegroup", suffix)

      # Build filter expression: vars in strata -> != "All", vars not in strata -> == "All"
      filter_parts <- character(0)
      for (v in filterable_vars) {
        if (v %in% s) {
          filter_parts <- c(filter_parts, paste0(v, " != 'All'"))
        } else {
          filter_parts <- c(filter_parts, paste0(v, " == 'All'"))
        }
      }
      filter_str <- paste(filter_parts, collapse = " & ")
      filter_expr <- if (length(filter_parts) > 0) parse(text = filter_str)[[1]] else quote(TRUE)

      if (standardised) {
        # Determine what was standardised by: all filterable vars not in the strata
        # Convert agegrp20 to age for standardisation naming
        s_for_std <- gsub("agegrp20", "age", s)
        all_std_vars <- gsub("agegrp20", "age", filterable_vars)
        standardised_vars <- setdiff(all_std_vars, s_for_std)
        std_suffix <- paste(standardised_vars, collapse = "-")
        list(strata = outstrata, suffix = suffix, filter_expr = filter_expr, std = std_suffix)
      } else {
        list(strata = outstrata, suffix = suffix, filter_expr = filter_expr)
      }
    })
  }

  # Helper to detect filterable vars: columns with "All" values from groupingsets
  detect_filterable_vars <- function(dt) {
    candidates <- setdiff(names(dt), c("mc", "scenario", "year"))
    candidates[vapply(candidates, function(v) {
      is.character(dt[[v]]) && "All" %in% dt[[v]]
    }, logical(1))]
  }

  # ---- Non-standardised (xps20) ----
  xps_path <- file.path(output_dir, "xps", "xps20")
  if (dir.exists(xps_path)) {
    xps_tab <- CKutils::read_parquet_dt(xps_path)
    # Convert year from short format (19) to full format (2019)
    xps_tab[, year := year + 2000L]

    # Detect which columns have "All" marginals from groupingsets
    filt_vars <- detect_filterable_vars(xps_tab)
    strata_configs <- make_xps_strata_configs(strata_ons, filt_vars, standardised = FALSE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- xps_tab[eval(cfg$filter_expr)]
      if (nrow(d) == 0) next
      d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
      d <- melt(d, id.vars = outstrata)
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "exposure", scales::percent(prbl, prefix = "xps_mean_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("exposures by ", cfg$suffix, " (not standardised).csv")))
      rm(d)  # Cleanup after each output
    }
    rm(xps_tab)  # Cleanup before next section
    gc(verbose = FALSE)
  }

  # ---- Standardised (xps5 / ESP) ----
  xps_path <- file.path(output_dir, "xps", "xps5")
  if (dir.exists(xps_path)) {
    xps_tab <- CKutils::read_parquet_dt(xps_path)
    # Convert year from short format (19) to full format (2019)
    xps_tab[, year := year + 2000L]

    # Detect which columns have "All" marginals from groupingsets
    filt_vars <- detect_filterable_vars(xps_tab)
    strata_configs <- make_xps_strata_configs(strata_esp, filt_vars, standardised = TRUE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- xps_tab[eval(cfg$filter_expr)]
      if (nrow(d) == 0) next
      d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
      d <- melt(d, id.vars = outstrata)
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "exposure", scales::percent(prbl, prefix = "xps_mean_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("exposures by ", cfg$suffix, " (", cfg$std, " standardised).csv")))
      rm(d)  # Cleanup after each output
    }
    rm(xps_tab)  # Final cleanup
  }

  invisible(NULL)
})
