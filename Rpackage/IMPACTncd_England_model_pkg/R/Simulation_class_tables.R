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
Simulation$set("public", "export_tables", function(
    baseline_year_for_change_outputs = 2019L,
    prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
    comparator_scenario = "sc0",
    two_agegrps = FALSE
) {
  tables_subdir <- if (two_agegrps) "tables2agegrps" else "tables"
  tables_dir <- private$output_dir(tables_subdir)
  private$create_new_folder(tables_dir)

  private$export_main_tables(
    prbl = prbl,
    baseline_year = baseline_year_for_change_outputs,
    output_dir = private$output_dir(),
    tables_dir = tables_dir,
    comparator_scenario = comparator_scenario,
    two_agegrps = two_agegrps
  )
  gc(verbose = FALSE)  # Cleanup after main tables

  private$export_all_cause_mrtl_tables(
    prbl = prbl,
    summaries_dir = private$output_dir("summaries"),
    tables_dir = tables_dir
  )
  gc(verbose = FALSE)  # Cleanup after mortality tables

  private$export_disease_characteristics_tables(
    prbl = prbl,
    summaries_dir = private$output_dir("summaries"),
    tables_dir = tables_dir
  )
  gc(verbose = FALSE)  # Cleanup after disease tables

  private$export_xps_tables(
    prbl = prbl,
    output_dir = private$output_dir(),
    tables_dir = tables_dir
  )
  gc(verbose = FALSE)  # Final cleanup

  invisible(self)
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
    two_agegrps = FALSE
) {
  # String mappings for file paths and column patterns (from process_out_Bradford.R)
  str0 <- c(
    "prvl" = "prvl", "prvl_change" = "prvl",
    "incd" = "incd", "incd_change" = "incd",
    "ftlt" = "dis_mrtl", "ftlt_change" = "dis_mrtl",
    "mrtl" = "mrtl", "mrtl_change" = "mrtl",
    "dis_mrtl" = "dis_mrtl", "dis_mrtl_change" = "dis_mrtl",
    "qalys" = "qalys", "net_qalys" = "qalys",
    "costs" = "costs", "net_costs" = "costs",
    "cypp" = "prvl", "cpp" = "incd", "dpp" = "mrtl",
    "pop" = "prvl"
  )

  # Column patterns for grep
  str2 <- c(
    "prvl" = "_prvl$|^popsize$", "prvl_change" = "_prvl$|^popsize$",
    "incd" = "_incd$|^popsize$", "incd_change" = "_incd$|^popsize$",
    "ftlt" = "_deaths$|_prvl$", "ftlt_change" = "_deaths$|_prvl$",
    "mrtl" = "_mrtl$|^popsize$", "mrtl_change" = "_mrtl$|^popsize$",
    "dis_mrtl" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
    "dis_mrtl_change" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
    "qalys" = "^EQ5D5L$|^HUI3$", "net_qalys" = "^EQ5D5L$|^HUI3$",
    "costs" = "_costs$", "net_costs" = "_costs$",
    "cypp" = "_prvl$", "cpp" = "_incd$", "dpp" = "_mrtl$",
    "pop" = "^popsize$"
  )

  # Output column name prefixes
  str3 <- c(
    "prvl" = "prvl_rate_", "prvl_change" = "prct_change_",
    "incd" = "incd_rate_", "incd_change" = "prct_change_",
    "ftlt" = "ftlt_rate_", "ftlt_change" = "ftlt_rate_",
    "mrtl" = "mrtl_rate_", "mrtl_change" = "mrtl_change_",
    "dis_mrtl" = "disease_mrtl_rate_", "dis_mrtl_change" = "disease_mrtl_change_",
    "qalys" = "qalys_", "net_qalys" = "net_qalys_",
    "costs" = "costs_", "net_costs" = "net_costs_",
    "cypp" = "cypp_", "cpp" = "cpp_", "dpp" = "dpp_",
    "pop" = "pop_size_"
  )

  # Output file descriptions
  str4 <- c(
    "prvl" = "prevalence by ", "prvl_change" = "prevalence change by ",
    "incd" = "incidence by ", "incd_change" = "incidence change by ",
    "ftlt" = "case fatality by ", "ftlt_change" = "case fatality change by ",
    "mrtl" = "all-cause mortality by ", "mrtl_change" = "all-cause mortality change by ",
    "dis_mrtl" = "disease-specific mortality by ",
    "dis_mrtl_change" = "disease-specific mortality change by ",
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
      # QALYs processing
      d <- tt[, .("EQ5D5L" = sum(EQ5D5L), "HUI3" = sum(HUI3)), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "scale", value.name = "QALYs")
      setkeyv(d, c(x[x != "year"], "scale", "year"))
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
      # Net QALYs (intervention - baseline)
      d <- tt[, .("EQ5D5L" = sum(EQ5D5L), "HUI3" = sum(HUI3)), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "scale", value.name = "QALYs")
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
      d <- tt[, lapply(.SD, sum), .SDcols = patterns("_costs$"), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "costs_type", value.name = "costs")
      d[, cumulative := cumsum(costs), keyby = c(setdiff(x, "year"), "costs_type")]
      d <- melt(d, id.vars = c(x, "costs_type"), variable.name = "type")
      d[type == "cumulative", type := "costs_cuml"]
      setkey(d, "type", "costs_type")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
             keyby = eval(setdiff(c(x, "costs_type"), "mc"))]
      setnames(d, c(setdiff(c(x, "costs_type"), "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(c(x, "costs_type"), "mc")))
      setcolorder(d, setdiff(c(x, "costs_type"), "mc"))

    } else if (grepl("^net_costs$", what)) {
      # Net costs (intervention - baseline)
      d <- tt[, lapply(.SD, sum), .SDcols = patterns("_costs$"), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "costs_type", value.name = "value")
      d_sc0 <- d[scenario == comparator_scenario & year >= comparison_starting_year][, scenario := NULL]
      d <- d[scenario != comparator_scenario & year >= comparison_starting_year][
        d_sc0, on = c(setdiff(x, "scenario"), "costs_type"), net_costs := value - i.value]
      d[, value := NULL]
      setkeyv(d, c(x[x != "year"], "costs_type", "year"))
      d[, cumulative := cumsum(net_costs), keyby = c(setdiff(x, "year"), "costs_type")]
      d <- melt(d, id.vars = c(x, "costs_type"), variable.name = "type")
      d[type == "cumulative", type := "net_costs_cuml"]
      setkey(d, "type", "costs_type")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
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

      if (grepl("_change$", what)) {
        # Relative change from baseline year
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year"), "variable"), value := value / i.value]
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
    two_agegrps = FALSE
) {
  if (self$design$sim_prm$logs) {
    message("Generating main summary tables...")
  }

  # String mappings for source datasets
  str0 <- c(
    "prvl" = "prvl", "prvl_change" = "prvl",
    "incd" = "incd", "incd_change" = "incd",
    "ftlt" = "dis_mrtl", "ftlt_change" = "dis_mrtl",
    "mrtl" = "mrtl", "mrtl_change" = "mrtl",
    "dis_mrtl" = "dis_mrtl", "dis_mrtl_change" = "dis_mrtl",
    "qalys" = "qalys", "net_qalys" = "qalys",
    "costs" = "costs", "net_costs" = "costs",
    "cypp" = "prvl", "cpp" = "incd", "dpp" = "mrtl",
    "pop" = "prvl"
  )
  str1 <- c("ons" = "scaled_up", "esp" = "esp")

  # Group metrics by source dataset for efficient memory usage
  source_to_metrics <- list(
    prvl = c("prvl", "prvl_change", "cypp", "pop"),
    incd = c("incd", "incd_change", "cpp"),
    dis_mrtl = c("ftlt", "ftlt_change", "dis_mrtl", "dis_mrtl_change"),
    mrtl = c("mrtl", "mrtl_change", "dpp"),
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

      # Also load prvl for ftlt metrics (case fatality denominator)
      prvl_for_ftlt <- NULL
      if (source_name == "dis_mrtl") {
        prvl_for_ftlt <- private$read_summary_dataset("prvl", pop_key)
      }

      # Process each metric that uses this source
      for (what in metrics_for_source) {
        # Skip pop for esp
        if (what == "pop" && pop_name == "esp") next
        if (grepl("_age", what) && pop_name == "esp") next

        # Define strata based on population
        if (pop_name == "ons") {
          if (two_agegrps) {
            strata <- list(c("year", "agegrp"), c("year", "agegrp", "sex"))
          } else {
            strata <- list("year", c("year", "sex"), c("year", "agegrp"),
                           c("year", "agegrp", "sex"))
          }
        } else {
          strata <- list("year", c("year", "sex"))
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
          two_agegrps = two_agegrps
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
Simulation$set("private", "export_all_cause_mrtl_tables", function(prbl, summaries_dir, tables_dir) {
  if (self$design$sim_prm$logs) {
    message("Generating all-cause mortality by disease tables...")
  }

  # Load datasets once (avoid duplicate reads)
  tt_scaled <- private$read_summary_dataset("all_cause_mrtl_by_dis", "scaled_up")
  pp_scaled <- private$read_summary_dataset("prvl", "scaled_up")
  tt_esp <- private$read_summary_dataset("all_cause_mrtl_by_dis", "esp")

  # ---- Non-standardised with disease denominator ----
  if (!is.null(tt_scaled)) {
    strata_configs <- list(
      list(strata = c("mc", "year", "scenario"), suffix = "year"),
      list(strata = c("mc", "year", "sex", "scenario"), suffix = "year-sex"),
      list(strata = c("mc", "year", "agegrp", "sex", "scenario"), suffix = "year-agegroup-sex")
    )

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
    strata_configs <- list(
      list(strata = c("mc", "year", "scenario"), suffix = "year"),
      list(strata = c("mc", "year", "sex", "scenario"), suffix = "year-sex"),
      list(strata = c("mc", "year", "agegrp", "sex", "scenario"), suffix = "year-agegroup-sex")
    )

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
    strata_configs <- list(
      list(strata = c("mc", "year", "scenario"), suffix = "year", std = "age-sex"),
      list(strata = c("mc", "year", "sex", "scenario"), suffix = "year-sex", std = "age")
    )

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
Simulation$set("private", "export_disease_characteristics_tables", function(prbl, summaries_dir, tables_dir) {
  if (self$design$sim_prm$logs) {
    message("Generating disease characteristics tables...")
  }

  tt <- private$read_summary_dataset("dis_characteristics", "scaled_up")
  if (is.null(tt)) return(invisible(NULL))

  # Type conversion if needed
  if ("mean_cms_count_cms1st_cont" %in% names(tt)) {
    tt[, mean_cms_count_cms1st_cont := as.numeric(mean_cms_count_cms1st_cont)]
  }

  # Extract case counts for weighting
  d1 <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|^cases_")]
  d1 <- melt(d1, id.vars = c("mc", "year", "scenario", "sex"))
  d1 <- unique(d1, by = c("mc", "year", "scenario", "sex", "variable"))
  d1[, disease := gsub("^cases_", "", variable)]
  d1[, variable := NULL]

  # Extract characteristics columns
  char_patterns <- "mc|scenario|year|sex|^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_"
  tt <- tt[, .SD, .SDcols = patterns(char_patterns)]

  if ("mean_cms_count_cmsmm1" %in% names(tt)) {
    tt[, mean_cms_count_cmsmm1 := as.double(mean_cms_count_cmsmm1)]
  }

  tt <- melt(tt, id.vars = c("mc", "year", "scenario", "sex"))
  tt[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
  tt[d1, on = c("mc", "year", "scenario", "sex", "disease"), cases := i.value]

  # Process each strata
  strata_configs <- list(
    list(strata = c("mc", "year", "scenario"), suffix = "year"),
    list(strata = c("mc", "year", "sex", "scenario"), suffix = "year-sex")
  )

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
Simulation$set("private", "export_xps_tables", function(prbl, output_dir, tables_dir) {
  if (self$design$sim_prm$logs) {
    message("Generating exposure tables...")
  }

  # ---- Non-standardised (xps20) ----
  xps_path <- file.path(output_dir, "xps", "xps20")
  if (dir.exists(xps_path)) {
    xps_tab <- CKutils::read_parquet_dt(xps_path)

    strata_configs <- list(
      list(
        filter_expr = quote(sex != "All" & agegrp20 != "All"),
        strata = c("mc", "year", "agegrp20", "sex", "scenario"),
        suffix = "year-agegroup-sex"
      ),
      list(
        filter_expr = quote(sex == "All" & agegrp20 == "All"),
        strata = c("mc", "year", "scenario"),
        suffix = "year"
      ),
      list(
        filter_expr = quote(sex == "All" & agegrp20 != "All"),
        strata = c("mc", "year", "agegrp20", "scenario"),
        suffix = "year-agegroup"
      ),
      list(
        filter_expr = quote(sex != "All"),
        strata = c("mc", "year", "sex", "scenario"),
        suffix = "year-sex"
      )
    )

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

    strata_configs <- list(
      list(
        filter_expr = quote(sex != "All"),
        strata = c("mc", "year", "sex", "scenario"),
        suffix = "year-sex",
        std = "age"
      ),
      list(
        filter_expr = quote(sex == "All"),
        strata = c("mc", "year", "scenario"),
        suffix = "year",
        std = "age-sex"
      )
    )

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
