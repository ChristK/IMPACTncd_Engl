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


# calc_equity_slope_indices ----
# Absolute and relative equity slope indices (analogues of the Slope Index of
# Inequality, SII, and Relative Index of Inequality, RII) for a benefit
# distributed across ordered socioeconomic groups (DIMD deciles).
#
# Statistical framework (following Moreno-Betancur et al. 2015, Epidemiology;
# Wagstaff et al. 1991; and, for cumulative benefits, Cookson et al. 2026,
# Pharmacoeconomics): each socioeconomic group is assigned a population-weighted
# midpoint cumulative rank (a *ridit* score) in [0, 1], ordered from LEAST
# deprived (r -> 0) to MOST deprived (r -> 1). The absolute equity index is the
# population-weighted OLS slope of the per-capita benefit on the ridit score,
# i.e. the predicted per-capita difference between the extreme ends of the
# distribution (r = 1 minus r = 0) -- exactly the SII. The closed-form weighted
# slope beta = Sxy / Sxx is identical to
# `coef(lm(y ~ ridit, weights = N))["ridit"]` but is computed by pure
# data.table aggregation so it can be evaluated over millions of
# (mc, scenario, year, ...) groups without a per-group lm() call.
#
# Sign convention: all supported metrics (CPP, CYPP, DPP, net QALYs) are
# "goods", so a POSITIVE index means the most deprived gain more (pro-poor,
# inequality-reducing); a negative index is pro-rich (inequality-widening).
#
# @param d  A data.table with (at least) columns `rank` (integer deprivation
#   rank, 1 = MOST deprived .. K = LEAST deprived, i.e. `as.integer(dimd)`),
#   `B` (cumulative benefit for that group) and `N` (reference population for
#   that group, used both as the regression weight and the ridit denominator),
#   plus the grouping columns named in `by`.
# @param by  Character vector of columns identifying one gradient per group
#   (must include the Monte-Carlo id `mc`; DIMD is NOT in `by` -- it is the
#   gradient axis consumed into the index).
# @return A data.table with one row per `by` group and columns AEI_total
#   (SII rescaled to the whole reference population = modelled total benefit gap
#   between the extreme ends), AEI_per100k (SII expressed per 100,000 people),
#   REI_rel (SII / population-weighted mean; NA when the mean is 0) and
#   RII_ratio (fitted benefit at the most-deprived end / least-deprived end;
#   NA when the fitted line crosses zero over [0, 1]).
calc_equity_slope_indices <- function(d, by) {
  d <- copy(d)
  # Order least -> most deprived within each group (rank descending) so the
  # cumulative population share runs from the bottom to the top of the ridit.
  setorderv(d, c(by, "rank"), order = c(rep(1L, length(by)), -1L))
  d[, .totN := sum(N), by = by]
  d[, .csum := cumsum(N), by = by]
  d[, ridit := (.csum - N / 2) / .totN]

  agg <- d[, .(
    Sw   = sum(N),
    Swx  = sum(N * ridit),
    Swy  = sum(B),
    Swxx = sum(N * ridit * ridit),
    Swxy = sum(ridit * B)
  ), by = by]

  agg[, xbar := Swx / Sw]
  agg[, ybar := Swy / Sw]                       # population-weighted mean benefit
  agg[, Sxx := Swxx - Sw * xbar * xbar]
  agg[, Sxy := Swxy - Sw * xbar * ybar]
  agg[, beta := fifelse(Sxx > 0, Sxy / Sxx, NA_real_)]
  agg[, fit0 := ybar - beta * xbar]             # fitted per-capita benefit at r = 0 (least deprived)
  agg[, fit1 := ybar + beta * (1 - xbar)]       # fitted per-capita benefit at r = 1 (most deprived)

  agg[, AEI_total   := beta * Sw]
  agg[, AEI_per100k := beta * 1e5]
  agg[, REI_rel     := fifelse(ybar != 0, beta / ybar, NA_real_)]
  agg[, RII_ratio   := fifelse(!is.na(beta) & fit0 != 0 & sign(fit0) == sign(fit1),
                               fit1 / fit0, NA_real_)]

  agg[, .SD, .SDcols = c(by, "AEI_total", "AEI_per100k", "REI_rel", "RII_ratio")]
}


# Deprivation axes the model reports on, with their canonical factor levels
# (repeated across the package, e.g. Simulation_class.R, Disease_class.R,
# ExposureEffect_class.R). Order matters: as.integer() on these gives the
# deprivation rank 1 = most deprived .. K = least deprived.
.deprivation_vars <- c("dimd", "qimd")
.deprivation_levels <- list(
  dimd = c("1 most deprived", as.character(2:9), "10 least deprived"),
  qimd = c("1 most deprived", as.character(2:4), "5 least deprived")
)

# Standard decile -> quintile collapse (deciles 1,2 -> quintile 1, and so on),
# matching `Simulation_class.R` run_sim() and `Disease_class.R`.
.dimd_to_qimd <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L)


# deprivation_rank ----
# Deprivation rank (1 = most deprived .. K = least deprived) for a `dimd`/`qimd`
# column. factor() silently returns NA for labels outside `levels`, and callers
# use the rank only to order or collapse, so an unexpected level set would yield
# a plausible-but-wrong result rather than any error (decile and quintile labels
# even share four of five labels). Fail loudly instead.
deprivation_rank <- function(x, var) {
  levs <- .deprivation_levels[[var]]
  r <- as.integer(factor(as.character(x), levels = levs))
  if (anyNA(r)) {
    stop(sprintf(
      paste0("unexpected `%s` level(s) in the summaries: %s. Expected exactly: ",
             "%s. The deprivation rank drives the ordering and the decile -> ",
             "quintile collapse, so an unrecognised level would silently ",
             "produce a wrong result."),
      var,
      paste(sQuote(unique(as.character(x)[is.na(r)])), collapse = ", "),
      paste(sQuote(levs), collapse = ", ")
    ), call. = FALSE)
  }
  r
}


# add_qimd_from_dimd ----
# Adds a `qimd` (IMD quintile) column derived from `dimd` (IMD decile), in
# place, so ANY table family can be stratified by quintiles without `qimd`
# having been named in `strata_for_output`.
#
# This is a relabelling only -- the actual collapsing is then done by whatever
# group-by the caller already performs, which is what makes it exact: summing
# two deciles' counts (and popsize) reproduces precisely what aggregating by
# quintile at summary time would have given, and rates built as
# summed-numerator / summed-denominator follow automatically. The converse is
# impossible: deciles cannot be recovered from quintiles.
#
# No-op (returning FALSE) when `dt` is NULL, already has `qimd`, or has no
# `dimd` to derive from. Adding the column is harmless when unused: every
# aggregation in this file selects value columns by anchored pattern, so `qimd`
# is only ever carried into a result when it is named in the strata.
add_qimd_from_dimd <- function(dt) {
  if (is.null(dt) || "qimd" %in% names(dt) || !"dimd" %in% names(dt)) {
    return(invisible(FALSE))
  }
  lv <- .deprivation_levels$qimd
  dt[, qimd := factor(lv[.dimd_to_qimd][deprivation_rank(dimd, "dimd")],
                      levels = lv)]
  invisible(TRUE)
}


# Columns the equity tables manage themselves, so they cannot also be named as
# an output stratum: `mc` is the Monte-Carlo axis the indices are quantiled
# over, and `scenario` is always present because every index is a contrast
# against the comparator.
.equity_reserved_vars <- c("mc", "scenario")


# validate_equity_strata ----
# Each entry of the equity strata list is one output table. Any variable carried
# by the summaries (`strata_for_output` in the sim_design YAML) can be an output
# stratum -- the gradient is then fit *within* each stratum -- except that
# `dimd`/`qimd` are not stratifications at all: they select the deprivation axis
# the gradient is fit over. Whether a stratum variable actually exists can only
# be judged against the summaries, so that check lives in export_equity_tables();
# this validates the structural rules that hold regardless of the data.
validate_equity_strata <- function(strata) {
  if (!is.list(strata)) {
    stop("`strata$equity` must be a list of character vectors.", call. = FALSE)
  }
  for (i in seq_along(strata)) {
    s <- strata[[i]]
    lbl <- sprintf("`strata$equity[[%d]]`", i)
    if (!is.character(s) || length(s) == 0L) {
      stop(lbl, " must be a non-empty character vector.", call. = FALSE)
    }
    if (anyDuplicated(s)) {
      stop(lbl, " has duplicated entries: ",
           paste(unique(s[duplicated(s)]), collapse = ", "), ".", call. = FALSE)
    }
    bad <- intersect(s, .equity_reserved_vars)
    if (length(bad) > 0L) {
      stop(lbl, " names reserved column(s): ", paste(bad, collapse = ", "),
           ". `mc` is the Monte-Carlo axis the indices are quantiled over and",
           " `scenario` is always included, so neither can be an output",
           " stratum.", call. = FALSE)
    }
    if (!"year" %in% s) {
      stop(lbl, " must include \"year\". The benefit is accumulated over years",
           " and each output row reports the benefit cumulated to that year,",
           " against that year's reference population -- so `year` defines the",
           " rows rather than merely subsetting them.", call. = FALSE)
    }
  }
  invisible(TRUE)
}


# build_equity_plans ----
# Expands an equity strata list into one plan per output table. `dimd`/`qimd`
# in an entry select the gradient axis; an entry naming both yields one table
# per axis. An entry naming neither falls back to `dimd`, which keeps the
# historical defaults -- and their filenames -- unchanged.
build_equity_plans <- function(strata) {
  unlist(lapply(strata, function(s) {
    out_vars <- setdiff(s, .deprivation_vars)   # keeps the caller's order
    grads <- intersect(.deprivation_vars, s)    # canonical dimd, qimd order
    implicit <- length(grads) == 0L
    if (implicit) grads <- "dimd"
    lapply(grads, function(g) list(
      out_vars = out_vars,
      gradient = g,
      # The gradient token is echoed in the filename only when the caller wrote
      # it, so `list("year", c("year", "sex"))` still produces exactly the
      # filenames it always has. `agegrp` is spelled `agegroup` in the suffix,
      # matching every other table (see make_strata_configs()).
      suffix = gsub(
        "agegrp", "agegroup",
        paste(if (implicit) out_vars else c(out_vars, g), collapse = "-")
      )
    ))
  }), recursive = FALSE)
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
#   - equity: List of character vectors for equity slope-index tables. `dimd`
#     and `qimd` mean something different here -- see below -- and the list is
#     validated.
#
#   Valid stratification variables:
#   - year: Simulation year (always required)
#   - sex: Sex (men/women)
#   - agegrp: Age groups (5-year bands)
#   - dimd: IMD deciles (10 levels: "1 most deprived" to "10 least deprived")
#   - agegrp20: 20-year age groups (xps tables only)
#   - qimd: IMD quintiles (5 levels)
#
#   `qimd` works in EVERY table family and does NOT need to be in
#   `strata_for_output`: whenever the summaries carry `dimd` it is derived by
#   the standard decile-pair collapse (see add_qimd_from_dimd()), which is
#   exact. The converse does not hold -- `dimd` must really be in the
#   summaries, and is unavailable in the xps tables, which are written
#   quintile-based by export_xps().
#
#   For `equity`, any variable the summaries carry can be an output stratum
#   (the gradient is then fit within it), and `year` is required. The
#   exception is `dimd`/`qimd`, which do NOT stratify -- they select the
#   deprivation gradient the index is fit over. Naming both in one entry writes
#   one table per gradient; naming neither falls back to `dimd`.
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
#
# # Every table family by IMD quintiles. `qimd` need not be in
# # strata_for_output -- it is derived from dimd whenever the summaries carry it
# sim$export_tables(strata = list(
#   ons          = list("year", c("year", "qimd"), c("year", "agegrp", "sex", "qimd")),
#   esp          = list("year", c("year", "qimd")),
#   mrtl_ons     = list("year", c("year", "qimd")),
#   mrtl_esp     = list("year", c("year", "qimd")),
#   disease_char = list("year", c("year", "qimd")),
#   equity       = list(c("year", "qimd"))
# ))
#
# # Equity indices over IMD quintiles instead of deciles, plus both gradients
# # by year, and a gradient fit within each age group
# sim$export_tables(strata = list(
#   equity = list(c("year", "dimd", "qimd"), c("year", "sex", "qimd"),
#                 c("year", "agegrp"))
# ))
# # -> equity cpp slope index by year-dimd (not standardised).csv
# #    equity cpp slope index by year-qimd (not standardised).csv
# #    equity cpp slope index by year-sex-qimd (not standardised).csv
# #    equity cpp slope index by year-agegroup (not standardised).csv
#' @description
#' Export summary tables for policy analysis.
#'
#' Builds main, all-cause mortality, disease-characteristics, exposure,
#' (optionally) cost-effectiveness and (optionally) equity slope-index tables
#' from the per-summary files produced by `$export_summaries()`. Outputs are
#' written to `output_dir/tables/` (or `output_dir/tables2agegrps/` when
#' `two_agegrps = TRUE`).
#'
#' When `equity = TRUE`, equity slope-index tables are written as
#' `equity <metric> slope index by <strata> (not standardised).csv` for each of
#' CPP, CYPP, DPP and net QALYs. These summarise how each cumulative benefit is
#' distributed across deprivation groups -- DIMD deciles by default, or IMD
#' quintiles (`qimd`) if asked for via `strata$equity` -- using absolute and
#' relative analogues of the Slope Index of Inequality and Relative Index of
#' Inequality (see the `equity` argument).
#'
#' When `cea = TRUE`, cost-effectiveness tables are written as
#' `cost-effectiveness by <strata> (<perspective>-<scale>) (not standardised).csv`,
#' containing the cumulative discounted incremental QALYs (`dQALYs_cuml`) and
#' costs (`dCosts_cuml`), the incremental cost-effectiveness ratio (`ICER`), and
#' the net monetary benefit (`NMB_at_wtp_<threshold>`) at each willingness-to-pay
#' threshold, for each QALY scale (EQ5D5L) and from two perspectives:
#' \emph{societal} (uses `total_cost`) and \emph{healthcare} (uses
#' `healthcare_cost`). User-defined `*_costs` columns are always added to the
#' societal perspective and, if named in `custom_costs_in_healthcare`, to the
#' healthcare perspective as well.
#'
#' @param baseline_year_for_change_outputs Integer. Reference year used for
#'   computing change-from-baseline columns and the start year for the
#'   cumulative cost-effectiveness calculations. Two-digit values (e.g. `19`)
#'   are auto-promoted to four-digit (`2019`).
#' @param prbl Numeric vector of probability levels for output quantiles
#'   (median plus uncertainty bounds). Default
#'   `c(0.5, 0.025, 0.975, 0.1, 0.9)`.
#' @param comparator_scenario Character. Name of the scenario used as the
#'   comparator when computing differences between scenarios.
#' @param two_agegrps Logical. If `TRUE`, uses a coarser two-age-group
#'   stratification; otherwise uses the standard age groups.
#' @param strata Optional named list overriding the default stratification
#'   configuration. See examples in the file header for shape; passed
#'   through `private$build_strata_config()`.
#'
#'   Any element may be stratified by `"qimd"` (IMD quintiles) instead of
#'   `"dimd"` (deciles), and `qimd` does **not** need to be listed in
#'   `strata_for_output` in the sim_design YAML: whenever the loaded summaries
#'   carry `dimd`, a `qimd` column is derived from it by the standard
#'   decile-pair collapse (deciles 1-2 -> quintile 1, and so on). That is exact
#'   rather than approximate, because every aggregation here groups first and
#'   then sums, so summing two deciles reproduces precisely what aggregating by
#'   quintile at summary time would have given -- and rates, which are formed as
#'   summed-numerator / summed-denominator after the group-by, follow
#'   automatically. The converse is impossible: deciles cannot be recovered from
#'   quintiles, so a `dimd` stratum requires `dimd` in the summaries, and is
#'   unavailable in the `xps_*` elements because `$export_xps()` writes those
#'   datasets quintile-based already.
#'
#'   In the `equity` element only, `"dimd"`/`"qimd"` select the deprivation
#'   gradient rather than a stratification, and `"year"` is required -- see the
#'   `equity` argument.
#' @param multicore Logical. If `TRUE`, runs table-building tasks in
#'   parallel with single-threaded workers; otherwise runs sequentially.
#' @param cea Logical. If `TRUE` (default), also build the cost-effectiveness
#'   (ICER / NMB) tables from the `qalys` and `costs` summaries.
#' @param wtp Numeric vector of willingness-to-pay thresholds (currency per
#'   QALY) at which the net monetary benefit is computed. Default
#'   `c(20000, 30000)` (the NICE thresholds, GBP/QALY).
#' @param qaly_discount_rate,cost_discount_rate Numeric. Annual discount rates
#'   (percent) applied to QALYs and costs. Default `3.5` each (UK NICE
#'   guidance). These arguments are the single source of truth for discounting;
#'   there is no `discounting` block in `sim_design.yaml`. Discounting is applied
#'   **only** to the cost-effectiveness tables - the prevalence, incidence,
#'   mortality, QALY and cost tables are reported undiscounted.
#' @param discount_from_year Integer. First year from which present values are
#'   discounted. When `NULL` (default) it defaults to
#'   `baseline_year_for_change_outputs`.
#' @param custom_costs_in_healthcare Character vector of user-defined `*_costs`
#'   column names to add to the healthcare perspective (in addition to the
#'   always-present `healthcare_cost`). `NULL`/`FALSE` (default) adds none;
#'   `TRUE` adds all user-defined cost columns. User-defined cost columns are
#'   always included in the societal perspective regardless of this argument.
#' @param equity Logical. If `TRUE` (default), also build the equity
#'   slope-index tables (absolute and relative analogues of the Slope Index of
#'   Inequality / Relative Index of Inequality) for the cumulative CPP, CYPP,
#'   DPP and net-QALYs benefits, distributed across deprivation groups.
#'   Written as `equity <metric> slope index by <strata> (not standardised).csv`
#'   with a `type` column giving four indices: `AEI_total` (absolute equity
#'   index rescaled to the whole population = modelled total benefit gap between
#'   the most- and least-deprived ends), `AEI_per100k` (the same slope per
#'   100,000 people), `REI_rel` (slope divided by the population-weighted mean
#'   benefit) and `RII_ratio` (fitted benefit at the most-deprived end divided
#'   by the least-deprived end). By the pro-poor sign convention a positive
#'   index means the most deprived gain more (inequality-reducing). The
#'   deprivation gradient is fit per Monte-Carlo iteration and quantiled across
#'   iterations, following the SII/RII framework (Moreno-Betancur et al. 2015;
#'   Wagstaff et al. 1991) and, for cumulative benefits, Cookson et al. 2026.
#'   A gradient requires at least two deprivation groups, so scenario-years that
#'   cover only one (e.g. interventions targeted at a single decile) have an
#'   undefined index and are omitted; the omission is reported when `logs` are
#'   enabled in the design.
#'
#'   **Output strata (`strata$equity`).** Each entry of the list is one output
#'   table, and any variable the summaries carry (`strata_for_output` in the
#'   sim_design YAML: `agegrp`, `sex`, `ethnicity`, `sha`, ...) can be an output
#'   stratum -- the gradient is then fit *within* each stratum, exactly as
#'   filtering the summaries to that stratum and fitting without it would. Every
#'   entry must name `"year"`, because the benefit is accumulated over years and
#'   each row reports the benefit cumulated to that year against that year's
#'   reference population, so `year` defines the rows rather than subsetting
#'   them. `"mc"` and `"scenario"` are reserved. Variables the summaries do not
#'   carry raise a `warning()` and that table is skipped.
#'
#'   **Choosing the gradient axis (`dimd` vs `qimd`).** The gradient is fit over
#'   DIMD deciles by default. Name `"dimd"` and/or `"qimd"` in an entry of
#'   `strata$equity` to choose explicitly; unlike every other variable these do
#'   *not* stratify, because the deprivation axis is consumed into the index. An
#'   entry naming both writes one table per axis, and an entry naming neither
#'   falls back to `dimd` so the built-in defaults keep producing exactly the
#'   files they always have. The chosen axis is echoed in a `gradient` column in
#'   every file, and (when the caller wrote it) in the filename suffix:
#'   `list(c("year", "sex", "qimd"))` writes
#'   `equity cpp slope index by year-sex-qimd (not standardised).csv`. As in
#'   every other table family, `agegrp` is spelled `agegroup` in the suffix.
#'   Quintiles give a coarser, lower-variance gradient with less power to detect
#'   departures from linearity; deciles are usually preferable when available.
#'   `qimd` does **not** need to be in `strata_for_output` -- if the summaries
#'   carry `dimd` it is derived by the standard decile-pair collapse, which is
#'   exact because event counts and populations are additive. The converse does
#'   not hold: a `dimd` gradient requires `dimd` in the summaries. If the
#'   requested axis is unavailable, or the summaries carry no deprivation column
#'   at all, a `warning()` is raised (not a `logs`-gated message) naming the
#'   tables that were skipped.
#' @param equity_ridit_reference Character, one of `"comparator"` (default) or
#'   `"scenario"`. Chooses the population used to build the ridit (deprivation)
#'   ranks and the population weights for the equity slope-index regression --
#'   the reference against which each decile's share of the population is
#'   measured. `"comparator"` uses the comparator scenario's decile populations,
#'   so the deprivation ranking is *identical* across the scenarios being
#'   compared and across years; `"scenario"` uses each intervention scenario's
#'   own decile populations, reflecting the socioeconomic composition that
#'   scenario actually produces.
#'
#'   **Which to choose (Renard et al. 2019).** Renard et al. showed that the SII
#'   and RII can move purely because the socioeconomic *composition* of the
#'   population shifts, even when group-specific rates are unchanged -- so a
#'   changing denominator population can masquerade as a change in inequality.
#'   For the usual question here -- *how does an intervention redistribute
#'   benefits across a fixed deprivation structure?* -- keep the default
#'   `"comparator"`: a common, fixed reference makes between-scenario and
#'   between-year comparisons clean and attributes differences to the benefits,
#'   not to composition. Choose `"scenario"` only when the intervention itself is
#'   expected to change the deprivation composition (e.g. large differential
#'   survival shifting who is alive in each decile) and you specifically want the
#'   index to reflect each scenario's own realised population; interpret
#'   between-scenario differences with the Renard caveat in mind. In this model
#'   DIMD deciles are structurally near-fixed, so the two options usually differ
#'   only slightly.
#' @return The `Simulation` object, invisibly.
Simulation$set("public", "export_tables", function(
    baseline_year_for_change_outputs = 2019L,
    prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
    comparator_scenario = "sc0",
    two_agegrps = FALSE,
    strata = NULL,
    multicore = TRUE,
    cea = TRUE,
    wtp = c(20000, 30000),
    qaly_discount_rate = 3.5,
    cost_discount_rate = 3.5,
    discount_from_year = NULL,
    custom_costs_in_healthcare = NULL,
    equity = TRUE,
    equity_ridit_reference = c("comparator", "scenario")
) {
  equity_ridit_reference <- match.arg(equity_ridit_reference)

  # Ensure baseline year is in full format (e.g. 2019, not 19)
  # Data is converted to full year format in export_main_tables()
  if (baseline_year_for_change_outputs <= 100) {
    baseline_year_for_change_outputs <- baseline_year_for_change_outputs + 2000L
  }

  # Discounting is applied ONLY in the cost-effectiveness (CEA) tables, and is
  # controlled solely by the arguments below (there is intentionally no
  # `discounting` block in sim_design.yaml). When `discount_from_year` is NULL
  # it defaults to the baseline/reference year inside export_cea_tables(). The
  # main prevalence / incidence / mortality / QALY / cost tables are reported
  # undiscounted.

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

  # Cost-effectiveness (ICER / NMB) tables, built from qalys + costs summaries.
  if (cea) {
    tasks[[length(tasks) + 1L]] <- list(
      id = 5L,
      type = "cea",
      prbl = prbl,
      summaries_dir = private$output_dir("summaries"),
      tables_dir = tables_dir,
      comparator_scenario = comparator_scenario,
      baseline_year = baseline_year_for_change_outputs,
      wtp = wtp,
      qaly_discount_rate = qaly_discount_rate,
      cost_discount_rate = cost_discount_rate,
      discount_from_year = discount_from_year,
      custom_costs_in_healthcare = custom_costs_in_healthcare,
      strata = strata_cfg$ons
    )
  }

  # Equity slope-index (SII / RII analogue) tables, built from the prvl / incd /
  # mrtl / qalys summaries. The deprivation gradient is fit across DIMD deciles.
  if (equity) {
    tasks[[length(tasks) + 1L]] <- list(
      id = 6L,
      type = "equity",
      prbl = prbl,
      summaries_dir = private$output_dir("summaries"),
      tables_dir = tables_dir,
      comparator_scenario = comparator_scenario,
      baseline_year = baseline_year_for_change_outputs,
      ridit_reference = equity_ridit_reference,
      strata = strata_cfg$equity,
      two_agegrps = two_agegrps
    )
  }

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
    ),
    "cea" = private$export_cea_tables(
      prbl = task$prbl,
      summaries_dir = task$summaries_dir,
      tables_dir = task$tables_dir,
      comparator_scenario = task$comparator_scenario,
      baseline_year = task$baseline_year,
      wtp = task$wtp,
      qaly_discount_rate = task$qaly_discount_rate,
      cost_discount_rate = task$cost_discount_rate,
      discount_from_year = task$discount_from_year,
      custom_costs_in_healthcare = task$custom_costs_in_healthcare,
      strata = task$strata
    ),
    "equity" = private$export_equity_tables(
      prbl = task$prbl,
      summaries_dir = task$summaries_dir,
      tables_dir = task$tables_dir,
      comparator_scenario = task$comparator_scenario,
      baseline_year = task$baseline_year,
      ridit_reference = task$ridit_reference,
      strata = task$strata,
      two_agegrps = task$two_agegrps
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
      ),
      # Equity slope-index output strata. DIMD deciles are always the gradient
      # axis (consumed into the index), so it never appears here.
      equity = list(
        "year",
        c("year", "sex")
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
      ),
      # Equity slope-index output strata. DIMD deciles are always the gradient
      # axis (consumed into the index), so it never appears here.
      equity = list(
        "year",
        c("year", "sex")
      )
    )
  }

  # Merge user-provided strata with defaults (user overrides defaults)
  result <- defaults
  if (!is.null(user_strata)) {
    for (name in names(user_strata)) {
      if (name %in% names(defaults)) {
        result[[name]] <- user_strata[[name]]
      }
    }
  }

  # Equity strata accept a much smaller vocabulary than the rest, and silently
  # ignoring an unsupported variable produced a mislabelled duplicate table.
  # Validate here so it fails fast, before any table-building work starts.
  validate_equity_strata(result$equity)

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
    two_agegrps = FALSE
) {
  # NOTE: the main tables (qalys, costs, net_qalys, net_costs, etc.) are
  # reported UNDISCOUNTED. Discounting is applied only in the cost-effectiveness
  # tables (export_cea_tables), controlled by export_tables() arguments.
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
    "costs" = "_cost$|_costs$|^economic_output$", "net_costs" = "_cost$|_costs$|^economic_output$",
    "contd" = "_contd$", "contd_change_relative" = "_contd$", "contd_change_absolute" = "_contd$",
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
    "contd" = "contd_mean_", "contd_change_relative" = "contd_change_relative_", "contd_change_absolute" = "contd_abs_change_",
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
    "contd" = "continuous outcome by ",
    "contd_change_relative" = "continuous outcome relative change by ",
    "contd_change_absolute" = "continuous outcome absolute change by ",
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
      d <- tt[, lapply(.SD, sum), .SDcols = patterns("_cost$|_costs$|^economic_output$"), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "costs_type", value.name = "costs")
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
      d <- tt[, lapply(.SD, sum), .SDcols = patterns("_cost$|_costs$|^economic_output$"), keyby = eval(x)]
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
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = TRUE),
             keyby = eval(setdiff(c(x, "costs_type"), "mc"))]
      x <- c(x, "costs_type")
      setnames(d, c(setdiff(x, "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(x, "mc")))
      setcolorder(d, setdiff(x, "mc"))

    } else if (grepl("^contd", what)) {
      # User-defined continuous outcomes (*_contd columns). The contd summaries
      # already hold a population-weighted mean per stratum, so collapsing to
      # coarser strata re-weights by popsize (NOT a sum, and no division by
      # popsize).
      contd_cols <- grep("_contd$", names(tt), value = TRUE)
      d <- tt[, lapply(.SD, function(v) weighted.mean(v, popsize, na.rm = TRUE)),
              .SDcols = contd_cols, keyby = eval(x)]
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

      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable), rounding = FALSE),
             keyby = eval(setdiff(x, "mc"))]
      setnames(d, c(setdiff(x, "mc"), "outcome", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, setdiff(x, "mc"))
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
    costs = c("costs", "net_costs"),
    # User-defined continuous outcomes (*_contd columns). Skipped silently when
    # the contd summary was not produced (read_summary_dataset returns NULL).
    contd = c("contd", "contd_change_relative", "contd_change_absolute")
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
      # Allow `qimd` strata even when the summaries only carry `dimd`. Every
      # aggregation below selects value columns by anchored pattern, so the
      # extra column is inert unless it is named in the strata.
      add_qimd_from_dimd(tt_base)

      # Also load prvl for ftlt metrics (case fatality denominator)
      prvl_for_ftlt <- NULL
      if (source_name == "dis_mrtl") {
        prvl_for_ftlt <- private$read_summary_dataset("prvl", pop_key)
        if (!is.null(prvl_for_ftlt)) {
          prvl_for_ftlt[, year := year + 2000L]
          add_qimd_from_dimd(prvl_for_ftlt)
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

        # Generate tables. Main tables are reported UNDISCOUNTED; discounting is
        # applied only in the cost-effectiveness tables (export_cea_tables).
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
        # Determine what was standardised by (variables NOT in strata).
        # `qimd` retains the deprivation axis just as `dimd` does (coarser),
        # so a qimd stratum must not be reported as standardised over dimd.
        possible_vars <- c("age", "sex", "dimd")
        s_axes <- if ("qimd" %in% s) c(s, "dimd") else s
        standardised_vars <- setdiff(possible_vars, s_axes)
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
  # Allow `qimd` strata even when the summaries only carry `dimd`
  add_qimd_from_dimd(tt_scaled)
  add_qimd_from_dimd(pp_scaled)
  add_qimd_from_dimd(tt_esp)

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
  # Allow `qimd` strata even when the summaries only carry `dimd`
  add_qimd_from_dimd(tt)

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

    # Identify exposure columns: _curr_xps columns plus any extra from exposures_for_output
    xps_cols <- grep("_curr_xps$", names(xps_tab), value = TRUE)
    extra_in_data <- intersect(self$design$sim_prm$exposures_for_output, names(xps_tab))
    extra_in_data <- extra_in_data[
      vapply(extra_in_data, function(col) is.numeric(xps_tab[[col]]), logical(1))
    ]
    xps_cols <- unique(c(xps_cols, extra_in_data))

    # Detect which columns have "All" marginals from groupingsets
    filt_vars <- detect_filterable_vars(xps_tab)
    strata_configs <- make_xps_strata_configs(strata_ons, filt_vars, standardised = FALSE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- xps_tab[eval(cfg$filter_expr)]
      if (nrow(d) == 0) next
      d <- d[, lapply(.SD, mean), .SDcols = xps_cols, keyby = eval(outstrata)]
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

    # Identify exposure columns: _curr_xps columns plus any extra from exposures_for_output
    xps_cols <- grep("_curr_xps$", names(xps_tab), value = TRUE)
    extra_in_data <- intersect(self$design$sim_prm$exposures_for_output, names(xps_tab))
    extra_in_data <- extra_in_data[
      vapply(extra_in_data, function(col) is.numeric(xps_tab[[col]]), logical(1))
    ]
    xps_cols <- unique(c(xps_cols, extra_in_data))

    # Detect which columns have "All" marginals from groupingsets
    filt_vars <- detect_filterable_vars(xps_tab)
    strata_configs <- make_xps_strata_configs(strata_esp, filt_vars, standardised = TRUE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- xps_tab[eval(cfg$filter_expr)]
      if (nrow(d) == 0) next
      d <- d[, lapply(.SD, mean), .SDcols = xps_cols, keyby = eval(outstrata)]
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


# export_cea_tables ----
# Generate cost-effectiveness (ICER / NMB) tables from the qalys and costs
# summaries. For each stratum, perspective (societal / healthcare) and QALY
# scale (EQ5D5L; HUI3 if present) it computes, per Monte-Carlo iteration, the
# cumulative discounted incremental QALYs and costs versus the comparator
# scenario, then the ICER and the net monetary benefit (NMB) at each
# willingness-to-pay threshold, and quantiles those across iterations.
#
# Cost perspectives:
#   societal   = total_cost      (+ all user *_costs columns)
#                total_cost already nets economic_output out of
#                healthcare + socialcare + informalcare costs.
#   healthcare = healthcare_cost (+ the user *_costs columns named in
#                custom_costs_in_healthcare) -- direct treatment costs only.
#
# Discounting: PV = FV / (1 + rate/100)^max(0, year - discount_from_year),
# with separate rates for QALYs and costs. Actual (scaled_up) population only.
Simulation$set("private", "export_cea_tables", function(
    prbl,
    summaries_dir,
    tables_dir,
    comparator_scenario = "sc0",
    baseline_year = 2019L,
    wtp = c(20000, 30000),
    qaly_discount_rate = 3.5,
    cost_discount_rate = 3.5,
    discount_from_year = NULL,
    custom_costs_in_healthcare = NULL,
    strata = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating cost-effectiveness (ICER/NMB) tables...")
  }

  if (is.null(discount_from_year)) discount_from_year <- baseline_year

  qalys <- private$read_summary_dataset("qalys", "scaled_up")
  costs <- private$read_summary_dataset("costs", "scaled_up")
  if (is.null(qalys) || is.null(costs)) {
    if (self$design$sim_prm$logs) {
      message("  qalys or costs summary missing; skipping CEA tables")
    }
    return(invisible(NULL))
  }

  # Summaries store year in short format (e.g. 19); promote to full (2019) so
  # baseline-year filtering and discounting use the same scale as the args.
  qalys[, year := year + 2000L]
  costs[, year := year + 2000L]
  # Allow `qimd` strata even when the summaries only carry `dimd`
  add_qimd_from_dimd(qalys)
  add_qimd_from_dimd(costs)

  # Need at least one intervention scenario to compare against the comparator
  non_comparator <- setdiff(unique(qalys$scenario), comparator_scenario)
  if (length(non_comparator) == 0L) {
    if (self$design$sim_prm$logs) {
      message("  no intervention scenarios (only '", comparator_scenario,
              "' found); skipping CEA tables")
    }
    return(invisible(NULL))
  }

  # QALY scales actually present in the summary (England: EQ5D5L only, but
  # HUI3 is honoured if a future summary provides it).
  scales_avail <- intersect(c("EQ5D5L", "HUI3"), names(qalys))
  if (length(scales_avail) == 0L) {
    if (self$design$sim_prm$logs) {
      message("  no EQ5D5L/HUI3 columns in qalys summary; skipping CEA tables")
    }
    return(invisible(NULL))
  }

  # User-defined cost columns use the plural `_costs$` suffix; the built-in cost
  # columns use the singular `_cost` suffix, so this grep matches only the
  # user-defined columns.
  custom_cost_cols <- grep("_costs$", names(costs), value = TRUE)

  # Resolve which custom cost columns to add to the healthcare perspective.
  # `custom_costs_in_healthcare` accepts a character vector of (custom) cost
  # column names to include there, in addition to the always-present
  # healthcare_cost. For convenience a logical is also honoured:
  # NULL/FALSE -> none (default), TRUE -> all user-defined custom cost columns.
  if (is.null(custom_costs_in_healthcare) ||
      isFALSE(custom_costs_in_healthcare)) {
    healthcare_custom_cols <- character(0)
  } else if (isTRUE(custom_costs_in_healthcare)) {
    healthcare_custom_cols <- custom_cost_cols
  } else {
    requested <- as.character(custom_costs_in_healthcare)
    healthcare_custom_cols <- intersect(requested, custom_cost_cols)
    unknown <- setdiff(requested, custom_cost_cols)
    if (length(unknown) > 0L && self$design$sim_prm$logs) {
      message("  custom_costs_in_healthcare: ignoring name(s) not matching a ",
              "user-defined cost column: ", paste(unknown, collapse = ", "))
    }
  }

  perspective_cols <- list(
    societal = c("total_cost", custom_cost_cols),
    healthcare = c("healthcare_cost", healthcare_custom_cols)
  )

  # WTP -> NMB column-name labels (plain integer form, e.g. 20000)
  wtp_labels <- paste0(
    "NMB_at_wtp_",
    vapply(wtp, function(w) format(w, scientific = FALSE, trim = TRUE,
                                   big.mark = ""), character(1))
  )

  disc <- function(v, year, rate) {
    v / (1 + rate / 100)^pmax(0, year - discount_from_year)
  }

  for (s in strata) {
    x <- c("mc", "scenario", s) # s always contains "year"

    for (persp in names(perspective_cols)) {
      pcols <- intersect(perspective_cols[[persp]], names(costs))
      if (!any(grepl("^(total_cost|healthcare_cost)$", pcols))) {
        if (self$design$sim_prm$logs) {
          message("  ", persp, ": required cost column missing; skipping")
        }
        next
      }

      # Aggregate (discounted) costs for this perspective, once per stratum
      cc <- copy(costs)
      cc[, .cost := Reduce(`+`, .SD), .SDcols = pcols]
      cc <- cc[, .(C = sum(.cost)), keyby = eval(x)]
      cc[, C := disc(C, year, cost_discount_rate)]

      for (scale in scales_avail) {
        # Aggregate (discounted) QALYs for this scale
        qq <- qalys[, .(Q = sum(get(scale))), keyby = eval(x)]
        qq[, Q := disc(Q, year, qaly_discount_rate)]

        d <- merge(qq, cc, by = x, all = TRUE)
        d[is.na(Q), Q := 0][is.na(C), C := 0]

        # Incremental vs comparator (intervention - comparator), from baseline
        cmp <- d[scenario == comparator_scenario & year >= baseline_year][
          , scenario := NULL]
        d <- d[scenario != comparator_scenario & year >= baseline_year]
        if (nrow(d) == 0L || nrow(cmp) == 0L) next
        d[cmp, on = setdiff(x, "scenario"), `:=`(dQ = Q - i.Q, dC = C - i.C)]
        d <- d[!is.na(dQ) & !is.na(dC)]
        if (nrow(d) == 0L) next

        # Cumulative over year within (mc, scenario, other strata)
        setkeyv(d, c(setdiff(x, "year"), "year"))
        d[, `:=`(dQALYs_cuml = cumsum(dQ), dCosts_cuml = cumsum(dC)),
          by = setdiff(x, "year")]

        # ICER and NMB at each WTP
        d[, ICER := fifelse(dQALYs_cuml == 0, NA_real_,
                            dCosts_cuml / dQALYs_cuml)]
        for (i in seq_along(wtp)) {
          set(d, NULL, wtp_labels[i], wtp[i] * d$dQALYs_cuml - d$dCosts_cuml)
        }

        metric_cols <- c("dCosts_cuml", "dQALYs_cuml", "ICER", wtp_labels)
        dm <- melt(d, id.vars = x, measure.vars = metric_cols,
                   variable.name = "type", value.name = "value")
        # Drop non-finite draws (e.g. ICER when dQALYs_cuml == 0) so the
        # quantile is taken over the finite Monte-Carlo iterations only.
        dm <- dm[is.finite(value)]
        if (nrow(dm) == 0L) next

        setkey(dm, "type")
        out <- dm[, safe_fquantile_byid(value, prbl, id = as.character(type),
                                        rounding = FALSE),
                  keyby = eval(setdiff(x, "mc"))]
        setnames(out, c(setdiff(x, "mc"), "type",
                        scales::percent(prbl, prefix = "value_")))
        setkeyv(out, c("type", setdiff(x, "mc")))
        setcolorder(out, setdiff(x, "mc"))

        suffix <- paste(setdiff(x, c("mc", "scenario")), collapse = "-")
        fwrite(out, file.path(
          tables_dir,
          paste0("cost-effectiveness by ", suffix,
                 " (", persp, "-", scale, ") (not standardised).csv")
        ))
        rm(d, dm, out, qq)
      }
      rm(cc)
    }
  }

  rm(qalys, costs)
  invisible(NULL)
})


# export_equity_tables ----
# Generate equity slope-index tables (absolute and relative analogues of the
# Slope Index of Inequality / Relative Index of Inequality) for the cumulative
# policy benefits CPP, CYPP, DPP and net QALYs, distributed across DIMD
# deprivation deciles. See calc_equity_slope_indices() for the statistical
# framework and sign convention.
#
# For each metric and each requested output stratum (e.g. "year" or
# c("year", "sex")) the deprivation gradient is fit PER Monte-Carlo iteration
# across the 10 DIMD deciles, and the four indices are then quantiled across the
# Monte-Carlo iterations (the MC iterations ARE the uncertainty), mirroring
# every other table. Uses the actual (scaled_up) population only -- "not
# standardised" -- because age-standardising across deciles would defeat the
# purpose of describing the real deprivation distribution.
Simulation$set("private", "export_equity_tables", function(
    prbl,
    summaries_dir,
    tables_dir,
    comparator_scenario = "sc0",
    baseline_year = 2019L,
    ridit_reference = "comparator",
    strata = NULL,
    two_agegrps = FALSE
) {
  ridit_reference <- match.arg(ridit_reference, c("comparator", "scenario"))
  if (self$design$sim_prm$logs) {
    message("Generating equity slope-index tables (ridit reference: ",
            ridit_reference, ")...")
  }

  if (is.null(strata)) strata <- list("year", c("year", "sex"))
  validate_equity_strata(strata)
  # One plan == one output table: the output stratification, the deprivation
  # gradient axis it is fit over, and the filename suffix.
  plans <- build_equity_plans(strata)
  wanted_grads <- unique(vapply(plans, `[[`, character(1), "gradient"))

  # Metric -> source dataset, value-column pattern, suffix to strip for the
  # disease label, and benefit direction. "prevented" = comparator - intervention
  # (CPP/CYPP/DPP), "gained" = intervention - comparator (net QALYs).
  metric_cfg <- list(
    cpp       = list(source = "incd",  pattern = "_incd$",   strip = "_incd$",  dir = "prevented"),
    cypp      = list(source = "prvl",  pattern = "_prvl$",   strip = "_prvl$",  dir = "prevented"),
    dpp       = list(source = "mrtl",  pattern = "_mrtl$",   strip = "_mrtl$",  dir = "prevented"),
    net_qalys = list(source = "qalys", pattern = "^EQ5D5L$", strip = NA,        dir = "gained")
  )

  for (metric in names(metric_cfg)) {
    cfg <- metric_cfg[[metric]]

    tt <- private$read_summary_dataset(cfg$source, "scaled_up")
    if (is.null(tt)) {
      if (self$design$sim_prm$logs) {
        message("  ", metric, ": ", cfg$source, " summary missing; skipping")
      }
      next
    }
    # A missing deprivation column means NO equity table at all, so say so
    # unconditionally rather than under `logs` -- an otherwise silent no-op is
    # indistinguishable from the feature having run.
    if (!any(.deprivation_vars %in% names(tt))) {
      warning(sprintf(
        paste0("equity tables: the `%s` summary has no deprivation column ",
               "(looked for %s), so no equity %s tables were written. Add ",
               "`dimd` to `strata_for_output` in the sim_design YAML and re-run ",
               "$export_summaries()."),
        cfg$source, paste(.deprivation_vars, collapse = "/"), metric),
        call. = FALSE, immediate. = TRUE)
      rm(tt); next
    }
    valcols <- grep(cfg$pattern, names(tt), value = TRUE)
    if (length(valcols) == 0L) {
      if (self$design$sim_prm$logs) {
        message("  ", metric, ": no value columns matching ", cfg$pattern, "; skipping")
      }
      rm(tt); next
    }

    # Promote short year (19) to full (2019) to match baseline_year.
    tt[, year := year + 2000L]

    # Match export_main_tables(): under two_agegrps, `agegrp` means the coarse
    # 30-64 / 65-99 split. Relevant now that agegrp can be an equity stratum.
    if (two_agegrps && "agegrp" %in% names(tt)) {
      tt[agegrp %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59",
                       "60-64"), agegrp := "30-64"]
      tt[agegrp %in% c("65-69", "70-74", "75-79", "80-84", "85-89", "90-94",
                       "95-99"), agegrp := "65-99"]
    }

    # The index is a benefit-vs-comparator contrast, so both the comparator and
    # at least one intervention scenario must be present.
    scns <- unique(tt$scenario)
    if (!comparator_scenario %in% scns) {
      warning(sprintf(
        paste0("equity tables: comparator scenario '%s' is not in the `%s` ",
               "summary (found: %s), so no equity %s tables were written."),
        comparator_scenario, cfg$source, paste(scns, collapse = ", "), metric),
        call. = FALSE, immediate. = TRUE)
      rm(tt); next
    }
    non_comparator <- setdiff(scns, comparator_scenario)
    if (length(non_comparator) == 0L) {
      if (self$design$sim_prm$logs) {
        message("  ", metric, ": no intervention scenarios (only '",
                comparator_scenario, "'); skipping")
      }
      rm(tt); next
    }

    # Derive `qimd` on demand so a quintile gradient does not require `qimd` in
    # `strata_for_output` (see add_qimd_from_dimd()).
    if ("qimd" %in% wanted_grads && add_qimd_from_dimd(tt) &&
        self$design$sim_prm$logs) {
      message("  ", metric, ": derived qimd quintiles from dimd deciles")
    }

    for (p in plans) {                        # out_vars always contains "year"
      grad <- p$gradient
      if (!grad %in% names(tt)) {
        warning(sprintf(
          paste0("equity tables: gradient '%s' was requested but the `%s` ",
                 "summary has no `%s` column and it cannot be derived from the ",
                 "columns present (%s); skipping 'equity %s slope index by %s'."),
          grad, cfg$source, grad,
          paste(intersect(.deprivation_vars, names(tt)), collapse = ", "),
          metric, p$suffix),
          call. = FALSE, immediate. = TRUE)
        next
      }
      out_vars <- p$out_vars                    # always contains "year"
      missing_vars <- setdiff(out_vars, names(tt))
      if (length(missing_vars) > 0L) {
        warning(sprintf(
          paste0("equity tables: stratification variable(s) %s are not in the ",
                 "`%s` summary (available: %s); skipping 'equity %s slope index ",
                 "by %s'. Add them to `strata_for_output` in the sim_design YAML ",
                 "and re-run $export_summaries()."),
          paste(sQuote(missing_vars), collapse = ", "), cfg$source,
          paste(setdiff(names(tt), c(valcols, "popsize")), collapse = ", "),
          metric, p$suffix),
          call. = FALSE, immediate. = TRUE)
        next
      }
      keep <- c("mc", "scenario", out_vars, grad)

      # Collapse every summary variable NOT named in this stratum (agegrp, sex,
      # ...); sum the disease value columns and the deprivation-group
      # population. When `grad` is `qimd` this same aggregation is what
      # collapses the deciles.
      agg <- tt[, c(lapply(.SD, sum), .(popsize = sum(popsize))),
                .SDcols = valcols, by = keep]

      long <- melt(agg, id.vars = c(keep, "popsize"),
                   measure.vars = valcols, variable.name = "disease",
                   value.name = "val")
      long[, disease := as.character(disease)]
      if (!is.na(cfg$strip)) long[, disease := gsub(cfg$strip, "", disease)]

      # Benefit versus comparator, from the baseline year onwards.
      idcols <- setdiff(c(keep, "disease"), "scenario")   # mc, <out_vars>, <grad>, disease
      cmp <- long[scenario == comparator_scenario & year >= baseline_year]
      cmp[, scenario := NULL]
      setnames(cmp, c("val", "popsize"), c("val_cmp", "pop_cmp"))
      d <- long[scenario != comparator_scenario & year >= baseline_year]
      if (nrow(d) == 0L || nrow(cmp) == 0L) next

      if (cfg$dir == "prevented") {
        d[cmp, on = idcols, `:=`(B_year = i.val_cmp - val, N_cmp = i.pop_cmp)]
      } else {
        d[cmp, on = idcols, `:=`(B_year = val - i.val_cmp, N_cmp = i.pop_cmp)]
      }
      # Reference population for the ridit ranking, the population weights of the
      # weighted regression and the per-capita denominator. "comparator" uses
      # the comparator scenario's population so the deprivation ranking is
      # IDENTICAL across the compared scenarios and years (recommended for
      # between-scenario comparison: it prevents intervention-induced shifts in
      # the socioeconomic composition from being read as changes in inequality,
      # per Renard et al. 2019). "scenario" uses each intervention scenario's own
      # population (the composition it actually produces). `popsize` here is the
      # intervention scenario's own deprivation-group population.
      d[, N := if (ridit_reference == "comparator") N_cmp else popsize]
      d <- d[is.finite(B_year) & is.finite(N) & N > 0]
      if (nrow(d) == 0L) next

      # Cumulative benefit over year within
      # (mc, scenario, <out_vars except year>, <grad>, disease).
      by_cum <- c("mc", "scenario", setdiff(out_vars, "year"), grad, "disease")
      setkeyv(d, c(by_cum, "year"))
      d[, B := cumsum(B_year), by = by_cum]

      # All-diseases summed row for CPP/CYPP (plain event-count sum; carries
      # comorbidity multiplicity -- not unique persons). DPP is already
      # all-cause; net QALYs is already a total. The CMS multimorbidity
      # indicators (cms1st_cont, cmsmm*, cmscs*) are derived metrics rather than
      # incident diseases, so they are kept as their own rows but excluded from
      # the sum.
      if (metric %in% c("cpp", "cypp")) {
        by_tot <- c("mc", "scenario", out_vars, grad)
        tot <- d[!grepl("^cms", disease),
                 .(disease = "all_diseases_sum", B = sum(B),
                   B_year = sum(B_year), N = N[1L]), by = by_tot]
        d <- rbind(d, tot, fill = TRUE)
      }

      set(d, j = "rank", value = deprivation_rank(d[[grad]], grad))

      # A deprivation gradient needs at least two groups. Scenario-years
      # covering a single group (e.g. interventions targeted at one decile, or
      # a scenario whose later years only touch one decile) have an undefined
      # slope index and are omitted downstream (the closed-form slope returns
      # NA and is dropped at the finite-value filter). Report the omission so it
      # is transparent rather than silent.
      if (self$design$sim_prm$logs) {
        cov_grp <- c("scenario", out_vars)
        ndec <- unique(d[, c(cov_grp, grad), with = FALSE])[, .(nd = .N),
                                                            by = cov_grp]
        insuff <- ndec[nd < 2L]
        if (nrow(insuff) > 0L) {
          message("  ", metric, " by ", p$suffix, ": omitting ",
                  nrow(insuff), " output group(s) with <2 ", grad,
                  " groups (undefined gradient); scenario(s): ",
                  paste(sort(unique(insuff$scenario)), collapse = ", "))
        }
      }

      by_grp <- c("mc", "scenario", out_vars, "disease")
      idx <- calc_equity_slope_indices(
        d[, .SD, .SDcols = c(by_grp, "rank", "B", "N")], by = by_grp)

      m <- melt(idx, id.vars = by_grp,
                measure.vars = c("AEI_total", "AEI_per100k", "REI_rel", "RII_ratio"),
                variable.name = "type", value.name = "value")
      # Quantile over the finite Monte-Carlo draws only (RII_ratio / REI_rel can
      # be NA for degenerate gradients).
      m <- m[is.finite(value)]
      if (nrow(m) == 0L) next

      setkey(m, "type")
      out <- m[, safe_fquantile_byid(value, prbl, id = as.character(type),
                                     rounding = FALSE),
               keyby = eval(setdiff(by_grp, "mc"))]
      setnames(out, c(setdiff(by_grp, "mc"), "type",
                      scales::percent(prbl, prefix = "equity_")))
      setkeyv(out, c(setdiff(by_grp, "mc"), "type"))
      # The filename suffix reads as "stratified by", but the gradient token in
      # it means "gradient over"; the `gradient` column makes each file
      # self-describing regardless of how it was named.
      out[, gradient := grad]
      setcolorder(out, c("gradient", setdiff(by_grp, "mc"), "type"))

      fwrite(out, file.path(
        tables_dir,
        paste0("equity ", metric, " slope index by ", p$suffix,
               " (not standardised).csv")
      ))
      rm(agg, long, cmp, d, idx, m, out)
    }

    rm(tt)
    gc(verbose = FALSE)
  }

  invisible(NULL)
})
