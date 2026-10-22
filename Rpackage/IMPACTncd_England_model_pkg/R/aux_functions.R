## IMPACTncdEngland is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngland: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
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


`:=` = function(...)
  NULL # due to NSE notes in R CMD check

#' Null coalescing operator
#'
#' @description Returns the left-hand side if it is not NULL, otherwise returns the right-hand side.
#' @param x Left-hand side value.
#' @param y Right-hand side value (default).
#' @return x if not NULL, otherwise y.
#' @examples
#' NULL %||% "default"
#' "value" %||% "default"
#' @export
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.onUnload <- function(libpath) {
  library.dynam.unload("IMPACTncdEngland", libpath)
}



# NOTE implementation below suffers by race conditions although it looks safer.
# fwrite_safe <- function(x,
#                         file,
#                         append = TRUE,
#                         threat_safe = append,
#                         ...) {
#   if (append) {
#     if (file.exists(file)) {
#       col_names_disk <- names(fread(file, nrows = 0))
#       col_names_file <- names(x)
#       col_names <- outersect(col_names_disk, col_names_file)
#       if (length(col_names) > 0)
#         x[, (col_names) := NA]
#       setcolorder(x, col_names_disk)
#     }
#   }
#
#   # create threat-safe mechanism
#   flock <- paste0(file, ".lock")
#
#
#   while (file.exists(flock)) {
#     Sys.sleep(runif(1))
#   }
#
#     file.create(flock)
#     fwrite(x, file, append, ...)
#     on.exit(if (file.exists(flock)) file.remove(flock))
#
# }

# inflate ----
#' Inflate Values Based on Percentage Rate, Year, and Baseline Year
#' @description
#' This function inflates values based on a given percentage rate, year, and baseline year. It calculates the inflated value using the formula: \code{x * (1 + percentage_rate / 100) ^ (year - baseline_year)}.
#'
#' @param x A numeric vector or matrix containing the values to be inflated.
#' @param percentage_rate A numeric value indicating the percentage rate for inflation.
#' @param year The target year for inflation.
#' @param baseline_year The baseline year to calculate the inflation from.
#'
#' @return A numeric vector or matrix containing the inflated values.
#'
#' @export
inflate <- function(x, percentage_rate, year, baseline_year) {
  x * (1 + percentage_rate / 100) ^ (year - baseline_year)
}
# inflate(1000, 3, 2011:2020, 2013)

# deflate ----
#' Deflate Values Based on Percentage Rate, Year, and Baseline Year
#' @description
#' This function deflates values based on a given percentage rate, year, and baseline year. It calculates the deflated value using the formula: \code{x * (1 - percentage_rate / 100) ^ (year - baseline_year)}.
#'
#' @param x A numeric vector or matrix containing the values to be deflated.
#' @param percentage_rate A numeric value indicating the percentage rate for deflation.
#' @param year The target year for deflation.
#' @param baseline_year The baseline year to calculate the deflation from.
#'
#' @return A numeric vector or matrix containing the deflated values.
#'
#' @export
deflate <- function(x, percentage_rate, year, baseline_year) {
  x * (1 - percentage_rate / 100) ^ (year - baseline_year)
}


# plot_cor ----
# Necessary aux functions
# plots mean x by y stratified by z
#' Plot Weighted Correlation Between Two Variables
#' @description
#' This function generates a scatter plot with a weighted correlation line between two variables based on provided data.
#'
#' @param x A character string specifying the variable for the x-axis.
#' @param y A character string specifying the variable for the y-axis.
#' @param z A character string specifying the variable used for faceting the plot.
#' @param wt A character string specifying the variable used for weighting the correlation.
#' @param dt A data table containing the variables specified by \code{x}, \code{y}, \code{z}, and \code{wt}.
#'
#' @return A plot object representing the scatter plot with a weighted correlation line.
#'
#' @export
plot_cor <- function(x, y, z, wt, dt) {
  xlab  <- toupper(x)
  ylab  <- toupper(y)
  title <- paste0(ylab, "~", xlab, "|", toupper(z))

  dt[, .(y = weighted.mean(get(y), get(wt), na.rm = TRUE)),
     keyby = .(x = round(get(x)), z = get(z))][, qplot(
    x,
    y,
    data = na.omit(.SD),
    facets = as.formula("~z"),
    geom = "smooth",
    main = title,
    xlab = xlab,
    ylab = ylab
  )]
}

# distr_best_fit ----
#' helper func for gamlss::fitDistr
#' @description
#' This function fits distribution parameters using the Generalized Additive Model for Location, Scale, and Shape (GAMLSS). It can perform fitting based on either minimum prediction global deviance (if \code{pred} is \code{TRUE}) or Bayesian Information Criterion (BIC) for model selection.
#'
#' @param dt A data table containing the variables specified by \code{var} and \code{wt}.
#' @param var A character string specifying the variable for which distribution parameters are to be fit.
#' @param wt A character string specifying the variable used for weighting the distribution fit.
#' @param distr_family A character string specifying the distribution family to fit (e.g., "GA" for Gamma).
#' @param distr_extra Additional parameters for the distribution family. Defaults to \code{NULL}.
#' @param pred Logical, indicating whether to perform fitting based on minimum prediction global deviance. Default is \code{FALSE}.
#' @param seed An optional seed for reproducibility when \code{pred} is \code{TRUE}.
#' @param trace Logical, indicating whether to print trace information during fitting. Default is \code{TRUE}.
#'
#' @return A fitted distribution object based on GAMLSS.
#'
#' @export
distr_best_fit <-
  function(dt,
           var,
           wt,
           distr_family,
           distr_extra = NULL,
           pred = FALSE,
           seed = NULL,
           trace = TRUE) {
    if (pred) {
      print("Selection based on minimum prediction global deviance")
      if (!is.null(seed))
        set.seed(seed)
      lns <- sample(nrow(dt), round(nrow(dt) * 0.8))
      dt_trn   <- dt[lns,] # train dataset
      dt_crv   <- dt[!lns,]  # cross-validation dataset
      marg_distr <- gamlss::fitDistPred(
        dt_trn[[var]],
        type = distr_family,
        weights = dt_trn[[wt]],
        extra = distr_extra,
        try.gamlss = TRUE,
        trace = trace,
        newdata = dt_crv[[var]]
      )
    } else {
      print("Selection based on BIC")
      marg_distr <-
        gamlss::fitDist(
          dt[[var]],
          log(nrow(dt)),
          type = distr_family,
          weights = dt[[wt]],
          extra = distr_extra,
          try.gamlss = TRUE,
          trace = trace
        )
    }
    marg_distr
  }

# distr_validation ----
#' Distribution Validation Diagnostics
#' @description
#' This function performs distribution validation diagnostics based on a fitted distribution model. It generates validation plots and diagnostic information to assess the goodness of fit.
#'
#' @param marg_distr A fitted distribution object obtained from the \code{distr_best_fit} function.
#' @param dt A data table containing the variables used in the fitting process.
#' @param title A character string specifying the title for the validation plots.
#' @param discrete Logical, indicating whether the variable is discrete. Default is \code{FALSE}.
#' @param smooth A numeric value specifying the smoothing parameter for the validation plots. Default is \code{0.35}.
#'
#' @return A list containing validation plots and diagnostic information.
#'
#' @export
distr_validation <-
  function(marg_distr,
           dt,
           title,
           discrete = FALSE,
           smooth = 0.35) {
    params <- vector("list", length(marg_distr$parameters) + 1L)
    names(params) <- c(marg_distr$parameters, "p")
    for (i in marg_distr$parameters)
      params[[i]] <- get(i, marg_distr)

    params$p <- seq(0.001, 0.999, length.out = nrow(dt))
    distr_nam <- marg_distr$family[[1]]
    y <- do.call(paste0("q", distr_nam), params)
    y <- y[between(y, dt[, min(var)], dt[, max(var)])]
    y_wt <- rep(1 / length(y), length(y))
    # validation plots
    # see http://www.csss.washington.edu/files/working-papers/2002/wp27.pdf
    out <- reldist_diagnostics(
      dt$var,
      y,
      dt[, wt / sum(wt)],
      y_wt,
      main = title,
      discrete = discrete,
      smooth = smooth
    )
    out
  }

# centile_predictAll ----
#' Predict Centiles from a GAMLSS Model for Multiple Observations
#' @description
#' This function predicts centiles from a Generalized Additive Model for Location, Scale, and Shape (GAMLSS) model for multiple observations.
#'
#' @param gamlss_model A GAMLSS model object created using the \code{gamlss} package.
#' @param orig_data The original data used to fit the GAMLSS model.
#' @param newdata The new data for which centiles are to be predicted.
#' @param cent A numeric value specifying the centile to predict. Default is \code{0.5} (median).
#'
#' @return A numeric vector containing the predicted centiles for each observation in \code{newdata}.
#'
#' @export
centile_predictAll <-
  function(gamlss_model, orig_data, newdata, cent = 0.5) {
    stopifnot("gamlss" %in% class(gamlss_model))
    tt <- predictAll(gamlss_model, newdata, data = orig_data)
    tt$p <- cent
    distr_nam <- attr(tt, "family")[[1]]
    out <- do.call(paste0("q", distr_nam), tt)
    return(out)
  }

# mean_predictAll ----
#' Predict Mean Values from a GAMLSS Model for Multiple Observations
#' @description
#' This function predicts mean values from a Generalized Additive Model for Location, Scale, and Shape (GAMLSS) model for multiple observations.
#'
#' @param gamlss_model A GAMLSS model object created using the \code{gamlss} package.
#' @param orig_data The original data used to fit the GAMLSS model.
#' @param newdata The new data for which predictions are to be made.
#'
#' @return A numeric vector containing the predicted mean values for each observation in \code{newdata}.
#'
#' @export
mean_predictAll <- function(gamlss_model, orig_data, newdata) {
  stopifnot("gamlss" %in% class(gamlss_model))
  tt <- predictAll(gamlss_model, newdata, data = orig_data)
  setDT(tt)
  tt[, n := 1e5]
  setcolorder(tt, "n")
  distr_nam <- attr(tt, "family")[[1]]
  out <- apply(tt, 1, function(x) {
    mean(do.call(paste0("r", distr_nam), as.list(x)))
  })
  return(out)
}

# plot_synthpop_val ----
#' Plot Synthpop Validation Density and Cumulative Density Functions
#' @description
#' This function generates and saves density and cumulative density function plots for synthpop validation based on the provided data.
#'
#' @param dt A data table containing the variables specified by \code{x}, \code{grp}, and \code{wt}.
#' @param x A variable to be plotted on the x-axis.
#' @param grp A grouping variable for faceting the plots.
#' @param wt A variable used for weighting the plots.
#' @param title A character string specifying the title for the plots.
#' @param x_label A character string specifying the label for the x-axis.
#' @param standardised_to_grp A logical vector indicating whether to standardize weights within each group. Default is \code{c(FALSE, TRUE)}.
#' @param print_to_screen A logical vector indicating whether to print the plots to the screen. Default is \code{c(FALSE, TRUE)}.
#'
#' @return None (saves density and cumulative density function plots to specified directory).
#'
#' @export
plot_synthpop_val <-
  function(dt,
           x,
           grp,
           wt,
           title,
           x_label,
           standardised_to_grp = c(FALSE, TRUE),
           print_to_screen = c(FALSE, TRUE)) {
    if (standardised_to_grp) {
      dt[, weight := get(wt) / sum(get(wt)), by = c("type", grp)]
    } else {
      dt[, weight := get(wt) / sum(get(wt)), by = c("type")]
    }

    if ("year" %in% grp)
      dt[, year := year + 2000L]
    x <- enquo(x)
    p <- ggplot(dt, aes(
      !!x,
      colour = type,
      weight = weight,
      linetype = type
    )) +
      geom_density() +
      facet_grid(grp) +
      xlab(x_label) +
      ggtitle(title)
    if (print_to_screen)
      print(p)
    suppressWarnings(
      cowplot::ggsave2(
        paste0(gsub(" ", "_", title), "_density.png"),
        p,
        width = 16,
        height = 9,
        units = "cm",
        scale = 2,
        dpi = 300,
        path = "./validation/synthpop_models"
      )
    )

    p <- ggplot(dt, aes(
      !!x,
      colour = type,
      weight = weight,
      linetype = type
    )) +
      stat_ecdf() +
      facet_grid(grp) +
      xlab(x_label) +
      ggtitle(title)
    if (print_to_screen)
      print(p)

    if ("year" %in% grp)
      dt[, year := year - 2000L]

    suppressWarnings(
      cowplot::ggsave2(
        paste0(gsub(" ", "_", title), "_cdf.png"),
        p,
        width = 16,
        height = 9,
        units = "cm",
        scale = 2,
        dpi = 300,
        path = "./validation/synthpop_models"
      )
    )
  }





# get_ons_mrtl ----
# get observed mortality from ONS data
# mortality_by_agegrp20_for_validation.fst; mortality_by_agegrp5_for_validation.fst;
#' Get ONS (Office for National Statistics) Mortality Data
#' @description
#' This function retrieves mortality data from the Office for National Statistics (ONS) for a specified disease, type, and age group width.
#'
#' @param disease A character string specifying the disease for which mortality data is requested.
#' @param type A character vector specifying the type of mortality data to retrieve. Options are "rate" for mortality rates and "absolute" for absolute mortality counts. Default is "rate".
#' @param agegrp_width An integer vector specifying the width of age groups in years. Options are 5 or 20. Default is 5.
#'
#' @return A data table containing mortality data for the specified disease, type, and age group width.
#'
#' @export
get_ons_mrtl <-
  function(disease,
           type = c("rate", "absolute"),
           agegrp_width = c(5L, 20L)) {
    if (type == "rate")
      prefix <- "Mx_"
    if (type == "absolute")
      prefix <- "deaths_"

    colnams <-
      c(
        "year",
        "sex",
        "qimd",
        paste0("agegrp", agegrp_width),
        paste0(prefix, disease),
        "pops"
      )
    if (agegrp_width == 5)
      tt <-
      read_fst(
        "./ONS_data/mortality_by_agegrp5_for_validation.fst",
        columns = colnams,
        as.data.table = TRUE
      )
    if (agegrp_width == 20)
      tt <-
      read_fst(
        "./ONS_data/mortality_by_agegrp20_for_validation.fst",
        columns = colnams,
        as.data.table = TRUE
      )
    tt
  }
# get_ons_mrtl("chd", "rate", 20)

# get_ons_incd ----
# get observed cancer incidence from ONS data
# cancer_incd_by_agegrp5_for_validation.fst; cancer_incd_by_agegrp5_for_validation.fst;
#' Get ONS (Office for National Statistics) Cancer Incidence Data
#' @description
#' This function retrieves cancer incidence data from the Office for National Statistics (ONS) for a specified disease, type, and age group width.
#'
#' @param disease A character string specifying the cancer type for which incidence data is requested.
#' @param type A character vector specifying the type of incidence data to retrieve. Options are "rate" for incidence rates and "absolute" for absolute incidence counts.
#' @param agegrp_width An integer specifying the width of age groups in years. Options are 5 or 20. Default is 20.
#'
#' @return A data table containing cancer incidence data for the specified disease, type, and age group width.
#'
#' @export
get_ons_incd <-
  function(disease,
           type = c("rate", "absolute"),
           agegrp_width = 20L) {
    if (type == "rate")
      prefix <- "rate_"
    if (type == "absolute")
      prefix <- "cases_"

    colnams <-
      c(
        "year",
        "sex",
        "qimd",
        paste0("agegrp", agegrp_width),
        paste0(prefix, disease),
        "pops"
      )
    if (agegrp_width == 5)
      tt <-
      read_fst(
        "./ONS_data/cancer_incd_by_agegrp5_for_validation.fst",
        columns = colnams,
        as.data.table = TRUE
      )
    if (agegrp_width == 20)
      tt <-
      read_fst(
        "./ONS_data/cancer_incd_by_agegrp20_for_validation.fst",
        columns = colnams,
        as.data.table = TRUE
      )
    tt
  }


# generate_rns ----
#' Generate Random Numbers for Monte Carlo iteration
#' @description
#' This function generates random numbers for a Monte Carlo iteration using the dqRNG package. It sets the seed based on the simulation run and updates specified columns in the input data table with random numbers.
#'
#' @param mc A numeric value representing the Monte Carlo iteration run.
#' @param dt A data table containing the columns to be updated with random numbers.
#' @param colnams A character vector specifying the names of columns to be updated with random numbers.
#'
#' @return The input data table with specified columns updated with random numbers.
#'
#' @export
generate_rns <- function(mc, dt, colnams) {
  dqRNGkind("pcg64") # dqRNGkind("Xoroshiro128+") ~10% faster
  SEED <- 4719349L # sample(1e7, 1)
  set.seed(SEED + mc)
  dqset.seed(SEED, mc)

  nrows <- nrow(dt)
  for (nam in colnams)
    set(dt, NULL, nam, dqrunif(nrows))
  invisible(dt)
}



# output_dir ----
#' Generate Output Directory Path
#' @description
#' This function generates the path to the output directory by combining the simulation parameter's output directory
#' and the provided suffix (if any).
#'
#' @param x A character string representing the suffix to be appended to the output directory path (default: "").
#'
#' @details
#' The function uses the `file.path` function to combine the simulation parameter's output directory with the provided
#' suffix (if any). It is typically used to organize output files in a structured manner.
#'
#' @return A character string representing the path to the output directory.
#'
#' @export
output_dir <- function(x = "") {
  file.path(design$sim_prm$output_dir, x)
}


# weighted_quantile ----
#' Weighted empirical quantiles
#'
#' @description
#' Quantiles of `x` under the sampling weights `w`, i.e. quantiles of the
#' population `x` represents rather than of the simulated sample. Used by
#' `Simulation$export_xps()` to describe exposure *distributions* alongside
#' their weighted means, so a scenario that shifts a whole distribution can be
#' told apart from one that merely moves its mean.
#'
#' @param x Numeric vector of values.
#' @param w Numeric vector of non-negative weights, the same length as `x`.
#' @param probs Numeric vector of probabilities, each strictly between 0 and 1.
#'
#' @details
#' This is the **inverse-CDF (lower) quantile**: the smallest observed `x`
#' whose cumulative weight reaches `p`. No interpolation is done, deliberately
#' -- several model exposures are discrete (`active_days` is 0-7, smoking
#' status is a small set of levels, `smok_packyrs` is integer-valued), and
#' interpolating between adjacent observations would publish values the
#' variable cannot take. For a continuous exposure with many simulants per
#' stratum the difference from an interpolating definition is far below the
#' Monte Carlo noise.
#'
#' Elements with a missing value, a missing weight or a non-positive weight are
#' dropped; if nothing is left, `NA_real_` is returned for every probability.
#'
#' **Not additive.** Unlike a weighted mean or a count, quantiles cannot be
#' combined: the median of two strata is not a function of their medians, and
#' averaging per-iteration quantiles across Monte Carlo iterations does not give
#' the quantile of the pooled draw. Summarise them across iterations the way the
#' table exports do (percentiles of the per-iteration values), and read each
#' stratum on its own.
#'
#' @return A numeric vector the same length as `probs`.
#'
#' @examples
#' # With equal weights this agrees with the type-1 (inverse CDF) definition
#' x <- c(1, 2, 3, 4, 5)
#' weighted_quantile(x, rep(1, 5), c(0.25, 0.5, 0.75))
#' stats::quantile(x, c(0.25, 0.5, 0.75), type = 1, names = FALSE)
#'
#' # Weighting shifts the median toward the heavily-weighted values
#' weighted_quantile(x, c(1, 1, 1, 50, 50), 0.5)
#'
#' @export
weighted_quantile <- function(x, w, probs) {
  wsel <- !is.na(w) & w > 0
  if (all(wsel)) {
    wsel <- NULL
    W <- sum(w)
  } else {
    W <- sum(w[wsel])
  }
  .weighted_quantile_pre(x, w, wsel, W, probs)
}


# .weighted_quantile_pre ----
# weighted_quantile() with the weight cleaning already done by the caller.
#
# The single implementation of the quantile itself; weighted_quantile() is a
# thin wrapper that computes `wsel`/`W` and delegates here, so there is no
# second copy of the definition to drift.
#
# It exists because export_xps() calls this once per exposure per group -- ~30
# times over the SAME weight vector. Recomputing `!is.na(w) & w > 0`,
# re-subsetting `w` and re-summing it 30 times was ~20% of the quantile cost.
# The caller hoists that out of its per-exposure loop and passes the result:
#
#   wsel  NULL when every weight is usable (the common case, no subsetting at
#         all), otherwise the logical mask.
#   W     sum of the usable weights. Recomputed here only if `x` has NAs, since
#         then the usable rows are a subset of the ones `W` was summed over.
#
# @noRd
.weighted_quantile_pre <- function(x, w, wsel, W, probs) {
  if (!is.null(wsel)) {
    x <- x[wsel]
    w <- w[wsel]
  }
  if (anyNA(x)) {
    ok <- !is.na(x)
    x <- x[ok]
    w <- w[ok]
    W <- sum(w)
  }
  if (!length(x)) return(rep(NA_real_, length(probs)))
  x <- as.double(x)
  # method = "radix" explicitly: order()'s "auto" picks radix for doubles
  # anyway, but saying so avoids the dispatch on every one of these calls.
  o <- order(x, method = "radix")
  xs <- x[o]
  cw <- cumsum(as.double(w[o]))
  # findInterval(..., left.open = TRUE) counts the cumulative weights strictly
  # below p, so +1 lands on the first one that reaches it. pmin() guards the
  # top end against floating-point drift in cw[.N] / sum(w).
  xs[pmin(
    findInterval(probs * W, cw, left.open = TRUE) + 1L,
    length(xs)
  )]
}

# with_original_dimd ----
# Evaluate `expr` with `pop$dimd` temporarily reset to the pre-scenario
# deprivation snapshotted in `pop$dimd_orig`, restoring the scenario's `dimd`
# on the way out (including on error, via on.exit).
#
# WHY THIS EXISTS
# ---------------
# `Disease$gen_parf()` recomputes nothing from the scenario's exposures. It
# reads the cached PARF table from `simulation/parf/` and joins it onto
# `sp$pop` with `on = .NATURAL` -- and `dimd` is one of those natural keys. The
# payload is
#     p0 = mu       * (1 - parf)          background ("unexposed") incidence
#     m0 = mu_mrtl  * (1 - parf_mrtl)     background disease mortality
# where `mu` is the OBSERVED, decile-specific rate. So p0/m0 carry the part of
# the deprivation gradient the modelled risk factors do not explain, and which
# decile's version a simulant receives is decided by the value of `dimd` at the
# moment gen_parf() runs.
#
# gen_parf() runs AFTER the primary-prevention scenario (so that
# `set_init_prvl()` sees the scenario world -- see run_sim()). A scenario that
# rewrites `dimd` to model deprivation change would therefore hand every
# simulant the target decile's background risk as an unavoidable side effect of
# the call order, moving the non-modelled gradient whether or not that was
# intended. The measured size of that side effect is large: for COPD the
# most-deprived decile's p0 is ~2.6x the least-deprived's, and the residual
# all-cause mortality m0 ~1.4x.
#
# Which of the two is wanted is a modelling choice, so it is stated by the
# scenario (`update_primary_prevention_scn(fn, parf_deprivation = ...)`) and
# applied here, rather than inferred from call order.
#
# The `dimd` axis is exactly the non-modelled channel: every exposure
# distribution table and the whole smoking sub-model are keyed on `qimd`
# instead, so a scenario that rewrites both keeps its modelled effect in full
# under either setting.
#
# @param pop A data.table (`sp$pop`), modified by reference.
# @param parf_deprivation `"original"` to swap, anything else to do nothing.
# @param expr Expression evaluated with the chosen `dimd` in place.
# @return `TRUE` if a swap happened, `FALSE` otherwise, invisibly.
# @noRd
with_original_dimd <- function(pop, parf_deprivation, expr) {
  swapped <- identical(parf_deprivation, "original") &&
    all(c("dimd", "dimd_orig") %chin% names(pop)) &&
    !identical(pop[["dimd"]], pop[["dimd_orig"]])

  if (swapped) {
    # `dimd_scn` is deliberately a third column rather than a saved vector:
    # `expr` may add or drop rows' worth of work but never reorders sp$pop, and
    # keeping it in the table means the restore is a plain column swap that
    # cannot silently mis-align.
    data.table::set(pop, NULL, "dimd_scn", pop[["dimd"]])
    data.table::set(pop, NULL, "dimd", pop[["dimd_orig"]])
    on.exit(
      {
        data.table::set(pop, NULL, "dimd", pop[["dimd_scn"]])
        data.table::set(pop, NULL, "dimd_scn", NULL)
      },
      add = TRUE
    )
  }

  force(expr)
  invisible(swapped)
}

# select_wrswor_by ----
# Weighted random sampling without replacement within groups, driven by one
# caller-supplied exponential variate per ROW, so that a group's selection
# depends on that group's rows alone (common-random-number safe).
#
# Within every group of `dt` defined by the columns in `by`, selects `size`
# rows with probability proportional to `weight`, without replacement, and
# returns their row numbers. Each row's key is `rexp / weight`; within each
# group the `size` rows with the smallest keys are kept. This is Algorithm A of
# Efraimidis & Spirakis (2006) with log-transformed keys -- i.e.
# `wrswoR::sample_int_rank()` vectorised over groups -- so the selected set has
# exactly the law of successive sampling: the same as
# `sample.int(n, size, prob = weight)` and as `wrswoR::sample_int_expj()`.
#
# WHY THIS EXISTS (GitHub issue #8)
# ---------------------------------
# `Disease$set_init_prvl()` used to call `wrswoR::sample_int_expj()` once per
# stratum with `by =`. All strata then shared ONE base-R stream, visited in
# order of first appearance, and `sample_int_expj()` consumes a data-dependent
# number of variates (its exponential jumps depend on the weights and on
# `size`). Any change to one stratum -- a weight, a row moving in or out, a
# different `size` -- shifted the stream for every stratum visited after it.
# Synthpop pids are ordered by age in the initial year, so the strata of future
# entrants come before every `init_year` stratum, and a scenario that touched
# only entrants from `sc_year` on re-drew WHICH simulants were prevalent in
# `init_year`: the baseline/scenario pairing (common random numbers) was lost
# before the scenario did anything.
#
# Here every row owns its variate. The contract the caller must keep is that
# `rexp` is drawn at a stream position that does not depend on the data (e.g.
# `dqrexp(nrow(dt))` straight after a `dqset.seed()`), and that the same
# simulant-year sits at the same row of `dt` in every scenario. Perturbing any
# other group, including moving rows into or out of it, then cannot change a
# group's selection.
#
# Semantics: reproduces `sample_int_expj()` wherever its result is well
# defined (positive, finite weights), and fails loudly wherever it errors or
# silently misbehaves. `sample_int_expj()` never inspects the weights: with NA,
# negative or Inf weights it always selects that row, and it returns
# zero-weight rows, without a warning.
# * Only rows with `size > 0` are candidates; the `weight`/`rexp` of other rows
#   are ignored (not validated).
# * `size` must be non-missing, finite, whole and non-negative on every row. It
#   must be constant among the rows of a group that have `size > 0` (rows with
#   `size == 0` are never candidates) and no larger than their number.
# * Candidate weights must be finite, non-negative and not NA. A zero weight is
#   never selected, as in `sample.int()`; needing one is an error.
# * NA in a `by` column forms a group of its own, as in data.table's `by =`.
# * Ties between keys (probability zero) are broken by row order; keys are
#   sorted on all 64 bits whatever the user's numeric rounding setting.
#
# @param dt A data.table. Not modified.
# @param by Character vector of grouping (stratum) column names.
# @param size Name of the column holding the number of rows to select from
#   each group (constant within a group).
# @param weight Name of the column holding the selection weights.
# @param rexp Numeric vector of `nrow(dt)` independent standard exponential
#   variates, one per row of `dt`, drawn at a data-independent stream position.
# @return Integer vector of the selected row numbers of `dt`, ascending;
#   exactly `size` rows per group.
# @references Efraimidis, P.S. & Spirakis, P.G. (2006). Weighted random
#   sampling with a reservoir. Information Processing Letters 97(5), 181-185.
#   doi:10.1016/j.ipl.2005.11.003
# @noRd
select_wrswor_by <- function(dt, by, size, weight, rexp) {
  # ---- argument checks ----
  if (!data.table::is.data.table(dt)) {
    stop("'dt' must be a data.table.", call. = FALSE)
  }
  if (!is.character(by) || length(by) == 0L || anyNA(by)) {
    stop("'by' must be a non-empty character vector.", call. = FALSE)
  }
  for (.arg in list(size, weight)) {
    if (!is.character(.arg) || length(.arg) != 1L || is.na(.arg)) {
      stop(
        "'size' and 'weight' must each be a single column name.",
        call. = FALSE
      )
    }
  }
  missing_cols <- setdiff(c(by, size, weight), names(dt))
  if (length(missing_cols) > 0L) {
    stop(
      "Column(s) not found in 'dt': ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  # Internal column names of the working table; they must not shadow `by`.
  internal <- c(".wrs_size", ".wrs_pos", ".wrs_key", ".wrs_row")
  if (any(internal %chin% by)) {
    stop(
      "'by' must not contain the reserved names ",
      paste(internal, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  if (!is.numeric(rexp) || length(rexp) != nrow(dt)) {
    stop(
      "'rexp' must be a numeric vector with one variate per row of 'dt' ",
      "(length ", nrow(dt), ", got ", length(rexp), ").",
      call. = FALSE
    )
  }

  sz <- dt[[size]]
  if (!is.numeric(sz)) {
    stop("Column '", size, "' must be numeric.", call. = FALSE)
  }
  if (anyNA(sz)) {
    stop("'size' must not be NA (column '", size, "').", call. = FALSE)
  }
  if (any(sz < 0)) {
    stop("'size' must be non-negative (column '", size, "').", call. = FALSE)
  }
  if (!is.integer(sz) && (any(!is.finite(sz)) || any(sz != trunc(sz)))) {
    stop(
      "'size' must be finite whole numbers (column '", size, "').",
      call. = FALSE
    )
  }

  # ---- candidate rows ----
  cand <- which(sz > 0)
  if (length(cand) == 0L) {
    return(integer(0L))
  }

  w <- dt[[weight]][cand]
  e <- rexp[cand]
  if (!is.numeric(w) || anyNA(w) || any(!is.finite(w)) || any(w < 0)) {
    stop(
      "'weight' (column '", weight, "') must be finite, non-negative and ",
      "not NA in every group with size > 0.",
      call. = FALSE
    )
  }
  if (anyNA(e) || any(!is.finite(e)) || any(e < 0)) {
    stop(
      "'rexp' must be finite, non-negative and not NA on every row of a ",
      "group with size > 0.",
      call. = FALSE
    )
  }

  # Only the stratum columns of the candidate rows are copied, never `dt`.
  wk <- dt[cand, by, with = FALSE]
  # Efraimidis-Spirakis key: smaller is selected first. A zero weight gets an
  # infinite key explicitly (0/0 would be NaN when the variate is 0).
  keys <- e / w
  keys[w == 0] <- Inf
  data.table::set(wk, j = ".wrs_key", value = keys)
  data.table::set(wk, j = ".wrs_size", value = as.integer(sz[cand]))
  data.table::set(wk, j = ".wrs_pos", value = w > 0)
  data.table::set(wk, j = ".wrs_row", value = cand)

  # ---- sort: groups contiguous, keys ascending within each group ----
  # Exact 64-bit ordering of the keys (data.table rounds the last bytes of
  # doubles when sorting if the user has changed this global setting).
  old_rounding <- data.table::getNumericRounding()
  if (old_rounding != 0L) {
    data.table::setNumericRounding(0L)
    on.exit(data.table::setNumericRounding(old_rounding), add = TRUE)
  }
  # na.last = TRUE is explicit: setorderv() defaults to na.last = FALSE, which
  # would put NA keys FIRST, i.e. select them. No key can be NA here (weights
  # and variates are validated above); NA stratum values only need to stay
  # together, which either setting guarantees.
  data.table::setorderv(wk, c(by, ".wrs_key"), na.last = TRUE)

  # ---- per-group checks and selection (plain vectors, no NSE) ----
  g <- data.table::rleidv(wk, cols = by) # contiguous group id, 1..G
  n_g <- tabulate(g)
  s <- wk[[".wrs_size"]]
  start_g <- cumsum(n_g) - n_g + 1L # first sorted row of each group
  s_g <- s[start_g]
  if (any(s != s_g[g])) {
    stop("'size' must be constant within each group of 'by'.", call. = FALSE)
  }
  if (any(s_g > n_g)) {
    stop("cannot take a sample larger than the population", call. = FALSE)
  }
  if (any(s_g > tabulate(g[wk[[".wrs_pos"]]], nbins = length(n_g)))) {
    stop("too few positive weights", call. = FALSE)
  }
  # Within-group rank of the key; keep the `size` smallest.
  rk <- seq_along(g) - start_g[g] + 1L
  sel <- wk[[".wrs_row"]][rk <= s]
  sort.int(sel, method = "radix")
}
