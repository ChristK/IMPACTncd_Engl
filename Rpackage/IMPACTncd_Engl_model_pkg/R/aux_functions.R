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


`:=` = function(...)
  NULL # due to NSE notes in R CMD check

.onUnload <- function(libpath) {
  library.dynam.unload("IMPACTncdEngl", libpath)
}

# fwrite_safe ----
#' Ensures that when fwrite appends file colnames of file to be written, match
#' those already in the file
#' @description This function writes a data table to a file with safe appending. It checks if the file already exists and, if appending is enabled, aligns column names and appends the data while preserving existing columns.
#'
#' @param x A data table to be written to the file.
#' @param file A character string specifying the file path to write the data table to.
#' @param append Logical, indicating whether to append to an existing file if it exists. Default is \code{TRUE}.
#' @param ... Additional arguments passed to the \code{\link[data.table]{fwrite}} function.
#'
#' @return None (invisibly returns the result of \code{\link[data.table]{fwrite}}).
#'
#' @export
fwrite_safe <- function(x,
                        file = "",
                        append = TRUE,
                        ...) {
  if (append) {
    if (file.exists(file)) {
      col_names_disk <- names(fread(file, nrows = 0))
      col_names_file <- names(x)
      col_names <- outersect(col_names_disk, col_names_file)
      if (length(col_names) > 0)
        x[, (col_names) := NA]
      setcolorder(x, col_names_disk)
    }
  }
  fwrite(x, file, append, ...)
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

# validate_gamlss_tbl ----
#' Validate GAMLSS Model against Observed Data in a Data Table
#' @description
#' This function validates a Generalized Additive Model for Location, Scale, and Shape (GAMLSS) model against observed data in a data table. It generates simulated data from the GAMLSS model and combines it with the observed data for comparison.
#'
#' @param dt A data table containing the observed data.
#' @param gamlss_tbl A data table containing the GAMLSS model parameters.
#' @param mc An integer specifying the number of Monte Carlo iteration for data validation. Default is 10.
#' @param colname A character string specifying the column name for the distributional parameter in the output.
#' @param distr_nam A character string specifying the distribution name to be used in the validation.
#'
#' @return A data table containing both observed and simulated data for model validation.
#'
#' @export
validate_gamlss_tbl <-
  function(dt,
           gamlss_tbl,
           mc = 10L,
           colname,
           distr_nam = distr_nam) {
    stopifnot(is.data.table(dt), is.data.table(gamlss_tbl), mc >= 1)
    nam_var <- intersect(names(dt), names(gamlss_tbl))
    nam_param <- setdiff(names(gamlss_tbl), nam_var)
    dt[, age := as.integer(age)]
    x <- copy(dt)
    x[, `:=`(type, "Observed")]
    z <- copy(dt)
    z[gamlss_tbl, (nam_param) := mget(nam_param), on = nam_var]
    if (z[is.na(mu), .N] > 0)
      stop("NAs produced in join")
    z[, `:=`(type, "Modelled")]
    z <- rbindlist(rep(list(z), mc))
    z[, p := dqrunif(.N, 0, 0.999)]
    z[, (colname) := do.call((distr_nam), .SD), .SDcols = c("p", nam_param)]
    z[, c("p", (nam_param)) := NULL]
    out <- rbind(x, z, use.names = TRUE, fill = TRUE)
  }

# dt <- dt[age >= 20L]
# gamlss_tbl <- copy(smok_incid_model_tbl)
# colname <- "smok_incid"
# z[, .SD, .SDcols = c("p", nam_param)]
# z[, do.call("qBI", .SD), .SDcols = c("p", nam_param)]

# shift_bypid ----
#' Shift Values by PID (Participant ID) for Integer, Logical, or Numeric Vectors
#' @description
#' This function shifts values by participant ID (PID) for integer, logical, or numeric vectors. It can handle different data types and supports shifting by a specified lag.
#'
#' @param x A vector of integers, logicals, or numerics to be shifted.
#' @param lag An integer specifying the lag for shifting values. A positive lag indicates shifting forward, and a negative lag indicates shifting backward.
#' @param id A vector representing participant IDs associated with each observation in \code{x}.
#' @param replace A value to replace shifted values beyond vector boundaries. Default is \code{NA}.
#'
#' @return A vector with values shifted by the specified lag for each participant ID.
#'
#' @export
shift_bypid <-
  function(x, lag, id, replace = NA) {
    if (lag == 0L) return(x)
    if (typeof(x) == "integer") {
      return(shift_bypidInt(x, lag, replace, id))
    } else if (typeof(x) == "logical") {
      return(shift_bypidBool(x, lag, replace, id))
    } else if (typeof(x) == "double") {
      return(shift_bypidNum(x, lag, replace, id))
    } else
      stop("type of x not supported")
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

# generate_corr_unifs ----
#' Given a correlation matrix (Pearson), produces a matrix of correlated uniforms
#' @description
#' This function generates correlated uniform random variables based on a given correlation matrix using the dqrnorm and pnorm functions.
#'
#' @param n The number of samples to generate.
#' @param M The correlation matrix.
#'
#' @return A matrix of correlated uniform random variables.
#'
#' @export
generate_corr_unifs <- function(n, M) {
  # generate normals, check correlations
  # from http://comisef.wikidot.com/tutorial:correlateduniformvariates
  stopifnot(is.matrix(M))
  # Check that matrix is semi-positive definite
  # NOTE next line crashes frequently!! see
  # https://stat.ethz.ch/pipermail/r-help/2006-March/102703.html for a
  # workaround and https://stat.ethz.ch/pipermail/r-help/2006-March/102647.html
  # for some explanation
  # stopifnot(min(eigen(M, only.values = TRUE)$values) >= 0)

  M_original <- M


  # adjust correlations for uniforms
  for (i in seq_len(dim(M)[[1L]])) {
    for (j in seq_len(dim(M)[[2L]])) {
      if (i != j) {
        M[i, j] <- 2 * sin(pi * M[i, j] / 6)
        M[j, i] <- 2 * sin(pi * M[j, i] / 6)
      }
    }
  }

  X <- matrix(dqrnorm(n * dim(M)[[2]]), n)
  colnames(X) <- colnames(M)

  # induce correlation, check correlations
  Y <- pnorm(X %*% chol(M))

  # message(paste0("Mean square error is: ", signif(sum((cor(Y) - M_original) ^
  # 2), 3)))
  return(Y)
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

#' Function for Linux
#' @return A named numeric vector with five non-negative elements `1min`,
#'   `5min`, and `15min` average CPU utilisation, used RAM, and available RAM.
#' The first values represent estimates of the CPU load during the last
#' minute, the last five minutes, and the last fifteen minutes \[1\]. An idle
#' system have values close to zero, and a heavily loaded system have values
#' near one`.
#'
#' @details
#' This function works only Unix-like system with \file{/proc/loadavg}. It is
#' heavily based on parallely::cpuLoad
#' (\url{https://github.com/HenrikBengtsson/parallelly})
#'
#' @references
#' 1. Linux Load Averages: Solving the Mystery,
#'    Brendan Gregg's Blog, 2017-08-08,
#'    \url{http://www.brendangregg.com/blog/2017-08-08/linux-load-averages.html}
#'
#' @keywords internal
#' @export
sysLoad <- function() {
  if (file.exists("/proc/loadavg")) {
    res <- readLines("/proc/loadavg", n = 1L)
    res <- strsplit(res, split=" ", fixed = TRUE)[[1]]
    res <- as.numeric(res[1:3])
    res <- signif(100 * res/parallel::detectCores(), 2)
    x <- system2('free', args = '-m', stdout = TRUE) # only for linux
    x <- strsplit(x[2], " +")[[1]][3:4]
    x <- round(as.numeric(x)/1024, 2)
    res <- as.numeric(c(res, x))
  } else {
    res <- rep(NA_real_, times = 5L)
  }
  names(res) <- c("1minAvgCPU(%)", "5minAvgCPU(%)", "15minAvgCPU(%)", "UsedRAM(Gb)", "FreeRAM(Gb)")
  res
}

