% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_functions.R
\name{distr_best_fit}
\alias{distr_best_fit}
\title{helper func for gamlss::fitDistr}
\usage{
distr_best_fit(
  dt,
  var,
  wt,
  distr_family,
  distr_extra = NULL,
  pred = FALSE,
  seed = NULL,
  trace = TRUE
)
}
\arguments{
\item{dt}{A data table containing the variables specified by \code{var} and \code{wt}.}

\item{var}{A character string specifying the variable for which distribution parameters are to be fit.}

\item{wt}{A character string specifying the variable used for weighting the distribution fit.}

\item{distr_family}{A character string specifying the distribution family to fit (e.g., "GA" for Gamma).}

\item{distr_extra}{Additional parameters for the distribution family. Defaults to \code{NULL}.}

\item{pred}{Logical, indicating whether to perform fitting based on minimum prediction global deviance. Default is \code{FALSE}.}

\item{seed}{An optional seed for reproducibility when \code{pred} is \code{TRUE}.}

\item{trace}{Logical, indicating whether to print trace information during fitting. Default is \code{TRUE}.}
}
\value{
A fitted distribution object based on GAMLSS.
}
\description{
This function fits distribution parameters using the Generalized Additive Model for Location, Scale, and Shape (GAMLSS). It can perform fitting based on either minimum prediction global deviance (if \code{pred} is \code{TRUE}) or Bayesian Information Criterion (BIC) for model selection.
}
