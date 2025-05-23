% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lookup_dt.R
\name{starts_from_1}
\alias{starts_from_1}
\title{Adjust Integer or Factor Columns to Start from 1}
\usage{
starts_from_1(tbl, on, i, min_lookup, cardinality)
}
\arguments{
\item{tbl}{A data frame containing the columns to be adjusted.}

\item{on}{A character vector specifying the column names in the data frame to be adjusted.}

\item{i}{An integer specifying the index of the column to be adjusted within the 'on' vector.}

\item{min_lookup}{A list containing the minimum values for adjustment corresponding to each column in 'on'.}

\item{cardinality}{A list containing the cardinality of each column in 'on'.}
}
\value{
A vector representing the adjusted column starting from 1.
}
\description{
This function adjusts integer or factor columns in a data frame to start from 1.
It is particularly useful when working with tables and joining operations that
involve integer or factor keys starting from a value other than 1.
}
