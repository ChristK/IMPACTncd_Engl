% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_functions.R
\name{fwrite_safe}
\alias{fwrite_safe}
\title{Ensures that when fwrite appends file colnames of file to be written, match
those already in the file}
\usage{
fwrite_safe(x, file = "", append = TRUE, ...)
}
\arguments{
\item{x}{A data table to be written to the file.}

\item{file}{A character string specifying the file path to write the data table to.}

\item{append}{Logical, indicating whether to append to an existing file if it exists. Default is \code{TRUE}.}

\item{...}{Additional arguments passed to the \code{\link[data.table]{fwrite}} function.}
}
\value{
None (invisibly returns the result of \code{\link[data.table]{fwrite}}).
}
\description{
This function writes a data table to a file with safe appending. It checks if the file already exists and, if appending is enabled, aligns column names and appends the data while preserving existing columns.
}
