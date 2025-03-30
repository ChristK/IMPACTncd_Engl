## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
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

# file.remove(list.files("./output/", full.names = TRUE, recursive = TRUE))
# file.remove(list.files("./Rpackage/IMPACTncd_Engl_model_pkg/src", full.names = TRUE, recursive = TRUE, pattern = "\\.o$|\\.so$"))

# If segfault from C stack overflow see
# https://github.com/Rdatatable/data.table/issues/1967

#' @section CRAN Mirror Selection:
#' This block sets the CRAN repository mirror only if:
#'   - The code is **not running inside a Docker container**, and
#'   - The CRAN repository option is either unset or set to the default placeholder ("@CRAN@").
#'
#' This prevents unnecessary or problematic CRAN mirror selection inside Docker,
#' where the environment is usually pre-configured or isolated from user input.
#'
#' The check `file.exists("/.dockerenv")` is a common way to detect Docker containers.

if (!file.exists("/.dockerenv")) {
  repos <- getOption("repos")
  if (is.null(repos) || repos["CRAN"] == "@CRAN@") {
    chooseCRANmirror(ind = 1)
  }
}

cat("Initialising IMPACTncd_Engl model...\n\n")

# Ensure 'pak' is installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Ensure 'CKutils' is installed from GitHub if missing
if (!requireNamespace("CKutils", quietly = TRUE)) {
  remotes::install_github("ChristK/CKutils", upgrade = "never", force = TRUE)
}
library(CKutils)

# Set development mode flag
dev_mode <- TRUE  # Set to FALSE for production

# Environment-specific options
options(rgl.useNULL = TRUE) # suppress error by demography in rstudio server
if (dev_mode) {
  options(future.fork.enable = TRUE) # enable for development only
  options(future.globals.maxSize = +Inf)
  options(future.rng.onMisuse = "ignore") # Remove false warning
}
options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

# install missing packages
pkg_list <- readLines("docker_setup/r-packages.txt", warn = FALSE)
pkg_list <- trimws(pkg_list)
pkg_list <- pkg_list[nzchar(pkg_list) & !grepl("^#", pkg_list)]
dependencies(pkg_list)
rm(pkg_list)


installLocalPackageIfChanged(pkg_path = "./Rpackage/IMPACTncd_Engl_model_pkg/",
                            snapshot_path = "./Rpackage/.IMPACTncd_Engl_model_pkg_snapshot.rds")
library(IMPACTncdEngl)

####################################################################################
## Fix: Added this to update .Random.seed / seed only for Windows.
## Run these two lines of code shown below for the first time and then we clear the work
## space using `rm(list = ls(all = TRUE))` after which we run these two lines again we get the error mentioned
## so to reproduce it remove this line of code `runif(1)`
# source("global.R")
# IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")
## Error generated otherwise (without the code below): Error in get(".Random.seed", .GlobalEnv) :
## object '.Random.seed' not found
####################################################################################
invisible(runif(1))

# roxygen2::roxygenise("./Rpackage/IMPACTncd_Engl_model_pkg/", clean = TRUE)
  # remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/", build_vignettes = TRUE, force = TRUE, upgrade = "never")