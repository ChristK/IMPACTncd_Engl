## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
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

# Set development mode flag
dev_mode <- FALSE # Set to FALSE for production

# Ensure a CRAN mirror is set
repos <- getOption("repos")
if (is.null(repos) || repos["CRAN"] == "@CRAN@") {
  # Set default CRAN mirror if not already set
  # Prefer Posit Public Package Manager for Linux binaries if on Linux
  if (Sys.info()["sysname"] == "Linux") {
    # Check if we are on RHEL/CentOS/Rocky/Alma
    if (file.exists("/etc/redhat-release")) {
      options(
        repos = c(
          CRAN = "https://packagemanager.posit.co/cran/__linux__/rhel9/latest"
        )
      )
    } else {
      options(repos = c(CRAN = "https://cloud.r-project.org"))
    }
  } else {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
}

# Define and ensure the user library path exists and is writable
# This is crucial when running in environments (like Docker) where the default
# system library might not be writable by the current user.
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  # Provide a default user library path if R_LIBS_USER is not set
  # Format: ~/R/<platform>/<major>.<minor>
  user_lib <- file.path(
    Sys.getenv("HOME"),
    "R",
    paste0(R.version$platform, "-library"),
    paste0(R.version$major, ".", substr(R.version$minor, 1, 1))
  )
}

# Create the directory if it doesn't exist
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}

# Check if the user library is writable
if (file.access(user_lib, mode = 2) != 0) {
  # mode = 2 checks for write permission
  stop(
    "User library path is not writable:",
    user_lib,
    ". Please check permissions or set the R_LIBS_USER environment variable to a writable path."
  )
}


# Add the user library to the library paths if not already present
# Prepending ensures it's the default location for installations
if (!user_lib %in% .libPaths()) {
  .libPaths(c(user_lib, .libPaths()))
}

if (dev_mode) cat("Using library path:", user_lib, "\n")

# --- End: User Library Path Configuration ---

cat("Initialising IMPACTncd_England model...\n\n")

# Ensure 'remotes' is installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Ensure 'CKutils' is installed from GitHub if missing
remotes::install_github(
  "ChristK/CKutils",
  upgrade = "never",
  force = FALSE,
  quiet = !dev_mode
)
if (dev_mode) {
  library(CKutils)
} else {
  suppressPackageStartupMessages(library(CKutils))
}


# Environment-specific options
options(rgl.useNULL = TRUE) # suppress error by demography in rstudio server
# if (dev_mode) {
#   options(future.fork.enable = TRUE) # enable for development only
#   options(future.globals.maxSize = +Inf)
#   options(future.rng.onMisuse = "ignore") # Remove false warning
# }
options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

# Install missing packages listed in r-packages.txt
# Assumes the working directory is the project root
pkg_list_file <- "./docker_setup/r-packages.txt"
if (file.exists(pkg_list_file)) {
  pkg_list <- readLines(pkg_list_file, warn = FALSE)
  pkg_list <- trimws(pkg_list)
  # Filter out empty lines and comments
  pkg_list <- pkg_list[nzchar(pkg_list) & !grepl("^#", pkg_list)]

  # Filter out packages that are already installed to avoid "in use" errors on Windows
  # and to respect the update = FALSE intent more strictly
  pkg_list <- pkg_list[!pkg_list %in% rownames(installed.packages())]

  if (length(pkg_list) > 0) {
    # update = FALSE prevents updating already installed packages
    CKutils::dependencies(pkg_list, update = FALSE)
  }
  rm(pkg_list, pkg_list_file) # Clean up
} else {
  warning("r-packages.txt not found at: ", pkg_list_file)
}

# Install the local R package if its source code has changed
# Uses a snapshot file to track changes
# Assumes the working directory is the project root
# Use withCallingHandlers to suppress "package in use" warnings on Windows
suppressPackageStartupMessages(withCallingHandlers(
  {
    installLocalPackageIfChanged(
      pkg_path = "./Rpackage/IMPACTncd_England_model_pkg/",
      snapshot_path = "./Rpackage/.IMPACTncd_England_model_pkg_snapshot.rds",
      debug = dev_mode # Set to TRUE for debug builds with -O0, FALSE for production builds with -O2
    )
  },
  warning = function(w) {
    if (grepl("is in use and will not be installed", conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  }
))

if (dev_mode) {
  library(IMPACTncdEngland)
  library(gamlss.dist) # necessary for ditributions in Exposure class
} else {
  suppressPackageStartupMessages(library(IMPACTncdEngland))
  suppressPackageStartupMessages(library(gamlss.dist))
}

data.table::setDTthreads(
  threads = 1L,
  restore_after_fork = FALSE
)
fst::threads_fst(
  nr_of_threads = 1L,
  reset_after_fork = FALSE
)
arrow::set_cpu_count(1L) # limit Arrow's internal threading

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