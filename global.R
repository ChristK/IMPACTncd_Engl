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

# Set a CRAN mirror if not already configured.
# On RHEL-based Linux, use Posit Package Manager for binary packages;
# otherwise fall back to the generic cloud mirror.

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
rm(repos)

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
    "User library path is not writable: ",
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
rm(user_lib)

cat("Initialising IMPACTncd_England model...\n\n")

# Prefer pre-built binary packages even when a newer source-only version is
# on CRAN (affects Windows/macOS, where type = "both" is the default).
# Compiling the newer source of a dependency mixes it with binaries built
# against the previous ABI, producing packages whose DLLs fail to load
# (e.g. 2026-07-23: RcppParallel 6.0.0 switched to oneTBB while the
# stringfish/qs2 binaries were still built against 5.x, which shipped
# tbb.dll). End users may also lack the build tools for source installs.
options(install.packages.compile.from.source = "never")

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
options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

# Increase download timeout for large packages (e.g. duckdb) on CI
options(timeout = max(300, getOption("timeout")))

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
    # update = FALSE prevents updating already installed packages. CKutils
    # (>= 3cffe8e) re-checks installed status per package, so it no longer
    # reinstalls a dependency that was loaded as a side effect of an earlier
    # one -- the previous Windows "package in use" warning muffle is gone.
    tryCatch(
      CKutils::dependencies(pkg_list, update = FALSE),
      error = function(e) {
        # During CRAN dependency transitions, pre-built binaries can be ABI
        # mismatched for a few days (e.g. 2026-07-23: the stringfish/qs2
        # binaries were built against RcppParallel 5.x, which shipped tbb.dll,
        # while RcppParallel 6.0.0 switched to oneTBB). The package then
        # installs fine but its DLL fails to load ("unable to load shared
        # object"). Recover by recompiling the affected packages from source
        # so they link against the versions actually installed (requires
        # Rtools on Windows). Install real sources from a CRAN mirror: Posit
        # Package Manager serves Linux binaries disguised as source tarballs,
        # which would reinstall the same broken build despite type = "source".
        msg <- conditionMessage(e)
        if (!grepl("unable to load shared object", msg)) {
          stop(e)
        }
        # R quotes names with ASCII or Unicode directional quotes (U+2018/19)
        # depending on the locale; accept both. The \u escapes keep this file
        # ASCII-only.
        lq <- "['\u2018]"
        rq <- "['\u2019]"
        nq <- "[^'\u2018\u2019]+" # a run free of quote characters
        # The package whose shared object failed to load, from its path
        # (<lib>/<pkg>/libs[/x64]/<dll>) ...
        dll_dir <- dirname(sub(
          paste0(".*unable to load shared object ", lq, "(", nq, ")", rq, ".*"),
          "\\1", msg
        ))
        if (basename(dll_dir) %in% c("x64", "i386")) {
          dll_dir <- dirname(dll_dir)
        }
        broken <- basename(dirname(dll_dir))
        # ... and the package whose namespace failed as a consequence
        if (grepl("namespace load failed for", msg)) {
          broken <- unique(c(broken, sub(
            paste0(".*namespace load failed for ", lq, "(", nq, ")", rq, ".*"),
            "\\1", msg
          )))
        }
        message(
          "Binary package(s) failed to load; reinstalling from source: ",
          paste(broken, collapse = ", ")
        )
        install.packages(
          broken,
          type = "source",
          repos = c(CRAN = "https://cloud.r-project.org")
        )
        CKutils::dependencies(pkg_list, update = FALSE) # retry once
      }
    )
  }
  rm(pkg_list, pkg_list_file) # Clean up
} else {
  warning("r-packages.txt not found at: ", pkg_list_file)
  rm(pkg_list_file)
}

# Install the local R package if its source code has changed
# Uses a snapshot file to track changes
# Assumes the working directory is the project root
suppressPackageStartupMessages(
  installLocalPackageIfChanged(
    pkg_path = "./Rpackage/IMPACTncd_England_model_pkg/",
    snapshot_path = "./Rpackage/.IMPACTncd_England_model_pkg_snapshot.rds",
    debug = dev_mode # Set to TRUE for debug builds with -O0, FALSE for production builds with -O2
  )
)

if (dev_mode) {
  library(IMPACTncdEngland)
  library(gamlss.dist) # necessary for distributions in Exposure class
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

# Initialise .Random.seed in .GlobalEnv. Without this, a fresh session
# (or one cleared with rm(list = ls(all = TRUE))) fails with:
#   Error in get(".Random.seed", .GlobalEnv): object '.Random.seed' not found
invisible(runif(1))
