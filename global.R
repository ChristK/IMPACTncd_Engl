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

# Install CKutils from GitHub, or update it when its default branch has moved:
# the install is skipped while the installed commit (RemoteSha) still matches
# the remote one, so this contacts GitHub on every run.
#
# Some corporate networks block HTTPS to github.com but allow SSH. There,
# CKutils is cloned over SSH with the system git and installed from that clone.
# remotes::install_git() does both and -- unlike a manual clone followed by
# install_local() -- records the commit it installed, against which the next
# run compares `git ls-remote`. SSH needs git on the PATH and an SSH key
# registered with a GitHub account: GitHub has no anonymous SSH access, even
# to public repos.
#
# The route is set by the environment variable IMPACTNCD_CKUTILS_INSTALL
# (e.g. in ~/.Renviron):
#   auto  (default) HTTPS, falling back to SSH if that fails;
#   https HTTPS only;
#   ssh   SSH only, so an SSH-only network does not wait on HTTPS every run.
# If every route fails but CKutils is installed, the installed version is used
# with a warning, so a session can still start offline.
# Tests: auxil/test_global_ckutils_install.R
install_ckutils <- function(
  quiet = TRUE,
  mode = Sys.getenv("IMPACTNCD_CKUTILS_INSTALL", "auto"),
  ssh_url = "git@github.com:ChristK/CKutils.git"
) {
  mode <- tolower(trimws(mode))
  if (!nzchar(mode)) {
    mode <- "auto" # the variable is set, but empty
  }
  if (!mode %in% c("auto", "https", "ssh")) {
    stop(
      "IMPACTNCD_CKUTILS_INSTALL must be \"auto\", \"https\" or \"ssh\", not \"",
      mode,
      "\".",
      call. = FALSE
    )
  }

  routes <- list(
    https = function() {
      remotes::install_github(
        "ChristK/CKutils",
        upgrade = "never",
        force = FALSE,
        quiet = quiet
      )
    },
    ssh = function() {
      if (!nzchar(Sys.which("git"))) {
        stop("git was not found on the PATH.", call. = FALSE)
      }
      # Compare commits here: remotes cannot tell whether an SSH install from
      # GitHub is current. It reads the package name with `git archive
      # --remote`, which GitHub refuses, so it would re-clone and recompile
      # on every run. This also keeps git's own error (e.g. a refused key),
      # which remotes discards.
      remote_head <- suppressWarnings(system2(
        "git",
        c("ls-remote", ssh_url, "HEAD"),
        stdout = TRUE,
        stderr = TRUE
      ))
      if (!is.null(attr(remote_head, "status"))) {
        stop(paste(remote_head, collapse = "\n"), call. = FALSE)
      }
      remote_head <- sub("\t.*", "", grep("\tHEAD$", remote_head, value = TRUE))
      installed <- suppressWarnings(
        utils::packageDescription("CKutils", fields = "RemoteSha")
      )
      if (length(remote_head) == 1L && identical(remote_head, installed)) {
        if (!quiet) {
          message(
            "Skipping install of 'CKutils' over SSH, the SHA1 (",
            substr(remote_head, 1L, 8L),
            ") has not changed since last install."
          )
        }
        return("CKutils")
      }
      remotes::install_git(
        ssh_url,
        git = "external", # git2r is often built without SSH support
        upgrade = "never",
        force = FALSE,
        quiet = quiet
      )
    }
  )
  if (mode != "auto") {
    routes <- routes[mode]
  }

  failures <- character()
  for (route in names(routes)) {
    res <- tryCatch(routes[[route]](), error = identity)
    # remotes returns NA, rather than raising an error, when the download or
    # clone itself fails -- silently when quiet = TRUE -- so test for both.
    if (!inherits(res, "error") && !anyNA(res)) {
      return(invisible(route))
    }
    failures[[route]] <- if (inherits(res, "error")) {
      conditionMessage(res)
    } else {
      "the download failed."
    }
    if (route == "https" && length(routes) > 1L) {
      message(
        "Installing CKutils from GitHub over HTTPS failed; trying SSH (",
        ssh_url,
        "). Set IMPACTNCD_CKUTILS_INSTALL=ssh to go straight to SSH."
      )
    }
  }

  failures <- paste0(
    "  [",
    names(failures),
    "] ",
    gsub("\n", "\n    ", failures, fixed = TRUE),
    collapse = "\n"
  )
  if (nzchar(system.file(package = "CKutils"))) {
    warning(
      "Could not install or update CKutils from GitHub; using the installed ",
      "version ",
      format(utils::packageVersion("CKutils")),
      ".\n",
      failures,
      call. = FALSE
    )
    return(invisible("installed"))
  }
  stop(
    "CKutils is not installed and could not be installed from GitHub.\n",
    failures,
    "\nHTTPS needs access to github.com and api.github.com. SSH needs git on ",
    "the PATH and an SSH key registered with GitHub (test with: ssh -T ",
    "git@github.com). On an SSH-only network, set ",
    "IMPACTNCD_CKUTILS_INSTALL=ssh in ~/.Renviron.",
    call. = FALSE
  )
}
install_ckutils(quiet = !dev_mode)
rm(install_ckutils)
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
#
# Defensive unload BEFORE the reinstall check: if IMPACTncdEngland is already
# loaded in this session (e.g. when re-sourcing global.R interactively), unload
# it now, while its on-disk files still match what is in memory. Reinstalling a
# loaded compiled package and then using it in the same session corrupts the
# lazy-load database, because R cannot hot-swap a loaded package -- and it is
# then too late to unload (unloadNamespace() runs .onUnload, which reads the
# already-replaced .rdb). The clean unload therefore MUST happen before the
# reinstall, never after. (installLocalPackageIfChanged itself no longer loads
# the namespace during its check; this guard covers the case where an earlier
# statement in the session did.)
if (isNamespaceLoaded("IMPACTncdEngland")) {
  message("Unloading already-loaded IMPACTncdEngland before reinstall check.")
  try(
    detach("package:IMPACTncdEngland", unload = TRUE, character.only = TRUE),
    silent = TRUE
  )
  try(
    library.dynam.unload(
      "IMPACTncdEngland",
      system.file(package = "IMPACTncdEngland")
    ),
    silent = TRUE
  )
  try(unloadNamespace("IMPACTncdEngland"), silent = TRUE)
}
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
