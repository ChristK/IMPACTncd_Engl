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

setOptions_for_repo <- function() {
  chooseCRANmirror(ind = 1)
  repos <- getOption("repos")
}
setOptions_for_repo()

cat("Initialising IMPACTncd_Engl model...\n\n")
if (!nzchar(system.file(package = "CKutils"))) {
  if (!nzchar(system.file(package = "pak"))) install.packages("pak")
  pak::pkg_install("ChristK/CKutils", upgrade = FALSE, ask = FALSE)
  # force = T
}

library(CKutils)
options(rgl.useNULL = TRUE) # suppress error by demography in rstudio server
options(future.fork.enable = TRUE) # TODO remove for production
options(future.globals.maxSize = +Inf)
options(future.rng.onMisuse = "ignore") # Remove false warning
options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

CKutils::dependencies(yaml::read_yaml("./dependencies.yaml")) # install missing packages

#' @description Detach library package.
detach_package <- function(pkg, character.only = FALSE) {
  if (!character.only) pkg <- deparse(substitute(pkg))
  search_item <- paste("package", pkg, sep = ":")
  while (search_item %in% search()) {
    detach(search_item, unload = TRUE, character.only = TRUE)
  }
}

#' @description Re/install the IMPACTncd_Engl package from local directory.
#' @param sIMPACTncdPackageDirPath string, IMPACTncd_Engl package directory path.
InstallIMPACTncdPackage <- function(sIMPACTncdPackageDirPath) {
  if (!nzchar(system.file(package = "pak"))) install.packages("pak")
  if (nzchar(system.file(package = "roxygen2"))) {
    roxygen2::roxygenise(sIMPACTncdPackageDirPath, clean = TRUE)
  } # update package exports (and docs) if necessary
  detach_package(IMPACTncdEngl)

  # ensure full install by removing intermediate files
  file.remove(list.files(sIMPACTncdPackageDirPath,
    pattern = ".o$|.dll&|.so&", recursive = TRUE,
    full.names = TRUE
  ))
  ## Find a solution to build vignettes while installing the package using pak
  # pak::local_install(sIMPACTncdPackageDirPath,
  #                    upgrade = FALSE,
  #                    ask = FALSE)
  remotes::install_local(sIMPACTncdPackageDirPath,
    build_vignettes = TRUE, force = TRUE, upgrade = "never"
  )
  # build_vignettes = T, force = T
}

#' @description Re/install the IMPACTncd_Engl package if not installed or if local package files have changed.
InstallIMPACTncdPackageOnChange <- function() {
  # load previously saved snapshot of package files
  sIMPACTncdPackageSnapshotFilePath <- "./Rpackage/.IMPACTncd_Engl_model_pkg_snapshot.qs"
  IMPACTncdPackageSnapshot <- NULL
  if (file.exists(sIMPACTncdPackageSnapshotFilePath)) {
    IMPACTncdPackageSnapshot <- changedFiles(qread(sIMPACTncdPackageSnapshotFilePath))
  }

  if (!nzchar(system.file(package = "IMPACTncdEngl")) || is.null(IMPACTncdPackageSnapshot) ||
    any(
      nzchar(IMPACTncdPackageSnapshot$added), nzchar(IMPACTncdPackageSnapshot$deleted),
      nzchar(IMPACTncdPackageSnapshot$changed)
    )) {
    # re/install IMPACTncd_Engl package and update snapshot
    sIMPACTncdPackageDirPath <- "./Rpackage/IMPACTncd_Engl_model_pkg/"
    InstallIMPACTncdPackage(sIMPACTncdPackageDirPath)
    if (!is.null(IMPACTncdPackageSnapshot)) file.remove(sIMPACTncdPackageSnapshotFilePath)
    IMPACTncdPackageSnapshot <- fileSnapshot(sIMPACTncdPackageDirPath,
      timestamp = NULL,
      md5sum = TRUE, recursive = TRUE
    )
    qsave(IMPACTncdPackageSnapshot, sIMPACTncdPackageSnapshotFilePath)
  }
}

# if(interactive())
InstallIMPACTncdPackageOnChange()

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
runif(1)

# roxygen2::roxygenise("./Rpackage/IMPACTncd_Engl_model_pkg/", clean = TRUE)
  # remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/", build_vignettes = TRUE, force = TRUE, upgrade = "never")