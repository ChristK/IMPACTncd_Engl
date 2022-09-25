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

cat("Initialising IMPACTncd_Engl model...\n\n")
if (interactive() && !nzchar(system.file(package = "CKutils"))) {
  if (!nzchar(system.file(package = "remotes"))) install.packages("remotes")
  remotes::install_github("ChristK/CKutils", force = TRUE, upgrade = "never")
}

library(CKutils)
options(rgl.useNULL = TRUE)  # suppress error by demography in rstudio server
options(future.fork.enable = TRUE) # TODO remove for production
options(future.rng.onMisuse = "ignore") # Remove false warning
options(datatable.verbose = FALSE)
options(datatable.showProgress = FALSE)

dependencies(yaml::read_yaml("./dependencies.yaml"))

if (interactive()) {
  snfile <- "./Rpackage/.IMPACTncd_Engl_model_pkg_snapshot.qs"
  if (file.exists(snfile)) snapshot <- changedFiles(qread(snfile))


  if (!nzchar(system.file(package = "IMPACTncdEngl")) ||
      !file.exists(snfile) || any(nzchar(snapshot$added),
        nzchar(snapshot$deleted),
        nzchar(snapshot$changed))) {
    if (!nzchar(system.file(package = "remotes")))
      install.packages("remotes")
    if (nzchar(system.file(package = "roxygen2")))
      roxygen2::roxygenise("./Rpackage/IMPACTncd_Engl_model_pkg/", clean = TRUE)
    detach_package <- function(pkg, character.only = FALSE)
    {
      if(!character.only)
      {
        pkg <- deparse(substitute(pkg))
      }
      search_item <- paste("package", pkg, sep = ":")
      while(search_item %in% search())
      {
        detach(search_item, unload = TRUE, character.only = TRUE)
      }
    }
    detach_package(IMPACTncdEngl)
    file.remove(list.files("./Rpackage/IMPACTncd_Engl_model_pkg/", pattern = ".o$|.dll&|.so&", recursive = TRUE, full.names = TRUE))
    remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/",
      force = TRUE,
      upgrade = "never")

    if (file.exists(snfile)) file.remove(snfile)
    qsave(
      fileSnapshot(
        "./Rpackage/IMPACTncd_Engl_model_pkg/",
        timestamp = NULL,
        md5sum = TRUE,
        recursive = TRUE
      ),
      snfile
    )
  }
}
library(IMPACTncdEngl)
