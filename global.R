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
# file.remove(list.files("./Rpackage/IMPACTncd_Engl_model_pkg/src/", full.names = TRUE, recursive = TRUE, pattern = ".o$"))

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


  if (!nzchar(system.file(package = "IMPACTncdEnglmisc")) ||
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
    detach_package(IMPACTncdEnglmisc)
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
      "./Rpackage/.IMPACTncd_Engl_model_pkg_snapshot.qs"
    )
  }
}
library(IMPACTncdEnglmisc)



design <- Design$new("./inputs/sim_design.yaml")
setDTthreads(threads = design$sim_prm$clusternumber, restore_after_fork = NULL)
threads_fst(nr_of_threads = design$sim_prm$clusternumber, reset_after_fork = NULL)
registerDoFuture()





strata_for_gui <- c("mc", "friendly_name", design$sim_prm$strata_for_output)

def_col <- viridis(16, option = "D")
def_col_small <- def_col[c(1, 9, 5, 3, 7, 2, 8, 4, 6)]



def_sym <- plotly::schema(F)$traces$scatter$attributes$marker$symbol$values

def_sym <- sort(grep("^[0-9]+$", def_sym, value = TRUE, invert = TRUE))

def_sym <- unique(c("circle-dot", "square", "diamond", "cross", "triangle-up",
  "pentagon", "star", "hexagon-open-dot", "triangle-down", def_sym))

# Produce scenario tabs 2-9 using scenario 1 as a template if not exist ------
for (i in 2:9) {
  if (!file.exists(paste0("ui/scenario", i, "_tab.R")) ||
      !identical(file.size("ui/scenario1_tab.R"), file.size(paste0("ui/scenario", i, "_tab.R")))) {
    tt <- readLines(file.path("ui", "scenario1_tab.R"))
    tt <- gsub("Scenario 1",
               paste0("Scenario ", i),
               gsub(
                 "input.level == 1",
                 paste0("input.level == ", i),
                 gsub("_sc1",  paste0("_sc", i),
                      gsub("def_col_small\\[\\[1\\]\\]",
                           paste0("def_col_small\\[\\[", i, "\\]\\]"),
                      gsub("def_sym\\[\\[1\\]\\]",
                           paste0("def_sym\\[\\[", i, "\\]\\]"), tt)
                      ))))
    writeLines(tt, paste0("ui/scenario", i, "_tab.R"))
  }
}



# load localities structure
localities_indx <- read_fst("./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst", as.data.table = TRUE)
localitities_list <- list(
  "Country" = list("England"),
  "Region" = sort(as.character((localities_indx[, unique(RGN11NM)]))),
  "Local Authority" = sort(as.character(localities_indx[, unique(LAD17NM)]))
)
rm(localities_indx)



# sum(gtools::rdirichlet(1, rep(1, 180))) # to simulate a 180 length vector that sums to 1
# extract_uptake_table(setDT(readRDS("./output/input_parameters.rds")), 2)[]
# To write when running from shiny server or docker see https://groups.google.com/forum/#!topic/shiny-discuss/srWETT6uL-I


# Table related ----
# default global search value
if (!exists("default_search_tbl_summary")) default_search_tbl_summary <- ""

# default column search values
if (!exists("default_search_columns_tbl_summary")) default_search_columns_tbl_summary <- NULL


# RR ----
# Create a named list of Exposure objects for the files in ./inputs/RR
fl <- list.files(path = "./inputs/RR/", pattern = ".csvy$", full.names = TRUE)
RR <- lapply(fl, Exposure$new)
names(RR) <- sapply(RR, function(x) x$get_name())
invisible(lapply(RR, function(x)
  x$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)))
# NOTE smooth cannot be exported to Design for now, because the first time
# this parameter changes we need logic to overwrite unsmoothed files
rm(fl)

# Generate diseases ----
diseases <- lapply(design$sim_prm$diseases, function(x) {
  do.call(Disease$new, c(x, "design" = design))
})
names(diseases) <- sapply(design$sim_prm$diseases, `[[`, "name")


# Parameters ----
# Create a named list of Exposure objects for the files in ./inputs/RR
# fl <- list.files(path = "./inputs/other_parameters/", pattern = ".csvy$")
# fl <- gsub(".csvy$", "", fl)
# fl2 <- strsplit(fl, "~", TRUE)
# fl2 <- lapply(fl2, function(x) {
#   tt <- x[2]
#   if (length(x) > 2) tt <- paste(x[-1], collapse = "_")
#   return(c(x[1], tt))
# })
#
# prmtrs <- mapply(Parameter$new, lapply(fl2, `[[`, 1), lapply(fl2, `[[`, 2))
# names(prmtrs) <- fl

# Read the files in ./inputs/other_parameters
# pfl <- paste0("./inputs/other_parameters/", fl, ".csvy")
# for (i in seq_along(pfl)) {
#   prmtrs[[i]]$read_exps_prm(pfl[[i]], design)
# }

# Generate stochastic parameterss
# for (i in seq_along(prmtrs)) {
#   prmtrs[[i]]$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)
#   # NOTE smooth cannot be exported to Design for now, because the first time
#   # this parameter changes we need logic to overwrite unsmoothed files
# }
# rm(fl, fl2, pfl)

# dockerfile_object <- containerit::dockerfile()
# dockerfile_object
# containerit::write(dockerfile_object, file = "./dockerfile")

