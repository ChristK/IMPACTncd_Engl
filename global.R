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

design <- Design$new("./inputs/sim_design.yaml")

setDTthreads(threads = design$sim_prm$clusternumber, restore_after_fork = NULL)
threads_fst(nr_of_threads = design$sim_prm$clusternumber, reset_after_fork = NULL)
# registerDoFuture()
plan(multicore, workers = design$sim_prm$clusternumber)


# RR ----
# Create a named list of Exposure objects for the files in ./inputs/RR
fl <- list.files(path = "./inputs/RR", pattern = ".csvy$", full.names = TRUE)
RR <- future_lapply(fl, Exposure$new, future.seed = 950480304L)
names(RR) <- sapply(RR, function(x) x$get_name())
invisible(future_lapply(RR, function(x) {
  x$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)
},
future.seed = 627524136L))
# NOTE smooth cannot be exported to Design for now, because the first time
# this parameter changes we need logic to overwrite unsmoothed files
rm(fl)

# Generate diseases ----
diseases <- lapply(design$sim_prm$diseases, function(x) {
  x[["design_"]] <- design
  x[["RR"]] <- RR
  do.call(Disease$new, x)
})
names(diseases) <- sapply(design$sim_prm$diseases, `[[`, "name")




# Plot causality structure
# ds <- unlist(strsplit(names(RR), "~"))
# ds[grep("^smok_", ds)] <- "smoking"
# ds <- gsub("_prvl$", "", ds)
#
# ds1 <- ds[as.logical(seq_along(ds) %% 2)]
# ds2 <- ds[!as.logical(seq_along(ds) %% 2)]
# ds <- unique(data.table(ds1, ds2))
#
# g <- make_graph(unlist(transpose(ds)), directed = TRUE)
# plot(g, vertex.shape = "none", edge.arrow.size = .3,
#     vertex.label.font = 2, vertex.label.color = "gray40",edge.arrow.width = .7,
#     vertex.label.cex = .7, edge.color = "gray85", layout = layout_components)
rm(RR)

# TODO move to class Simulation
mk_scenario_init2 <- function(scenario_name, diseases_, sp, design_) {
  # scenario_suffix_for_pop <- paste0("_", scenario_name) # TODO get suffix from design
  scenario_suffix_for_pop <- scenario_name
  list(
    "exposures"          = design_$sim_prm$exposures,
    "scenarios"          = design_$sim_prm$scenarios, # to be generated programmatically
    "scenario"           = scenario_name,
    "kismet"             = design_$sim_prm$kismet, # If TRUE random numbers are the same for each scenario.
    "init_year"          = design_$sim_prm$init_year,
    "pids"               = "pid",
    "years"              = "year",
    "ages"               = "age",
    "ageL"               = design_$sim_prm$ageL,
    "all_cause_mrtl"     = paste0("all_cause_mrtl", scenario_suffix_for_pop),
    "strata_for_outputs" = c("pid", "year", "age", "sex", "dimd"),
    "diseases"           = lapply(diseases_, function(x) x$to_cpp(sp, design_))
  )
}

run_sim <- function(mc, sp, diseases, design) {
  sp <- SynthPop$new(mc, design)
  lapply(diseases, function(x) {
    x$gen_parf(sp, design)$
      set_init_prvl(sp, design)$
      set_rr(sp, design)$
      set_incd_prb(sp, design)$
      set_dgns_prb(sp, design)$
      set_mrtl_prb(sp, design)
  })
  l <- mk_scenario_init2("", diseases, sp, design)
  simcpp(sp$pop, l, sp$mc)

  sp$update_pop_weights()

}
