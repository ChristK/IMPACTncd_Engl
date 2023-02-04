source("./global.R")
design <- Design$new("./inputs/sim_design.yaml")
# RR ----
# Create a named list of Exposure objects for the files in ./inputs/RR
fl <- list.files(path = "./inputs/RR", pattern = ".csvy$", full.names = TRUE)
# RR <- lapply(fl, Exposure$new, design)
# names(RR) <- sapply(RR, function(x) x$get_name())
# lapply(RR, function(x) {
#     x$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)
# })
RR <- future_lapply(fl, Exposure$new, design,future.seed = 950480304L)
names(RR) <- sapply(RR, function(x) x$get_name())
invisible(future_lapply(RR, function(x) {
  x$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)
},
future.seed = 627524136L))
# NOTE smooth cannot be exported to Design for now, because the first time
# this parameter changes we need logic to overwrite unsmoothed files
rm(fl)
#
# Generate diseases ----
diseases <- lapply(design$sim_prm$diseases, function(x) {
  print(x$name)
  x[["design_"]] <- design
  x[["RR"]] <- RR
  do.call(Disease$new, x)
})
names(diseases) <- sapply(design$sim_prm$diseases, `[[`, "name")

lapply(diseases, function(x) {
  print(x)
  x$gen_parf_files(design)
})
