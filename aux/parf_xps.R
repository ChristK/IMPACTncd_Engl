source("./global.R")
design <- Design$new("./aux/sim_design_parf.yaml")

sp <- SynthPop$new(1L, design)
path <- file.path(design$sim_prm$output_dir, "parf")
if (!dir.exists(path)) {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

fl <-
  list.files(path = "./inputs/RR",
             pattern = ".csvy$",
             full.names = TRUE)
RR <-
  future_lapply(fl, Exposure$new, design, future.seed = 950480304L)
names(RR) <- sapply(RR, function(x)
  x$get_name())

it <- c("", unique(sapply(RR, `[[`, "name")), "smoking")
it <- it[!it %in% c("met", "statin_px", "af_prvl", "t2dm_prvl")]
it <- it[grep("^smok_", it, invert = TRUE)]

for (i in it) {
  print(i)
  design$sim_prm$ignore_xps <- i

  filenam <-
    file.path(path, paste0("parf_", design$sim_prm$ignore_xps, ".csv"))

  # RR ----
  # Create a named list of Exposure objects for the files in ./inputs/RR
  RR <-
    future_lapply(fl, Exposure$new, design, future.seed = 950480304L)
  names(RR) <- sapply(RR, function(x)
    x$get_name())
  invisible(future_lapply(RR, function(x) {
    x$gen_stochastic_effect(design, overwrite = FALSE, smooth = FALSE)
  },
  future.seed = 627524136L))
  # # NOTE smooth cannot be exported to Design for now, because the first time
  # # this parameter changes we need logic to overwrite unsmoothed files
  #
  # # Generate diseases ----
  diseases <- lapply(design$sim_prm$diseases, function(x) {
    x[["design_"]] <- design
    x[["RR"]] <- RR
    do.call(Disease$new, x)
  })
  names(diseases) <- sapply(design$sim_prm$diseases, `[[`, "name")

  lapply(diseases, function(x) {
    print(x$name)

    x$gen_parf(
      sp = sp,
      design_ = design,
      diseases_ = diseases,
      popsize = 100,
      check = design$sim_prm$logs,
      keep_intermediate_file = TRUE
    )

    if (sum(dim(x$get_parf())) > 0) {
      parf_dt <-
        cbind(sp$pop[, .(wt_immrtl, age, sex, dimd, ethnicity, sha, year)],
              "parf" = x$get_parf("parf")$parf)[year == design$sim_prm$init_year]
      parf_dt <-
        parf_dt[!is.na(parf), .(parf = unique(parf), pop_size = sum(wt_immrtl)),
                keyby = .(age, sex, dimd, ethnicity, sha)]
      parf_dt[, `:=`(disease = x$name, mc = sp$mc)] # not sp$mc_aggr
      fwrite_safe(parf_dt, filenam)
    }
  })

}

tt  <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_.csv"))
tpa <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_active_days.csv"))
tal <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_alcohol.csv"))
tbm <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_bmi.csv"))
tfr <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_fruit.csv"))
tbp <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_sbp.csv"))
tch <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_tchol.csv"))
tve <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_veg.csv"))
tsm <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_smoking.csv"))
tet <- fread(paste0(design$sim_prm$output_dir, "/parf/parf_ets.csv"))

tt[tpa, on = .(age, sex, dimd, ethnicity, sha, disease), parf_pa  := parf - i.parf]
tt[tal, on = .(age, sex, dimd, ethnicity, sha, disease), parf_alc := parf - i.parf]
tt[tbm, on = .(age, sex, dimd, ethnicity, sha, disease), parf_bmi := parf - i.parf]
tt[tfr, on = .(age, sex, dimd, ethnicity, sha, disease), parf_fru := parf - i.parf]
tt[tbp, on = .(age, sex, dimd, ethnicity, sha, disease), parf_sbp := parf - i.parf]
tt[tch, on = .(age, sex, dimd, ethnicity, sha, disease), parf_tch := parf - i.parf]
tt[tve, on = .(age, sex, dimd, ethnicity, sha, disease), parf_veg := parf - i.parf]
tt[tsm, on = .(age, sex, dimd, ethnicity, sha, disease), parf_smo := parf - i.parf]
tt[tet, on = .(age, sex, dimd, ethnicity, sha, disease), parf_ets := parf - i.parf]

summary(tt)
View(tt[parf_fru < 0])

nam <- grep("^parf_", names(tt), value = TRUE) # remove negative values. Might not be entirely correct
tt[, (nam) := lapply(.SD, function(x) {
  x[x < 0] <- 0
  x
}), .SDcols = nam]

fwrite(tt, paste0(design$sim_prm$output_dir, "/parf/parf_final.csv"))

tt[, parf_correction := parf/Reduce(`+`, .SD), .SDcols = nam]
fwrite(tt, paste0(design$sim_prm$output_dir, "/parf/parf_final_corrected.csv"))

# tt[, .(weighted.mean(parf, pop_size), weighted.mean(parf_smo, pop_size), weighted.mean(parf_ets, pop_size)), keyby = .(disease)]
# tt[, .(weighted.mean(parf, pop_size), weighted.mean(parf_smo * parf_correction, pop_size), weighted.mean(parf_ets * parf_correction, pop_size)), keyby = .(disease)]

