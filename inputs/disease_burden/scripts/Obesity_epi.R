library(data.table) # for fast data manipulation
library(fst) # Fast way to save and load data.tables
library(gamlss)
library(qs)

disnm <- "Obesity" # disease name
disnm2 <- "obesity" # disease name for saving
overwrite_incd <- TRUE
overwrite_prvl <- TRUE
overwrite_ftlt <- TRUE
overwrite_dur  <- TRUE
overwrite_pred <- FALSE

# disease list
# "Anxiety_Depression"            "Asthma"
# "Atrial Fibrillation"           "CHD"
# "COPD"                          "Chronic Kidney Disease"
# "Dementia"                      "Heart failure"
# "Hypertension"                  "Obesity"
# "Other cancers"                 "Primary Malignancy_Breast"
# "Primary Malignancy_Colorectal" "Primary Malignancy_Lung"
# "Primary Malignancy_Prostate"   "Stroke"
# "Type 2 Diabetes Mellitus"

strata <- c("year", "age", "sex", "dimd", "sha", "ethnicity")
strata_ftlt <- c("year", "age", "sex", "dimd")
source(paste0("/mnt/", Sys.info()[["user"]], "/UoL/CPRD2021/epi_models/scripts/aux_fn.R"))


# Duration ====
if (overwrite_dur ||
    !file.exists(output_path(paste0(disnm2, "_dur.qs")))) {

  dt <- harmonise(read_fst(input_path("panel_short_prev_2018_years.fst"),
                           as.data.table = TRUE)[gender != "I"])[between(age, 20, 100) &
                                                                   get(paste0(disnm, "_years")) >= 2L,
                                                                 .SD, .SDcols = c(paste0(disnm, "_years"), strata)]

  setnames(dt, paste0(disnm, "_years"), "dur")
  dt[, dur := dur - 2L] # during the sim will add 2
  #marg_distr <- fitDist(
  #  dt$dur,
  #  log(nrow(dt)),
  #  type = "count", # "realplus",
  #  try.gamlss = TRUE,
  #  trace = TRUE
  #)
  #head(marg_distr$fits)

  # workHORSEmisc::distr_validation(marg_distr, dt[between(dur, 0, 50), .(var = dur, wt = 1)],
  #                  expression(bold(duration ~ (years))), discrete = TRUE)

#  distr_nam <- names(marg_distr$fits[1]) # pick appropriately and note here ZANBI

  dur_model <- gamlss(
    dur ~ pb(age) + pcat(sex) + pcat(dimd) + pcat(ethnicity),
    ~pb(age) + pcat(sex) + pcat(dimd),
    ~pb(age),
    family = "ZINBI",
    data = dt,
    method = mixed(20, 100)
  )

  qsave(dur_model, output_path(paste0(disnm2, "_dur.qs")), "archive")
  print(paste0(disnm, "_dur model saved!"))

  trms <- all.vars(formula(dur_model))[-1] # -1 excludes dependent var
  newdata <-
    CJ(
      age = 20:100,
      sex = levels(dt$sex),
      dimd = levels(dt$dimd),
      ethnicity = levels(dt$ethnicity)
    )
  newdata <- split(newdata, by = "dimd")
  newdata <-
    # assignment necessary! Copies of data.tables are happening
    lapply(newdata, function(x)
      x[, c("mu", "sigma", "nu") := predictAll(dur_model, .SD, data = dt), .SDcols = trms])
  newdata <- rbindlist(newdata)
  newdata[, dimd := factor(dimd, as.character(1:10))]
  setkeyv(newdata, c("age", "sex", "dimd", "ethnicity"))
  write_fst(newdata, output_path(paste0(disnm2, "_dur.fst")), 100L)
  print(paste0(disnm, "_dur table saved!"))
  rm(dt, dur_model, newdata, trms)
}
# Incidence ====
if (overwrite_incd ||
    !file.exists(output_path(paste0(disnm2, "_incd.qs")))) {
  dt <- harmonise(read_fst(input_path("panel_short_inc.fst"),
                           as.data.table = TRUE)[gender != "I"])[between(age, 20, 100) &
                                                                   get(disnm) < 2L &
                                                                   year < 2020, .SD, .SDcols = c(disnm, strata)][, .(incd = sum(get(disnm), na.rm = TRUE), n = .N),
                                                                                                                 keyby = strata]
  dt[, no_incd := n - incd]
  dt[, year := year - 2000]
  y <- cbind(dt$incd, dt$no_incd)
  dt[, c("incd", "no_incd", "n") := NULL]
  mod_max <- gamlss(
    y ~ (
      log(year) + pb(age) + pcat(sex) + pcat(dimd) + pcat(sha) + pcat(ethnicity)
    ) ^ 2,
    family = BI(),
    data = dt,
    method = mixed(20, 100)
  )
  validate_plots(dt, y, mod_max, "_incd", disnm, strata)
  qsave(mod_max, output_path(paste0(disnm2, "_incd.qs")), "archive")
  print(paste0(disnm, "_incd model saved!"))

  trms <- all.vars(formula(mod_max))[-1] # -1 excludes dependent var
  newdata <-
    CJ(
      age = 20:100,
      year = 3:100,
      sex = levels(dt$sex),
      dimd = levels(dt$dimd),
      ethnicity = levels(dt$ethnicity),
      sha = levels(dt$sha)
    )
  newdata <- split(newdata, by = "sha")
  newdata <-
    # assignment necessary! Copies of data.tables are happening
    lapply(newdata, function(x)
      x[, c("mu") := predictAll(mod_max, .SD, data = dt), .SDcols = trms])
  newdata <- rbindlist(newdata)
  newdata[, dimd := factor(dimd, as.character(1:10))]
  setkeyv(newdata, strata)
  write_fst(newdata, output_path(paste0(disnm2, "_incd.fst")), 100L)
  print(paste0(disnm, "_incd table saved!"))
  rm(dt, mod_max, newdata, trms)
}

# Prevalence ====
if (overwrite_prvl ||
    !file.exists(output_path(paste0(disnm2, "_prvl.qs")))) {
  dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                           as.data.table = TRUE)[gender != "I"])[between(age, 20, 100) &
                                                                   get(disnm) <= 2L &
                                                                   year < 2020, .SD, .SDcols = c(disnm, strata)][, .(prvl = sum(get(disnm) == 2L, na.rm = TRUE), n = .N),
                                                                                                                 keyby = strata]
  dt[, no_prvl := n - prvl]
  dt[, year := year - 2000]
  y <- cbind(dt$prvl, dt$no_prvl)
  dt[, c("prvl", "no_prvl", "n") := NULL]
  mod_max <- gamlss(
    y ~ (
      log(year) + pb(age) + pcat(sex) + pcat(dimd) + pcat(sha) + pcat(ethnicity)
    ) ^ 2,
    family = BI(),
    data = dt,
    method = mixed(20, 100)
  )
  validate_plots(dt, y, mod_max, "_prvl", disnm, strata)
  qsave(mod_max, output_path(paste0(disnm2, "_prvl.qs")), "archive")
  print(paste0(disnm, "_prvl model saved!"))

  trms <- all.vars(formula(mod_max))[-1] # -1 excludes dependent var
  newdata <-
    CJ(
      age = 20:100,
      year = 3:100,
      sex = levels(dt$sex),
      dimd = levels(dt$dimd),
      ethnicity = levels(dt$ethnicity),
      sha = levels(dt$sha)
    )
  newdata <- split(newdata, by = "dimd")
  newdata <-
    # assignment necessary! Copies of data.tables are happening
    lapply(newdata, function(x)
      x[, c("mu") := predictAll(mod_max, .SD, data = dt), .SDcols = trms])
  newdata <- rbindlist(newdata)
  newdata[, dimd := factor(dimd, as.character(1:10))]
  setkeyv(newdata, strata)
  write_fst(newdata, output_path(paste0(disnm2, "_prvl.fst")), 100L)
  print(paste0(disnm, "_prvl table saved!"))
  rm(dt, mod_max, newdata, trms)
}

# Case Fatality 1st year ====
if (overwrite_ftlt ||
    !file.exists(output_path(paste0(disnm2, "_ftlt1.qs"))) ||
    !file.exists(output_path(paste0(disnm2, "_ftlt2.qs")))) {
  dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                           as.data.table = TRUE)[gender != "I"]
                  )[between(age, 20, 100) & get(disnm) == 1L &
                     year < 2020, .SD, .SDcols = c(disnm, strata_ftlt, "death_cause")
      ][, .(ftlt1 = sum(death_cause == disnm, na.rm = TRUE),
                                 n = .N), keyby = strata_ftlt]
  dt[, no_ftlt1 := n - ftlt1]
  dt[, year := year - 2000]
  y <- cbind(dt$ftlt1, dt$no_ftlt1)
  dt[, c("ftlt1", "no_ftlt1", "n") := NULL]
  mod_max <- gamlss(
    y ~ (log(year) + pb(age) + pcat(sex) + dimd) ^ 2,
    family = BI(),
    data = dt,
    method = mixed(20, 100)
  )
  validate_plots(dt, y, mod_max, "_ftlt1", disnm, strata_ftlt)
  qsave(mod_max, output_path(paste0(disnm2, "_ftlt1.qs")), "archive")
  print(paste0(disnm, "_ftlt1 model saved!"))

  trms <- all.vars(formula(mod_max))[-1] # -1 excludes dependent var
  newdata1 <-
    CJ(
      age = 20:100,
      year = 3:100,
      sex = levels(dt$sex),
      dimd = levels(dt$dimd)
    )
  newdata1 <- split(newdata1, by = "dimd")
  newdata1 <-
    # assignment necessary! Copies of data.tables are happening
    lapply(newdata1, function(x)
      x[, c("mu1") := predictAll(mod_max, .SD, data = dt), .SDcols = trms])
  newdata1 <- rbindlist(newdata1)

  # 2+ years case fatality
  dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                           as.data.table = TRUE)[gender != "I"]
  )[between(age, 20, 100) & get(disnm) == 2L &
      year < 2020, .SD, .SDcols = c(disnm, strata_ftlt, "death_cause")
  ][, .(ftlt2 = sum(death_cause == disnm, na.rm = TRUE),
        n = .N), keyby = strata_ftlt]
  dt[, no_ftlt2 := n - ftlt2]
  dt[, year := year - 2000]
  y <- cbind(dt$ftlt2, dt$no_ftlt2)
  dt[, c("ftlt2", "no_ftlt2", "n") := NULL]
  mod_max <- gamlss(
    y ~ (log(year) + pb(age) + pcat(sex) + dimd) ^ 2,
    family = BI(),
    data = dt,
    method = mixed(20, 100)
  )
  validate_plots(dt, y, mod_max, "_ftlt2", disnm, strata_ftlt)
  qsave(mod_max, output_path(paste0(disnm2, "_ftlt2.qs")), "archive")
  print(paste0(disnm, "_ftlt2 model saved!"))

  trms <- all.vars(formula(mod_max))[-1] # -1 excludes dependent var
  newdata2 <-
    CJ(
      age = 20:100,
      year = 3:100,
      sex = levels(dt$sex),
      dimd = levels(dt$dimd)
    )
  newdata2 <- split(newdata2, by = "dimd")
  newdata2 <-
    # assignment necessary! Copies of data.tables are happening
    lapply(newdata2, function(x)
      x[, c("mu2") := predictAll(mod_max, .SD, data = dt), .SDcols = trms])
  newdata2 <- rbindlist(newdata2)

  stopifnot(nrow(newdata1) == nrow(newdata2))
  stopifnot(identical(newdata1[, .SD, .SDcols = -"mu1"],
                      newdata2[, .SD, .SDcols = -"mu2"]))

  newdata1[newdata2, on = .NATURAL, mu2 := i.mu2]
  newdata1[, dimd := factor(dimd, as.character(1:10))]
  setkeyv(newdata1, strata_ftlt)
  write_fst(newdata1, output_path(paste0(disnm2, "_ftlt.fst")), 100L)
  print(paste0(disnm, "_ftlt model saved!"))
  rm(dt, mod_max, newdata1, newdata2, trms)
}



# Longer predictions ----
if (overwrite_pred) {
  template <-
    CJ(
      age = 20:100,
      year = 3:100,
      sex = factor(c("men", "women")),
      dimd = factor(1:10),
      ethnicity = factor(
        c(
          "white",
          "indian",
          "pakistani",
          "bangladeshi",
          "other asian",
          "black caribbean",
          "black african",
          "chinese",
          "other"
        )
      ),
      sha = factor(
        c(
          "North East",
          "North West",
          "Yorkshire and the Humber",
          "East Midlands",
          "West Midlands",
          "East of England",
          "London",
          "South East Coast",
          "South Central",
          "South West"
        )
      )
    )


  for (i in c("_incd", "_prvl")) {

    if (i == "_incd") {
      dt <- harmonise(read_fst(input_path("panel_short_inc.fst"),
                               as.data.table = TRUE)[gender != "I"])[between(age, 20, 100) &
                                                                       get(disnm) < 2L &
                                                                       year < 2020, .SD, .SDcols = c(disnm, strata)][, .(incd = sum(get(disnm), na.rm = TRUE), n = .N),
                                                                                                                     keyby = strata]
    }

    if (i == "_prvl") {
      dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                               as.data.table = TRUE)[gender != "I"])[between(age, 20, 100) &
                                                                       get(disnm) <= 2L &
                                                                       year < 2020, .SD, .SDcols = c(disnm, strata)][, .(prvl = sum(get(disnm) == 2L, na.rm = TRUE), n = .N),
                                                                                                                     keyby = strata]
    }
    dt[, year := year - 2000]


    mod_max <- qread(output_path(paste0(disnm2, i, ".qs")))
    trms <-
      all.vars(formula(mod_max))[-1] # -1 excludes dependent var
    newdata <- copy(template)
    newdata <- split(newdata, by = "sha")
    newdata <- lapply(newdata, function(x)
      x[, c("mu") := predictAll(mod_max, .SD, data = dt),
        .SDcols = trms])
    newdata <- rbindlist(newdata)
    newdata[, dimd := factor(dimd, as.character(1:10))]
    setkeyv(newdata, strata)
    write_fst(newdata,
              output_path(paste0(disnm2, i, ".fst")), 100L)
    print(paste0(disnm, " ", i, " table saved!"))
  }
}




