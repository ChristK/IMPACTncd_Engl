library(data.table) # for fast data manipulation
library(fst) # Fast way to save and load data.tables
library(gamlss)
library(qs)

disnm <- "Type 2 Diabetes Mellitus" # disease name
overwrite_incd <- FALSE
overwrite_prvl <- FALSE
overwrite_ftlt <- FALSE

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
source("./inputs/disease_burden/scripts/aux_fn.R")

# Incidence ====
if (overwrite_incd ||
    !file.exists(output_path(paste0(disnm, "_incd.qs")))) {
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
      pb(year) + pb(age) + pcat(sex) + pcat(dimd) + pcat(sha) + pcat(ethnicity)
    ) ^ 2,
    family = BI(),
    data = dt,
    method = mixed(20, 100)
  )
  validate_plots(dt, y, mod_max, "_incd", disnm, strata)
  
  qsave(mod_max, output_path(paste0(disnm, "_incd.qs")), "archive")
  print(paste0(disnm, "_incd model saved!"))

  
  
  
  trms <- all.vars(formula(mod_max))[-1] # -1 excludes dependent var
  newdata <-
    CJ(
      age = 20:100,
      year = unique(dt$year),
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
  write_fst(newdata, output_path(paste0(disnm, "_incd.fst")), 100L)
  print(paste0(disnm, "_incd table saved!"))
  rm(dt, mod_max, newdata, trms)
}





# Prevalence ====
if (overwrite_prvl ||
    !file.exists(output_path(paste0(disnm, "_prvl.qs")))) {
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
      pb(year) + pb(age) + pcat(sex) + pcat(dimd) + pcat(sha) + pcat(ethnicity)
    ) ^ 2,
    family = BI(),
    data = dt,
    method = mixed(20, 100)
  )
  validate_plots(dt, y, mod_max, "_prvl", disnm, strata)
  qsave(mod_max, output_path(paste0(disnm, "_prvl.qs")), "archive")
  print(paste0(disnm, "_prvl model saved!"))

  trms <- all.vars(formula(mod_max))[-1] # -1 excludes dependent var
  newdata <-
    CJ(
      age = 20:100,
      year = unique(dt$year),
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
  write_fst(newdata, output_path(paste0(disnm, "_prvl.fst")), 100L)
  print(paste0(disnm, "_prvl table saved!"))
  rm(dt, mod_max, newdata, trms)
}

# Case Fatality all years ====

if (overwrite_ftlt ||
    !file.exists(output_path(paste0(disnm, "_ftlt.qs")))) 
  
  {
  
  dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                           as.data.table = TRUE)[gender != "I"]
                  )[between(age, 20, 100) & get(disnm) > 0L &
                     year < 2020, .SD, .SDcols = c(disnm, strata_ftlt, "death_cause")
      ][, .(ftlt = sum(death_cause == disnm, na.rm = TRUE),
                                 n = .N), keyby = strata_ftlt]
  
  dt[, no_ftlt := n - ftlt]
  dt[, year := year - 2000]
  
  y <- cbind(dt$ftlt, dt$no_ftlt)
  dt[, c("ftlt", "no_ftlt", "n") := NULL]
  
  mod_max <- gamlss(
    y ~ (log(year) + pb(age) + pcat(sex) + dimd) ^ 2,
    family = BI(),
    data = dt,
    method = mixed(20, 100)
  )
  
  validate_plots(dt, y, mod_max, "_ftlt", disnm, strata_ftlt)
  qsave(mod_max, output_path(paste0(disnm, "_ftlt.qs")), "archive")
  print(paste0(disnm, "_ftlt model saved!"))

  trms <- all.vars(formula(mod_max))[-1] # -1 excludes dependent var
  newdata1 <-
    CJ(
      age = 20:100,
      year = 8:100,
      sex = levels(dt$sex),
      dimd = levels(dt$dimd)
    )
  newdata1 <- split(newdata1, by = "dimd")
  
  newdata1 <-
    # assignment necessary! Copies of data.tables are happening
    lapply(newdata1, function(x)
      x[, c("mu") := predictAll(mod_max, .SD, data = dt), .SDcols = trms])
  newdata1 <- rbindlist(newdata1)


  newdata1[, dimd := factor(dimd, as.character(1:10))]
  setkeyv(newdata1, strata_ftlt)
  write_fst(newdata1, output_path(paste0(disnm, "_ftlt.fst")), 100L)
  print(paste0(disnm, "_ftlt model saved!"))
  rm(dt, mod_max, newdata1, trms)
}
