library(data.table) # for fast data manipulation
library(fst) # Fast way to save and load data.tables
library(gamlss)
library(qs)

disnm <- "nonmodelled" # disease name
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
strata_ftlt <- c("year", "age", "sex", "dimd", "sha", "ethnicity")
source(paste0("/mnt/", Sys.info()[["user"]], "/UoL/CPRD2021/epi_models/scripts/aux_fn.R"))


# Case Fatality 1st year ====
if (overwrite_ftlt ||
    # !file.exists(output_path(paste0(disnm, "_ftlt1.qs"))) ||
    !file.exists(output_path(paste0(disnm, "_ftlt2.qs")))) {


  # 2+ years case fatality
  dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                           as.data.table = TRUE)[gender != "I"]
  )[between(age, 20, 100)  &
      year < 2020, .SD, .SDcols = c(strata_ftlt, "death_cause", "diedinyr", "diedinCPRD")
  ][, .(ftlt2 = sum(diedinyr == 1L & is.na(death_cause), na.rm = TRUE),
        n = .N), keyby = strata_ftlt]
  dt[, no_ftlt2 := n - ftlt2]
  dt[, year := year - 2000]
  y <- cbind(dt$ftlt2, dt$no_ftlt2)
  dt[, c("ftlt2", "no_ftlt2", "n") := NULL]
  mod_max <- gamlss(
    y ~ (log(year) + pb(age) + pcat(sex) + dimd + pcat(sha) + pcat(ethnicity)) ^ 2,
    family = BI(),
    data = dt,
    method = mixed(20, 100)
  )
  validate_plots(dt, y, mod_max, "_ftlt2", disnm, strata_ftlt)
  qsave(mod_max, output_path(paste0(disnm, "_ftlt2.qs")), "archive")
  print(paste0(disnm, "_ftlt2 model saved!"))

  trms <- all.vars(formula(mod_max))[-1] # -1 excludes dependent var
  newdata2 <-
    CJ(
      age = 20:100,
      year = 8:100,
      sex = levels(dt$sex),
      dimd = levels(dt$dimd),
      ethnicity = levels(dt$ethnicity),
      sha = levels(dt$sha)
    )
  newdata2 <- split(newdata2, by = "dimd")
  newdata2 <-
    # assignment necessary! Copies of data.tables are happening
    lapply(newdata2, function(x)
      x[, c("mu2") := predictAll(mod_max, .SD, data = dt), .SDcols = trms])
  newdata2 <- rbindlist(newdata2)


  newdata2[, dimd := factor(dimd, as.character(1:10))]
  setkeyv(newdata2, strata_ftlt)
  write_fst(newdata2, output_path(paste0(disnm, "_ftlt.fst")), 100L)
  print(paste0(disnm, "_ftlt model saved!"))
  rm(dt, mod_max, newdata1, newdata2, trms)
}
