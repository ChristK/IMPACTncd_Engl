## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
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


setwd("/home/ckypridMy_ModelsIMPACTncd_Engl")
# For ages 20 to 90
univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- TRUE
seed                <- 43L


if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("qs", "fst", "MASS", "splines", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multiprocess)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}

dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90), .(
  education, age, agegrp10, sex, qimd, ethnicity, sha, wt_int, year)]
)
dt[, age := scale(age, 54.52, 15.28)]
dt[, education := ordered(education)]
set.seed(seed)

# formula from qs::qread("./secure_data/lifecourse_models/education_model.qs")
education_model <- polr(
  education ~ ns(age, 4) + qimd + year + sha + sex +
    ethnicity + ns(age, 4):sex + year:sex + ns(age, 4):qimd,
  weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = TRUE
)

education_model$data <- copy(dt)

qsave(education_model, "./secure_data/lifecourse_models/education_model_for_imputation.qs", preset = "high")
print("Model saved")

trms <- all.vars(formula(education_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:14, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd),
              sha = unique(dt$sha), ethnicity = unique(dt$ethnicity))
newdata[, age := scale(age_int, 54.52, 15.28)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("ed", 1:7)) := data.table(rowCumsums(predict(education_model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata)
newdata[, age := age_int]
newdata[, c("age_int", "ed7") := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/education_table_for_imputation.fst", 100L)

print("Table saved")

if (plots) {

  }


