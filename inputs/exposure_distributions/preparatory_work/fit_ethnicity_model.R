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



setwd("/home/ckyprid/My Models/IMPACTncd_Engl")

# For ages 20 to 90
# Only used to impute HSE_ts during HSE_ts correlation structure extraction

univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- TRUE
seed                <- 1L


if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
if (!require(IMPACTncdEnglmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/")
  library(IMPACTncdEnglmisc)
}
dependencies(c("qs", "fst", "MASS", "nnet", "splines", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multiprocess)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}

dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90), .(
  ethnicity, age, agegrp10, sex, qimd, sha, wt_int, year)]
)
# dt[, age := scale(age, 54.52, 15.28)]
set.seed(seed)

if (univariate_analysis) {
  # age
  age_scaled <- 20:90
  dt[, .(ethnicity_mean = wtd.mean(as.integer(ethnicity), weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, ethnicity_mean, ylim = c(1, 3))]

  m_age1 <- multinom(
    ethnicity ~ age,
    weights = dt$wt_int,
    data = dt,
    Hess = TRUE
  )
  tt <- predict(m_age1, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "blue1")

  m_age2 <- multinom(
    ethnicity ~ poly(logb(age, 100), 3),
    weights = dt$wt_int,
    data = dt,
    Hess = TRUE
  )
  tt <- predict(m_age2, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "red1")

  m_age3 <- multinom(
    ethnicity ~ logb(age, 100),
    weights = dt$wt_int,
    data = dt,
    Hess = TRUE
  )
  tt <- predict(m_age3, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "green1")

  setDT(BIC(m_age1, m_age2, m_age3), keep.rownames = TRUE, key = "BIC")[] # m_age3

  # year
  dt[, .(ethnicity_mean = wtd.mean(as.integer(ethnicity), weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, ethnicity_mean, xlim = c(3, 30), ylim = c(1, 9))]

  m_year1 <- multinom(
    ethnicity ~ year,
    weights = dt$wt_int,
    data = dt,
    Hess = TRUE
  )
  tt <- predict(m_year1 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "blue1")

  m_year2 <- multinom(
    ethnicity ~ logb(year, 20),
    weights = dt$wt_int,
    data = dt,
    Hess = TRUE
  )
  tt <- predict(m_year2 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "red1")

  setDT(BIC(m_year1, m_year2), keep.rownames = TRUE, key = "BIC")[] # m_year2
}

ethnicity_model <- multinom(
  ethnicity ~ logb(year, 10) +  poly(logb(age, 100), 3) + qimd + sex + sha,
  weights = dt$wt_int,
  data = dt,
  maxit = 1e3,
  Hess = TRUE
)

ethnicity_model$data <- copy(dt)
qsave(ethnicity_model, "./secure_data/lifecourse_models/ethnicity_model.qs", preset = "high")
print("Model saved")

trms <- all.vars(formula(ethnicity_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:14, age = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd),
              sha = unique(dt$sha))

newdata[, (levels(dt$ethnicity)) := data.table(matrixStats::rowCumsums(predict(ethnicity_model, type = "p", newdata = .SD))), .SDcols = trms]
# newdata[, c("other") := NULL]
write_fst(newdata, "./inputs/exposure_distributions/ethnicity_table.fst", 100L)

print("Table saved")

if (plots) {
 # only used for imputation of ~300 cases

}



