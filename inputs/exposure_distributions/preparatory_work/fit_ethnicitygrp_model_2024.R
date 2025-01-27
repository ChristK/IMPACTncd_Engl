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



#setwd("/home/ckyprid/My_Models/IMPACTncd_Engl")

# For ages 20 to 90
# Only used to impute HSE_ts during HSE_ts correlation structure extraction

library(IMPACTncdEngl)
library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk
library(nnet) # for multinomial

# Set some variables we will use later in automation
univariable_analysis <- FALSE
plots <- FALSE
seed                <- 1L



# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)


ds<- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90), .(
  ethnicity_grp, age, agegrp10, sex, qimd, sha, urban_rural, wt_int, year)]
)
# ds[, age := scale(age, 54.52, 15.28)]
set.seed(seed)

if (univariable_analysis) {
  # age
  age_scaled <- 20:90
  ds[, .(ethnicity_mean = wtd.mean(as.integer(ethnicity_grp), weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, ethnicity_mean, ylim = c(1, 3))]

  m_age1 <- multinom(
    ethnicity_grp ~ age,
    weights = ds$wt_int,
    data = ds,
    Hess = TRUE
  )
  # tt <- predict(m_age1, type = "p", newdata = data.frame(age = age_scaled))
  # lines(age_scaled, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "blue1")

  m_age2 <- multinom(
    ethnicity_grp ~ poly(logb(age, 100), 3),
    weights = ds$wt_int,
    data = ds,
    Hess = TRUE
  )
  # tt <- predict(m_age2, type = "p", newdata = data.frame(age = age_scaled))
  # lines(age_scaled, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "red1")
  # 
  m_age3 <- multinom(
    ethnicity_grp ~ logb(age, 100),
    weights = ds$wt_int,
    data = ds,
    Hess = TRUE
  )
  # tt <- predict(m_age3, type = "p", newdata = data.frame(age = age_scaled))
  # lines(age_scaled, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "green1")

  setDT(BIC(m_age1, m_age2, m_age3), keep.rownames = TRUE, key = "BIC")[] # m_age2

  # year
  ds[, .(ethnicity_mean = wtd.mean(as.integer(ethnicity_grp), weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, ethnicity_mean, xlim = c(3, 30), ylim = c(1, 9))]

  m_year1 <- multinom(
    ethnicity_grp ~ year,
    weights = ds$wt_int,
    data = ds,
    Hess = TRUE
  )
  # tt <- predict(m_year1 , type = "p", newdata = data.frame(year = 3:30))
  # lines(3:30, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "blue1")

  m_year2 <- multinom(
    ethnicity_grp ~ logb(year, 10),
    weights = ds$wt_int,
    data = ds,
    Hess = TRUE
  )
  # tt <- predict(m_year2 , type = "p", newdata = data.frame(year = 3:30))
  # lines(3:30, apply(tt, 1, function(x) mean(sample(9, 1e5, TRUE, x))), col = "red1")

  setDT(BIC(m_year1, m_year2), keep.rownames = TRUE, key = "BIC")[] # m_year2
}

ethnicity_grp_model <- multinom(
  ethnicity_grp ~ logb(year, 10) +  poly(logb(age, 100), 3) + qimd + sex + sha + urban_rural,
  weights = ds$wt_int,
  data = ds,
  maxit = 1e3,
  Hess = TRUE
)

ethnicity_grp_model$data <- copy(ds)
qsave(ethnicity_grp_model, "./secure_data/lifecourse_models/ethnicity_grp_model.qs", preset = "high")
print("Model saved")

trms <- all.vars(formula(ethnicity_grp_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:19,
              age = 20:90, 
              sex = unique(ds$sex),
              qimd = unique(ds$qimd),
              sha = unique(ds$sha),
              urban_rural = unique(ds$urban_rural))

newdata[, (levels(ds$ethnicity_grp)) := data.table(matrixStats::rowCumsums(predict(ethnicity_grp_model, type = "p", newdata = .SD))), .SDcols = trms]
# newdata[, c("other") := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/ethnicity_grp_table.fst", 100L)

print("Table saved")

if (plots) {
 # only used for imputation of ~300 cases

}



