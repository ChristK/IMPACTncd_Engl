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

#Final distribution chosen: MN4

library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk
library(IMPACTncdEngl)
library(CKutils)
setDTthreads(5)

# Set some variables we will use later in automation
distributions <- FALSE
univariable_analysis <- FALSE
fit_model <- TRUE
diagnostics <- TRUE
plots <- TRUE
seed                <- 43L


# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)


ds <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90),
  .(smok_status, year, age, agegrp10, sex, qimd,
    ethnicity_grp, sha, urban_rural, wt_int)])
ds[, smok_status := as.integer(smok_status)]
ds[smok_status == 2, smok_status := 1]
# ds[, smok_status := factor(smok_status, 
#                               levels = c("1", "3", "4"))]
set.seed(seed)




distr_nam <- "MN4" #multinomial with 4 parameters 
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis.

if (univariable_analysis) {
  ds[, .(smok_status_median = wtd.quantile(smok_status, weight = wt_int)), keyby = .(age)
    ][, scatter.smooth(age, smok_status_median)]

  m_age1 <- gamlss(
    smok_status ~ pb(age),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")

  m_age2 <- gamlss(
    smok_status ~ poly(age, 3),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")

  m_age3 <- gamlss(
    smok_status ~ cs(age),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")
  GAIC(m_age1, m_age2, m_age3, k = log(nrow(ds))) # BIC m_age3 but not much more thatn m_age1
  GAIC(m_age1, m_age2, m_age3, k = 2) # AIC m_age1
  centiles(m_age1, xvar = ds$age)
  centiles(m_age2, xvar = ds$age)

  # ds[, .(smok_status_median = wtd.quantile(smok_status, weight = wt_int)), keyby = .(year)
  #   ][, scatter.smooth(year, smok_status_median, xlim = c(3, 40), ylim = c(0, 5))]
  # 
  m_year1 <- gamlss(
    smok_status ~ year,
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "blue1")

  m_year2 <- gamlss(
    smok_status ~ log(year),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "red1")

  m_year3 <- gamlss(
    smok_status ~ log(year - 2),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "green1")

  m_year4 <- gamlss(
    smok_status ~ log(year + 10),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "green1")
  
  GAIC(m_year1, m_year2, m_year3, m_year4, k = log(nrow(ds))) # BIC m_year1
  GAIC(m_year1, m_year2, m_year3, m_year4, k = 2) # AIC m_year1
  centiles(m_year1, xvar = ds$age)
  centiles(m_year2, xvar = ds$age)

}



if(fit_model){
# Setting a minimum model - we know we want these things
mod_min <- gamlss(
  smok_status ~ year + pb(age) + pcat(qimd) +
    sex  + pcat(ethnicity_grp) + pcat(sha) + urban_rural,
  family = distr_nam,
  weights = ds$wt_int,
  data = ds,
  method = mixed(5, 100),
  control = con1 # for speed but looses accuracy if not 1e-3
)
qsave(mod_min, "./inputs/exposure_distributions/AH_test/smok_status_mod_min.qs", preset = "high")
warnings() # check for warnings
mod_min <- qread( "./inputs/exposure_distributions/AH_test/smok_status_mod_min.qs")

# Stepwise model selection using a Generalized Akaike Information Criterion
smok_status_modelA <- stepGAICAll.A(
  mod_min,
  scope = list(
    lower = ~ year + pb(age) + pcat(qimd) +
      sex  + pcat(ethnicity_grp) + pcat(sha) + urban_rural,
    upper = ~ (year + pb(age) + pcat(qimd) +
                 sex  + pcat(ethnicity_grp) + pcat(sha) + urban_rural)^2
  ),
  sigma.scope = list(
    lower = ~1,
    upper = ~ year + pb(age) + pcat(qimd) +
      sex  + pcat(ethnicity_grp) + pcat(sha) + urban_rural
  ),
  nu.scope = list(
    lower = ~1,
    upper = ~ year + pb(age) + pcat(qimd) +
      sex  + pcat(ethnicity_grp) + pcat(sha) + urban_rural
  ),
  parallel = "multicore",
  ncpus = 5L,
  weights = wt_int,
  trace = TRUE  # This is so you can see all the models it has tried to fit
)
warnings()

smok_status_modelA

smok_status_modelA <- update(smok_status_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy

smok_status_modelA$call


qsave(smok_status_modelA, "./inputs/exposure_distributions/AH_test/smok_status_model.qs", preset = "high")
#smok_status_modelA <- qread("./inputs/exposure_distributions/AH_test/smok_status_modelAwt_July.qs")


GAIC.table(smok_status_modelA,
           #smok_status_modelB,
           mod_min)

smok_status_modelA$hsedata <- copy(ds)
qsave(smok_status_modelA,"./secure_data/lifecourse_models/smok_status_model.qs" )

}



smok_status_modelbest <- qread( "./inputs/exposure_distributions/AH_test/smok_status_model.qs")
# # Code to create a table with predictions
trms <- all.vars(formula(smok_status_modelbest))[-1] # -1 excludes dependent var
newdata <- CJ(
  year = 3:50,
  age = 20:90,
  sex = unique(ds$sex),
  qimd = unique(ds$qimd),
  sha = unique(ds$sha),
  ethnicity_grp = unique(ds$ethnicity_grp),
  urban_rural = unique(ds$urban_rural)
)
# newdata holds all the possible combinations of predictors

# This is to be able to parallelise
newdata <- split(newdata, by = "year")
newdata <- # assignment necessary! Copies of data.tables are happening
  lapply(
    newdata,
    function(x) {
      x[, (smok_status_modelbest$parameters) := predictAll(smok_status_modelbest, .SD, data = ds), .SDcols = trms]
    }
  )
newdata <- rbindlist(newdata) # bind the chunks back
View(head(newdata, 1e3))

write_fst(newdata, path = "./inputs/exposure_distributions/AH_test/smok_status_table.fst", compress = 100L) # This is what the model use as an input
# newdata <- read_fst("./inputs/exposure_distributions/AH_test/smok_status_table.fst", as.data.table = T)
# 
# print("Table saved")

if (diagnostics) {
  #smok_status_model <- qread("./smok_status_modelbest.qs")
  smok_status_model <- qread("./inputs/exposure_distributions/AH_test/smok_status_model.qs")
  
  plot(smok_status_model) 
  
  wp(smok_status_model) # detrended QQ-plot
  wp(resid = resid(smok_status_model)) # equivalen to the one above
  wp(resid = resid(smok_status_model), ylim.all = 80 * sqrt(1 / length(resid(smok_status_model))))
  tt <- wp(smok_status_model, xvar = ds$age, n.inter = 6)
  wp(smok_status_model, xvar = ds$year, n.inter = 10) 
  wp(resid = resid(smok_status_model), xvar = ~ ds$qimd)
  
  dtop(smok_status_model, xvar = ds$age)
  
  rqres.plot(smok_status_model)
  rqres.plot(smok_status_model, type = "QQ")
}

if (plots) {
  #smok_status_model <- qread("./smok_status_modelbest.qs")
  smok_status_model <- qread("./inputs/exposure_distributions/AH_test/smok_status_model.qs")
  dir.create("./inputs/exposure_distributions/AH_test/validation_synthpop_models", FALSE)
  
  #source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  
  
  xlab_nam <- expression(bold(smok_status))
  smok_status_model_tbl <- read_fst("./inputs/exposure_distributions/AH_test/smok_status_table.fst", as.data.table = TRUE)
  zz <-
    validate_gamlss_tbl(ds[age > 19], smok_status_model_tbl, 50, "smok_status", paste0("q", smok_status_model$family[1]))[
      between(smok_status, quantile(smok_status, 0.01), quantile(smok_status, 0.99))
    ]
  zz[, weight := wt_int / sum(wt_int), by = "type"]
  
  png(
    "./inputs/exposure_distributions/AH_test/validation_synthpop_models/smok_status_rel_dist.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", smok_status],
                      zz[type == "Modelled", smok_status],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      100
  )
  dev.off()
  
  plot_synthpop_val(zz, smok_status, "agegrp10", "wt_int", "smok_status by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, "year", "wt_int", "smok_status by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, "qimd", "wt_int", "smok_status by QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, "sha", "wt_int", "smok_status by SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, "ethnicity_grp", "wt_int", "smok_status by ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, "urban_rural", "wt_int", "smok_status by rurality", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, c("year", "agegrp10"), "wt_int", "smok_status by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, c("year", "qimd"), "wt_int", "smok_status by year and QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, c("year", "sha"), "wt_int", "BMI by year and SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, c("year", "ethnicity_grp"), "wt_int", "smok_status by year and ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, smok_status, c("year", "urban_rural"), "wt_int", "smok_status by year and smoking", xlab_nam, FALSE, FALSE)
}
