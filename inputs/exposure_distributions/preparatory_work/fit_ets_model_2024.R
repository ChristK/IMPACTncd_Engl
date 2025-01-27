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

# Final distribution chosen: BI 

print("ets")

library(IMPACTncdEngl)
library(data.table)
#library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk
library(future.apply)
library(CKutils)
source("./inputs/exposure_distributions/preparatory_work/fix_stepGAICall.A.R")


#options(future.fork.enable = FALSE)
plan("multicore")

# Set some variables we will use later in automation
distributions <- FALSE
univariable_analysis <- FALSE
fit_model <- FALSE
update <- FALSE
diagnostics <- FALSE
plots <- TRUE
seed                <- 43L





# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)

set.seed(seed)

ds <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90) & year != 7 & year != 8, #year 7 is half with and half without smokefree policy. Year 8 was influencial
                     .(ets, year, age, agegrp10, sex, qimd, ethnicity_grp, sha, wt_int, urban_rural, smok_status)])

ds[, ets := as.integer(as.character(ets))]
# Calculate smoking prevalence to be used as externality
tt1 <- ds[smok_status == "4",  sum(wt_int), by = .(year, sha)]
tt2 <- ds[,  sum(wt_int), by = .(year, sha)]
tt1[tt2, smok_prev := V1/i.V1, on = c("year", "sha")]
ds[tt1, smok_prev := i.smok_prev, on = c("year", "sha")]
rm(tt1, tt2)
ds[, legislation := year > 7]
set.seed(seed)


ds[, `:=`(
  smok_status = factor(smok_status), # 1 = never smoker; 2 = occasional; 3 = ex-; 4 = current
  sex = factor(sex),
  qimd = factor(qimd),
  sha = factor(sha), 
  urban_rural = factor(urban_rural), 
  ethnicity_grp = factor(ethnicity_grp) # 5 groups (more granular not available)
)]

# Treat occasional smokers as never smokers 
ds[smok_status == "2", smok_status := "1"]
summary(ds)

if(distributions){
marg_distr <- fitDist(ets,
                      k = 2, # Default is for AIC. use log(length(ds$ets)) for BIC
                      type = "binom", # Other options "realline", "realplus", "real0to1", "counts", "binom"
                      try.gamlss = TRUE, # Better but slow
                      extra = NULL,
                      data = ds,
                      trace = TRUE, 
                      weights = ds$wt_int
)

# Check the best 10 distr - look at their properties, how many parameters? 
head(marg_distr$fits)
# BI       BB     ZIBI      DBI     ZABI     ZIBB 
# 152296.4 152298.4 152298.4 152298.4 152298.4 152300.4 
distr_nam <- names(marg_distr$fits[1]) # BI

}

distr_nam <- "BI"
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.


if (univariable_analysis) {
  age_scaled <- scale(18:90, 50.8, 17.4)
  age <- c(18:90)
  ds[, .(ets_mean = wtd.mean(ets, weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, ets_mean)]

  m_age1 <- gamlss(
    ets ~ pb(age),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(age, mean_predictAll(m_age1, ds, data.frame("age" = age)), col = "blue1")

  m_age2 <- gamlss(
    ets ~ log(age+3),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(age, mean_predictAll(m_age2, ds, data.frame("age" = age)), col = "red1")

  m_age3 <- gamlss(
    ets ~ cs(age),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(age, mean_predictAll(m_age3, ds, data.frame("age" = age)), col = "green1")
  GAIC.table(m_age1, m_age2, m_age3, k = log(nrow(ds))) # BIC m_age3
  GAIC.table(m_age1, m_age2, m_age3, k = 2) # AIC m_age1

  ds[, .(ets_mean = wtd.mean(ets, weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, ets_mean, xlim = c(3, 40), ylim = c(0, 1))]

  m_year1 <- gamlss(
    ets ~ year * (year>7),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(3:40, mean_predictAll(m_year1, ds, data.frame("year" = 3:40)), col = "blue1")

  m_year2 <- gamlss(
    ets ~ log(year) * (year>7),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year2, ds, data.frame("year" = 3:40)), col = "red1")

  m_year3 <- gamlss(
    ets ~ log(year - 2),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year3, ds, data.frame("year" = 3:40)), col = "green1")

  GAIC.table(m_year1, m_year2, m_year3, k = log(nrow(ds))) # BIC m_year1
  GAIC.table(m_year1, m_year2, m_year3, k = 2) # AIC m_year1
  centiles(m_year1, xvar = ds$age)
  centiles(m_year2, xvar = ds$age)
  
  
  
  # Smoking status 
  m_sm1 <- gamlss(
    ets ~ smok_status,
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  summary(m_sm1)
  
  # Using a smoother 
  m_sm2 <- gamlss(
    ets ~ pcat(smok_status, method = "GAIC"),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  summary(m_sm2)
  
  GAIC.table(m_sm1, m_sm2) #m_sm1 is best
  
  #Looking at the levels pcat has reduced smok_status to, and how they 
  #correspond to the original levels - 
  f2 <- getSmo(m_sm2)$factor
  levels(f2)
  table(f2,ds$smok_status)

}

if(fit_model){



# Setting a minimum model - we know we want these things


mod_min <- gamlss(
  ets ~  year * legislation + pb(age) + pcat(qimd),#+ pb(age) + pcat(qimd) +
   # sex + smok_status + smok_prev, # + pcat(ethnicity_grp) + pcat(sha), #+ urban_rural + smok_prev,
  family = distr_nam,
  weights = ds$wt_int,
  data = ds,
  method = mixed(5, 100),
  control = con1 # for speed but looses accuracy if not 1e-3
)
qsave(mod_min, "./inputs/exposure_distributions/AH_test/ets_mod_min.qs", preset = "high")
warnings() # check for warnings
mod_min <- qread( "./inputs/exposure_distributions/AH_test/ets_mod_min.qs")

# Stepwise model selection using a Generalized Akaike Information Criterion

ets_modelA <- stepGAICAll.A(
  mod_min,
  scope = list(
    lower = ~ year * legislation + pb(age) + pcat(qimd)   ,
    upper = ~ (year * legislation + pb(age) + pcat(qimd) +
     sex + smok_status  + pcat(ethnicity_grp) + pcat(sha) +  urban_rural + smok_prev)^3
  ),
  parallel = "multicore",
  ncpus = 12L, #,
  weights = wt_int,
 # trace = F  # This is so you can see all the models it has tried to fit
)
warnings()

ets_modelA

ets_modelA <- update(ets_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy

ets_modelA$call
# gamlss(formula = ets ~ year + legislation + pb(age) + pcat(qimd) +
# smok_status + sex + pcat(ethnicity_grp) + pcat(sha) + urban_rural +
#   year:legislation + legislation:smok_status + legislation:sex +
#   pb(age):sex + year:pb(age) + pb(age):urban_rural + smok_status:urban_rural +
#   legislation:urban_rural + sex:urban_rural + pb(age):smok_status +
#   smok_status:sex + legislation:smok_status:urban_rural + legislation:smok_status:sex +
#   pb(age):smok_status:sex, family = distr_nam, data = ds, method = mixed(5,
#                                                                          100), control = gamlss.control(c.crit = 0.001), trace = FALSE)

qsave(ets_modelA, "./inputs/exposure_distributions/AH_test/ets_modelA.qs", preset = "high")

}


if(update){
ets_model <- qread("./inputs/exposure_distributions/AH_test/ets_modelA.qs")

# tt <- chooseDist(ets_model,
#                  type = "binom",
#                  trace = TRUE, data = ds,
#                  parallel = "multicore", ncpus = 15L
# )
# 
# which(tt == min(tt[, 1], na.rm = TRUE), arr.ind = TRUE) # Best distr based on AIC
# which(tt == min(tt[, 3], na.rm = TRUE), arr.ind = TRUE) # Best distr based on BIC
# BI is the best for both 

# 
# 
# 
ets_model$hsedata <- copy(ds)
# 
qsave(ets_model, "./secure_data/lifecourse_models/ets_model.qs", preset = "high")
# 
# print("Model saved")
# 
trms <- all.vars(formula(ets_model))[-1] # -1 excludes dependent var
newdata <- CJ(
  year = 3:50,
  age = 20:90,
  sex = unique(ds$sex),
  qimd = unique(ds$qimd),
  sha = unique(ds$sha),
  ethnicity_grp = unique(ds$ethnicity_grp),
  urban_rural = unique(ds$urban_rural), 
  smok_status = unique(ds$smok_status)
)
newdata[, legislation := year > 7]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(ets_model, .SD, data = ds), .SDcols = trms], future.seed = TRUE)
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)


kc <- sort(setdiff(names(newdata), c("mu")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/AH_test/ets_table.fst", 100L)

print("Table saved")

}
# Diagnostics -------------------------------------------------------------
if (diagnostics) {
  ets_model <- qread("./inputs/exposure_distributions/AH_test/ets_modelA.qs")
  ets_model$data <- NULL
  wp(ets_model)
  #wp(ets_model, xvar = ~ ds$age) this doesn't work 

  plot(ets_model)
}

# Plots -------------------------------------------------------------
if (plots) {
  xlab_nam <- expression(bold(ETS))
  ets_model_tbl <- read_fst("./inputs/exposure_distributions/AH_test/ets_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(ds[age >= 20L], ets_model_tbl, 50, "ets", paste0("q", distr_nam))

  # [between(ets, quantile(ets, 0.01), quantile(ets, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_int / sum(wt_int), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/ETS_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", ets],
                        zz[type == "Modelled", ets],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        2/3, TRUE)
    dev.off()
  })
  future(plot_synthpop_val(zz, ets, "agegrp10", "wt_int",
                           "ETS by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, "year", "wt_int",
                           "ETS by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, "qimd", "wt_int",
                           "ETS by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, "sha", "wt_int",
                           "ETS by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, "ethnicity_grp", "wt_int",
                           "ETS by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, "urban_rural", "wt_int",
                           "ETS by rurality", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("year", "agegrp10"), "wt_int",
                           "ETS by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("year", "qimd"), "wt_int",
                           "ETS by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("year", "sha"), "wt_int",
                           "ETS by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("year", "ethnicity_grp"), "wt_int",
                           "ETS by year and ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("year", "urban_rural"), "wt_int",
                           "ETS by year and rurality", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("smok_status"), "wt_int",
                           "ETS by smoking status", xlab_nam, FALSE, FALSE))

}

