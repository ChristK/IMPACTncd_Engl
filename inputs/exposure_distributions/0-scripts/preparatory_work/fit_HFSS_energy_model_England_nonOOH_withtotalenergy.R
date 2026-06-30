## workHORSE is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - workHORSE: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## workHORSE is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

rm(list=ls())
gc()

library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk
library(haven) #have to install this... as the data were from STATA

#setDTthreads(1L, restore_after_fork = FALSE)
#threads_fst(1L, reset_after_fork = FALSE)

# Set some variables we will use later in automation
diagnostics <- TRUE
plots <- TRUE

setwd("/mnt/storage_fast4/Salt_England/stata13_se")
ds_ <- qread("HFSS_personlevel.qs")
# Load the dataset #15655 obs

#Here we exclude potential under-reporting, using Black Method, performed by Zoe Colombet,
#file from Zoe
under_report <- readRDS("underreporters.rds") #this is for the UK, and our current data for England >= 20 years
table(under_report$supp_scho_155)
ds_$seriali <- as.character(ds_$seriali)
ds_ <- ds_[under_report, on= "seriali"]

#Adding energy total to the prediction
energy_total <- readRDS("energy_diet_consump.rds")
setDT(energy_total)
energy_total <- energy_total [, .(seriali, energy_total = Energykcal)]
energy_total$seriali <- as.character(energy_total$seriali)
ds_ <- ds_[energy_total, on= "seriali"]

test <- ds_[, diff:= (Energykcal-energy_total)]
summary(test$diff) 
test <- test[diff>0,] #no energy HFSS > total_energy

#prop energy from HFSS
test <- ds_[, prop:= (Energykcal/energy_total)]
mean(ds_$prop)
weighted.mean(ds_$prop, w = ds_$s_weight_diary) #46%

ds <- ds_[Age>=20 & country==101 & supp_scho_155==0, .(  # will not include those < 20 years old
  seriali,
  age = as.numeric(Age),
  age_group=ifelse(Age >= 20 & Age <= 29, "20-29",
                   ifelse(Age >= 30 & Age <= 39, "30-39",
                   ifelse(Age >= 40 & Age <= 49, "40-49",
                   ifelse(Age >= 50 & Age <= 59, "50-59",
                   ifelse(Age >= 60 & Age <= 69, "60-69",
                   ifelse(Age >= 70 & Age <= 79, "70-79",
                   ifelse(Age >= 80, "80+", NA))))))),
  sex = as.factor(Sex), #1 Men 2 Women
  SES = as.factor(IMD), #1 least deprived, 5 most deprived
  bmi = bmi,
  bmi_group = ifelse(bmi <18.5, "<18.5", 
                     ifelse (bmi >=18.5 & bmi <25, "18.5- <25",
                     ifelse (bmi >=25 & bmi <30, "25- < 30",
                     ifelse (bmi >=30, ">=30",NA)))),
  energy = Energykcal,
  #food_weight = TotalGrams,
  energy_total = energy_total,
  energy_total_group = ifelse(energy_total <=1500, "<=1500", ">1500"),
  year= as.integer(surveyyr+8), # survey years 1 - 11, representing 2009 - 2019
  weight= s_weight_diary,
  #ethnic = as.factor(ethnic_5), #1 White, 2 Mixed ethnic group, 3 Black or Black British, 4 Asian or asian British, 5 Other
  #ethnic2 = as.factor(ifelse(ethnic_5==5, 2, ethnic_5)), #for IMPACT-NCD model, Mixed and Other were combined together
  #gor= as.factor(gor), #1 North East, 2 North West, 3 Yorkshire and the Humber, 4 East Midlands, 5 West Midlands, 6 East of England, 7 London, 8 South East, 9 South West
  OOH= OOH
)]

summary(ds)
ds[, energy:= energy*(1-OOH)] #here we estimated energy from non-OOH
weighted.mean(ds$energy, ds$weight)
#823.6717

ds <- ds[complete.cases(ds), ]
sum(is.na(ds)) # 3925 obs
#View(ds)

ds <- ds[year > 12, ] #selecting year from 2013 to 2019 because of this is predicted by GAMLSS MTL #2539 obs
summary(ds)

# First identify what distribution is best fit for energy consumption (the left hand side variable)
#dev.off()
plot(density(ds$energy))
summary(ds$energy)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 34.48  553.34  752.78  793.88  994.52 2894.25

# keep energy < 4000, it is likely overestimate
# ds <- ds[energy<4000,] # 3923 obs
# plot(density(ds$energy))
# summary(ds$energy)

# Let's see what distribution work best for this. We will use the function
# fitDist from gamlss. Please have a read of the documentation for this function
marg_distr <- fitDist(energy,
                      k = 2, # Default is for AIC. use log(length(ds$energy)) for BIC
                      type = "realAll", # Other options "realline", "realplus", "real0to1", "counts", "binom"
                      #note: real0to1 if var between 0 and 1
                      #counts: for count variables
                      #binom: for binomial with only 0 and 1
                      try.gamlss = TRUE, # Better but slow
                      extra = NULL,
                      data = ds,
                      trace = TRUE
)

# Check the best 10 distr
head(marg_distr$fits, 10)
#for the current data, 10 best distributions
# BCTo      BCT    BCPEo     BCPE      ST5     JSUo      JSU      GB2     SEP4     BCCG 
# 36650.14 36650.14 36651.76 36651.76 36653.61 36654.17 36654.17 36654.22 36657.36 36660.67

#class(marg_distr)
#str(marg_distr)
#saveRDS(marg_distr, file="marg_distr.rds") #it is better to save the object as the R studio is sometime error
#qsave(marg_distr, "marg_distr.qs", preset = "high")
#View(marg_distr)
#marg_distr <- qread("marg_distr.qs")

# lets plot some of the best three distributions
m1 <- gamlss(energy ~ 1,                          
             family = names(marg_distr$fits[1]), 
             # weights = ds$s_weight,
             method = mixed(20, 40),
             data = ds
)
#saveRDS(m1, file="m1.rds")
#qsave(m1, "m1.qs", preset = "high")
#m1 <- qread("m1.qs")


# Let's create a helper function
sample_from_gamlss <- function(m, sample_size = 1e4) {
  # m is a gamlss model
  stopifnot(inherits(m, "gamlss")) # Stop if model is not gamlss class
  stopifnot(unique(sapply(
    lapply(
      m$parameters,
      function(x) {
        as.numeric(coef(m, x))
      }
    ), length
  )) == 1L) # stop if predictors are present
  parm <- lapply(m$parameters, function(x) {
    fitted(m, x)[1]
  })
  names(parm) <- m$parameters
  y <-
    do.call(paste0("r", m$family[1]), c("n" = sample_size, parm)) # Sample from the model
}

plot(density(ds$energy))
lines(density(sample_from_gamlss(m1)), col = alpha("red", 0.4))

m2 <- gamlss(energy ~ 1, family = names(marg_distr$fits[2]), data = ds, method = mixed(20, 40))
#saveRDS(m2, file="m2.rds")  ### Some warnings for SEP1 ###Warning messages:
                             ### 1: In log(0.5 + w * exp((1 - (1/t1)) * log(t1) - lgamma(1/t1) -  ... : NaNs produced
#qsave(m2, "m2.qs", preset = "high")
#m2 <- qread("m2.qs")
lines(density(sample_from_gamlss(m2)), col = alpha("blue", 0.4))

m3 <- gamlss(energy ~ 1, family = names(marg_distr$fits[3]), data = ds, method = mixed(20, 40))
#saveRDS(m3, file="m3.rds")
#qsave(m3, "m3.qs", preset = "high")
#m3 <- qread("m3.qs")
lines(density(sample_from_gamlss(m3)), col = alpha("green", 0.4))


# An alternative plot
#reldist(sample_from_gamlss(m1), ds$sbp, method = "bgk", bar = TRUE, show = "effect")

# What would the plot looked like if the distribution was off
#reldist(sample_from_gamlss(m1) + 1, ds$sbp, method = "bgk", bar = TRUE, show = "effect")
#reldist(sample_from_gamlss(m1) - 1, ds$sbp, method = "bgk", bar = TRUE, show = "effect")

#Let's look at the AIC for these 3 models (the lower the AIC, the best)
GAIC.table(m1, m2, m3)
# minimum GAIC(k= 2 ) model: m2 
# minimum GAIC(k= 3.84 ) model: m2 
# minimum GAIC(k= 7.84 ) model: m2 
# df      k=2   k=3.84   k=7.84
# m1  4 36650.14 36657.50 36673.50
# m2  4 36650.14 36657.50 36673.50
# m3  4 36651.76 36659.12 36675.12

distr_nam <- names(marg_distr$fits[2]) # better fit based on the plot and GAIC, BCT
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3

distr_nam

# Second, identify what distribution is best fit for predictors (the right hand side variable)
# Predictors
# Age
ds[, .(energy_median = wtd.quantile(energy)), keyby = .(age)][, scatter.smooth(age, energy_median)]
#ds[, .(energy_median = wtd.quantile(energy, weight = s_weight)), keyby = .(age)][, scatter.smooth(age, energy_median)]

#linear
m_age0 <- gamlss(
  energy ~ age,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20) #few iteration to have a quick results
)
#saveRDS(m_age0, file="m_age0.rds")
#qsave(m_age0, "m_age0.qs", preset = "high")
lines(centiles.pred(m_age0, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "purple")

#penalised B-spline
m_age1 <- gamlss(
  energy ~ pb(age),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(50, 100)
)
#saveRDS(m_age1, file="m_age1.rds")
#qsave(m_age1, "m_age1.qs", preset = "high")
lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")

#polynomial
m_age2 <- gamlss(
  energy ~ poly(age, 3),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(50, 100)
)
#saveRDS(m_age2, file="m_age2.rds")
#qsave(m_age2, "m_age2.qs", preset = "high")
lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")

#cubic-spline
#m_age3 <- gamlss(
#  energy ~ cs(age),
#  family = distr_nam,
#  # weights = ds$s_weight,
#  data = ds,
#  method = mixed(20, 20)
#)
#lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")

GAIC.table(m_age0, m_age1, m_age2)
# minimum GAIC(k= 2 ) model: m_age1 
# minimum GAIC(k= 3.84 ) model: m_age1 
# minimum GAIC(k= 7.84 ) model: m_age0 
# df      k=2   k=3.84   k=7.84
# m_age0 5.000000 36615.61 36624.81 36644.81
# m_age1 6.397918 36611.22 36622.99 36648.58
# m_age2 7.000000 36612.87 36625.75 36653.75 

#centiles(m_age1, xvar = ds$age)
#centiles(m_age2, xvar = ds$age)


# Year
ds[, .(energy_median = wtd.quantile(energy)), keyby = .(year)][, scatter.smooth(year, energy_median, xlim = c(min(ds$year),max(ds$year)))]
#ds_summary <- ds[, .(energy_median = wtd.quantile(energy)), keyby = .(year)]
#plot(ds_summary$year, ds_summary$energy_median, main = "energy Median vs. Year")

m_year1 <- gamlss(
  energy ~ year,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)
#qsave(m_year1, "m_year1.qs", preset = "high")
lines(centiles.pred(m_year1, xname = "year", xvalues = min(ds$year):max(ds$year), cent = 50, data = ds), col = "blue1")

m_year2 <- gamlss(
  energy ~ log(year),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)
#qsave(m_year2, "m_year2.qs", preset = "high")
lines(centiles.pred(m_year2, xname = "year", xvalues = min(ds$year):max(ds$year), cent = 50, data = ds), col = "red1")

m_year3 <- gamlss(
  energy ~ log(year + 100),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)
#qsave(m_year3, "m_year3.qs", preset = "high")
lines(centiles.pred(m_year3, xname = "year", xvalues = min(ds$year):max(ds$year), cent = 50, data = ds), col = "green1")

GAIC.table(m_year1, m_year2, m_year3)
# minimum GAIC(k= 2 ) model: m_year3 
# minimum GAIC(k= 3.84 ) model: m_year3 
# minimum GAIC(k= 7.84 ) model: m_year3 
# df      k=2   k=3.84   k=7.84
# m_year1  5 36647.35 36656.55 36676.55
# m_year2  5 36647.39 36656.59 36676.59
# m_year3  5 36647.35 36656.55 36676.55
#centiles(m_year1, xvar = ds$year)
#centiles(m_year2, xvar = ds$year)
# based on GAIC.table, is m_year3 (log + 100), but we select year1 (linear)

# bmi
ds[, .(energy_median = wtd.quantile(energy)), keyby = .(bmi)][, scatter.smooth(bmi, energy_median)]

m_bmi1 <- gamlss(
  energy ~ bmi,
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50,20) #few iteration to have a quick results
)
lines(centiles.pred(m_bmi1, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "purple")

#for this one, using log, be carefull that your bmi are not 0 because log will not work!!!
m_bmi2 <- gamlss(
  energy ~ log(bmi),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20)
)
lines(centiles.pred(m_bmi2, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "red1")

m_bmi3 <- gamlss(
  energy ~ log(bmi+100),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20)
)
lines(centiles.pred(m_bmi3, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "green1")

m_bmi4 <- gamlss(
  energy ~ poly(bmi, 3),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20)
)
lines(centiles.pred(m_bmi4, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "blue")

m_bmi5 <- gamlss(
  energy ~ pb(bmi),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20) 
)
lines(centiles.pred(m_bmi5, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "yellow")

GAIC.table(m_bmi1, m_bmi2, m_bmi3, m_bmi4, m_bmi5)
# minimum GAIC(k= 2 ) model: m_bmi5 
# minimum GAIC(k= 3.84 ) model: m_bmi5 
# minimum GAIC(k= 7.84 ) model: m_bmi2 
# df      k=2   k=3.84   k=7.84
# m_bmi1 5.000000 36648.95 36658.15 36678.15
# m_bmi2 5.000000 36647.59 36656.79 36676.79
# m_bmi3 5.000000 36648.67 36657.87 36677.87
# m_bmi4 7.000000 36644.96 36657.84 36685.84
# m_bmi5 7.732506 36641.37 36655.60 36686.53
# based on GAIC.table, here we will pick m_bmi5 (pb) 

# energy_total
ds[, .(energy_median = wtd.quantile(energy)), keyby = .(energy_total)][, scatter.smooth(energy_total, energy_median)]

m_energy_total1 <- gamlss(
  energy ~ energy_total,
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50,20) #few iteration to have a quick results
)
lines(centiles.pred(m_energy_total1, xname = "energy_total", xvalues = min(ds$energy_total):max(ds$energy_total), cent = 50, data = ds), col = "purple")

#for this one, using log, be carefull that your bmi are not 0 because log will not work!!!
m_energy_total2 <- gamlss(
  energy ~ log(energy_total),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20)
)
lines(centiles.pred(m_energy_total2, xname = "energy_total", xvalues = min(ds$energy_total):max(ds$energy_total), cent = 50, data = ds), col = "red1")

m_energy_total3 <- gamlss(
  energy ~ log(energy_total+100),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20)
)
lines(centiles.pred(m_energy_total3, xname = "energy_total", xvalues = min(ds$energy_total):max(ds$energy_total), cent = 50, data = ds), col = "green1")

m_energy_total4 <- gamlss(
  energy ~ poly(energy_total, 3),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20)
)
lines(centiles.pred(m_energy_total4, xname = "energy_total", xvalues = min(ds$energy_total):max(ds$energy_total), cent = 50, data = ds), col = "blue")

m_energy_total5 <- gamlss(
  energy ~ pb(energy_total),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20) 
)
lines(centiles.pred(m_energy_total5, xname = "energy_total", xvalues = min(ds$energy_total):max(ds$energy_total), cent = 50, data = ds), col = "yellow")

GAIC.table(m_energy_total1, m_energy_total2, m_energy_total3, m_energy_total4, m_energy_total5)
# minimum GAIC(k= 2 ) model: m_energy_total5 
# minimum GAIC(k= 3.84 ) model: m_energy_total5 
# minimum GAIC(k= 7.84 ) model: m_energy_total5 
# df      k=2   k=3.84   k=7.84
# m_energy_total1  5 34631.01 34640.21 34660.21
# m_energy_total2  5 34828.98 34838.18 34858.18
# m_energy_total3  5 34804.97 34814.17 34834.17
# m_energy_total4  7 34632.20 34645.08 34673.08
# m_energy_total5  5 34631.01 34640.21 34660.21
# based on GAIC.table, here we will pick m_energy_total5 (pb) 

#for factors, not much to do. only one smoother
#SES
m_ses1 <- gamlss(
  energy ~ SES,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

m_ses2 <- gamlss(
  energy ~ pcat(SES), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

GAIC.table(m_ses1, m_ses2)


#Sex
m_sex1 <- gamlss(
  energy ~ sex,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

m_sex2 <- gamlss(
  energy ~ pcat(sex), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

GAIC.table(m_sex1, m_sex2)


# #Ethnicity
# m_eth1 <- gamlss(
#   energy ~ ethnic2,
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# m_eth2 <- gamlss(
#   energy ~ pcat(ethnic2), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# GAIC.table(m_eth1, m_eth2)
# 
# 
# #GOR
# m_gor1 <- gamlss(
#   energy ~ gor,
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# m_gor2 <- gamlss(
#   energy ~ pcat(gor), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# GAIC.table(m_gor1, m_gor2)


##### 
# Theory driven model selection check ->?BCT, there are FOUR parameters to be included, the FOURTH will not be used.
# This takes a long time to run
energy_model <- gamlss(
  energy ~ year + poly(age,3) + pb(bmi) + pb(energy_total) + sex + SES, #poly is used here for age for a better fit and reduce computation
  ~ year + poly(age,3) + pb(bmi) + pb(energy_total) + sex + SES, #interaction is not added due to negative mu s and limited data
  ~ year + poly(age,3) + pb(energy_total),
  ~ year + poly(age,3) + pb(energy_total),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 500),
  control = gamlss.control(c.crit = 1e-3)
)

qsave(energy_model, "HFSS_energy_model_BCT_nonOOH_byenergytotal.qs", preset = "high")
energy_model <- qread("HFSS_energy_model_BCT_nonOOH_byenergytotal.qs")

energy_model2 <- update(energy_model, family="BCTo") #we update using BCTo because BCT resulted in some negative mu's

GAIC.table(energy_model, energy_model2) #model 2 is better
# minimum GAIC(k= 2 ) model: energy_model2 
# minimum GAIC(k= 3.84 ) model: energy_model 
# minimum GAIC(k= 7.84 ) model: energy_model 
# df      k=2   k=3.84   k=7.84
# energy_model  42.41505 34476.49 34554.53 34724.19
# energy_model2 46.04365 34471.84 34556.56 34740.74

qsave(energy_model2, "HFSS_energy_model_BCTo_nonOOH_byenergytotal.qs", preset = "high")
model_final <- qread("HFSS_energy_model_BCTo_nonOOH_byenergytotal.qs")

# Code to create a table with predictions
trms <- all.vars(formula(model_final))[-1] # -1 excludes dependent var

#creating newdata which holds all the possible combinations of predictors // you can leave it like that, it's fine
newdata <- CJ(
  year = 3:50, 
  age = 20:99, 
  sex = unique(factor(ds$sex)), 
  SES = unique(ds$SES),
  #policy = unique(factor(ds$policy)),
  #gor = unique(ds$gor),
  #ethnic2 = unique(factor(ds$ethnic2)),
  bmi = seq(14L,70L, 1L), #min(ds$bmi) #integer is needed for IMPACT-NCD fitting #I just created until 100 which later can be adjusted to our needs.
  energy_total = seq(0L, 5000L, 100L)
)

# This is to be able to parallelise
newdata <- split(newdata, by = "year")
newdata <- # assignment necessary! Copies of data.tables are happening
  lapply(
    newdata,
    function(x) {
      x[, (model_final$parameters) := predictAll(model_final, .SD, data = ds), .SDcols = trms]
    }
  )
newdata <- rbindlist(newdata) # bind the chunks back
#View(head(newdata, 1e3))
summary(newdata) #BCT creates negative mu, lets try with BCTo

write_fst(newdata, path = "./HFSS_energy_table_BCTo_nonOOH_byenergytotal.fst", compress = 100L) #I don't save this, the final table was adjusted for its use in IMPACT-NCD
qsave(model_final, "./HFSS_energy_modelfinal_BCTo_nonOOH_byenergytotal.qs", preset = "high")

print("Table saved")


###LETS SEE THE DIAGNOSTICS
#IF YOU CAN SAVE THE PLOT IN A FILE, IT WILL BE GREAT

if (diagnostics) {
  energy_model_final <- qread("./HFSS_energy_modelfinal_BCTo_nonOOH_byenergytotal.qs")
  
  plot(energy_model_final)
  
  wp(energy_model_final) # detrended QQ-plot
  wp(resid = resid(energy_model_final)) # equivalent to the one above
  wp(resid = resid(energy_model_final), ylim.all = 1000 * sqrt(1 / length(resid(energy_model_final))))
  tt <- wp(energy_model_final, xvar = ds$age, n.inter = 6)
  wp(energy_model_final, xvar = ds$year, n.inter = 10)
  wp(resid = resid(energy_model_final), xvar = ~ ds$SES)
  
  dtop(energy_model_final, xvar = ds$age) # here we want the blue line to be horizontal
  
  rqres.plot(energy_model_final)
  rqres.plot(energy_model_final, type = "QQ")
}

if (plots) {
  energy_model_final <- qread("./HFSS_energy_modelfinal_BCTo_nonOOH_byenergytotal.qs")
  dir.create("./validation/synthpop_models", recursive = TRUE)
  
  source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  
  
  xlab_nam <- expression(bold(energy ~ (kcal)))
  energy_model_final_tbl <- read_fst("./HFSS_energy_table_BCTo_nonOOH_byenergytotal.fst", as.data.table = TRUE)
  ds$bmisaved <- ds$bmi 
  ds$energy_total_saved <- ds$energy_total
  
  #I had issue with merging because of the BMI, so let's do a trick
  ds$bmi <- as.integer(round(ds$bmi,0))
  ds$energy_total <- as.integer(round(ds$energy_total / 100) * 100)

  #energy_model_final_tbl$bmi <- energy_model_final_tbl$bmi
  #energy_model_final_tbl[, bmi := as.integer(bmi)] 
  
  zz <-
    validate_gamlss_tbl(ds[age > 19 & age<=90 & bmi>=14 & bmi<=70 & energy_total>=0 & energy_total<=5000], energy_model_final_tbl, 50, "energy", paste0("q", energy_model_final$family[1]))[
      between(energy, quantile(energy, 0.01), quantile(energy, 0.99))
    ]
  
  #summary(zz$bmi)
  zz[, weight := weight / sum(weight), by = "type"]
  #summary(zz)
  #summary(energy_model_final_tbl)
  png(
    "./HFSS_energy_rel_dist_HFSS_BCTo_nonOOH_byenergytotal.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", energy],
                      zz[type == "Modelled", energy],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      100
  )
  dev.off()
  
  plot_synthpop_val(zz, energy, "age_group", "weight", "energy by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, energy, "sex", "weight", "energy by sex", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, energy, "year", "weight", "energy by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, energy, "SES", "weight", "energy by SES", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, energy, "gor", "weight", "energy by region", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, energy, "ethnic2", "weight", "energy by ethnic", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, energy, "bmi_group", "weight", "energy by bmi", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, energy, "energy_total_group", "weight", "energy by energy total", xlab_nam, FALSE, FALSE)
  
  plot_synthpop_val(zz, energy, c("year", "age_group"), "weight", "energy by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, energy, c("year", "sex"), "weight", "energy by year and sex", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, energy, c("year", "SES"), "weight", "energy by year and SES", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, energy, c("year", "gor"), "weight", "energy by year and region", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, energy, c("year", "ethnic2"), "weight", "energy by year and ethnic", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, energy, c("year", "bmi_group"), "weight", "energy by year and bmi", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, energy, c("year", "energy_total_group"), "weight", "energy by year and energy total", xlab_nam, FALSE, FALSE)
}

### Editing the data to merge with IMPACTncd model

#svy <- read_fst("/mnt/storage_fast4/Salt_England/IMPACTncd_Engl/synthpop_20/hf_real/synthpop_99ba702352a4ba3f55419baabacfb658_1.fst")
#summary(svy$ethnicity) #in the synthetic population, there are 9 ethic groups:
                        #white, indian, pakistani, bangladeshi, other asian, black caribbean, black african, chinese, other
                        #in the IMPACT-NCD, mixed and other were combined
#summary(svy$sha) #in the synthetic population, North East, North West, Yorkshire and the Humber,
                  #East Midlands, West Midlands, East of England, London, South East Coast,
                  #South Central, South West

energy_model_final_tbl <- read_fst("./HFSS_energy_table_BCTo_nonOOH_byenergytotal.fst", as.data.table = TRUE)
HFSS_energy_table_BCTo_nonOOH_byenergytotal <- energy_model_final_tbl  [, .(year=year,
                                                                    age=age,
                                                                    sex=as.factor(fifelse(sex==1, "men", 
                                                                                  fifelse(sex==2, "women", NA))),
                                                                    qimd=as.factor(fifelse(SES==1, "1 most deprived",
                                                                                   fifelse(SES==2, "2",
                                                                                   fifelse(SES==3, "3",
                                                                                   fifelse(SES==4, "4",
                                                                                   fifelse(SES==5, "5 least deprived", NA)))))),
                                                                    bmi = bmi,
                                                                    energy_total = energy_total,
                                                                    mu=mu,
                                                                    sigma=sigma,
                                                                    nu=nu,
                                                                    tau=tau)]
summary(HFSS_energy_table_BCTo_nonOOH_byenergytotal)
write_fst(HFSS_energy_table_BCTo_nonOOH_byenergytotal, path = "./HFSS_energy_table_final_BCTo_nonOOH_byenergytotal.fst", compress = 100L) # This is what my models use as an input

print("Table saved")
