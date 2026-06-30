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

ds <- ds_[Age>=20 & country==101 & supp_scho_155==0, .(  # will not include those < 20 years old
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
  food_weight = TotalGrams,  
  energy = Energykcal, 
  energy_group = ifelse(Energykcal <1500, "<1500", 
                 ifelse(Energykcal >=1500, ">=1500",NA)), 
  year= as.integer(surveyyr+8), # survey years 1 - 11, representing 2009 - 2019
  weight= s_weight_diary,
  #ethnic = as.factor(ethnic_5), #1 White, 2 Mixed ethnic group, 3 Black or Black British, 4 Asian or asian British, 5 Other
  #ethnic2 = as.factor(ifelse(ethnic_5==5, 2, ethnic_5)), #for IMPACT-NCD model, Mixed and Other were combined together
  #gor= as.factor(gor), #1 North East, 2 North West, 3 Yorkshire and the Humber, 4 East Midlands, 5 West Midlands, 6 East of England, 7 London, 8 South East, 9 South West
  OOH= OOH
)]

summary(ds)
ds[, food_weight:= food_weight*(1-OOH)] #here we estimated energy from non-OOH
ds[, energy:= energy*(1-OOH)]
weighted.mean(ds$food_weight, ds$weight)
#336.0511
weighted.mean(ds$energy, ds$weight)
#823.6717

ds <- ds[complete.cases(ds), ]
sum(is.na(ds)) # 3925 obs
#View(ds)

#ds <- ds[year > 12, ] #selecting year from 2013 to 2019 because of this is predicted by GAMLSS MTL #2539 obs
summary(ds) #error in fitting the model if data are few

# # creating a variable in which before and after the introduction of TLS
# ds [, policy:= as.factor(ifelse(year<13, 1, 2))]
# summary(ds)

# First identify what distribution is best fit for food_weight consumption (the left hand side variable)
#dev.off()
plot(density(ds$food_weight))
summary(ds$food_weight)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 6.429  204.906  295.861  328.056  412.560 1987.325


# Let's see what distribution work best for this. We will use the function
# fitDist from gamlss. Please have a read of the documentation for this function
marg_distr <- fitDist(food_weight,
                      k = 2, # Default is for AIC. use log(length(ds$food_weight)) for BIC
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
# BCT     BCTo      ST5     JSUo      JSU      GB2     BCPE    BCPEo     SEP4 
# 50661.12 50661.12 50663.65 50667.70 50667.70 50667.99 50673.23 50673.23 50679.81 
# ST1 
# 50683.91 

#class(marg_distr)
#str(marg_distr)
#saveRDS(marg_distr, file="marg_distr.rds") #it is better to save the object as the R studio is sometime error
#qsave(marg_distr, "marg_distr.qs", preset = "high")
#View(marg_distr)
#marg_distr <- qread("marg_distr.qs")

# lets plot some of the best three distributions
m1 <- gamlss(food_weight ~ 1,                          
             family = names(marg_distr$fits[1]), 
             # weights = ds$s_weight,
             method = mixed(50, 100),
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

plot(density(ds$food_weight))
lines(density(sample_from_gamlss(m1)), col = alpha("red", 0.4))

m2 <- gamlss(food_weight ~ 1, family = names(marg_distr$fits[2]), data = ds, method = mixed(20, 40))
#saveRDS(m2, file="m2.rds")  ### Some warnings for SEP1 ###Warning messages:
                             ### 1: In log(0.5 + w * exp((1 - (1/t1)) * log(t1) - lgamma(1/t1) -  ... : NaNs produced
#qsave(m2, "m2.qs", preset = "high")
#m2 <- qread("m2.qs")
lines(density(sample_from_gamlss(m2)), col = alpha("blue", 0.4))

m3 <- gamlss(food_weight ~ 1, family = names(marg_distr$fits[3]), data = ds, method = mixed(20, 40))
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
# minimum GAIC(k= 8.28 ) model: m2 
# df      k=2   k=3.84   k=8.28
# m1  4 50661.12 50668.48 50686.24
# m2  4 50661.12 50668.48 50686.24
# m3  4 50667.37 50674.73 50692.49

distr_nam <- names(marg_distr$fits[2]) # better fit based on the plot, BCTo
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3

distr_nam

# Second, identify what distribution is best fit for predictors (the right hand side variable)
# Predictors
# Age
ds[, .(food_weight_median = wtd.quantile(food_weight)), keyby = .(age)][, scatter.smooth(age, food_weight_median)]
#ds[, .(food_weight_median = wtd.quantile(food_weight, weight = s_weight)), keyby = .(age)][, scatter.smooth(age, food_weight_median)]

#linear
m_age0 <- gamlss(
  food_weight ~ age,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(50, 100) #few iteration to have a quick results
)
#saveRDS(m_age0, file="m_age0.rds")
#qsave(m_age0, "m_age0.qs", preset = "high")
lines(centiles.pred(m_age0, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "purple")

#penalised B-spline
m_age1 <- gamlss(
  food_weight ~ pb(age),
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
  food_weight ~ poly(age, 3),
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
#  food_weight ~ cs(age),
#  family = distr_nam,
#  # weights = ds$s_weight,
#  data = ds,
#  method = mixed(20, 20)
#)
#lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")

GAIC.table(m_age0, m_age1, m_age2)
# minimum GAIC(k= 2 ) model: m_age1 
# minimum GAIC(k= 3.84 ) model: m_age1 
# minimum GAIC(k= 8.28 ) model: m_age1 
# df      k=2   k=3.84   k=8.28
# m_age0 5.000000 50636.99 50646.19 50668.39
# m_age1 6.863966 50624.85 50637.48 50667.96
# m_age2 7.000000 50626.37 50639.25 50670.33
# based on GAIC.table, here we will pick m_age1 (pb) 

#centiles(m_age1, xvar = ds$age)
#centiles(m_age2, xvar = ds$age)


# Year
ds[, .(food_weight_median = wtd.quantile(food_weight)), keyby = .(year)][, scatter.smooth(year, food_weight_median, xlim = c(min(ds$year),max(ds$year)))]
#ds_summary <- ds[, .(food_weight_median = wtd.quantile(food_weight)), keyby = .(year)]
#plot(ds_summary$year, ds_summary$food_weight_median, main = "food_weight Median vs. Year")

m_year1 <- gamlss(
  food_weight ~ year,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(50, 100)
)
#qsave(m_year1, "m_year1.qs", preset = "high")
lines(centiles.pred(m_year1, xname = "year", xvalues = min(ds$year):max(ds$year), cent = 50, data = ds), col = "blue1")

m_year2 <- gamlss(
  food_weight ~ log(year),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(50, 100)
)
#qsave(m_year2, "m_year2.qs", preset = "high")
lines(centiles.pred(m_year2, xname = "year", xvalues = min(ds$year):max(ds$year), cent = 50, data = ds), col = "red1")

m_year3 <- gamlss(
  food_weight ~ log(year + 100),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(50, 100)
)
#qsave(m_year3, "m_year3.qs", preset = "high")
lines(centiles.pred(m_year3, xname = "year", xvalues = min(ds$year):max(ds$year), cent = 50, data = ds), col = "green1")

GAIC.table(m_year1, m_year2, m_year3)
# minimum GAIC(k= 2 ) model: m_year1 
# minimum GAIC(k= 3.84 ) model: m_year1 
# minimum GAIC(k= 8.28 ) model: m_year1 
# df      k=2   k=3.84   k=8.28
# m_year1  5 50615.13 50624.33 50646.53
# m_year2  5 50615.80 50625.00 50647.20
# m_year3  5 50615.14 50624.34 50646.54
# based on GAIC.table, is m_year1 (linear)

# bmi
ds[, .(food_weight_median = wtd.quantile(food_weight)), keyby = .(bmi)][, scatter.smooth(bmi, food_weight_median)]

m_bmi1 <- gamlss(
  food_weight ~ bmi,
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50,100) #few iteration to have a quick results
)
lines(centiles.pred(m_bmi1, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "purple")

#for this one, using log, be carefull that your bmi are not 0 because log will not work!!!
m_bmi2 <- gamlss(
  food_weight ~ log(bmi),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 100)
)
lines(centiles.pred(m_bmi2, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "red1")

m_bmi3 <- gamlss(
  food_weight ~ log(bmi+100),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 100)
)
lines(centiles.pred(m_bmi3, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "green1")

m_bmi4 <- gamlss(
  food_weight ~ poly(bmi, 3),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 100)
)
lines(centiles.pred(m_bmi4, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "blue")

m_bmi5 <- gamlss(
  food_weight ~ pb(bmi),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 100)
)
lines(centiles.pred(m_bmi4, xname = "bmi", xvalues = min(ds$bmi):max(ds$bmi), cent = 50, data = ds), col = "yellow")

GAIC.table(m_bmi1, m_bmi2, m_bmi3, m_bmi4, m_bmi5)
# minimum GAIC(k= 2 ) model: m_bmi5 
# minimum GAIC(k= 3.84 ) model: m_bmi2 
# minimum GAIC(k= 8.28 ) model: m_bmi2 
# df      k=2   k=3.84   k=8.28
# m_bmi1 5.000000 50653.54 50662.74 50684.94
# m_bmi2 5.000000 50651.75 50660.95 50683.15
# m_bmi3 5.000000 50653.15 50662.35 50684.55
# m_bmi4 7.000000 50650.86 50663.74 50694.82
# m_bmi5 7.121762 50649.28 50662.38 50694.00
# based on GAIC.table, here we will pick is pb(bmi)


# Energy
ds[, .(food_weight_median = wtd.quantile(food_weight)), keyby = .(energy)][, scatter.smooth(energy, food_weight_median)]

#linear
m_energy1 <- gamlss(
  food_weight ~ energy,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(50, 100) #few iteration to have a quick results
)

lines(centiles.pred(m_energy1, xname = "energy", xvalues = min(ds$energy):max(ds$energy), cent = 50, data = ds), col = "purple")

#log + 100
m_energy2 <- gamlss(
  food_weight ~ log(energy+100),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 100)
)
lines(centiles.pred(m_energy2, xname = "energy", xvalues = min(ds$energy):max(ds$energy), cent = 50, data = ds), col = "red1")

#penalised B-spline
m_energy3 <- gamlss(
  food_weight ~ pb(energy),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(50, 100)
)

lines(centiles.pred(m_energy3, xname = "energy", xvalues = min(ds$energy):max(ds$energy), cent = 50, data = ds), col = "blue1")

#polynomial
m_energy4 <- gamlss(
  food_weight ~ poly(energy, 3),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(50, 100)
)

lines(centiles.pred(m_energy4, xname = "energy", xvalues = min(ds$energy):max(ds$energy), cent = 50, data = ds), col = "red1")


GAIC.table(m_energy1, m_energy2, m_energy3, m_energy4)
# minimum GAIC(k= 2 ) model: m_energy3 
# minimum GAIC(k= 3.84 ) model: m_energy3 
# minimum GAIC(k= 8.28 ) model: m_energy3 
# df      k=2   k=3.84   k=8.28
# m_energy1  5.00000 46978.13 46987.33 47009.53
# m_energy2  5.00000 46086.35 46095.55 46117.75
# m_energy3 18.10698 45974.20 46007.51 46087.91
# m_energy4  7.00000 46165.23 46178.11 46209.19
# We select m_energy3 


#for factors, not much to do. only one smoother
#SES
m_ses1 <- gamlss(
  food_weight ~ SES,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

m_ses2 <- gamlss(
  food_weight ~ pcat(SES), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

GAIC.table(m_ses1, m_ses2)


#Sex
m_sex1 <- gamlss(
  food_weight ~ sex,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

m_sex2 <- gamlss(
  food_weight ~ pcat(sex), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

GAIC.table(m_sex1, m_sex2)


# #Ethnicity
# m_eth1 <- gamlss(
#   food_weight ~ ethnic2,
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# m_eth2 <- gamlss(
#   food_weight ~ pcat(ethnic2), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
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
#   food_weight ~ gor,
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# m_gor2 <- gamlss(
#   food_weight ~ pcat(gor), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# GAIC.table(m_gor1, m_gor2)

##### 
# Theory driven model selection check ->?BCTo, there are FOUR parameters to be included, the FOURTH will not be used.
# This takes a long time to run
# Simplifying the model to reduce computational burden in generating synthetic population
food_weight_model <- gamlss( 
  food_weight ~ year*pb(energy)*SES, 
  ~ year + pb(energy) + SES, #error when including BMI, and energy was predicted by BMI, therefore BMI was not included
  ~ year + pb(energy),
  ~ year + pb(energy),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(500, 1000),
  control = gamlss.control(c.crit = 1e-3)
)

qsave(food_weight_model, "HFSS_food_weight_model_BCTo_nonOOH_byenergyHFSS.qs", preset = "high")
food_weight_model <- qread("HFSS_food_weight_model_BCTo_nonOOH_byenergyHFSS.qs")

model_final <- food_weight_model
model_final

# Code to create a table with predictions
trms <- all.vars(formula(model_final))[-1] # -1 excludes dependent var

#creating newdata which holds all the possible combinations of predictors // you can leave it like that, it's fine
newdata <- CJ(
  year = 3:50, 
  #age = 20:99, 
  #sex = unique(factor(ds$sex)), 
  SES = unique(ds$SES),
  #policy = unique(factor(ds$policy)),
  #gor = unique(ds$gor),
  #ethnic2 = unique(factor(ds$ethnic2)),
  energy = seq(0L,3000L, 10L) #min(ds$bmi) #integer is needed for IMPACT-NCD fitting #I just created until 100 which later can be adjusted to our needs.
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
summary(newdata) 
                

write_fst(newdata, path = "./HFSS_food_weight_table_BCTo_nonOOH_byenergyHFSS.fst", compress = 100L) #I don't save this, the final table was adjusted for its use in IMPACT-NCD
qsave(model_final, "./HFSS_food_weight_modelfinal_BCTo_nonOOH_byenergyHFSS.qs", preset = "high")

print("Table saved")


###LETS SEE THE DIAGNOSTICS
#IF YOU CAN SAVE THE PLOT IN A FILE, IT WILL BE GREAT

if (diagnostics) {
  food_weight_model_final <- qread("./HFSS_food_weight_modelfinal_BCTo_nonOOH_byenergyHFSS.qs")
  
  plot(food_weight_model_final)
  
  wp(food_weight_model_final) # detrended QQ-plot
  wp(resid = resid(food_weight_model_final)) # equivalent to the one above
  wp(resid = resid(food_weight_model_final), ylim.all = 1000 * sqrt(1 / length(resid(food_weight_model_final))))
  tt <- wp(food_weight_model_final, xvar = ds$age, n.inter = 6)
  wp(food_weight_model_final, xvar = ds$year, n.inter = 10)
  wp(resid = resid(food_weight_model_final), xvar = ~ ds$SES)
  
  dtop(food_weight_model_final, xvar = ds$age) # here we want the blue line to be horizontal
  
  rqres.plot(food_weight_model_final)
  rqres.plot(food_weight_model_final, type = "QQ")
}

if (plots) {
  food_weight_model_final <- qread("./HFSS_food_weight_modelfinal_BCTo_nonOOH_byenergyHFSS.qs")
  dir.create("./validation/synthpop_models", recursive = TRUE)
  
  source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  
  
  xlab_nam <- expression(bold(food_weight ~ (grams)))
  food_weight_model_final_tbl <- read_fst("./HFSS_food_weight_table_BCTo_nonOOH_byenergyHFSS.fst", as.data.table = TRUE)
  ds$energysaved <- ds$energy 
  
  #I had issue with merging because of the BMI, so let's do a trick
  ds$energy <- as.integer(round(ds$energy / 10) * 10)

  #food_weight_model_final_tbl$bmi <- food_weight_model_final_tbl$bmi
  #food_weight_model_final_tbl[, bmi := as.integer(bmi)] 
  
  zz <-
    validate_gamlss_tbl(ds[age > 19 & age<=90 & energy>=10L & energy<=3000L], food_weight_model_final_tbl, 50, "food_weight", paste0("q", food_weight_model_final$family[1]))[
      between(food_weight, quantile(food_weight, 0.01), quantile(food_weight, 0.99))
    ]
  
  #summary(zz$bmi)
  zz[, weight := weight / sum(weight), by = "type"]
  #summary(zz)
  #summary(food_weight_model_final_tbl)
  png(
    "./HFSS_food_weight_rel_dist_HFSS_BCTo_nonOOH_byenergyHFSS.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", food_weight],
                      zz[type == "Modelled", food_weight],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      100
  )
  dev.off()
  
  plot_synthpop_val(zz, food_weight, "age_group", "weight", "food_weight by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "sex", "weight", "food_weight by sex", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "year", "weight", "food_weight by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "SES", "weight", "food_weight by SES", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, food_weight, "gor", "weight", "food_weight by region", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, food_weight, "ethnic2", "weight", "food_weight by ethnic", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "bmi_group", "weight", "food_weight by bmi", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "energy_group", "weight", "food_weight by energy", xlab_nam, FALSE, FALSE)
  
  plot_synthpop_val(zz, food_weight, c("year", "age_group"), "weight", "food_weight by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, c("year", "sex"), "weight", "food_weight by year and sex", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, c("year", "SES"), "weight", "food_weight by year and SES", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, food_weight, c("year", "gor"), "weight", "food_weight by year and region", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, food_weight, c("year", "ethnic2"), "weight", "food_weight by year and ethnic", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, c("year", "bmi_group"), "weight", "food_weight by year and bmi", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, c("year", "energy_group"), "weight", "food_weight by year and energy", xlab_nam, FALSE, FALSE)
}


### Editing the data to merge with IMPACTncd model

#svy <- read_fst("/mnt/storage_fast4/Salt_England/IMPACTncd_Engl/synthpop_20/hf_real/synthpop_99ba702352a4ba3f55419baabacfb658_1.fst")
#summary(svy$ethnicity) #in the synthetic population, there are 9 ethic groups:
                        #white, indian, pakistani, bangladeshi, other asian, black caribbean, black african, chinese, other
                        #in the IMPACT-NCD, mixed and other were combined
#summary(svy$sha) #in the synthetic population, North East, North West, Yorkshire and the Humber,
                  #East Midlands, West Midlands, East of England, London, South East Coast,
                  #South Central, South West

food_weight_model_final_tbl <- read_fst("./HFSS_food_weight_table_BCTo_nonOOH_byenergyHFSS.fst", as.data.table = TRUE)
HFSS_food_weight_table_BCTo_nonOOH_byenergyHFSS <- food_weight_model_final_tbl  [, .(year=year,
                                                                    #age=age,
                                                                    #sex=as.factor(fifelse(sex==1, "men", 
                                                                    #               fifelse(sex==2, "women", NA))),
                                                                    qimd=as.factor(fifelse(SES==1, "1 most deprived",
                                                                                   fifelse(SES==2, "2",
                                                                                   fifelse(SES==3, "3",
                                                                                   fifelse(SES==4, "4",
                                                                                   fifelse(SES==5, "5 least deprived", NA)))))),
                                                                    energy_HFSS = energy, 
                                                                    mu=mu,
                                                                    sigma=sigma,
                                                                    nu=nu,
                                                                    tau=tau)]
summary(HFSS_food_weight_table_BCTo_nonOOH_byenergyHFSS)
str(HFSS_food_weight_table_BCTo_nonOOH_byenergyHFSS)
write_fst(HFSS_food_weight_table_BCTo_nonOOH_byenergyHFSS, path = "./HFSS_food_weight_table_final_BCTo_nonOOH_byenergyHFSS.fst", compress = 100L) # This is what my models use as an input

print("Table saved")

food_weight_model_final_tbl[, mm := qBCTo(0.5, mu, sigma, nu, tau)] #based on this, the trends do no look good
food_weight_model_final_tbl[, mean(mm), keyby = year]



tt <- chooseDist(food_weight_model,
                 type = "realplus",
                 trace = TRUE, data = ds,
                 parallel = "multicore", ncpus = 15L # the result is exGAUS, the model remains the same --> model 2
)
tt

# minimum GAIC(k= 2 ) family: exGAUS 
# minimum GAIC(k= 3.84 ) family: exGAUS 
# minimum GAIC(k= 8.28 ) family: exGAUS 
# > tt
# 2     3.84     8.28
# EXP      52676.53 52732.71 52868.28
# GA       46458.47 46530.23 46703.37
# IG             NA       NA       NA
# LOGNO    46122.83 46198.06 46379.60
# LOGNO2   46122.83 46198.06 46379.60
# WEI      47586.93 47646.92 47791.67
# WEI2     47594.83 47652.54 47791.80
# WEI3     47586.92 47646.89 47791.60
# IGAMMA         NA       NA       NA
# PARETO2  52693.35 52762.42 52929.09
# PARETO2o 52741.00 52790.68 52910.56
# GP       52693.35 52762.42 52929.09
# BCCG     45785.02 45849.85 46006.27
# BCCGo          NA       NA       NA
# exGAUS   45759.62 45824.64 45981.52
# GG       45792.20 45875.48 46076.45
# GIG            NA       NA       NA
# LNO            NA       NA       NA
# BCTo     45815.99 45910.16 46137.39
# BCT      45813.27 45903.96 46122.78
# BCPEo    45794.20 45883.84 46100.16
# BCPE     45775.46 45845.68 46015.13
# GB2      45814.51 45904.92 46123.10

#exGAUS takes time to fit, we selected another distribution --> BCPE

food_weight_model2 <- update(food_weight_model, family="BCPE") 

GAIC.table(food_weight_model, food_weight_model2) #model 2 is better

model_final <- food_weight_model2
model_final

# Code to create a table with predictions
trms <- all.vars(formula(model_final))[-1] # -1 excludes dependent var

#creating newdata which holds all the possible combinations of predictors // you can leave it like that, it's fine
newdata <- CJ(
  year = 3:50, 
  #age = 20:99, 
  #sex = unique(factor(ds$sex)), 
  SES = unique(ds$SES),
  #policy = unique(factor(ds$policy)),
  #gor = unique(ds$gor),
  #ethnic2 = unique(factor(ds$ethnic2)),
  energy = seq(0L,3000L, 10L) #min(ds$bmi) #integer is needed for IMPACT-NCD fitting #I just created until 100 which later can be adjusted to our needs.
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
summary(newdata) #mu is negative, therefore BCPE cannot be used.
                

# Therefore, we changed the formula and change the distribution using GG (faster distribution)
food_weight_model2 <- gamlss( 
  food_weight ~ year + pb(energy) + SES, 
  ~ year + pb(energy) + SES, #error when including BMI, and energy was predicted by BMI, therefore BMI was not included
  ~ year + pb(energy),
  ~ year + pb(energy),
  family = "GG",
  #weights = ds$weight,
  data = ds,
  method = mixed(500, 1000),
  control = gamlss.control(c.crit = 1e-3)
)

qsave(food_weight_model2, "HFSS_food_weight_model_GG_nonOOH_byenergyHFSS.qs", preset = "high")
food_weight_model2 <- qread("HFSS_food_weight_model_GG_nonOOH_byenergyHFSS.qs")

model_final <- food_weight_model2
model_final

# Code to create a table with predictions
trms <- all.vars(formula(model_final))[-1] # -1 excludes dependent var

#creating newdata which holds all the possible combinations of predictors // you can leave it like that, it's fine
newdata <- CJ(
  year = 3:50, 
  #age = 20:99, 
  #sex = unique(factor(ds$sex)), 
  SES = unique(ds$SES),
  #policy = unique(factor(ds$policy)),
  #gor = unique(ds$gor),
  #ethnic2 = unique(factor(ds$ethnic2)),
  energy = seq(0L,3000L, 10L) #min(ds$bmi) #integer is needed for IMPACT-NCD fitting #I just created until 100 which later can be adjusted to our needs.
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
summary(newdata) 

write_fst(newdata, path = "./HFSS_food_weight_table_GG_nonOOH_byenergyHFSS.fst", compress = 100L) #I don't save this, the final table was adjusted for its use in IMPACT-NCD
qsave(model_final, "./HFSS_food_weight_modelfinal_GG_nonOOH_byenergyHFSS.qs", preset = "high")

newdata[, mm := qGG(0.5, mu, sigma, nu)]
newdata[, mean(mm), keyby = year]

print("Table saved")


###LETS SEE THE DIAGNOSTICS
#IF YOU CAN SAVE THE PLOT IN A FILE, IT WILL BE GREAT

if (diagnostics) {
  food_weight_model_final <- qread("./HFSS_food_weight_modelfinal_GG_nonOOH_byenergyHFSS.qs")
  
  plot(food_weight_model_final)
  
  wp(food_weight_model_final) # detrended QQ-plot
  wp(resid = resid(food_weight_model_final)) # equivalent to the one above
  wp(resid = resid(food_weight_model_final), ylim.all = 1000 * sqrt(1 / length(resid(food_weight_model_final))))
  tt <- wp(food_weight_model_final, xvar = ds$age, n.inter = 6)
  wp(food_weight_model_final, xvar = ds$year, n.inter = 10)
  wp(resid = resid(food_weight_model_final), xvar = ~ ds$SES)
  
  dtop(food_weight_model_final, xvar = ds$age) # here we want the blue line to be horizontal
  
  rqres.plot(food_weight_model_final)
  rqres.plot(food_weight_model_final, type = "QQ")
}

if (plots) {
  food_weight_model_final <- qread("./HFSS_food_weight_modelfinal_GG_nonOOH_byenergyHFSS.qs")
  dir.create("./validation/synthpop_models", recursive = TRUE)
  
  source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  
  
  xlab_nam <- expression(bold(food_weight ~ (grams)))
  food_weight_model_final_tbl <- read_fst("./HFSS_food_weight_table_GG_nonOOH_byenergyHFSS.fst", as.data.table = TRUE)
  ds$energysaved <- ds$energy 
  
  #I had issue with merging because of the BMI, so let's do a trick
  ds$energy <- as.integer(round(ds$energy / 10) * 10)

  #food_weight_model_final_tbl$bmi <- food_weight_model_final_tbl$bmi
  #food_weight_model_final_tbl[, bmi := as.integer(bmi)] 
  
  zz <-
    validate_gamlss_tbl(ds[age > 19 & age<=90 & energy>=10L & energy<=3000L], food_weight_model_final_tbl, 50, "food_weight", paste0("q", food_weight_model_final$family[1]))[
      between(food_weight, quantile(food_weight, 0.01), quantile(food_weight, 0.99))
    ]
  
  #summary(zz$bmi)
  zz[, weight := weight / sum(weight), by = "type"]
  #summary(zz)
  #summary(food_weight_model_final_tbl)
  png(
    "./HFSS_food_weight_rel_dist_HFSS_GG_nonOOH_byenergyHFSS.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", food_weight],
                      zz[type == "Modelled", food_weight],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      100
  )
  dev.off()
  
  plot_synthpop_val(zz, food_weight, "age_group", "weight", "food_weight by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "sex", "weight", "food_weight by sex", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "year", "weight", "food_weight by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "SES", "weight", "food_weight by SES", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, food_weight, "gor", "weight", "food_weight by region", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, food_weight, "ethnic2", "weight", "food_weight by ethnic", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "bmi_group", "weight", "food_weight by bmi", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, "energy_group", "weight", "food_weight by energy", xlab_nam, FALSE, FALSE)
  
  plot_synthpop_val(zz, food_weight, c("year", "age_group"), "weight", "food_weight by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, c("year", "sex"), "weight", "food_weight by year and sex", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, c("year", "SES"), "weight", "food_weight by year and SES", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, food_weight, c("year", "gor"), "weight", "food_weight by year and region", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, food_weight, c("year", "ethnic2"), "weight", "food_weight by year and ethnic", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, c("year", "bmi_group"), "weight", "food_weight by year and bmi", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, food_weight, c("year", "energy_group"), "weight", "food_weight by year and energy", xlab_nam, FALSE, FALSE)
}


### Editing the data to merge with IMPACTncd model

#svy <- read_fst("/mnt/storage_fast4/Salt_England/IMPACTncd_Engl/synthpop_20/hf_real/synthpop_99ba702352a4ba3f55419baabacfb658_1.fst")
#summary(svy$ethnicity) #in the synthetic population, there are 9 ethic groups:
                        #white, indian, pakistani, bangladeshi, other asian, black caribbean, black african, chinese, other
                        #in the IMPACT-NCD, mixed and other were combined
#summary(svy$sha) #in the synthetic population, North East, North West, Yorkshire and the Humber,
                  #East Midlands, West Midlands, East of England, London, South East Coast,
                  #South Central, South West

food_weight_model_final_tbl <- read_fst("./HFSS_food_weight_table_GG_nonOOH_byenergyHFSS.fst", as.data.table = TRUE)
HFSS_food_weight_table_GG_nonOOH_byenergyHFSS <- food_weight_model_final_tbl  [, .(year=year,
                                                                    #age=age,
                                                                    #sex=as.factor(fifelse(sex==1, "men", 
                                                                    #               fifelse(sex==2, "women", NA))),
                                                                    qimd=as.factor(fifelse(SES==1, "1 most deprived",
                                                                                   fifelse(SES==2, "2",
                                                                                   fifelse(SES==3, "3",
                                                                                   fifelse(SES==4, "4",
                                                                                   fifelse(SES==5, "5 least deprived", NA)))))),
                                                                    hfss_energy = energy, 
                                                                    mu=mu,
                                                                    sigma=sigma,
                                                                    nu=nu)]
                                                                    #tau=tau)]
summary(HFSS_food_weight_table_GG_nonOOH_byenergyHFSS)
str(HFSS_food_weight_table_GG_nonOOH_byenergyHFSS)
write_fst(HFSS_food_weight_table_GG_nonOOH_byenergyHFSS, path = "./HFSS_food_weight_table_final_GG_nonOOH_byenergyHFSS.fst", compress = 100L) # This is what my models use as an input

print("Table saved")
