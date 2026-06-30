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

#Adding salt total to the prediction
salt_total <- readRDS("salt_diet_consump.rds")
setDT(salt_total)
salt_total <- salt_total [, .(seriali, salt_total = sodium_diet_mg/393.4)] #to make consistent with salt calculation for HFSS
salt_total$seriali <- as.character(salt_total$seriali)
ds_ <- ds_[salt_total, on= "seriali"]

test <- ds_[, diff:= (salt-salt_total)]
summary(test$diff) 
test <- test[diff>0,] #no salt HFSS > total_salt

#prop salt from HFSS
test <- ds_[, prop:= (salt/salt_total)]
mean(ds_$prop)
weighted.mean(ds_$prop, w = ds_$s_weight_diary) #58%

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
  salt = salt,
  #food_weight = TotalGrams,
  salt_total = salt_total,
  salt_total_group = fifelse(salt_total<=10, "<=10", ">10"),
  year= as.integer(surveyyr+8), # survey years 1 - 11, representing 2009 - 2019
  weight= s_weight_diary,
  #ethnic = as.factor(ethnic_5), #1 White, 2 Mixed ethnic group, 3 Black or Black British, 4 Asian or asian British, 5 Other
  #ethnic2 = as.factor(ifelse(ethnic_5==5, 2, ethnic_5)), #for IMPACT-NCD model, Mixed and Other were combined together
  #gor= as.factor(gor), #1 North East, 2 North West, 3 Yorkshire and the Humber, 4 East Midlands, 5 West Midlands, 6 East of England, 7 London, 8 South East, 9 South West
  OOH= OOH
)]

summary(ds)
ds[, salt:= salt*(1-OOH)] #here we estimated salt from non-OOH
#ds[, food_weight:= food_weight*(1-OOH)] #here we estimated salt from non-OOH
weighted.mean(ds$salt, ds$weight)
#3.240978

ds <- ds[complete.cases(ds), ]
sum(is.na(ds)) # 3925 obs
#View(ds)

ds <- ds[year > 12, ] #selecting year from 2013 to 2019 because of this is predicted by GAMLSS MTL #2539 obs
summary(ds)

# # creating a variable in which before and after the introduction of TLS
# ds [, policy:= as.factor(ifelse(year<13, 1, 2))]
# summary(ds)

# First identify what distribution is best fit for salt consumption (the left hand side variable)
#dev.off()
plot(density(ds$salt))
summary(ds$salt)
# This is a continuous distribution that cannot be 0 
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.001112  2.006059  2.781615  3.025152  3.762830 21.361274



# Let's see what distribution work best for this. We will use the function
# fitDist from gamlss. Please have a read of the documentation for this function
marg_distr <- fitDist(salt,
                      k = 2, # Default is for AIC. use log(length(ds$salt)) for BIC
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
# BCTo      BCT      ST5      JSU     JSUo    BCPEo     BCPE       RG   exGAUS 
# 8700.380 8700.380 8709.893 8716.543 8716.543 8718.368 8718.368 8721.143 8721.523 
# ST1 
# 8721.553 

#class(marg_distr)
#str(marg_distr)
#saveRDS(marg_distr, file="marg_distr.rds") #it is better to save the object as the R studio is sometime error
#qsave(marg_distr, "marg_distr.qs", preset = "high")
#View(marg_distr)
#marg_distr <- qread("marg_distr.qs")

# lets plot some of the best three distributions
m1 <- gamlss(salt ~ 1,
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

plot(density(ds$salt), lwd = 1, ylim = c(0, 0.40))
lines(density(sample_from_gamlss(m1)), col = alpha("red", 0.4))

m2 <- gamlss(salt ~ 1, family = names(marg_distr$fits[2]), data = ds, method = mixed(20, 40))
#saveRDS(m2, file="m2.rds")
#qsave(m2, "m2.qs", preset = "high")
#m2 <- qread("m2.qs")
lines(density(sample_from_gamlss(m2)), col = alpha("blue", 0.4))

m3 <- gamlss(salt ~ 1, family = names(marg_distr$fits[3]), data = ds, method = mixed(20, 40))
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
# minimum GAIC(k= 2 ) model: m1 
# minimum GAIC(k= 3.84 ) model: m1 
# minimum GAIC(k= 7.84 ) model: m1 
# df      k=2   k=3.84   k=7.84
# m1  4 8700.380 8707.740 8723.740
# m2  4 8700.380 8707.740 8723.740
# m3  4 8710.287 8717.647 8733.647

distr_nam <- names(marg_distr$fits[1]) # based on GAIC.table, and better fit based on the plot, BCTo
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3

distr_nam

# Second, identify what distribution is best fit for predictors (the right hand side variable)
# Predictors
# Age
ds[, .(salt_median = wtd.quantile(salt)), keyby = .(age)][, scatter.smooth(age, salt_median)]
#ds[, .(salt_median = wtd.quantile(salt, weight = s_weight)), keyby = .(age)][, scatter.smooth(age, salt_median)]

#linear
m_age0 <- gamlss(
  salt ~ age,
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
  salt ~ pb(age),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)
#saveRDS(m_age1, file="m_age1.rds")
#qsave(m_age1, "m_age1.qs", preset = "high")
lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")

#polynomial
m_age2 <- gamlss(
  salt ~ poly(age, 3),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)
#saveRDS(m_age2, file="m_age2.rds")
#qsave(m_age2, "m_age2.qs", preset = "high")
lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")

#cubic-spline
#m_age3 <- gamlss(
#  salt ~ cs(age),
#  family = distr_nam,
#  # weights = ds$s_weight,
#  data = ds,
#  method = mixed(20, 20)
#)
#lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")

GAIC.table(m_age0, m_age1, m_age2)
# minimum GAIC(k= 2 ) model: m_age2 
# minimum GAIC(k= 3.84 ) model: m_age0 
# minimum GAIC(k= 7.84 ) model: m_age0 
# df      k=2   k=3.84   k=7.84
# m_age0 5.000000 8626.459 8635.659 8655.659
# m_age1 5.018104 8626.471 8635.704 8655.776
# m_age2 7.000000 8623.338 8636.218 8664.218
# based on GAIC.table, here we will pick m_age2 (poly) 

#centiles(m_age1, xvar = ds$age)
#centiles(m_age2, xvar = ds$age)


# Year
ds[, .(salt_median = wtd.quantile(salt)), keyby = .(year)][, scatter.smooth(year, salt_median, xlim = c(min(ds$year),max(ds$year)))]
#ds_summary <- ds[, .(salt_median = wtd.quantile(salt)), keyby = .(year)]
#plot(ds_summary$year, ds_summary$salt_median, main = "salt Median vs. Year")

m_year1 <- gamlss(
  salt ~ year,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)
#qsave(m_year1, "m_year1.qs", preset = "high")
lines(centiles.pred(m_year1, xname = "year", xvalues = min(ds$year):max(ds$year), cent = 50, data = ds), col = "blue1")

m_year2 <- gamlss(
  salt ~ log(year),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)
#qsave(m_year2, "m_year2.qs", preset = "high")
lines(centiles.pred(m_year2, xname = "year", xvalues = min(ds$year):max(ds$year), cent = 50, data = ds), col = "red1")

m_year3 <- gamlss(
  salt ~ log(year + 100),
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)
#qsave(m_year3, "m_year3.qs", preset = "high")
lines(centiles.pred(m_year3, xname = "year", xvalues = min(ds$year):max(ds$year), cent = 50, data = ds), col = "green1")

GAIC.table(m_year1, m_year2, m_year3)
# minimum GAIC(k= 2 ) model: m_year1 
# minimum GAIC(k= 3.84 ) model: m_year1 
# minimum GAIC(k= 7.84 ) model: m_year1 
# df      k=2   k=3.84   k=7.84
# m_year1  5 8699.485 8708.685 8728.685
# m_year2  5 8699.980 8709.180 8729.180
# m_year3  5 8699.554 8708.754 8728.754
# based on GAIC.table, here we will pick m_year1 (linear) 

#centiles(m_year1, xvar = ds$year)
#centiles(m_year2, xvar = ds$year)

# salt_total
ds[, .(salt_median = wtd.quantile(salt)), keyby = .(salt_total)][, scatter.smooth(salt_total, salt_median)]

m_salt_total1 <- gamlss(
  salt ~ salt_total,
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50,20) #few iteration to have a quick results
)
lines(centiles.pred(m_salt_total1, xname = "salt_total", xvalues = min(ds$salt_total):max(ds$salt_total), cent = 50, data = ds), col = "purple")

#for this one, using log, be carefull that your bmi are not 0 because log will not work!!!
m_salt_total2 <- gamlss(
  salt ~ log(salt_total),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20)
)
lines(centiles.pred(m_salt_total2, xname = "salt_total", xvalues = min(ds$salt_total):max(ds$salt_total), cent = 50, data = ds), col = "red1")

m_salt_total3 <- gamlss(
  salt ~ log(salt_total+100),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20)
)
lines(centiles.pred(m_salt_total3, xname = "salt_total", xvalues = min(ds$salt_total):max(ds$salt_total), cent = 50, data = ds), col = "green1")

m_salt_total4 <- gamlss(
  salt ~ poly(salt_total, 3),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20)
)
lines(centiles.pred(m_salt_total4, xname = "salt_total", xvalues = min(ds$salt_total):max(ds$salt_total), cent = 50, data = ds), col = "blue")

m_salt_total5 <- gamlss(
  salt ~ pb(salt_total),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(50, 20) 
)
lines(centiles.pred(m_salt_total5, xname = "salt_total", xvalues = min(ds$salt_total):max(ds$salt_total), cent = 50, data = ds), col = "yellow")

GAIC.table(m_salt_total1, m_salt_total2, m_salt_total3, m_salt_total4, m_salt_total5)
# minimum GAIC(k= 2 ) model: m_salt_total5 
# minimum GAIC(k= 3.84 ) model: m_salt_total2 
# minimum GAIC(k= 7.84 ) model: m_salt_total2 
# df      k=2   k=3.84   k=7.84
# m_salt_total1  5.00000 5409.962 5419.162 5439.162
# m_salt_total2  5.00000 5076.265 5085.465 5105.465
# m_salt_total3  5.00000 5379.637 5388.837 5408.837
# m_salt_total4  7.00000 5117.297 5130.177 5158.177
# m_salt_total5 11.82278 5074.790 5096.544 5143.835
# based on GAIC.table, here we will pick m_salt_total5 (pb) 

#for factors, not much to do. only one smoother
#SES
m_ses1 <- gamlss(
  salt ~ SES,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

m_ses2 <- gamlss(
  salt ~ pcat(SES), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

GAIC.table(m_ses1, m_ses2)


#Sex
m_sex1 <- gamlss(
  salt ~ sex,
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

m_sex2 <- gamlss(
  salt ~ pcat(sex), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
  family = distr_nam,
  # weights = ds$s_weight,
  data = ds,
  method = mixed(20, 20)
)

GAIC.table(m_sex1, m_sex2)


# #Ethnicity
# m_eth1 <- gamlss(
#   salt ~ ethnic,
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# m_eth2 <- gamlss(
#   salt ~ pcat(ethnic), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# GAIC.table(m_eth1, m_eth2) ### There is an error for using pcat for ethnicity


# #GOR
# m_gor1 <- gamlss(
#   salt ~ gor,
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# m_gor2 <- gamlss(
#   salt ~ pcat(gor), #pcat is a smoother for categorical variables, allow to see if its better to collapse different levels
#   family = distr_nam,
#   # weights = ds$s_weight,
#   data = ds,
#   method = mixed(20, 20)
# )
# 
# GAIC.table(m_gor1, m_gor2)


#####
# Theory driven model selection check ->?BCTo, there are four parameters to be included
# This takes a long time to run
salt_model <- gamlss(
  salt ~ year + poly(age,3) + pb(salt_total) + sex + SES, #this formula is better fit
  ~ year + poly(age,3) + sex + SES,      
  ~ year + poly(age,3),
  ~ year + poly(age,3),
  family = distr_nam,
  #weights = ds$weight,
  data = ds,
  method = mixed(100, 500),
  control = gamlss.control(c.crit = 1e-3)
)

qsave(salt_model, "HFSS_salt_model_BCTo_nonOOH_bysalttotal.qs", preset = "high") #the model is converged

model_final <- salt_model

# Code to create a table with predictions
trms <- all.vars(formula(model_final))[-1] # -1 excludes dependent var

#creating newdata which holds all the possible combinations of predictors // you can leave it like that, it's fine
newdata <- CJ(
  year = 3:50, 
  age = 20:99, 
  sex = unique(factor(ds$sex)), 
  SES = unique(factor(ds$SES)),
  salt_total = seq(0L,35L, 1L)
  #policy = unique(factor(ds$policy)),
  #gor = unique(ds$gor),
  #ethnic2 = unique(factor(ds$ethnic2)),
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

write_fst(newdata, path = "./HFSS_salt_table_BCTo_nonOOH_bysalttotal.fst", compress = 100L) # This is what my models use as an input
qsave(model_final, "./HFSS_salt_modelfinal_BCTo_nonOOH_bysalttotal.qs", preset = "high")

print("Table saved")


###LETS SEE THE DIAGNOSTICS
#IF YOU CAN SAVE THE PLOT IN A FILE, IT WILL BE GREAT

if (diagnostics) {
  salt_model_final <- qread("./HFSS_salt_modelfinal_BCTo_nonOOH_bysalttotal.qs")
  
  plot(salt_model_final)
  
  wp(salt_model_final) # detrended QQ-plot
  wp(resid = resid(salt_model_final)) # equivalent to the one above
  wp(resid = resid(salt_model_final), ylim.all = 1000 * sqrt(1 / length(resid(salt_model_final))))
  tt <- wp(salt_model_final, xvar = ds$age, n.inter = 6)
  wp(salt_model_final, xvar = ds$year, n.inter = 10)
  wp(resid = resid(salt_model_final), xvar = ~ ds$SES)
  
  dtop(salt_model_final, xvar = ds$age) # here we want the blue line to be horizontal
  
  rqres.plot(salt_model_final)
  rqres.plot(salt_model_final, type = "QQ")
}

if (plots) {
  salt_model_final <- qread("./HFSS_salt_modelfinal_BCTo_nonOOH_bysalttotal.qs")
  dir.create("./validation/synthpop_models", recursive = TRUE)
  
  source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  
  
  xlab_nam <- expression(bold(salt ~ (gram)))
  salt_model_final_tbl <- read_fst("./HFSS_salt_table_BCTo_nonOOH_bysalttotal.fst", as.data.table = TRUE)
  ds$salt_total_saved <- ds$salt_total
  
  ds$salt_total <- as.integer(round(ds$salt_total,0))
  
  zz <-
    validate_gamlss_tbl(ds[age > 19 & age<=90 & salt_total>=0 & salt_total<=35], salt_model_final_tbl, 50, "salt", paste0("q", salt_model_final$family[1]))[
      between(salt, quantile(salt, 0.01), quantile(salt, 0.99))
    ]
  
  zz[, weight := weight / sum(weight), by = "type"]
  #summary(zz)
  #summary(salt_model_final_tbl)
  png(
    "./HFSS_salt_rel_dist_HFSS_BCTo_nonOOH_bysalttotal.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", salt],
                      zz[type == "Modelled", salt],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      100
  )
  dev.off()
  
  plot_synthpop_val(zz, salt, "age_group", "weight", "salt by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, salt, "sex", "weight", "salt by sex", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, salt, "year", "weight", "salt by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, salt, "SES", "weight", "salt by SES", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, salt, "gor", "weight", "salt by region", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, salt, "ethnic2", "weight", "salt by ethnic", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, salt, "salt_total_group", "weight", "salt by salt total", xlab_nam, FALSE, FALSE)
  
  plot_synthpop_val(zz, salt, c("year", "age_group"), "weight", "salt by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, salt, c("year", "sex"), "weight", "salt by year and sex", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, salt, c("year", "SES"), "weight", "salt by year and SES", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, salt, c("year", "gor"), "weight", "salt by year and region", xlab_nam, FALSE, FALSE)
  #plot_synthpop_val(zz, salt, c("year", "ethnic2"), "weight", "salt by year and ethnic", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, salt, c("year", "salt_total_group"), "weight", "salt by year and salt total", xlab_nam, FALSE, FALSE)
}

### Editing the data to merge with IMPACTncd model

#svy <- read_fst("/mnt/storage_fast4/Salt_England/IMPACTncd_Engl/synthpop_20/hf_real/synthpop_99ba702352a4ba3f55419baabacfb658_1.fst")
#summary(svy$ethnicity) #in the synthetic population, there are 9 ethic groups:
                        #white, indian, pakistani, bangladeshi, other asian, black caribbean, black african, chinese, other
                        #in the IMPACT-NCD, mixed and other were combined
#summary(svy$sha) #in the synthetic population, North East, North West, Yorkshire and the Humber,
                  #East Midlands, West Midlands, East of England, London, South East Coast,
                  #South Central, South West

salt_model_final_tbl <- read_fst("./HFSS_salt_table_BCTo_nonOOH_bysalttotal.fst", as.data.table = TRUE)
HFSS_salt_table_BCTo_nonOOH_bysalttotal <- salt_model_final_tbl  [, .(year=year,
                                                                    age=age,
                                                                    sex=as.factor(fifelse(sex==1, "men", 
                                                                                  fifelse(sex==2, "women", NA))),
                                                                    qimd=as.factor(fifelse(SES==1, "1 most deprived",
                                                                                   fifelse(SES==2, "2",
                                                                                   fifelse(SES==3, "3",
                                                                                   fifelse(SES==4, "4",
                                                                                   fifelse(SES==5, "5 least deprived", NA)))))),
                                                                    salt_total = salt_total,
                                                                    mu=mu,
                                                                    sigma=sigma,
                                                                    nu=nu,
                                                                    tau=tau)]
summary(HFSS_salt_table_BCTo_nonOOH_bysalttotal)
str(HFSS_salt_table_BCTo_nonOOH_bysalttotal)
write_fst(HFSS_salt_table_BCTo_nonOOH_bysalttotal, path = "./HFSS_salt_table_final_BCTo_nonOOH_bysalttotal.fst", compress = 100L) # This is what my models use as an input

print("Table saved")

