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

# Final distribution chosen: ST5 

print("tchol")

library(IMPACTncdEngl)
library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk

# Set some variables we will use later in automation
distributions <- FALSE
univariable_analysis <- FALSE
runmodel <- FALSE
update <- TRUE
diagnostics <- FALSE
plots <- FALSE
seed                <- 43L



# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)

# Complete case analysis for the variables fitting in the model, including 
# a valid 'nurse weight' (this is calculated by HSE to account for non-response rates) 
ds <- na.omit(HSE_ts[wt_blood > 0 & age >= 20,
                     .(tchol, year, age, agegrp10, sex, qimd, ethnicity_grp, sha, smok_status, wt_blood, urban_rural, bmi)])
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
  # First identify what distribution is best fit for tchol (the left hand side
  # variable)
  plot(density(ds$tchol))
  summary(ds$tchol)
  # This is a continuous distribution that cannot be 0
  
  # Let's see what distribution work best for this. We will use the function
  # fitDist from gamlss - this fits all the relevant distributions as specified 
  # by 'type' (see the documentation) and then orders them by some criteria (AIC is the default) 
  # Please have a read of the documentation for this function
  marg_distr <- fitDist(tchol,
                        k = 2, # Default is for AIC. use log(length(ds$tchol)) for BIC
                        type = "realAll", # Other options "realline", "realplus", "real0to1", "counts", "binom"
                        try.gamlss = TRUE, # Better but slow
                        extra = NULL,
                        data = ds,
                        trace = TRUE, 
                        weights = ds$wt_blood
  )
  
  # Check the best 10 distr - look at their properties, how many parameters? 
  head(marg_distr$fits, 10)
  # ST5      BCT     BCTo     JSUo      JSU      GB2     EGB2    BCPEo     BCPE    BCCGo 
  # 166129.3 166129.8 166129.8 166131.1 166131.1 166131.5 166132.2 166154.2 166154.2 166155.5 
  # 
  # This function does "distribution validation diagnostics based on a fitted distribution model"
  # It is from the IMPACRncdEngl package, so if you don't have this loaded it won't work.
  # It also takes a long time
  distr_validation(marg_distr, ds[between(tchol, 0, 15), .(var = tchol, wt = wt_blood)],
                   expression(bold(tchol ~ (mg/dL))))
  
  # lets plot some of them, starting with the best AIC 'marg_distr$fits[1]'
  m1 <- gamlss(tchol ~ 1,
               family = names(marg_distr$fits[1]),
               weights = ds$wt_blood,
               # method = mixed(20, 20)
               data = ds
  )
  
  # Let's create a helper function to plot the fits from different distributions 
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
  
  plot(density(ds$tchol), lwd = 3, ylim = c(0, 0.4))
  lines(density(sample_from_gamlss(m1)), col = alpha("red", 0.4))
  
  # The 2nd best 
  m2 <- gamlss(tchol ~ 1, family = names(marg_distr$fits[2]), data = ds)
  lines(density(sample_from_gamlss(m2)), col = alpha("blue", 0.4))
  
  # The 4th best
  m3 <- gamlss(tchol ~ 1, family = names(marg_distr$fits[4]), data = ds)
  lines(density(sample_from_gamlss(m3)), col = alpha("green", 0.4))
  # they all look terrible 
  
  # An alternative plot (skip these if the reldist package doesn't work)
  reldist(sample_from_gamlss(m1), ds$tchol, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m2), ds$tchol, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m3), ds$tchol, method = "bgk", bar = TRUE, show = "effect")
  

  
  
  distr_nam <- names(marg_distr$fits[1]) # I will pick ST5
  
}

distr_nam <- "ST5"

print(distr_nam)
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3

if (univariable_analysis) {
  
# Age
ds[, .(tchol_median = wtd.quantile(tchol, weight = wt_blood)), keyby = .(age)][, scatter.smooth(age, tchol_median)]

m_age0 <- gamlss(
  tchol ~ age,
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)  # Number of iterations using Rigby and Stasinopoulos algorith (default), 
  # followed by number using Cole and Green algorithm 
)
lines(centiles.pred(m_age0, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "purple")

m_age1 <- gamlss(
  tchol ~ pb(age), # penalised beta spline - an automated way for finding the splines & knots
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")

m_age2 <- gamlss(
  tchol ~ poly(age, degree = 3),#polynomial
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")

m_age3 <- gamlss(
  tchol ~ cs(age),
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")
GAIC.table(m_age0, m_age1, m_age2, m_age3)

#This is another way of looking at the different models 
centiles(m_age1, xvar = ds$age)
centiles(m_age2, xvar = ds$age)

# m_age1 - pb spline is the best


ds[, .(tchol_median = wtd.quantile(tchol, weight = wt_blood)), keyby = .(year)][
  , scatter.smooth(year, tchol_median, xlim = c(3, 40), ylim = c(4, 6))]


# Year
m_year0 <- gamlss(
  tchol ~ year,
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_year0, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "blue1")

m_year1 <- gamlss(
  tchol ~ log(year),
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "red1")

m_year2 <- gamlss(
  tchol ~ log(year + 10),
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "green1")

m_year3 <- gamlss(
  tchol ~ log(year + 100),
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "pink")

GAIC.table(m_year0, m_year1, m_year2, m_year3) #m_year3 is the best 
centiles(m_year0, xvar = ds$year)
centiles(m_year1, xvar = ds$year)
centiles(m_year2, xvar = ds$year)
centiles(m_year3, xvar = ds$year)


# Smoking status 
m_sm1 <- gamlss(
  tchol ~ smok_status,
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
summary(m_sm1)

# Using a smoother 
m_sm2 <- gamlss(
  tchol ~ pcat(smok_status, method = "GAIC"),
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
summary(m_sm2)

GAIC.table(m_sm1, m_sm2) #without pcat  is better

#Looking at the levels pcat has reduced smok_status to, and how they 
#correspond to the original levels
f2 <- getSmo(m_sm2)$factor
levels(f2)
table(f2,ds$smok_status)

#In theory you can plot this, but it seems to take an incredibly long tine 
# plotLambda(tchol,factor=smok_status,data=ds,family=JSU,along=seq(-10,2,.2))
# abline(v=log(getSmo(m_sm2)$lambda),col='gray',lwd=2.5)

# BMI
ds[, .(bmi_median = wtd.quantile(bmi, weight = wt_blood)), keyby = .(bmi)][, scatter.smooth(bmi, bmi_median)]

m_bmi0 <- gamlss(
  tchol ~ bmi,
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_bmi0, xname = "bmi", xvalues = 20:90, cent = 50, data = ds), col = "purple")

m_bmi1 <- gamlss(
  tchol ~ pb(bmi), 
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_bmi1, xname = "bmi", xvalues = 20:90, cent = 50, data = ds), col = "blue1")

m_bmi2 <- gamlss(
  tchol ~ poly(bmi, 3),
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_bmi2, xname = "bmi", xvalues = 20:90, cent = 50, data = ds), col = "red1")

m_bmi3 <- gamlss(
  tchol ~ cs(bmi),
  family = distr_nam,
  weights = ds$wt_blood,
  data = ds,
  method = mixed(100, 100)
)
lines(centiles.pred(m_bmi3, xname = "bmi", xvalues = 20:90, cent = 50, data = ds), col = "green1")
GAIC.table(m_bmi0, m_bmi1, m_bmi2, m_bmi3) #bmi3 has v slightly better aic than pb spline, but pb spline better bic
centiles(m_bmi1, xvar = ds$bmi)
centiles(m_bmi2, xvar = ds$bmi)

# m_bmi1 - using pb spline 


}


if(runmodel){

# mod_min <- gamlss(
#   tchol ~ log(year + 10) + pb(age) + pcat(qimd) +
#     sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural + pb(bmi) ,
#   family = distr_nam,
#   weights = ds$wt_blood,
#   data = ds,
#   method = mixed(5, 100),
#   control = con1 # for speed but looses accuracy if not 1e-3
# )
# # qsave(mod_min, "./inputs/exposure_distributions/AH_test/tchol_mod_min.qs", preset = "high")
# # warnings() # check for warnings
# 
# mod_min2 <- gamlss(
#   tchol ~ log(year + 10) + pb(age) + pcat(qimd) +
#     sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural + pb(bmi) ,
#   family = "BCT",
#   weights = ds$wt_blood,
#   data = ds,
#   method = mixed(5, 100),
#   control = con1 # for speed but looses accuracy if not 1e-3
# )
# # qsave(mod_min2, "./inputs/exposure_distributions/AH_test/tchol_mod_min_BCT.qs", preset = "high")
# # warnings() # check for warnings
# GAIC.table(mod_min, mod_min2)
# minimum GAIC(k= 2 ) model: mod_min2 
# minimum GAIC(k= 3.84 ) model: mod_min2 
# minimum GAIC(k= 10.61 ) model: mod_min2 
# df      k=2   k=3.84  k=10.61
# mod_min  29.15275 116244.5 116298.2 116495.6
# mod_min2 29.78201 115761.8 115816.6 116018.2


# mod_min3 <- gamlss(
#   tchol ~ log(year + 10) + pb(age) + pcat(qimd) +
#     sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural + pb(bmi) ,
#   family = "BCTo",
#   weights = ds$wt_blood,
#   data = ds,
#   method = mixed(5, 100),
#   control = con1 # for speed but looses accuracy if not 1e-3
# )
#  qsave(mod_min3, "./inputs/exposure_distributions/AH_test/tchol_mod_min_BCTo.qs", preset = "high")

 # mod_min4 <- gamlss(
 #   tchol ~ log(year + 10) + pb(age) + pcat(qimd) +
 #     sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural  ,
 #   family = distr_nam,
 #   weights = ds$wt_blood,
 #   data = ds,
 #   method = mixed(5, 100),
 #   control = con1 # for speed but looses accuracy if not 1e-3
 # )
 #  qsave(mod_min4, "./inputs/exposure_distributions/AH_test/tchol_mod_min_nobmi.qs", preset = "high")
 # warnings() # check for warnings
 
  # mod_min5 <- gamlss(
  #   tchol ~ log(year + 100) + pb(age) + pcat(qimd) +
  #     sex  + pcat(ethnicity_grp) + pcat(sha) + urban_rural  ,
  #   family = distr_nam,
  #   weights = ds$wt_blood,
  #   data = ds,
  #   method = mixed(100, 100),
  #   control = con1 # for speed but looses accuracy if not 1e-3
  # )
  # qsave(mod_min5, "./inputs/exposure_distributions/AH_test/tchol_mod_min_nobmi_nosmk.qs", preset = "high")
  # warnings() # check for warnings
  
  mod_min <- gamlss(
    tchol ~ log(year + 100)   ,
    family = distr_nam,
    weights = ds$wt_blood,
    data = ds,
    method = mixed(100, 100),
    control = con1 # for speed but looses accuracy if not 1e-3
  )
  qsave(mod_min, "./inputs/exposure_distributions/AH_test/tchol_mod_min.qs", preset = "high")
  warnings() # check for warnings

  # Stepwise model selection using a Generalized Akaike Information Criterion
tchol_modelA <- stepGAICAll.A(
  mod_min,
  scope = list(
    lower = ~ log(year + 100) ,
    upper = ~ (log(year + 100) + pb(age) + pcat(qimd) + 
                 sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural  + pb(bmi))^2
  ),
  sigma.scope = list(
    lower = ~1,
    upper = ~ log(year + 100) + pb(age) + pcat(qimd) + 
      sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural + pb(bmi)
  ),
  nu.scope = list(
    lower = ~1,
    upper = ~ log(year + 100) + pb(age) + pcat(qimd) + 
      sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural + pb(bmi)
  ),
  tau.scope = list(
    lower = ~1,
    upper = ~ log(year + 100) + pb(age) + pcat(qimd) + 
      sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural + pb(bmi)
  ),
  parallel = "multicore",
  ncpus = 12L,
  weights = wt_blood, 
  trace = TRUE  # This is so you can see all the models it has tried to fit 
)
warnings()

tchol_modelA

tchol_modelA <- update(tchol_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy

tchol_modelA$call 

qsave(tchol_modelA, "./inputs/exposure_distributions/AH_test/tchol_model.qs", preset = "high")

# 
 GAIC.table(tchol_modelA, 
#            #tchol_modelB,
            mod_min)

# Double check that the distribution is still a good one
tt <- chooseDist(tchol_modelA,
                 type = "realplus",
                 trace = TRUE, data = ds,
                 parallel = "multicore", ncpus = 15L
)
qsave(tt, "./inputs/exposure_distributions/AH_test/new_distr_tchol.qs", preset = "high")
# # tt <- qread("./inputs/exposure_distributions/AH_test/new_distr.qs")
# 
which(tt == min(tt[, 1], na.rm = TRUE), arr.ind = TRUE) # Best distr based on AIC #GB2
which(tt == min(tt[, 3], na.rm = TRUE), arr.ind = TRUE) # Best distr based on BIC #GB2


}

if (update){
tchol_modelA <- qread("./inputs/exposure_distributions/AH_test/tchol_model.qs")
  
tchol_modelbest <- update(tchol_modelA, family = "GB2", method = mixed(50,400))
warnings()
GAIC.table(tchol_modelbest, tchol_modelA)
# Actually tchol_modelA is better 
tchol_modelbest <- copy(tchol_modelA)
qsave(tchol_modelbest, "./inputs/exposure_distributions/AH_test/tchol_modelbest.qs", preset = "high")
tchol_modelbest$hsedata <- copy(ds)
qsave(tchol_modelbest,"./secure_data/lifecourse_models/tchol_model.qs" )



# Code to create a table with predictions
trms <- all.vars(formula(tchol_modelbest))[-1] # -1 excludes dependent var
newdata <- CJ(
  year = 3:50, 
  age = 20:90, 
  sex = unique(ds$sex), 
  qimd = unique(ds$qimd),
  sha = unique(ds$sha), 
  smok_status = unique(ds$smok_status), 
  ethnicity_grp = unique(ds$ethnicity_grp),
  urban_rural = unique(ds$urban_rural),
  bmi = seq(18, 50, 1) 
)
# newdata holds all the possible combinations of predictors

# This is to be able to parallelise
newdata <- split(newdata, by = "year")
newdata <- # assignment necessary! Copies of data.tables are happening
  lapply(
    newdata,
    function(x) {
      x[, (tchol_modelbest$parameters) := predictAll(tchol_modelbest, .SD, data = ds), .SDcols = trms]
    }
  )
newdata <- rbindlist(newdata) # bind the chunks back
View(head(newdata, 1e3))

write_fst(newdata, path = "./inputs/exposure_distributions/AH_test/tchol_table.fst", compress = 100L) # This is what the model use as an input
# newdata <- read_fst("./inputs/exposure_distributions/AH_test/tchol_table.fst", as.data.table = T)

print("Table saved")
}

if (diagnostics) {
  tchol_model <- qread("./inputs/exposure_distributions/AH_test/tchol_modelbest.qs")
  #tchol_model <- qread("./inputs/exposure_distributions/AH_test/tchol_modelA.qs")
  
  plot(tchol_model)
  
  wp(tchol_model) # detrended QQ-plot
  wp(resid = resid(tchol_model)) # equivalen to the one above
  wp(resid = resid(tchol_model), ylim.all = 80 * sqrt(1 / length(resid(tchol_model))))
  tt <- wp(tchol_model, xvar = ds$age, n.inter = 6)
  wp(tchol_model, xvar = ds$year, n.inter = 10)
  wp(resid = resid(tchol_model), xvar = ~ ds$qimd)
  
  dtop(tchol_model, xvar = ds$age)
  
  rqres.plot(tchol_model)
  rqres.plot(tchol_model, type = "QQ")
}

if (plots) {
  #tchol_model <- qread("./tchol_modelbest.qs")
  tchol_model <- qread("./inputs/exposure_distributions/AH_test/tchol_modelbest.qs")
  dir.create("./inputs/exposure_distributions/AH_test/validation_synthpop_models", FALSE)
  
  #source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  
  # Restricting BMI to predicted values
  ds[, bmi := round(bmi, digits = 0L)]
  ds <- ds[bmi >= 18 & bmi <= 50]
  
  xlab_nam <- expression(bold(tchol ~ (mmHg)))
  tchol_model_tbl <- read_fst("./inputs/exposure_distributions/AH_test/tchol_table.fst", as.data.table = TRUE)
  zz <-
    validate_gamlss_tbl(ds[age > 19 & bmi <=50 & bmi >= 18], tchol_model_tbl, 50, "tchol", paste0("q", tchol_model$family[1]))[
      between(tchol, quantile(tchol, 0.01), quantile(tchol, 0.99))
    ]
  zz[, weight := wt_blood / sum(wt_blood), by = "type"]
  
  png(
    "./inputs/exposure_distributions/AH_test/validation_synthpop_models/tchol_rel_dist.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", tchol],
                      zz[type == "Modelled", tchol],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      100
  )
  dev.off()
  
  plot_synthpop_val(zz, tchol, "agegrp10", "wt_blood", "Total cholesterol by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, "year", "wt_blood", "Total cholesterol by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, "qimd", "wt_blood", "Total cholesterol by QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, "sha", "wt_blood", "Total cholesterol by SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, "ethnicity_grp", "wt_blood", "Total cholesterol by ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, "urban_rural", "wt_blood", "Total cholesterol by rurality", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, "smok_status", "wt_blood", "Total cholesterol by smoking status", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, c("year", "agegrp10"), "wt_blood", "Total cholesterol by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, c("year", "qimd"), "wt_blood", "Total cholesterol by year and QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, c("year", "sha"), "wt_blood", "Total cholesterol by year and SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, c("year", "ethnicity_grp"), "wt_blood", "Total cholesterol by year and ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, tchol, c("year", "urban_rural"), "wt_blood", "Total cholesterol by year and rurality", xlab_nam, FALSE, FALSE)
}
