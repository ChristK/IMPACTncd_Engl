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

# Final distribution chosen: GB2 

library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk
library(IMPACTncdEngl)
library(CKutils)

# Set some variables we will use later in automation
distributions <- FALSE
univariable_analysis <- FALSE
fit_model <- FALSE
update <- TRUE
diagnostics <- TRUE
plots <- TRUE
seed                <- 43L


# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)

# Complete case analysis for the variables fitting in the model, including 
# a valid 'nurse weight' (this is calculated by HSE to account for non-response rates) 
ds <- na.omit(HSE_ts[wt_nurse > 0 & age >= 20,
                     .(sbp, year, age, agegrp10, sex, qimd, ethnicity_grp, sha, smok_status, wt_nurse, urban_rural, bmi)])
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
# First identify what distribution is best fit for SBP (the left hand side
# variable)
plot(density(ds$sbp))
summary(ds$sbp)
# This is a continuous distribution that cannot be 0

# Let's see what distribution work best for this. We will use the function
# fitDist from gamlss - this fits all the relevant distributions as specified 
# by 'type' (see the documentation) and then orders them by some criteria (AIC is the default) 
# Please have a read of the documentation for this function
marg_distr <- fitDist(sbp,
  k = 2, # Default is for AIC. use log(length(ds$sbp)) for BIC
  type = "realAll", # Other options "realline", "realplus", "real0to1", "counts", "binom"
  try.gamlss = TRUE, # Better but slow
  extra = NULL,
  data = ds,
  trace = TRUE, 
  weights = ds$wt_nurse
)

# Check the best 10 distr - look at their properties, how many parameters? 
head(marg_distr$fits, 10)
#   JSU     JSUo      ST5      GB2     SEP4     BCTo      BCT    BCPEo     BCPE      ST1 
# 635248.9 635248.9 635262.5 635263.3 635268.7 635273.8 635273.8 635279.9 635279.9 635312.7 
# 
# This function does "distribution validation diagnostics based on a fitted distribution model"
# It is from the IMPACRncdEngl package, so if you don't have this loaded it won't work.
# It also takes a long time
distr_validation(marg_distr, ds[between(sbp, 70, 240), .(var = sbp, wt = wt_nurse)],
                 expression(bold(SBP ~ (mmHg))))

# lets plot some of them, starting with the best AIC 'marg_distr$fits[1]'
m1 <- gamlss(sbp ~ 1,
  family = names(marg_distr$fits[1]),
  weights = ds$wt_nurse,
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

plot(density(ds$sbp), lwd = 3, ylim = c(0, 0.03))
lines(density(sample_from_gamlss(m1)), col = alpha("red", 0.4))

# The 4th best 
m2 <- gamlss(sbp ~ 1, family = names(marg_distr$fits[4]), data = ds)
lines(density(sample_from_gamlss(m2)), col = alpha("blue", 0.4))

# The 20th best
m3 <- gamlss(sbp ~ 1, family = names(marg_distr$fits[20]), data = ds)
lines(density(sample_from_gamlss(m3)), col = alpha("green", 0.4))
# they all look terrible 

# An alternative plot (skip these if the reldist package doesn't work)
reldist(sample_from_gamlss(m1), ds$sbp, method = "bgk", bar = TRUE, show = "effect")
reldist(sample_from_gamlss(m3), ds$sbp, method = "bgk", bar = TRUE, show = "effect")

# What would the plot looked like if the distribution was off
reldist(sample_from_gamlss(m1) + 1, ds$sbp, method = "bgk", bar = TRUE, show = "effect")
reldist(sample_from_gamlss(m1) - 1, ds$sbp, method = "bgk", bar = TRUE, show = "effect")


distr_nam <- names(marg_distr$fits[1]) # I will pick JSU

}

distr_nam <- "JSU"

print(distr_nam)
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3


if (univariable_analysis) {

  # Wtd.quantile function -----
  # If you cannot get the reldist package to load, try using this 
  wtd.quantile <- function (x, q = 0.5, na.rm = FALSE, weight = NULL) 
  {
    if (mode(x) != "numeric") 
      stop("need numeric data")
    if (!length(weight)) {
      quantile(x = x, probs = q, na.rm = na.rm)
    }
    else {
      x <- as.vector(x)
      wnas <- is.na(x)
      if (sum(wnas) > 0) {
        if (na.rm) {
          x <- x[!wnas]
          weight <- weight[!wnas]
        }
        else {
          return(NA)
        }
      }
      o <- order(x)
      n <- length(weight)
      order <- 1 + (n - 1) * q
      low <- pmax(floor(order), 1)
      high <- pmin(low + 1, n)
      order <- order%%1
      allq <- approx(x = cumsum(weight[o])/sum(weight), y = x[o], 
                     xout = c(low, high)/n, method = "constant", f = 1, 
                     rule = 2)$y
      k <- length(q)
      (1 - order) * allq[1:k] + order * allq[-(1:k)]
    }
  }

  # ----
    # Age
  ds[, .(sbp_median = wtd.quantile(sbp, weight = wt_nurse)), keyby = .(age)][, scatter.smooth(age, sbp_median)]

  m_age0 <- gamlss(
    sbp ~ age,
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)  # Number of iterations using Rigby and Stasinopoulos algorith (default), 
    # followed by number using Cole and Green algorithm 
  )
  lines(centiles.pred(m_age0, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "purple")

  m_age1 <- gamlss(
    sbp ~ pb(age), # penalised beta spline - an automated way for finding the splines & knots
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")

  m_age2 <- gamlss(
    sbp ~ poly(age, degree = 3),#polynomial
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")

  m_age3 <- gamlss(
    sbp ~ cs(age),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")
  GAIC.table(m_age0, m_age1, m_age2, m_age3)
  
  #This is another way of looking at the different models 
  centiles(m_age1, xvar = ds$age)
  centiles(m_age2, xvar = ds$age)

  # m_age1 - pb spline is the best
  
  
  ds[, .(sbp_median = wtd.quantile(sbp, weight = wt_nurse)), keyby = .(year)][, scatter.smooth(year, sbp_median, xlim = c(3, 40), ylim = c(110, 130))]

  # Year
  m_year1 <- gamlss(
    sbp ~ year,
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "blue1")

  m_year2 <- gamlss(
    sbp ~ log(year),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "red1")

  m_year3 <- gamlss(
    sbp ~ log(year + 10),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "green1")
  
  m_year4 <- gamlss(
    sbp ~ log(year + 100),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year4, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "pink")

  GAIC.table(m_year1, m_year2, m_year3, m_year4) #m_year4 is the best 
  centiles(m_year1, xvar = ds$year)
  centiles(m_year2, xvar = ds$year)
  centiles(m_year3, xvar = ds$year)
  
  
  # Smoking status 
  m_sm1 <- gamlss(
    sbp ~ smok_status,
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  summary(m_sm1)
  
  # Using a smoother 
  m_sm2 <- gamlss(
    sbp ~ pcat(smok_status, method = "GAIC"),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  summary(m_sm2)
  
  GAIC.table(m_sm1, m_sm2)
  #without is better
  
  #Looking at the levels pcat has reduced smok_status to, and how they 
  #correspond to the original levels
  f2 <- getSmo(m_sm2)$factor
  levels(f2)
  table(f2,ds$smok_status)
  
  #In theory you can plot this, but it seems to take an incredibly long tine 
  # plotLambda(sbp,factor=smok_status,data=ds,family=JSU,along=seq(-10,2,.2))
  # abline(v=log(getSmo(m_sm2)$lambda),col='gray',lwd=2.5)
  
  # BMI
  ds[, .(bmi_median = wtd.quantile(bmi, weight = wt_nurse)), keyby = .(bmi)][, scatter.smooth(bmi, bmi_median)]
  
  m_bmi0 <- gamlss(
    sbp ~ bmi,
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_bmi0, xname = "bmi", xvalues = 20:90, cent = 50, data = ds), col = "purple")
  
  m_bmi1 <- gamlss(
    sbp ~ pb(bmi), 
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_bmi1, xname = "bmi", xvalues = 20:90, cent = 50, data = ds), col = "blue1")
  
  m_bmi2 <- gamlss(
    sbp ~ poly(bmi, 3),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_bmi2, xname = "bmi", xvalues = 20:90, cent = 50, data = ds), col = "red1")
  
  m_bmi3 <- gamlss(
    sbp ~ cs(bmi),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_bmi3, xname = "bmi", xvalues = 20:90, cent = 50, data = ds), col = "green1")
  GAIC.table(m_bmi0, m_bmi1, m_bmi2, m_bmi3)
  centiles(m_bmi1, xvar = ds$bmi)
  centiles(m_bmi2, xvar = ds$bmi)
  
  # m_bmi1 - pb spline is the best
  
  
  
  
}

if(fit_model){
# Setting a minimum model - we know we want these things - ideally we would have BMI included, but the model doesn't fit
# mod_min <- gamlss(
#   sbp ~ log(year + 100) + pb(age) + pcat(qimd) +
#     sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural + pb(bmi) ,
#   family = distr_nam,
#   weights = ds$wt_nurse,
#   data = ds,
#   method = mixed(5, 100),
#   control = con1 # for speed but looses accuracy if not 1e-3
# )
# qsave(mod_min, "./inputs/exposure_distributions/AH_test/sbp_mod_min.qs", preset = "high")
# warnings()

# Setting a minimum model - we know we want these things 
mod_min <- gamlss(
  sbp ~ log(year + 100) + pb(age) + pcat(qimd) +
    sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural,
  family = distr_nam,
  weights = ds$wt_nurse,
  data = ds,
  method = mixed(5, 100),
  control = con1 # for speed but looses accuracy if not 1e-3
)
qsave(mod_min, "./inputs/exposure_distributions/AH_test/sbp_mod_min_nobmi.qs", preset = "high")
warnings() # check for warnings
# mod_min <-qread( "./inputs/exposure_distributions/AH_test/sbp_mod_min_nobmi.qs")

# Stepwise model selection using a Generalized Akaike Information Criterion
sbp_modelA <- stepGAICAll.A(
  mod_min,
  scope = list(
    lower = ~ log(year + 100) + pb(age) + pcat(qimd) + 
      sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural,
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
  weights = wt_nurse, 
  trace = TRUE  # This is so you can see all the models it has tried to fit 
)
warnings()

sbp_modelA

sbp_modelA <- update(sbp_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy

sbp_modelA$call 


qsave(sbp_modelA, "./inputs/exposure_distributions/AH_test/sbp_modelAwt_Aug.qs", preset = "high")
#sbp_modelA <- qread("./inputs/exposure_distributions/AH_test/sbp_modelAwt_July.qs")


GAIC.table(sbp_modelA, 
           #sbp_modelB,
           mod_min)

# Double check that the distribution is still a good one
tt <- chooseDist(sbp_modelA,
  type = "realplus",
  trace = TRUE, data = ds,
  parallel = "multicore", ncpus = 15L
)
 qsave(tt, "./inputs/exposure_distributions/AH_test/new_distr2.qs", preset = "high")
# # tt <- qread("./inputs/exposure_distributions/AH_test/new_distr.qs")
# 
which(tt == min(tt[, 1], na.rm = TRUE), arr.ind = TRUE) # Best distr based on AIC
which(tt == min(tt[, 3], na.rm = TRUE), arr.ind = TRUE) # Best distr based on BIC
# best is GB2 

}
if(update){
  sbp_modelA <- qread("./inputs/exposure_distributions/AH_test/sbp_modelAwt_Aug.qs")
 sbp_modelbest <- update(sbp_modelA, family = "GB2") 
 GAIC.table(sbp_modelbest, sbp_modelA)
 qsave(sbp_modelbest, "./inputs/exposure_distributions/sbp_modelbest.qs", preset = "high")
 
 sbp_modelbest$hsedata <- copy(ds)
 
 qsave(sbp_modelbest,"./secure_data/lifecourse_models/sbp_model.qs" )
 
# sbp_modelbest <- qread("./inputs/exposure_distributions/sbp_modelbest.qs")
# 
# qsave(sbp_modelbest, "./sbp_modelbest.qs", preset = "high")

# Code to create a table with predictions
trms <- all.vars(formula(sbp_modelbest))[-1] # -1 excludes dependent var
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
      x[, (sbp_modelbest$parameters) := predictAll(sbp_modelbest, .SD, data = ds), .SDcols = trms]
    }
  )
newdata <- rbindlist(newdata) # bind the chunks back
View(head(newdata, 1e3))

write_fst(newdata, path = "./inputs/exposure_distributions/sbp_table.fst", compress = 100L) # This is what the model use as an input
# newdata <- read_fst("./inputs/exposure_distributions/sbp_table.fst", as.data.table = T)

print("Table saved")
}

if (diagnostics) {
  sbp_model <- qread("./inputs/exposure_distributions/sbp_modelbest.qs")

  plot(sbp_model)

  wp(sbp_model) # detrended QQ-plot
  wp(resid = resid(sbp_model)) # equivalen to the one above
  wp(resid = resid(sbp_model), ylim.all = 80 * sqrt(1 / length(resid(sbp_model))))
  tt <- wp(sbp_model, xvar = ds$age, n.inter = 6)
  wp(sbp_model, xvar = ds$year, n.inter = 10)
  wp(resid = resid(sbp_model), xvar = ~ ds$qimd)

  dtop(sbp_model, xvar = ds$age)

  rqres.plot(sbp_model)
  rqres.plot(sbp_model, type = "QQ")
}

if (plots) {
  #sbp_model <- qread("./sbp_modelbest.qs")
  sbp_model <- qread("./inputs/exposure_distributions/sbp_modelbest.qs")
 # dir.create("./inputs/exposure_distributions/AH_test/validation_synthpop_models", FALSE)

  #source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())

  # Restricting BMI to predicted values
  ds[, bmi := round(bmi, digits = 0L)]
  ds <- ds[bmi >= 18 & bmi <= 50]

  xlab_nam <- expression(bold(SBP ~ (mmHg)))
  sbp_model_tbl <- read_fst("./inputs/exposure_distributions/sbp_table.fst", as.data.table = TRUE)
  zz <-
    validate_gamlss_tbl(ds[age > 19], sbp_model_tbl, 50, "sbp", paste0("q", sbp_model$family[1]))[
      between(sbp, quantile(sbp, 0.01), quantile(sbp, 0.99))
    ]
  zz[, weight := wt_nurse / sum(wt_nurse), by = "type"]

  png(
    "./inputs/exposure_distributions/AH_test/validation_synthpop_models/SBP_rel_dist.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", sbp],
    zz[type == "Modelled", sbp],
    zz[type == "Observed", weight],
    zz[type == "Modelled", weight],
    main = xlab_nam,
    100
  )
  dev.off()

  plot_synthpop_val(zz, sbp, "agegrp10", "wt_nurse", "SBP by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, "year", "wt_nurse", "SBP by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, "qimd", "wt_nurse", "SBP by QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, "sha", "wt_nurse", "SBP by SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, "ethnicity_grp", "wt_nurse", "SBP by ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, "smok_status", "wt_nurse", "SBP by smoking status", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, c("year", "agegrp10"), "wt_nurse", "SBP by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, c("year", "qimd"), "wt_nurse", "SBP by year and QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, c("year", "sha"), "wt_nurse", "SBP by year and SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, sbp, c("year", "ethnicity_grp"), "wt_nurse", "SBP by year and ethnicity", xlab_nam, FALSE, FALSE)
}
