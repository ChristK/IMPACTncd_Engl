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
diagnostics <- TRUE
plots <- TRUE
seed                <- 43L

print("bmi")

# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)

# Complete case analysis for the variables fitting in the model, including 
# a valid 'nurse weight' (this is calculated by HSE to account for non-response rates) 
ds <- na.omit(HSE_ts[wt_nurse > 0 & age >= 20,
                     .(year, age, agegrp10, sex, qimd, ethnicity_grp, sha, smok_status, wt_nurse, urban_rural, bmi)])
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
  plot(density(ds$bmi))
  summary(ds$bmi)
  # This is a continuous distribution that cannot be 0
  
  # Let's see what distribution work best for this. We will use the function
  # fitDist from gamlss - this fits all the relevant distributions as specified 
  # by 'type' (see the documentation) and then orders them by some criteria (AIC is the default) 
  # Please have a read of the documentation for this function
  marg_distr <- fitDist(bmi,
                        k = 2, # Default is for AIC. use log(length(ds$bmi)) for BIC
                        type = "realAll", # Other options "realline", "realplus", "real0to1", "counts", "binom"
                        try.gamlss = TRUE, # Better but slow
                        extra = NULL,
                        data = ds,
                        trace = TRUE, 
                        weights = ds$wt_nurse
  )
  
  # Check the best 10 distr - look at their properties, how many parameters? 
  head(marg_distr$fits, 10)
  #     EGB2      JSU     JSUo      ST5      GB2     SEP4    BCPEo     BCPE     BCTo      BCT 
  #   529258.7 529264.9 529264.9 529270.6 529272.3 529280.7 529281.3 529281.3 529285.7 529285.7 
  # 
  # This function does "distribution validation diagnostics based on a fitted distribution model"
  # It is from the IMPACRncdEngl package, so if you don't have this loaded it won't work.
  # It also takes a long time
  distr_validation(marg_distr, ds[between(bmi, 10, 85), .(var = bmi, wt = wt_nurse)],
                   expression(bold(BMI)))
  
  # lets plot some of them, starting with the best AIC 'marg_distr$fits[1]'
  m1 <- gamlss(bmi ~ 1,
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
  
  plot(density(ds$bmi), lwd = 3, ylim = c(0, 0.1))
  # The 1st - plots weird
  lines(density(sample_from_gamlss(m1)), col = alpha("red", 0.4))

  # The 2nd best 
  m2 <- gamlss(bmi ~ 1, family = names(marg_distr$fits[2]), data = ds)
  lines(density(sample_from_gamlss(m2)), col = alpha("pink", 0.4))
  
  # The 4th best 
  m3 <- gamlss(bmi ~ 1, family = names(marg_distr$fits[4]), data = ds)
  lines(density(sample_from_gamlss(m2)), col = alpha("blue", 0.4))
  
  # The 20th best
  m4 <- gamlss(bmi ~ 1, family = names(marg_distr$fits[20]), data = ds)
  lines(density(sample_from_gamlss(m3)), col = alpha("green", 0.4))
  # they all look terrible 
  
  # An alternative plot (skip these if the reldist package doesn't work)
  reldist(sample_from_gamlss(m1), ds$bmi, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m2), ds$bmi, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m3), ds$bmi, method = "bgk", bar = TRUE, show = "effect")
  
  # What would the plot looked like if the distribution was off
  reldist(sample_from_gamlss(m1) + 1, ds$bmi, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m1) - 1, ds$bmi, method = "bgk", bar = TRUE, show = "effect")
  
  
  distr_nam <- names(marg_distr$fits[2]) # I will pick JSU (the second)
  
}

distr_nam <- "JSU"

print(distr_nam)
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3


if (univariable_analysis) {
  
  # Wtd.quantile function -----
  # If you cannot get the reldist package to load, try using this 
  # wtd.quantile <- function (x, q = 0.5, na.rm = FALSE, weight = NULL)
  # {
  #   if (mode(x) != "numeric")
  #     stop("need numeric data")
  #   if (!length(weight)) {
  #     quantile(x = x, probs = q, na.rm = na.rm)
  #   }
  #   else {
  #     x <- as.vector(x)
  #     wnas <- is.na(x)
  #     if (sum(wnas) > 0) {
  #       if (na.rm) {
  #         x <- x[!wnas]
  #         weight <- weight[!wnas]
  #       }
  #       else {
  #         return(NA)
  #       }
  #     }
  #     o <- order(x)
  #     n <- length(weight)
  #     order <- 1 + (n - 1) * q
  #     low <- pmax(floor(order), 1)
  #     high <- pmin(low + 1, n)
  #     order <- order%%1
  #     allq <- approx(x = cumsum(weight[o])/sum(weight), y = x[o],
  #                    xout = c(low, high)/n, method = "constant", f = 1,
  #                    rule = 2)$y
  #     k <- length(q)
  #     (1 - order) * allq[1:k] + order * allq[-(1:k)]
  #   }
  # }
  
  # ----
  # Age
  ds[, .(bmi_median = wtd.quantile(bmi, weight = wt_nurse)), keyby = .(age)][, scatter.smooth(age, bmi_median)]
  
  m_age0 <- gamlss(
    bmi ~ age,
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)  # Number of iterations using Rigby and Stasinopoulos algorith (default), 
    # followed by number using Cole and Green algorithm 
  )
  lines(centiles.pred(m_age0, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "purple")
  
  m_age1 <- gamlss(
    bmi ~ pb(age), # penalised beta spline - an automated way for finding the splines & knots
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")
  
  m_age2 <- gamlss(
    bmi ~ poly(age, degree = 3),#polynomial
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")
  
  m_age3 <- gamlss(
    bmi ~ cs(age),
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
  
  
  ds[, .(bmi_median = wtd.quantile(bmi, weight = wt_nurse)), keyby = .(year)][, scatter.smooth(year, bmi_median, xlim = c(3, 40), ylim = c(20, 30))]
  
  # Year
  m_year1 <- gamlss(
    bmi ~ year,
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "blue1")
  
  m_year2 <- gamlss(
    bmi ~ log(year),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "red1")
  
  m_year3 <- gamlss(
    bmi ~ log(year + 10),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "green1")
  
  m_year4 <- gamlss(
    bmi ~ log(year + 100),
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year4, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "pink")
  
  GAIC.table(m_year1, m_year2, m_year3, m_year4) #m_year2 is the best 
  centiles(m_year1, xvar = ds$year)
  centiles(m_year2, xvar = ds$year)
  centiles(m_year3, xvar = ds$year)
  
  
  # Smoking status 
  m_sm1 <- gamlss(
    bmi ~ smok_status,
    family = distr_nam,
    weights = ds$wt_nurse,
    data = ds,
    method = mixed(50, 20)
  )
  summary(m_sm1)
  
  # Using a smoother 
  m_sm2 <- gamlss(
    bmi ~ pcat(smok_status, method = "GAIC"),
    family = distr_nam,
    weights = ds$wt_nurse,
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
  
  #In theory you can plot this, but it seems to take an incredibly long tine 
  # plotLambda(sbp,factor=smok_status,data=ds,family=JSU,along=seq(-10,2,.2))
  # abline(v=log(getSmo(m_sm2)$lambda),col='gray',lwd=2.5)
  
 
  
  
  
}


run <- FALSE 
if(run){
# Setting a minimum model - we know we want these things
mod_min <- gamlss(
  bmi ~ log(year) + pb(age) + pcat(qimd) +
    sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural,
  family = distr_nam,
  weights = ds$wt_nurse,
  data = ds,
  method = mixed(5, 100),
  control = con1 # for speed but looses accuracy if not 1e-3
)
qsave(mod_min, "./inputs/exposure_distributions/AH_test/bmi_mod_min.qs", preset = "high")
warnings() # check for warnings
mod_min <- qread( "./inputs/exposure_distributions/AH_test/bmi_mod_min.qs")

# Stepwise model selection using a Generalized Akaike Information Criterion
bmi_modelA <- stepGAICAll.A(
  mod_min,
  scope = list(
    lower = ~ bmi ~ log(year) + pb(age) + pcat(qimd) +
    sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural,
    upper = ~ (log(year) + pb(age) + pcat(qimd) +
                 sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural)^2
  ),
  sigma.scope = list(
    lower = ~1,
    upper = ~ log(year) + pb(age) + pcat(qimd) +
      sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural
  ),
  nu.scope = list(
    lower = ~1,
    upper = ~ log(year) + pb(age) + pcat(qimd) +
      sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural 
  ),
  tau.scope = list(
    lower = ~1,
    upper = ~ log(year) + pb(age) + pcat(qimd) +
      sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural 
  ),
  parallel = "multicore",
  ncpus = 12L,
  weights = wt_nurse,
  trace = TRUE  # This is so you can see all the models it has tried to fit
)
warnings()

bmi_modelA

bmi_modelA <- update(bmi_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy

bmi_modelA$call


qsave(bmi_modelA, "./inputs/exposure_distributions/AH_test/bmi_modelAwt_July.qs", preset = "high")
#bmi_modelA <- qread("./inputs/exposure_distributions/AH_test/bmi_modelAwt_July.qs")


GAIC.table(bmi_modelA,
           #bmi_modelB,
           mod_min)

# Double check that the distribution is still a good one
tt <- chooseDist(bmi_modelA,
                 type = "realplus",
                 trace = TRUE, data = ds,
                 parallel = "multicore", ncpus = 15L
)
#qsave(tt, "./inputs/exposure_distributions/AH_test/new_distr2.qs", preset = "high")
# # tt <- qread("./inputs/exposure_distributions/AH_test/new_distr.qs")
#
which(tt == min(tt[, 1], na.rm = TRUE), arr.ind = TRUE) # Best distr based on AIC
which(tt == min(tt[, 3], na.rm = TRUE), arr.ind = TRUE) # Best distr based on BIC
# 
bmi_modelbest <- update(bmi_modelA, family = "GB2") 
GAIC.table(bmi_modelbest, bmi_modelA, mod_min )
}
distr_nam <- "GB2"

bmi_modelbest <- qread("./inputs/exposure_distributions/AH_test/bmi_model.qs")
bmi_modelbest <- update(bmi_modelbest, method = mixed(100, 400))
warnings()
qsave(bmi_modelbest, "./inputs/exposure_distributions/AH_test/bmi_modelbest.qs", preset = "high")

bmi_modelbest$hsedata <- copy(ds)
qsave(bmi_modelbest,"./secure_data/lifecourse_models/bmi_model.qs" )


# # 
# 
# # Code to create a table with predictions
trms <- all.vars(formula(bmi_modelbest))[-1] # -1 excludes dependent var
newdata <- CJ(
  year = 3:50,
  age = 20:90,
  sex = unique(ds$sex),
  qimd = unique(ds$qimd),
  sha = unique(ds$sha),
  smok_status = unique(ds$smok_status),
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
      x[, (bmi_modelbest$parameters) := predictAll(bmi_modelbest, .SD, data = ds), .SDcols = trms]
    }
  )
newdata <- rbindlist(newdata) # bind the chunks back
View(head(newdata, 1e3))

write_fst(newdata, path = "./inputs/exposure_distributions/AH_test/bmi_table.fst", compress = 100L) # This is what the model use as an input
# newdata <- read_fst("./inputs/exposure_distributions/AH_test/bmi_table.fst", as.data.table = T)
# 
# print("Table saved")

if (diagnostics) {
  #bmi_model <- qread("./bmi_modelbest.qs")
  bmi_model <- qread("./inputs/exposure_distributions/AH_test/bmi_modelbest.qs")
  
  plot(bmi_model) 
  
  wp(bmi_model) # detrended QQ-plot
  wp(resid = resid(bmi_model)) # equivalen to the one above
  wp(resid = resid(bmi_model), ylim.all = 80 * sqrt(1 / length(resid(bmi_model))))
  tt <- wp(bmi_model, xvar = ds$age, n.inter = 6)
  wp(bmi_model, xvar = ds$year, n.inter = 10) 
  wp(resid = resid(bmi_model), xvar = ~ ds$qimd)
  
  dtop(bmi_model, xvar = ds$age)
  
  rqres.plot(bmi_model)
  rqres.plot(bmi_model, type = "QQ")
}

if (plots) {
  #bmi_model <- qread("./bmi_modelbest.qs")
  bmi_model <- qread("./inputs/exposure_distributions/AH_test/bmi_modelbest.qs")
  dir.create("./inputs/exposure_distributions/AH_test/validation_synthpop_models", FALSE)
  
  #source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  
  
  xlab_nam <- expression(bold(BMI))
  bmi_model_tbl <- read_fst("./inputs/exposure_distributions/AH_test/bmi_table.fst", as.data.table = TRUE)
  zz <-
    validate_gamlss_tbl(ds[age > 19], bmi_model_tbl, 50, "bmi", paste0("q", bmi_model$family[1]))[
      between(bmi, quantile(bmi, 0.01), quantile(bmi, 0.99))
    ]
  zz[, weight := wt_nurse / sum(wt_nurse), by = "type"]
  
  png(
    "./inputs/exposure_distributions/AH_test/validation_synthpop_models/bmi_rel_dist.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", bmi],
                      zz[type == "Modelled", bmi],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      100
  )
  dev.off()
  
  plot_synthpop_val(zz, bmi, "agegrp10", "wt_nurse", "BMI by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, "year", "wt_nurse", "BMI by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, "qimd", "wt_nurse", "bmi by QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, "sha", "wt_nurse", "BMI by SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, "ethnicity_grp", "wt_nurse", "BMI by ethinicty", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, "smok_status", "wt_nurse", "BMI by smoking status", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, "urban_rural", "wt_nurse", "BMI by rurality", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, c("year", "agegrp10"), "wt_nurse", "BMI by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, c("year", "qimd"), "wt_nurse", "BMI by year and QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, c("year", "sha"), "wt_nurse", "BMI by year and SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, c("year", "ethnicity_grp"), "wt_nurse", "BMI by year and ethinicty", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, bmi, c("year", "urban_rural"), "wt_nurse", "BMI by year and rurality", xlab_nam, FALSE, FALSE)
  
  
}
