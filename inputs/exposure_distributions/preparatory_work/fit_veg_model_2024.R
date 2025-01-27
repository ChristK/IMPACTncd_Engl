
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
fit_model <- TRUE
update <- FALSE
diagnostics <- FALSE
plots <- FALSE
seed                <- 43L



# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)

ds <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90) & vegpor <= 30 ,
                     .(vegpor, year, age, agegrp10, sex, qimd, ethnicity_grp, sha, wt_int, urban_rural)])

set.seed(seed)
# lns <- sample(nrow(ds), round(nrow(ds) * 0.8))
# ds_trn   <- ds[lns] # train dataset
# ds_crv   <- ds[!lns]  # cross-validation dataset




ds[, vegpor := round(vegpor)]

if(distributions){
  # First identify what distribution is best fit for vegpor (the left hand side
  # variable)
  # For count variables, need to round 
  ds[, vegpor := round(vegpor)]
  
  plot(hist(ds$vegpor, breaks = 100 ))
  plot(density(ds$vegpor))
  summary(ds$vegpor)
  # This is a count distribution, but 0-inflated
  
  # Let's see what distribution work best for this. We will use the function
  # fitDist from gamlss - this fits all the relevant distributions as specified 
  # by 'type' (see the documentation) and then orders them by some criteria (AIC is the default) 
  # Please have a read of the documentation for this function
  
  # Veg are count variables 
  marg_distr <- fitDist(vegpor,
                        k = 2, # Default is for AIC. use log(length(ds$vegpor)) for BIC
                        type = "count", # Other options "realline", "realplus", "real0to1", "counts", "binom"
                        try.gamlss = TRUE, # Better but slow
                        extra = NULL, # could put the extra distributions here + ZAIG
                        data = ds,
                        trace = TRUE, 
                        weights = ds$wt_int
  )
  
  # Check the best 10 distr - look at their properties, how many parameters? 
  head(marg_distr$fits, 10)
  #         LG     ZIPF      DEL ZASICHEL    ZAPIG    ZABNB       SI   SICHEL ZISICHEL    ZANBI 
  #    264690.7 290395.5 373445.1 374086.8 374595.0 374662.9 374855.7 374855.7 374857.7 374934.9 
  # need to check whether they crashed 
  

  
  # This function does "distribution validation diagnostics based on a fitted distribution model"
  # It is from the IMPACRncdEngl package, so if you don't have this loaded it won't work.
  # It also takes a long time
  # distr_validation(marg_distr, ds[between(vegpor, 0, 15), .(var = vegpor, wt = wt_int)],
  #                  expression(bold(vegpor)))
  
  # lets plot some of them, starting with the best AIC 'marg_distr$fits[1]'
  m1 <- gamlss(vegpor ~ 1,
               family = names(marg_distr$fits[1]),
               weights = ds$wt_int,
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
  
  plot(hist(ds$vegpor, breaks = 100), xlim = c(0,15))
  plot(density(ds$vegpor), lwd = 3, ylim = c(0, 0.7))
  lines(density(sample_from_gamlss(m_zaga)), col = alpha("red", 0.4))
  
  # The 2nd best 
  m2 <- gamlss(vegpor ~ 1, family = names(marg_distr$fits[2]), data = ds)
  lines(density(sample_from_gamlss(m2)), col = alpha("blue", 0.4))
  
  # The 3rd best
  m3 <- gamlss(vegpor ~ 1, family = names(marg_distr$fits[3]), data = ds)
  lines(density(sample_from_gamlss(m3)), col = alpha("green", 0.4))
  hist(sample_from_gamlss(m3))
  # they all look terrible 
  
  # The 4th best
  m4 <- gamlss(vegpor ~ 1, family = names(marg_distr$fits[4]), data = ds)
  lines(density(sample_from_gamlss(m4)), col = alpha("purple", 0.4))

  #  ZANBI as it's faster (number 10)
  m5 <- gamlss(vegpor ~ 1, family = names(marg_distr$fits[10]), data = ds)
  lines(density(sample_from_gamlss(m5)), col = alpha("purple", 0.4))
  hist(sample_from_gamlss(m5))
  # they all look terrible 
  
  # An alternative plot (skip these if the reldist package doesn't work)
  reldist(sample_from_gamlss(m1), ds$vegpor, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m2), ds$vegpor, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m3), ds$vegpor, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m4), ds$vegpor, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m_zaga), ds$vegpor, method = "bgk", bar = TRUE, show = "effect")
  
  
  
  
  distr_nam <- names(marg_distr$fits[10]) # I'm going for ZANBI as it's fast  
  
}

distr_nam <- "ZANBI"

print(distr_nam)
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3

# Univariate analysis ----------------------------------------
if (univariable_analysis) {
  ds[, .(vegpor_median = wtd.quantile(vegpor, weight = wt_int)), keyby = .(age)
  ][, scatter.smooth(age, vegpor_median)]
  
  
  m_age0 <- gamlss(
    vegpor ~ age,
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)  # Number of iterations using Rigby and Stasinopoulos algorith (default), 
    # followed by number using Cole and Green algorithm 
  )
  lines(centiles.pred(m_age0, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "purple")
  
  m_age1 <- gamlss(
    vegpor ~ pb(age), # penalised beta spline - an automated way for finding the splines & knots
    family = distr_nam,
    weights = ds$wt_int,
    data = ds#, #,
   # method = mixed(5, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")
  
  m_age2 <- gamlss(
    vegpor ~ poly(age, degree = 3),#polynomial
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")
  
  # m_age3 <- gamlss(
  #   vegpor ~ cs(age),
  #   family = distr_nam,
  #   weights = ds$wt_int,
  #   data = ds,
  #   method = mixed(20, 20)
  # )
  # lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")
#  GAIC.table(m_age0, m_age1, m_age2, m_age3)
  GAIC.table(m_age0, m_age1, m_age2) #m_age1 the best pb()
  print(GAIC.table(m_age0, m_age1, m_age2))
  
  #This is another way of looking at the different models 
  centiles(m_age1, xvar = ds$age)
  centiles(m_age2, xvar = ds$age)
  
  # m_age1 - pb spline is the best
  
  
  ds[, .(vegpor_median = wtd.quantile(vegpor, weight = wt_int)), keyby = .(year)][
    , scatter.smooth(year, vegpor_median, xlim = c(3, 40))]
  
  
  # Year
  m_year0 <- gamlss(
    vegpor ~ year,
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year0, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "blue1")
  
  m_year1 <- gamlss(
    vegpor ~ log(year),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "red1")
  
  m_year2 <- gamlss(
    vegpor ~ log(year + 10),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "green1")
  
  m_year3 <- gamlss( 
    vegpor ~ log(year + 100),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "pink")
  
  GAIC.table(m_year0, m_year1, m_year2, m_year3) #m_year3 is the best 
  centiles(m_year0, xvar = ds$year)
  centiles(m_year1, xvar = ds$year)
  centiles(m_year2, xvar = ds$year)
  centiles(m_year3, xvar = ds$year)
  
}

if(fit_model){
  
  ds[, vegpor := round(vegpor)]
  
 #  mod_min <- gamlss(
 #    vegpor ~ log(year + 100)   + pb(age) + pcat(qimd) +
 #      sex  + pcat(ethnicity_grp) + pcat(sha) + urban_rural  ,
 #    family = distr_nam,
 #    weights = ds$wt_int,
 #    data = ds,
 #    method = RS(400),
 #    control = con1 # for speed but looses accuracy if not 1e-3
 #  )
 #  qsave(mod_min, "./inputs/exposure_distributions/AH_test/vegpor_mod_min.qs", preset = "high")
 #   warnings() # check for warnings
 # 
 # 
 #   # Stepwise model selection using a Generalized Akaike Information Criterion
 # # mod_min <- qread( "./inputs/exposure_distributions/AH_test/vegpor_mod_min.qs")
 #  vegpor_modelA <- stepGAICAll.A(
 #    mod_min,
 #    scope = list(
 #      lower = ~ log(year + 100) +  pb(age) + pcat(qimd) +
 #        sex +  pcat(ethnicity_grp) + pcat(sha) + urban_rural,
 #      upper = ~ (year + pb(age) + pcat(qimd) +
 #                   sex +  pcat(ethnicity_grp) + pcat(sha) + urban_rural  )^2
 #    ),
 #    sigma.scope = list(
 #      lower = ~1,
 #      upper = ~ year + pb(age) + pcat(qimd) +
 #        sex + pcat(ethnicity_grp) + pcat(sha) + urban_rural
 #    ),
 #    nu.scope = list(
 #      lower = ~1,
 #      upper = ~ year + pb(age) + pcat(qimd) +
 #        sex +  pcat(ethnicity_grp) + pcat(sha) + urban_rural
 #    ),
 #    parallel = "multicore",
 #    ncpus = 12L,
 #    weights = wt_int,
 #    trace = TRUE  # This is so you can see all the models it has tried to fit
 #  )
 #  warnings()
 # 
 #  vegpor_modelA
 # 
 #  vegpor_modelA <- update(vegpor_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy
 # 
 #  vegpor_modelA$call
 # 
 #  qsave(vegpor_modelA, "./inputs/exposure_distributions/AH_test/vegpor_model.qs", preset = "high")
 #   # 
 #   # 
 #   GAIC.table(vegpor_modelA,
 #              #vegpor_modelB,
 #              mod_min)
 #   # 
 #   # # # Double check that the distribution is still a good one
 #   tt <- chooseDist(vegpor_modelA,
 #                    type = "count",
 #                    trace = TRUE, data = ds,
 #                    parallel = "multicore", ncpus = 15L
 #   )
 #   qsave(tt, "./inputs/exposure_distributions/AH_test/new_distr.qs", preset = "high")
 #   # # tt <- qread("./inputs/exposure_distributions/AH_test/new_distr.qs")
 #   #
 #   which(tt == min(tt[, 1], na.rm = TRUE), arr.ind = TRUE) # Best distr based on AIC
 #   which(tt == min(tt[, 3], na.rm = TRUE), arr.ind = TRUE) # Best distr based on BIC

 #  This is saying ZAP. 2 parameter distribution
   # 
   
   
   ds[, vegpor := round(vegpor)]
   distr_nam <- "ZAP"
   
   mod_min <- gamlss(
     vegpor ~ log(year + 100) , #  + pb(age) + pcat(qimd) + sex   ,
     family = distr_nam,
     weights = ds$wt_int,
     data = ds,
     method = RS(10),
     control = con1 # for speed but looses accuracy if not 1e-3
   )
   warnings()
   
   
   # Stepwise model selection using a Generalized Akaike Information Criterion
   
   vegpor_modelB <- stepGAICAll.A(
     mod_min,
     scope = list(
       lower = ~ log(year + 100),
       upper = ~ (year + pb(age) + pcat(qimd) +
                    sex +  pcat(ethnicity_grp) + pcat(sha) + urban_rural  )^2
     ),
     sigma.scope = list(
       lower = ~1,
       upper = ~ year + pb(age) + pcat(qimd) +
         sex + pcat(ethnicity_grp) + pcat(sha) + urban_rural
     ),
     parallel = "multicore",
     ncpus = 12L,
     weights = wt_int,
     trace = TRUE  # This is so you can see all the models it has tried to fit
   )
   warnings()
   
   vegpor_modelB
   
   vegpor_modelB <- update(vegpor_modelB, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy
   
   
   qsave(vegpor_modelB, "./vegpor_modelZAP2.qs", preset = "high")
   
vegpor_modelA <- qread("./inputs/exposure_distributions/AH_test/vegpor_model.qs")
vegpor_modelZAP1 <- qread("./inputs/exposure_distributions/AH_test/vegpor_model_ZAP.qs")

GAIC.table(vegpor_modelA,
           vegpor_modelZAP1  ,
           vegpor_modelB)

     # # # Double check that the distribution is still a good one
     tt <- chooseDist(vegpor_modelB,
                      type = "count",
                      trace = TRUE, data = ds,
                      parallel = "multicore", ncpus = 15L
     )
     qsave(tt, "./inputs/exposure_distributions/AH_test/new_distr_ZAP.qs", preset = "high")
     # # tt <- qread("./inputs/exposure_distributions/AH_test/new_distr.qs")
     #
     which(tt == min(tt[, 1], na.rm = TRUE), arr.ind = TRUE) # Best distr based on AIC
     which(tt == min(tt[, 3], na.rm = TRUE), arr.ind = TRUE) # Best distr based on BIC


}


if (update){
  # vegpor_modelA <- qread("./inputs/exposure_distributions/AH_test/vegpor_model.qs")
  # # 
  # vegpor_modelbest <- update(vegpor_modelA, family = "ZAP")
  # # warnings()
  # GAIC.table(vegpor_modelbest, vegpor_modelA)
  # # warnings()
  # qsave(vegpor_modelbest, "./inputs/exposure_distributions/AH_test/vegpor_model_ZAP.qs", preset = "high")
  # #vegpor_modelbest <- qread("./vegpor_modelbest.qs")
  # 
  
  distr_nam <- "ZAP"
vegpor_modelbest <- qread("./inputs/exposure_distributions/AH_test/vegpor_model_ZAP.qs")
vegpor_modelbest <- update(vegpor_modelbest, method = RS(600))
warnings()
qsave(vegpor_modelbest, "./inputs/exposure_distributions/AH_test/vegpor_model_ZAP600.qs")

  #vegpor_modelbest <- copy(vegpor_modelA )
  # Code to create a table with predictions
  trms <- all.vars(formula(vegpor_modelbest))[-1] # -1 excludes dependent var
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
        x[, (vegpor_modelbest$parameters) := predictAll(vegpor_modelbest, .SD, data = ds), .SDcols = trms]
      }
    )
  newdata <- rbindlist(newdata) # bind the chunks back
  #View(head(newdata, 1e3))
  setattr(newdata, "distribution", distr_nam)
  
  write_fst(newdata, path = "./inputs/exposure_distributions/AH_test/vegpor_table_ZAP.fst", compress = 100L) # This is what the model use as an input
  # newdata <- read_fst("./inputs/exposure_distributions/AH_test/vegpor_table.fst", as.data.table = T)

  print("Table saved")

  
}

if (diagnostics) {
  #vegpor_model <- qread("./inputs/exposure_distributions/AH_test/vegpor_modelbest.qs")
 # vegpor_model <- qread("./inputs/exposure_distributions/AH_test/vegpor_model.qs")
 vegpor_model <- qread("./inputs/exposure_distributions/AH_test/vegpor_model_ZAP.qs")
  
  plot(vegpor_model)
  
  wp(vegpor_model) # detrended QQ-plot
  wp(resid = resid(vegpor_model)) # equivalen to the one above
  wp(resid = resid(vegpor_model), ylim.all = 80 * sqrt(1 / length(resid(vegpor_model))))
  tt <- wp(vegpor_model, xvar = ds$age, n.inter = 6)
  wp(vegpor_model, xvar = ds$year, n.inter = 10)
  wp(resid = resid(vegpor_model), xvar = ~ ds$qimd)
  
  dtop(vegpor_model, xvar = ds$age)
  
  rqres.plot(vegpor_model)
  rqres.plot(vegpor_model, type = "QQ")
}

if (plots) {
  #vegpor_model <- qread("./vegpor_modelbest.qs")
  #vegpor_model <- qread("./inputs/exposure_distributions/AH_test/vegpor_modelbest.qs")
 # dir.create("./inputs/exposure_distributions/AH_test/validation_synthpop_models", FALSE)
  
  #source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  
  
  xlab_nam <- expression(bold(vegpor ~ (mmHg)))
  vegpor_model_tbl <- read_fst("./inputs/exposure_distributions/AH_test/vegpor_table_ZAP.fst", as.data.table = TRUE)
  zz <-
    validate_gamlss_tbl(ds[age > 19], vegpor_model_tbl, 50, "vegpor", paste0("q", vegpor_model$family[1]))[
      between(vegpor, quantile(vegpor, 0.01), quantile(vegpor, 0.99))
    ]
  zz[, weight := wt_int / sum(wt_int), by = "type"]
  
  png(
    "./inputs/exposure_distributions/AH_test/validation_synthpop_models/vegpor_rel_dist_ZAP.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", vegpor],
                      zz[type == "Modelled", vegpor],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      100
  )
  dev.off()
  
  plot_synthpop_val(zz, vegpor, "agegrp10", "wt_int", "Veg portions by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, "year", "wt_int", "Veg portions  by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, "qimd", "wt_int", "Veg portions  by QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, "sha", "wt_int", "Veg portions  by SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, "ethnicity_grp", "wt_int", "Veg portions  by ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, "urban_rural", "wt_int", "Veg portions  by rurality", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, c("year", "agegrp10"), "wt_int", "Veg portions  by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, c("year", "qimd"), "wt_int", "Veg portions  by year and QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, c("year", "sha"), "wt_int", "Veg portions  by year and SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, c("year", "ethnicity_grp"), "wt_int", "Veg portions  by year and ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, vegpor, c("year", "urban_rural"), "wt_int", "Veg portions  by year and rurality", xlab_nam, FALSE, FALSE)
}
