# Final distribution chosen: ZINBI 

print("fruit")

library(IMPACTncdEngl)
library(CKutils)
library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk

# Set some variables we will use later in automation
distributions <- FALSE
univariable_analysis <- FALSE
fit_model <- FALSE
update <- FALSE
diagnostics <- TRUE
plots <- TRUE
seed                <- 43L



# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)

ds <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90) & frtpor <= 30 ,
                     .(frtpor, year, age, agegrp10, sex, qimd, ethnicity_grp, sha, wt_int, urban_rural)])

set.seed(seed)
# lns <- sample(nrow(ds), round(nrow(ds) * 0.8))
# ds_trn   <- ds[lns] # train dataset
# ds_crv   <- ds[!lns]  # cross-validation dataset


ZA_models <- FALSE
if(ZA_models){
  plot(hist(ds$frtpor, breaks = 100 ))
  
  # Many of the count distributions are very slow 
  # Could try zero adjusted gamma 
  m_zaga <- gamlss(frtpor ~ 1,
                   family = "ZAGA",
                   weights = ds$wt_int,
                   # method = mixed(20, 20)
                   data = ds
  )
  # AIC = 450080 
  
  
  m_zaig <- gamlss(frtpor ~ 1,
                   family = "ZAIG",
                   weights = ds$wt_int,
                   # method = mixed(20, 20)
                   data = ds
  )
  # AIC = 446701 
  
  GAIC.table(m_zaga, m_zaig)
  # minimum GAIC(k= 2 ) model: m_zaig 
  # minimum GAIC(k= 3.84 ) model: m_zaig 
  # minimum GAIC(k= 11.69 ) model: m_zaig 
  # df      k=2   k=3.84  k=11.69
  # m_zaga  3 450079.9 450085.4 450109.0
  # m_zaig  3 446701.4 446706.9 446730.4
  # 
  
  
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
  
  plot(hist(ds$frtpor, breaks = 100 ), xlim = c(0,15))
  plot(hist(sample_from_gamlss(m_zaga), breaks = 100), xlim = c(0,15))
  plot(hist(sample_from_gamlss(m_zaig), breaks = 100), xlim = c(0,15))
  
  
}




if(distributions){
  # First identify what distribution is best fit for frtpor (the left hand side
  # variable)
  # For count variables, need to round 
  ds[, frtpor := round(frtpor)]
  
  plot(hist(ds$frtpor, breaks = 100 ))
  plot(density(ds$frtpor))
  summary(ds$frtpor)
  # This is a count distribution, but 0-inflated
  
  # Let's see what distribution work best for this. We will use the function
  # fitDist from gamlss - this fits all the relevant distributions as specified 
  # by 'type' (see the documentation) and then orders them by some criteria (AIC is the default) 
  # Please have a read of the documentation for this function
  
  # Fruit (& veg) are count variables 
  marg_distr <- fitDist(frtpor,
                        k = 2, # Default is for AIC. use log(length(ds$frtpor)) for BIC
                        type = "count", # Other options "realline", "realplus", "real0to1", "counts", "binom"
                        try.gamlss = TRUE, # Better but slow
                        extra = NULL, # could put the extra distributions here + ZAIG
                        data = ds,
                        trace = TRUE, 
                        weights = ds$wt_int
  )
  
  # Check the best 10 distr - look at their properties, how many parameters? 
  head(marg_distr$fits, 10)
  #    LG     ZIPF ZASICHEL ZISICHEL    ZIPIG    ZAPIG    ZIBNB    ZABNB    ZINBI    ZANBI 
  # 342835.5 379777.3 460866.7 460866.7 461365.2 461365.2 461429.1 461429.1 461828.1 461828.1 # 
  # need to check whether they crashed 
  

  
  # This function does "distribution validation diagnostics based on a fitted distribution model"
  # It is from the IMPACRncdEngl package, so if you don't have this loaded it won't work.
  # It also takes a long time
  # distr_validation(marg_distr, ds[between(frtpor, 0, 15), .(var = frtpor, wt = wt_int)],
  #                  expression(bold(frtpor)))
  
  # lets plot some of them, starting with the best AIC 'marg_distr$fits[1]'
  m1 <- gamlss(frtpor ~ 1,
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
  
  plot(hist(ds$frtpor, breaks = 100), xlim = c(0,15))
  plot(density(ds$frtpor), lwd = 3, ylim = c(0, 0.7))
  lines(density(sample_from_gamlss(m_zaga)), col = alpha("red", 0.4))
  
  # The 2nd best 
  m2 <- gamlss(frtpor ~ 1, family = names(marg_distr$fits[2]), data = ds)
  lines(density(sample_from_gamlss(m2)), col = alpha("blue", 0.4))
  
  # The 3rd best
  m3 <- gamlss(frtpor ~ 1, family = names(marg_distr$fits[3]), data = ds)
  lines(density(sample_from_gamlss(m3)), col = alpha("green", 0.4))
  hist(sample_from_gamlss(m3))
  # they all look terrible 
  
  # The 4th best
  m4 <- gamlss(frtpor ~ 1, family = names(marg_distr$fits[4]), data = ds)
  lines(density(sample_from_gamlss(m4)), col = alpha("purple", 0.4))

  #  ZINBI as it's faster (number 9)
  m5 <- gamlss(frtpor ~ 1, family = names(marg_distr$fits[9]), data = ds)
  lines(density(sample_from_gamlss(m5)), col = alpha("purple", 0.4))
  hist(sample_from_gamlss(m5))
  # they all look terrible 
  
  # An alternative plot (skip these if the reldist package doesn't work)
  reldist(sample_from_gamlss(m1), ds$frtpor, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m2), ds$frtpor, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m3), ds$frtpor, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m4), ds$frtpor, method = "bgk", bar = TRUE, show = "effect")
  reldist(sample_from_gamlss(m_zaga), ds$frtpor, method = "bgk", bar = TRUE, show = "effect")
  
  
  
  
  distr_nam <- names(marg_distr$fits[9]) # I'm going for ZINBI as it's fast  
  
}

distr_nam <- "ZAIG"

print(distr_nam)
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis, i.e. 1e-2. the default is 1e-3

# Univariate analysis ----------------------------------------
if (univariable_analysis) {
  ds[, .(frtpor_median = wtd.quantile(frtpor, weight = wt_int)), keyby = .(age)
  ][, scatter.smooth(age, frtpor_median)]
  
  
  m_age0 <- gamlss(
    frtpor ~ age,
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)  # Number of iterations using Rigby and Stasinopoulos algorith (default), 
    # followed by number using Cole and Green algorithm 
  )
  lines(centiles.pred(m_age0, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "purple")
  
  m_age1 <- gamlss(
    frtpor ~ pb(age), # penalised beta spline - an automated way for finding the splines & knots
    family = distr_nam,
    weights = ds$wt_int,
    data = ds#, #,
   # method = mixed(5, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")
  
  m_age2 <- gamlss(
    frtpor ~ poly(age, degree = 3),#polynomial
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")
  
  m_age3 <- gamlss(
    frtpor ~ cs(age),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")
#  GAIC.table(m_age0, m_age1, m_age2, m_age3)
  GAIC.table(m_age0, m_age1, m_age2)
  
  #This is another way of looking at the different models 
  centiles(m_age1, xvar = ds$age)
  centiles(m_age2, xvar = ds$age)
  
  # m_age1 - pb spline is the best
  
  
  ds[, .(frtpor_median = wtd.quantile(frtpor, weight = wt_int)), keyby = .(year)][
    , scatter.smooth(year, frtpor_median, xlim = c(3, 40))]
  
  
  # Year
  m_year0 <- gamlss(
    frtpor ~ year,
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year0, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "blue1")
  
  m_year1 <- gamlss(
    frtpor ~ log(year),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "red1")
  
  m_year2 <- gamlss(
    frtpor ~ log(year + 10),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "green1")
  
  m_year3 <- gamlss(
    frtpor ~ log(year + 100),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "pink")
  
  GAIC.table(m_year0, m_year1, m_year2, m_year3) #m_year0 
  centiles(m_year0, xvar = ds$year)
  centiles(m_year1, xvar = ds$year)
  centiles(m_year2, xvar = ds$year)
  centiles(m_year3, xvar = ds$year)
  
}

if(fit_model){
  
  ds[, frtpor := round(frtpor)]
  
  # mod_min <- gamlss(
  #   frtpor ~ year + pb(age) + pcat(qimd) +
  #     sex  + pcat(ethnicity_grp) + pcat(sha) + urban_rural  ,
  #   family = distr_nam,
  #   weights = ds$wt_int,
  #   data = ds,
  #   method = mixed(5, 100),
  #   control = con1 # for speed but looses accuracy if not 1e-3
  # )
  # qsave(mod_min, "./inputs/exposure_distributions/AH_test/frtpor_mod_min_zinbi.qs", preset = "high")
  #  warnings() # check for warnings
  
  mod_min <- gamlss(
    frtpor ~ year + pb(age) + pcat(qimd) +
      sex  + pcat(ethnicity_grp) + pcat(sha) + urban_rural  ,
    family = "ZAIG",
    weights = ds$wt_int,
    data = ds,
    #method = mixed(5, 100),
    control = con1 # for speed but looses accuracy if not 1e-3
  )
  qsave(mod_min, "./inputs/exposure_distributions/AH_test/frtpor_mod_min_zaig.qs", preset = "high")
   warnings() # check for warnings

   # Stepwise model selection using a Generalized Akaike Information Criterion
  # mod_min <- qread( "./inputs/exposure_distributions/AH_test/frtpor_mod_min_zinbi.qs")
   frtpor_modelA <- stepGAICAll.A(
     mod_min,
     scope = list(
       lower = ~ year +  + pb(age) + pcat(qimd) + 
         sex +  pcat(ethnicity_grp) + pcat(sha) + urban_rural,
       upper = ~ (year + pb(age) + pcat(qimd) + 
                    sex +  pcat(ethnicity_grp) + pcat(sha) + urban_rural  )^2
     ),
     sigma.scope = list(
       lower = ~1,
       upper = ~ year + pb(age) + pcat(qimd) + 
         sex + pcat(ethnicity_grp) + pcat(sha) + urban_rural
     ),
     nu.scope = list(
       lower = ~1,
       upper = ~ year + pb(age) + pcat(qimd) + 
         sex +  pcat(ethnicity_grp) + pcat(sha) + urban_rural 
     ),
     parallel = "multicore",
     ncpus = 12L,
     weights = wt_int, 
     trace = TRUE  # This is so you can see all the models it has tried to fit 
   )
   warnings()
   
   frtpor_modelA
   
   frtpor_modelA <- update(frtpor_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy
   
   frtpor_modelA$call 
   
   qsave(frtpor_modelA, "./inputs/exposure_distributions/AH_test/frtpor_model_ZAIG.qs", preset = "high")

   
   GAIC.table(frtpor_modelA, 
              #frtpor_modelB,
              mod_min)
   
   # # Double check that the distribution is still a good one
   # tt <- chooseDist(frtpor_modelA,
   #                  type = "realplus",
   #                  trace = TRUE, data = ds,
   #                  parallel = "multicore", ncpus = 15L
   # )
   # qsave(tt, "./inputs/exposure_distributions/AH_test/new_distr.qs", preset = "high")
   # # # tt <- qread("./inputs/exposure_distributions/AH_test/new_distr.qs")
   # # 
   # which(tt == min(tt[, 1], na.rm = TRUE), arr.ind = TRUE) # Best distr based on AIC
   # which(tt == min(tt[, 3], na.rm = TRUE), arr.ind = TRUE) # Best distr based on BIC
   # 
   #This is saying pareto type 2. 2 parameter distribution
   

}


if (update){
  frtpor_model <- qread("./inputs/exposure_distributions/AH_test/frtpor_model_ZINBI.qs")
   frtpor_model <- update(frtpor_model, method = mixed(100, 400))
   qsave(frtpor_model, "./inputs/exposure_distributions/AH_test/frtpor_model_ZINBI.qs")
   frtpor_model$hsedata <- copy(ds)
   qsave(frtpor_model, "./secure_data/lifecourse_models/frtpor_model.qs")
   
  # 
  # frtpor_modelbest <- update(frtpor_modelA, family = "ZAGA")
  # # warnings()
  # GAIC.table(frtpor_modelbest, frtpor_modelA)
  # # warnings()
  # qsave(frtpor_modelbest, "./inputs/exposure_distributions/AH_test/frtpor_model_ZAIG.qs", preset = "high")
  # #frtpor_modelbest <- qread("./frtpor_modelbest.qs")
  # 
  # #qsave(frtpor_modelbest, "./frtpor_modelbest.qs", preset = "high")
  #  
  # frtpor_modelbest <- copy(frtpor_modelA )
   # Code to create a table with predictions
  trms <- all.vars(formula(frtpor_model))[-1] # -1 excludes dependent var
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
        x[, (frtpor_model$parameters) := predictAll(frtpor_model, .SD, data = ds), .SDcols = trms]
      }
    )
  newdata <- rbindlist(newdata) # bind the chunks back
  #View(head(newdata, 1e3))

  write_fst(newdata, path = "./inputs/exposure_distributions/AH_test/frtpor_table.fst", compress = 100L) # This is what the model use as an input
  # newdata <- read_fst("./inputs/exposure_distributions/AH_test/frtpor_table.fst", as.data.table = T)

  print("Table saved")

  
}

if (diagnostics) {
  #frtpor_model <- qread("./inputs/exposure_distributions/AH_test/frtpor_modelbest.qs")
  #frtpor_model <- qread("./inputs/exposure_distributions/AH_test/frtpor_model_ZAIG.qs")
  frtpor_model <- qread("./inputs/exposure_distributions/AH_test/frtpor_model_ZINBI.qs")
  
  plot(frtpor_model)
  
  wp(frtpor_model) # detrended QQ-plot
  wp(resid = resid(frtpor_model)) # equivalen to the one above
  wp(resid = resid(frtpor_model), ylim.all = 80 * sqrt(1 / length(resid(frtpor_model))))
  tt <- wp(frtpor_model, xvar = ds$age, n.inter = 6)
  wp(frtpor_model, xvar = ds$year, n.inter = 10)
  wp(resid = resid(frtpor_model), xvar = ~ ds$qimd)
  
  dtop(frtpor_model, xvar = ds$age)
  
  rqres.plot(frtpor_model)
  rqres.plot(frtpor_model, type = "QQ")
}

if (plots) {
  #frtpor_model <- qread("./frtpor_modelbest.qs")
  #frtpor_model <- qread("./inputs/exposure_distributions/AH_test/frtpor_modelbest.qs")
  dir.create("./inputs/exposure_distributions/AH_test/validation_synthpop_models", FALSE)
  
  #source("misc_functions.R")
  library(dqrng)
  library(ggplot2)
  library(cowplot)
  theme_set(theme_cowplot())
  
  
  xlab_nam <- expression(bold(frtpor ~ (mmHg)))
  frtpor_model_tbl <- read_fst("./inputs/exposure_distributions/AH_test/frtpor_table.fst", as.data.table = TRUE)
  zz <-
    validate_gamlss_tbl(ds[age > 19], frtpor_model_tbl, 50, "frtpor", paste0("q", frtpor_model$family[1]))[
      between(frtpor, quantile(frtpor, 0.01), quantile(frtpor, 0.99))
    ]
  zz[, weight := wt_int / sum(wt_int), by = "type"]
  
  png(
    "./inputs/exposure_distributions/AH_test/validation_synthpop_models/frtpor_rel_dist.png",
    3840,
    2160,
    pointsize = 48
  )
  reldist_diagnostics(zz[type == "Observed", frtpor],
                      zz[type == "Modelled", frtpor],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = xlab_nam,
                      100
  )
  dev.off()
  
  plot_synthpop_val(zz, frtpor, "agegrp10", "wt_int", "Fruit portions by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, "year", "wt_int", "Fruit portions  by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, "qimd", "wt_int", "Fruit portions  by QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, "sha", "wt_int", "Fruit portions  by SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, "ethnicity_grp", "wt_int", "Fruit portions  by ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, "urban_rural", "wt_int", "Fruit portions  by rurality", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, c("year", "agegrp10"), "wt_int", "Fruit portions  by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, c("year", "qimd"), "wt_int", "Fruit portions  by year and QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, c("year", "sha"), "wt_int", "Fruit portions  by year and SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, c("year", "ethnicity_grp"), "wt_int", "Fruit portions  by year and ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, frtpor, c("year", "urban_rural"), "wt_int", "Fruit portions  by year and rurality", xlab_nam, FALSE, FALSE)
}
