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


# Final distribution chosen: ZINBI 

library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk
library(IMPACTncdEngl)
library(CKutils)
library(future)
library(future.apply)

# Set some variables we will use later in automation
distributions <- FALSE
univariable_analysis <- FALSE
diagnostics <- TRUE
plots <- TRUE
seed                <- 43L


# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)


ds <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90),
                     .(totalwu, year, age, agegrp10, sex, qimd, ethnicity_grp, sha, urban_rural,smok_status, wt_int)])
ds[, totalwu := round(totalwu)] # convert to counts to enable distr that are 0 inflated
# Treat occasional smokers as never smokers 
ds[smok_status == "2", smok_status := "1"]
set.seed(seed)

if(distributions){
marg_distr <- distr_best_fit(ds, 
                               "totalwu", 
                               "wt_int", 
                               "counts")
head(marg_distr$fits)
# LG     ZIPF    ZABNB    ZIBNB    ZINBI    ZANBI 
# 359449.4 387306.9 430439.0 430439.0 431174.9 431174.9 
#going for ZINBI bc it's relatively fast

# distr_validation(marg_distr, ds[between(totalwu, 0, 100), .(var = totalwu, wt = wt_int)],
#                  expression(bold(Alcohol ~ (g/d))), TRUE)
 m_zinbi <- gamlssML(ds$totalwu, ZINBI, ds$wt_int)
distr_validation(m_zinbi, ds[between(totalwu, 0, 100), .(var = totalwu, wt = wt_int)],
                 expression(bold(Alcohol ~ (g/d)~ZINBI)), TRUE)

distr_nam <- names(marg_distr$fits[5]) # No 5ZINBI is fast and a good compromise
}
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis.
distr_nam <- "ZINBI"
# Univariable analysis ----------------------------------------
if (univariable_analysis) {
  ds[, .(alcohol_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, alcohol_median)]

  m_age1 <- gamlss(
    totalwu ~ pb(age),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "blue1")

  m_age2 <- gamlss(
    totalwu ~ poly(age, 3),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "red1")

  m_age3 <- gamlss(
    totalwu ~ cs(age, 2),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = 20:90, cent = 50, data = ds), col = "green1")
  GAIC(m_age1, m_age2, m_age3, k = log(nrow(ds))) # BIC m_age1
  GAIC(m_age1, m_age2, m_age3, k = 2) # AIC m_age1
  centiles(m_age1, xvar = ds$age)
  centiles(m_age2, xvar = ds$age)

  ds[, .(alcohol_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, alcohol_median, xlim = c(3, 40), ylim = c(0, 10))]

  m_year1 <- gamlss(
    totalwu ~ pb(year),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "blue1")

  m_year2 <- gamlss(
    totalwu ~ log(year),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "red1")

  m_year3 <- gamlss(
    totalwu ~ log(year - 2),
    family = distr_nam,
    weights = ds$wt_int,
    data = ds,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = ds), col = "green1")

  GAIC(m_year1, m_year2, m_year3, k = log(nrow(ds))) # BIC m_year2
  GAIC(m_year1, m_year2, m_year3, k = 2) # AIC m_year2
  centiles(m_year1, xvar = ds$age)
  centiles(m_year2, xvar = ds$age)

}


run <- FALSE
if(run){
# Setting a minimum model - we know we want these things
mod_min <- gamlss(
  totalwu ~ log(year) + pb(age) + pcat(qimd) +
    sex + smok_status + pcat(ethnicity_grp) + pcat(sha) + urban_rural,
  family = distr_nam,
  weights = ds$wt_int,
  data = ds,
  method = mixed(5, 100),
  control = con1 # for speed but looses accuracy if not 1e-3
)
qsave(mod_min, "./inputs/exposure_distributions/AH_test/alcohol_mod_min.qs", preset = "high")
warnings() # check for warnings
mod_min <- qread( "./inputs/exposure_distributions/AH_test/alcohol_mod_min.qs")

# Stepwise model selection using a Generalized Akaike Information Criterion
alcohol_modelA <- stepGAICAll.A(
  mod_min,
  scope = list(
    lower =  ~ log(year) + pb(age) + pcat(qimd) +
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
  ) ,
  parallel = "multicore",
  ncpus = 12L,
  weights = wt_int,
  trace = TRUE  # This is so you can see all the models it has tried to fit
)
warnings()

alcohol_modelA

alcohol_modelA <- update(alcohol_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy

alcohol_modelA$call


qsave(alcohol_modelA, "./inputs/exposure_distributions/AH_test/alcohol_model.qs", preset = "high")


GAIC.table(alcohol_modelA,
           #alcohol_modelB,
           mod_min)

# Double check that the distribution is still a good one
tt <- chooseDist(alcohol_modelA,
                 type = "count",
                 trace = TRUE, data = ds,
                 parallel = "multicore", ncpus = 15L
)
qsave(tt, "./inputs/exposure_distributions/AH_test/new_distr_alc.qs", preset = "high")
#tt <- qread("./inputs/exposure_distributions/AH_test/new_distr.qs")

which(tt == min(tt[, 1], na.rm = TRUE), arr.ind = TRUE) # Best distr based on AIC
which(tt == min(tt[, 3], na.rm = TRUE), arr.ind = TRUE) # Best distr based on BIC

#ZISCHEL is better, but v slow so sticking with ZINBI 
# Did try ZINBF but was terrible

# alcohol_model <- qread( "./inputs/exposure_distributions/AH_test/alcohol_model.qs")
# alcohol_model <- update(alcohol_model, method = mixed(50, 400))
# warnings()
# qsave(alcohol_model_best,  "./inputs/exposure_distributions/AH_test/alcohol_model.qs", preset = "high")
# GAIC.table(alcohol_model, alcohol_model_best)
alcohol_model <- qread( "./inputs/exposure_distributions/AH_test/alcohol_model.qs")
alcohol_model$hsedata <- copy(ds)

qsave(alcohol_model,"./secure_data/lifecourse_models/alcohol_model.qs" )

# 

trms <- all.vars(formula(alcohol_model))[-1] # -1 excludes dependent var
newdata <- CJ(  year = 3:50,
                age = 20:90,
                sex = unique(ds$sex),
                qimd = unique(ds$qimd),
                sha = unique(ds$sha),
                smok_status = unique(ds$smok_status),
                ethnicity_grp = unique(ds$ethnicity_grp),
                urban_rural = unique(ds$urban_rural)
)

newdata <- split(newdata, by = "year")
distr_nam <- "ZINBI"
newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu") := predictAll(alcohol_model, .SD, data = alcohol_model$data), .SDcols = trms], future.seed = TRUE)
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", "ZINBI")


kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/AH_test/alcohol_table.fst", 100)

print("Table saved")

tt <- copy(ds)
tt[, .N, by = .(agegrp10, sex, qimd)][, table(N)]
tt[, totalwu_rank := pctl_rank(totalwu, "random"), by = .(agegrp10, sex, qimd)]
tt[newdata, on=.NATURAL, nomatch = NULL, `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu)]
tt[totalwu_rank == 0, totalwu_rank := tt[totalwu_rank > 0, min(totalwu_rank)]]
tt[totalwu_rank == 1, totalwu_rank := tt[totalwu_rank < 1, max(totalwu_rank)]]

tt[, totalwu := qBCPEo(totalwu_rank, mu, sigma, nu)]


tt[, type := "Modelled"]
ds[, type := "Observed"]
zz <- rbind(tt, ds, fill = TRUE)
zz[totalwu > 50, totalwu := 50]
zz[totalwu < 16, totalwu := 16]
zz <- zz[between(totalwu, 17, 49)]
zz[, weight := wt_int / sum(wt_int), by = .(type)]
reldist_diagnostics(zz[type == "Observed", totalwu],
                    zz[type == "Modelled", totalwu],
                    zz[type == "Observed", weight],
                    zz[type == "Modelled", weight],
                    main = expression(bold(Alcohol ~ (g / d))),
                    100)
# 
}
if (diagnostics) {
  alcohol_model <- qread("./inputs/exposure_distributions/AH_test/alcohol_model.qs")


  wp(alcohol_model)
  wp(alcohol_model, xvar = age)

  plot(alcohol_model)
}

if (plots) {
  xlab_nam <- expression(bold(Alcohol ~ (g / d)))
  alcohol_model_tbl <- read_fst("./inputs/exposure_distributions/AH_test/alcohol_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(ds, alcohol_model_tbl, 50, "totalwu", paste0("q", distr_nam)
                        )[between(totalwu, quantile(totalwu, 0.01), quantile(totalwu, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_int / sum(wt_int), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/Alcohol_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", totalwu],
                        zz[type == "Modelled", totalwu],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        100)
    dev.off()
  })
  plot_synthpop_val(zz, totalwu, "agegrp10", "wt_int", "Alcohol by agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, "year", "wt_int", "Alcohol by year", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, "qimd", "wt_int", "Alcohol by QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, "sha", "wt_int", "Alcohol by SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, "ethnicity_grp", "wt_int", "Alcohol by ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, "urban_rural", "wt_int", "Alcohol by rurality", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, "smok_status", "wt_int", "Alcohol by smoking status", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, c("year", "agegrp10"), "wt_int", "Alcohol by year and agegroup", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, c("year", "qimd"), "wt_int", "Alcohol by year and QIMD", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, c("year", "sha"), "wt_int", "Alcohol by year and SHA", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, c("year", "ethnicity_grp"), "wt_int", "Alcohol by year and ethnicity", xlab_nam, FALSE, FALSE)
  plot_synthpop_val(zz, totalwu, c("year", "urban_rural"), "wt_int", "Alcohol by year and rurality", xlab_nam, FALSE, FALSE)
}

