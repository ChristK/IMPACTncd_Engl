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


setwd("/home/ckyprid/My Models/IMPACTncd_Engl/")
univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- TRUE
seed                <- 37L


if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("Rcpp", "dqrng", "fst", "qs", "gamlss", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multiprocess)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}
source("./preparatory_work/aux_functions.R", local = TRUE)
sourceCpp("./preparatory_work/DPO_distribution.cpp", cacheDir = "./.CppCache/")

dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90) & smok_status %in% 2:3, .(
  smok_quit_yrs, smok_status, year, age, agegrp10, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, smok_quit_yrs := as.integer(round(smok_quit_yrs))]
dt[, age := scale(age, 56.2, 17.0)]

set.seed(seed)

if (file.exists("./preparatory_work/marginal_distr_smok_quit_yrs.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_smok_quit_yrs.qs")
} else {
  marg_distr <- distr_best_fit(dt, "smok_quit_yrs", "wt_int", "counts")
  qsave(marg_distr, "./preparatory_work/marginal_distr_smok_quit_yrs.qs")
}
head(marg_distr$fits)
distr_nam <- names(marg_distr$fits[1])
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.

# Univariate analysis ----------------------------------------
if (univariate_analysis) {
  age_scaled <- scale(20:90, 56.2, 17.0)
  dt[, .(smok_quit_yrs_median = wtd.quantile(smok_quit_yrs, weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, smok_quit_yrs_median)]

  m_age1 <- gamlss(
    smok_quit_yrs ~ pb(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "blue1")

  m_age2 <- gamlss(
    smok_quit_yrs ~ exp(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "red1")

  m_age3 <- gamlss(
    smok_quit_yrs ~ age,
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "green1")
  GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))) # BIC m_age1
  GAIC(m_age1, m_age2, m_age3, k = 2) # AIC m_age1. I inteniolly chose this for consistency with orher models

  dt[, .(smok_quit_yrs_median = wtd.quantile(smok_quit_yrs, weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, smok_quit_yrs_median, xlim = c(3, 40), ylim = c(0, 20))]

  m_year1 <- gamlss(
    smok_quit_yrs ~ pb(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "blue1")

  m_year2 <- gamlss(
    smok_quit_yrs ~ log(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "red1")

  m_year3 <- gamlss(
    smok_quit_yrs ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "green1")

  GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))) # BIC m_year3
  GAIC(m_year1, m_year2, m_year3, k = 2) # AIC m_year3.
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

}

smok_quit_yrs_model <- gamlss(
  smok_quit_yrs ~ log(year - 2) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha)),
  ~log(year - 2) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha),
  family = distr_nam,
  weights = dt$wt_int,
  data = dt,
  method = mixed(400, 400),
  control = con1
)
# I intentionally exclude smok_status from above because dt has only levels 2 and 3 of the factor.
# If I did prediction could be messy

smok_quit_yrs_model$data <- copy(dt)

qsave(smok_quit_yrs_model, "./secure_data/lifecourse_models/smok_quit_yrs_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(smok_quit_yrs_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd),
              ethnicity = unique(dt$ethnicity),
              sha = unique(dt$sha))
newdata[, age := scale(age_int, 56.2, 17.0)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma") := predictAll(smok_quit_yrs_model, .SD, data = smok_quit_yrs_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]
write_fst(newdata, "./inputs/exposure_distributions/smok_quit_yrs_table.fst", 100L)

print("Table saved")



if (diagnostics) {
  smok_quit_yrs_model <- qread("./secure_data/lifecourse_models/smok_quit_yrs_model.qs")

  wp(smok_quit_yrs_model)
  wp(smok_quit_yrs_model, xvar = age)

  plot(smok_quit_yrs_model)
}

if (plots) {
  xlab_nam <- expression(bold(Years~since~cessation))
  dt[, age := age * 17 + 56.2] # descale
  smok_quit_yrs_model_tbl <- read_fst("./inputs/exposure_distributions/smok_quit_yrs_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt, smok_quit_yrs_model_tbl, 50, "smok_quit_yrs", paste0("my_q", distr_nam)
    )[between(smok_quit_yrs, quantile(smok_quit_yrs, 0.01), quantile(smok_quit_yrs, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_int / sum(wt_int), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/Smok_quit_yrs_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", smok_quit_yrs],
                        zz[type == "Modelled", smok_quit_yrs],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        100)
    dev.off()
  })
  future(plot_synthpop_val(zz, smok_quit_yrs, "agegrp10", "wt_int", "Years since cessation by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_quit_yrs, "year", "wt_int", "Years since cessation by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_quit_yrs, "qimd", "wt_int", "Years since cessation by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_quit_yrs, "sha", "wt_int", "Years since cessation by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_quit_yrs, "ethnicity", "wt_int", "Years since cessation by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_quit_yrs, "smok_status", "wt_int", "Years since cessation by smoking status", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_quit_yrs, c("year", "agegrp10"), "wt_int", "Years since cessation by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_quit_yrs, c("year", "qimd"), "wt_int", "Years since cessation by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_quit_yrs, c("year", "sha"), "wt_int", "Years since cessation by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_quit_yrs, c("year", "ethnicity"), "wt_int", "Years since cessation by year and ethnicity", xlab_nam, FALSE, FALSE))
}

