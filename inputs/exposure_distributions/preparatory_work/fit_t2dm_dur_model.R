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
seed                <- 16L


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
if (!require(IMPACTncdEnglmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/")
  library(IMPACTncdEnglmisc)
}

dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90),
                     .(dm_dur, year, age, agegrp10, sex, qimd, bmi,
                       ethnicity, sha, wt_int)])
dt[, age := scale(age, 62.5, 13.3)]
set.seed(seed)

if (file.exists("./preparatory_work/marginal_distr_dm_dur.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_dm_dur.qs")
} else {
  marg_distr <- distr_best_fit(dt, "dm_dur", "wt_int", "counts")
  qsave(marg_distr, "./preparatory_work/marginal_distr_dm_dur.qs")
}
head(marg_distr$fits)
distr_nam <- names(marg_distr$fits[1]) # GPO
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.
distr_validation(marg_distr, dt[between(dm_dur, 0, 60),
                                .(var = dm_dur, wt = wt_int)],
                 expression(bold(dm_dur)), TRUE)


if (univariate_analysis) {
  age_scaled <- scale(18:90, 62.5, 13.3)
  dt[, .(dm_dur_mean = wtd.mean(dm_dur, weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, dm_dur_mean)]

  m_age1 <- gamlss(
    dm_dur ~ pb(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age1, dt, data.frame("age" = age_scaled)), col = "blue1")

  m_age2 <- gamlss(
    dm_dur ~ age,
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age2, dt, data.frame("age" = age_scaled)), col = "red1")

  m_age3 <- gamlss(
    dm_dur ~ cs(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age3, dt, data.frame("age" = age_scaled)), col = "green1")
  setDT(GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_age2
  setDT(GAIC(m_age1, m_age2, m_age3, k = 2), TRUE)[, setorder(.SD, AIC)] # AIC m_age1

  dt[, .(dm_dur_mean = wtd.mean(dm_dur, weight = wt_int)), keyby = .(year)
     ][, plot(year, dm_dur_mean, xlim = c(3, 40), ylim = c(0, 60))]

  m_year1 <- gamlss(
    dm_dur ~ year,
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(50, 20)
  )
  lines(3:40, mean_predictAll(m_year1, dt, data.frame("year" = 3:40)), col = "blue1")


  m_year2 <- gamlss(
    dm_dur ~ log(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year2, dt, data.frame("year" = 3:40)), col = "red1")

  m_year3 <- gamlss(
    dm_dur ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year3, dt, data.frame("year" = 3:40)), col = "green1")

  setDT(GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_year32
  setDT(GAIC(m_year1, m_year2, m_year3, k = 2), TRUE)[, setorder(.SD, -AIC)] # AIC m_year3 = m_year2
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

}


dm_dur_model <- gamlss(
  dm_dur ~ age * pcat(qimd) * sex,
  sigma.fo = ~ age,
  family = distr_nam,
  weights = dt$wt_int,
  data = dt,
  method = mixed(400, 400),
  control = con1
)

dm_dur_model <- stepGAIC(dm_dur_model, direction = "backward")

dm_dur_model$data <- copy(dt)

qsave(dm_dur_model, "./secure_data/lifecourse_models/dm_dur_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(dm_dur_model))[-1] # -1 excludes dependent var
newdata <- CJ(age_int = 20:90, sex = unique(dt$sex))
newdata[, age := scale(age_int, 62.5, 13.3)]

newdata <- split(newdata, by = "sex")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma") := predictAll(dm_dur_model, .SD, data = dm_dur_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]
write_fst(newdata, "./inputs/exposure_distributions/dm_dur_table.fst", 100L)

print("Table saved")


# Diagnostics -------------------------------------------------------------
if (diagnostics) {
  dm_dur_model <- qread("./secure_data/lifecourse_models/dm_dur_model.qs")
  dm_dur_model$data <- NULL
  wp(dm_dur_model)
  wp(dm_dur_model, xvar = ~dt$age)

  plot(dm_dur_model)
}

# Plots -------------------------------------------------------------
if (plots) {
  xlab_nam <- expression(bold(T2DM~(dgn)))
  dt[, age := age * 13.3 + 62.5] # descale
  dt[, bmi := round(clamp(dt$bmi, 18, 50), -1)]
  dm_dur_model_tbl <- read_fst("./inputs/exposure_distributions/dm_dur_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt[age >= 20L], dm_dur_model_tbl, 10, "dm_dur", paste0("q", distr_nam))

  # [between(dm_dur, quantile(dm_dur, 0.01), quantile(dm_dur, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_int / sum(wt_int), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/T2dm_dur_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", dm_dur],
                        zz[type == "Modelled", dm_dur],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        2/3, TRUE)
    dev.off()
  })
  future(plot_synthpop_val(zz, dm_dur, "agegrp10", "wt_int",
                           "T2DM (dgn) by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dur, "year", "wt_int",
                           "T2DM (dgn) by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dur, "qimd", "wt_int",
                           "T2DM (dgn) by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dur, "sha", "wt_int",
                           "T2DM (dgn) by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dur, "ethnicity", "wt_int",
                           "T2DM (dgn) by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dur, c("year", "agegrp10"), "wt_int",
                           "T2DM (dgn) by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dur, c("year", "qimd"), "wt_int",
                           "T2DM (dgn) by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dur, c("year", "sha"), "wt_int",
                           "T2DM (dgn) by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dur, c("year", "ethnicity"), "wt_int",
                           "T2DM (dgn) by year and ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dur, "bmi", "wt_int",
                           "T2DM (dgn) by bmi", xlab_nam, FALSE, FALSE))

}

