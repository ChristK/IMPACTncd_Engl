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


setwd("/home/ckyprid/My_Models/IMPACTncd_Engl/")
univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- TRUE
seed                <- 21L


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
sourceCpp("./preparatory_work/aux_functions.cpp", cacheDir = "./.CppCache/")


dt <- na.omit(HSE_ts[wt_blood > 0 & between(age, 20, 90) & year >= 9,
                     .(statin_px, tchol, year, age, agegrp10, sex, qimd,
                       ethnicity, sha, wt_blood)])
dt[, statin_px := as.integer(as.character(statin_px))]
dt[, age := scale(age, 52.4, 16.6)]
set.seed(seed)

if (file.exists("./preparatory_work/marginal_distr_statin_px.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_statin_px.qs")
} else {
  marg_distr <- distr_best_fit(dt, "statin_px", "wt_blood", "binom")
  qsave(marg_distr, "./preparatory_work/marginal_distr_statin_px.qs")
}
head(marg_distr$fits)
distr_nam <- names(marg_distr$fits[1])
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.


if (univariate_analysis) {
  age_scaled <- scale(20:90, 52.4, 16.6)
  dt[, .(statin_px_mean = wtd.mean(statin_px, weight = wt_blood)), keyby = .(age)
     ][, scatter.smooth(age, statin_px_mean)]

  m_age1 <- gamlss(
    statin_px ~ pb(age),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age1, dt, data.frame("age" = age_scaled)), col = "blue1")

  m_age2 <- gamlss(
    statin_px ~ log(age+3),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age2, dt, data.frame("age" = age_scaled)), col = "red1")

  m_age3 <- gamlss(
    statin_px ~ cs(age),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age3, dt, data.frame("age" = age_scaled)), col = "green1")
  setDT(GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_age3
  setDT(GAIC(m_age1, m_age2, m_age3, k = 2), TRUE)[, setorder(.SD, AIC)] # AIC m_age1

  dt[, .(statin_px_mean = wtd.mean(statin_px, weight = wt_blood)), keyby = .(year)
     ][, plot(year, statin_px_mean, xlim = c(3, 40), ylim = c(0, 1))]

  m_year1 <- gamlss(
    statin_px ~ year,
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(50, 20)
  )
  lines(3:40, mean_predictAll(m_year1, dt, data.frame("year" = 3:40)), col = "blue1")

  m_year2 <- gamlss(
    statin_px ~ log(year),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(50, 100)
  )
  lines(3:40, mean_predictAll(m_year2, dt, data.frame("year" = 3:40)), col = "red1")

  m_year3 <- gamlss(
    statin_px ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(50, 100)
  )
  lines(3:40, mean_predictAll(m_year3, dt, data.frame("year" = 3:40)), col = "green1")

  setDT(GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_year32
  setDT(GAIC(m_year1, m_year2, m_year3, k = 2), TRUE)[, setorder(.SD, -AIC)] # AIC m_year3 = m_year2
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

  dt[, .(statin_px_mean = wtd.mean(statin_px, weight = wt_blood)), keyby = .(tchol = round(tchol))
     ][, scatter.smooth(tchol, statin_px_mean)]

  m_tchol1 <- gamlss(
    statin_px ~ pb(tchol),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(2:12, mean_predictAll(m_tchol1, dt, data.frame("tchol" = 2:12)), col = "blue1")
}


statin_px_model <- gamlss(
  statin_px ~ log(year) * pb(age) * pcat(qimd) * (pb(tchol) + sex +  pcat(ethnicity) + pcat(sha)),
  family = distr_nam,
  weights = dt$wt_blood,
  data = dt,
  method = mixed(400, 400),
  control = con1
)

# I have to do backward elimination in the hope to reduce the number of predictors.
# otherwise the table has almost a billion rows
statin_px_model <- stepGAIC(statin_px_model, direction = "backward", parallel = "multicore", ncpus = 20L)

statin_px_model$data <- copy(dt)

qsave(statin_px_model, "./secure_data/lifecourse_models/statin_px_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(statin_px_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd),
              ethnicity = unique(dt$ethnicity), sha = unique(dt$sha), tchol = seq(2, 12, 1))
newdata[, age := scale(age_int, 52.4, 16.6)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(statin_px_model, .SD, data = statin_px_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/statin_px_table.fst", 100L)

print("Table saved")


# Diagnostics -------------------------------------------------------------
if (diagnostics) {
  statin_px_model <- qread("./secure_data/lifecourse_models/statin_px_model.qs")
  statin_px_model$data <- NULL
  wp(statin_px_model)
  wp(statin_px_model, xvar = ~dt$age)

  plot(statin_px_model)
}

# Plots -------------------------------------------------------------
if (plots) {
  xlab_nam <- expression(bold(Statins~px))
  dt[, age := age * 16.6 + 52.4] # descale
  statin_px_model_tbl <- read_fst("./inputs/exposure_distributions/statin_px_table.fst", as.data.table =  TRUE)
  dt[, tchol:= round(clamp(tchol, 2, 12), 0)]
  zz <-
    validate_gamlss_tbl(dt[age >= 20L], statin_px_model_tbl, 10, "statin_px", paste0("q", distr_nam))

  # [between(statin_px, quantile(statin_px, 0.01), quantile(statin_px, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_blood / sum(wt_blood), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/statin_px_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", statin_px],
                        zz[type == "Modelled", statin_px],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        2/3, TRUE)
    dev.off()
  })
  future(plot_synthpop_val(zz, statin_px, "agegrp10", "wt_blood",
                           "Statins px by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, statin_px, "year", "wt_blood",
                           "Statins px by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, statin_px, "qimd", "wt_blood",
                           "Statins px by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, statin_px, "sha", "wt_blood",
                           "Statins px by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, statin_px, "ethnicity", "wt_blood",
                           "Statins px by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, statin_px, c("year", "agegrp10"), "wt_blood",
                           "Statins px by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, statin_px, c("year", "qimd"), "wt_blood",
                           "Statins px by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, statin_px, c("year", "sha"), "wt_blood",
                           "Statins px by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, statin_px, c("year", "ethnicity"), "wt_blood",
                           "Statins px by year and ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, statin_px, "sbp", "wt_blood",
                           "Statins px by sbp", xlab_nam, FALSE, FALSE))

}

