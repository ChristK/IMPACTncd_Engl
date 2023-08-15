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


setwd("/home/ckypridMy_ModelsIMPACTncd_Engl/")
univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- TRUE
seed                <- 28L


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


dt <- na.omit(HSE_ts[wt_nurse > 0 & between(age, 20, 90) & year >= 12, # methodology changed in 2012
                     .(bp_med, sbp, year, age, agegrp10, sex, qimd,
                       ethnicity, sha, wt_nurse)])
dt[, bp_med := as.integer(as.character(bp_med))]
dt[, age := scale(age, 53.2, 17.2)]
set.seed(seed)

if (file.exists("./preparatory_work/marginal_distr_bp_med.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_bp_med.qs")
} else {
  marg_distr <- distr_best_fit(dt, "bp_med", "wt_nurse", "binom")
  qsave(marg_distr, "./preparatory_work/marginal_distr_bp_med.qs")
}
head(marg_distr$fits)
distr_nam <- names(marg_distr$fits[1])
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.


if (univariate_analysis) {
  age_scaled <- scale(18:90, 53.2, 17.2)
  dt[, .(bp_med_mean = wtd.mean(bp_med, weight = wt_nurse)), keyby = .(age)
     ][, scatter.smooth(age, bp_med_mean)]

  m_age1 <- gamlss(
    bp_med ~ pb(age),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age1, dt, data.frame("age" = age_scaled)), col = "blue1")

  m_age2 <- gamlss(
    bp_med ~ log(age+3),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age2, dt, data.frame("age" = age_scaled)), col = "red1")

  m_age3 <- gamlss(
    bp_med ~ cs(age),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age3, dt, data.frame("age" = age_scaled)), col = "green1")
  setDT(GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_age3
  setDT(GAIC(m_age1, m_age2, m_age3, k = 2), TRUE)[, setorder(.SD, AIC)] # AIC m_age1

  dt[, .(bp_med_mean = wtd.mean(bp_med, weight = wt_nurse)), keyby = .(year)
     ][, plot(year, bp_med_mean, xlim = c(3, 40), ylim = c(0, 1))]

  m_year1 <- gamlss(
    bp_med ~ year,
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(50, 20)
  )
  lines(3:40, mean_predictAll(m_year1, dt, data.frame("year" = 3:40)), col = "blue1")

  m_year2 <- gamlss(
    bp_med ~ log(year),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year2, dt, data.frame("year" = 3:40)), col = "red1")

  m_year3 <- gamlss(
    bp_med ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year3, dt, data.frame("year" = 3:40)), col = "green1")

  setDT(GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_year32  setDT(GAIC(m_year1, m_year2, m_year3, k = 2), TRUE)[, setorder(.SD, -AIC)] # AIC m_year3 = m_year2
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

  dt[, .(bp_med_mean = wtd.mean(bp_med, weight = wt_nurse)), keyby = .(sbp = round(sbp))
     ][, scatter.smooth(sbp, bp_med_mean)]

  m_sbp1 <- gamlss(
    bp_med ~ pb(sbp),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(20, 20)
  )
  lines(70:233, mean_predictAll(m_sbp1, dt, data.frame("sbp" = 70:233)), col = "blue1")


}


bp_med_model <- gamlss(
  bp_med ~ year * pb(age) * pcat(qimd) * (pb(sbp) + sex +  pcat(ethnicity) + pcat(sha)),
  family = distr_nam,
  weights = dt$wt_nurse,
  data = dt,
  method = mixed(400, 400),
  control = con1
)

# I have to do backward elimination in the hope to reduce the number of predictors.
# otherwise the table has almost a billion rows
bp_med_model <- stepGAIC(bp_med_model, direction = "backward", parallel = "multicore", ncpus = 20L)

bp_med_model$data <- copy(dt)

qsave(bp_med_model, "./secure_data/lifecourse_models/bp_med_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(bp_med_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd),
              ethnicity = unique(dt$ethnicity), sha = unique(dt$sha), sbp = seq(110, 200, 10))
newdata[, age := scale(age_int, 53.2, 17.2)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(bp_med_model, .SD, data = bp_med_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]
write_fst(newdata, "./inputs/exposure_distributions/bp_med_table.fst", 100L)

print("Table saved")


# Diagnostics -------------------------------------------------------------
if (diagnostics) {
  bp_med_model <- qread("./secure_data/lifecourse_models/bp_med_model.qs")
  bp_med_model$data <- NULL
  wp(bp_med_model)
  wp(bp_med_model, xvar = ~dt$age)

  plot(bp_med_model)
}

# Plots -------------------------------------------------------------
if (plots) {
  xlab_nam <- expression(bold(bp~med))
  dt[, age := age * 17.2 + 53.2] # descale
  bp_med_model_tbl <- read_fst("./inputs/exposure_distributions/bp_med_table.fst", as.data.table =  TRUE)
  dt[, sbp:= round(clamp(sbp, 110, 200), -1)]
  zz <-
    validate_gamlss_tbl(dt[age >= 20L], bp_med_model_tbl, 10, "bp_med", paste0("q", distr_nam))

  # [between(bp_med, quantile(bp_med, 0.01), quantile(bp_med, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_nurse / sum(wt_nurse), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/BP_med_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", bp_med],
                        zz[type == "Modelled", bp_med],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        2/3, TRUE)
    dev.off()
  })
  future(plot_synthpop_val(zz, bp_med, "agegrp10", "wt_nurse",
                           "BP med by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bp_med, "year", "wt_nurse",
                           "BP med by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bp_med, "qimd", "wt_nurse",
                           "BP med by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bp_med, "sha", "wt_nurse",
                           "BP med by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bp_med, "ethnicity", "wt_nurse",
                           "BP med by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bp_med, c("year", "agegrp10"), "wt_nurse",
                           "BP med by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bp_med, c("year", "qimd"), "wt_nurse",
                           "BP med by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bp_med, c("year", "sha"), "wt_nurse",
                           "BP med by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bp_med, c("year", "ethnicity"), "wt_nurse",
                           "BP med by year and ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bp_med, "sbp", "wt_nurse",
                           "BP med by sbp", xlab_nam, FALSE, FALSE))

}

