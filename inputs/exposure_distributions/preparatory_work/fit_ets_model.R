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
seed                <- 30L


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
# sourceCpp("./preparatory_work/MN_distribution.cpp", cacheDir = "./.CppCache/")


dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90) & year != 7 & year != 8, #year 7 is half with and half without smokefree policy. Year 8 was influencial
                     .(ets, smok_status, year, age, agegrp10, sex, qimd,
                       ethnicity, sha, wt_int)])
dt[, ets := as.integer(as.character(ets))]
# Calculate smoking prevalence to be used as externality
tt1 <- dt[smok_status == "4",  sum(wt_int), by = .(year, sha)]
tt2 <- dt[,  sum(wt_int), by = .(year, sha)]
tt1[tt2, smok_prev := V1/i.V1, on = c("year", "sha")]
dt[tt1, smok_prev := i.smok_prev, on = c("year", "sha")]
rm(tt1, tt2)
dt[, age := scale(age, 50.8, 17.4)]
set.seed(seed)

if (file.exists("./preparatory_work/marginal_distr_ets.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_ets.qs")
} else {
  marg_distr <- distr_best_fit(dt, "ets", "wt_int", "binom")
  qsave(marg_distr, "./preparatory_work/marginal_distr_ets.qs")
}
head(marg_distr$fits)
distr_nam <- names(marg_distr$fits[1])
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.


if (univariate_analysis) {
  age_scaled <- scale(18:90, 50.8, 17.4)
  dt[, .(ets_mean = wtd.mean(ets, weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, ets_mean)]

  m_age1 <- gamlss(
    ets ~ pb(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age1, dt, data.frame("age" = age_scaled)), col = "blue1")

  m_age2 <- gamlss(
    ets ~ log(age+3),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age2, dt, data.frame("age" = age_scaled)), col = "red1")

  m_age3 <- gamlss(
    ets ~ cs(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age3, dt, data.frame("age" = age_scaled)), col = "green1")
  GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))) # BIC m_age3
  GAIC(m_age1, m_age2, m_age3, k = 2) # AIC m_age1

  dt[, .(ets_mean = wtd.mean(ets, weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, ets_mean, xlim = c(3, 40), ylim = c(0, 1))]

  m_year1 <- gamlss(
    ets ~ year * (year>7),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(50, 20)
  )
  lines(3:40, mean_predictAll(m_year1, dt, data.frame("year" = 3:40)), col = "blue1")

  m_year2 <- gamlss(
    ets ~ log(year) * (year>7),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year2, dt, data.frame("year" = 3:40)), col = "red1")

  m_year3 <- gamlss(
    ets ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year3, dt, data.frame("year" = 3:40)), col = "green1")

  GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))) # BIC m_year1
  GAIC(m_year1, m_year2, m_year3, k = 2) # AIC m_year1
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

}


ets_model <- gamlss(
  ets ~ year * (year>7) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha) + pcat(smok_status)),
  family = distr_nam,
  weights = dt$wt_int,
  data = dt,
  method = mixed(400, 400),
  control = con1
)

# smok_prev not significant

ets_model$data <- copy(dt)

qsave(ets_model, "./secure_data/lifecourse_models/ets_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(ets_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd),
              ethnicity = unique(dt$ethnicity), sha = unique(dt$sha),
              smok_status = unique(dt$smok_status))
newdata[, age := scale(age_int, 50.8, 17.4)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(ets_model, .SD, data = ets_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/ets_table.fst", 100L)

print("Table saved")


# Diagnostics -------------------------------------------------------------
if (diagnostics) {
  ets_model <- qread("./secure_data/lifecourse_models/ets_model.qs")
  ets_model$data <- NULL
  wp(ets_model)
  wp(ets_model, xvar = ~dt$age)

  plot(ets_model)
}

# Plots -------------------------------------------------------------
if (plots) {
  xlab_nam <- expression(bold(ETS))
  dt[, age := age * 17.4 + 50.8] # descale
  ets_model_tbl <- read_fst("./inputs/exposure_distributions/ets_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt[age >= 20L], ets_model_tbl, 50, "ets", paste0("q", distr_nam))

  # [between(ets, quantile(ets, 0.01), quantile(ets, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_int / sum(wt_int), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/ETS_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", ets],
                        zz[type == "Modelled", ets],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        2/3, TRUE)
    dev.off()
  })
  future(plot_synthpop_val(zz, ets, "agegrp10", "wt_int",
                           "ETS by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, "year", "wt_int",
                           "ETS by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, "qimd", "wt_int",
                           "ETS by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, "sha", "wt_int",
                           "ETS by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, "ethnicity", "wt_int",
                           "ETS by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("year", "agegrp10"), "wt_int",
                           "ETS by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("year", "qimd"), "wt_int",
                           "ETS by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("year", "sha"), "wt_int",
                           "ETS by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("year", "ethnicity"), "wt_int",
                           "ETS by year and ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ets, c("smok_status"), "wt_int",
                           "ETS by smoking status", xlab_nam, FALSE, FALSE))

}

