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
seed                <- 39L


if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("Rcpp", "dqrng", "fst", "qs", "gamlss", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multicore)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}
library(IMPACTncdEnglmisc)

dt <- HSE_ts[wt_int > 0 & between(age, 16, 90),
                     .(smok_status, smok_init_age, year, age, agegrp10, sex, qimd,
                       ethnicity, sha, wt_int)]
# dt[smok_init_age >= (age - 1), .N]
dt[smok_status == "1", event := 0L] # only never smokers in the denominator
dt[smok_init_age >= (age - 1) & smok_status == "4", event := 1L] # only smokers for up to a year in the denominator
dt[, smok_init_age := NULL]
dt <- na.omit(dt)
# dt[, age := scale(age, 51.7, 17.4)]
set.seed(seed)

if (file.exists("./preparatory_work/marginal_distr_smok_incid.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_smok_incid.qs")
} else {
  marg_distr <- distr_best_fit(dt, "event", "wt_int", "binom")
  qsave(marg_distr, "./preparatory_work/marginal_distr_smok_incid.qs")
}
head(marg_distr$fits)
distr_nam <- "BI"
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.


if (univariate_analysis) {
  # age_scaled <- scale(18:90, 51.7, 17.4)
  dt[, .(event_mean = wtd.mean(event, weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, event_mean)]

  m_age1 <- gamlss(
    event ~ pb(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age1, dt, data.frame("age" = age_scaled)), col = "blue1")

  m_age2 <- gamlss(
    event ~ exp(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age2, dt, data.frame("age" = age_scaled)), col = "red1")

  m_age3 <- gamlss(
    event ~ cs(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age3, dt, data.frame("age" = age_scaled)), col = "green1")
  GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))) # BIC m_age1
  GAIC(m_age1, m_age2, m_age3, k = 2) # AIC m_age1

  dt[, .(event_mean = wtd.mean(event, weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, event_mean, xlim = c(3, 40), ylim = c(0, 0.01))]

  m_year1 <- gamlss(
    event ~ pb(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(50, 20)
  )
  lines(3:40, mean_predictAll(m_year1, dt, data.frame("year" = 3:40)), col = "blue1")

  m_year2 <- gamlss(
    event ~ log(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year2, dt, data.frame("year" = 3:40)), col = "red1")

  m_year3 <- gamlss(
    event ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year3, dt, data.frame("year" = 3:40)), col = "green1")

  GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))) # BIC m_year3
  GAIC(m_year1, m_year2, m_year3, k = 2) # AIC m_year3
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

}

smok_incid_model <- gamlss(
  event ~ log(year - 2L) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha)),
  family = distr_nam,
  weights = dt$wt_int,
  data = dt,
  method = mixed(400, 400),
  control = con1
)

smok_incid_model$data <- copy(dt)

qsave(smok_incid_model, "./secure_data/lifecourse_models/smok_incid_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(smok_incid_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd), ethnicity = unique(dt$ethnicity),
              sha = unique(dt$sha))
newdata[, age := scale(age_int, 51.7, 17.4)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(smok_incid_model, .SD, data = smok_incid_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]
write_fst(newdata, "./inputs/exposure_distributions/smok_incid_table.fst", 100L)

print("Table saved")


# Diagnostics -------------------------------------------------------------
if (diagnostics) {
  smok_incid_model <- qread("./secure_data/lifecourse_models/smok_incid_model.qs")
  smok_incid_model$data <- NULL
  wp(smok_incid_model)
  wp(smok_incid_model, xvar = ~dt$age)

  plot(smok_incid_model)
}

# Plots -------------------------------------------------------------
if (plots) {
  xlab_nam <- expression(bold(Smoking~Incidence))
  dt[, age := age * 17.4 + 51.7] # descale
  smok_incid_model_tbl <- read_fst("./inputs/exposure_distributions/smok_incid_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt[age >= 20L], smok_incid_model_tbl, 50, "event", paste0("q", distr_nam))

  # [between(smok_incid, quantile(smok_incid, 0.01), quantile(smok_incid, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_int / sum(wt_int), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/Smoking_incid_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", event],
                        zz[type == "Modelled", event],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        2/3, TRUE)
    dev.off()
  })
  future(plot_synthpop_val(zz, event, "agegrp10", "wt_int",
                           "Smoking incid by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, event, "sex", "wt_int",
    "Smoking incid by sex", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, event, "year", "wt_int",
                           "Smoking incid by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, event, "qimd", "wt_int",
                           "Smoking incid by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, event, "sha", "wt_int",
                           "Smoking incid by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, event, "ethnicity", "wt_int",
                           "Smoking incid by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, event, c("year", "agegrp10"), "wt_int",
                           "Smoking incid by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, event, c("year", "qimd"), "wt_int",
                           "Smoking incid by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, event, c("year", "sha"), "wt_int",
                           "Smoking incid by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, event, c("year", "ethnicity"), "wt_int",
                           "Smoking incid by year and ethnicity", xlab_nam, FALSE, FALSE))
}

