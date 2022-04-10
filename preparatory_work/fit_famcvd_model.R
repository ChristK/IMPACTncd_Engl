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
# sourceCpp("./preparatory_work/MN_distribution.cpp", cacheDir = "./.CppCache/")


dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90),
                     .(famcvd, year, age, agegrp10, sex, qimd,
                       ethnicity, sha, wt_int)])
dt[, famcvd := as.integer(as.character(famcvd))]
dt[, age := scale(age, 47.3, 15.6)]
set.seed(seed)

if (file.exists("./preparatory_work/marginal_distr_famcvd.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_famcvd.qs")
} else {
  marg_distr <- distr_best_fit(dt, "famcvd", "wt_int", "binom")
  qsave(marg_distr, "./preparatory_work/marginal_distr_famcvd.qs")
}
head(marg_distr$fits)
distr_nam <- names(marg_distr$fits[1])
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.


if (univariate_analysis) {
  age_scaled <- scale(18:90, 47.3, 15.6)
  dt[, .(famcvd_mean = wtd.mean(famcvd, weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, famcvd_mean)]

  m_age1 <- gamlss(
    famcvd ~ pb(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age1, dt, data.frame("age" = age_scaled)), col = "blue1")

  m_age2 <- gamlss(
    famcvd ~ log(age+3),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age2, dt, data.frame("age" = age_scaled)), col = "red1")

  m_age3 <- gamlss(
    famcvd ~ cs(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age3, dt, data.frame("age" = age_scaled)), col = "green1")
  setDT(GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_age1
  setDT(GAIC(m_age1, m_age2, m_age3, k = 2), TRUE)[, setorder(.SD, AIC)] # AIC m_age1

  dt[, .(famcvd_mean = wtd.mean(famcvd, weight = wt_int)), keyby = .(year)
     ][, plot(year, famcvd_mean, xlim = c(3, 40), ylim = c(0, 1))]

  m_year1 <- gamlss(
    famcvd ~ year,
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(50, 20)
  )
  lines(3:40, mean_predictAll(m_year1, dt, data.frame("year" = 3:40)), col = "blue1")

  m_year2 <- gamlss(
    famcvd ~ log(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year2, dt, data.frame("year" = 3:40)), col = "red1")

  m_year3 <- gamlss(
    famcvd ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year3, dt, data.frame("year" = 3:40)), col = "green1")

  setDT(GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_year32  setDT(GAIC(m_year1, m_year2, m_year3, k = 2), TRUE)[, setorder(.SD, -AIC)] # AIC m_year3 = m_year2
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

  dt[, .(famcvd_mean = wtd.mean(famcvd, weight = wt_int)), keyby = .(sbp = round(sbp))
     ][, scatter.smooth(sbp, famcvd_mean)]

  m_sbp1 <- gamlss(
    famcvd ~ pb(sbp),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(70:233, mean_predictAll(m_sbp1, dt, data.frame("sbp" = 70:233)), col = "blue1")


}


famcvd_model <- gamlss(
  famcvd ~ pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha)),
  family = distr_nam,
  weights = dt$wt_int,
  data = dt,
  method = mixed(400, 400),
  control = con1
)

famcvd_model <- stepGAIC(famcvd_model, direction = "backward")

famcvd_model$data <- copy(dt)

qsave(famcvd_model, "./secure_data/lifecourse_models/famcvd_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(famcvd_model))[-1] # -1 excludes dependent var
newdata <- CJ( age_int = 20:90, qimd = unique(dt$qimd),
              ethnicity = unique(dt$ethnicity), sha = unique(dt$sha))
newdata[, age := scale(age_int, 47.3, 15.6)]

newdata <- split(newdata, by = "sha")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(famcvd_model, .SD, data = famcvd_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]
write_fst(newdata, "./inputs/exposure_distributions/famcvd_table.fst", 100L)

print("Table saved")


# Diagnostics -------------------------------------------------------------
if (diagnostics) {
  famcvd_model <- qread("./secure_data/lifecourse_models/famcvd_model.qs")
  famcvd_model$data <- NULL
  wp(famcvd_model)
  wp(famcvd_model, xvar = ~dt$age)

  plot(famcvd_model)
}

# Plots -------------------------------------------------------------
if (plots) {
  xlab_nam <- expression(bold(family~cvd))
  dt[, age := age * 15.6 + 47.3] # descale
  famcvd_model_tbl <- read_fst("./inputs/exposure_distributions/famcvd_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt[age >= 20L], famcvd_model_tbl, 10, "famcvd", paste0("q", distr_nam))

  # [between(famcvd, quantile(famcvd, 0.01), quantile(famcvd, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_int / sum(wt_int), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/Famcvd_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", famcvd],
                        zz[type == "Modelled", famcvd],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        2/3, TRUE)
    dev.off()
  })
  future(plot_synthpop_val(zz, famcvd, "agegrp10", "wt_int",
                           "Family CVD by agegroup", xlab_nam, FALSE, FALSE))
  # future(plot_synthpop_val(zz, famcvd, "year", "wt_int",
  #                          "Family CVD (diagnosed) by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, famcvd, "qimd", "wt_int",
                           "Family CVD (diagnosed) by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, famcvd, "sha", "wt_int",
                           "Family CVD (diagnosed) by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, famcvd, "ethnicity", "wt_int",
                           "Family CVD (diagnosed) by ethnicity", xlab_nam, FALSE, FALSE))
  # future(plot_synthpop_val(zz, famcvd, c("year", "agegrp10"), "wt_int",
  #                          "Family CVD (diagnosed) by year and agegroup", xlab_nam, FALSE, FALSE))
  # future(plot_synthpop_val(zz, famcvd, c("year", "qimd"), "wt_int",
  #                          "Family CVD (diagnosed) by year and QIMD", xlab_nam, FALSE, FALSE))
  # future(plot_synthpop_val(zz, famcvd, c("year", "sha"), "wt_int",
  #                          "Family CVD (diagnosed) by year and SHA", xlab_nam, FALSE, FALSE))
  # future(plot_synthpop_val(zz, famcvd, c("year", "ethnicity"), "wt_int",
  #                          "Family CVD (diagnosed) by year and ethnicity", xlab_nam, FALSE, FALSE))
  # future(plot_synthpop_val(zz, famcvd, "sbp", "wt_nurse",
  #                          "bp med by sbp", xlab_nam, FALSE, FALSE))

}

