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
seed                <- 41L



if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("Rcpp", "fst", "qs", "gamlss", "reldist", "future", "future.apply", "data.table"))
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
                     .(totalwu, year, age, agegrp10, sex, qimd, ethnicity, sha, smok_status, wt_int)])
dt[, totalwu := round(totalwu)] # convert to counts to enable distr that are 0 inflated
dt[, age := scale(age, 51.3, 17.5)]
set.seed(seed)
# lns <- sample(nrow(dt), round(nrow(dt) * 0.8))
# dt_trn   <- dt[lns] # train dataset
# dt_crv   <- dt[!lns]  # cross-validation dataset

if (file.exists("./preparatory_work/marginal_distr_alcohol.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_alcohol.qs")
} else {
  marg_distr <- distr_best_fit(dt, "totalwu", "wt_int", "counts")
  qsave(marg_distr, "./preparatory_work/marginal_distr_alcohol.qs")
}
head(marg_distr$fits)
# distr_validation(marg_distr, dt[between(totalwu, 0, 100), .(var = totalwu, wt = wt_int)],
#                  expression(bold(Alcohol ~ (g/d))), TRUE)
# m11 <- gamlssML(dt$totalwu, ZINBI, dt$wt_int)
# distr_validation(m11, dt[between(totalwu, 0, 100), .(var = totalwu, wt = wt_int)],
#                  expression(bold(Alcohol ~ (g/d)~ZINBI)), TRUE)

distr_nam <- names(marg_distr$fits[3]) # 1 is ZIBNB and is very slow. No 3, ZINBI is fast and a good compromise
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.

# Univariate analysis ----------------------------------------
if (univariate_analysis) {
  age_scaled <- scale(20:90, 51.3, 17.5)
  dt[, .(alcohol_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, alcohol_median)]

  m_age1 <- gamlss(
    totalwu ~ pb(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "blue1")

  m_age2 <- gamlss(
    totalwu ~ poly(age, 3),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "red1")

  m_age3 <- gamlss(
    totalwu ~ cs(age, 2),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "green1")
  GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))) # BIC m_age1
  GAIC(m_age1, m_age2, m_age3, k = 2) # AIC m_age1
  centiles(m_age1, xvar = dt$age)
  centiles(m_age2, xvar = dt$age)

  dt[, .(alcohol_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, alcohol_median, xlim = c(3, 40), ylim = c(0, 10))]

  m_year1 <- gamlss(
    totalwu ~ pb(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "blue1")

  m_year2 <- gamlss(
    totalwu ~ log(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "red1")

  m_year3 <- gamlss(
    totalwu ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "green1")

  GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))) # BIC m_year3
  GAIC(m_year1, m_year2, m_year3, k = 2) # AIC m_year3
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

}

# mod_min <- gamlss(
#   totalwu ~ 1,
#   family = distr_nam,
#   weights = dt_trn$wt_int,
#   data = dt_trn,
#   method = mixed(5, 100)
# )
#
# alcohol_model <- stepTGDAll.A(
#   mod_min,
#   scope = list(
#     lower =  ~ 1,
#     upper =  ~ ga( ~ s(log(year), age, by = sex)) + (
#       log(year) + pb(age) + sex + pcat(qimd) +
#        pcat(ethnicity) + pcat(sha)
#     ) ^ 2
#   ),
#   sigma.scope = list(
#     lower =  ~ 1,
#     upper =  ~ ga( ~ s(log(year), age, by = sex)) + (
#       log(year) +  pb(age) + sex + pcat(qimd) + pcat(sha)) ^ 2
#   ),
#   nu.scope = list(
#     lower =  ~ 1,
#     upper =  ~ ga( ~ s(log(year), age, by = sex)) + (
#       log(year) +  pb(age) + sex + pcat(qimd)) ^ 2
#     ),
#   newdata = dt_crv,
#   parallel = "multicore",
#   ncpus = 16L
# )

con1 <- gamlss.control(c.crit = 1e-2) # increase for faster exploratory analysis.
alcohol_model <- gamlss(
  totalwu ~ log(year) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha) + pcat(smok_status)),
  ~log(year) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha) + pcat(smok_status),
  ~log(year) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha) + pcat(smok_status),
  family = distr_nam,
  weights = dt$wt_int,
  data = dt,
  method = mixed(100, 100),
  control = con1
)

# alcohol_model <- update(alcohol_model, method = mixed(200, 200))
# tt <- chooseDist(alcohol_model, type = "realplus", trace = TRUE, data = dt)
alcohol_model$data <- copy(dt)

qsave(alcohol_model, "./secure_data/lifecourse_models/alcohol_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(alcohol_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd), ethnicity = unique(dt$ethnicity),
              sha = unique(dt$sha), smok_status = unique(dt$smok_status))
newdata[, age := scale(age_int, 51.3, 17.5)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu") := predictAll(alcohol_model, .SD, data = alcohol_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/alcohol_table.fst", 100)

print("Table saved")

# tt <- copy(dt)
# tt[, .N, by = .(agegrp10, sex, qimd)][, table(N)]
# tt[, totalwu_rank := pctl_rank(totalwu, "random"), by = .(agegrp10, sex, qimd)]
# tt[, age := as.integer(round(age * 16.8 + 50.8))]
# tt[newdata, on=.NATURAL, nomatch = NULL, `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu, tau = i.tau)]
# tt[totalwu_rank == 0, totalwu_rank := tt[totalwu_rank > 0, min(totalwu_rank)]]
# tt[totalwu_rank == 1, totalwu_rank := tt[totalwu_rank < 1, max(totalwu_rank)]]
#
# tt[, totalwu := qBCPEo(totalwu_rank, mu, sigma, nu, tau)]
#
#
# tt[, type := "Modelled"]
# dt[, type := "Observed"]
# zz <- rbind(tt, dt, fill = TRUE)
# zz[totalwu > 50, totalwu := 50]
# zz[totalwu < 16, totalwu := 16]
# zz <- zz[between(totalwu, 17, 49)]
# zz[, weight := wt_int / sum(wt_int), by = .(type)]
# reldist_diagnostics(zz[type == "Observed", totalwu],
#                     zz[type == "Modelled", totalwu],
#                     zz[type == "Observed", weight],
#                     zz[type == "Modelled", weight],
#                     main = expression(bold(Alcohol ~ (g / d))),
#                     100)

if (diagnostics) {
  alcohol_model <- qread("./secure_data/lifecourse_models/alcohol_model.qs")

  wp(alcohol_model)
  wp(alcohol_model, xvar = age)

  plot(alcohol_model)
}

if (plots) {
  xlab_nam <- expression(bold(Alcohol ~ (g / d)))
  dt[, age := age * 17.5 + 51.3] # descale
  alcohol_model_tbl <- read_fst("./inputs/exposure_distributions/alcohol_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt, alcohol_model_tbl, 50, "totalwu", paste0("q", distr_nam)
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
  future(plot_synthpop_val(zz, totalwu, "agegrp10", "wt_int", "Alcohol by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, totalwu, "year", "wt_int", "Alcohol by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, totalwu, "qimd", "wt_int", "Alcohol by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, totalwu, "sha", "wt_int", "Alcohol by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, totalwu, "ethnicity", "wt_int", "Alcohol by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, totalwu, "smok_status", "wt_int", "Alcohol by smoking status", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, totalwu, c("year", "agegrp10"), "wt_int", "Alcohol by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, totalwu, c("year", "qimd"), "wt_int", "Alcohol by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, totalwu, c("year", "sha"), "wt_int", "Alcohol by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, totalwu, c("year", "ethnicity"), "wt_int", "Alcohol by year and ethnicity", xlab_nam, FALSE, FALSE))
}

