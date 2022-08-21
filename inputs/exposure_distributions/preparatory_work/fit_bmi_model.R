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
seed                <- 42L



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
sourceCpp("./preparatory_work/BCPEo_distribution.cpp", cacheDir = "./.CppCache/")

dt <- na.omit(HSE_ts[wt_nurse > 0 & between(age, 20, 90),
                     .(bmi, year, age, agegrp10, sex, qimd, ethnicity, sha, smok_status, wt_nurse)])
dt[, age := scale(age, 50.8, 16.8)]
set.seed(seed)
# lns <- sample(nrow(dt), round(nrow(dt) * 0.8))
# dt_trn   <- dt[lns] # train dataset
# dt_crv   <- dt[!lns]  # cross-validation dataset

if (file.exists("./preparatory_work/marginal_distr_bmi.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_bmi.qs")
} else {
  marg_distr <- distr_best_fit(dt, "bmi", "wt_nurse", "realplus")
  qsave(marg_distr, "./preparatory_work/marginal_distr_bmi.qs")
}
head(marg_distr$fits)
# distr_validation(marg_distr, dt[between(bmi, 16, 50), .(var = bmi, wt = wt_nurse)],
#                  expression(bold(BMI ~ (kg/m^{2}))))

distr_nam <- names(marg_distr$fits[2]) # With the first option (GB2) final model doesn't converge
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.

if (univariate_analysis) {
  age_scaled <- scale(20:90, 50.8, 16.8)
  dt[, .(bmi_median = wtd.quantile(bmi, weight = wt_nurse)), keyby = .(age)
     ][, scatter.smooth(age, bmi_median)]

  m_age1 <- gamlss(
    bmi ~ pb(age),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "blue1")

  m_age2 <- gamlss(
    bmi ~ poly(age, 2),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "red1")

  m_age3 <- gamlss(
    bmi ~ cs(age, 2),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "green1")
  GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))) # BIC m_age2
  GAIC(m_age1, m_age2, m_age3, k = 2) # AIC m_age1
  centiles(m_age1, xvar = dt$age)
  centiles(m_age2, xvar = dt$age)

  dt[, .(bmi_median = wtd.quantile(bmi, weight = wt_nurse)), keyby = .(year)
     ][, scatter.smooth(year, bmi_median, xlim = c(3, 40), ylim = c(26, 28))]

  m_year1 <- gamlss(
    bmi ~ pb(year),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "blue1")

  m_year2 <- gamlss(
    bmi ~ log(year),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "red1")

  m_year3 <- gamlss(
    bmi ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_nurse,
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
#   bmi ~ 1,
#   family = distr_nam,
#   weights = dt_trn$wt_nurse,
#   data = dt_trn,
#   method = mixed(5, 100)
# )
#
# bmi_model <- stepTGDAll.A(
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

bmi_model <- gamlss(
  bmi ~ log(year - 2) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha) + pcat(smok_status)),
  ~log(year - 2) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha) + pcat(smok_status),
  ~log(year - 2) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha) + pcat(smok_status),
  ~log(year - 2) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha) + pcat(smok_status),
  family = distr_nam,
  weights = dt$wt_nurse,
  data = dt,
  method = mixed(50, 50),
  control = con1
)

# bmi_model <- update(bmi_model, control = con1)
# tt <- chooseDist(bmi_model, type = "realplus", trace = TRUE, data = dt)
bmi_model$data <- copy(dt)

qsave(bmi_model, "./secure_data/lifecourse_models/bmi_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(bmi_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd), ethnicity = unique(dt$ethnicity),
              sha = unique(dt$sha), smok_status = unique(dt$smok_status))
newdata[, age := scale(age_int, 50.8, 16.8)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(bmi_model, .SD, data = bmi_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/bmi_table.fst", 100L)

print("Table saved")

# tt <- copy(dt)
# tt[, .N, by = .(agegrp10, sex, qimd)][, table(N)]
# tt[, bmi_rank := pctl_rank(bmi, "random"), by = .(agegrp10, sex, qimd)]
# tt[, age := as.integer(round(age * 16.8 + 50.8))]
# tt[newdata, on=.NATURAL, nomatch = NULL, `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu, tau = i.tau)]
# tt[bmi_rank == 0, bmi_rank := tt[bmi_rank > 0, min(bmi_rank)]]
# tt[bmi_rank == 1, bmi_rank := tt[bmi_rank < 1, max(bmi_rank)]]
#
# tt[, bmi := qBCPEo(bmi_rank, mu, sigma, nu, tau)]
#
#
# tt[, type := "Modelled"]
# dt[, type := "Observed"]
# zz <- rbind(tt, dt, fill = TRUE)
# zz[bmi > 50, bmi := 50]
# zz[bmi < 16, bmi := 16]
# zz <- zz[between(bmi, 17, 49)]
# zz[, weight := wt_nurse / sum(wt_nurse), by = .(type)]
# reldist_diagnostics(zz[type == "Observed", bmi],
#                     zz[type == "Modelled", bmi],
#                     zz[type == "Observed", weight],
#                     zz[type == "Modelled", weight],
#                     main = expression(bold(BMI ~ (kg / m ^ {
#                       2
#                     }))),
#                     100)

if (diagnostics) {
  bmi_model <- qread("./secure_data/lifecourse_models/bmi_model.qs")

  wp(bmi_model)
  wp(bmi_model, xvar = age)

  plot(bmi_model)
}

if (plots) {
  xlab_nam <- expression(bold(BMI ~ (kg / m^{2})))
  bmi_model <- qread("./secure_data/lifecourse_models/bmi_model.qs")

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())
  zz <-
    validate_gamlss(dt, bmi_model, 50, bmi_model$data)[between(bmi, quantile(bmi, 0.01), quantile(bmi, 0.99))]
  # [between(bmi, 16, 50)]


  future({
    zz[, weight := wt_nurse / sum(wt_nurse), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/BMI_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", bmi],
                        zz[type == "Modelled", bmi],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        100)
    dev.off()
  })
  future(plot_synthpop_val(zz, bmi, "agegrp10", "wt_nurse", "BMI by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bmi, "year", "wt_nurse", "BMI by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bmi, "qimd", "wt_nurse", "BMI by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bmi, "sha", "wt_nurse", "BMI by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bmi, "ethnicity", "wt_nurse", "BMI by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bmi, "smok_status", "wt_nurse", "BMI by smoking status", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bmi, c("year", "agegrp10"), "wt_nurse", "BMI by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bmi, c("year", "qimd"), "wt_nurse", "BMI by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bmi, c("year", "sha"), "wt_nurse", "BMI by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, bmi, c("year", "ethnicity"), "wt_nurse", "BMI by year and ethnicity", xlab_nam, FALSE, FALSE))
}

