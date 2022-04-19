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
plots               <- FALSE
seed                <- 43L



if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("Rcpp", "dqrng", "fst", "qs", "gamlss", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multiprocess)

if (file.exists("./secure_data/HSE_ts.fst")) {
  HSE_ts <- read_fst("./secure_data/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./secure_data/preprocess_HSE.R", local = TRUE)
}
source("./secure_data/aux_functions.R", local = TRUE)
sourceCpp("./secure_data/BCPEo_distribution.cpp", cacheDir = "./.CppCache/")

dt <- na.omit(HSE_ts[wt_nurse > 0 & between(age, 20, 90),
                     .(sbp, year, age, agegrp10, sex, qimd, ethnicity, sha, smok_status, wt_nurse)])
dt[, age := scale(age, 52.1, 17.1)]
set.seed(seed)


if (file.exists("./secure_data/marginal_distr_sbp.qs")) {
  marg_distr <- qread("./secure_data/marginal_distr_sbp.qs")
} else {
  marg_distr <- distr_best_fit(dt, "sbp", "wt_nurse", "realplus")
  qsave(marg_distr, "./secure_data/marginal_distr_sbp.qs")
}
head(marg_distr$fits)
# distr_validation(marg_distr, dt[between(sbp, 16, 50), .(var = sbp, wt = wt_nurse)],
#                  expression(bold(SBP ~ (mmHg))))

distr_nam <- names(marg_distr$fits[2]) # For convenience as I already have BCPEo in C++
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.

if (univariate_analysis) {
  age_scaled <- scale(20:90, 52.1, 17.1)
  dt[, .(sbp_median = wtd.quantile(sbp, weight = wt_nurse)), keyby = .(age)
     ][, scatter.smooth(age, sbp_median)]

  m_age1 <- gamlss(
    sbp ~ pb(age),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "blue1")

  m_age2 <- gamlss(
    sbp ~ poly(age, 3),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "red1")

  m_age3 <- gamlss(
    sbp ~ cs(age),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "green1")
  GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))) # BIC m_age2
  GAIC(m_age1, m_age2, m_age3, k = 2) # AIC m_age1
  centiles(m_age1, xvar = dt$age)
  centiles(m_age2, xvar = dt$age)

  dt[, .(sbp_median = wtd.quantile(sbp, weight = wt_nurse)), keyby = .(year)
     ][, scatter.smooth(year, sbp_median, xlim = c(3, 40), ylim = c(110, 130))]

  m_year1 <- gamlss(
    sbp ~ year,
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "blue1")

  m_year2 <- gamlss(
    sbp ~ log(year),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "red1")

  m_year3 <- gamlss(
    sbp ~ log(year + 100),
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "green1")

  GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))) # BIC m_year3 = m_year1
  GAIC(m_year1, m_year2, m_year3, k = 2) # AIC m_year3 = m_year1
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

  m_sm1 <- gamlss(
    sbp ~ smok_status,
    family = distr_nam,
    weights = dt$wt_nurse,
    data = dt,
    method = mixed(50, 20)
  )
  summary(m_sm1)
  # (Intercept)  ***
  # smok_status2
  # smok_status3 ***
  # smok_status4 ***
}

# mod_min <- gamlss(
#   sbp ~ 1,
#   family = distr_nam,
#   weights = dt_trn$wt_nurse,
#   data = dt_trn,
#   method = mixed(5, 100)
# )
#
# sbp_model <- stepTGDAll.A(
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

sbp_model <- gamlss(
  sbp ~ log(year + 100) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha) + pcat(smok_status)),
  ~log(year + 100) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha) + pcat(smok_status),
  ~log(year + 100) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha) + pcat(smok_status),
  ~log(year + 100) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha) + pcat(smok_status),
  family = distr_nam,
  weights = dt$wt_nurse,
  data = dt,
  method = mixed(400, 400),
  control = con1
)

# sbp_model <- update(sbp_model, control = con1)
# tt <- chooseDist(sbp_model, type = "realplus", trace = TRUE, data = sbp_model$dt)
sbp_model$data <- copy(dt)

qsave(sbp_model, "./secure_data/lifecourse_models/sbp_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(sbp_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd), ethnicity = unique(dt$ethnicity),
              sha = unique(dt$sha), smok_status = unique(dt$smok_status))
newdata[, age := scale(age_int, 52.1, 17.1)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(sbp_model, .SD, data = sbp_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]
write_fst(newdata, "./inputs/exposure_distributions/sbp_table.fst", 100L)

print("Table saved")

# tt <- copy(dt)
# tt[, .N, by = .(agegrp10, sex, qimd)][, table(N)]
# tt[, sbp_rank := pctl_rank(sbp, "random"), by = .(agegrp10, sex, qimd)]
# tt[, age := as.integer(round(age * 16.8 + 50.8))]
# tt[newdata, on=.NATURAL, nomatch = NULL, `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu, tau = i.tau)]
# tt[sbp_rank == 0, sbp_rank := tt[sbp_rank > 0, min(sbp_rank)]]
# tt[sbp_rank == 1, sbp_rank := tt[sbp_rank < 1, max(sbp_rank)]]
#
# tt[, sbp := qBCPEo(sbp_rank, mu, sigma, nu, tau)]
#
#
# tt[, type := "Modelled"]
# dt[, type := "Observed"]
# zz <- rbind(tt, dt, fill = TRUE)
# zz[sbp > 50, sbp := 50]
# zz[sbp < 16, sbp := 16]
# zz <- zz[between(sbp, 17, 49)]
# zz[, weight := wt_nurse / sum(wt_nurse), by = .(type)]
# reldist_diagnostics(zz[type == "Observed", sbp],
#                     zz[type == "Modelled", sbp],
#                     zz[type == "Observed", weight],
#                     zz[type == "Modelled", weight],
#                     main = expression(bold(SBP ~ (kg / m ^ {
#                       2
#                     }))),
#                     100)

if (diagnostics) {
  sbp_model <- qread("./secure_data/lifecourse_models/sbp_model.qs")

  wp(sbp_model)
  wp(sbp_model, xvar = age)

  plot(sbp_model)
}

if (plots) {
  xlab_nam <- expression(bold(SBP ~ (mmHg)))
  dt[, age := age * 17.1 + 52.1] # descale
  sbp_model_tbl <- read_fst("./inputs/exposure_distributions/sbp_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt, sbp_model_tbl, 50, "sbp", paste0("my_q", distr_nam)
    )[between(sbp, quantile(sbp, 0.01), quantile(sbp, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_nurse / sum(wt_nurse), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/SBP_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", sbp],
                        zz[type == "Modelled", sbp],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        100)
    dev.off()
  })
  future(plot_synthpop_val(zz, sbp, "agegrp10", "wt_nurse", "SBP by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, sbp, "year", "wt_nurse", "SBP by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, sbp, "qimd", "wt_nurse", "SBP by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, sbp, "sha", "wt_nurse", "SBP by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, sbp, "ethnicity", "wt_nurse", "SBP by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, sbp, "smok_status", "wt_nurse", "SBP by smoking status", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, sbp, c("year", "agegrp10"), "wt_nurse", "SBP by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, sbp, c("year", "qimd"), "wt_nurse", "SBP by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, sbp, c("year", "sha"), "wt_nurse", "SBP by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, sbp, c("year", "ethnicity"), "wt_nurse", "SBP by year and ethnicity", xlab_nam, FALSE, FALSE))
}

