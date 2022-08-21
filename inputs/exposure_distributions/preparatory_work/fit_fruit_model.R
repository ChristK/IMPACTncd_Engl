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
plots               <- FALSE
seed                <- 32L



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
sourceCpp("./preparatory_work/SICHEL_distribution.cpp", cacheDir = "./.CppCache/")


dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90) & frtpor <= 30 & year < 14,
                     .(frtpor, year, age, agegrp10, sex, qimd, ethnicity, sha, wt_int)])
dt[, frtpor := round(frtpor)]
dt[, age := scale(age, 50.6, 17.4)]
set.seed(seed)
# lns <- sample(nrow(dt), round(nrow(dt) * 0.8))
# dt_trn   <- dt[lns] # train dataset
# dt_crv   <- dt[!lns]  # cross-validation dataset

if (file.exists("./preparatory_work/marginal_distr_frtpor.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_frtpor.qs")
} else {
  marg_distr <- distr_best_fit(dt, "frtpor", "wt_int", "count")
  qsave(marg_distr, "./preparatory_work/marginal_distr_frtpor.qs")
}
head(marg_distr$fits)
# distr_validation(marg_distr, dt[between(frtpor, 0, 30), .(var = frtpor, wt = wt_int)],
#                  expression(bold(frtpor ~ (portions/d))), TRUE)

distr_nam <- names(marg_distr$fits[1]) # ZISICHEL
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.

# Univariate analysis ----------------------------------------
if (univariate_analysis) {
  age_scaled <- scale(20:90, 50.6, 17.4)
  dt[, .(frtpor_median = wtd.quantile(frtpor, weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, frtpor_median)]

  m_age1 <- gamlss(
    frtpor ~ pb(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "blue1")

  m_age2 <- gamlss(
    frtpor ~ poly(age, 3),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "red1")

  m_age3 <- gamlss(
    frtpor ~ cs(age, 2),
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

  dt[, .(frtpor_median = wtd.quantile(frtpor, weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, frtpor_median, xlim = c(3, 40), ylim = c(0, 3))]

  m_year1 <- gamlss(
    frtpor ~ year,
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = RS(50)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "blue1")

  m_year2 <- gamlss(
    frtpor ~ log(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = RS(50)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "red1")

  m_year3 <- gamlss(
    frtpor ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = RS(50)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "green1")

  GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))) # BIC m_year1
  GAIC(m_year1, m_year2, m_year3, k = 2) # AIC m_year1
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

}


frtpor_model <- gamlss(
  frtpor ~ year *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha)),
  ~year + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha),
  ~year + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha),
  ~year + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha) ,
  family = distr_nam,
  weights = dt$wt_int,
  data = dt,
  method = RS(800),
  control = con1
)

# frtpor_model <- update(frtpor_model, method = mixed(200, 200))
# tt <- chooseDist(frtpor_model, type = "realplus", trace = TRUE, data = dt)
frtpor_model$data <- copy(dt)

qsave(frtpor_model, "./secure_data/lifecourse_models/frtpor_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(frtpor_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd), ethnicity = unique(dt$ethnicity),
              sha = unique(dt$sha))
newdata[, age := scale(age_int, 50.6, 17.4)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(frtpor_model, .SD, data = frtpor_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/frtpor_table.fst", 100L)

print("Table saved")

# tt <- copy(dt)
# tt[, .N, by = .(agegrp10, sex, qimd)][, table(N)]
# tt[, frtpor_rank := pctl_rank(frtpor, "random"), by = .(agegrp10, sex, qimd)]
# tt[, age := as.integer(round(age * 16.8 + 50.8))]
# tt[newdata, on=.NATURAL, nomatch = NULL, `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu, tau = i.tau)]
# tt[frtpor_rank == 0, frtpor_rank := tt[frtpor_rank > 0, min(frtpor_rank)]]
# tt[frtpor_rank == 1, frtpor_rank := tt[frtpor_rank < 1, max(frtpor_rank)]]
#
# tt[, frtpor := qBCPEo(frtpor_rank, mu, sigma, nu, tau)]
#
#
# tt[, type := "Modelled"]
# dt[, type := "Observed"]
# zz <- rbind(tt, dt, fill = TRUE)
# zz[frtpor > 50, frtpor := 50]
# zz[frtpor < 16, frtpor := 16]
# zz <- zz[between(frtpor, 17, 49)]
# zz[, weight := wt_int / sum(wt_int), by = .(type)]
# reldist_diagnostics(zz[type == "Observed", frtpor],
#                     zz[type == "Modelled", frtpor],
#                     zz[type == "Observed", weight],
#                     zz[type == "Modelled", weight],
#                     main = expression(bold(frtpor ~ (g / d))),
#                     100)

if (diagnostics) {
  frtpor_model <- qread("./secure_data/lifecourse_models/frtpor_model.qs")

  wp(frtpor_model)
  wp(frtpor_model, xvar = age)

  plot(frtpor_model)
}

if (plots) {
  xlab_nam <- expression(bold(Fruit ~ intake ~ (portions/d)))
  dt[, age := age * 17.4 + 50.6] # descale
  frtpor_model_tbl <- read_fst("./inputs/exposure_distributions/frtpor_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt, frtpor_model_tbl, 10, "frtpor", paste0("my_q", distr_nam)
                        )[between(frtpor, quantile(frtpor, 0.01), quantile(frtpor, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_int / sum(wt_int), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/Fruit_intake_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", frtpor],
                        zz[type == "Modelled", frtpor],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        100)
    dev.off()
  })
  future(plot_synthpop_val(zz, frtpor, "agegrp10", "wt_int", "Fruit_intake by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, frtpor, "year", "wt_int", "Fruit_intake by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, frtpor, "qimd", "wt_int", "Fruit_intake by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, frtpor, "sha", "wt_int", "Fruit_intake by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, frtpor, "ethnicity", "wt_int", "Fruit_intake by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, frtpor, c("year", "agegrp10"), "wt_int", "Fruit_intake by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, frtpor, c("year", "qimd"), "wt_int", "Fruit_intake by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, frtpor, c("year", "sha"), "wt_int", "Fruit_intake by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, frtpor, c("year", "ethnicity"), "wt_int", "Fruit_intake by year and ethnicity", xlab_nam, FALSE, FALSE))
}

