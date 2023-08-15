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
seed                <- 44L



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
sourceCpp("./preparatory_work/BCT_distribution.cpp", cacheDir = "./.CppCache/")

dt <- na.omit(HSE_ts[wt_blood > 0 & between(age, 20, 90),
                     .(tchol, year, age, agegrp10, sex, qimd, ethnicity, sha, wt_blood)])
dt[, age := scale(age, 52, 16.6)]
set.seed(seed)


if (file.exists("./preparatory_work/marginal_distr_tchol.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_tchol.qs")
} else {
  marg_distr <- distr_best_fit(dt, "tchol", "wt_blood", "realplus")
  qsave(marg_distr, "./preparatory_work/marginal_distr_tchol.qs")
}
head(marg_distr$fits)
# distr_validation(marg_distr, dt[between(tchol, 16, 50), .(var = tchol, wt = wt_blood)],
#                  expression(bold(TCHOL ~ (mmHg))))

distr_nam <- names(marg_distr$fits[1]) # BCT
con1 <- gamlss.control(c.crit = 1e-2) # increase for faster exploratory analysis.

if (univariate_analysis) {
  age_scaled <- scale(20:90, 52, 16.6)
  dt[, .(tchol_median = wtd.quantile(tchol, weight = wt_blood)), keyby = .(age)
     ][, scatter.smooth(age, tchol_median)]

  m_age1 <- gamlss(
    tchol ~ pb(age),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "blue1")

  m_age2 <- gamlss(
    tchol ~ poly(age, 3),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "red1")

  m_age3 <- gamlss(
    tchol ~ cs(age),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age3, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "green1")
  GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))) # BIC m_age3
  GAIC(m_age1, m_age2, m_age3, k = 2) # AIC m_age1
  centiles(m_age1, xvar = dt$age)
  centiles(m_age2, xvar = dt$age)

  dt[, .(tchol_median = wtd.quantile(tchol, weight = wt_blood)), keyby = .(year)
     ][, scatter.smooth(year, tchol_median, xlim = c(3, 40), ylim = c(4, 6))]

  m_year1 <- gamlss(
    tchol ~ year,
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "blue1")

  m_year2 <- gamlss(
    tchol ~ log(year),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "red1")

  m_year3 <- gamlss(
    tchol ~ log(year + 100),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "green1")

  GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))) # BIC m_year3 = m_year1
  GAIC(m_year1, m_year2, m_year3, k = 2) # AIC m_year3 = m_year1
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)
}


tchol_model <- gamlss(
  tchol ~ log(year + 100) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha)),
  ~log(year + 100) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha),
  ~log(year + 100) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha),
  ~log(year + 100) + pb(age) + sex + pcat(qimd) +
    pcat(ethnicity) + pcat(sha),
  family = distr_nam,
  weights = dt$wt_blood,
  data = dt,
  method = RS(400),
  control = con1
)

# tchol_model <- update(tchol_model, control = con1, family = "BCT")
# tt <- chooseDist(tchol_model, type = "realplus", trace = TRUE, data = tchol_model$dt)
tchol_model$data <- copy(dt)

qsave(tchol_model, "./secure_data/lifecourse_models/tchol_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(tchol_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd), ethnicity = unique(dt$ethnicity),
              sha = unique(dt$sha))
newdata[, age := scale(age_int, 52, 16.6)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu", "tau") := predictAll(tchol_model, .SD, data = tchol_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/tchol_table.fst", 100L)

print("Table saved")

# tt <- copy(dt)
# tt[, .N, by = .(agegrp10, sex, qimd)][, table(N)]
# tt[, tchol_rank := pctl_rank(tchol, "random"), by = .(agegrp10, sex, qimd)]
# tt[, age := as.integer(round(age * 16.8 + 50.8))]
# tt[newdata, on=.NATURAL, nomatch = NULL, `:=` (mu = i.mu, sigma = i.sigma, nu = i.nu, tau = i.tau)]
# tt[tchol_rank == 0, tchol_rank := tt[tchol_rank > 0, min(tchol_rank)]]
# tt[tchol_rank == 1, tchol_rank := tt[tchol_rank < 1, max(tchol_rank)]]
#
# tt[, tchol := qBCPEo(tchol_rank, mu, sigma, nu, tau)]
#
#
# tt[, type := "Modelled"]
# dt[, type := "Observed"]
# zz <- rbind(tt, dt, fill = TRUE)
# zz[tchol > 50, tchol := 50]
# zz[tchol < 16, tchol := 16]
# zz <- zz[between(tchol, 17, 49)]
# zz[, weight := wt_blood / sum(wt_blood), by = .(type)]
# reldist_diagnostics(zz[type == "Observed", tchol],
#                     zz[type == "Modelled", tchol],
#                     zz[type == "Observed", weight],
#                     zz[type == "Modelled", weight],
#                     main = expression(bold(TCHOL ~ (kg / m ^ {
#                       2
#                     }))),
#                     100)

if (diagnostics) {
  tchol_model <- qread("./secure_data/lifecourse_models/tchol_model.qs")

  wp(tchol_model)
  wp(tchol_model, xvar = age)

  plot(tchol_model)
}

if (plots) {
  xlab_nam <- expression(bold(TCHOL ~ (mmol/l)))
  dt[, age := age * 16.6 + 52] # descale
  tchol_model_tbl <- read_fst("./inputs/exposure_distributions/tchol_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt, tchol_model_tbl, 50, "tchol", paste0("my_q", distr_nam)
    )[between(tchol, quantile(tchol, 0.01), quantile(tchol, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_blood / sum(wt_blood), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/TCHOL_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", tchol],
                        zz[type == "Modelled", tchol],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        100)
    dev.off()
  })
  future(plot_synthpop_val(zz, tchol, "agegrp10", "wt_blood", "TCHOL by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, tchol, "year", "wt_blood", "TCHOL by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, tchol, "qimd", "wt_blood", "TCHOL by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, tchol, "sha", "wt_blood", "TCHOL by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, tchol, "ethnicity", "wt_blood", "TCHOL by ethnicity", xlab_nam, FALSE, FALSE))
  # future(plot_synthpop_val(zz, tchol, "smok_status", "wt_blood", "TCHOL by smoking status", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, tchol, c("year", "agegrp10"), "wt_blood", "TCHOL by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, tchol, c("year", "qimd"), "wt_blood", "TCHOL by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, tchol, c("year", "sha"), "wt_blood", "TCHOL by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, tchol, c("year", "ethnicity"), "wt_blood", "TCHOL by year and ethnicity", xlab_nam, FALSE, FALSE))
}

