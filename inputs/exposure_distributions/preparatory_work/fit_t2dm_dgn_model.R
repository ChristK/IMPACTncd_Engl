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
seed                <- 26L


if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}

if (!require(IMPACTncdEnglmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/")
  library(IMPACTncdEnglmisc)
}

dependencies(c("Rcpp", "dqrng", "fst", "qs", "gamlss", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multiprocess)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}
# sourceCpp("./preparatory_work/MN_distribution.cpp", cacheDir = "./.CppCache/")


dt <- na.omit(HSE_ts[wt_blood > 0 & between(age, 20, 90) & !between(year, 4, 5),
                     .(dm, dm_dgn, year, age, agegrp10, sex, qimd, bmi,
                       ethnicity, sha, wt_blood)])
dt[, dm := as.integer(as.character(dm))]
dt[, dm_dgn := as.integer(as.character(dm_dgn))]
dt <- dt[dm == 1L, ]
dt[, table(dm_dgn)]
dt[, age := scale(age, 62.5, 13.3)]
set.seed(seed)

if (file.exists("./preparatory_work/marginal_distr_dm_dgn.qs")) {
  marg_distr <- qread("./preparatory_work/marginal_distr_dm_dgn.qs")
} else {
  marg_distr <- distr_best_fit(dt, "dm_dgn", "wt_blood", "binom")
  qsave(marg_distr, "./preparatory_work/marginal_distr_dm_dgn.qs")
}
head(marg_distr$fits)
distr_nam <- names(marg_distr$fits[1])
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.


if (univariate_analysis) {
  age_scaled <- scale(18:90, 62.5, 13.3)
  dt[, .(dm_dgn_mean = wtd.mean(dm_dgn, weight = wt_blood)), keyby = .(age)
     ][, scatter.smooth(age, dm_dgn_mean)]

  m_age1 <- gamlss(
    dm_dgn ~ pb(age),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age1, dt, data.frame("age" = age_scaled)), col = "blue1")

  m_age2 <- gamlss(
    dm_dgn ~ age,
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age2, dt, data.frame("age" = age_scaled)), col = "red1")

  m_age3 <- gamlss(
    dm_dgn ~ cs(age),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(age_scaled, mean_predictAll(m_age3, dt, data.frame("age" = age_scaled)), col = "green1")
  setDT(GAIC(m_age1, m_age2, m_age3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_age1
  setDT(GAIC(m_age1, m_age2, m_age3, k = 2), TRUE)[, setorder(.SD, AIC)] # AIC m_age1

  dt[, .(dm_dgn_mean = wtd.mean(dm_dgn, weight = wt_blood)), keyby = .(year)
     ][, plot(year, dm_dgn_mean, xlim = c(3, 40), ylim = c(0, 0.9))]

  m_year1 <- gamlss(
    dm_dgn ~ year,
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(50, 20)
  )
  lines(3:40, mean_predictAll(m_year1, dt, data.frame("year" = 3:40)), col = "blue1")

  m_year2 <- gamlss(
    dm_dgn ~ log(year),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year2, dt, data.frame("year" = 3:40)), col = "red1")

  m_year3 <- gamlss(
    dm_dgn ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(5, 100)
  )
  lines(3:40, mean_predictAll(m_year3, dt, data.frame("year" = 3:40)), col = "green1")

  setDT(GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))), TRUE)[, setorder(.SD, AIC)] # BIC m_year32
  setDT(GAIC(m_year1, m_year2, m_year3, k = 2), TRUE)[, setorder(.SD, -AIC)] # AIC m_year3 = m_year2
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

  dt[, .(dm_dgn_mean = wtd.mean(dm_dgn, weight = wt_blood)), keyby = .(bmi = round(bmi))
     ][, scatter.smooth(bmi, dm_dgn_mean)]

  m_bmi1 <- gamlss(
    dm_dgn ~ pb(bmi),
    family = distr_nam,
    weights = dt$wt_blood,
    data = dt,
    method = mixed(20, 20)
  )
  lines(15:70, mean_predictAll(m_bmi1, dt, data.frame("bmi" = 15:70)), col = "blue1")


}


dm_dgn_model <- gamlss(
  dm_dgn ~ year * age * pcat(qimd) * (bmi + sex +  pcat(ethnicity) + pcat(sha)),
  family = distr_nam,
  weights = dt$wt_blood,
  data = dt,
  method = mixed(400, 400),
  control = con1
)

dm_dgn_model <- stepGAIC(dm_dgn_model, direction = "backward")

dm_dgn_model$data <- copy(dt)

qsave(dm_dgn_model, "./secure_data/lifecourse_models/dm_dgn_model.qs", preset = "high")

print("Model saved")

trms <- all.vars(formula(dm_dgn_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, qimd = unique(dt$qimd), sex = unique(dt$sex), sha = unique(dt$sha),
              ethnicity = unique(dt$ethnicity), bmi = unique(round(clamp(dt$bmi, 18, 50), -1)))
newdata[, age := scale(age_int, 62.5, 13.3)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu") := predictAll(dm_dgn_model, .SD, data = dm_dgn_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]

# NOTE dt has no chinese ethnicity. I will assume chinese have the same probability with white ethnicity
tt <- newdata[ethnicity == "white"]
tt[, ethnicity := "chinese"]
newdata <- rbind(newdata, tt)
write_fst(newdata, "./inputs/exposure_distributions/dm_dgn_table.fst", 100L)

print("Table saved")


# Diagnostics -------------------------------------------------------------
if (diagnostics) {
  dm_dgn_model <- qread("./secure_data/lifecourse_models/dm_dgn_model.qs")
  dm_dgn_model$data <- NULL
  wp(dm_dgn_model)
  wp(dm_dgn_model, xvar = ~dt$age)

  plot(dm_dgn_model)
}

# Plots -------------------------------------------------------------
if (plots) {
  xlab_nam <- expression(bold(T2DM~(dgn)))
  dt[, age := age * 13.3 + 62.5] # descale
  dt[, bmi := round(clamp(dt$bmi, 18, 50), -1)]
  dm_dgn_model_tbl <- read_fst("./inputs/exposure_distributions/dm_dgn_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt[age >= 20L], dm_dgn_model_tbl, 10, "dm_dgn", paste0("q", distr_nam))

  # [between(dm_dgn, quantile(dm_dgn, 0.01), quantile(dm_dgn, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_blood / sum(wt_blood), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/T2DM_dgn_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", dm_dgn],
                        zz[type == "Modelled", dm_dgn],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        2/3, TRUE)
    dev.off()
  })
  future(plot_synthpop_val(zz, dm_dgn, "agegrp10", "wt_blood",
                           "T2DM (dgn) by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dgn, "year", "wt_blood",
                           "T2DM (dgn) by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dgn, "qimd", "wt_blood",
                           "T2DM (dgn) by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dgn, "sha", "wt_blood",
                           "T2DM (dgn) by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dgn, "ethnicity", "wt_blood",
                           "T2DM (dgn) by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dgn, c("year", "agegrp10"), "wt_blood",
                           "T2DM (dgn) by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dgn, c("year", "qimd"), "wt_blood",
                           "T2DM (dgn) by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dgn, c("year", "sha"), "wt_blood",
                           "T2DM (dgn) by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dgn, c("year", "ethnicity"), "wt_blood",
                           "T2DM (dgn) by year and ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, dm_dgn, "bmi", "wt_blood",
                           "T2DM (dgn) by bmi", xlab_nam, FALSE, FALSE))

}

