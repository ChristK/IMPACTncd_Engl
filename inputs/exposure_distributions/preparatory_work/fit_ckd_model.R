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
# For ages 20 to 90
univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- TRUE
seed                <- 20L


if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("qs", "fst", "MASS", "splines", "matrixStats", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE)
plan(multiprocess)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}
if (!require(IMPACTncdEnglmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  roxygen2::roxygenise("./Rpackage/IMPACTncd_Engl_model_pkg/")
  remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/", force = TRUE)
  library(IMPACTncdEnglmisc)
}

dt <- na.omit(HSE_ts[wt_blood > 0 & between(age, 20, 90), .(
  ckd, year, age, agegrp10, sex, qimd, ethnicity, sha, wt_blood)]
)
dt[, age := scale(age, 52.2, 16.4)]
dt[, ckd := ordered(ckd)]
set.seed(seed)

if (univariate_analysis) {
  # age
  age_scaled <- scale(20:90, 52.2, 16.4)
  dt[, .(ckd_mean = wtd.mean(as.integer(ckd), weight = wt_blood)), keyby = .(age)
     ][, scatter.smooth(age, ckd_mean, ylim = c(0, 7))]

  m_age1 <- polr(
    ckd ~ ns(age, 4),
    weights = dt$wt_blood,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age1, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "blue1")

  m_age2 <- polr(
    ckd ~ bs(age, 3),
    weights = dt$wt_blood,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age2, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "red1")

  m_age3 <- polr(
    ckd ~ bs(age, 2),
    weights = dt$wt_blood,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age3, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "green1")

  setDT(BIC(m_age1, m_age2, m_age3), keep.rownames = TRUE, key = "BIC")[] # m_age3
}


ckd_model <- polr(
  ckd ~ bs(age, 3) *  qimd + (sex + sha + ethnicity),
  weights = dt$wt_blood,
  data = dt,
  method = "logistic",
  Hess = TRUE
) ## Better AIC but worst BAIC

ckd_model2 <- stepAIC(ckd_model,
                             direction = "back",
                             k = log(nrow(dt))
)

AIC(ckd_model, ckd_model2)
BIC(ckd_model, ckd_model2)

ckd_model <- ckd_model2
ckd_model$data <- copy(dt)

qsave(ckd_model, "./secure_data/lifecourse_models/ckd_model.qs", preset = "high")
print("Model saved")

trms <- all.vars(formula(ckd_model))[-1] # -1 excludes dependent var
newdata <- CJ(age_int = 20:90) #
newdata[, age := scale(age_int, 52.2, 16.4)]
newdata[, c(paste0("ckd", 0:4)) := data.table(rowCumsums(predict(ckd_model, type = "p", newdata = .SD))), .SDcols = trms]

newdata[, age := age_int]
newdata[, c("age_int", "ckd4") := NULL]
write_fst(newdata, "./inputs/exposure_distributions/ckd_table.fst", 100L)

print("Table saved")

if (plots) {
  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())
  xlab_nam <- expression(bold(CKD))

  tbl <- read_fst("./inputs/exposure_distributions/ckd_table.fst", as.data.table = TRUE)

  ckd_model <- qread("./secure_data/lifecourse_models/ckd_model.qs")

  zz <- clone_dt(ckd_model$data, 10)
  zz[, ckd := NULL]
  zz[, age := round(age*16.4+52.2)]
  zz[, rank_pa := runif(.N)]
  nam <- intersect(names(zz), names(tbl))
  zz[tbl, ckd := (rank_pa > ckd0) + (rank_pa > ckd1) + (rank_pa > ckd2) +
        (rank_pa > ckd3) + 1L,
      on = nam]
  zz[, `:=` (
    type = "Modelled",
    ckd = factor(
      ckd,
      levels = 1:5,
      labels = 1:5,
      ordered = TRUE
    ),
    .id = NULL
  )]
  zz[, rank_pa := NULL]
  zz <- rbind(zz, ckd_model$data[, type := "Observed"])
  zz[, ckd := as.integer(as.character(ckd))]


  future({
    dir.create("./validation/synthpop_models", FALSE)
    zz[, weight := wt_blood/sum(wt_blood), by = type]
    png(
      "./validation/synthpop_models/CKD_rel_dist.png",
      3840,
      2160,
      pointsize = 48
    )
    reldist_diagnostics(zz[type == "Observed", ckd],
                        zz[type == "Modelled", ckd],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        discrete = TRUE)
    dev.off()
  })

  future(plot_synthpop_val(zz, ckd, "agegrp10", "wt_blood", "CKD by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ckd, "sex", "wt_blood", "CKD by sex", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ckd, "qimd", "wt_blood", "CKD by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ckd, "sha", "wt_blood", "CKD by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, ckd, "ethnicity", "wt_blood", "CKD by ethnicity", xlab_nam, FALSE, FALSE))
}



