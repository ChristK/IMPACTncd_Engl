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
# For ages 30 to 90
univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- TRUE
seed                <- 43L


if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("qs", "fst", "MASS", "splines", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multiprocess)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}
source("./preparatory_work/aux_functions.R", local = TRUE)

dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 30, 90), .(
  education, age, agegrp10, sex, qimd, ethnicity, sha, wt_int, year)]
)
dt[, age := scale(age, 54.52, 15.28)]
dt[, education := ordered(education)]
set.seed(seed)

if (univariate_analysis) {
  # age
  age_scaled <- scale(30:90, 54.52, 15.28)
  dt[, .(education_mean = wtd.mean(as.integer(education), weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, education_mean, ylim = c(1, 7))]

  m_age1 <- polr(
    education ~ age,
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age1, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(7, 1e5, TRUE, x))), col = "blue1")

  m_age2 <- polr(
    education ~ poly(age, 2),
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age2, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(7, 1e5, TRUE, x))), col = "red1")

  m_age3 <- polr(
    education ~ ns(age, 4),
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age3, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(7, 1e5, TRUE, x))), col = "green1")

  setDT(BIC(m_age1, m_age2, m_age3), keep.rownames = TRUE, key = "BIC")[] # m_age3

  # year
  dt[, .(education_mean = wtd.mean(as.integer(education), weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, education_mean, xlim = c(3, 30), ylim = c(1, 7))]

  m_year1 <- polr(
    education ~ year,
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_year1 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(7, 1e5, TRUE, x))), col = "blue1")

  m_year2 <- polr(
    education ~ log(year),
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_year2 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(7, 1e5, TRUE, x))), col = "red1")

  setDT(BIC(m_year1, m_year2), keep.rownames = TRUE, key = "BIC")[] # m_year1

  # sex
  dt[, .(education_mean = wtd.mean(as.integer(education), weight = wt_int)), keyby = .(sex)
     ][, scatter.smooth(sex, education_mean, ylim = c(1, 7))]

}


education_model <- polr(
  education ~ year +  ns(age, 4) + qimd + sex + sha + ethnicity,
  weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = TRUE
)

mod_min <- polr(
  education ~ 1,
  weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = T
)

education_model2 <- stepAIC(mod_min,
                             scope = list(
                               lower = ~ 1,
                               upper = ~ (
                                 year + ns(age, 4) + sex + qimd +
                                   ethnicity + sha
                               ) ^ 3
                             ),
                             direction = "both",
                             k = log(nrow(dt))
)
AIC(education_model, education_model2)
AIC(education_model, education_model2, k = log(nrow(dt)))

education_model <- education_model2
education_model$data <- copy(dt)

qsave(education_model, "./secure_data/lifecourse_models/education_model.qs", preset = "high")
print("Model saved")

trms <- all.vars(formula(education_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 30:90, sex = unique(dt$sex), qimd = unique(dt$qimd),
              sha = unique(dt$sha), ethnicity = unique(dt$ethnicity))
newdata[, age := scale(age_int, 54.52, 15.28)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("ed", 1:7)) := data.table(rowCumsums(predict(education_model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata)
newdata[, age := age_int]
newdata[, c("age_int", "ed7") := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/education_table.fst", 100L)

print("Table saved")

if (plots) {
  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())
  xlab_nam <- expression(bold(Education ~ ("1 = Higher")))
  tbl <- read_fst("./inputs/exposure_distributions/education_table.fst", as.data.table = TRUE)

  education_model <- qread("./secure_data/lifecourse_models/education_model.qs")

  zz <- clone_dt(education_model$data, 10)
  zz[, education := NULL]
  zz[, age := round(age * 15.28 + 54.52)]
  zz[, rank_edu := runif(.N)]
  nam <- intersect(names(zz), names(tbl))
  zz[tbl, education :=  (rank_edu > ed1) + (rank_edu > ed2) +
       (rank_edu > ed3) + (rank_edu > ed4) + (rank_edu > ed5) + (rank_edu > ed6) + 1L,
     on = nam]
  zz[, `:=` (type = "Modelled", education = factor(
    education,
    levels = 1:7,
    labels = c(
      "NVQ4/NVQ5/Degree or equiv",
      "Higher ed below degree",
      "NVQ3/GCE A Level equiv",
      "NVQ2/GCE O Level equiv",
      "NVQ1/CSE other grade equiv",
      "Foreign/other",
      "No qualification"
    ), ordered = TRUE), .id = NULL)]
  zz[, rank_edu := NULL]
  zz <- rbind(zz, education_model$data[, type := "Observed"])
  zz[, education := as.integer(education) ]

  future({
    dir.create("./validation/synthpop_models", FALSE)
    zz[, weight := wt_int/sum(wt_int), by = type]
    png(
      "./validation/synthpop_models/Education_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", education],
                        zz[type == "Modelled", education],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = expression(bold(Education ~ ("1 = Higher"))),
                        discrete = TRUE)
    dev.off()
  })

  future(plot_synthpop_val(zz, education, "agegrp10", "wt_int", "Education by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, education, "year", "wt_int", "Education by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, education, "qimd", "wt_int", "Education by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, education, "sha", "wt_int", "Education by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, education, "ethnicity", "wt_int", "Education by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, education, c("year", "agegrp10"), "wt_int", "Education by year and agegroup", xlab_nam, FALSE, FALSE))
}



