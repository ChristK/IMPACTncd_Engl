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
seed                <- 44L


if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("qs", "fst", "MASS", "splines", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE)
plan(multiprocess)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}
source("./preparatory_work/aux_functions.R", local = TRUE)

dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90), .(
  income, year, age, agegrp10, sex, qimd, education, ethnicity, sha, wt_int)]
)
dt[, age := scale(age, 49.96, 16.98)]
dt[, income := ordered(income)]
set.seed(seed)

if (univariate_analysis) {
  # age
  age_scaled <- scale(20:90, 49.96, 16.98)
  dt[, .(income_mean = wtd.mean(as.integer(income), weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, income_mean, ylim = c(1, 5))]

  m_age1 <- polr(
    income ~ ns(age, 8),
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age1, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "blue1")

  m_age2 <- polr(
    income ~ poly(age, 6),
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age2, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "red1")

  m_age3 <- polr(
    income ~ bs(age, 8),
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age3, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "green1")

  setDT(BIC(m_age1, m_age2, m_age3), keep.rownames = TRUE, key = "BIC")[] # m_age1

  # year (This is likely meaningless as we project quintiles)
  dt[, .(income_mean = wtd.mean(as.integer(income), weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, income_mean, xlim = c(3, 30), ylim = c(1, 7))]

  m_year1 <- polr(
    income ~ year,
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_year1 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "blue1")
}


income_model <- polr(
  income ~ ns(age, 8) + qimd + sex + sha + education + ethnicity,
  weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = TRUE
)

mod_min <- polr(
  income ~ 1,
  weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = T
)

income_model2 <- stepAIC(mod_min,
                          scope = list(
                            lower = ~ 1,
                            upper = ~ (
                              ns(age, 8) +
                                sex + qimd + education +
                                ethnicity + sha
                            ) ^ 2
                          ),
                          direction = "both",
                          k = log(nrow(dt))
                          )

AIC(income_model, income_model2)
AIC(income_model, income_model2, k = log(nrow(dt)))

income_model <- income_model2
income_model$data <- copy(dt)

qsave(income_model, "./secure_data/lifecourse_models/income_model.qs", preset = "high")
print("Model saved")

trms <- all.vars(formula(income_model))[-1] # -1 excludes dependent var
newdata <- CJ(age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd),
              sha = unique(dt$sha), education = unique(dt$education), ethnicity = unique(dt$ethnicity))
newdata[, age := scale(age_int, 49.96, 16.98)]

newdata <- split(newdata, by = "education")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("inc", 1:5)) := data.table(rowCumsums(predict(income_model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata)
newdata[, age := age_int]
newdata[, c("age_int", "inc5") := NULL]
write_fst(newdata, "./inputs/exposure_distributions/income_table.fst", 100L)

print("Table saved")

if (plots) {
  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())
  xlab_nam <- expression(bold(Income ~ ("5 = Higher")))
  tbl <- read_fst("./inputs/exposure_distributions/income_table.fst", as.data.table = TRUE)
  income_model <- qread("./secure_data/lifecourse_models/income_model.qs")

  zz <- clone_dt(income_model$data, 10)
  zz[, income := NULL]
  zz[, age := round(age * 16.98 + 49.96)]
  zz[, rank_inc := runif(.N)]
  nam <- intersect(names(zz), names(tbl))
  zz[tbl, income :=  (rank_inc > inc1) + (rank_inc > inc2) +
       (rank_inc > inc3) + (rank_inc > inc4) + 1L,
     on = nam]

  zz[, `:=` (
    type = "Modelled",
    income = factor(
      income,
      levels = 1:5,
      labels = c("1 Highest", "2", "3", "4", "5 Lowest"),
      ordered = TRUE
    ),
    .id = NULL
  )]
  zz[, rank_inc := NULL]
  zz <- rbind(zz, income_model$data[, type := "Observed"])
  zz[, income := as.integer(income) ]


  future({
    dir.create("./validation/synthpop_models", FALSE)
    zz[, weight := wt_int/sum(wt_int), by = type]
    png(
      "./validation/synthpop_models/Income_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", income],
                        zz[type == "Modelled", income],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = expression(bold(Income ~ ("5 = Higher"))),
                        discrete = TRUE)
    dev.off()
  })

  future(plot_synthpop_val(zz, income, "agegrp10", "wt_int", "Income by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, income, "year", "wt_int", "Income by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, income, "qimd", "wt_int", "Income by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, income, "sha", "wt_int", "Income by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, income, "ethnicity", "wt_int", "Income by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, income, "education", "wt_int", "Income by education", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, income, c("education", "agegrp10"), "wt_int", "Income by education and agegroup", xlab_nam, FALSE, FALSE))
}



