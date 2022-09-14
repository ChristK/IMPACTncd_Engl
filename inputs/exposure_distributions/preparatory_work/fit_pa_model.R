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
# For ages 20 to 90
univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- TRUE
seed                <- 30L


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
source("./preparatory_work/aux_functions.R", local = TRUE)

dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90) & year > 5, .(
  active_days, year, age, agegrp10, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, age := scale(age, 49.6, 17)]
dt[, active_days := ordered(round(active_days))]
set.seed(seed)

if (univariate_analysis) {
  # age
  age_scaled <- scale(20:90, 49.6, 17)
  dt[, .(active_days_mean = wtd.mean(as.integer(active_days), weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, active_days_mean, ylim = c(0, 7))]

  m_age1 <- polr(
    active_days ~ ns(age, 4),
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age1, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "blue1")

  m_age2 <- polr(
    active_days ~ poly(age, 6),
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age2, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "red1")

  m_age3 <- polr(
    active_days ~ bs(age, 4),
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age3, type = "p", newdata = data.frame(age = age_scaled))
  lines(age_scaled, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "green1")

  setDT(BIC(m_age1, m_age2, m_age3), keep.rownames = TRUE, key = "BIC")[] # m_age3

  # year (This is likely meaningless as we project quintiles)
  dt[, .(active_days_mean = wtd.mean(as.integer(active_days), weight = wt_int)), keyby = .(year)
     ][, plot(year, active_days_mean, xlim = c(3, 30), ylim = c(1, 7))]

  m_year1 <- polr(
    active_days ~ year,
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_year1 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "blue1")

  m_year2 <- polr(
    active_days ~ year,
    weights = dt$wt_int,
    data = dt,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_year2 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "red1")

  setDT(BIC(m_year1, m_year2), keep.rownames = TRUE, key = "BIC")[] # equal I will accept the linear

}


active_days_model <- polr(
  active_days ~ year * bs(age, 4) *  qimd + (sex + sha + ethnicity),
  weights = dt$wt_int,
  data = dt,
  method = "logistic",
  Hess = TRUE
) ## Better AIC but worst BAIC

# mod_min <- polr(
#   active_days ~ 1,
#   weights = dt$wt_int,
#   data = dt,
#   method = "logistic",
#   Hess = T
# )
#
# active_days_model2 <- stepAIC(mod_min,
#                              scope = list(
#                                lower = ~ 1,
#                                upper = ~ (
#                                  year + bs(age, 4) + sex + qimd +
#                                    ethnicity + sha
#                                ) ^ 3
#                              ),
#                              direction = "both",
#                              k = log(nrow(dt))
# )


active_days_model$data <- copy(dt)

qsave(active_days_model, "./secure_data/lifecourse_models/active_days_model.qs", preset = "high")
print("Model saved")

trms <- all.vars(formula(active_days_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd),
              sha = unique(dt$sha), ethnicity = unique(dt$ethnicity)) #
newdata[, age := scale(age_int, 49.6, 17)]

newdata <- split(newdata, by = "ethnicity")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("pa", 0:7)) := data.table(rowCumsums(predict(active_days_model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata)
newdata[, age := age_int]
newdata[, c("age_int", "pa7") := NULL]
write_fst(newdata, "./inputs/exposure_distributions/active_days_table.fst", 100L)

print("Table saved")

if (plots) {
  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())
  xlab_nam <- expression(bold(Active ~ days ~ (d/week)))

  tbl <- read_fst("./inputs/exposure_distributions/active_days_table.fst", as.data.table = TRUE)

  active_days_model <- qread("./secure_data/lifecourse_models/active_days_model.qs")

  zz <- clone_dt(active_days_model$data, 10)
  zz[, active_days := NULL]
  zz[, age := round(age*17+49.6)]
  zz[, rank_pa := runif(.N)]
  nam <- intersect(names(zz), names(tbl))
  zz[tbl, active_days := (rank_pa > pa0) + (rank_pa > pa1) + (rank_pa > pa2) +
        (rank_pa > pa3) + (rank_pa > pa4) + (rank_pa > pa5) + (rank_pa > pa6) + 1L,
      on = nam]
  zz[, `:=` (
    type = "Modelled",
    active_days = factor(
      active_days,
      levels = 0:7,
      labels = 0:7,
      ordered = TRUE
    ),
    .id = NULL
  )]
  zz[, rank_pa := NULL]
  zz <- rbind(zz, active_days_model$data[, type := "Observed"])
  zz[, active_days := as.integer(as.character(active_days))]


  future({
    dir.create("./validation/synthpop_models", FALSE)
    zz[, weight := wt_int/sum(wt_int), by = type]
    png(
      "./validation/synthpop_models/Active_days_rel_dist.png",
      3840,
      2160,
      pointsize = 48

    )
    reldist_diagnostics(zz[type == "Observed", active_days],
                        zz[type == "Modelled", active_days],
                        zz[type == "Observed", weight],
                        zz[type == "Modelled", weight],
                        main = xlab_nam,
                        discrete = TRUE)
    dev.off()
  })

  future(plot_synthpop_val(zz, active_days, "agegrp10", "wt_int", "Active days by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "year", "wt_int", "Active days by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "sex", "wt_int", "Active days by sex", xlab_nam, FALSE, FALSE))

  future(plot_synthpop_val(zz, active_days, "qimd", "wt_int", "Active days by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "sha", "wt_int", "Active days by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "ethnicity", "wt_int", "Active days by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, c("year", "agegrp10"), "wt_int", "Active days by year and agegroup", xlab_nam, FALSE, FALSE))
}



