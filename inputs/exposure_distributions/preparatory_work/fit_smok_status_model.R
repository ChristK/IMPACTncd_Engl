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
seed                <- 40L


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
# source("./preparatory_work/aux_functions.R", local = TRUE)
# sourceCpp("./preparatory_work/MN_distribution.cpp", cacheDir = "./.CppCache/")


dt <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90),
  .(smok_status, year, age, agegrp10, sex, qimd,
    ethnicity, sha, wt_int)])
dt[, smok_status := as.integer(smok_status)]
dt[, age := scale(age, 50.7, 17.4)]
set.seed(seed)

distr_nam <- "MN4"
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.

if (univariate_analysis) {
  age_scaled <- scale(20:90, 50.7, 17.4)
  dt[, .(smok_status_median = wtd.quantile(smok_status, weight = wt_int)), keyby = .(age)
    ][, scatter.smooth(age, smok_status_median)]

  m_age1 <- gamlss(
    smok_status ~ pb(age),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age1, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "blue1")

  m_age2 <- gamlss(
    smok_status ~ poly(age, 3),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(20, 20)
  )
  lines(centiles.pred(m_age2, xname = "age", xvalues = age_scaled, cent = 50, data = dt), col = "red1")

  m_age3 <- gamlss(
    smok_status ~ cs(age),
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

  dt[, .(smok_status_median = wtd.quantile(smok_status, weight = wt_int)), keyby = .(year)
    ][, scatter.smooth(year, smok_status_median, xlim = c(3, 40), ylim = c(0, 5))]

  m_year1 <- gamlss(
    smok_status ~ pb(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(50, 20)
  )
  lines(centiles.pred(m_year1, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "blue1")

  m_year2 <- gamlss(
    smok_status ~ log(year),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year2, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "red1")

  m_year3 <- gamlss(
    smok_status ~ log(year - 2),
    family = distr_nam,
    weights = dt$wt_int,
    data = dt,
    method = mixed(5, 100)
  )
  lines(centiles.pred(m_year3, xname = "year", xvalues = 3:40, cent = 50, data = dt), col = "green1")

  GAIC(m_year1, m_year2, m_year3, k = log(nrow(dt))) # BIC m_year2
  GAIC(m_year1, m_year2, m_year3, k = 2) # AIC m_year2
  centiles(m_year1, xvar = dt$age)
  centiles(m_year2, xvar = dt$age)

}


smok_status_model <- gamlss(
  smok_status ~ log(year) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha)),
  ~log(year) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha)),
  ~log(year) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha)),
  ~log(year) *  pb(age) * pcat(qimd) * (sex +  pcat(ethnicity) + pcat(sha)),
  family = distr_nam,
  weights = dt$wt_int,
  data = dt,
  method = mixed(400, 400),
  control = con1
)
# smok_status_model2 <-
#   stepGAIC(
#     smok_status_model,
#     direction = "backward",
#     parallel = "multicore",
#     ncpus = 20L
#   )

# smok_status_model <- update(smok_status_model, method = mixed(200, 200))
# tt <- chooseDist(smok_status_model, type = "realplus", trace = TRUE, data = dt)
smok_status_model$data <- copy(dt)

qsave(smok_status_model, "./secure_data/lifecourse_models/smok_status_model.qs", preset = "high")
# smok_status_model <- qread("./secure_data/lifecourse_models/smok_status_model.qs")
print("Model saved")

trms <- all.vars(formula(smok_status_model))[-1] # -1 excludes dependent var
newdata <- CJ(year = 3:50, age_int = 20:90, sex = unique(dt$sex), qimd = unique(dt$qimd), ethnicity = unique(dt$ethnicity),
  sha = unique(dt$sha))
newdata[, age := scale(age_int, 50.7, 17.4)]

newdata <- split(newdata, by = "year")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c("mu", "sigma", "nu") := predictAll(smok_status_model, .SD, data = smok_status_model$data), .SDcols = trms])
newdata <- rbindlist(newdata)
setattr(newdata, "distribution", distr_nam)
newdata[, age := age_int]
newdata[, age_int := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

write_fst(newdata, "./inputs/exposure_distributions/smok_status_table.fst", 100L)

print("Table saved")

if (diagnostics) {
  smok_status_model <- qread("./secure_data/lifecourse_models/smok_status_model.qs")
  smok_status_model$data <- NULL
  wp(smok_status_model)
  wp(smok_status_model, xvar = ~smok_status_model$data$age)

  plot(smok_status_model)
}

if (plots) {
  xlab_nam <- expression(bold(Smoking~Status))
  dt[, age := age * 17.4 + 50.7] # descale
  smok_status_model_tbl <- read_fst("./inputs/exposure_distributions/smok_status_table.fst", as.data.table =  TRUE)
  zz <-
    validate_gamlss_tbl(dt, smok_status_model_tbl, 50, "smok_status", paste0("my_q", distr_nam))[between(smok_status, quantile(smok_status, 0.01), quantile(smok_status, 0.99))]

  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())

  future({
    zz[, weight := wt_int / sum(wt_int), by = .(type)]
    dir.create("./validation/synthpop_models", FALSE)
    png(
      "./validation/synthpop_models/Smoking_status_rel_dist.png",
      3840,
      2160,
      pointsize = 48
    )
    reldist_diagnostics(zz[type == "Observed", smok_status],
      zz[type == "Modelled", smok_status],
      zz[type == "Observed", weight],
      zz[type == "Modelled", weight],
      main = xlab_nam,
      2/3, TRUE)
    dev.off()
  })
  future(plot_synthpop_val(zz, smok_status, "agegrp10", "wt_int",
    "Smoking status by agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_status, "sex", "wt_int",
    "Smoking status by sex", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_status, "year", "wt_int",
    "Smoking status by year", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_status, "qimd", "wt_int",
    "Smoking status by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_status, "sha", "wt_int",
    "Smoking status by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_status, "ethnicity", "wt_int",
    "Smoking status by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_status, c("year", "agegrp10"), "wt_int",
    "Smoking status by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_status, c("year", "qimd"), "wt_int",
    "Smoking status by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_status, c("year", "sha"), "wt_int",
    "Smoking status by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, smok_status, c("year", "ethnicity"), "wt_int",
    "Smoking status by year and ethnicity", xlab_nam, FALSE, FALSE))

  # zz[!is.na(smok_status), smok_never := fifelse(smok_status == 1, 1, 0)]
  # zz[!is.na(smok_status), smok_active := fifelse(smok_status == 4, 1, 0)]
  # tt <- zz[, mean(smok_never), keyby = .(year, sex, type)]
  # gg <- ggplot(
  #   tt,
  #   aes(
  #     x = year,
  #     y = V1,
  #     col = type,
  #     fill = type
  #   )
  # ) +
  #   geom_point(size = 1,
  #     alpha = 5 / 5,
  #     show.legend = FALSE) +
  #   geom_line(size = 1, alpha = 5 / 5) +
  #   # geom_ribbon(alpha = 1 / 5,
  #   #   linetype = 0,
  #   #   show.legend = FALSE) +
  #   scale_x_continuous(name = "Year") +
  #   scale_y_continuous(name = "proportion") +
  #   # ggtitle(filename) +
  #   scale_color_viridis_d() +
  #   scale_fill_viridis_d() +
  #   theme(
  #     legend.position="bottom",
  #     legend.title = element_blank(),
  #     strip.text.y = element_text(angle = 0)) +
  #   facet_grid(sex ~ .)
}

