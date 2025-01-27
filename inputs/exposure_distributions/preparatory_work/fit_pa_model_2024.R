## workHORSE is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - workHORSE: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## workHORSE is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

# Final distribution chosen: polr 


library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk
library(IMPACTncdEngl)
library(CKutils)
library(MASS)
library(future.apply)

# Set some variables we will use later in automation
distributions <- FALSE
univariable_analysis <- FALSE
diagnostics <- FALSE
plots <- FALSE
seed                <- 43L


# Read in the cleaned time series of HSE data (2003-2019)
HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)


ds <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90) & year > 5, .(
  active_days, year, age, agegrp10, sex, qimd, ethnicity_grp, sha, wt_int, urban_rural)]
)

ds[, active_days := ordered(round(active_days))]
ds[, `:=`(
  sex = factor(sex),
  qimd = factor(qimd),
  sha = factor(sha), 
  urban_rural = factor(urban_rural), 
  ethnicity_grp = factor(ethnicity_grp) # 5 groups (more granular not available)
)]



summary(ds)

set.seed(seed)

if (univariable_analysis) {
  # age
  ds[, .(active_days_mean = wtd.mean(as.integer(active_days), weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, active_days_mean, ylim = c(0, 7))]

  age <- c(18:90)
  
  m_age0 <- polr(
    active_days ~ pb(age),
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age0, type = "p", newdata = data.frame(age = age))
  lines(age, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "pink")
  
    m_age1 <- polr(
    active_days ~ ns(age, 4),
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age1, type = "p", newdata = data.frame(age = age))
  lines(age, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "blue1")

  m_age2 <- polr(
    active_days ~ poly(age, 6),
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age2, type = "p", newdata = data.frame(age = age))
  lines(age, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "red1")

  m_age3 <- polr(
    active_days ~ bs(age, 4),
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age3, type = "p", newdata = data.frame(age = age))
  lines(age, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "green1")

#  AIC(m_age0, m_age1, m_age2,m_age3)
  BIC(m_age0, m_age1, m_age2, m_age3) #m_age3

  # year (This is likely meaningless as we project quintiles)
  ds[, .(active_days_mean = wtd.mean(as.integer(active_days), weight = wt_int)), keyby = .(year)
     ][, plot(year, active_days_mean, xlim = c(3, 30), ylim = c(1, 7))]

  m_year1 <- polr(
    active_days ~ year,
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_year1 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "blue1")

  m_year2 <- polr(
    active_days ~ year,
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_year2 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(8, 1e5, TRUE, x))), col = "red1")

  BIC(m_year1, m_year2)  # equal I will accept the linear

}


# active_days_model <- polr(
#   active_days ~ year * bs(age, 4) *  qimd + (sex + sha + ethnicity),
#   weights = ds$wt_int,
#   data = ds,
#   method = "logistic",
#   Hess = TRUE
# ) ## Better AIC but worst BAIC

mod_min <- polr(
  active_days ~ year + bs(age, 4) +  qimd + sex + sha + ethnicity_grp + urban_rural,
  weights = ds$wt_int,
  data = ds,
  method = "logistic",
  Hess = T
)
warnings() # check for warnings
qsave(mod_min, "./inputs/exposure_distributions/AH_test/pa_mod_min.qs", preset = "high")
mod_min <- qread( "./inputs/exposure_distributions/AH_test/pa_mod_min.qs")


pa_modelA <- stepAIC(
  mod_min,
  scope = list(
    lower = ~ year + bs(age, 4) +  qimd + sex + sha + ethnicity_grp + urban_rural,
    upper = ~ (year + bs(age, 4) +  qimd + sex + sha + ethnicity_grp + urban_rural)^2
  ),
  direction = "both",
  k = log(nrow(ds)),
  trace = TRUE  # This is so you can see all the models it has tried to fit
)
warnings()

pa_modelA

pa_modelA <- update(pa_modelA, control = gamlss.control(c.crit = 1e-3)) # back to default accuracy

pa_modelA$call


qsave(pa_modelA, "./inputs/exposure_distributions/AH_test/active_days_model.qs", preset = "high")
pa_modelA$hsedata <- copy(ds)
qsave(pa_modelA, "./secure_data/lifecourse_models/active_days_model.qs")

print("Model saved")

trms <- all.vars(formula(active_days_model))[-1] # -1 excludes dependent var
newdata <- CJ(  year = 3:50,
                age = 20:90,
                sex = unique(ds$sex),
                qimd = unique(ds$qimd),
                sha = unique(ds$sha),
                ethnicity_grp = unique(ds$ethnicity_grp),
                urban_rural = unique(ds$urban_rural)
)

newdata <- split(newdata, by = "ethnicity_grp")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("pa", 0:7)) := data.table(rowCumsums(predict(active_days_model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata)
write_fst(newdata, "./inputs/exposure_distributions/AH_test/active_days_table.fst", 100L)

print("Table saved")

if (plots) {
  dependencies(c("ggplot2", "cowplot"))
  theme_set(theme_cowplot())
  xlab_nam <- expression(bold(Active ~ days ~ (d/week)))

  tbl <- read_fst("./inputs/exposure_distributions/AH_test/active_days_table.fst", as.data.table = TRUE)

  active_days_model <- qread("./inputs/exposure_distributions/AH_test/pa_modelA.qs")

  zz <- as.data.table(sapply(ds, rep.int, times=10, simplify = F ))
  #zz <- clone_ds(active_days_model$data, 10)
  zz[, active_days := NULL]
  zz[, rank_pa := runif(.N)]
  nam <- intersect(names(zz), names(tbl))
  zz[tbl, active_days := (rank_pa > pa0) + (rank_pa > pa1) + (rank_pa > pa2) +
        (rank_pa > pa3) + (rank_pa > pa4) + (rank_pa > pa5) + (rank_pa > pa6) + 1L,
      on = nam]
  zz[, `:=` (
    type = "Modelled",
    active_days = factor(
      active_days,
      levels = 1:8,
      labels = 0:7,
      ordered = TRUE
    ),
    .id = NULL
  )]
  zz[, rank_pa := NULL]
  zz <- rbind(zz, ds[, type := "Observed"])
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
  future(plot_synthpop_val(zz, active_days, "qimd", "wt_int", "Active days by QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "sex", "wt_int", "Active days by sex", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "sha", "wt_int", "Active days by SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "ethnicity_grp", "wt_int", "Active days by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, "urban_rural", "wt_int", "Active days by rurality", xlab_nam, FALSE, FALSE))
  
  future(plot_synthpop_val(zz, active_days, c("year", "agegrp10"), "wt_int", "Active days by year and agegroup", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, c("year", "qimd"), "wt_int", "Active days by year and QIMD", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, c("year", "sha"), "wt_int", "Active days by year and SHA", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, c("year", "ethnicity_grp"), "wt_int", "Active days by year and ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, active_days, c("year", "urban_rural"), "wt_int", "Active days by year and rurality", xlab_nam, FALSE, FALSE))
  
  

}



