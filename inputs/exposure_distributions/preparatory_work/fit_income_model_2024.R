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

library(IMPACTncdEngl)
library(data.table)
library(gamlss)
library(scales) # for alpha
library(reldist)
library(qs) # save R objects to disk
library(fst) # save R tabular objects to disk
library(MASS) 
library(splines)
library(future.apply)
library(CKutils)

# Set some variables we will use later in automation
univariable_analysis <- FALSE
plots <- TRUE
seed                <- 44L



HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)

ds <- na.omit(HSE_ts[wt_int > 0 & between(age, 20, 90), .(
  income, year, age, agegrp10, sex, qimd, education, ethnicity_grp, sha,urban_rural, wt_int)]
)
ds[, income := ordered(income)]
set.seed(seed)

if (univariable_analysis) {
  # age
  ds[, .(income_mean = wtd.mean(as.integer(income), weight = wt_int)), keyby = .(age)
     ][, scatter.smooth(age, income_mean, ylim = c(1, 5))]

  m_age1 <- polr(
    income ~ ns(age, 8),
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age1, type = "p", newdata = data.frame(age = 20:90))
  lines(20:90, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "blue1")

  m_age2 <- polr(
    income ~ poly(age, 6),
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age2, type = "p", newdata = data.frame(age = 20:90))
  lines(20:90, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "red1")

  m_age3 <- polr(
    income ~ bs(age, 8),
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_age3, type = "p", newdata = data.frame(age = 20:90))
  lines(20:90, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "green1")

  setDT(BIC(m_age1, m_age2, m_age3), keep.rownames = TRUE, key = "BIC")[] # m_age1

  # year (This is likely meaningless as we project quintiles)
  ds[, .(income_mean = wtd.mean(as.integer(income), weight = wt_int)), keyby = .(year)
     ][, scatter.smooth(year, income_mean, xlim = c(3, 30), ylim = c(1, 7))]

  m_year1 <- polr(
    income ~ year,
    weights = ds$wt_int,
    data = ds,
    method = "logistic",
    Hess = TRUE
  )
  tt <- predict(m_year1 , type = "p", newdata = data.frame(year = 3:30))
  lines(3:30, apply(tt, 1, function(x) mean(sample(5, 1e5, TRUE, x))), col = "blue1")
}


income_model <- polr(
  income ~ ns(age, 8) + qimd + sex + sha + education + ethnicity_grp + urban_rural,
  weights = ds$wt_int,
  data = ds,
  method = "logistic",
  Hess = TRUE
)

mod_min <- polr(
  income ~ 1,
  weights = ds$wt_int,
  data = ds,
  method = "logistic",
  Hess = T
)

income_model2 <- stepAIC(mod_min,
                          scope = list(
                            lower = ~ 1,
                            upper = ~ (
                              ns(age, 8) +
                                sex + qimd + education +
                                ethnicity_grp + sha + urban_rural
                            ) ^ 2
                          ),
                          direction = "both",
                          k = log(nrow(ds))
                          )

AIC(income_model, income_model2)
AIC(income_model, income_model2, k = log(nrow(ds)))

income_model <- income_model2
income_model$data <- copy(ds)

qsave(income_model, "./secure_data/lifecourse_models/income_model.qs", preset = "high")
print("Model saved")

trms <- all.vars(formula(income_model))[-1] # -1 excludes dependent var
newdata <- CJ(age= 20:90, 
              sex = unique(ds$sex), 
              qimd = unique(ds$qimd),
              sha = unique(ds$sha),
              education = unique(ds$education), 
              ethnicity_grp = unique(ds$ethnicity_grp), 
              urban_rural = unique(ds$urban_rural))

newdata <- split(newdata, by = "education")

newdata <- # assignment necessary! Copies of data.tables are happening
  future_lapply(newdata, function(x)
    x[, c(paste0("inc", 1:5)) := data.table(rowCumsums(predict(income_model, type = "p", newdata = .SD))), .SDcols = trms],
    future.packages = c("MASS", "splines", "matrixStats"))
newdata <- rbindlist(newdata)
newdata[, c( "inc5") := NULL]

kc <- sort(setdiff(names(newdata), c("mu", "sigma", "nu", "tau")))
kc <- kc[order(match(kc, "year"))]
setcolorder(newdata, kc)
setkeyv(newdata, kc)

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
  future(plot_synthpop_val(zz, income, "ethnicity_grp", "wt_int", "Income by ethnicity", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, income, "education", "wt_int", "Income by education", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, income, "urban_rural", "wt_int", "Income by rurality", xlab_nam, FALSE, FALSE))
  future(plot_synthpop_val(zz, income, c("education", "agegrp10"), "wt_int", "Income by education and agegroup", xlab_nam, FALSE, FALSE))
}



