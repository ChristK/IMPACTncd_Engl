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



# preample ---------------------------------------------------
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(c("ggplot2", "cowplot", "fst", "gamlss", "gamlss.tr", "gamlss.add",
               "future", "future.apply", "reldist", "data.table", "minerva"), FALSE, FALSE, TRUE, FALSE)
theme_set(theme_cowplot())

options(
  future.fork.enable = TRUE, # enable fork in Rstudio
  future.globals.maxSize = +Inf,
  future.globals.onReference = "ignore",
  "scipen" = 999,
  "digits" = 4
)
plan(multiprocess, workers = 20L)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}
source("./preparatory_work/aux_functions.R", local = TRUE)


# Correlation exploratory analysis ----
# Using Maximal Information Coefficient (MIC) to account for linear and non-linear correlations
# https://www.freecodecamp.org/news/how-machines-make-predictions-finding-correlations-in-complex-data-dfd9f0d87889/
# I will use bootstrap to calculate CI for MIC but most importantly to account for
# the weighted data. Note that because I use pairwise.complete.obs the sample size in each iteration
# for each paiwise MIC is different because of missing values. Also note that this method
# does not account for partial correlations i.e. cor(x,y|z). Furthermore, sample size
# for bootstrap was set to 1e4 for efficiency but there are vars that ~90% are missing
HSE_ts[, smok_status_cont := as.integer(smok_status)]
HSE_ts[, income_cont := as.integer(income)]
HSE_ts[, education_cont := as.integer(education)]

if (file.exists("./preparatory_work/correlation_structure.csv")) {
  correlation_structure <- fread("./preparatory_work/correlation_structure.csv")
  correlation_structure_post_2010 <- fread("./preparatory_work/correlation_structure_post_2010.csv")
} else {
  source("./preparatory_work/maximal_information_coefficient.R", local = TRUE)
}

# Education exploration ---------------------
correlation_structure[, .SD, .SDcols = c("rn", "education_cont")
                      ][, setorderv(.SD, "education_cont", order = -1)]


# Income exploration ---------------------
correlation_structure[, .SD, .SDcols = c("rn", "income_cont")
                      ][, setorderv(.SD, "income_cont", order = -1)]
plot_cor("education_cont", "income_cont", "agegrp20", "wt_int", HSE_ts)

# Smoking status exploration ---------------------
correlation_structure[, .SD, .SDcols = c("rn", "smok_status_cont")
                      ][, setorderv(.SD, "smok_status_cont", order = -1)]
plot_cor("smok_status_cont", "bmi", "agegrp20", "wt_int", HSE_ts)


# Alcohol exploration (based on average 7d intake, totalwu) ---------------------
correlation_structure_post_2010[, .SD, .SDcols = c("rn", "totalwu")
                      ][, setorderv(.SD, "totalwu", order = -1)]

plot_cor("alcohol", "totalwu", "agegrp20", "wt_int", HSE_ts)
plot_cor("smok_status_cont", "totalwu", "agegrp20", "wt_int", HSE_ts)
plot_cor("bmi", "totalwu", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("income_cont", "totalwu", "agegrp20", "wt_int", HSE_ts)
plot_cor("education_cont", "totalwu", "agegrp20", "wt_int", HSE_ts)
plot_cor("education_cont", "totalwu", "income_cont", "wt_int", HSE_ts)
plot_cor("income_cont", "totalwu", "education_cont", "wt_int", HSE_ts)
plot_cor("sbp", "totalwu", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("tchol", "totalwu", "agegrp20", "wt_blood", HSE_ts)
plot_cor("hba1c", "totalwu", "agegrp20", "wt_blood", HSE_ts)
plot_cor("frtpor", "totalwu", "agegrp20", "wt_int", HSE_ts)
plot_cor("vegpor", "totalwu", "agegrp20", "wt_int", HSE_ts)

# Fruit exploration ---------------------
correlation_structure[, .SD, .SDcols = c("rn", "frtpor")
                                ][, setorderv(.SD, "frtpor", order = -1)]

plot_cor("bmi", "frtpor", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("alcohol", "frtpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("sbp", "frtpor", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("smok_status_cont", "frtpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("vegpor", "frtpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("education_cont", "frtpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("education_cont", "frtpor", "income_cont", "wt_int", HSE_ts)
plot_cor("income_cont", "frtpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("income_cont", "frtpor", "education_cont", "wt_int", HSE_ts)
plot_cor("active_days", "frtpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("hba1c", "frtpor", "agegrp20", "wt_blood", HSE_ts)
plot_cor("tchol", "frtpor", "agegrp20", "wt_blood", HSE_ts)

# Veg exploration ---------------------
correlation_structure[, .SD, .SDcols = c("rn", "vegpor")
                      ][, setorderv(.SD, "vegpor", order = -1)]

plot_cor("bmi", "vegpor", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("frtpor", "vegpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("alcohol", "vegpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("sbp", "vegpor", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("education_cont", "vegpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("education_cont", "vegpor", "income_cont", "wt_int", HSE_ts)
plot_cor("active_days", "vegpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("tchol", "vegpor", "agegrp20", "wt_blood", HSE_ts)
plot_cor("hba1c", "vegpor", "agegrp20", "wt_blood", HSE_ts)
plot_cor("income_cont", "vegpor", "agegrp20", "wt_int", HSE_ts)
plot_cor("income_cont", "vegpor", "education_cont", "wt_int", HSE_ts)
plot_cor("smok_status_cont", "vegpor", "agegrp20", "wt_int", HSE_ts)

# Active_days exploration ---------------------
correlation_structure[, .SD, .SDcols = c("rn", "active_days")
                      ][, setorderv(.SD, "active_days", order = -1)]

plot_cor("bmi", "active_days", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("sbp", "active_days", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("alcohol", "active_days", "agegrp20", "wt_int", HSE_ts)
plot_cor("education_cont", "active_days", "agegrp20", "wt_int", HSE_ts)
plot_cor("education_cont", "active_days", "income_cont", "wt_int", HSE_ts)
plot_cor("hba1c", "active_days", "agegrp20", "wt_blood", HSE_ts)
plot_cor("income_cont", "active_days", "agegrp20", "wt_int", HSE_ts)
plot_cor("income_cont", "active_days", "education_cont", "wt_int", HSE_ts)
plot_cor("frtpor", "active_days", "agegrp20", "wt_int", HSE_ts)
plot_cor("tchol", "active_days", "agegrp20", "wt_blood", HSE_ts)
plot_cor("vegpor", "active_days", "agegrp20", "wt_int", HSE_ts)
plot_cor("smok_status_cont", "active_days", "agegrp20", "wt_int", HSE_ts)

# BMI exploration & GAMLSS model ---------------------
correlation_structure[, .SD, .SDcols = c("rn", "bmi")][, setorderv(.SD, "bmi", order = -1)]

plot_cor("sbp", "bmi", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("hba1c", "bmi", "agegrp20", "wt_blood", HSE_ts)
plot_cor("tchol", "bmi", "agegrp20", "wt_blood", HSE_ts)
plot_cor("active_days", "bmi", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("frtpor", "bmi", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("alcohol", "bmi", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("vegpor", "bmi", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("smok_status_cont", "bmi", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("income_cont", "bmi", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("education_cont", "bmi", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("education_cont", "bmi", "income_cont", "wt_nurse", HSE_ts)
plot_cor("income_cont", "bmi", "education_cont", "wt_nurse", HSE_ts)

# SBP exploration & GAMLSS model ---------------------
correlation_structure[, .SD, .SDcols = c("rn", "sbp")][, setorderv(.SD, "sbp", order = -1)]

plot_cor("bmi", "sbp", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("tchol", "sbp", "agegrp20", "wt_blood", HSE_ts)
plot_cor("hba1c", "sbp", "agegrp20", "wt_blood", HSE_ts)
plot_cor("active_days", "sbp", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("education_cont", "sbp", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("education_cont", "sbp", "income_cont", "wt_nurse", HSE_ts)
plot_cor("income_cont", "sbp", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("income_cont", "sbp", "education_cont", "wt_nurse", HSE_ts)
plot_cor("alcohol", "sbp", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("frtpor", "sbp", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("smok_status_cont", "sbp", "agegrp20", "wt_nurse", HSE_ts)
plot_cor("vegpor", "sbp", "agegrp20", "wt_nurse", HSE_ts)

# Tchol exploration & GAMLSS model ---------------------
correlation_structure[, .SD, .SDcols = c("rn", "tchol")][, setorderv(.SD, "tchol", order = -1)]

plot_cor("bmi", "tchol", "agegrp20", "wt_blood", HSE_ts)
plot_cor("sbp", "tchol", "agegrp20", "wt_blood", HSE_ts)
plot_cor("alcohol", "tchol", "agegrp20", "wt_blood", HSE_ts)
plot_cor("hba1c", "tchol", "agegrp20", "wt_blood", HSE_ts)
plot_cor("education_cont", "tchol", "agegrp20", "wt_blood", HSE_ts)
plot_cor("education_cont", "tchol", "income_cont", "wt_blood", HSE_ts)
plot_cor("frtpor", "tchol", "agegrp20", "wt_blood", HSE_ts)
plot_cor("active_days", "tchol", "agegrp20", "wt_blood", HSE_ts)
plot_cor("income_cont", "tchol", "agegrp20", "wt_blood", HSE_ts)
plot_cor("income_cont", "tchol", "education_cont", "wt_blood", HSE_ts)
plot_cor("smok_status_cont", "tchol", "agegrp20", "wt_blood", HSE_ts)
plot_cor("vegpor", "tchol", "agegrp20", "wt_blood", HSE_ts)





dt <- na.omit(HSE_ts[between(age, 20, 90), .(
  bmi, year, age, agegrp10, sex, qimd, ethnicity, sha, wt_nurse)]
  )
marg_distr <- distr_best_fit(dt, "bmi", "wt_nurse", "realplus")

# resolved(futureOf(marg_distr))
# while (!resolved(futureOf(marg_distr))) {print("Not yet");Sys.sleep(5)}
# start_val <- gamlssML(bmi, GB2_bmi, weights = dt$wt_nurse,
# data = dt)

head(marg_distr$fits)
distr_validation(marg_distr, dt[between(bmi, 16, 50), .(var = bmi, wt = wt_nurse)],
                 expression(bold(BMI ~ (kg / m ^ {2}))))
#           Summary statistics   Measure Lower 95% CI Upper 95% CI p-value
# 1:    Overall change entropy  0.000958           NA           NA      NA
# 2:     Median effect entropy  0.000195           NA           NA      NA
# 3:      Shape effect entropy  0.000758           NA           NA      NA
# 4: Median polarization index  0.000565     -0.00551      0.00664   0.428
# 5:  Lower polarization index  0.003244     -0.00890      0.01538   0.300
# 6:  Upper polarization index -0.002138     -0.01429      0.01001   0.365

dt <- na.omit(HSE_ts[between(age, 20, 90), .(
  bmi, year, age, agegrp10, sex, qimd, ethnicity, sha, year, sbp, tchol, smok_status, wt_nurse, wt_blood)]
)
m1 <- gamlss(bmi ~ age + sex + qimd + ethnicity + sha, data = dt, weights = dt$wt_nurse,
       family = names(marg_distr$fits[1]), method = mixed(5, 100))
summary(m1)
m2 <- gamlss(bmi ~ age + sex + qimd + ethnicity + sha + sbp + tchol + smok_status,
             data = dt, weights = dt$wt_blood,
             family = names(marg_distr$fits[1]), method = mixed(5, 100))
summary(m2, robust = TRUE, save = TRUE)



con1 <- gamlss.control(c.crit = 1e-1) # reduced for faster exploratory analysis.
# nC <- floor(detectCores()/4)
m1 %<-% gamlss(bmi~pb(age) * sex, data = dt_trn, weights = dt_trn$wt_nurse,
             family = distr_nam, method = mixed(5, 100))
m2 %<-% gamlss(bmi~pvc(age, by = sex) + sex, data = dt_trn, weights = dt_trn$wt_nurse,
             family = distr_nam, method = mixed(5, 100))
m3 %<-% gamlss(bmi~pvc(age, by = sex) * sex, data = dt_trn, weights = dt_trn$wt_nurse,
             family = distr_nam, method = mixed(5, 100))
GAIC(m1, m2, m3, k = log(nrow(dt_trn))) # m2
GAIC(m1, m2, m3, k = 2) # m2

dt[, .(bmi_median = wtd.quantile(bmi, weight = wt_nurse)), keyby = .(age, sex)
   ][, scatter.smooth(age, bmi_median, col = sex)]
lines(20:90, centile_predictAll(m1, newdata = data.frame(age = 20:90, sex = 1)), col = "blue1")
lines(20:84, centile_predictAll(m1, newdata = data.frame(age = 20:84, sex = 2)), col = "blue4")
lines(20:84, centile_predictAll(m2, newdata = data.frame(age = 20:84, sex = 1)), col = "red1")
lines(20:84, centile_predictAll(m2, newdata = data.frame(age = 20:84, sex = 2)), col = "red4")
lines(20:84, centile_predictAll(m3, newdata = data.frame(age = 20:84, sex = 1)), col = "green1")
lines(20:84, centile_predictAll(m3, newdata = data.frame(age = 20:84, sex = 2)), col = "green4")

m4 %<-% gamlss(bmi~log(year) + qimd, data = dt_trn, weights = dt_trn$wt_nurse,
             family = distr_nam, method = mixed(5, 100))
m4r %<-% gamlss(bmi~log(year) + pcat(qimd), data = dt_trn, weights = dt_trn$wt_nurse,
               family = distr_nam, method = mixed(5, 100))
m5 %<-% gamlss(bmi~log(year) * pcat(qimd), data = dt_trn, weights = dt_trn$wt_nurse,
             family = distr_nam, method = mixed(5, 100))
m6 %<-% gamlss(bmi~pvc(year,  by = qimd), data = dt_trn, weights = dt_trn$wt_nurse,
               family = distr_nam, method = mixed(5, 100))
GAIC(m4, m4r, m5, m6, k = log(nrow(dt_trn))) # m4
GAIC(m4, m4r, m5, m6, k = 2) # m6

dt[, .(bmi_median = wtd.quantile(bmi, weight = wt_nurse)), keyby = .(year, qimd)
   ][, scatter.smooth(year, bmi_median, col = qimd, xlim = c(3, 40), ylim = c(25, 30))]
lines(3:40, centile_predictAll(m4, newdata = data.frame(year = 3:40, qimd = 1)), col = "red")
lines(3:40, centile_predictAll(m4, newdata = data.frame(year = 3:40, qimd = 2)), col = "red1")
lines(3:40, centile_predictAll(m4, newdata = data.frame(year = 3:40, qimd = 3)), col = "red2")
lines(3:40, centile_predictAll(m4, newdata = data.frame(year = 3:40, qimd = 4)), col = "red3")
lines(3:40, centile_predictAll(m4, newdata = data.frame(year = 3:40, qimd = 5)), col = "red4")

lines(3:40, centile_predictAll(m4r, newdata = data.frame(year = 3:40, qimd = 1)), col = "purple")
lines(3:40, centile_predictAll(m4r, newdata = data.frame(year = 3:40, qimd = 2)), col = "purple1")
lines(3:40, centile_predictAll(m4r, newdata = data.frame(year = 3:40, qimd = 3)), col = "purple2")
lines(3:40, centile_predictAll(m4r, newdata = data.frame(year = 3:40, qimd = 4)), col = "purple3")
lines(3:40, centile_predictAll(m4r, newdata = data.frame(year = 3:40, qimd = 5)), col = "purple4")

lines(3:40, centile_predictAll(m5, newdata = data.frame(year = 3:40, qimd = 1)), col = "green")
lines(3:40, centile_predictAll(m5, newdata = data.frame(year = 3:40, qimd = 2)), col = "green1")
lines(3:40, centile_predictAll(m5, newdata = data.frame(year = 3:40, qimd = 3)), col = "green2")
lines(3:40, centile_predictAll(m5, newdata = data.frame(year = 3:40, qimd = 4)), col = "green3")
lines(3:40, centile_predictAll(m5, newdata = data.frame(year = 3:40, qimd = 5)), col = "green4")

lines(3:40, centile_predictAll(m6, newdata = data.frame(year = 3:40, qimd = 1)), col = "blue")
lines(3:40, centile_predictAll(m6, newdata = data.frame(year = 3:40, qimd = 2)), col = "blue1")
lines(3:40, centile_predictAll(m6, newdata = data.frame(year = 3:40, qimd = 3)), col = "blue2")
lines(3:40, centile_predictAll(m6, newdata = data.frame(year = 3:40, qimd = 4)), col = "blue3")
lines(3:40, centile_predictAll(m6, newdata = data.frame(year = 3:40, qimd = 5)), col = "blue4")

m7 %<-% gamlss(bmi~ethnicity, data = dt_trn, weights = dt_trn$wt_nurse,
               family = distr_nam, method = mixed(5, 100))
m8 %<-% gamlss(bmi~pcat(ethnicity), data = dt_trn, weights = dt_trn$wt_nurse,
               family = distr_nam, method = mixed(5, 100))
getSmo(m8)
coef(getSmo(m8))
fitted(getSmo(m8))[1:10]
plot(getSmo(m8)) #
# After the fit a new factor is created  this factor has the reduced levels
levels(getSmo(m8)$factor)
GAIC(m7, m8, k = log(nrow(dt_trn))) # m8
GAIC(m7, m8, k = 2) # m8

m9 %<-% gamlss(bmi~sha, data = dt_trn, weights = dt_trn$wt_nurse,
               family = distr_nam, method = mixed(5, 100))
m10 %<-% gamlss(bmi~random(sha), data = dt_trn, weights = dt_trn$wt_nurse,
               family = distr_nam, method = mixed(5, 100))
GAIC(m9, m10, k = log(nrow(dt_trn))) # m9
GAIC(m9, m10, k = 2) # m9

m11 <- quote(gamlss(bmi~random(sha, df=p[1]),
                   sigma.fo= ~random(sha, df=p[2]),
                   nu.fo= ~random(sha, df=p[3]),
                  data = dt_trn, weights = dt_trn$wt_nurse, family=distr_nam,
                  method = mixed(5, 100),
                   control = gamlss.control(trace=FALSE)))
optimis2 %<-% find.hyper(m11,par=c(8,8,8), lower=c(0.001,0.001, 0.001),
           upper=c(10,10,10), pen=2)

mod1_min %<-% gamlss(bmi~log(year) + pb(age, by = sex) + sex + pcat(qimd) +
                      pcat(ethnicity) +  random(sha),
                     sigma.fo = ~log(year) + pb(age, by = sex) + sex + pcat(qimd) +
                         random(sha),
                     nu.fo = ~log(year) + pb(age, by = sex),
                       family = distr_nam, weights = dt_small$wt_nurse,
                       data = dt_small, method = mixed(5, 100))
mod2_min %<-% gamlss(bmi~1,
                     sigma.fo = ~1,
                     nu.fo = ~1,
                       family = distr_nam, weights = dt_small$wt_nurse,
                       data = dt_small, method = mixed(5, 100))
GAIC(mod1_min, mod1_max, k = log(nrow(dt_trn))) # m9
GAIC(mod1_min, mod1_max, k = 2) # m9

bmi_model %<-% stepTGDAll.A(
  mod2_min,
  scope = list(
    lower =  ~ 1,
    upper =  ~ (
      log(year) + pb(age, by = sex) + sex + pcat(qimd) +
        pcat(ethnicity) +  random(sha)
    ) ^ 2
  ),
  sigma.scope = list(
    lower =  ~ 1,
    upper =  ~ (
      log(year) + pb(age, by = sex) + sex + pcat(qimd) +
        pcat(ethnicity) +  random(sha)
    ) ^ 2
  ),
  nu.scope = list(lower =  ~ 1, upper =  ~ (log(year) + pb(age, by = sex) + pcat(qimd)) ^
                    2),
  newdata = dt_crv,
  parallel = "multicore",
  ncpus = 12L
)
resolved(futureOf(bmi_model))



mod1_small %<-% gamlss(bmi~pb(year) + pb(age) + sex + qimd + sha + pcat(ethnicity),
                     family = distr_nam, weights = dt_small$wt_nurse,
                     data = dt_small, method = mixed(5, 100))
mod2_small %<-% gamlss(bmi~pb(year) * (pb(age) + sex + qimd + sha + ethnicity),
                     family = distr_nam, weights = dt_small$wt_nurse,
                     data = dt_small, method = mixed(5, 100))
mod3_small %<-% gamlss(bmi~pb(year) * (pb(age) + sex + qimd + sha),
                       family = distr_nam, weights = dt_small$wt_nurse,
                       data = dt_small, method = mixed(5, 100))
mod4_small %<-% gamlss(bmi~pvc(year, by = qimd) * (pvc(age, by = sex) + sex +
                                                     pcat(ethnicity) + random(sha)) + qimd,
                       family = distr_nam, weights = dt_small$wt_nurse,
                       data = dt_small, method = mixed(5, 100))
GAIC(mod1_small, mod2_small, mod3_small, mod4_small,
     k = log(nrow(dt_small))) # mod1 has lower BIC
GAIC(mod1_small, mod2_small, mod3_small, mod4_small,
     k = 2) # mod4 has lower AIC

plot(mod4_small, par = newpar)
wp(mod4_small)
wp(mod2_small, dt_small$age)
# mean should be ~0, var ~1, skew ~0, kurt = ~3

mod41_small %<-% update(mod4_small,
                        sigma.fo = ~ pvc(year, by = qimd) * pvc(age, by = sex))
mod42_small %<-% update(mod4_small,
                        sigma.fo = ~ ~pvc(year, by = qimd) * (pvc(age, by = sex)
                                                              + sex + ethnicity) + qimd)
GAIC(mod4_small, mod41_small, mod42_small) # mod22 has lower AIC
plot(mod42_small, par = newpar)
wp(mod42_small)

mod222_small %<-% update(mod22_small, nu.fo = ~ pb(age, df = 4))
GAIC(mod22_small, mod222_small) # mod222 has lower AIC
plot(mod222_small, par = newpar)
wp(mod222_small)

mod2222_small %<-% stepGAIC(mod222_small)
GAIC(mod222_small, mod2222_small) # mod222 has lowerAIC
wp(mod2222_small)
term.plot(mod2222_small, se=TRUE, partial=TRUE)
zz <- validate_gamlss(dt_crv, mod2222_small, 50, dt_small)[between(bmi, 17, 45)]
reldist_diagnostics(zz[type == "Observed", bmi],
                    zz[type == "Modelled", bmi],
                    zz[type == "Observed", wt_nurse/sum(wt_nurse)],
                    zz[type == "Modelled", wt_nurse/sum(wt_nurse)],
                    main = expression(bold(BMI~(kg/m^{2}))),
                    100)
yy <- zz[sha == "2", ]
reldist_diagnostics(yy[type == "Observed", bmi],
                    yy[type == "Modelled", bmi],
                    yy[type == "Observed", wt_nurse/sum(wt_nurse)],
                    yy[type == "Modelled", wt_nurse/sum(wt_nurse)],
                    main = expression(bold(BMI~(kg/m^{2}))),
                    100)
bmi_model %<-%  update(mod2222_small, data = dt_trn, weights = dt_trn$wt_nurse)

  bmi ~ pb(year) * (pb(age) + sex + qimd + sha + ethnicity),
  sigma.fo = ~ pb(year) * (pb(age) + sex + qimd + sha),
  nu.fo = ~ pb(age, df = 4),
  family = distr_nam,
  weights = dt_trn$wt_nurse,
  data = dt_trn,
  method = mixed(5, 100)
)
bmi_model %<-%  gamlss(
  bmi ~ pb(year) * (pb(age) + sex + qimd + sha + ethnicity),
  sigma.fo = ~ pb(year) * (pb(age) + sex + qimd + sha),
  nu.fo = ~ pb(age, df = 4),
  family = distr_nam,
  weights = dt_trn$wt_nurse,
  data = dt_trn,
  method = mixed(5, 100)
)

bmi_model2 %<-% stepGAIC(bmi_model)

zz <- validate_gamlss(dt_crv, bmi_model, 50, dt_trn)[between(bmi, 17, 45)]
reldist_diagnostics(zz[type == "Observed", bmi],
                    zz[type == "Modelled", bmi],
                    zz[type == "Observed", wt_nurse/sum(wt_nurse)],
                    zz[type == "Modelled", wt_nurse/sum(wt_nurse)],
                    main = expression(bold(BMI~(kg/m^{2}))),
                    100)
plot(bmi_model)
wp(bmi_model)



# GAMLSS alcohol model ---------------------
# TODO Build a model to predict `totalwu` from `alcohol`
dt <- na.omit(HSE_ts[alcohol < 200, .(
alcohol, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, .(alcohol_median = wtd.quantile(alcohol, weight = wt_int)), keyby = .(year, sex)
   ][, plot(year, alcohol_median, col = sex, xlim = c(2, 15), ylim = c(0, 50))]
dt <- na.omit(HSE_ts[alcohol < 200 & totalwu < 200, .(
  alcohol, totalwu, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, .(totalwu_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(year, sex)
   ][, points(year, totalwu_median, col = sex)]

dt[, plot(alcohol, totalwu)]

plot(density(dt$alcohol, weights = dt$wt_int))
lines(density(dt$totalwu, weights = dt$wt_int), col = "red")

HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- na.omit(HSE_ts[totalwu < 100 & wt_int > 0, .(
  totalwu, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, age := scale(age, 50, 16.8)]
dt[, totalwu := round(totalwu)]
set.seed(43)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
marg_distr <- fitDistPred(dt_trn$totalwu, type = "realplus", extra = c("plogSSTZadj", "pGGZadj"),
                          weights = dt_trn$wt_int,
                            try.gamlss = TRUE, trace = TRUE, newdata = dt_crv$totalwu,
                            parallel = "multicore", ncpus = 12L)
marg_distr$fits # BNB family does not converge
# ZASICHEL    ZANBI ZISICHEL    ZINBF    ZINBF    ZINBI


params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
# params$n <- 1e4
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y <- y[y < 100]
y_wt <- rep(1/length(y), length(y))
reldist_diagnostics(dt[, totalwu], y, dt[, wt_int/sum(wt_int)], y_wt,
                    main = expression(bold(Alcohol~(g/d))), discrete = TRUE)

centile_predictAll <- function(gamlss_model, newdata, cent = 0.5) {
  stopifnot("gamlss" %in% class(gamlss_model))
  tt <- predictAll(gamlss_model, newdata)
  tt$p <- 0.5
  dist_nam <- attr(tt, "family")[[1]]
  out <- do.call(paste0("q", dist_nam), tt)
  return(out)
}

con1 <- gamlss.control(c.crit = 1e-2) # reduced for faster exploratory analysis.
setMKLthreads(1L)
distr_nam <- "ZINBF"
dt_trn[, age_scale := scale(age)]
m1 %<-% gamlss(totalwu~pvc(age_scale, by = sex), data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100))
m2 %<-% gamlss(totalwu~pvc(age_scale, by = sex) + sex, data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100))
m3 %<-% gamlss(totalwu~pvc(log(age), by = sex), data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100))
GAIC(m1, m2, m3, k = log(nrow(dt_trn))) # m1
GAIC(m1, m2, m3, k = 2) # m2

dt_trn[, .(totalwu_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(age, sex)
   ][, scatter.smooth(age, totalwu_median, col = sex)]
lines(20:84, centile_predictAll(m1, newdata = data.frame(age_scale = scale(20:84,
                                                                           attr(dt_trn$age_scale, "scaled:center"),
                                                                           attr(dt_trn$age_scale, "scaled:scale")), sex = 1)), col = "blue1")
lines(20:84, centile_predictAll(m1, newdata = data.frame(age_scale = scale(20:84,
                                                                           attr(dt_trn$age_scale, "scaled:center"),
                                                                           attr(dt_trn$age_scale, "scaled:scale")), sex = 2)), col = "blue4")
lines(20:84, centile_predictAll(m2, newdata = data.frame(age_scale = scale(20:84,
                                                                           attr(dt_trn$age_scale, "scaled:center"),
                                                                           attr(dt_trn$age_scale, "scaled:scale")), sex = 1)), col = "red1")
lines(20:84, centile_predictAll(m2, newdata = data.frame(age_scale = scale(20:84,
                                                                           attr(dt_trn$age_scale, "scaled:center"),
                                                                           attr(dt_trn$age_scale, "scaled:scale")), sex = 2)), col = "red4")
lines(20:84, centile_predictAll(m3, newdata = data.frame(age = 20:84, sex = 1)), col = "green1")
lines(20:84, centile_predictAll(m3, newdata = data.frame(age = 20:84, sex = 2)), col = "green4")

m4 %<-% gamlss(totalwu~year, data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100))
m5 %<-% gamlss(totalwu~log(year), data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100))
GAIC(m4, m5, k = log(nrow(dt_trn)))
GAIC(m4, m5, k = 2)

dt_trn[, .(totalwu_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(year, sex)
       ][, scatter.smooth(year, totalwu_median, col = sex)]
lines(3:14, centile_predictAll(m4, newdata = data.frame(year = 3:14)), col = "green1")
lines(3:14, centile_predictAll(m5, newdata = data.frame(year = 3:14)), col = "red")

m6 %<-% gamlss(totalwu~pcat(ethnicity, Lp = 1), data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100))
m7 %<-% gamlss(totalwu~pcat(ethnicity), data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100))
GAIC(m6, m7, k = log(nrow(dt_trn))) # m7
GAIC(m6, m7, k = 2) # m7
dt_trn[, .(totalwu_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(ethnicity, sex)
       ][, scatter.smooth(ethnicity, totalwu_median, col = sex)]
lines(1:9, centile_predictAll(m6, newdata = data.frame(ethnicity = 1:9)), col = "green1")
lines(1:9, centile_predictAll(m7, newdata = data.frame(ethnicity = 1:9)), col = "red")

m8 %<-% gamlss(totalwu~pcat(sha, Lp = 1), data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100))
m9 %<-% gamlss(totalwu~pcat(sha), data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100))
GAIC(m8, m9, k = log(nrow(dt_trn)))
GAIC(m8, m9, k = 2)
dt_trn[, .(totalwu_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(sha, sex)
       ][, scatter.smooth(sha, totalwu_median, col = sex)]
lines(1:10, centile_predictAll(m8, newdata = data.frame(sha = 1:10)), col = "green1")
lines(1:10, centile_predictAll(m9, newdata = data.frame(sha = 1:10)), col = "red")

con1 <- gamlss.control(c.crit = 1e-1)
m10 <- gamlss(totalwu~pcat(qimd, Lp = 1, start = 2.35), data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100), control = con1)
m11 %<-% gamlss(totalwu~pcat(qimd), data = dt_trn, weights = dt_trn$wt_int,
               family = distr_nam, method = mixed(5, 100), control = con1)
GAIC(m10, m11, k = log(nrow(dt_trn))) # m10
GAIC(m10, m11, k = 2) # m10
dt_trn[, .(totalwu_median = wtd.quantile(totalwu, weight = wt_int)), keyby = .(qimd, sex)
       ][, scatter.smooth(qimd, totalwu_median, col = sex)]
lines(1:5, centile_predictAll(m10, newdata = data.frame(qimd = 1:5)), col = "green1")
lines(1:5, centile_predictAll(m11, newdata = data.frame(qimd = 1:5)), col = "red")

m12 %<-% gamlss(totalwu~pcat(qimd, Lp = 1, start = 2.35) + pvc(age, by = sex), data = dt_trn, weights = dt_trn$wt_int,
                family = distr_nam, method = mixed(5, 100), control = con1)
resolved(futureOf(m12))

m1 %<-% gamlss(
  totalwu ~ ga( ~ s(log(year), age), method="REML"),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
m2 %<-% gamlss(
  totalwu ~ ga( ~ te(log(year), age, k=3), method="REML"),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
m3 %<-% gamlss(
  totalwu ~ ga( ~ ti(log(year), age, k=3), method="REML") + log(year) *  age,
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
m4 %<-% gamlss(
  totalwu ~ ga( ~ ti(log(year), age, k=3) , method="REML")+ log(year) *  pb(age),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
m5 %<-% gamlss(
  totalwu ~ ga( ~ s(log(year), age, by = sex)) + sex,
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
m7 %<-% gamlss(
  totalwu ~ log(year) * pb(age, by = sex) + sex,
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
GAIC(m1, m2, m3, m4, m5, m7, k = 0) # m5
GAIC(m1, m2, m3, m4, m5, m7, k = log(nrow(dt_trn)))

# Ordinal physical activity model ---------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- na.omit(HSE_ts[year > 5 & wt_int > 0, .(
  active_days, bmi, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
) # years before 2006 were excluded for more realistic projections
dt[, .(bmi_median = median(bmi)), keyby = active_days][, scatter.smooth(.SD)]
plot(density(dt$active_days, weights = dt$wt_int))
setkey(dt[, sum(wt_int), keyby = active_days], V1)[]

dt[, age := scale(age, 50, 16.8)]
dt[, bmi := scale(bmi, 27.4, 5.2)]
dt[, active_days := active_days/7]

set.seed(44)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset


set.seed(44)
lns <- sample(nrow(dt), nrow(dt) * 0.9)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
setMKLthreads(12L)

dt[, .(active_days_mean = wtd.mean(active_days, weight = wt_int)), keyby = .(age)
   ][, scatter.smooth(age, active_days_mean)]
lines(scale(20:84, 50, 16.8), predict(glm(active_days~age, weights = wt_int, data = dt),
                     newdata = data.frame(age = scale(20:84, 50, 16.8))), col = "red")
lines(scale(20:84, 50, 16.8), predict(glm(active_days~ns(age, 4), weights = wt_int, data = dt),
                     newdata = data.frame(age = scale(20:84, 50, 16.8))), col = "blue")

dt[, .(active_days_mean = wtd.mean(active_days, weight = wt_int)), keyby = .(year)
   ][, plot(year, active_days_mean, ylim = c(0, 5), xlim = c(6, 40))]
lines(5:40, predict(glm(active_days~year, weights = wt_int, data = dt),
                     newdata = data.frame(year = 5:40)), col = "red")
lines(5:40, predict(glm(active_days~log(year), weights = wt_int, data = dt),
                    newdata = data.frame(year = 5:40)),  col = "blue")

dt[, .(active_days_mean = wtd.mean(active_days, weight = wt_int)), keyby = .(bmi = round(bmi))
   ][, plot(bmi, active_days_mean, ylim = c(0, 5), xlim = c(-2.2, 5.4))]
lines(scale(16:55, 27.4, 5.2), predict(glm(active_days~bmi, weights = wt_int, data = dt),
                    newdata = data.frame(bmi = scale(16:55, 27.4, 5.2))), col = "red")
lines(scale(16:55, 27.4, 5.2), predict(glm(active_days~log(bmi), weights = wt_int, data = dt),
                    newdata = data.frame(bmi = scale(16:55, 27.4, 5.2))),  col = "green")
lines(scale(16:55, 27.4, 5.2), predict(glm(active_days~ns(bmi, 4), weights = wt_int, data = dt),
                     newdata = data.frame(bmi = scale(16:55, 27.4, 5.2))),  col = "blue")

dt[, .(active_days_mean = wtd.mean(active_days, weight = wt_int)), keyby = .(sha)
   ][, plot(sha, active_days_mean, ylim = c(0, 5))]
lines(1:10, predict(glm(active_days~sha, weights = wt_int, data = dt),
                     newdata = data.frame(sha = as.factor(1:10))), col = "red")
lines(1:10, predict(glm(active_days~pcat(sha, 8), weights = wt_int, data = dt),
                    newdata = data.frame(sha = as.factor(1:10))), col = "green")
lines(1:10, predict(glm(active_days~random(sha), weights = wt_int, data = dt),
                    newdata = data.frame(sha = as.factor(1:10))), col = "blue")

dt[, .(active_days_mean = wtd.mean(active_days, weight = wt_int)), keyby = .(ethnicity)
   ][, plot(ethnicity, active_days_mean, ylim = c(0, 5))]
lines(1:9, predict(glm(active_days~ethnicity, weights = wt_int, data = dt),
                    newdata = data.frame(ethnicity = as.factor(1:9))), col = "red")
lines(1:9, predict(glm(active_days~pcat(ethnicity, 8), weights = wt_int, data = dt),
                    newdata = data.frame(ethnicity = as.factor(1:9))), col = "green")
lines(1:9, predict(glm(active_days~random(ethnicity), weights = wt_int, data = dt),
                    newdata = data.frame(ethnicity = as.factor(1:9))), col = "blue")

dt[, .(active_days_mean = wtd.mean(active_days, weight = wt_int)), keyby = .(qimd)
   ][, plot(qimd, active_days_mean, ylim = c(0, 5))]
lines(1:5, predict(glm(active_days~qimd, weights = wt_int, data = dt),
                    newdata = data.frame(qimd = as.factor(1:5))), col = "red")
lines(1:5, predict(glm(active_days~pcat(qimd, 4), weights = wt_int, data = dt),
                    newdata = data.frame(qimd = as.factor(1:5))), col = "green")
lines(1:5, predict(glm(active_days~random(qimd), weights = wt_int, data = dt),
                    newdata = data.frame(qimd = as.factor(1:5))), col = "blue")

# GAMLSS ETS model ---------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)

dt <- na.omit(HSE_ts[wt_int > 0, .(
  ets, smok_status,
  agegrp, year, age, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, age := scale(age, 50, 16.8)]
# Calculate smoking prevalence to be used as externality
dt[, smok_prev := prop_if(smok_status == "4"), by = .(year, qimd, sha)]

set.seed(44)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
diagnostics <- FALSE

centile_predictAll <- function(gamlss_model, newdata, cent = 0.5) {
  stopifnot("gamlss" %in% class(gamlss_model))
  tt <- predictAll(gamlss_model, newdata)
  tt$p <- 0.5
  if ("y" %in% names(tt)) tt$y <- NULL
  dist_nam <- attr(tt, "family")[[1]]
  if (dist_nam %in% .gamlss.bi.list) {
    return(tt$mu)
  } else {
    out <- do.call(paste0("q", dist_nam), tt)
  }
}

dt[, .(ets_mean = wtd.mean(as.integer(as.character(ets)), weight = wt_int)), keyby = .(year)
   ][, plot(year, ets_mean, xlim = c(0, 35), ylim = c(0, 0.6))]
dt[qimd == 1, .(ets_mean = wtd.mean(as.integer(as.character(ets)), weight = wt_int)), keyby = .(year)
   ][, points(year, ets_mean, xlim = c(0, 35), ylim = c(0, 0.6), col = "blue")]
dt[qimd == 5, .(ets_mean = wtd.mean(as.integer(as.character(ets)), weight = wt_int)), keyby = .(year)
   ][, points(year, ets_mean, xlim = c(0, 35), ylim = c(0, 0.6), col = "red")]
lines(0:35, predict(glm(ets~exp(year) * (year >= 8), weights = wt_int, data = dt, family = binomial()),
                    newdata = data.frame(year = 0:35), type = "response"), col = "blue")
lines(0:35, predict(glm(ets~year * (year >= 6), weights = wt_int, data = dt, family = binomial()),
                    newdata = data.frame(year = 0:35), type = "response"), col = "green")
lines(0:35, predict(glm(ets~year * (year >= 8) * (year >= 6), weights = wt_int, data = dt,
                        family = binomial(link = "identity")),
                    newdata = data.frame(year = 0:35), type = "response"), col = "red")
summary(glm(ets~year * (year >= 8) * (year >= 6), weights = wt_int, data = dt, family = binomial(link = "identity")))

dt2 <- dt[, .(ets_prev = prop_if(ets == "1"), smok_prev = mean(smok_prev),
              wt_int = sum(wt_int)), by = .(year, qimd, sha)]
dt2[, scatter.smooth(smok_prev, ets_prev)]
lines(seq(0, 1, by = 0.01), predict(glm(ets_prev~log(smok_prev),
                                        weights = wt_int, data = dt2[smok_prev > 0], family = gaussian()),
                                    newdata = data.frame(smok_prev = seq(0, 1, by = 0.01)), type = "response"), col = "red")
lines(seq(0, 1, by = 0.01), predict(glm(ets_prev~poly(smok_prev, 2),
                                        weights = wt_int, data = dt2, family = gaussian()),
                                    newdata = data.frame(smok_prev = seq(0, 1, by = 0.01)), type = "response"), col = "blue")

distr_nam <- "BI"
setMKLthreads(1L)
m1 <- gamlss(
  ets ~ year * (year > 7) * (year > 5),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
m2 <- gamlss(
  ets ~ pb(year) * qimd,
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
dt[, .(ets_mean = wtd.mean(as.integer(as.character(ets)), weight = wt_int)), keyby = .(year)
   ][, plot(year, ets_mean, xlim = c(0, 35), ylim = c(0, 1))]
lines(1:30, predictAll(m1, newdata = data.frame(year = 1:30))$mu, col = "yellow")
lines(1:30, predictAll(m2, newdata = data.frame(year = 1:30, qimd = 1))$mu, col = "blue")
lines(1:30, predictAll(m2, newdata = data.frame(year = 1:30, qimd = 2))$mu, col = "green")
lines(1:30, predictAll(m2, newdata = data.frame(year = 1:30, qimd = 3))$mu, col = "black")
lines(1:30, predictAll(m2, newdata = data.frame(year = 1:30, qimd = 4))$mu, col = "purple")
lines(1:30, predictAll(m2, newdata = data.frame(year = 1:30, qimd = 5))$mu, col = "red")
GAIC(m1, m2, k = 2) # m2 better
GAIC(m1, m2, k = log(nrow(dt_trn))) # m2 better

m3 <- gamlss(
  ets ~ pb(smok_prev, df = 3) * pcat(qimd),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
dt[, .(ets_mean = wtd.mean(as.integer(as.character(ets)), weight = wt_int), smok_prev), keyby = .(round(smok_prev, 2))
   ][, plot(smok_prev, ets_mean, xlim = c(0, 1), ylim = c(0, 1))]
lines(seq(0, 1, by = 0.01), predictAll(m3, newdata = data.frame(smok_prev = seq(0, 1, by = 0.01)))$mu, col = "blue")
lines(seq(0, 1, by = 0.01), predictAll(m3, newdata = data.frame(smok_prev = seq(0, 1, by = 0.01), qimd = 1))$mu, col = "blue")
lines(seq(0, 1, by = 0.01), predictAll(m3, newdata = data.frame(smok_prev = seq(0, 1, by = 0.01), qimd = 2))$mu, col = "green")
lines(seq(0, 1, by = 0.01), predictAll(m3, newdata = data.frame(smok_prev = seq(0, 1, by = 0.01), qimd = 3))$mu, col = "black")
lines(seq(0, 1, by = 0.01), predictAll(m3, newdata = data.frame(smok_prev = seq(0, 1, by = 0.01), qimd = 4))$mu, col = "purple")
lines(seq(0, 1, by = 0.01), predictAll(m3, newdata = data.frame(smok_prev = seq(0, 1, by = 0.01), qimd = 5))$mu, col = "red")

# GAMLSS smoking initiation model ---------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- HSE_ts[wt_int > 0, .(
  wt_int, smok_init_age, smok_status,
  agegrp, year, age, sex, qimd, ethnicity, sha, smok_dur, smok_quit_yrs)]
dt[, wt_int := wt_int*10000/sum(wt_int), by = year] # make pop of each hse =10000
dt[smok_status == 4 & between(age - smok_init_age, 0, 3),
   `:=`(year = year - (age - smok_init_age))]
dt[smok_status == 4 & between(age - smok_init_age, 0, 3),
   `:=`(age = smok_init_age, smok_incid = 1L)]
dt[smok_status %in% c("2", "3") & between(smok_quit_yrs + smok_dur, 0, 3),
   `:=`(
     age  = age  - smok_quit_yrs - smok_dur,
     year = year - smok_quit_yrs - smok_dur,
     smok_incid = 1L
   )]
dt[smok_status == "1", smok_incid := 0L]

dt <- na.omit(dt[between(year, 3, 11) & age > 17, .(
  wt_int, smok_incid, agegrp, year, age, sex, qimd, ethnicity, sha)])
# dt[, age := scale(age, 50, 16.8)]

set.seed(45)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset

distr_nam <- "BI"
setMKLthreads(16L)
m1 <- gamlss(
  smok_incid ~ log(year),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
tt <- dt[, .(prop_if(smok_incid == 1)), by = year]
tt[, pr := predictAll(m1, newdata = as.data.frame(year))$mu]
tt[, plot(year, V1)]
tt[, lines(year, pr)]

m2 <- gamlss(
  smok_incid ~ I(1/age),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)

m3 <- gamlss(
  smok_incid ~ log(age),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
GAIC(m2, m3, k = 2)
GAIC(m2, m3, k = log(nrow(dt_trn)))

tt <- dt[, .(prop_if(smok_incid == 1)), keyby = .(age)]
tt[, pr := predictAll(m3, newdata = data.frame(age))$mu]
tt[, plot(age, V1)]
tt[, lines(age, pr)]

m4 <- gamlss(
  smok_incid ~ I(1/age) + log(year),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)

GAIC(m4, m3, k = 2)
GAIC(m4, m3, k = log(nrow(dt_trn)))


tt <- dt[, .(prop_if(smok_incid == 1)), keyby = .(age, year)]
tt[, pr := predictAll(m4, newdata = data.frame(age, year))$mu]
tt[, plot(age, V1)]
tt[, lines(age, pr)]


# GAMLSS smoking cessation model ---------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- HSE_ts[wt_int > 0, .(
  wt_int, smok_quit_yrs, smok_status,
  agegrp, year, age, sex, qimd, ethnicity, sha)]
dt[smok_status == "4", smok_cess := 0L]
dt[smok_quit_yrs == 1, smok_cess := 1L]

dt <- na.omit(dt[, .(
  wt_int, smok_cess, agegrp, year, age, sex, qimd, ethnicity, sha)])
# dt[, age := scale(age, 50, 16.8)]

set.seed(46)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset

distr_nam <- "BI"
setMKLthreads(16L)
m1 <- gamlss(
  smok_cess ~ log(year),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
tt <- dt[, .(prop_if(smok_cess == 1)), by = year]
tt[, pr := predictAll(m1, newdata = as.data.frame(year))$mu]

tt[, plot(year, V1)]
tt[, lines(year, pr)]

m2 <- gamlss(
  smok_cess ~ log(age),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)

m3 <- gamlss(
  smok_cess ~ pb(age),
  data = dt_trn,
  family = distr_nam,
  weights = dt_trn$wt_int,
  method = mixed(5, 100)
)
GAIC(m2, m3, k = 2)
GAIC(m2, m3, k = log(nrow(dt_trn)))

tt <- dt[, .(prop_if(smok_cess == 1)), keyby = .(age)]
tt[, pr := predictAll(m3, newdata = data.frame(age))$mu]
tt[, plot(age, V1)]
tt[, lines(age, pr)]


# Relapse method ----------------------------------------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- HSE_ts[wt_int > 0, .(
  wt_int, smok_quit_yrs, smok_status,
  agegrp, year, age, sex, qimd, ethnicity, sha)]
dt[, wt_int := wt_int*10000/sum(wt_int), by = year] # make pop of each hse =10000

dt[smok_quit_yrs < 20, sum(wt_int), keyby = smok_quit_yrs][, .(V1/head(V1, 1))][, plot(V1)]
tt <- dt[smok_quit_yrs < 20, sum(wt_int), keyby = .(smok_quit_yrs, sex, qimd)]
setkey(tt, smok_quit_yrs)
tt[, pr := V1/head(V1, 1), by = .(sex, qimd)]
tt[qimd == "2" & sex == "1", plot(smok_quit_yrs, pr, ylim = c(0, max(pr) * 1.1))]
tt[smok_quit_yrs > 0 & qimd == "2" & sex == "1",
   lines(1:19, predict(glm(pr~log(smok_quit_yrs))))]
tt[smok_quit_yrs > 0,
   pr := predict(glm(pr~log(smok_quit_yrs))), by = .(sex, qimd)]
tt[qimd == "1" & sex == "1", plot(smok_quit_yrs, pr, ylim = c(0, max(pr) * 1.1))]
tt[qimd == "1" & sex == "1", lines(smok_quit_yrs, pr, ylim = c(0, max(pr) * 1.1))]
tt[qimd == "5" & sex == "1", lines(smok_quit_yrs, pr, ylim = c(0, max(pr) * 1.1), col = "red")]
tt[qimd == "3" & sex == "1", lines(smok_quit_yrs, pr, ylim = c(0, max(pr) * 1.1), col = "blue")]
tt[qimd == "2" & sex == "1", lines(smok_quit_yrs, pr, ylim = c(0, max(pr) * 1.1), col = "green")]
tt[qimd == "4" & sex == "1", lines(smok_quit_yrs, pr, ylim = c(0, max(pr) * 1.1), col = "purple")]
tt[, pr := pr/shift(pr), by = .(sex, qimd)]
saveRDS(smok_relapse_model, "./secure_data/lifecourse_models/smok_relapse_model.rds")

# GAMLSS smok_cig_curr ----------------------------------------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- na.omit(HSE_ts[wt_int > 0 & smok_status == "4" & age > 15, .(
  smok_cig_curr, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, smok_cig_curr := round(smok_cig_curr/5)] # final result need to multiply by 5
dt[, age := scale(age, 50, 16.8)]
set.seed(47)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
diagnostics <- FALSE

marg_distr <- fitDistPred(dt_trn$smok_cig_curr, type = "counts",
                          weights = dt_trn$wt_int,
                          try.gamlss = TRUE, trace = TRUE, newdata = dt_crv$smok_cig_curr,
                          parallel = "multicore", ncpus = 12L)
marg_distr$fits
#  DEL ZASICHEL ZISICHEL   SICHEL


params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
# params$n <- 1e4
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
reldist_diagnostics(dt[, smok_cig_curr], y, dt[, wt_int/sum(wt_int)], y_wt,
                    main = expression(bold(Number~of~Cigarettes~per~day)), discrete = TRUE)

con1 <- gamlss.control(c.crit = 1e-2) # reduced for faster exploratory analysis.
setMKLthreads(1L)
distr_nam <- "DEL"

dt_trn[, .(smok_cig_curr_median = wtd.quantile(smok_cig_curr, weight = wt_int)), keyby = .(age, sex)
       ][, scatter.smooth(age, smok_cig_curr_median, col = sex)]


# GAMLSS smok_cig_ex ------------------------------------------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- na.omit(HSE_ts[wt_int > 0 & smok_status == "3" & age > 15, .(
  smok_cig_ex, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, smok_cig_ex := round(smok_cig_ex/5)] # final result need to multiply by 5
dt[, age := scale(age, 50, 16.8)]
set.seed(48)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
diagnostics <- FALSE

marg_distr <- fitDistPred(dt_trn$smok_cig_ex, type = "counts",
                          weights = dt_trn$wt_int,
                          try.gamlss = TRUE, trace = TRUE, newdata = dt_crv$smok_cig_ex,
                          parallel = "multicore", ncpus = 12L)
marg_distr$fits
#  ZASICHEL      DEL    ZAPIG


params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
# params$n <- 1e4
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
reldist_diagnostics(dt[, smok_cig_ex], y, dt[, wt_int/sum(wt_int)], y_wt,
                    main = expression(bold(Number~of~Cigarettes~per~day)), discrete = TRUE)

con1 <- gamlss.control(c.crit = 1e-2) # reduced for faster exploratory analysis.
setMKLthreads(1L)
distr_nam <- "ZASICHEL"

dt_trn[, .(smok_cig_curr_median = wtd.quantile(smok_cig_curr, weight = wt_int)), keyby = .(age, sex)
       ][, scatter.smooth(age, smok_cig_curr_median, col = sex)]

# GAMLSS smok_dur_ex ------------------------------------------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- na.omit(HSE_ts[wt_int > 0 & smok_status %in% c("2", "3") & age > 15, .(
  smok_dur_ex, smok_status, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, hist(smok_dur_ex, breaks = 100)]
dt[, smok_dur_ex := round(smok_dur_ex)]
dt[, age := scale(age, 50, 16.8)]
set.seed(49)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
diagnostics <- FALSE

marg_distr <- fitDistPred(dt_trn$smok_dur_ex, type = "counts",
                          weights = dt_trn$wt_int,
                          try.gamlss = TRUE, trace = TRUE, newdata = dt_crv$smok_dur_ex,
                          parallel = "multicore", ncpus = 12L)
marg_distr$fits
#   DPO    ZABNB    ZINBI


params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
# params$n <- 1e4
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
reldist_diagnostics(dt[, smok_dur_ex], y, dt[, wt_int/sum(wt_int)], y_wt,
                    main = expression(bold(Number~of~Years)), discrete = TRUE)

con1 <- gamlss.control(c.crit = 1e-2) # reduced for faster exploratory analysis.
setMKLthreads(1L)
distr_nam <- "DPO"

dt_trn[, .(smok_dur_ex_median = wtd.quantile(smok_dur_ex, weight = wt_int)), keyby = .(age, sex)
       ][, scatter.smooth(age, smok_dur_ex_median, col = sex)]
dt_trn[, .(smok_dur_ex_median = wtd.quantile(smok_dur_ex, weight = wt_int)), keyby = .(year)
       ][, scatter.smooth(year, smok_dur_ex_median)]

# GAMLSS smok_dur_curr ------------------------------------------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- na.omit(HSE_ts[wt_int > 0 & smok_status %in% c("4") & age > 15, .(
  smok_dur_curr = age - smok_init_age, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, hist(smok_dur_curr, breaks = 100)]
dt[, smok_dur_curr := round(smok_dur_curr)]
dt[, age := scale(age, 50, 16.8)]
set.seed(50)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
diagnostics <- FALSE

marg_distr <- fitDistPred(dt_trn$smok_dur_curr, type = "counts",
                          weights = dt_trn$wt_int,
                          try.gamlss = TRUE, trace = TRUE, newdata = dt_crv$smok_dur_curr,
                          parallel = "multicore", ncpus = 12L)
marg_distr$fits
#  DPO      NBF      NBI   SICHEL    ZINBI       SI


params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
# params$n <- 1e4
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
reldist_diagnostics(dt[, smok_dur_curr], y, dt[, wt_int/sum(wt_int)], y_wt,
                    main = expression(bold(Number~of~Years)), discrete = TRUE)

con1 <- gamlss.control(c.crit = 1e-2) # reduced for faster exploratory analysis.
setMKLthreads(1L)
distr_nam <- "DPO"

dt_trn[, .(smok_dur_curr_median = wtd.quantile(smok_dur_curr, weight = wt_int)), keyby = .(age, sex)
       ][, scatter.smooth(age, smok_dur_curr_median, col = sex)]
dt_trn[, .(smok_dur_curr_median = wtd.quantile(smok_dur_curr, weight = wt_int)), keyby = .(year)
       ][, scatter.smooth(year, smok_dur_curr_median)]

# GAMLSS smok_quit_yrs ------------------------------------------------------
HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- na.omit(HSE_ts[wt_int > 0 & smok_status %in% c("2", "3") & age > 15, .(
  smok_quit_yrs, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, hist(smok_quit_yrs, breaks = 100)]
dt[, smok_dur_curr := round(smok_dur_curr)]
dt[, age := scale(age, 50, 16.8)]
set.seed(50)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
diagnostics <- FALSE

marg_distr <- fitDistPred(dt_trn$smok_dur_curr, type = "counts",
                          weights = dt_trn$wt_int,
                          try.gamlss = TRUE, trace = TRUE, newdata = dt_crv$smok_dur_curr,
                          parallel = "multicore", ncpus = 12L)
marg_distr$fits
#  DPO      NBF      NBI   SICHEL    ZINBI       SI


params <- vector("list")
for (i in seq_along(marg_distr$parameters)) {
  nam <- marg_distr$parameters[[i]]
  params[[i]] <- get(nam, marg_distr)
  names(params)[i] <- nam
}
params$p <- seq(0.00001, 0.99999, 0.00001)
# params$n <- 1e4
distr_nam <- marg_distr$family[[1]]
y <- do.call(paste0("q", distr_nam), params)
y_wt <- rep(1/length(y), length(y))
reldist_diagnostics(dt[, smok_dur_curr], y, dt[, wt_int/sum(wt_int)], y_wt,
                    main = expression(bold(Number~of~Years)), discrete = TRUE)

con1 <- gamlss.control(c.crit = 1e-2) # reduced for faster exploratory analysis.
setMKLthreads(1L)
distr_nam <- "DPO"

dt_trn[, .(smok_dur_curr_median = wtd.quantile(smok_dur_curr, weight = wt_int)), keyby = .(age, sex)
       ][, scatter.smooth(age, smok_dur_curr_median, col = sex)]
dt_trn[, .(smok_dur_curr_median = wtd.quantile(smok_dur_curr, weight = wt_int)), keyby = .(year)
       ][, scatter.smooth(year, smok_dur_curr_median)]


# GAMLSS SBP model ---------------------
dt <- na.omit(HSE_ts[, .(sbp, year, age, agegrp, sex, qimd, ethnicity, sha, wt_nurse)])
marg_distr <- fitDist(sbp, data = dt, type = "realplus", weights = dt$wt_nurse,
                      try.gamlss = TRUE, trace = TRUE,
                      parallel = "multicore", ncpus = 12L)
# validation
plot(density(dt[, sbp], weights = dt[, wt_nurse/sum(wt_nurse)]))
lines(density(qGB2(seq(0.0001, 0.9999, 0.0001), marg_distr$mu, marg_distr$sigma,
     marg_distr$nu, marg_distr$tau)), col = "red")

set.seed(242)
lns <- sample(nrow(dt), nrow(dt) * 0.6)
dt_trn <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
dt_small <- dt_train[sample(.N, 10e3)]
con1 <- gamlss.control(c.crit = 1e-1) # reduced for exploratory analysis.

mod1_small <- gamlss(sbp~log(year) + age + sex + qimd,
                family = "GB2", weights = dt_small$wt_nurse,
                data = dt_small, method = mixed(20, 100))
mod1 <- gamlss(sbp~log(year) + age + sex + qimd, start.from = mod1_small,
               family = "GB2", weights = dt_train$wt_nurse,
               data = dt_train, method = mixed(5, 100), control = con1)
newpar <- par(mfrow = c(2,2), mar = par("mar") + c(0,1,0,0), col.axis = "blue4",
            col = "blue4",
            col.main = "blue4", col.lab = "blue4", pch = "+", cex = 0.45,
            cex.lab = 1.2, cex.axis = 1, cex.main = 1.2)
summary(mod1)
plot(mod1, par = newpar)
wp(mod1)
# mean should be ~0, var ~1, skew ~0, kurt = ~3
wp(mod1)
wp(mod1, dt_train$age, n.inter = 10)
Q.stats(mod1, xvar = dt_train$age, n.inter = 10)
zz <- validate_gamlss(dt_crv, mod1, 10, dt_train)
ggplot(zz, aes(sbp, colour = type)) +
  geom_density() +
  facet_wrap(.~agegrp, nrow = 3)

zz <- validate_gamlss(dt_crv, mod1, 1, dt_train)
Metrics::bias(zz[type == "Observed", sbp], zz[type == "Modelled", sbp])
Metrics::rmse(zz[type == "Observed", sbp], zz[type == "Modelled", sbp])
Metrics::mae(zz[type == "Observed", sbp], zz[type == "Modelled", sbp])

mod1_small <- update(mod1_small, sigma.fo = ~log(year) + age + sex + qimd)
dropterm(mod1_small, "log(year)", what = "sigma")
mod2 <- update(mod1, sigma.fo = ~log(year) + age + sex + qimd, start.from = mod1_small)
mod2 <- refit(mod2)
zz <- validate_gamlss(dt_crv, mod1, 10, dt_train)ggplot(zz, aes(sbp, colour = type)) +
  geom_density() +
  facet_wrap(.~agegrp, nrow = 2)
# Validation plot ---------------------------------------------



zz <- validate_gamlss(dt, mod1, 10)
ggplot(zz, aes(sbp, colour = type)) +
  geom_density() +
  facet_wrap(.~agegrp, nrow = 2)

gg <- ggplot(salt3, aes(salt, colour = type)) +
  geom_density() +
  facet_wrap(time~sex, ncol = 2) +
  ylim(0, 0.18) +
  scale_x_continuous("Salt (g/d)", breaks = c(0, 10, 20, 30), limits = c(0, 30)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(
          margin = margin(r = 24, unit = "pt")))
cowplot::ggsave("./Exposure/itsgamlss_lin_validation.png", gg, dpi = 600, width = 8, height = 6,
                units = "cm", scale = 3)

gg <- ggplot(salt4, aes(salt, colour = type)) +
  geom_density() +
  facet_wrap(time~sex, ncol = 2) +
  ylim(0, 0.18)+
  scale_x_continuous("Salt (g/d)", breaks = c(0, 10, 20, 30), limits = c(0, 30)) +
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(
          margin = margin(r = 24, unit = "pt")))
cowplot::ggsave("./Exposure/itsgamlss_log_validation.png", gg, dpi = 600, width = 8, height = 6,
                units = "cm", scale = 3)



mod1 <- gamlss(sbp~log(year) + pb(age) + sex + qimd,
                ~log(year) + pb(age) + sex + qimd,
                ~log(year) + pb(age) + sex + qimd,
                ~log(year) + pb(age) + sex + qimd,
                family = "GB2", weights = dt$wt_nurse,
                data = dt, method = mixed(5, 50))
distr_nam_tr_short <- c("GB2", "LOGNO") # selected based on marg_distr
t1 <- chooseDist(mod1, type = "extra", extra = distr_nam_tr_short) #, parallel = "multicore", ncpus = 12
getOrder(t1,3) # based on GAIC & BIC RGtr fits best

# Impute PA --------------------------
HSE.ts[, id := 1:.N]
HSE.ts[, year_ch := factor(year)]
A = copy(HSE.ts[is.na(a30to06m) == T & age >15 & is.na(agegroup)==F & is.na(sex)==F & is.na(qimd)==F, .(id, agegroup, sex, qimd, year_ch, wt_int, psu, cluster)]) # missing pa to be imputed

B = copy(HSE.ts[age >15 & is.na(a30to06m)==F & is.na(agegroup)==F & is.na(sex)==F & is.na(qimd)==F, .(id, a30to06m, agegroup, sex, year_ch, qimd, wt_int, psu, cluster)]) #

B[, a30to06m := ordered(a30to06m)]
A[, agegroup := ordered(agegroup)]
B[, agegroup := ordered(agegroup)]

A.srv <- svydesign(id=~psu, strata=~cluster, weights=~wt_int, nest=F, data=A, check.strata = T)
B.srv <- svydesign(id=~psu, strata=~cluster, weights=~wt_int, nest=F, data=B, check.strata = T)

out.hz <- harmonize.x(svy.A  = A.srv,
                      svy.B  = B.srv,
                      x.tot  = NULL,
                      form.x =~ year_ch + agegroup + sex + qimd, cal.method = "poststratify")
A.srv = out.hz$cal.A
B.srv = out.hz$cal.B

tt <- comb.samples(out.hz$cal.A,
                   out.hz$cal.B,
                   NULL,
                   "id", "a30to06m",
                   form.x = ~year_ch + agegroup + sex + qimd,
                   micro = T)

ttt <- cbind(A, tt$Z.A)
setnames(ttt, paste0("a30to06m", 1:8), paste0(0:7))

ttt <- melt(ttt, id=c("id", "agegroup", "sex", "qimd", "year_ch", "wt_int", "psu", "cluster"), variable.name = "a30to06m.imp")
setkey(ttt, id, value)
ttt[value < 0, value := 0 ]
ttt[value > 1, value := 1 ]
ttt <- ttt[, .SD[sample(.N, 1, prob = value)], by = id]
HSE.ts = copy(merge(HSE.ts, ttt[,.(id, a30to06m.imp)], by = "id", all.x = T))
HSE.ts[, a30to06m.imp := as.integer(as.character(a30to06m.imp))]
HSE.ts[is.na(a30to06m.imp), a30to06m.imp := a30to06m]

# Impute porftvg (only consider 2011 to impute in 2012)
A = copy(HSE.ts[is.na(porftvg)==T & age >15 & is.na(agegroup)==F & is.na(sex)==F & is.na(qimd)==F & year == 1, .(id, agegroup, sex, qimd, wt_int, psu, cluster)])

B = copy(HSE.ts[age >15 & is.na(porftvg)==F & is.na(agegroup)==F & is.na(sex)==F & is.na(qimd)==F & year == 0, .(id, porftvg, agegroup, sex, qimd, wt_int, psu, cluster)]) #

B[, porftvg := ordered(porftvg)]
A[, agegroup := ordered(agegroup)]
B[, agegroup := ordered(agegroup)]

A.srv <- svydesign(id=~psu, strata=~cluster, weights=~wt_int, nest=F, data=A, check.strata = T)
B.srv <- svydesign(id=~psu, strata=~cluster, weights=~wt_int, nest=F, data=B, check.strata = T)

out.hz <- harmonize.x(svy.A  = A.srv,
                      svy.B  = B.srv,
                      x.tot  = NULL,
                      form.x = ~-1+(agegroup+sex+qimd)^2,
                      cal.method = "linear")

tt <- comb.samples(out.hz$cal.A,
                   out.hz$cal.B,
                   NULL,
                   "id", "porftvg",
                   form.x=~-1+(agegroup+sex+qimd)^2,
                   micro = T)

ttt <- cbind(A, tt$Z.A)
setnames(ttt, paste0("porftvg", 1:9), paste0(0:8))

ttt <- melt(ttt, id=c("id", "agegroup", "sex", "qimd", "wt_int", "psu", "cluster"), variable.name = "porftvg.imp")
setkey(ttt, id, value)
ttt[value < 0, value := 0 ]
ttt[value > 1, value := 1 ]
ttt <- ttt[, .SD[sample(.N, 1, prob = value)], by = id]
HSE.ts = copy(merge(HSE.ts, ttt[,.(id, porftvg.imp)], by = "id", all.x = T))
HSE.ts[, porftvg.imp := as.integer(as.character(porftvg.imp))]
HSE.ts[is.na(porftvg.imp), porftvg.imp := porftvg]

#HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts[wt_int>0,], check.strata = T)
#HSE.ts.srv.nurse <- svydesign(id=~psu, strata =~cluster, weights = ~wt_nurse, nest=F, data=HSE.ts[wt_nurse>0,], check.strata = T)
#HSE.ts.srv.blood <- svydesign(id=~psu, strata =~cluster, weights = ~wt_blood, nest=F, data=HSE.ts[wt_blood>0,], check.strata = T)
#save(HSE.ts, file="./Lagtimes/HSE.ts.RData")
#load(file="./Models/IMPACTncd/Lagtimes/HSE.ts.RData")

# Build Models

# Salt prepare file  -------------------------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
#load(file="./Models/IMPACTncd Liverpool/Lagtimes/HSE.ts.RData")

# Salt
loadcmp(file="./salt.Rc")
#source(file="./Models/IMPACTncd/salt.R")
HSE.ts[bmival<16 & age>19, bmival := 16]
HSE.ts[bmival>50 & age>19, bmival := 50]
HSE.ts[is.na(wt_urine), wt_urine := 0]
HSE.ts[age>85, age:= 85]
# HSE.ts[age>15 & sex== "1", Na24 := Na24.men(.N, sodium, creatin, potass, bmival, age)]
# HSE.ts[age>15 & sex== "2", Na24 := Na24.women(.N, sodium, creatin, potass, bmival, age)]
# HSE.ts[, salt := Na24 * 58.5/1000]
# HSE.ts[, summary(salt)]

# ignore survey effects and match with 24h known distributions
HSE = copy(HSE.ts[between(age, 19, 64) & !is.na(qimd) & !is.na(sodium) & !is.na(bmival), ])
HSE <- rbindlist(sample(list(HSE), 500, T),idcol = T)
HSE[age>15 & sex== "1", Na24 := Na24.men(.N, sodium, creatin, potass, bmival, age)]
HSE[age>15 & sex== "2", Na24 := Na24.women(.N, sodium, creatin, potass, bmival, age)]
HSE[, salt := Na24 * 58.5/1000]
# year 2001. I will use HSE2003 sample but 24h distribution from 2001. Year will be -10 (and not -8)
#HSE[year == -8, year := -10]
tmp1 = copy(HSE[year == -8,])

tmp1[between(age, 19, 24) & sex == "1", grp := 1]
tmp1[between(age, 19, 24) & sex == "2", grp := 2]
tmp1[between(age, 25, 34) & sex == "1", grp := 3]
tmp1[between(age, 25, 34) & sex == "2", grp := 4]
tmp1[between(age, 35, 49) & sex == "1", grp := 5]
tmp1[between(age, 35, 49) & sex == "2", grp := 6]
tmp1[between(age, 50, 64) & sex == "1", grp := 7]
tmp1[between(age, 50, 64) & sex == "2", grp := 8]

tmp1[grp > 0, salt24h.rank := (frank(salt, na.last = F, ties.method="random")-1)/(.N - 1), by = grp]
setkey(tmp1, grp, salt24h.rank)

# define distributions from 24h urine from 2001 (table 4.1, 4.2)
# Men 19-24
p <- c(0,	0,	0.09,	0.33,	0.46,	0.6,	0.75,	0.96,	1,	0.02,	0.37,	0.81,	1,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 14.0, 15.8, 3, 6, 9, 15, 18, 10.6, 6, 16.6) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p == 1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=7, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt1 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 19-24
p <- c(0.04,	0.17,	0.34,	0.65,	0.7,	0.84,	0.88,	0.9,	0.66,	0.84,	0.92,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 14.0, 15.8, 9, 12, 18, 7.6, 1.7, 23.2) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=25, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.lnorm.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt2 <- get.lnorm.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


# Men 25-34
p <- c(0.06,	0.12,	0.27,	0.33,	0.48,	0.59,	0.7,	0.78,	0.05,	0.2,	0.34,	0.57,	0.73,	0.89,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 14.0, 15.8, 3, 6, 9, 12, 15, 18, 10.9, 2.2, 22.3) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=5, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt3 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 25-34
p <- c(0.05,	0.25,	0.41,	0.57,	0.74,	0.85,	0.91,	0.95,	0.06,	0.29,	0.59,	0.81,	0.92,	0.97,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 14.0, 15.8, 3, 6, 9, 12, 15, 18, 8, 1.9, 22.2) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=7, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.chisq.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt4 <- get.chisq.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Men 35-49
p <- c(0.02,	0.09,	0.19,	0.36,	0.51,	0.6,	0.74,	0.82,	0.02,	0.13,	0.39,	0.58,	0.8,	0.91,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 14.0, 15.8, 3, 6, 9, 12, 15, 18, 10.2, 2.4, 22.1) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=6, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt5 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 35-49
p <- c(0.05,	0.2,	0.43,	0.67,	0.8,	0.87,	0.92,	0.97,	0.05,	0.31,	0.68,	0.85,	0.96,	1,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 14.0, 15.8, 3, 6, 9, 12, 15, 18, 7.6, 2.6, 16.2) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=6, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt6 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Men 50-64 (50 - 100)
p <- c(0.08,	0.15,	0.27,	0.42,	0.54,	0.67,	0.77,	0.85,	0.05,	0.18,	0.42,	0.65,	0.83,	0.91,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 14.0, 15.8, 3, 6, 9, 12, 15, 18, 10.1, 2.1, 21.2) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=5, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt7 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 50-64 (50 - 100)
p <- c(0.12,	0.28,	0.5,	0.68,	0.84,	0.92,	0.94,	0.98,	0.07,	0.38,	0.69,	0.91,	0.96,	0.99,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 14.0, 15.8, 3, 6, 9, 12, 15, 18, 7, 2.3, 15.7) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=8, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt8 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.001, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


tmp <- data.table(
  salt24h = c( rtriang(tmp1[grp == 1, .N], tt1[1], tt1[2], tt1[3]),
               rlnorm (tmp1[grp == 2, .N], tt2[1], tt2[2]),
               rtriang(tmp1[grp == 3, .N], tt3[1], tt3[2], tt3[3]),
               rchisq (tmp1[grp == 4, .N], tt4[1]),
               rtriang(tmp1[grp == 5, .N], tt5[1], tt5[2], tt5[3]),
               rgamma (tmp1[grp == 6, .N], tt6[1]  , tt6[2]),
               rtriang(tmp1[grp == 7, .N], tt7[1], tt7[2], tt7[3]),
               rgamma (tmp1[grp == 8, .N], tt8[1]  , tt8[2])
  ),
  grp = c(rep(1, tmp1[grp == 1, .N]),
          rep(2, tmp1[grp == 2, .N]),
          rep(3, tmp1[grp == 3, .N]),
          rep(4, tmp1[grp == 4, .N]),
          rep(5, tmp1[grp == 5, .N]),
          rep(6, tmp1[grp == 6, .N]),
          rep(7, tmp1[grp == 7, .N]),
          rep(8, tmp1[grp == 8, .N])
  )
)

tmp[, salt24h.rank := (frank(salt24h, na.last = F, ties.method="random")-1)/(.N - 1), by = grp]
setkey(tmp, grp, salt24h.rank)
tmp1 <-   tmp[tmp1, roll = "nearest"]
tmp1[, `:=` (grp = NULL, salt24h.rank = NULL)]

# year 2006. for 24h urine participants were selected from HSE2005
tmp2 = copy(HSE[year == -5,])
tmp2[, year := -5]

tmp2[between(age, 19, 24) & sex == "1", grp := 1]
tmp2[between(age, 19, 24) & sex == "2", grp := 2]
tmp2[between(age, 25, 34) & sex == "1", grp := 3]
tmp2[between(age, 25, 34) & sex == "2", grp := 4]
tmp2[between(age, 35, 49) & sex == "1", grp := 5]
tmp2[between(age, 35, 49) & sex == "2", grp := 6]
tmp2[between(age, 50, 64) & sex == "1", grp := 7]
tmp2[between(age, 50, 64) & sex == "2", grp := 8]

tmp2[grp > 0, salt24h.rank := (frank(salt, na.last = F, ties.method="random")-1)/(.N - 1), by = grp]
setkey(tmp2, grp, salt24h.rank)

# define distributions from 24h urine from 2006
# Men 19-24
p <- c(0,	0,	0,	0.35,	0.5, 0.66,	0.94,	0,	0,	0.31,	0.59,	0.94,	0.5,	0.025)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 15, 18, 10.4, 7.5) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=10, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt1 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 19-24
p <- c(0,	0.15,	0.62,	1,	1,	1,	1,	0,	0.15,	1,	1,	1,	0.5,	0.025)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 15, 18, 6.6, 4.6) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=.95, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.norm.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt2 <- get.norm.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


# Men 25-34
p <- c(0,	0,	0.24,	0.44,	0.63,	0.73,	0.95,	0,	0.15,	0.49,	0.73,	0.95,	0.95,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 9.8, 5.2, 23.6) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=11, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.lnorm.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt3 <- get.lnorm.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 25-34
p <- c(0.02,	0.1,	0.31,	0.51,	0.77,	0.86,	0.98,	0.02,	0.16,	0.56,	0.84,	0.98,	1,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 8.7, 2.5, 16.8) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=6, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt4 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Men 35-49
p <- c(0,	0.04,	0.28,	0.47,	0.61,	0.68,	0.94,	0,	0.12,	0.47,	0.66,	0.92,	0.97,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 10.2, 2.4, 22.1) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=13, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt5 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 35-49
p <- c(0.02,	0.18,	0.41,	0.66,	0.82,	0.95,	0.98,	0.01,	0.28,	0.69,	0.95,	0.98,	1,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 7.6, 2.6, 16.2) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=6.25, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt6 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Men 50-64 (50 - 100)
p <- c(0.02,	0.07,	0.21,	0.39,	0.59,	0.72,	0.92,	0.02,	0.11,	0.41,	0.71,	0.89,	0.96,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 9.7, 2.7, 20.8) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=7.3, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt7 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 50-64 (50 - 100)
p <- c(0.1,	0.27,	0.59,	0.78,	0.91,	0.96,	0.99,	0.06,	0.43,	0.82,	0.95,	0.99,	0.1,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 6.4, 2, 14.4) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=11, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt8 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.001, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


tmp <- data.table(
  salt24h = c( rtriang(tmp2[grp == 1, .N], tt1[1], tt1[2], tt1[3]),
               rnorm  (tmp2[grp == 2, .N], tt2[1], tt2[2]),
               rlnorm (tmp2[grp == 3, .N], tt3[1], tt3[2]),
               rgamma (tmp2[grp == 4, .N], tt4[1], tt4[2]),
               rgamma (tmp2[grp == 5, .N], tt5[1], tt5[2]),
               rgamma (tmp2[grp == 6, .N], tt6[1], tt6[2]),
               rgamma (tmp2[grp == 7, .N], tt7[1], tt7[2]),
               rgamma (tmp2[grp == 8, .N], tt8[1], tt8[2])
  ),
  grp = c(rep(1, tmp2[grp == 1, .N]),
          rep(2, tmp2[grp == 2, .N]),
          rep(3, tmp2[grp == 3, .N]),
          rep(4, tmp2[grp == 4, .N]),
          rep(5, tmp2[grp == 5, .N]),
          rep(6, tmp2[grp == 6, .N]),
          rep(7, tmp2[grp == 7, .N]),
          rep(8, tmp2[grp == 8, .N])
  )
)

tmp[, salt24h.rank := (frank(salt24h, na.last = F, ties.method="random")-1)/(.N - 1), by = grp]
setkey(tmp, grp, salt24h.rank)
tmp2 <-   tmp[tmp2, roll = "nearest"]
tmp2[, `:=` (grp = NULL, salt24h.rank = NULL)]

# Year 2008. HSE 2007 will be used. Original 24h sample from UK (not only England)
tmp3 = copy(HSE[year == -2, ])
tmp3[, year := -3]

tmp3[between(age, 19, 24) & sex == "1", grp := 1]
tmp3[between(age, 19, 24) & sex == "2", grp := 2]
tmp3[between(age, 25, 34) & sex == "1", grp := 3]
tmp3[between(age, 25, 34) & sex == "2", grp := 4]
tmp3[between(age, 35, 49) & sex == "1", grp := 5]
tmp3[between(age, 35, 49) & sex == "2", grp := 6]
tmp3[between(age, 50, 64) & sex == "1", grp := 7]
tmp3[between(age, 50, 64) & sex == "2", grp := 8]

tmp3[grp > 0, salt24h.rank := (frank(salt, na.last = F, ties.method="random")-1)/(.N - 1), by = grp]
setkey(tmp3, grp, salt24h.rank)

# define distributions from 24h urine from 2006
# Men 19-24
p <- c(0,	0.14,	0.14,	0.25,	0.34,	0.57,	1,	0,	0.14,	0.25,	0.46,	1,	1,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 12.1, 3.7, 14.4) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=8, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt1 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 19-24
p <- c(0.07,	0.07,	0.39,	0.49,	0.49,	0.72,	0.91,	0,	0.21,	0.49,	0.72,	0.91,	0.91,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 11.3, 3.3, 18) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=15, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt2 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


# Men 25-34
p <- c(0.07,	0.07,	0.19,	0.37,	0.6,	0.72,	0.97,	0.05,	0.13,	0.43,	0.7,	0.87,	0.97,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 9.1, 2.5, 22.1) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=13.1, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt3 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 25-34
p <- c(0.01,	0.19,	0.49,	0.68,	0.79,	0.83,	0.98,	0,	0.28,	0.68,	0.82,	0.98,	0.99,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 7.1, 3.9, 14.8) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=9.5, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt4 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Men 35-49
p <- c(0.01,	0.1,	0.28,	0.43,	0.65,	0.78,	0.96,	0.01,	0.19,	0.47,	0.76,	0.95,	0.97,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 9.4, 4, 19.1) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=7, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt5 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 35-49
p <- c(0.04,	0.24,	0.52,	0.72,	0.86,	0.92,	0.99,	0.02,	0.38,	0.76,	0.92,	0.99,	0.99,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 7, 3, 14.3) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=6.25, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt6 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Men 50-64 (50 - 100)
p <- c(0.04,	0.12,	0.26,	0.46,	0.64,	0.8,	0.95,	0.02,	0.2,	0.48,	0.78,	0.93,	0.99,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 9.5, 3.2, 16.6) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=5, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt7 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 50-64 (50 - 100)
p <- c(0.1,	0.31,	0.55,	0.78,	0.88,	0.95,	0.99,	0.06,	0.39,	0.8,	0.95,	0.98,	1,	0.5,	0.025,	0.975)  # Known percentiles probabilities
e <- c(3.5, 5.3, 7.0, 8.8, 10.5, 12.3, 15.8, 3, 6, 9, 12, 15, 18, 6.7, 2.3, 13.9) # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=3.8, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt8 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.001, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


tmp <- data.table(
  salt24h = c( rtriang(tmp3[grp == 1, .N], tt1[1], tt1[2], tt1[3]),
               rtriang(tmp3[grp == 2, .N], tt2[1], tt2[2], tt2[3]),
               rgamma (tmp3[grp == 3, .N], tt3[1], tt3[2]),
               rtriang(tmp3[grp == 4, .N], tt4[1], tt4[2], tt4[3]),
               rgamma (tmp3[grp == 5, .N], tt5[1], tt5[2]),
               rgamma (tmp3[grp == 6, .N], tt6[1], tt6[2]),
               rtriang(tmp3[grp == 7, .N], tt7[1], tt7[2], tt7[3]),
               rgamma (tmp3[grp == 8, .N], tt8[1], tt8[2])
  ),
  grp = c(rep(1, tmp3[grp == 1, .N]),
          rep(2, tmp3[grp == 2, .N]),
          rep(3, tmp3[grp == 3, .N]),
          rep(4, tmp3[grp == 4, .N]),
          rep(5, tmp3[grp == 5, .N]),
          rep(6, tmp3[grp == 6, .N]),
          rep(7, tmp3[grp == 7, .N]),
          rep(8, tmp3[grp == 8, .N])
  )
)

tmp[, salt24h.rank := (frank(salt24h, na.last = F, ties.method="random")-1)/(.N - 1), by = grp]
setkey(tmp, grp, salt24h.rank)
tmp3 <- tmp[tmp3, roll = "nearest"]
tmp3[, `:=` (grp = NULL, salt24h.rank = NULL)]

# Year 2011. HSE 2010 & 12 wll be used. Note 10 has no sec gradiaent. an exception
tmp4 = copy(HSE[year == 1, ])
tmp4[, year := 0]

tmp4[between(age, 19, 34) & sex == "1", grp := 1]
tmp4[between(age, 19, 34) & sex == "2", grp := 2]
tmp4[between(age, 35, 49) & sex == "1", grp := 3]
tmp4[between(age, 35, 49) & sex == "2", grp := 4]
tmp4[between(age, 50, 64) & sex == "1", grp := 5]
tmp4[between(age, 50, 64) & sex == "2", grp := 6]

tmp4[grp > 0, salt24h.rank := (frank(salt, na.last = F, ties.method="random")-1)/(.N - 1), by = grp]
setkey(tmp4, grp, salt24h.rank)

# Men 19-34 (16-34)
p <- c(0.025, 0.15, 0.3, 0.44, 0.5, 0.61, 0.76, 0.95, 0.975, 0, 0.21,  0.72,  0.98)  # Known percentiles probabilities
e <- c(4.2  , 5.3 , 7  , 8.8 , 9.3, 10.5, 12.3, 15.8, 17.3 , 3, 6   ,  12  ,   18)  # salt from Sodium Survey England 2011 tables 9 and 10
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=4, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt1 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


# Women 19-34 (16-34)
p <- c(0.025, 0.12, 0.23, 0.5, 0.52, 0.77, 0.91, 0.94, 0.975, 0.05, 0.32, 0.81,  0.97, 1)
e <- c(2.8  , 3.5 , 5.3 , 7  , 7   , 8.8 , 10.5, 12.3, 15.2 , 3   , 6   , 9   ,  15  , 18)
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=7.8, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.lnorm.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt2 <- get.lnorm.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Men 35-49
p <- c(0.025, 0.06, 0.26, 0.39, 0.5, 0.6 , 0.77, 0.94, 0.975, 0.01, 0.11, 0.42, 0.89, 0.96)
e <- c(4.3  , 5.3 , 7   , 8.8 , 9.7, 10.5, 12.3, 15.8, 18.8 , 3   , 6   , 9   , 15  , 18)
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=5, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt3 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


# Women 35-49
p <- c(0.025, 0.3, 0.5, 0.66, 0.78, 0.9 , 0.975, 0.98, 1, 0.02, 0.48, 0.80, 0.94,  1)
e <- c(3.4  , 5.3, 6.1, 7   , 8.8 , 10.5, 12.1 , 12.3, 15.8, 3   , 6   , 9   , 12  , 18)
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=6, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.triang.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt4 <- get.triang.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


# Men 50-64 (50-80)
p <- c(0.025, 0.21, 0.38, 0.5, 0.62, 0.81, 0.9 , 0.97, 0.975, 0, 0.28, 0.64, 0.88, 0.96)
e <- c(3.1  , 5.3 , 7   , 7.8, 8.8 , 10.5, 12.3, 15.8, 18   , 3   , 6   , 9   , 12  , 15)
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=4, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.lnorm.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt5 <- get.lnorm.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))

# Women 50-64 (50-80)
p <- c(0.025, 0.09, 0.33, 0.5, 0.6, 0.82, 0.95, 0.96, 0.975, 0.05, 0.46, 0.85, 1   , 1)
e <- c(2.6  , 3.5 , 5.3 , 6.3, 7  , 8.8 , 10.5, 12.3, 12.7 , 3   , 6   , 9   , 15  , 18)
xx <- which(p == 0 | p ==1) # index for 0, 1
p <- sort(p[p != 0 & p != 1]) # remove 0, 1
ifelse (length(xx)==0, e <- sort(e), e <- sort(e[-xx])) # remove relevant values from e
#fit.perc(p, e, show.output=F, tolPlot=3, tolConv=0.01, fit.weights=1/(abs(0.5-p)+1))
#get.gamma.par(p, e, show.output=T, plot=T, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))
tt6 <- get.gamma.par(p, e, show.output=F, plot=F, tol=0.01, fit.weights=1/(abs(0.5-p)+1), scaleX=c(0.025,0.975))


tmp <- data.table(
  salt24h = c( rtriang(tmp4[grp == 1, .N], tt1[1], tt1[2], tt1[3]),
               rlnorm (tmp4[grp == 2, .N], tt2[1], tt2[2]),
               rgamma (tmp4[grp == 3, .N], tt3[1], tt3[2]),
               rtriang(tmp4[grp == 4, .N], tt4[1], tt4[2], tt4[3]),
               rlnorm (tmp4[grp == 5, .N], tt5[1], tt5[2]),
               rgamma (tmp4[grp == 6, .N], tt6[1], tt6[2])
  ),
  grp = c(rep(1, tmp4[grp == 1, .N]),
          rep(2, tmp4[grp == 2, .N]),
          rep(3, tmp4[grp == 3, .N]),
          rep(4, tmp4[grp == 4, .N]),
          rep(5, tmp4[grp == 5, .N]),
          rep(6, tmp4[grp == 6, .N])
  )
) # parameters from fit distr to 24h.R

tmp[, salt24h.rank := (frank(salt24h, na.last = F, ties.method="random")-1)/(.N - 1), by = grp]
setkey(tmp, grp, salt24h.rank)
tmp4 <-   tmp[tmp4, roll = "nearest"]
tmp4[, `:=` (grp = NULL, salt24h.rank = NULL)]

# combine all
HSE.salt <- setkey(rbind(tmp1, tmp2, tmp3, tmp4, use.names = T), age, sex, qimd, year)
rm(tmp1, tmp2, tmp3, tmp4, tmp, HSE, HSE.ts)

# explore salt correlations
tt <- HSE.salt[, mean(salt24h, na.rm=T), by = year]
#tt[, V2 := (V1-6)/(max(V1)-6)] # scale from 0 t0 1 with threshold 6
#plot((V2*tt[, max(V1)-6]+6)~year, tt, xlim=c(-10, 60), ylim=c(0, 10))
plot(tt, xlim=c(-10, 25), ylim=c(0, 12))
abline(h = 6)
lines(y=(predict(glm(I(V1^(1/3))~I(log(year+14)), family=gaussian(link="identity"), data=tt), data.frame(year=-10:60), type="response"))^(3), x=-10:60, col="red")

lines(y=predict(glm(V1~I(log(year+14)^-2)+I((year+14)^-2), family=gaussian(link="identity"), data=tt), data.frame(year=-10:60), type="response"), x=-10:60, col="green")
lines(y=6 + tt[, max(V1)-6]*(predict(glm(V2~year, family=binomial(logit), data=tt), data.frame(year=-10:600), type="response")), x=-10:600, col="purple")

lines(y=predict(glm(V1~I(log(year^(1/year))), family=gaussian, data=tt), data.frame(year=-10:60), type="response"), x=-10:60, col="green")
lines(y=predict(glm(V1~poly(1/(year+30), 2), family=gaussian(), data=tt), data.frame(year=-10:60), type="response"), x=-10:60, col="black")
lines(y=predict(glm(V1~poly(year, 2), family=gaussian(), data=tt), data.frame(year=-10:60), type="response"), x=-10:60, col="red")

lines(y=-6 + predict(glm(I((V1+6)/20)~year, family=binomial(logit), data=tt), data.frame(year=-10:60), type="response")*20, x=-10:60, col="blue")

tt <- HSE.salt[, mean(salt24h, na.rm=T), by = age]
scatter.smooth(tt, xlim=c(19, 100), ylim=c(6, 10))
lines(y=predict(glm(V1~age, family=gaussian(link="identity"), data=tt), data.frame(age=19:100), type="response"), x=19:100, col="green")
lines(y=predict(glm(V1~poly(log(age), 3), family=gaussian(link="identity"), data=tt), data.frame(age=19:100), type="response"), x=19:100, col="red")

tt <- HSE.salt[, mean(salt24h, na.rm=T), by = qimd]
tt[, qimd := as.numeric(as.character(qimd))]
plot(tt, xlim=c(0, 6), ylim=c(6, 10))
lines(y=predict(glm(V1~poly(qimd, 2), family=gaussian(link="identity"), data=tt), data.frame(qimd=1:5), type="response"), x=1:5, col="green")

tt <- HSE.salt[, mean(salt24h, na.rm=T), by = round(bmival)]
tt[, round := as.integer(as.character(round))]
plot(tt, xlim=c(20, 50), ylim=c(0, 12))
lines(y=predict(glm(V1~poly(round, 2), family=gaussian(link="identity"), data=tt), data.frame(round=20:50), type="response"), x=20:50, col="green")

# Salt fitting ------------------------------------------------------------
HSE.salt <- HSE.salt[, .(age, sex, qimd, year, .id, id, salt24h, wt_urine)]
HSE.salt[, hist(salt24h^(1/3))]
q1 <- glm(I(salt24h^(1/3))~ (I(log(year+9)) + poly(log(age), 3) + sex + qimd)^3, weights = wt_urine, family = gaussian(), data = HSE.salt[.id < 20])
summary(q1)
q11 <- stepAIC(q1) #AIC: 1,634,000
summary(q11)
rm(q1, q11)


q2 <- glm(I(salt24h^(-1/3))~ I(log(year+14)^-2)+I((year+14)^-2) + log(age) + sex + qimd + qimd:I(log(year+14)^-2), weights = wt_urine, family = gaussian(), data = HSE.salt)
AIC(q2)
rm(q2)

q3 <- glm(I(salt24h^(1/3))~ poly(year, 2) + log(age) + sex + qimd + qimd:poly(year, 2), weights = wt_urine, family = gaussian(), data = HSE.salt)
AIC(q3)
rm(q3)

summary(q12)
View(data.frame(list(a=coef(q11), b=coef(q12))))

q2 <- rq(I(salt24h^(1/3))~ log(year+14),
         tau = .5,
         data = HSE.salt,
         weights = wt_urine,
         method = "fn")
pr <- (predict(q2, data.frame(year= -10:20)))^3
tt <- HSE.salt[, quantileWt(salt24h, wt_urine,  0.5, na.rm=T), by = year]
plot(tt, xlim=c(-10, 20), ylim=c(0, 20))
abline(v=10, h=6)
lines(y=pr, x = -10:20)

salt.rq <- rq(I(salt24h^(1/3)) ~ (I(log(year + 14)) + poly(log(age), 3) + sex + qimd)^2,
              tau = c(0.01, 1:19/20, 0.99),
              data = HSE.salt[.id == 1, ], # Then coeff to be injected from all versions
              weights = wt_urine,
              method = "br")

# produce a model for each version of population. NOTE year+16 used. steeper decline than +14
# but works for init.year == 2006 to avoid log(0) (-5 -10(cancer.lag))
if (Sys.info()["sysname"] == "Linux") {
  registerDoParallel(6L) # used for forking. only linux
}
if (Sys.info()["sysname"] == "Windows") {
  cl <- makeCluster(6) # used for clustering. win compatible
  registerDoParallel(cl)
}

salt.rq.coef <- vector("list", 500)
salt.rq.coef <-
  foreach(jj = 1 : 500,
          .inorder = F,
          .verbose = T,
          .packages = c("data.table",
                        "dplyr",
                        "quantreg")
  ) %dopar% {
    md1 <- rq(I(salt24h^(1/3)) ~ (I(log(year + 14)) + poly(log(age), 3) + sex + qimd)^2,
              tau = c(0.01, 1:19/20, 0.99),
              data = HSE.salt[.id == jj, ],
              weights = wt_urine,
              method = "br")
    coef(md1)
  }
if (exists("cl")) stopCluster(cl)

save(salt.rq, file="./Lagtimes/salt.rq.rda")
save(salt.rq.coef, file="./Lagtimes/salt.rq.coef.rda")

#save(salt.rq.coef, file="./Lagtimes/salt.rq.rda")
#save(salt.rq.coef, file="./Models/IMPACTncd/Lagtimes/salt.rq.coef.rda")
apply(simplify2array(salt.rq.coef), 1:2, mean) # calculate mean per element
apply(simplify2array(salt.rq.coef), 1:2, quantile, c(0.025, 0.5, 0.975)) # 95%CI
View(apply(simplify2array(salt.rq.coef), 1:2, quantile, c(0.025, 0.5, 0.975)))
#salt.rq$coefficients <- sample(salt.rq.coef,1)[[1]] for the MC simulation
#salt.rq$coefficients <- apply(simplify2array(salt.rq.coef), 1:2, mean)
q2 <- salt.rq
# When I use 100 percentiles, give some wrong values. eg pct .4 gives salt 5 and
# pct .41 gives salt 4.9 (should be above 5). The -tiles I used seems the best option
salt.rq.old <- rq(salt24h ~ log(year+14) + log(age) + sex + qimd + log(age):sex,
                  tau = c(0.01, 1:19/20, 0.99),
                  data = HSE.salt,
                  weights = wt_urine,
                  method = "fn")
# log(year + 14). 14 was selected so median salt for an 85 yo woman of qimd=1 is
# above 3gr in 2111. See below. Other wise linear decrease is too steep and produces
# unrealistic results.
#
# Also the predictors selected seem stable for each iteration of HSE.salt.
# interaction terms were very unstable

# needed for plots
q2$coefficients <- apply(simplify2array(salt.rq.coef), 1:2, mean)[, 11]
pred.s <- cmpfun(function(year, age, sex, qimd) {
  if (is.factor(sex)==F) {
    sex <-  factor(sex,
                   levels = c(1,2),
                   ordered = F)
  }
  if (is.ordered(qimd)==F) {
    qimd <- factor(qimd,
                   levels = c(1,2,3,4,5),
                   ordered = T)
  }
  cc <- predict(q2,
                data.frame(year= year,
                           age = age,
                           sex = sex,
                           qimd = qimd
                )
  )

  return(cc^3)
}
)
q2$coefficients <- sample(salt.rq.coef, 1)[[1]][,11]
plot (pred.s(-10:30, 30:70, 1, 1), ylim=c(0,12))
abline(h=3)
lines(pred.s(-10:30, 30:70, 1, 3), col = "blue")
lines(pred.s(-10:30, 30:70, 1, 5), col = "red")
lines(pred.s(-10:30, 30:70, 1, 1), col = "green")
#save(salt.rq, file="./Lagtimes/salt.rq.rda")
#save(salt.rq.coef, file="./Lagtimes/salt.rq.coef.rda")

#save(salt.rq, file="./Models/IMPACTncd/Lagtimes/salt.rq.rda")
#save(salt.rq.coef, file="./Models/IMPACTncd/Lagtimes/salt.rq.coef.rda")

# PA -----------------------------------------------------
# PA (imputation is ignored, as it should for dependent variable)
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[, mean(a30to06m, na.rm = T), by = age][, plot(age, V1, col="red")]
HSE.nat[, mean(a30to06m, na.rm = T), by = age][, points(age, V1)]

HSE.ts[, mean(a30to06m, na.rm = T), by = year][, plot(year, V1, col="red")]
HSE.nat[, mean(a30to06m, na.rm = T), by = year][, points(year, V1)]

HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age>19 & year > -7 & wt_int>0 & is.na(a30to06m)== F & is.na(qimd)==F) # years before 2006 were excluded for more realistic projections
pp <- svyby(~a30to06m, by=~age, design=HSE.ts.srv.int, svymean)
scatter.smooth(pp, family = "gaussian", ylim=c(0,6))
lines(y=predict(svyglm(a30to06m~age + I(age^2) + I(age^3), family=quasipoisson(), design=HSE.ts.srv.int), data.frame(age=20:100), type="response"), x=20:100, col="red")

pp <- svyby(~a30to06m, by=~year, design=HSE.ts.srv.int, svymean)
scatter.smooth(pp, family = "gaussian", ylim=c(0,6), xlim=c(-10,1))
lines(y=predict(svyglm(a30to06m~year, family=quasipoisson(), design=HSE.ts.srv.int), data.frame(year=-10:50), type="response"), x=-10:50, col="red")
scatter.smooth(pp, family = "gaussian", ylim=c(0,6), xlim=c(-10,50))
lines(y=predict(svyglm(a30to06m~year, family=quasipoisson(), design=HSE.ts.srv.int), data.frame(year=-10:50), type="response"), x=-10:50, col="red")
lines(y=predict(svyglm(a30to06m~log(year+25), family=quasipoisson(), design=HSE.ts.srv.int), data.frame(year=-10:50), type="response"), x=-10:50, col="green")

HSE.ts[bmival<16 & age>19, bmival := 16]
HSE.ts[bmival>40 & age>19, bmival := 40]
HSE.ts[age>85, age := 85]
HSE.ts[, a30to06m := ordered(a30to06m)]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age>19 & year > -7 & wt_int>0 & is.na(a30to06m)== F & is.na(qimd)==F) # years before 2006 were excluded for more realistic projections
HSE.nat[, hist(a30to06m)]
HSE.ts[, hist(as.integer(a30to06m))]


tt = copy(HSE.ts[age > 19 &
                   year > -7 &
                   wt_int > 0 &
                   is.na(a30to06m)== F &
                   is.na(qimd) == F, ])
tt[, a30to06m := ordered(a30to06m)]


# fit ordinal multinomial regression model and ignore sampling design
# I need to scale count vars for stepaic to work
# also tried with rcs and ns and got worse predictions
ttt <- data.frame(a30to06m=tt[,a30to06m],scale(tt[,.(year,age)]), tt[,.(sex,qimd,wt_int)])
mod1 <- polr(a30to06m~(log(year+25) + age  + sex + qimd)^2 + I(age^2) + I(age^3),
             data = ttt,
             weights = wt_int,
             method = "logistic",
             Hess = T)
summary(mod1)
mod2 <- stepAIC(mod1)

mod2 <- # best model based on aic
  polr(
    a30to06m ~ log(year + 25) + age + sex + qimd +
      I(age^2) + I(age^3) + log(year + 25):qimd + age:sex + age:qimd,
    data = tt,
    weights = wt_int,
    method = "logistic",
    Hess = T)

pa.svylr <- #apply formula of mod2 to svy
  svyolr(
    a30to06m ~ log(year + 25) + age + sex + qimd +
      I(age^2) + I(age^3) + log(year + 25):qimd + age:sex + age:qimd,
    design = HSE.ts.srv.int,
    method = "logistic")

# copy parameters of svy model to polr model so predict can work
for (k in intersect(names(mod2),  names(pa.svylr))) mod2[k] <- pa.svylr[k]

pa.svylr <- mod2
pa.svylr$data <- NULL
pa.svylr$lp <- NULL
pa.svylr$fitted.values <- NULL
#save(pa.svylr, file="./Lagtimes/pa.svylr.rda")

# bmi model  --------------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.nat[, mean(bmival, na.rm = T), by = age][, plot(age, V1)]
HSE.ts[, mean(bmival, na.rm = T), by = age][, lines(age, V1, col="red")]
HSE.nat[, mean(bmival, na.rm = T), by = year][, plot(year, V1)]
HSE.ts[, mean(bmival, na.rm = T), by = year][, lines(year, V1, col="red")]
HSE.ts[bmival<16 & age>19, bmival := 16]
HSE.ts[bmival>40 & age>19, bmival := 40]
HSE.ts[age>85, age:= 85]
HSE.ts.srv.nurse <- svydesign(id = ~psu, strata = ~cluster,
                              weights = ~wt_nurse, nest = F,
                              data=HSE.ts, check.strata = T)
HSE.ts.srv.nurse <- subset(HSE.ts.srv.nurse, age>19 &
                             wt_nurse > 0 &
                             is.na(qimd) == F &
                             is.na(bmival) == F &
                             is.na(a30to06m.imp) == F)
dt <- subset(HSE.ts, age>19 & wt_nurse > 0 &
               is.na(qimd) == F & is.na(bmival) == F &
               is.na(a30to06m.imp)== F)
dt[between(age, 20, 84), hist(bmival)]
dt[between(age, 20, 84), plot(density(log(bmival), na.rm = T))]
dt[between(age, 20, 84), plot(sort(log(bmival)), pch = ".")]
dt[between(age, 20, 84), boxplot(log(bmival) ~ age)]
dt[between(age, 20, 84), boxplot(log(bmival) ~ year)]
dt[between(age, 20, 84), boxplot(log(bmival) ~ a30to06m.imp)]
dt[between(age, 20, 84), boxplot(log(bmival) ~ qimd)]

pp <- svyby(~bmival, by=~age, design=HSE.ts.srv.nurse, svymean)
scatter.smooth(pp, ylim=c(0,30), family="gaussian", span = 1/3)
lines(y=predict(svyglm(bmival~age + I(age^2), family=gaussian(link="inverse"), design=HSE.ts.srv.nurse), data.frame(age=20:100), type="response"), x=20:100, col="green")
lines(y=predict(svyglm(bmival~age + I(age^2), family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(age=20:100), type="response"), x=20:100, col="red")
lines(y=predict(svyglm(bmival~age + I(age^2), family=inverse.gaussian, design=subset(HSE.ts.srv.nurse, sex == 1)), data.frame(age=20:100), type="response"), x=20:100, col="blue")

scatter.smooth(svyby(~bmival, by=~year, design=HSE.ts.srv.nurse, svymean), ylim=c(25,30), xlim=c(-10, 20))
lines(y=predict(svyglm(bmival~year, design=HSE.ts.srv.nurse), data.frame(year=-10:50), type="response"), x=-10:50, col="red")
lines(y=predict(svyglm(bmival~year, family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(year=-10:50), type="response"), x=-10:50, col="green")
lines(y=predict(svyglm(bmival~log(year + 25), family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(year=-10:50), type="response"), x=-10:50, col="green")
lines(y=predict(svyglm(bmival~year, family=inverse.gaussian, design=HSE.ts.srv.nurse), data.frame(year=-10:50), type="response"), x=-10:50, col="blue")
lines(y=predict(svyglm(bmival~year, family=gaussian(link="inverse"), design=HSE.ts.srv.nurse), data.frame(year=-10:50), type="response"), x=-10:50, col="purple")

scatter.smooth(svyby(~bmival, by=~a30to06m.imp, design=HSE.ts.srv.nurse, svymean), ylim=c(25,30), xlim=c(0, 7))
lines(y=predict(svyglm(bmival~a30to06m.imp + I(a30to06m.imp^2) + I(a30to06m.imp^3), design=HSE.ts.srv.nurse, family=gaussian(link="log")), data.frame(a30to06m.imp=0:7), type="response"), x=0:7, col="red")
lines(y=predict(svyglm(bmival~a30to06m.imp, design=HSE.ts.srv.nurse, family=gaussian(link="log")), data.frame(a30to06m.imp=0:7), type="response"), x=0:7, col="blue")

scatter.smooth(svyby(~bmival, by=~porftvg, design=HSE.ts.srv.nurse, svymean), ylim=c(25,30), xlim=c(0, 8))

bmi.svylm <- svyglm(bmival~log(year + 25) + age + sex + qimd + a30to06m.imp + I(age^2),
                    method = "glm.fit2",
                    family=gaussian(link = "log"),
                    design = HSE.ts.srv.nurse)
anova(bmi.svylm) # fv are not significant. literature support this

bmi.svylm2 <- svyglm(bmival~(log(year + 25) + age + sex + qimd + a30to06m.imp)^2 +
                       I(age^2) + I(a30to06m.imp^2) + I(a30to06m.imp^3),
                     method = "glm.fit2",
                     family=gaussian(link = "log"),
                     design = HSE.ts.srv.nurse)
anova(bmi.svylm2)
anova(bmi.svylm, bmi.svylm2)
AIC(bmi.svylm, bmi.svylm2)

bmi.svylm3 <- stepAIC(bmi.svylm2)
anova(bmi.svylm3)

anova(bmi.svylm2, bmi.svylm3)
anova(bmi.svylm, bmi.svylm3) #The simplest model is not signifficantly different than the complex ones
AIC(bmi.svylm, bmi.svylm2, bmi.svylm3)


bmi.svylm$deviance/bmi.svylm$df.null
1 - bmi.svylm$deviance/bmi.svylm$null.deviance # R^2

bmi.svylm$data <- NULL
bmi.svylm$survey.design <- NULL
bmi.svylm$qr <- NULL
bmi.svylm$residuals <- NULL
bmi.svylm$y <- NULL
bmi.svylm$linear.predictors <- NULL
bmi.svylm$fitted.values <- NULL
bmi.svylm$effects <- NULL
bmi.svylm$weights <- NULL
bmi.svylm$prior.weights <- NULL
#save(bmi.svylm, file="./Lagtimes/bmi.svylm.rda")

# Alcohol model ------------------------------------------------------


# SBP model ------------------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.nat[, mean(omsysval, na.rm = T), by = age][, plot(age, V1)]
HSE.ts[, mean(omsysval, na.rm = T), by = age][, lines(age, V1, col="red")]
HSE.nat[, mean(omsysval, na.rm = T), by = year][, plot(year, V1)]
HSE.ts[, mean(omsysval, na.rm = T), by = year][, lines(year, V1, col="red")]
HSE.ts[bmival<16 & age>19, bmival := 16]
HSE.ts[bmival>40 & age>19, bmival := 40]
HSE.ts[age>85, age:= 85]
HSE.ts[, cigst2 := mapvalues(cigst1,  c(4:1 ), c(1,0,0,0))]
HSE.ts.srv.nurse <- svydesign(id=~psu, strata =~cluster, weights = ~wt_nurse, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.nurse <- subset(HSE.ts.srv.nurse,
                           age>19 & wt_nurse>0 & is.na(omsysval)== F &
                             is.na(qimd)== F & is.na(bmival)== F &
                             is.na(cigst2)== F &
                             is.na(a30to06m.imp) == F)
dt <- subset(HSE.ts, age>19 & wt_nurse > 0 &
               is.na(qimd) == F & is.na(bmival) == F &
               is.na(a30to06m.imp)== F)
dt[between(age, 20, 84), hist(omsysval)]
dt[between(age, 20, 84), plot(density(log(omsysval), na.rm = T))]

pp <- svyby(~omsysval, by=~age, design=HSE.ts.srv.nurse, svymean)
scatter.smooth(pp, ylim=c(100,150))
lines(y=predict(svyglm(omsysval~age + I(age^2), family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(age=20:100), type="response"), x=20:100, col="green")

bb <- svyby(~omsysval, by=~round(bmival), design=HSE.ts.srv.nurse, svymean)
scatter.smooth(bb, ylim=c(100,150))
lines(y=predict(svyglm(omsysval~ bmival + I(bmival^2), family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(bmival=10:80), type="response"), x=10:80, col="green")
lines(y=predict(svyglm(omsysval~ bmival + I(bmival^2) + I(bmival^3), family=gaussian, design=HSE.ts.srv.nurse), data.frame(bmival=10:80), type="response"), x=10:80, col="red")
lines(y=predict(svyglm(omsysval~ log(bmival) + I(log(bmival)^2), family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(bmival=10:80), type="response"), x=10:80, col="blue")

scatter.smooth(svyby(~omsysval, by=~year, design=HSE.ts.srv.nurse, svymean), family = "gaussian", ylim=c(120,140), xlim=c(-10, 1))

scatter.smooth(svyby(~omsysval, by=~year, design=HSE.ts.srv.nurse, svymean), ylim=c(90,140), xlim=c(-10, 20))
lines(y=predict(svyglm(omsysval~year, design=HSE.ts.srv.nurse), data.frame(year=-10:100), type="response"), x=-10:100, col="red")
lines(y=predict(svyglm(omsysval~I((year+50)^-1), family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(year=-10:100), type="response"), x=-10:100, col="red") # play with+50 to adjust limit
lines(y=predict(svyglm(omsysval~I(log(year+25)), family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(year=-10:100), type="response"), x=-10:100, col="blue")

scatter.smooth(svyby(~omsysval, by=~porftvg.imp, design=HSE.ts.srv.nurse, svymean), ylim=c(125,135), xlim=c(0, 8))
lines(y=predict(svyglm(omsysval~porftvg.imp, family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(porftvg.imp=0:8), type="response"), x=0:8, col="red")
lines(y=predict(svyglm(omsysval~porftvg.imp + I(porftvg.imp^2), family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(porftvg.imp=0:8), type="response"), x=0:8, col="green")

scatter.smooth(svyby(~omsysval, by=~cigst2, design=HSE.ts.srv.nurse, svymean), ylim=c(90,140), xlim=c(1, 4))
scatter.smooth(svyby(~omsysval, by=~a30to06m.imp, design=HSE.ts.srv.nurse, svymean), ylim=c(90,140), xlim=c(0, 7))
lines(y=predict(svyglm(omsysval~a30to06m.imp, family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(a30to06m.imp=0:7), type="response"), x=0:7, col="red")
lines(y=predict(svyglm(omsysval~a30to06m.imp + I(a30to06m.imp^2), family=gaussian(link="log"), design=HSE.ts.srv.nurse), data.frame(a30to06m.imp=0:7), type="response"), x=0:7, col="green")


sbp.svylm <- svyglm(omsysval ~ I(age^2) + log(year+25) +
                      age + sex + qimd + log(bmival) +
                      I(log(bmival)^2) + cigst2 + a30to06m.imp +
                      I(a30to06m.imp^2) + I(porftvg.imp^2),
                    family = gaussian(link="log"),
                    method = "glm.fit2",
                    design = HSE.ts.srv.nurse)
anova(sbp.svylm)

sbp.svylm1 <- stepAIC(sbp.svylm)
anova(sbp.svylm, sbp.svylm1)
anova(sbp.svylm1)

sbp.svylm2 <- svyglm(omsysval ~ (I(age^2) + log(year + 25) + sex +
                       qimd + log(bmival) + cigst2 + I(a30to06m.imp^2))^2 + I(log(bmival)^2),
                    family = gaussian(link="log"),
                    method = "glm.fit2",
                    design = HSE.ts.srv.nurse)
anova(sbp.svylm1, sbp.svylm2)

sbp.svylm3 <- stepAIC(sbp.svylm2)
anova(sbp.svylm2, sbp.svylm3)

AIC(sbp.svylm, sbp.svylm2, sbp.svylm3)

sbp.svylm <- sbp.svylm3
sbp.svylm$deviance/sbp.svylm$df.null
1 - sbp.svylm$deviance/sbp.svylm$null.deviance # R^2
# It seems like increase PA does not influence sbp
sbp.svylm$data <- NULL
sbp.svylm$survey.design <- NULL
sbp.svylm$qr <- NULL
sbp.svylm$residuals <- NULL
sbp.svylm$y <- NULL
sbp.svylm$linear.predictors <- NULL
sbp.svylm$fitted.values <- NULL
sbp.svylm$effects <- NULL
sbp.svylm$weights <- NULL
sbp.svylm$prior.weights <- NULL
# save(sbp.svylm, file="./Lagtimes/sbp.svylm.rda")


# chol model -----------------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.nat[, mean(cholval1, na.rm = T), by = age][, plot(age, V1)]
HSE.ts[, mean(cholval1, na.rm = T), by = age][, lines(age, V1, col="red")]
HSE.nat[, mean(cholval1, na.rm = T), by = year][, plot(year, V1)]
HSE.ts[, mean(cholval1, na.rm = T), by = year][, points(year, V1, col="red")]
HSE.ts[bmival<16 & age>19, bmival := 16]
HSE.ts[bmival>40 & age>19, bmival := 40]
HSE.ts[age>85, age:= 85]
HSE.ts.srv.blood <- svydesign(id=~psu, strata =~cluster, weights = ~wt_blood, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.blood <- subset(HSE.ts.srv.blood, age>19 & wt_blood>0 & is.na(cholval1)== F & is.na(qimd)== F & is.na(bmival)== F & is.na(a30to06m.imp)==F & is.na(porftvg.imp)==F) #
dt <- subset(HSE.ts, age>19 & wt_nurse > 0 &
               is.na(qimd) == F & is.na(cholval1) == F &
               is.na(a30to06m.imp)== F)
dt[between(age, 20, 84), hist(cholval1)]
dt[between(age, 20, 84), plot(density((cholval1), na.rm = T))]
dt[between(age, 20, 84), plot(density(log(cholval1), na.rm = T))]
dt[between(age, 20, 84), plot(sort(log(cholval1)), pch = ".")]


scatter.smooth(svyby(~cholval1, by=~age, design=HSE.ts.srv.blood, svymean), ylim=c(3,7))
lines(y=predict(svyglm(cholval1~age + I(age^2) + I(age^3), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(age=20:100), type="response"), x=20:100, col="red")
lines(y=predict(svyglm(cholval1~age + I(age^2), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(age=20:100), type="response"), x=20:100, col="green")
lines(y=predict(svyglm(cholval1~log(age) + I(log(age)^2), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(age=20:100), type="response"), x=20:100, col="blue")

scatter.smooth(svyby(~cholval1, by=~round(bmival), design=HSE.ts.srv.blood, svymean), ylim=c(3,7))
lines(y=predict(svyglm(cholval1~bmival + I(1/bmival), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(bmival=10:80), type="response"), x=10:80, col="red")
lines(y=predict(svyglm(cholval1~bmival + I(bmival^2), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(bmival=10:80), type="response"), x=10:80, col="green")


scatter.smooth(svyby(~cholval1, by=~year, design=HSE.ts.srv.blood, svymean), family="gaussian", ylim=c(4,6), xlim=c(-10, 3))

scatter.smooth(svyby(~cholval1, by=~year, design=HSE.ts.srv.blood, svymean), family="gaussian",ylim=c(4,6), xlim=c(-10, 20))

lines(y=predict(svyglm(cholval1~year, family = gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(year=-10:100), type="response"), x=-10:100, col="red")
lines(y=predict(svyglm(cholval1~I(log(year+25)), family = gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(year=-10:100), type="response"), x=-10:100, col="blue")

scatter.smooth(svyby(~cholval1, by=~porftvg.imp, design=subset(HSE.ts.srv.blood, sex == 1), svymean), family="gaussian", ylim=c(3,7), xlim=c(0, 8))

scatter.smooth(svyby(~cholval1, by=~a30to06m.imp, design=subset(HSE.ts.srv.blood, sex == 1), svymean), family="gaussian", ylim=c(3,7), xlim=c(0, 7))

chol.svylm <- svyglm(cholval1~I(log(year+25)) + age + sex + qimd + porftvg.imp + a30to06m.imp +
                       bmival + I(age^2) + I(bmival^2),
                     family = gaussian(link="log"),
                     method = "glm.fit2",
                     design = HSE.ts.srv.blood)
anova(chol.svylm)

chol.svylm1 <- stepAIC(chol.svylm)
anova(chol.svylm1)
anova(chol.svylm, chol.svylm1)

chol.svylm2 <- svyglm(cholval1~(I(log(year+25)) + age + sex + qimd + bmival)^2 +
                        I(age^2) + I(bmival^2),
                      family = gaussian(link="log"),
                      method = "glm.fit2",
                      design = HSE.ts.srv.blood)
anova(chol.svylm1, chol.svylm2)
anova(chol.svylm2)

chol.svylm3 <- stepAIC(chol.svylm2)

anova(chol.svylm2, chol.svylm3)
AIC(chol.svylm, chol.svylm1, chol.svylm2, chol.svylm3)
anova(chol.svylm3)
chol.svylm <- chol.svylm3

chol.svylm$deviance/chol.svylm$df.null
1-chol.svylm$deviance/chol.svylm$null.deviance


chol.svylm$data <- NULL
chol.svylm$survey.design <- NULL
chol.svylm$qr <- NULL
chol.svylm$residuals <- NULL
chol.svylm$y <- NULL
chol.svylm$linear.predictors <- NULL
chol.svylm$fitted.values <- NULL
chol.svylm$effects <- NULL
chol.svylm$weights <- NULL
chol.svylm$prior.weights <- NULL
#save(chol.svylm, file="./Lagtimes/chol.svylm.rda")

# HDL model -----------------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.nat[, mean(hdlval1, na.rm = T), by = age][, plot(age, V1)]
HSE.ts[, mean(hdlval1, na.rm = T), by = age][, lines(age, V1, col="red")]
HSE.nat[, mean(hdlval1, na.rm = T), by = year][, plot(year, V1)]
HSE.ts[, mean(hdlval1, na.rm = T), by = year][, points(year, V1, col="red")]
HSE.ts[bmival<19 & age>19, bmival := 19]
HSE.ts[bmival>35 & age>19, bmival := 35]
HSE.ts[age>85, age:= 85]
HSE.ts[, tc.to.hdl := cholval1/hdlval1]
HSE.ts[, cigst2 := mapvalues(cigst1,  c(4:1 ), c(1,0,0,0))]

HSE.ts.srv.blood <- svydesign(id=~psu, strata =~cluster, weights = ~wt_blood, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.blood <- subset(HSE.ts.srv.blood, age>19 & wt_blood>0 & !is.na(cholval1) & !is.na(qimd) & !is.na(bmival) &
                             !is.na(hdlval1) & !is.na(a30to06m.imp) & !is.na(cigst1))

HSE.ts[tc.to.hdl < 10, hist(log(tc.to.hdl))]
HSE.ts[, scatter.smooth(cholval1, tc.to.hdl)]
HSE.ts[, mean(tc.to.hdl, na.rm = T), by = year][, plot(year, V1, col="red")]

scatter.smooth(svyby(~tc.to.hdl, by=~age, design=HSE.ts.srv.blood, svymean), ylim=c(3,5))
lines(y=predict(svyglm(tc.to.hdl~age + I(age^2) + I(age^3), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(age=20:100), type="response"), x=20:100, col="red")
lines(y=predict(svyglm(tc.to.hdl~log(age) + I(log(age)^2), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(age=20:100), type="response"), x=20:100, col="green")

scatter.smooth(svyby(~tc.to.hdl, by=~round(bmival), design=HSE.ts.srv.blood, svymean), ylim=c(3,5))
lines(y=predict(svyglm(tc.to.hdl~bmival + I(bmival^2)+ I(bmival^3), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(bmival=10:80), type="response"), x=10:80, col="red")
lines(y=predict(svyglm(tc.to.hdl~ log(bmival) + I(log(bmival)^2), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(bmival=10:80), type="response"), x=10:80, col="green")


scatter.smooth(svyby(~tc.to.hdl, by=~round(cholval1), design=HSE.ts.srv.blood, svymean), ylim=c(0,20))
lines(y=predict(svyglm(tc.to.hdl~log(cholval1), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(cholval1=0:16), type="response"), x=0:16, col="red")
lines(y=predict(svyglm(tc.to.hdl~cholval1 + I(cholval1^2), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(cholval1=0:16), type="response"), x=0:16, col="blue")
lines(y=predict(svyglm(tc.to.hdl~ cholval1 + I(log(cholval1)), family=gaussian(link="log"), design=HSE.ts.srv.blood), data.frame(cholval1=0:16), type="response"), x=0:16, col="green")

scatter.smooth(svyby(~tc.to.hdl, by=~as.integer(qimd), design=subset(HSE.ts.srv.blood, sex == 1), svymean), family="gaussian", ylim=c(3,5), xlim=c(0, 8))

scatter.smooth(svyby(~tc.to.hdl, by=~as.integer(cigst2), design=subset(HSE.ts.srv.blood, sex == 1), svymean), family="gaussian", ylim=c(3,5), xlim=c(0, 8))

scatter.smooth(svyby(~tc.to.hdl, by=~porftvg.imp, design=subset(HSE.ts.srv.blood, sex == 1), svymean), family="gaussian", ylim=c(3,5), xlim=c(0, 8))

scatter.smooth(svyby(~tc.to.hdl, by=~a30to06m.imp, design=subset(HSE.ts.srv.blood, sex == 1), svymean), family="gaussian", ylim=c(3,5), xlim=c(0, 7))

tctohdl.svylm <- svyglm(tc.to.hdl~ log(cholval1) + age + sex +
                          qimd + log(bmival) + log(age) +
                          a30to06m.imp + porftvg.imp +
                          cigst2 + I(log(bmival)^2) +
                          I(log(age)^2),
                        family = gaussian(link="log"),
                        method = "glm.fit2",
                        design = HSE.ts.srv.blood)
anova(tctohdl.svylm)

tctohdl.svylm1 <- stepAIC(tctohdl.svylm)
anova(tctohdl.svylm1)

tctohdl.svylm2 <- svyglm(tc.to.hdl~ (log(cholval1) + age + sex +
                           qimd + log(bmival) + log(age) +
                           a30to06m.imp +
                           cigst2)^2 + I(log(age)^2),
                         family = gaussian(link="log"),
                         method = "glm.fit2",
                         design = HSE.ts.srv.blood)
anova(tctohdl.svylm1, tctohdl.svylm2)

tctohdl.svylm3 <- stepAIC(tctohdl.svylm2)
anova(tctohdl.svylm2, tctohdl.svylm3)
AIC(tctohdl.svylm1, tctohdl.svylm2, tctohdl.svylm3)

# tctohdl.svylm <- tctohdl.svylm3 # Model 3 computationally expensive for predict
tctohdl.svylm <- tctohdl.svylm1
tctohdl.svylm$deviance/tctohdl.svylm$df.null
1-tctohdl.svylm$deviance/tctohdl.svylm$null.deviance

tctohdl.svylm$data <- NULL
tctohdl.svylm$survey.design$variables <- NULL
#save(tctohdl.svylm, file="./Lagtimes/tctohdl.svylm.rda")


# diab model -----------------------------------------------
# diab model (no time trend as will consider diabetes is totally dependant to BMI and age)
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.nat[, mean(diabtotr, na.rm = T), by = age][, plot(age, V1)]
HSE.ts[, mean(diabtotr, na.rm = T), by = age][, lines(age, V1, col="red")]
HSE.nat[, mean(diabtotr, na.rm = T), by = year][, plot(year, V1)]
HSE.ts[, mean(diabtotr, na.rm = T), by = year][, points(year, V1, col="red")]
HSE.ts[diabtotr == 1, diabtotr :=0]
HSE.ts[diabtotr == 2, diabtotr :=1]
HSE.ts[bmival<16 & age>19, bmival := 16]
HSE.ts[bmival>40 & age>19, bmival := 40]
HSE.ts[age>85, age:= 85] # different than the usual 85
HSE.ts.srv.blood <- svydesign(id=~psu, strata =~cluster, weights = ~wt_blood, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.blood <- subset(HSE.ts.srv.blood, age>19 & wt_blood>0 & is.na(diabtotr)== F & is.na(qimd)== F & is.na(bmival)== F & is.na(a30to06m.imp)== F) #

scatter.smooth(svyby(~diabtotr, by=~age, design=HSE.ts.srv.blood, svymean), ylim=c(0,0.25), xlim=c(20,100))
lines(y=predict(svyglm(diabtotr~age + I(age^2), family=quasibinomial(link="logit"), design=HSE.ts.srv.blood), data.frame(age=20:100), type="response"), x=20:100, col="purple")
lines(y=predict(svyglm(diabtotr~log(age) , family=quasibinomial(link="logit"), design=HSE.ts.srv.blood), data.frame(age=20:100), type="response"), x=20:100, col="red")

scatter.smooth(svyby(~diabtotr, by=~round(bmival), design=HSE.ts.srv.blood, svymean), ylim=c(0,0.25))
lines(y=predict(svyglm(diabtotr~bmival + I(bmival^2), family=quasibinomial(link="logit"), design=HSE.ts.srv.blood), data.frame(bmival=15:60), type="response"), x=15:60, col="blue")
lines(y=predict(svyglm(diabtotr~bmival , family=quasibinomial(link="logit"), design=HSE.ts.srv.blood), data.frame(bmival=15:60), type="response"), x=15:60, col="green")

scatter.smooth(svyby(~diabtotr, by=~year, design=HSE.ts.srv.blood, svymean), ylim=c(0,0.25), xlim=c(-9, 1))

scatter.smooth(svyby(~diabtotr, by=~year, design=HSE.ts.srv.blood, svymean), ylim=c(0,0.75), xlim=c(-11, 60))
lines(y=predict(svyglm(diabtotr~year, family=gaussian(link="identity"), design=HSE.ts.srv.blood), data.frame(year=-11:60), type="response"), x=-11:60, col="red")
lines(y=predict(svyglm(diabtotr~exp(-year+52), family=quasibinomial(link="logit"), design=HSE.ts.srv.blood), data.frame(year=-11:60), type="response"), x=-11:60, col="blue")
lines(y=predict(svyglm(diabtotr~I(log(year+11)), family=quasibinomial(link="probit"), design=HSE.ts.srv.blood), data.frame(year=-11:60), type="response"), x=-11:60, col="green")

scatter.smooth(svyby(~diabtotr, by=~a30to06m.imp, design=HSE.ts.srv.blood, svymean), ylim=c(0,0.25), xlim=c(0, 7))
lines(y=predict(svyglm(diabtotr~a30to06m.imp, family=quasibinomial, design=HSE.ts.srv.blood), data.frame(a30to06m.imp=0:7), type="response"), x=0:7, col="red")
lines(y=predict(svyglm(diabtotr~a30to06m.imp + I(a30to06m.imp^2), family=quasibinomial, design=HSE.ts.srv.blood), data.frame(a30to06m.imp=0:7), type="response"), x=0:7, col="blue")

diab.svylr <- svyglm(diabtotr~age + sex + qimd + bmival + a30to06m.imp + I(age^2),
                     design = HSE.ts.srv.blood,
                     family=quasibinomial,
                     method = "glm.fit")
anova(diab.svylr)

diab.svylr2 <- svyglm(diabtotr~(age + sex + qimd + bmival +
                        a30to06m.imp)^2 + I(age^2),
                      design = HSE.ts.srv.blood,
                      family=quasibinomial,
                      method = "glm.fit")

anova(diab.svylr2)
anova(diab.svylr, diab.svylr2)

diab.svylr3 <- stepAIC(diab.svylr)
anova(diab.svylr, diab.svylr3)
AIC(diab.svylr, diab.svylr2, diab.svylr3)
anova(diab.svylr3)

diab.svylr <- diab.svylr3
diab.svylr$deviance/diab.svylr$df.null
1-diab.svylr$deviance/diab.svylr$null.deviance

diab.svylr$data <- NULL
diab.svylr$survey.design <- NULL
diab.svylr$qr <- NULL
diab.svylr$residuals <- NULL
diab.svylr$y <- NULL
diab.svylr$linear.predictors <- NULL
diab.svylr$fitted.values <- NULL
diab.svylr$effects <- NULL
diab.svylr$weights <- NULL
diab.svylr$prior.weights <- NULL
#save(diab.svylr, file="./Lagtimes/diab.svylr.rda")


# smoking start model ------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
# numbers are small to restrict analysis in sha=2
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[,wt_int := wt_int*1000/sum(wt_int), by = year] # make pop of each hse =1000
HSE.ts[,sum(wt_int), by=year]
HSE.ts[cigst1==4 & between(age-startsmk,0,3), `:=`(year = year - (age-startsmk))]
HSE.ts[cigst1==4 & between(age-startsmk,0,3), `:=`(age = startsmk, smok.incid = 1)]
HSE.ts[between(cigst1, 2, 3) & between(endsmoke + smokyrs,0,3), `:=`(age  = age  - endsmoke - smokyrs,
    year = year - endsmoke - smokyrs, smok.incid = 1)]
HSE.ts[cigst1 == 1, smok.incid := 0]
HSE.ts[age > 85, age := 85]
HSE.ts[, qimd := mapvalues(qimd, 1:5, c(1,2,2,2,3))]
HSE.ts[, smok.incid:= factor(smok.incid)]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster,
                            weights = ~wt_int, nest=F,
                            data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age>15 & wt_int>0 &
                           !is.na(qimd) &
                           !is.na(smok.incid) & age < 50 &
                           between(year, -9, 0))

HSE.ts[, as.list(table(smok.incid)), keyby=year]
#HSE.ts.srv.int <- subset(HSE.ts.srv.int,qimd=="5")
# this gives the probability of a never smoker to become a smoker
# I had to do it this way because of the way the question was asked in HSE
pp <- svyby(~smok.incid, by=~age, design=HSE.ts.srv.int, svymean)

scatter.smooth(y=pp[[3]], x=as.numeric(pp[[1]]), ylim=c(0, 0.2), family = "gaussian")

lines(y=predict(svyglm(smok.incid~age , family=quasibinomial(link="logit"), design=HSE.ts.srv.int), data.frame(age=16:100), type="response"), x=16:100, col="red")
lines(y=predict(svyglm(smok.incid~age+ I(age^4)+ I(age^3) + I(age^2) , family=quasibinomial(link="logit"), design=HSE.ts.srv.int), data.frame(age=16:100), type="response"), x=16:100, col="blue")

pp <- svyby(~smok.incid, by=~qimd, design=HSE.ts.srv.int, svymean)
scatter.smooth(pp[3], ylim=c(0, 0.1), family = "gaussian")
lines(y=predict(svyglm(smok.incid~qimd , family=quasibinomial(link="probit"), design=HSE.ts.srv.int), data.frame(qimd=ordered(1:2)), type="response"), x= 1:2, col="red")
lines(y=predict(svyglm(smok.incid~qimd , family=quasibinomial(link="logit"), design=HSE.ts.srv.int), data.frame(qimd=ordered(1:2)), type="response"), x= 1:2, col="blue")

pp <- svyby(~smok.incid, by=~year, design=HSE.ts.srv.int, svymean)
scatter.smooth(pp[[1]], pp[[3]], ylim=c(0, .1), family = "gaussian", xlim=c(-10, 10))
lines(y=predict(svyglm(smok.incid~year, family=quasibinomial(link="logit"), design=HSE.ts.srv.int), data.frame(year=-10:50), type="response"), x=-10:50, col="red")
lines(y=predict(svyglm(smok.incid~year , family=quasibinomial(link="probit"), design=subset(HSE.ts.srv.int, qimd==1)), data.frame(year=-10:50), type="response"), x=-10:50, col="blue")
lines(y=predict(svyglm(smok.incid~year , family=quasibinomial(link="probit"), design=subset(HSE.ts.srv.int, qimd==2)), data.frame(year=-10:50), type="response"), x=-10:50, col="green")
lines(y=predict(svyglm(smok.incid~year , family=quasibinomial(link="probit"), design=subset(HSE.ts.srv.int, qimd==3)), data.frame(year=-10:50), type="response"), x=-10:50, col="yellow")
lines(y=predict(svyglm(smok.incid~year , family=quasibinomial(link="probit"), design=subset(HSE.ts.srv.int, qimd==4)), data.frame(year=-10:50), type="response"), x=-10:50, col="gray")
lines(y=predict(svyglm(smok.incid~year , family=quasibinomial(link="probit"), design=subset(HSE.ts.srv.int, qimd==5)), data.frame(year=-10:50), type="response"), x=-10:50, col="pink")


smok.start.svylr <- svyglm(smok.incid~ year + sex + age + qimd + I(age^4)+ I(age^3) + I(age^2), design = HSE.ts.srv.int,
                           family=quasibinomial(link="logit"),
                           method = "glm.fit2")
anova(smok.start.svylr)

smok.start.svylr1 <- stepAIC(smok.start.svylr)
anova(smok.start.svylr1)
anova(smok.start.svylr, smok.start.svylr1)

smok.start.svylr2 <- svyglm(smok.incid~ (year + sex + age + qimd)^2 + I(age^4)+ I(age^3) + I(age^2), design = HSE.ts.srv.int,
                           family=quasibinomial(link="logit"),
                           method = "glm.fit2")
anova(smok.start.svylr2)

smok.start.svylr3 <- stepAIC(smok.start.svylr2)
anova(smok.start.svylr2, smok.start.svylr3)
anova(smok.start.svylr3)

smok.start.svylr <- smok.start.svylr3
smok.start.svylr$deviance/smok.start.svylr$df.null
1-smok.start.svylr$deviance/smok.start.svylr$null.deviance

smok.start.svylr$data <- NULL
smok.start.svylr$survey.design <- NULL
smok.start.svylr$qr <- NULL
smok.start.svylr$residuals <- NULL
smok.start.svylr$y <- NULL
smok.start.svylr$linear.predictors <- NULL
smok.start.svylr$fitted.values <- NULL
smok.start.svylr$effects <- NULL
smok.start.svylr$weights <- NULL
smok.start.svylr$prior.weights <- NULL
#save(smok.start.svylr, file="./Lagtimes/smok.start.svylr.rda")

# Smok cessation -------------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[endsmoke < 1, smok.cess:= 1]
HSE.ts[cigst1 == 4, smok.cess:= 0]
HSE.ts[, qimd := mapvalues(qimd, 1:5, c(1,2,2,2,3))]
HSE.ts[age > 85, age := 85]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age>15 & wt_int>0 & is.na(qimd)==F & is.na(smok.cess)==F)
HSE.ts[, as.list(table(smok.cess)), keyby=year]
HSE.ts[, as.list(table(smok.cess)), keyby=year][, `1`/`0`]

pp<- svyby(~smok.cess, by=~year, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(pp, ylim=c(0, 0.4), xlim=c(-10,20), family = "gaussian")
lines(y=predict(svyglm(smok.cess~year , family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(year=-10:20), type="response"), x=-10:20, col="red")
lines(y=predict(svyglm(smok.cess~year+I(year^2) , family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(year=-10:50), type="response"), x=-10:50, col="blue")

aa<- svyby(~smok.cess, by=~age, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(aa, ylim=c(0, 0.2), xlim=c(10,100), family = "gaussian")
lines(y=predict(svyglm(smok.cess~age , family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(age=10:100), type="response"), x=10:100, col="red")
lines(y=predict(svyglm(smok.cess~age + I(age^2) + I(age^3)+ I(age^4), family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(age=10:100), type="response"), x=10:100, col="blue")
lines(y=predict(svyglm(smok.cess~age + (I(age^2) + I(age^3) + I(age^4))^2, family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(age=10:100), type="response"), x=10:100, col="green")

smok.cess.svylr <- svyglm(factor(smok.cess)~ (year + age + sex + qimd)^2,
                          design = HSE.ts.srv.int,
                          family=quasibinomial(link="logit"),
                          method = "glm.fit2")
anova(smok.cess.svylr)

smok.cess.svylr1 <- stepAIC(smok.cess.svylr)
anova(smok.cess.svylr1)
anova(smok.cess.svylr, smok.cess.svylr1)


smok.cess.svylr <- smok.cess.svylr1
smok.cess.svylr$deviance/smok.cess.svylr$df.null

smok.cess.svylr$data <- NULL
smok.cess.svylr$survey.design <- NULL
smok.cess.svylr$qr <- NULL
smok.cess.svylr$residuals <- NULL
smok.cess.svylr$y <- NULL
smok.cess.svylr$linear.predictors <- NULL
smok.cess.svylr$fitted.values <- NULL
smok.cess.svylr$effects <- NULL
smok.cess.svylr$weights <- NULL
smok.cess.svylr$prior.weights <- NULL
#save(smok.cess.svylr, file="./Lagtimes/smok.cess.svylr.rda")


# smoking relapse ------------------------------------------------
# Assuming perfect representation of hse with no sampling error,
# no deaths
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[, sum(wt_int), by=year]
HSE.ts[, wt_int := wt_int*1000/sum(wt_int), by=year] # make pop of each hse =1000
HSE.ts[, sum(wt_int), by=year]
#HSE.ts[endsmoke>10, endsmoke := 11]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int,
                            nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age > 15 & wt_int > 0 &
                           !is.na(qimd) & endsmoke >= 0 & cigst1 == 3)
#xx<- as.data.table(prop.table(svytable(~endsmoke+year+sex+qimd, HSE.ts.srv.int),2))
# Above was tested and year is not significant in the regression (and also too complicated)
pp <- as.data.table(svytable(~sex+qimd, HSE.ts.srv.int))
pp[, sex := as.factor(sex)]
pp[, qimd := as.ordered(qimd)]
xx <- as.data.table(svytable(~endsmoke+sex+qimd, HSE.ts.srv.int))
xx[, endsmoke:=as.numeric(endsmoke)]
xx[, sex := as.factor(sex)]
xx[, qimd := as.ordered(qimd)]
xx <- merge(xx, pp, by = c("sex", "qimd"),  all.x = T)
setkey(xx, sex, qimd, endsmoke)
xx[, pct:= N.x/N.y]
# smoothing to adjust for round number bias
xx[sex==1 & qimd == 1, plot(N.x)]
xx[sex==1 & qimd == 1 & endsmoke<20, lines(shift(predict(glm(N.x ~ log(endsmoke+1))), type = "lead"))]

xx[endsmoke<10, sm := predict(glm(N.x ~ log(endsmoke + 1))), by = .(sex, qimd)]
xx[endsmoke == 0 & sm < N.x, sm := N.x]

xx[, failure := shift(sm) - sm, by=.(sex,qimd)]

xx[,pct:=sm/(sm+failure)]
xx[sex==1 & qimd == 5, plot(failure)]
xx[sex==1 & qimd == 5,
   lines(predict(glm(failure ~ poly(log(endsmoke+1), 2) + endsmoke)))]
xx <- xx[endsmoke<10, ]
smok.cess.success <- glm(cbind(round(sm), round(failure))~ (poly(log(endsmoke), 2) + endsmoke +
                                                              sex + qimd)^2, xx[endsmoke > 0], family="binomial"(link="logit"))
require(MASS)
mod2 <- stepAIC(smok.cess.success)
predict(smok.cess.success, data.frame(endsmoke = 1:9, sex = factor(1, levels = c(1:2)), qimd = ordered(1, levels = c(1:5))), type="response", se.fit=F)
predict(mod2, data.frame(endsmoke = 1:9, sex = factor(1, levels = c(1:2)), qimd = ordered(1, levels = c(1:5))), type="response", se.fit=F)

xx[sex==2 & qimd==5, scatter.smooth(endsmoke,pct, ylim=c(0,1))]
lines(x=c(0:20), predict(smok.cess.success, data.frame(endsmoke = 0:20, sex = factor(2, levels = c(1:2)), qimd = ordered(5, levels = c(1:5))), type="response"), col="red")
smok.cess.success.log <- mod2
smok.cess.success.log$data <- NULL

#save(smok.cess.success.log, file="./Lagtimes/smok.cess.success.log.rda")

load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[,sum(wt_int), by=year]
HSE.ts[,wt_int := wt_int*1000/sum(wt_int), by=year] # make pop of each hse =1000
HSE.ts[,sum(wt_int), by=year]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int,
                            nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age > 15 & wt_int > 0 &
                           !is.na(qimd) & endsmoke >= 0 & cigst1 == 3)
pp <- as.data.table(svytable(~sex+qimd, HSE.ts.srv.int))
pp[, sex := as.factor(sex)]
pp[, qimd := as.ordered(qimd)]
xx <- as.data.table(svytable(~endsmoke+sex+qimd, HSE.ts.srv.int))
xx[, endsmoke:=as.numeric(endsmoke)]
xx[, sex := as.factor(sex)]
xx[, qimd := as.ordered(qimd)]
xx <- merge(xx, pp, by = c("sex", "qimd"),  all.x = T)
setkey(xx, sex, qimd, endsmoke)
xx[, pct:= N.x/N.y]
# smoothing to edjust for round number bias
xx[sex==1 & qimd == 5, plot(N.x)]
xx[endsmoke > 0 & sex==1 & qimd == 5, lines(predict(glm(N.x ~ poly(endsmoke^-1, 1))))]
xx[endsmoke > 0 & sex==1 & qimd == 5, lines(predict(glm(N.x ~ log(endsmoke))))]

xx[, sm := N.x]
xx[endsmoke > 0, sm := shift(predict(glm(N.x ~ I(endsmoke^-1))), type = "lead"), by = .(sex, qimd)]

xx[, failure := shift(sm) - sm, by=.(sex,qimd)]

xx[,pct:=sm/(sm+failure)]
xx <- xx[endsmoke < 10, ]

smok.cess.success.parabola <- glm(cbind(sm, failure)~ I(endsmoke^-1) +
                                    endsmoke + sex + qimd, xx[endsmoke > 0], family="quasibinomial"(link="logit"))

predict(smok.cess.success.parabola, data.frame(endsmoke = 0:9, sex = factor(1, levels = c(1:2)), qimd = ordered(1, levels = c(1:5))), type="response", se.fit=F)

xx[sex==2 & qimd==1, scatter.smooth(endsmoke,pct, ylim=c(0,1))]
lines(x=c(0:20), predict(smok.cess.success.parabola, data.frame(endsmoke = 0:20, sex = factor(2, levels = c(1:2)), qimd = ordered(1, levels = c(1:5))), type="response"), col="red")

smok.cess.success.parabola$data <- NULL
#save(smok.cess.success.parabola, file="./Lagtimes/smok.cess.success.parabola.rda")

# smoking prevalence ---------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[, smok.active := 0]
HSE.ts[between(cigst1, 2, 3), smok.active := 1]
HSE.ts[age>85, age := 85]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age>15 & age <21 & wt_int>0 & !is.na(qimd) & !is.na(smok.active))

pp<- svyby(~smok.active, by=~year, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(pp, ylim=c(0, 0.4), xlim=c(-10,30), family = "gaussian")
lines(y=predict(svyglm(smok.active~year , family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(year=-10:50), type="response"), x=-10:50, col="blue")
lines(y=predict(svyglm(smok.active~year+I(year^2) , family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(year=-10:50), type="response"), x=-10:50, col="red")
lines(y=predict(svyglm(smok.active~log(year+25) , family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(year=-10:50), type="response"), x=-10:50, col="green")

aa<- svyby(~smok.active, by=~age, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(aa, ylim=c(0, 0.4), xlim=c(15,30), family = "gaussian")
lines(y=predict(svyglm(smok.active~ age , family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(age=10:100), type="response"), x=10:100, col="red")
lines(y=predict(svyglm(smok.active~log(age) + I(log(age)^2) + age, family=quasibinomial(link="logit"),
                       design=HSE.ts.srv.int),
                data.frame(age=10:100), type="response"), x=10:100, col="blue")
lines(y=predict(svyglm(smok.active~I(age^-1) + I(age^2), family=quasibinomial(link="logit"),
                       design=HSE.ts.srv.int),
                data.frame(age=10:100), type="response"), x=10:100, col="green")

pp<- svyby(~smok.active, by=~qimd, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(pp, ylim=c(0, 0.4), xlim=c(1,5), family = "gaussian")
lines(y=predict(svyglm(smok.active~ qimd , family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(qimd=ordered(1:5)), type="response"), x=1:5, col="red")

smok.active.svylr <- svyglm(smok.active~ year + age + sex + qimd,
                            design = HSE.ts.srv.int,
                            family=quasibinomial(link="logit"),
                            method = "glm.fit2")
anova(smok.active.svylr)

smok.active.svylr2 <- svyglm(smok.active~ (year + age+ sex + qimd)^2,
                             design = HSE.ts.srv.int,
                             family=quasibinomial(link="logit"),
                             method = "glm.fit2")
anova(smok.active.svylr2)

smok.active.svylr3 <- stepAIC(smok.active.svylr2)
anova(smok.active.svylr3)
anova(smok.active.svylr, smok.active.svylr3)
smok.active.svylr <- smok.active.svylr3
smok.active.svylr$deviance/smok.active.svylr$df.null

smok.active.svylr$data <- NULL
smok.active.svylr$survey.design <- NULL
smok.active.svylr$qr <- NULL
smok.active.svylr$residuals <- NULL
smok.active.svylr$y <- NULL
smok.active.svylr$linear.predictors <- NULL
smok.active.svylr$fitted.values <- NULL
smok.active.svylr$effects <- NULL
smok.active.svylr$weights <- NULL
smok.active.svylr$prior.weights <- NULL
#save(smok.active.svylr, file="./Lagtimes/smok.active.svylr.rda")

# ex-smoking prevalence ------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[, smok.active := 0]
HSE.ts[between(cigst1, 2, 3), smok.active := 1]
HSE.ts[age>85, age := 85]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age>15 & age <21 & wt_int>0 & !is.na(qimd) & !is.na(smok.active))

pp <- svyby(~smok.active, by=~year, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(pp, ylim=c(0, 0.25), xlim=c(-10,30), family = "gaussian")
lines(y=predict(svyglm(smok.active~year , family=quasibinomial(link="logit"),
                       design=HSE.ts.srv.int),
                data.frame(year=-10:50), type="response"), x=-10:50, col="blue")
lines(y=predict(svyglm(smok.active~year+I(year^2) , family=quasibinomial(link="logit"),
                       design=HSE.ts.srv.int),
                data.frame(year=-10:50), type="response"), x=-10:50, col="red")
lines(y=predict(svyglm(smok.active~log(year+25) , family=quasibinomial(link="logit"),
                       design=HSE.ts.srv.int),
                data.frame(year=-10:50), type="response"), x=-10:50, col="green")

aa<- svyby(~smok.active, by=~age, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(aa, ylim=c(0, 0.4), xlim=c(16,20), family = "gaussian")
lines(y=predict(svyglm(smok.active~ age , family=quasibinomial(link="logit"), design=HSE.ts.srv.int),
                data.frame(age=10:100), type="response"), x=10:100, col="red")
lines(y=predict(svyglm(smok.active~log(age) + I(log(age)^2) + age, family=quasibinomial(link="logit"),
                       design=HSE.ts.srv.int),
                data.frame(age=10:100), type="response"), x=10:100, col="blue")
lines(y=predict(svyglm(smok.active~I(age^-1) + I(age^2), family=quasibinomial(link="logit"),
                       design=HSE.ts.srv.int),
                data.frame(age=10:100), type="response"), x=10:100, col="green")

aa<- svyby(~smok.active, by=~qimd, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(aa, ylim=c(0, 0.4), xlim=c(1,5), family = "gaussian")

smok.exactive.svylr <- svyglm(smok.active~ year + age + sex + qimd,
                              design = HSE.ts.srv.int,
                              family=quasibinomial(link="logit"),
                              method = "glm.fit2")
anova(smok.exactive.svylr)

smok.exactive.svylr2 <- svyglm(smok.active~ (year + age + sex + qimd)^2,
                               design = HSE.ts.srv.int,
                               family=quasibinomial(link="logit"),
                               method = "glm.fit2")
anova(smok.exactive.svylr2)

smok.exactive.svylr3 <- stepAIC(smok.exactive.svylr2)
anova(smok.exactive.svylr3)
anova(smok.exactive.svylr, smok.exactive.svylr3)

smok.exactive.svylr <- smok.exactive.svylr3
smok.exactive.svylr$deviance/smok.exactive.svylr$df.null

smok.exactive.svylr$data <- NULL
smok.exactive.svylr$survey.design <- NULL
smok.exactive.svylr$qr <- NULL
smok.exactive.svylr$residuals <- NULL
smok.exactive.svylr$y <- NULL
smok.exactive.svylr$linear.predictors <- NULL
smok.exactive.svylr$fitted.values <- NULL
smok.exactive.svylr$effects <- NULL
smok.exactive.svylr$weights <- NULL
smok.exactive.svylr$prior.weights <- NULL
#save(smok.exactive.svylr, file="./Lagtimes/smok.exactive.svylr.rda")

# F&V ------------------------------------------------------
load(file = "./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.nat[, mean(porftvg, na.rm = T), by = age][, plot(age, V1)]
HSE.ts[, mean(porftvg, na.rm = T), by = age][, lines(age, V1, col="red")]
HSE.nat[, mean(porftvg, na.rm = T), by = year][, plot(year, V1)]
HSE.ts[, mean(porftvg, na.rm = T), by = year][, lines(year, V1, col="red")]
HSE.ts[bmival<16 & age>19, bmival := 16]
HSE.ts[bmival>40 & age>19, bmival := 40]
HSE.ts[age > 85, age := 85]
HSE.ts[porftvg>8, porftvg := 8L]
HSE.ts[, porftvg := ordered(porftvg)]

HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age>19 & wt_int>0 & is.na(porftvg)== F & is.na(qimd)==F & year != -8)

tt = copy(HSE.ts[age>19 & wt_int>0 & is.na(porftvg)== F & is.na(qimd)==F])

tt[, porftvg := ordered(porftvg)]

require(MASS)
# I need to scale count vars for stepaic to work
ttt <-data.frame(porftvg=tt[,porftvg],scale(tt[,.(year,age)]), tt[,.(sex,qimd,wt_int)])
mod1 <- polr(porftvg~(log(year + 25) + age + sex + qimd)^2 + I(age^2) ,
             data = ttt,
             weights = wt_int,
             method = "logistic",
             Hess = T)
summary(mod1)
mod2 <- stepAIC(mod1)

mod2 <- # best model based on aic
  polr(
    porftvg ~ (log(year + 25) + age + sex + qimd)^2 + I(age^2),
    data = tt,
    weights = wt_int,
    method = "logistic",
    Hess = T)

fv.svylr <- #apply formula of mod2 to svy
  svyolr(
    porftvg ~ (log(year + 25) + age + sex + qimd)^2 + I(age^2),
    design = HSE.ts.srv.int,
    method = "logistic")

# copy parameters of svy model to polr model so predict can work
for (k in intersect(names(mod2),  names(fv.svylr))) mod2[k] <- fv.svylr[k]

fv.svylr <- mod2

pp <- svyby(~porftvg, by=~age, design=HSE.ts.srv.int, svymean)
scatter.smooth(pp, family = "gaussian", ylim=c(0,6))
lines(y=predict(svyglm(porftvg~age + I(age^2), family=quasipoisson(), design=HSE.ts.srv.int), data.frame(age=20:100), type="response"), x=20:100, col="red")

pp <- svyby(~porftvg, by=~year, design=HSE.ts.srv.int, svymean)
scatter.smooth(pp, family = "gaussian", ylim=c(3,4.5), xlim=c(-10,1))
scatter.smooth(pp, family = "gaussian", ylim=c(3,4.5), xlim=c(-10,50))
lines(y=predict(svyglm(porftvg~year, family=quasipoisson(), design=HSE.ts.srv.int), data.frame(year=-10:50), type="response"), x=-10:50, col="red")

pp<- svyby(~as.integer(porftvg), by=~round(bmival), design=HSE.ts.srv.int, svymean)
scatter.smooth(pp, family = "gaussian", ylim=c(0,6))
lines(y=predict(svyglm(porftvg~bmival + I(bmival^2) + I(bmival^3), family=quasipoisson(), design=subset(HSE.ts.srv.int, sex == 1)), data.frame(bmival=10:60), type="response"), x=10:60, col="red")
lines(y=predict(svyglm(porftvg~bmival + I(bmival^2) , family=quasipoisson(), design=subset(HSE.ts.srv.int, sex == 1)), data.frame(bmival=10:60), type="response"), x=10:60, col="blue")

fv.svylr$deviance/fv.svylr$df.null
1-fv.svylr$deviance/fv.svylr$null.deviance

fv.svylr$data <- NULL
fv.svylr$lp <- NULL
fv.svylr$fitted.values <- NULL
#save(fv.svylr, file="./Lagtimes/fv.svylr.rda")

# FamCVD --------------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.nat[, mean(porftvg, na.rm = T), by = age][, plot(age, V1)]
HSE.ts[, mean(porftvg, na.rm = T), by = age][, lines(age, V1, col="red")]
HSE.ts[famcvd == 2,famcvd := 0]
HSE.ts[, famcvd := factor(famcvd)]
HSE.ts[, agegroup := agegroup.fn(age)]
HSE.ts[, table(famcvd)]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, wt_int>0 & !is.na(famcvd) & !is.na(qimd))

famcvd.svylr <- svyglm(famcvd~age + I(age^2) + qimd,
                       design = HSE.ts.srv.int,
                       family=quasibinomial,
                       method = "glm.fit2")
anova(famcvd.svylr)
1-famcvd.svylr$deviance/famcvd.svylr$null.deviance

famcvd.svylr$data <- NULL
famcvd.svylr$survey.design$variables <- NULL
#save(famcvd.svylr, file="./Lagtimes/famcvd.svylr.rda")

# AF prevalence  -----------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[iregdef == 2,iregdef := 0]
HSE.ts[, iregdef := factor(iregdef)]
HSE.ts[, agegroup := agegroup.fn(age)]
HSE.nat[, table(iregdef)]
HSE.ts[, table(iregdef)]
HSE.ts[, cigst2 := mapvalues(cigst1,  c(4:1 ), c(1,0,0,0))]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, wt_int>0 & !is.na(iregdef) & !is.na(qimd) & !is.na(cigst1))

af.svylr <- svyglm(iregdef~ (age + qimd + cigst2 + sex + omsysval)^2,
                   design = HSE.ts.srv.int,
                   family=quasibinomial,
                   method = "glm.fit2")
anova(af.svylr)

af.svylr1 <- stepAIC(af.svylr)
anova(af.svylr1)

af.svylr <- svyglm(iregdef~ age + omsysval,
                   design = HSE.ts.srv.int,
                   family=quasibinomial,
                   method = "glm.fit2")
anova(af.svylr)
af.svylr$data <- NULL
af.svylr$survey.design$variables <- NULL
#save(af.svylr, file="./Lagtimes/af.svylr.rda")

# Kidney disease prevalence  -----------------------------------
# kidfailgp
# 1 'Normal: eGFR 60+ ml/min/1.73m2 and normal albuminuria'
# 2 'Stage 1: eGFR 90+ ml/min/1.73m2 and micro- or macro-albuminuria'
# 3 'Stage 2: eGFR 60-89 ml/min/1.73m2 and micro- or macro-albuminuria'
# 4 'Stage 3a/3b: eGFR 30-59 ml/min/1.73m2 regardless of albuminuria'
# 5 'Stage 4/5: eGFR less than 30 ml/min/1.73m2 regardless of albuminuria' ONLY this level is considered in QRisk

load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
#HSE.ts <- HSE.ts[sha==2] # North West SHA. very few cases
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[, table(kidfailgp)]
HSE.ts[, kidfailgp := factor(kidfailgp)]
HSE.ts[, kidfailgp := mapvalues(kidfailgp, 1:5, c(0,0,0,0,1))]
HSE.ts[, agegroup := agegroup.fn(age)]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, wt_int>0 & !is.na(kidfailgp) & omsysval > 0)

kidfailgp.svylr <- svyglm(kidfailgp~ age + omsysval,
                        design = HSE.ts.srv.int,
                        family=quasibinomial,
                        method = "glm.fit2")
anova(kidfailgp.svylr)

kidfailgp.svylr$data <- NULL
kidfailgp.svylr$survey.design$variables <- NULL
#save(kidfailgp.svylr, file="./Lagtimes/kidfailgp.svylr.rda")

# BP medication  ---------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[, table(bpmedd)]
HSE.ts[, bpmedd := factor(bpmedd)]
HSE.ts[, agegroup := agegroup.fn(age)]
HSE.ts.srv.nurse <- svydesign(id=~psu, strata =~cluster, weights = ~wt_nurse, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.nurse <- subset(HSE.ts.srv.nurse, wt_nurse>0 & !is.na(bpmedd) & !is.na(qimd) & !is.na(sex) & !is.na(omsysval))
svytable(~bpmedd, subset(HSE.ts.srv.nurse, age > 30 & age < 84))

bpmed.svylr <- svyglm(bpmedd~ year + age + I(age^2) + sex + qimd + omsysval,
                      design = HSE.ts.srv.nurse,
                      family=quasibinomial,
                      method = "glm.fit2")
anova(bpmed.svylr)

bpmed.svylr1 <- stepAIC(bpmed.svylr)
anova(bpmed.svylr1)

bpmed.svylr2 <- svyglm(bpmedd~ (year + age + I(age^2) + sex + qimd + omsysval)^2,
                      design = HSE.ts.srv.nurse,
                      family=quasibinomial,
                      method = "glm.fit2")
bpmed.svylr3 <- stepAIC(bpmed.svylr2)
anova(bpmed.svylr3)
anova(bpmed.svylr1, bpmed.svylr3)
AIC(bpmed.svylr, bpmed.svylr1,bpmed.svylr2, bpmed.svylr3)

bpmed.svylr <- bpmed.svylr3
1-bpmed.svylr$deviance/bpmed.svylr$null.deviance

bpmed.svylr$data <- NULL
bpmed.svylr$survey.design <- NULL
bpmed.svylr$qr <- NULL
bpmed.svylr$residuals <- NULL
bpmed.svylr$y <- NULL
bpmed.svylr$linear.predictors <- NULL
bpmed.svylr$fitted.values <- NULL
bpmed.svylr$effects <- NULL
bpmed.svylr$weights <- NULL
bpmed.svylr$prior.weights <- NULL
#save(bpmed.svylr, file="./Lagtimes/bpmed.svylr.rda")

# Undiagnosed diab --------------------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[, diabtotr := factor(diabtotr - 1)] #everybody
HSE.ts[diabete2 == 2, diabete2 := 0] #only diagnosed
HSE.ts[diabtotr == "1" & diabete2 == "0", undiag.diab := 1L]
HSE.ts[diabtotr == "1" & diabete2 == "1", undiag.diab := 0L]
HSE.ts[, undiag.diab :=  factor(undiag.diab)]
HSE.ts[, diabete2 := factor(diabete2)]
HSE.ts[, agegroup := agegroup.fn(age)]
HSE.ts[, table(undiag.diab)]
HSE.ts.srv.blood <- svydesign(id=~psu, strata =~cluster, weights = ~wt_blood, nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.blood <- subset(HSE.ts.srv.blood, wt_blood > 0 & !is.na(qimd) & !is.na(sex) &
                             !is.na(undiag.diab)) #-7 is influential outliar
svytable(~undiag.diab, subset(HSE.ts.srv.blood, age > 30 & age < 84))
pp<- svyby(~undiag.diab, by=~year, design=HSE.ts.srv.blood, svymean)
scatter.smooth(x=pp$year, y=pp$undiag.diab1, family = "gaussian")

undiag.diab.svylr <- svyglm(undiag.diab ~ qimd, # no other factors were signif.
                            # nor interactions.
                            design = HSE.ts.srv.blood,
                            family=quasibinomial,
                            method = "glm.fit2")
anova(undiag.diab.svylr) # year non significant

undiag.diab.svylr$data <- NULL
undiag.diab.svylr$survey.design$variables <- NULL
#save(undiag.diab.svylr, file="./Lagtimes/undiag.diab.svylr.rda")

# smoking cigdyal ---------------------------------------------------------
load(file="./Lagtimes/HSE.ts.RData")
HSE.nat = copy(HSE.ts)
HSE.ts <- HSE.ts[sha==2] # North West SHA
setkey(HSE.nat, age)
setkey(HSE.ts, age)
HSE.ts[startsmk == 97, startsmk := NA]
HSE.ts[cigdyal == 97, cigdyal := NA]
HSE.ts[cigst1 == 4, smokyrs:= age - startsmk]
HSE.ts[smokyrs <0, smokyrs := 0L]
HSE.ts[age > 85, age := 85L]
HSE.ts[, year := as.integer(year)]
# HSE.ts[cigst1 == 3, `:=` (year = year - endsmoke,
#                           cigdyal = as.numeric(numsmok),
#                           cigst1 = 4L)]
HSE.ts[cigdyal>30, cigdyal := 30]
HSE.ts[cigst1 == 4 & cigdyal < 1, cigdyal := 1]
HSE.ts[cigst1 == "4", hist(cigdyal)]
HSE.ts[cigst1 == "4", cigdyalCat := cut(cigdyal, c(0, 5, 10, 20, 30), include.lowest = T,
                                        right = T, ordered_result = T)]
HSE.ts[cigst1 == "4", plot(cigdyalCat)]
HSE.ts[cigst1 == "4", mean(cigdyal, na.rm = T), keyby = .(cigdyalCat)]
HSE.ts[cigdyalCat == "[0,5]", hist(cigdyal)]
HSE.ts[cigdyalCat == "(5,10]", hist(cigdyal)]
HSE.ts[cigdyalCat == "(10,20]", hist(cigdyal)]
HSE.ts[cigdyalCat == "(20,30]", hist(cigdyal)]

#HSE.ts[, cigdyal := ordered(as.integer(cigdyal/4))]
HSE.ts.srv.int <- svydesign(id=~psu, strata =~cluster, weights = ~wt_int,
                            nest=F, data=HSE.ts, check.strata = T)
HSE.ts.srv.int <- subset(HSE.ts.srv.int, age > 15 & age < 85 & wt_int > 0 &
                           !is.na(smokyrs) &
                           !is.na(qimd) & !is.na(cigdyalCat) & cigst1 == 4)
tt = copy(HSE.ts[between(age, 16, 84) & wt_int > 0 &
                   !is.na(smokyrs) & year > -30 &
                   !is.na(qimd) & !is.na(cigdyalCat) & cigst1 == 4])

pp<- svyby(~cigdyal, by=~year, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(pp, ylim=c(0, 30), xlim=c(-40,50), family = "gaussian")
lines(y=predict(svyglm(cigdyal~I(year), design=HSE.ts.srv.int, family = "quasipoisson"),
                data.frame(year=-80:50), type="response"), x=-80:50, col="blue")
lines(y=predict(svyglm(cigdyal~log(year+50) , design=HSE.ts.srv.int, family = "quasipoisson"),
                data.frame(year=-30:50), type="response"), x=-30:50, col="green")

pp<- svyby(~cigdyal, by=~age, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(pp, ylim=c(0, 30), xlim=c(20,84), family = "gaussian")
lines(y=predict(svyglm(cigdyal~age, design=HSE.ts.srv.int, family = "quasipoisson"),
                data.frame(age=15:84), type="response"), x=15:84, col="blue")
lines(y=predict(svyglm(cigdyal~age+I(age^2), design=HSE.ts.srv.int, family = "quasipoisson"),
                data.frame(age=15:84), type="response"), x=15:84, col="red")
lines(y=predict(svyglm(cigdyal~log(age) + I(log(age)^2), design=HSE.ts.srv.int, family = "quasipoisson"),
                data.frame(age=15:84), type="response"), x=15:84, col="green")

pp<- svyby(~cigdyal, by=~smokyrs, design=HSE.ts.srv.int, svymean, na.rm=T)
scatter.smooth(pp, ylim=c(0, 20), xlim=c(0,84), family = "gaussian")
lines(y=predict(svyglm(cigdyal~smokyrs, design=HSE.ts.srv.int, family = "quasipoisson"),
                data.frame(smokyrs=0:84), type="response"), x=0:84, col="blue")
lines(y=predict(svyglm(cigdyal~smokyrs+I(smokyrs^2), design=HSE.ts.srv.int, family = "quasipoisson"),
                data.frame(smokyrs=0:84), type="response"), x=0:84, col="red")
lines(y=predict(svyglm(cigdyal~log(smokyrs+1) + I(log(smokyrs+1)^2), design=HSE.ts.srv.int, family = "quasipoisson"),
                data.frame(smokyrs=0:84), type="response"), x=0:84, col="green")

smok.cigdyal.svylm <- svyglm(cigdyal~ year + sex + qimd + smokyrs + I(smokyrs^2) + I(age^2),# age excluded due to colinearity
                             design = HSE.ts.srv.int,
                             family = "quasipoisson",
                             method = "glm.fit2")
anova(smok.cigdyal.svylm)

smok.cigdyal.svylm2 <- svyglm(cigdyal~ (year + sex + qimd + smokyrs)^2 + I(smokyrs^2) + I(age^2),
                              design = HSE.ts.srv.int,
                              family = "quasipoisson",
                              method = "glm.fit2")
anova(smok.cigdyal.svylm2)
anova(smok.cigdyal.svylm, smok.cigdyal.svylm2)

smok.cigdyal.svylm3 <-stepAIC(smok.cigdyal.svylm2)
anova(smok.cigdyal.svylm3)
anova(smok.cigdyal.svylm, smok.cigdyal.svylm3) # the simpler model not signifficantly different



smok.cigdyal.svylm$deviance/smok.cigdyal.svylm$df.null

smok.cigdyal.svylm$data <- NULL
smok.cigdyal.svylm$survey.design <- NULL
smok.cigdyal.svylm$qr <- NULL
smok.cigdyal.svylm$residuals <- NULL
smok.cigdyal.svylm$y <- NULL
smok.cigdyal.svylm$linear.predictors <- NULL
smok.cigdyal.svylm$fitted.values <- NULL
smok.cigdyal.svylm$effects <- NULL
smok.cigdyal.svylm$weights <- NULL
smok.cigdyal.svylm$prior.weights <- NULL
#save(smok.cigdyal.svylm, file="./Lagtimes/smok.cigdyal.svylm.rda")

ttt <- data.frame(cigdyalCat=tt[,cigdyalCat],scale(tt[,.(year, age, smokyrs)]), tt[,.(sex,qimd,wt_int)])
mod1 <- polr(cigdyalCat ~ year + sex + qimd + age + smokyrs +I(age^2),
             data = ttt,
             weights = wt_int,
             method = "logistic",
             Hess = T)
summary(mod1)
mod2 <- stepAIC(mod1)

cigdyal.svylr <- #apply formula of mod2 to svy
  svyolr(
    cigdyalCat ~ year + sex + qimd + age + smokyrs + I(age^2),
    design = HSE.ts.srv.int,
    method = "logistic")

# copy parameters of svy model to polr model so predict can work
for (k in intersect(names(mod2),  names(cigdyal.svylr))) mod2[k] <- cigdyal.svylr[k]

cigdyal.svylr <- mod2

cigdyal.svylr$data <- NULL
cigdyal.svylr$lp <- NULL
cigdyal.svylr$fitted.values <- NULL
#save(cigdyal.svylr, file="./Lagtimes/cigdyal.svylr.rda")

HSE.ts[cigst1 == 4, meansd(cigdyal)]
HSE.ts[cigst1 == 4, plot(density(cigdyal, na.rm = T))]
meansd(rnbinom(1e5, mu = 13.3, size = 2.9))
lines(density(rnbinom(1e5, mu = 13.3, size = 2.9)))
meansd(c(rnbinom(1e5, mu = 11.3, size = 2.9), rnbinom(1e5, mu = 15.3, size = 2.9)))
HSE.ts[cigst1 == 4 & qimd == "5", meansd(cigdyal)]
meansd(rnbinom(1e5, mu = 14.26, size = 2.9))

y1 <- HSE.ts[cigst1 == 4, density(cigdyal, na.rm = T)]
y2 <- density(rnbinom(1e5, mu = 13.3, size = 2.9))
auc(y1$x, y1$y, type = "spline", to = 30)
auc(y2$x, y2$y, type = "spline", to = 30)
