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


setwd("/home/ckyprid/My Models/IMPACTncd_Engl/")
univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- TRUE
seed                <- 37L


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

dt <- HSE_ts[wt_int > 0, .(
  wt_int, smok_quit_yrs, smok_status,
  agegrp10, year, age, sex, qimd, ethnicity, sha)]
dt[, wt_int := wt_int*10000/sum(wt_int), by = year] # make pop of each hse =10000

dt[smok_quit_yrs < 20, sum(wt_int), keyby = smok_quit_yrs][, .(V1/first(V1))][, plot(V1)]
tt <- dt[smok_quit_yrs < 20, sum(wt_int), keyby = .(smok_quit_yrs, sex, qimd)]
setkey(tt, smok_quit_yrs)
tt[, pr := V1/first(V1), by = .(sex, qimd)]
m1 <- glm(pr~log(smok_quit_yrs) + sex + qimd, data = tt[smok_quit_yrs > 0, ])
newdata <- CJ(smok_quit_yrs = 1:100, sex = unique(dt$sex), qimd = unique(dt$qimd))
newdata[, pr := predict(m1, .SD, type = "response")]
newdata[pr < 0, pr := 0]
newdata[, pr := pr/shift(pr, fill = 1), by = .(sex, qimd)]
newdata[is.na(pr), pr := 1]

# Calculate probability of relapse
newdata[, pr := 1 - pr]
# Assume whoever achieves 15 years of cessation is not at risk of relapse
newdata[smok_quit_yrs > 15, pr := 0]

write_fst(newdata, "./inputs/exposure_distributions/smok_relapse_table.fst", 100L)

print("Table saved")

