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
setwd("/home/ckyprid/My Models/IMPACTncd_Engl/")

if (!require(CKutils)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
dependencies(
  c(
    "fst",
    "future",
    "future.apply",
    "data.table",
    "minerva"
  ),
  FALSE,
  FALSE,
  TRUE,
  FALSE
)
options(future.fork.enable = TRUE) # enable fork in Rstudio

plan(multiprocess, workers = 20L)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <-
    read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
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
HSE_ts <- HSE_ts[wt_blood > 0, ]

important_varnames <-
  c(
    "bmi",
    "tchol",
    "sbp",
    "active_days",
    "alcohol",
    # "totalwu",
    "frtpor",
    "vegpor",
    "hba1c",
    "smok_status_cont",
    "income_cont",
    "education_cont"
  ) # Only numeric vars

ll <-
  future_replicate(1e3L, minerva::mine(setDT(HSE_ts[sample(.N, 2e4, TRUE, wt_blood), .SD, .SDcols = important_varnames]),
                                       use = 'pairwise.complete.obs'))
mic <-
  rbindlist(lapply(ll[1,], as.data.table, keep.rownames = TRUE))# a dt with the MIC statistic correlation matrices
setcolorder(mic, order(names(mic)))
mic <-
  mic[, lapply(.SD, quantile, c(0.5, 0.025, 0.975)), keyby = rn]
correlation_structure <- mic[, lapply(.SD, function(x) {
  paste0(
    format(x[[1]], TRUE, 1, 2),
    " (95% CI:",
    format(x[[2]], TRUE, 1, 2),
    "-",
    format(x[[3]], TRUE, 1, 2),
    ")"
  )
}), by = rn]
fwrite(correlation_structure,
       "./preparatory_work/correlation_structure.csv")

# Post 2010 for totalwu
HSE_ts <- HSE_ts[year > 10]

important_varnames <-
  c(
    "bmi",
    "tchol",
    "sbp",
    # "active_days",
    "alcohol",
    "totalwu",
    "frtpor",
    "vegpor",
    "hba1c",
    "smok_status_cont",
    "income_cont",
    "education_cont"
  ) # Only numeric vars


ll <-
  future_replicate(1e3L, minerva::mine(setDT(HSE_ts[sample(.N, 1e4, TRUE, wt_blood), .SD, .SDcols = important_varnames]),
                                       use = 'pairwise.complete.obs'))
mic <-
  rbindlist(lapply(ll[1,], as.data.table, keep.rownames = TRUE))# a dt with the MIC statistic correlation matrices
setcolorder(mic, order(names(mic)))
mic <-
  mic[, lapply(.SD, quantile, c(0.5, 0.025, 0.975)), keyby = rn]
correlation_structure <- mic[, lapply(.SD, function(x) {
  paste0(
    format(x[[1]], TRUE, 1, 2),
    " (95% CI:",
    format(x[[2]], TRUE, 1, 2),
    "-",
    format(x[[3]], TRUE, 1, 2),
    ")"
  )
}), by = rn]
fwrite(correlation_structure,
       "./preparatory_work/correlation_structure_post_2010.csv")
