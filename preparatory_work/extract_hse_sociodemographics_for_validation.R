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

# Extract sociodemographic weights from HSE
# To be used to directly adjust POP to HSE sociodemographics -----
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  roxygen2::roxygenise("../CKutils/") # TODO remove before deployment
  remotes::install_github("ChristK/CKutils")
  library(CKutils)
}
if (!require(IMPACTncdEnglmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  roxygen2::roxygenise("./Rpackage/IMPACTncd_Engl_model_pkg/") # TODO remove before deployment
  remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/", force = TRUE)
  library(IMPACTncdEnglmisc)
}
dependencies(
  c(
    # "gamlss", # only necesary when fitting the models
    # "gamlss.tr", # only necesary when fitting the models
    # "mc2d", # only necessary for generating fixed_mc
    "doParallel",
    "doRNG",
    "foreach",
    # "mc2d", # for rpert()
    "gamlss.dist", # For distr in prevalence.R
    "dqrng",
    "qs",
    "fst",
    "wrswoR",
    "ggplot2",
    "cowplot",
    "viridis",
    "dichromat",
    "promises", # needed for %...>%
    "future",
    "data.table"
  ), TRUE, FALSE, FALSE, FALSE
)

if (file.exists("./preparatory_work/HSE_ts.fst")) {
  HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./preparatory_work/preprocess_HSE.R", local = TRUE)
}

HSE <- HSE_ts[!is.na(age) & !is.na(sex) & !is.na(qimd)  & !is.na(ethnicity) &
    !is.na(sha) & between(age, design$ageL, design$ageH)]

HSE[, year := year + 2000L]
HSE[, agegrp20 := NULL]
to_agegrp(HSE, 20L, 89L, "age", "agegrp20", to_factor = TRUE)

hse_wt <- CJ( # ignore year for simplicity
  # year = unique(HSE$year),
  age = unique(HSE$age),
  sex = unique(HSE$sex),
  qimd = unique(HSE$qimd),
  ethnicity = unique(HSE$ethnicity),
  sha = unique(HSE$sha)
  )

tt <- HSE[, .(hse_wt = sum(wt_int)), keyby = .(age, sex, qimd, ethnicity, sha)]

absorb_dt(hse_wt, tt)
setnafill(hse_wt, "c", 0, cols = "hse_wt")
write_fst(hse_wt, "./synthpop/hse_sociodemographics.fst", 100)
