## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2022 University of Liverpool, Chris Kypridemos
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

library(demography)
library(data.table)
library(CKutils)
library(Rcpp)
library(fst)
library(future.apply)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
if (!require(IMPACTncdEnglmisc)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/")
  library(IMPACTncdEnglmisc)
}
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multicore, workers = 10L)

hor <- 100L # maximum simulation horizon

# load data ----
deaths <- fread("./inputs/mortality/mortality.csv")[year > 2000]
deaths <- deaths[, lapply(.SD, function(x) gsub(",", "", x))]
deaths <- melt(deaths, id.vars = c("year", "dimd", "sex", "cod"),
     variable.name = "agegroup", value.name = "deaths")
deaths[, `:=` (year = as.integer(year),
               dimd = factor(dimd, levels = as.character(1:10)),
               sex = factor(sex, levels = c("Male", "Female"), labels = c("men", "women")),
               cod = factor(cod),
               deaths = as.integer(deaths)
               )]

pop <- fread("./inputs/mortality/population.csv")
pop <- pop[, lapply(.SD, function(x) gsub(",", "", x))]
pop <- melt(pop, id.vars = c("year", "dimd", "sex"),
               variable.name = "agegroup", value.name = "pops")
pop[, `:=` (year = as.integer(year),
            dimd = factor(dimd, levels = as.character(1:10)),
            sex = factor(sex, levels = c("Male", "Female"), labels = c("men", "women")),
            pops = as.integer(pops)
)]

levels(deaths$cod)
deaths[, cod := fcase(
  cod == "02 Malignant neoplasm of colon rectosigmoid junction and rectum (i.e.colorectal)", "colorectal_ca",
  cod == "05 Malignant neoplasm of trachea bronchus and lung", "lunc_ca",
  cod == "06 Malignant neoplasm of prostate", "prostate_ca",
  cod == "07 Malignant neoplasm of breast", "breast_ca",
  cod %in% c("01 Malignant neoplasm of liver and intrahepatic bile ducts",
             "03 Malignant neoplasms of digestive organs except liver and intrahepatic bile ducts",
             "04 Malignant neoplasms of lymphoid haematopoietic and related tissue",
             "08 Cervical Cancer",
             "09 All other cancers (C00-D48)"), "other_ca",
  cod == "11 Diabetes" , "t2dm", # although it includes T1DM as well. Bias should be small
  cod == "12 Ischaemic heart disease", "chd",
  cod == "13 Stroke (cerebrovascular)", "stroke",
  cod == "17 Dementia and Alzheimers", "dementia",
  default = "nonmodelled"
)]
deaths <- deaths[, .(deaths = sum(deaths)), keyby = .(year, agegroup, sex, dimd, cod)]
deaths <- dcast(deaths, ... ~ cod, value.var = "deaths", drop = TRUE, fill = 0L)

absorb_dt(deaths, pop)
rm(pop)
deaths[, sum(nonmodelled)/sum(pops), keyby = year]
