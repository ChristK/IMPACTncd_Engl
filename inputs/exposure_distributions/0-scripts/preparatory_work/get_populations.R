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

library(CKutils)
dt <- fread("/home/ckyprid/pCloudDriveMy_ModelsIMPACTncd_Engl data/ward pops by age gender ethnicity from census2011.csv")
dt <- melt(dt, 1:3, value.name = "population_size")
dt[, ward_code := substr(as.character(variable), 1, 6)]
dt[, ward_name := substr(as.character(variable), 10, nchar(as.character(variable)))]
dt[, variable := NULL]
agegrp_name(0, 90, 5, FALSE)
replace_from_table(dt, "Age", dt[, sort(unique(Age))], agegrp_name(0, 90, 5, FALSE))

# From NOMIS ONS Census 2011 population by ward age sex ethnicity
dt <- fread("~My_ModelsIMPACTncd_Engl data/Data_AGE_ETHGEW_SEX_UNIT.csv")
dt <- melt(dt, 1:2, value.name = "population_size")
dt[variable %like% "Male", sex := "men"]
dt[variable %like% "Female", sex := "women"]
dt[, agegrp :=  rapply(regmatches(variable,regexec("Age : (.*?) years old", variable)), `[`,
                       classes = "ANY", deflt = NULL, how = "unlist", 2)]
dt[, sort(unique(agegrp))]
dt <- dt[!agegrp %in% c("0 to 15", "50 to 64")]
dt[agegrp %in% c("5 to 7", "8 to 9"), agegrp := "5 to 9"]
agegrp_name(0, 90, 5, FALSE)
replace_from_table(dt, "Age", c( "0 to 4", "5 to 9", "10 to 14", "20 to 24" "25 to 29" "30 to 34" "35 to 39" "40 to 44" "45 to 49"
                                 "50 to 54" "55 to 59" "60 to 64" "65 to 69" "70 to 74" "75 to 79" "80 to 84"
                                 "85 to 89" "90" ), agegrp_name(0, 90, 5, FALSE))



dt <- dt[, sum(population_size), by = .(sex, agegrp, ethnicity, GEO_CODE, GEO_LABEL)]
dt[, unique(variable)]
dt[variable %like% "16", unique(variable)]

dt[, ward_code := substr(as.character(variable), 1, 6)]
dt[, ward_name := substr(as.character(variable), 10, nchar(as.character(variable)))]
dt[, variable := NULL]
agegrp_name(0, 90, 5, FALSE)
replace_from_table(dt, "Age", dt[, sort(unique(Age))], agegrp_name(0, 90, 5, FALSE))

dt2 = dt[sample(.N, 100), ]
dt2[, agegrp :=  rapply(regmatches(variable,regexec("Age : (.*?) years old", variable)), `[`,
             classes = "ANY", deflt = NULL, how = "list", 2)]



test = dt[sample(.N, 10), variable]
pattern="Age : (.*?) years old"
result <- regmatches(test,regexec(pattern,test))
lapply(result, `[[`, 2)
result[[1]][2]
