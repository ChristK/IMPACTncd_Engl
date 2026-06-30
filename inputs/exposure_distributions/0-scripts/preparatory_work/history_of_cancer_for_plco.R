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




# History of cancer: from Maddams J, et al. Projections of cancer
# prevalence in the United  Kingdom, 2010-2040.
# Br J Cancer 2012; 107: 1195-1202 (table 5, scenario 1)

# This is needed for PLCO lung cancer calculation

dt <- CJ(age = 20:89,
  sex = c("men", "women"),
  year = 3:50)



dt[age < 45             &
    sex == "men"  & year < 10L, mu := 381 / 1e5]
dt[between(age, 45, 64) &
    sex == "men"  & year < 10L, mu := 2669 / 1e5]
dt[age > 64             &
    sex == "men"  & year < 10L, mu := 12656 / 1e5]
dt[age < 45             &
    sex == "women"  & year < 10L, mu := 525 / 1e5]
dt[between(age, 45, 64) &
    sex == "women"  & year < 10L, mu := 4952 / 1e5]
dt[age > 64             &
    sex == "women"  & year < 10L, mu := 12801 / 1e5]
dt[age < 45             &
    sex == "men" & between(year, 10L, 19L), mu := 391 / 1e5]
dt[between(age, 45, 64) &
    sex == "men" & between(year, 10L, 19L), mu := 3037 / 1e5]
dt[age > 64             &
    sex == "men" & between(year, 10L, 19L), mu := 15558 / 1e5]
dt[age < 45             &
    sex == "women" & between(year, 10L, 19L), mu := 567 / 1e5]
dt[between(age, 45, 64) &
    sex == "women" & between(year, 10L, 19L), mu := 5914 / 1e5]
dt[age > 64             &
    sex == "women" & between(year, 10L, 19L), mu := 15909 / 1e5]
dt[age < 45             &
    sex == "men" & between(year, 20L, 29L), mu := 424 / 1e5]
dt[between(age, 45, 64) &
    sex == "men" & between(year, 20L, 29L), mu := 3459 / 1e5]
dt[age > 64             &
    sex == "men" & between(year, 20L, 29L), mu := 18698 / 1e5]
dt[age < 45             &
    sex == "women" & between(year, 20L, 29L), mu := 655 / 1e5]
dt[between(age, 45, 64) &
    sex == "women" & between(year, 20L, 29L), mu := 7299 / 1e5]
dt[age > 64             &
    sex == "women" & between(year, 20L, 29L), mu := 19261 / 1e5]
dt[age < 45             &
    sex == "men" & year > 29L, mu := 431 / 1e5]
dt[between(age, 45, 64) &
    sex == "men" & year > 29L, mu := 3632 / 1e5]
dt[age > 64             &
    sex == "men" & year > 29L, mu := 23301 / 1e5]
dt[age < 45             &
    sex == "women" & year > 29L, mu := 690 / 1e5]
dt[between(age, 45, 64) &
    sex == "women" & year > 29L, mu := 8419 / 1e5]
dt[age > 64             &
    sex == "women" & year > 29L, mu := 24852 / 1e5]

write_fst(dt, "./inputs/exposure_distributions/history_of_cancer.fst", 100L)
