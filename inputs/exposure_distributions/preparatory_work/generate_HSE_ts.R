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



# preample ---------------------------------------------------
if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
  library(CKutils)
  }
dependencies(c("fst", "data.table"))

# Load single datasets -----------------------------------------
HSE_files <- sort(list.files("./secure_data/HSE_data", "\\.tab$", full.names = TRUE))
names(HSE_files) <- paste0("hse", 2003:2018)
HSE_files <- lapply(HSE_files, data.table::fread)
invisible(lapply(HSE_files, function(x) setnames(x, tolower(names(x)))))

HSE_files$hse2015 <- NULL # until proper harmonisation
HSE_files$hse2016 <- NULL # until proper harmonisation
HSE_files$hse2017 <- NULL # until proper harmonisation
HSE_files$hse2018 <- NULL # until proper harmonisation

# Statins
statin_px <- function(dt) {
  if ("medbi01" %in% names(dt)) {
    set(dt, NULL, "statin_px", 0L)
    dt[medbi01 == "21201" | medbi02 == "21201" | medbi03 == "21201" |
         medbi04 == "21201" | medbi05 == "21201" | medbi06 == "21201" |
         medbi07 == "21201" | medbi08 == "21201" | medbi09 == "21201" |
         medbi10 == "21201" | medbi11 == "21201" | medbi12 == "21201" |
         medbi13 == "21201" | medbi14 == "21201" | medbi15 == "21201" |
         medbi16 == "21201" | medbi17 == "21201" | medbi18 == "21201" |
         medbi19 == "21201" | medbi20 == "21201" | medbi21 == "21201" |
         medbi22 == "21201", statin_px := 1L ]
  }
}
invisible(lapply(HSE_files, statin_px))
HSE_files$hse2014[, table(statin_px)] # statin prescribed

statin_tkn <- function(dt) {
  if (("medbia1" %in% names(dt)) & ("medbi01" %in% names(dt))) {
    set(dt, NULL, "statin_tkn", 0L)
    dt[(medbi01 == "21201" & medbia1 == 1) | (medbi02 == "21201" & medbia2 == 1) |
       (medbi03 == "21201" & medbia3 == 1) | (medbi04 == "21201" & medbia4 == 1) |
         (medbi05 == "21201" & medbia5 == 1) | (medbi06 == "21201" & medbia6 == 1) |
         (medbi07 == "21201" & medbia7 == 1) | (medbi08 == "21201" & medbia8 == 1) |
         (medbi09 == "21201" & medbia9 == 1) | (medbi10 == "21201" & medbia10 == 1) |
         (medbi11 == "21201" & medbia11 == 1) | (medbi12 == "21201" & medbia12 == 1) |
         (medbi13 == "21201" & medbia13 == 1) | (medbi14 == "21201" & medbia14 == 1) |
         (medbi15 == "21201" & medbia15 == 1) | (medbi16 == "21201" & medbia16 == 1) |
         (medbi17 == "21201" & medbia17 == 1) | (medbi18 == "21201" & medbia18 == 1) |
         (medbi19 == "21201" & medbia19 == 1) | (medbi20 == "21201" & medbia20 == 1) |
         (medbi21 == "21201" & medbia21 == 1) | (medbi22 == "21201" & medbia22 == 1),
       statin_tkn := 1L ]
  }
}
invisible(lapply(HSE_files, statin_tkn))
HSE_files$hse2014[, table(statin_tkn)] # statin taken over the last 7 days

lapply(HSE_files, function(x) grep("diage", names(x), value = TRUE))
lapply(HSE_files, function(x) grep("diabtot", names(x), value = TRUE))

# 2015
# HSE_files$hse2015 <-
#   HSE_files$hse2015[, .(
#     age35g, sex, qimd, bmival, cholval13, omsysval, cigst1, startsmk, endsmokg,
#   numsmok, smokyrsg, cigdyal,  origin2, hdlval13, bpmedd2, diabete2, sha,
#   vegpor, frtpor15,  potass, creatin, wt_int,wt_nurse, wt_blood,
#   wt_urine, psu, cluster, glyhbval, diabtotr, alcohol = d7unitwg * 8,
#   totalwu = totalwu * 8 / 7, xpsmok, eqv5, topqual3, bp1, statina,
#   statin_px, statin_tkn, diage)]
# setnames(HSE_files$hse2015, c("age90", "origin3", "sodiumval", "cholval13",
#                               "hdlval12", "bpmedd2"),
#          c("age", "origin", "sodium", "cholval1", "hdlval1", "bpmedd"))
# HSE_files$hse2015[, `:=`(year = 15L, a30to06 = NA, sodiumval = NA)] # year 2000 = year 0
# HSE_files$hse2015[, sha := as.integer(sha)]
# replace_from_table(HSE_files$hse2015, "origin", c(1:2, 4:18), c(1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 7L, 6L, 9L, 9L, 9L))
#

# 2014
# code bp1 (diagnosed htn), validation to be applied in 2013
# HSE_files$hse2014[, bp2 := docbp]
# HSE_files$hse2014[docbp == -1L, bp2 := 2]
# HSE_files$hse2014[sex == 2L & othbp == 2L, bp2 := 2]
# HSE_files$hse2014[docbp == -9L | pregbp == -9L | othbp == -9L, bp2 := -9]
# HSE_files$hse2014[docbp == -8L | pregbp == -8L | othbp == -8L, bp2 := -8]
# HSE_files$hse2014[, table(bp1 == bp2)]

HSE_files$hse2014 <- HSE_files$hse2014[, .(
  age90, sex, qimd, bmival, cholval12, omsysval, cigst1, startsmk, endsmoke,
  numsmok, smokyrs, cigdyal,  origin3, hdlval12, bpmedd2, diabete2, sha,
  vegpor, frtpor, sodiumval, potass, creatin, wt_int,wt_nurse, wt_blood,
  wt_urine, psu, cluster, glyhbval, diabtotr, alcohol = d7unitwg * 8,
  totalwu = totalwu * 8 / 7, xpsmok, eqv5, topqual3, bp1, statina,
  statin_px, statin_tkn, diage)]
setnames(HSE_files$hse2014, c("age90", "origin3", "sodiumval", "cholval12",
                              "hdlval12", "bpmedd2"),
         c("age", "origin", "sodium", "cholval1", "hdlval1", "bpmedd"))
HSE_files$hse2014[, `:=`(year = 14L, a30to06 = NA)] # year 2000 = year 0
HSE_files$hse2014[, sha := as.integer(sha)]
replace_from_table(HSE_files$hse2014, "origin", c(1:2, 4:18), c(1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 7L, 6L, 9L, 9L, 9L))

# 2013
# code bp1 (diagnosed htn)
HSE_files$hse2013[, bp1 := docbp]
HSE_files$hse2013[docbp == -1L, bp1 := 2]
HSE_files$hse2013[sex == 2L & othbp == 2L, bp1 := 2]
HSE_files$hse2013[docbp == -9L | pregbp == -9L | othbp == -9L, bp1 := -9]
HSE_files$hse2013[docbp == -8L | pregbp == -8L | othbp == -8L, bp1 := -8]

HSE_files$hse2013 <- HSE_files$hse2013[, .(
  wt_int, wt_nurse, wt_blood, psu, cluster, age, sex, qimd,
  bmival, cholval12, omsysval, diabtotr, cigst1, startsmk, endsmoke, numsmok, bp1,
  smokyrs, cigdyal,  origin, hdlval12, bpmedd2, diabete2, sha, vegpor, frtpor,
  glyhbval, alcohol = d7unitwg * 8, totalwu = totalwu * 8 / 7, xpsmok, eqv5,
  topqual3, statina, statin_px, statin_tkn, diage)]
HSE_files$hse2013[, `:=`(year = 13L, a30to06 = NA, sodium = NA, potass = NA,
                         creatin = NA,  wt_urine = NA)]
setnames(HSE_files$hse2013, c("cholval12", "hdlval12", "bpmedd2"),
         c("cholval1", "hdlval1", "bpmedd"))
replace_from_table(HSE_files$hse2013, "origin", 1:18, c(1L, 1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 6L, 7L, 9L, 9L, 9L))

HSE_files$hse2013[, sha := as.integer(sha)]
HSE_files$hse2013[, length(unique(psu)), by = cluster][, table(V1)] # stratas can have many psus
HSE_files$hse2013[, length(unique(cluster)), by = psu][, table(V1)] # but here one psu has 2 stratas
HSE_files$hse2013[, length(unique(cluster)), by = psu][V1 == 2, ] #2131015
HSE_files$hse2013[psu == 2131015 , cluster]
HSE_files$hse2013[cluster == 213235, psu]
HSE_files$hse2013[psu == 2131015 & cluster == 213235, psu := 2131340]

# 2012
HSE_files$hse2012 <- HSE_files$hse2012[, .(
  wt_int, wt_nurse, wt_blood, psu, cluster, age, sex, qimd, bmival, cholval12,
  omsysval, diabtotr, cigst1, startsmk, endsmoke, numsmok, smokyrs, cigdyal,
  a30to06, sodiumval, potass, creatin, wt_urine, origin, hdlval12, bpmedd2, bp1,
  diabete2, sha, glyhbval, alcohol = d7unitwg * 8, totalwu = totalwu * 8 / 7,
  xpsmok, eqv5, topqual3, statina, statin_px, statin_tkn, diage)]
HSE_files$hse2012[, `:=`(year = 12L, vegpor = NA, frtpor = NA)]
setnames(HSE_files$hse2012, c("cholval12", "hdlval12", "bpmedd2", "sodiumval"),
         c("cholval1", "hdlval1", "bpmedd", "sodium"))
replace_from_table(HSE_files$hse2012, "origin", 1:18, c(1L, 1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 6L, 7L, 9L, 9L, 9L))

# 2011
HSE_files$hse2011 <- HSE_files$hse2011[, .(
  wt_int, wt_nurse, wt_blood, psu, cluster, age, sex, qimd,
  bmival, cholval1, omsysval, diabtotr, cigst1, startsmk, endsmoke, vegpor,
  frtpor, numsmok, smokyrs, cigdyal, origin, hdlval1, iregdef, bpmedd, bp1,
  diabete2, sha, glyhbval, alcohol = d7unitwg * 8, totalwu = totalwu * 8 / 7,
  xpsmok, eqv5, topqual3, statina, statin_px, statin_tkn, diage = ageinfo1)]
HSE_files$hse2011[, `:=`(year = 11L, a30to06 = NA, sodium = NA, potass = NA,
                         creatin = NA, wt_urine = 0)]
HSE_files$hse2011[, sha := as.integer(substr(sha, nchar(sha) - 2L, nchar(sha)))]
replace_from_table(HSE_files$hse2011, "origin", 1:18, c(1L, 1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 6L, 7L, 9L, 9L, 9L))

# 2010
HSE_files$hse2010[cholflag == 1L, `:=`(cholval1 = cholval1 + 0.1, hdlval1 = hdlval1 - 0.1)]
HSE_files$hse2010 <- HSE_files$hse2010[samptype == 1, .(
  wt_int, wt_nurse, wt_blood, psu, cluster, age, sex, imd2007, bmival, cholval1,
  omsysval, cigst1, startsmk, endsmoke, vegpor, frtpor, numsmok, bp1,
  smokyrs, cigdyal, sodival, potass, creatin, wt_urine, origin, hdlval1,
  kidfailgp, kiddiag, bpmedd, diabete2, sha, glyhbval, alcohol = d7unitwg * 8,
  xpsmok, eqv5, topqual3, statina, statin_px, statin_tkn, diage)]
setnames(HSE_files$hse2010, c("imd2007", "sodival"), c("qimd", "sodium"))
HSE_files$hse2010[is.na(wt_nurse), wt_nurse := 0]
HSE_files$hse2010[is.na(wt_blood), wt_blood := 0]
HSE_files$hse2010[diabete2 == 2, diabtotr := 1]
HSE_files$hse2010[diabete2 == 1 | glyhbval > 6.5, diabtotr := 2]
HSE_files$hse2010[, `:=`(year = 10L, a30to06 = NA)]
replace_from_table(HSE_files$hse2010, "sha", paste0("Q3", 0:9), paste0(1:10), "sha2")
HSE_files$hse2010[, sha := as.integer(sha2)]
HSE_files$hse2010[, sha2 := NULL]
replace_from_table(HSE_files$hse2010, "origin", 1:16, c(1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 5L, 6L, 7L, 9L, 8L, 9L))

# 2009
HSE_files$hse2009[, `:=` (cholval1 = cholval1 + 0.1, hdlval1 = hdlval1 - 0.1)]
HSE_files$hse2009 <- HSE_files$hse2009[samptype == 1, .(
  wt_int, wt_nurse, wt_blood, psu, cluster, age, sex, imd2007, bmival, cholval1,
  omsysval, cigst1, startsmk, endsmoke, vegpor, frtpor, numsmok, bp1,
  smokyrs, cigdyal, sodium, potass, creatin, wt_urine, origin, hdlval1, bpmedd,
  diabete2, sha, glyhbval, alcohol = d7unitwg * 8, xpsmok, eqv5, topqual3,
  statina, statin_px, statin_tkn, diage)]
setnames(HSE_files$hse2009, "imd2007" , "qimd")
HSE_files$hse2009[diabete2 == 2, diabtotr := 1]
HSE_files$hse2009[diabete2 == 1 | glyhbval > 6.5, diabtotr := 2]
HSE_files$hse2009[, `:=`(year = 9L, a30to06 = NA)]
replace_from_table(HSE_files$hse2009, "sha", paste0("Q3", 0:9), paste0(1:10), "sha2")
HSE_files$hse2009[, `:=` (sha = as.integer(sha2), sha2 = NULL)]
replace_from_table(HSE_files$hse2009, "origin", 1:16, c(1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 5L, 6L, 7L, 9L, 8L, 9L))

# 2008
HSE_files$hse2008[, `:=` (cholval1 = cholval1 + 0.1, hdlval1 = hdlval1 - 0.1)]
HSE_files$hse2008 <- HSE_files$hse2008[samptype == 1, .(
  wt_int, wt_nurse, wt_blood, psu, cluster, age, sex, qimd, bmival, cholval1,
  omsysval,  cigst1, startsmk, endsmoke, vegpor, frtpor, numsmok, smokyrs,
  cigdyal, a30to06, origin, hdlval1, bpmedd, sha, glyhbval, alcohol = d7unitwg * 8,
  xpsmok, eqv5, topqual3, statina, statin_px, statin_tkn, diage = NA)]
HSE_files$hse2008[, `:=`(year = 8L, diabtotr = NA, diabete2 = NA, sodium = NA, potass = NA,
                         creatin = NA, wt_urine = 0, bp1 = NA)]
replace_from_table(HSE_files$hse2008, "origin", 1:16, c(1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 5L, 6L, 7L, 9L, 8L, 9L))

# 2007
HSE_files$hse2007[ethinda  == 1,                  ethnicity := 1L] # white
HSE_files$hse2007[indcult1 == 1 | indcult2 == 1,  ethnicity := 2L] # indian
HSE_files$hse2007[indcult3 == 1,                  ethnicity := 3L] # pakistani
HSE_files$hse2007[indcult4 == 1,                  ethnicity := 4L] # bangladesh
HSE_files$hse2007[indcult5 == 1,                  ethnicity := 5L] # other asian
HSE_files$hse2007[blacult1 == 1,                  ethnicity := 6L] # black caribbean
HSE_files$hse2007[blacult2 == 1,                  ethnicity := 7L] # black african
HSE_files$hse2007[othcult1 == 1,                  ethnicity := 8L] # chinese
HSE_files$hse2007[ethinda > 0 & is.na(ethnicity), ethnicity := 9L] # other
HSE_files$hse2007 <- HSE_files$hse2007[samptype == 1, .(
  wt_int, wt_nurse,  area, cluster, age, sex, imd2007, bmival, omsysval, cigst1,
  startsmk, endsmoke, vegpor, frtpor, numsmok, smokyrs, cigdyal, sodium,
  potass, creatin, bpmedd, newsha, ethnicity, alcohol = d7unitwg * 8,
  xpsmok = xpsm, eqv5, topqual3, statina, statin_px, statin_tkn, diage = NA)]
setnames(HSE_files$hse2007, c("imd2007", "area", "newsha"), c("qimd", "psu", "sha"))
HSE_files$hse2007[, `:=`(year = 7L, diabtotr = NA, diabete2 = NA, wt_blood = 1,
                         cholval1 = NA, a30to06 = NA, wt_urine = wt_nurse, bp1 = NA,
                         glyhbval = NA)]

# 2006
HSE_files$hse2006[ethinda  == 1,                  ethnicity := 1L] # white
HSE_files$hse2006[indcult1 == 1 | indcult2 == 1,  ethnicity := 2L] # indian
HSE_files$hse2006[indcult3 == 1,                  ethnicity := 3L] # pakistani
HSE_files$hse2006[indcult4 == 1,                  ethnicity := 4L] # bangladesh
HSE_files$hse2006[indcult5 == 1,                  ethnicity := 5L] # other asian
HSE_files$hse2006[blacult1 == 1,                  ethnicity := 6L] # black caribbean
HSE_files$hse2006[blacult2 == 1,                  ethnicity := 7L] # black african
HSE_files$hse2006[othcult1 == 1,                  ethnicity := 8L] # chinese
HSE_files$hse2006[ethinda > 0 & is.na(ethnicity), ethnicity := 9L] # other
HSE_files$hse2006[, `:=` (cholval1 = cholval1 + 0.1, hdlval1 = hdlval1 - 0.1)]
HSE_files$hse2006 <- HSE_files$hse2006[samptype != 3, .(
  wt_int, wt_nurse, wt_blood, psu, cluster, age, sex, imd2004, bmival, cholval1,
  omsysval, cigst1, startsmk, endsmoke, vegpor, frtpor, numsmok, bp1,
  smokyrs, cigdyal, a30to06, sodium, potass, creatin, famcvd, hdlval1, bpmedd,
  diabete2, newsha, glyhbval, ethnicity, alcohol = drevunit * 8, xpsmok = xpsm,
  eqv5, topqual3, statina, statin_px, statin_tkn, diage = NA)]
setnames(HSE_files$hse2006, c("imd2004", "newsha"), c("qimd", "sha"))
HSE_files$hse2006[diabete2 == 2, diabtotr := 1]
HSE_files$hse2006[diabete2 == 1 | glyhbval > 6.5, diabtotr := 2]
HSE_files$hse2006[, `:=`(year = 6L, wt_urine = wt_nurse)]

# 2005
HSE_files$hse2005[ethinda  == 1,                  ethnicity := 1L] # white
HSE_files$hse2005[indcult1 == 1 | indcult2 == 1,  ethnicity := 2L] # indian
HSE_files$hse2005[indcult3 == 1,                  ethnicity := 3L] # pakistani
HSE_files$hse2005[indcult4 == 1,                  ethnicity := 4L] # bangladesh
HSE_files$hse2005[indcult5 == 1,                  ethnicity := 5L] # other asian
HSE_files$hse2005[blacult1 == 1,                  ethnicity := 6L] # black caribbean
HSE_files$hse2005[blacult2 == 1,                  ethnicity := 7L] # black african
HSE_files$hse2005[othcult1 == 1,                  ethnicity := 8L] # chinese
HSE_files$hse2005[ethinda > 0 & is.na(ethnicity), ethnicity := 9L] # other
HSE_files$hse2005[, `:=` (cholval1 = cholval1 + 0.1, hdlval1 = hdlval1 - 0.1)]
HSE_files$hse2005 <- HSE_files$hse2005[samptype == 1, .(
  wt_int, wt_nurse, wt_bldel, area, cluster, age, sex,  imd2004, bmival, cholval1,
  omsysval, cigst1, startsmk, endsmoke, vegpor, frtpor, numsmok, bp1,
  smokyrs, cigdyal, sodium, potass, creatin, hdlval1, bpmedd, diabete2, newsha,
  glyhbval, ethnicity, alcohol = d7unit * 8, xpsmok = xpsm,
  eqv5, topqual3, statina, statin_px, statin_tkn, diage = NA)]
setnames(HSE_files$hse2005, c("imd2004", "area", "wt_bldel", "newsha"),
         c("qimd", "psu", "wt_blood", "sha"))
HSE_files$hse2005[diabete2 == 2, diabtotr := 1]
HSE_files$hse2005[diabete2 == 1 | glyhbval > 6.5, diabtotr := 2]
HSE_files$hse2005[, `:=`(year = 5L, a30to06 = NA, wt_urine = wt_nurse)]

# 2004
replace_from_table(HSE_files$hse2004, "dmethn04",
                   1:9, c(6L, 7L, 2L, 3L, 4L, 8L, 1L, 1L, 9L), "ethnicity")
HSE_files$hse2004[indclt6 == 1 & dmethn04 == 9, ethnicity := 5L] # other asian
HSE_files$hse2004[, `:=` (cholval1 = cholval1 + 0.1, hdlval1 = hdlval1 - 0.1)]
HSE_files$hse2004[samptype == 7 & eqvinc >= 0, eqv5 :=
                    as.integer(cut(eqvinc, quantile(eqvinc, c(0, 0.2, 0.4, 0.6, 0.8, 1)), 1:5))]
HSE_files$hse2004 <- HSE_files$hse2004[samptype == 7, .(
  wt_int, area, cluster, age, sex, imd2004, bmival, cholval1, bp1,
  omsysval, cigst1, startsmk, endsmoke, vegpor, frtpor, numsmok,
  smokyrs, cigdyal, adtot30, sodium, potass, creatin, hdlval1, bpmedd, diabete2,
  sha, glyhbval, ethnicity, alcohol = d7unit * 8, xpsmok = xpsm,
  eqv5, topqual3, diage = NA)]
setnames(HSE_files$hse2004, c("imd2004", "area", "adtot30"), c("qimd", "psu", "a30to06"))
replace_from_table(HSE_files$hse2004, "sha", paste0("Q", sprintf("%02.0f", 1:28)),
                   as.character(c(6, 6, 6, 7, 7, 7, 7, 7, 1, 1, 3, 3, 2, 2, 2, 9, 9, 8,
                                  8, 10, 10, 10, 3, 4, 4, 5, 5, 5)), "sha2")
HSE_files$hse2004[, `:=` (sha = as.integer(sha2), sha2 = NULL)]
HSE_files$hse2004[diabete2 == 2, diabtotr := 1]
HSE_files$hse2004[diabete2 == 1 | glyhbval > 6.5, diabtotr := 2]
HSE_files$hse2004[, `:=`(year = 4L, wt_urine = wt_int, wt_blood = wt_int, wt_nurse = wt_int)]

# 2003
HSE_files$hse2003[ethnici  == 1,                  ethnicity := 1L] # white
HSE_files$hse2003[indcult1 == 1 | indcult2 == 1,  ethnicity := 2L] # indian
HSE_files$hse2003[indcult3 == 1,                  ethnicity := 3L] # pakistani
HSE_files$hse2003[indcult4 == 1,                  ethnicity := 4L] # bangladesh
HSE_files$hse2003[indcult5 == 1,                  ethnicity := 5L] # other asian
HSE_files$hse2003[blacult1 == 1,                  ethnicity := 6L] # black caribbean
HSE_files$hse2003[blacult2 == 1,                  ethnicity := 7L] # black african
HSE_files$hse2003[othcult1 == 1,                  ethnicity := 8L] # chinese
HSE_files$hse2003[ethnici > 0 & is.na(ethnicity), ethnicity := 9L] # other
HSE_files$hse2003[, `:=` (cholval1 = cholval1 + 0.1, hdlval1 = hdlval1 - 0.1)]
HSE_files$hse2003 <- HSE_files$hse2003[, .(
  int_wt, nurse_wt, blood_wt, area, cluster, age, sex, imd2004, bmival,
  cholval1, omsysval, cigst1, startsmk, endsmoke, vegpor, frtpor,bp1,
  numsmok, smokyrs, cigdyal, adtot30, sodium, potass, creatin, hdlval1,
  bpmedd, diabete2, sha, glyhbval, ethnicity, alcohol = d7unit * 8, xpsmok = xpsm,
  eqv5, topqual3, diage = NA)]
setnames(HSE_files$hse2003,
         c("imd2004", "area", "int_wt", "blood_wt", "nurse_wt","adtot30"),
         c("qimd", "psu", "wt_int", "wt_blood", "wt_nurse", "a30to06"))
HSE_files$hse2003[diabete2 == 2, diabtotr := 1]
HSE_files$hse2003[diabete2 == 1 | glyhbval > 6.5, diabtotr := 2]
replace_from_table(HSE_files$hse2003, "sha", paste0("Q", sprintf("%02.0f", 1:28)),
                   as.character(c(6, 6, 6, 7, 7, 7, 7, 7, 1, 1, 3, 3, 2, 2, 2, 9, 9, 8,
                                  8, 10, 10, 10, 3, 4, 4, 5, 5, 5)), "sha2")
HSE_files$hse2003[, `:=` (sha = as.integer(sha2), sha2 = NULL)]
HSE_files$hse2003[, `:=`(year = 3L, wt_urine = wt_nurse)]


HSE_ts <- rbindlist(HSE_files,use.names = TRUE, fill = TRUE)
HSE_ts[bmival    < 0, bmival    := NA]
HSE_ts[cholval1  < 0, cholval1  := NA]
HSE_ts[omsysval  < 0, omsysval  := NA]
HSE_ts[cigst1    < 0, cigst1    := NA]
HSE_ts[endsmoke  < 0, endsmoke  := NA]
HSE_ts[numsmok   < 0, numsmok   := NA]
HSE_ts[smokyrs   < 0, smokyrs   := NA]
HSE_ts[cigdyal   < 0, cigdyal   := NA]
HSE_ts[vegpor    < 0, vegpor    := NA]
HSE_ts[frtpor    < 0, frtpor    := NA]
HSE_ts[a30to06   < 0, a30to06   := NA]
HSE_ts[sodium    < 0, sodium    := NA]
HSE_ts[potass    < 0, potass    := NA]
HSE_ts[creatin   < 0, creatin   := NA]
HSE_ts[hdlval1   < 0, hdlval1   := NA]
HSE_ts[bpmedd    < 0, bpmedd    := NA]
HSE_ts[bp1       < 0, bp1       := NA]
HSE_ts[statina   < 0, statina   := NA]
HSE_ts[glyhbval  < 0, glyhbval  := NA]
HSE_ts[diabtotr  < 0, diabtotr  := NA]
HSE_ts[diabete2  < 0, diabete2  := NA]
HSE_ts[famcvd    < 0, famcvd    := NA]
HSE_ts[origin    < 0, origin    := NA]
HSE_ts[kidfailgp < 0, kidfailgp := NA]
HSE_ts[kiddiag   < 0, kiddiag   := NA]
HSE_ts[iregdef   < 0, iregdef   := NA]
HSE_ts[alcohol   < 0, alcohol   := NA]
HSE_ts[totalwu   < 0, totalwu   := NA]
HSE_ts[xpsmok   < 0, xpsmok   := NA] # 97 = more than 97 hours (for some years)
HSE_ts[eqv5      < 0, eqv5      := NA]
HSE_ts[topqual3  < 0, topqual3  := NA]
HSE_ts[diage     < 0, diage     := NA]
HSE_ts[is.na(ethnicity), ethnicity   := origin]
HSE_ts[(startsmk   <= 0) | (startsmk   == 97), startsmk   := NA] # 97 = never smoked regularly
HSE_ts[, dm_dur := age - diage]
HSE_ts[, diage  := NULL]


HSE_ts[, `:=` (
  sex = factor(sex, levels = 1:2, labels = c("men", "women")),
  qimd = factor(
    6L - qimd,
    levels = 1:5,
    labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
  ),
  cigst1 = factor(cigst1),
  bpmedd = factor(bpmedd),
  diabete2 = factor(diabete2, levels = 2:1, labels = 0:1),
  sha = factor(
    sha,
    levels = 1:10,
    labels = c(
      "North East",
      "North West",
      "Yorkshire and the Humber",
      "East Midlands",
      "West Midlands",
      "East of England",
      "London",
      "South East Coast",
      "South Central",
      "South West"
    )
  ),
  diabtotr = factor(diabtotr, levels = 1:2, labels = 0:1),
  bp1 = factor(bp1, levels = 2:1, labels = 0:1),
  statina = factor(statina, levels = 2:1, labels = 0:1),
  eqv5 = factor(
    eqv5,
    levels = 1:5,
    labels = c(
      "1 Highest",
      "2", "3", "4",
      "5 Lowest"
    )),
  topqual3 = factor(
    topqual3,
    levels = 1:7,
    labels = c(
      "NVQ4/NVQ5/Degree or equiv",
      "Higher ed below degree",
      "NVQ3/GCE A Level equiv",
      "NVQ2/GCE O Level equiv",
      "NVQ1/CSE other grade equiv",
      "Foreign/other",
      "No qualification"
    )),
  a30to06 = a30to06 / 4,
  origin = NULL,
  famcvd = factor(famcvd, levels = 2:1, labels = 0:1),
  ethnicity = factor(
    ethnicity,
    levels = 1:9,
    labels = c(
      "white",
      "indian",
      "pakistani",
      "bangladeshi",
      "other asian",
      "black caribbean",
      "black african",
      "chinese",
      "other"
    )
  ),
  kidfailgp = factor(kidfailgp),
  kiddiag = factor(kiddiag, levels = 2:1, labels = 0:1),
  iregdef = factor(iregdef, levels = 2:1, labels = 0:1),
  xpsmok = as.factor(as.integer(xpsmok > 0))
)]


HSE_ts[(omsysval >=140 | bpmedd == "1") & bp1 == "1", htn_dgn := 1L] # diagnosed hypertensives
HSE_ts[(omsysval >=140 | bpmedd == "1") & bp1 == "0", htn_dgn := 0L] # undiagnosed hypertensives
HSE_ts[, htn_dgn := factor(htn_dgn)]
# HSE_ts[, totalwu := as.integer(round(totalwu))]
# HSE_ts[, alcohol := as.integer(round(alcohol))]
# HSE_ts[, endsmoke := as.integer(round(endsmoke))]
# HSE_ts[, numsmok := as.integer(round(numsmok))]
# HSE_ts[, smokyrs := as.integer(round(smokyrs))]
# HSE_ts[, cigdyal := as.integer(round(cigdyal))]
# HSE_ts[, a30to06 := as.integer(round(a30to06))]
# HSE_ts[, frtpor := as.integer(round(frtpor))]
# HSE_ts[, vegpor := as.integer(round(vegpor))]

HSE_ts <- HSE_ts[between(age, 15, 100), ]


to_agegrp(HSE_ts, 20L, 85L, "age", "agegrp20", to_factor = TRUE)
to_agegrp(HSE_ts, 10L, 85L, "age", "agegrp10", to_factor = TRUE)
to_agegrp(HSE_ts,  5L, 85L, "age", "agegrp5" , to_factor = TRUE)
setnames(HSE_ts,
         c("bmival", "cholval1", "omsysval", "cigst1", "startsmk",
           "endsmoke", "numsmok", "smokyrs", "cigdyal", "a30to06",
           "hdlval1", "diabtotr", "diabete2", "iregdef", "kidfailgp", "glyhbval",
           "bpmedd", "xpsmok", "eqv5", "topqual3", "kiddiag", "statina"),
         c("bmi", "tchol", "sbp", "smok_status", "smok_init_age",
           "smok_quit_yrs", "smok_cig_ex", "smok_dur_ex", "smok_cig_curr",
           "active_days", "hdl", "dm", "dm_dgn", "af", "ckd", "hba1c", "bp_med", "ets",
           "income", "education", "ckd_dgn", "statin_otc"))
setcolorder(HSE_ts, sort(copy(names(HSE_ts))))
# HSE_ts[, qimd := ordered(qimd)]
# summary(HSE_ts)

# Standartise weights to the number of participants every year
HSE_ts[, wt_int   := sum(wt_int > 0)   * wt_int/sum(wt_int),     keyby = year]
HSE_ts[, wt_nurse := sum(wt_nurse > 0) * wt_nurse/sum(wt_nurse), keyby = year]
HSE_ts[, wt_blood := sum(wt_blood > 0) * wt_blood/sum(wt_blood), keyby = year]
# HSE_ts[, .(sum(wt_blood), .N), keyby = year]

rm(HSE_files)
write_fst(HSE_ts, "./preparatory_work/HSE_ts.fst")

# HSE_ts[, lapply(.SD, function(x) sum(is.na(x))), keyby = year]
#
# tt <- missForest::missForest(HSE_ts)

# EQV5: (D) Equivalised Income Quintiles
# 5 Highest Quintile (>£44,318)
# 4 Second highest Quintile (>£27,637 <=£44,318)
# 3 Middle Quintile (>£19,180 <=£27,637)
# 2 Second lowest Quintile (>£12,118 <=£19,180)
# 1 Lowest Quintile (<=£12,118)

# TOPQUAL3: (D) Highest Educational Qualification
# 1 NVQ4/NVQ5/Degree or equiv
# 2 Higher ed below degree
# 3 NVQ3/GCE A Level equiv
# 4 NVQ2/GCE O Level equiv
# 5 NVQ1/CSE other grade equiv
# 6 Foreign/other
# 7 No qualification

# CKD: Chronic kidney disease stage
# 1 Normal: eGFR 60+ ml/min/1.73m2 and normal albuminuria
# 2 Stage 1: eGFR 90+ ml/min/1.73m2 and micro- or macro-albuminuria
# 3 Stage 2: eGFR 60-89 ml/min/1.73m2 and micro- or macro-albuminuria
# 4 Stage 3a/3b: eGFR 30-59 ml/min/1.73m2 regardless of albuminuria
# 5 Stage 4/5: eGFR less than 30 ml/min/1.73m2 regardless of albuminuria

# totalwu is the grams of ethanol per day, based on average weekly consumption
# alcohol is the grams of ethanol per day, based on the day with highest alcohol
# intake over a week of observation

# NOTE method to prescribed medication changed since 2012. Since then they record
# Actual use rather than prescription for all categories of medicine.
