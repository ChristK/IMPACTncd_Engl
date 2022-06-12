## adapted the generate_HSE_ts.R to clean the 2015 - 2018 data
## Mar 2022
## documents on the dataset:
## 2016: smoking, social care, adult physical activity
## new topics for 2017: chronic pain, VVD - long module, social care provision, EQ5D-5L, IPAQ, 
## interview content removed from 2017: cardiovascular disease

library(fst)
library(data.table)

# Load single datasets -----------------------------------------
HSE_files <- sort(list.files("./data", "\\.tab$", full.names = TRUE))
names(HSE_files) <- paste0("hse", 2013:2018)
HSE_files <- lapply(HSE_files, data.table::fread)
invisible(lapply(HSE_files, function(x) setnames(x, tolower(names(x)))))

original = read_fst("HSE_ts.fst", as.data.table = TRUE)
# 
# # Statins
# # TODO: don't have below
# statin_px <- function(dt) {
#   if ("medbi01" %in% names(dt)) {
#     set(dt, NULL, "statin_px", 0L)
#     dt[medbi01 == "21201" | medbi02 == "21201" | medbi03 == "21201" |
#          medbi04 == "21201" | medbi05 == "21201" | medbi06 == "21201" |
#          medbi07 == "21201" | medbi08 == "21201" | medbi09 == "21201" |
#          medbi10 == "21201" | medbi11 == "21201" | medbi12 == "21201" |
#          medbi13 == "21201" | medbi14 == "21201" | medbi15 == "21201" |
#          medbi16 == "21201" | medbi17 == "21201" | medbi18 == "21201" |
#          medbi19 == "21201" | medbi20 == "21201" | medbi21 == "21201" |
#          medbi22 == "21201", statin_px := 1L ]
#   }
# }
# invisible(lapply(HSE_files, statin_px))
# HSE_files$hse2018[, table(statin_px)] # statin prescribed
# # TODO: can't run: this name not in 2018 file
# 
# statin_tkn <- function(dt) {
#   if (("medbia1" %in% names(dt)) & ("medbi01" %in% names(dt))) {
#     set(dt, NULL, "statin_tkn", 0L)
#     dt[(medbi01 == "21201" & medbia1 == 1) | (medbi02 == "21201" & medbia2 == 1) |
#        (medbi03 == "21201" & medbia3 == 1) | (medbi04 == "21201" & medbia4 == 1) |
#          (medbi05 == "21201" & medbia5 == 1) | (medbi06 == "21201" & medbia6 == 1) |
#          (medbi07 == "21201" & medbia7 == 1) | (medbi08 == "21201" & medbia8 == 1) |
#          (medbi09 == "21201" & medbia9 == 1) | (medbi10 == "21201" & medbia10 == 1) |
#          (medbi11 == "21201" & medbia11 == 1) | (medbi12 == "21201" & medbia12 == 1) |
#          (medbi13 == "21201" & medbia13 == 1) | (medbi14 == "21201" & medbia14 == 1) |
#          (medbi15 == "21201" & medbia15 == 1) | (medbi16 == "21201" & medbia16 == 1) |
#          (medbi17 == "21201" & medbia17 == 1) | (medbi18 == "21201" & medbia18 == 1) |
#          (medbi19 == "21201" & medbia19 == 1) | (medbi20 == "21201" & medbia20 == 1) |
#          (medbi21 == "21201" & medbia21 == 1) | (medbi22 == "21201" & medbia22 == 1),
#        statin_tkn := 1L ]
#   }
# }
# invisible(lapply(HSE_files, statin_tkn))
# HSE_files$hse2018[, table(statin_tkn)] # statin taken over the last 7 days


# 2015
# code bp1 (diagnosed htn), validation to be applied in 2013
# HSE_files$hse2014[, bp2 := docbp]
# HSE_files$hse2014[docbp == -1L, bp2 := 2]
# HSE_files$hse2014[sex == 2L & othbp == 2L, bp2 := 2]
# HSE_files$hse2014[docbp == -9L | pregbp == -9L | othbp == -9L, bp2 := -9]
# HSE_files$hse2014[docbp == -8L | pregbp == -8L | othbp == -8L, bp2 := -8]
# # HSE_files$hse2014[, table(bp1 == bp2)]
# 
# lapply(HSE_files, function(x) grep("diage", names(x), value = TRUE))
# lapply(HSE_files, function(x) grep("diabtot", names(x), value = TRUE))

HSE_files$hse2015 <- HSE_files$hse2015[, .(
  seriala, age35g, sex, qimd, bmival, 
  #cholval12, # revised to cholval3 
  omsysval, cigst1, startsmk, 
  #endsmoke, # revised to ENDSMOKG # measured as grouped
  numsmok, 
  #smokyrs, # revised to SMOKYRSG # measured as grouped
  cigdyal,  
  #origin3, #revised to origin2 (grouped)
  #hdlval12, #revised to halval13
  bpmedd2, diabete2, sha,
  vegpor,
  #frtpor, # revised to porfrt
  #sodiumval, # not exist
  #potass, # can't find
  #creatin, # can't find
  wt_int,wt_nurse, wt_blood,
  #wt_urine, # can't find
  psu, #revised to psu_scr 
  cluster, 
  #glyhbval,# revised to GLYHBVALA
  diabtotr, alcohol = d7unitwg * 8,
  totalwu = totalwu * 8 / 7, xpsmok, eqv5, topqual3, bp1, 
  #statina, # can't find
  #statin_px, # can't find
  #statin_tkn, # can't find
  diage, # can only find diageg
  frtpor15, # revised from frtpor
  cholval13,
  hdlval13,
  glyhbvala,
  origin2, # grouped category
  smokyrsg, # measured as grouped
  endsmokg # measured as grouped
  )] # remove a few

setnames(HSE_files$hse2015, c( 
                              "origin2", 
                              #"sodiumval", 
                              "cholval13",
                              "hdlval13", 
                              "bpmedd2",
                              "frtpor15", "glyhbvala", "seriala"),
         c( "origin", 
           #"sodium", 
           "cholval1", 
           "hdlval1", 
           "bpmedd",
           "frtpor", "glyhbval", "id"))
HSE_files$hse2015[, `:=`(year = 15L, a30to06 = NA,
                         sodium = NA, # not exist
                         potass = NA, # can't find
                         creatin = NA, # can't find
                         wt_urine = NA, # can't find
                         statina = NA, # can't find
                         statin_px = NA, # can't find
                         statin_tkn = NA, # can't find
                         kiddiag = NA,
                         famcvd = NA
                         )] # year 2000 = year 0 # TODO: what is the a30to06
HSE_files$hse2015[, sha := as.integer(sha)]
#replace_from_table(HSE_files$hse2015, "origin", 
#                   c(1:2, 4:18), c(1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 7L, 6L, 9L, 9L, 9L))
# TODO: how to recode origin? unique(origin):  1  2  4  3  5 -9 -8

# 2016
# code bp1 (diagnosed htn)
# TODO: A30HS06: Number of days heavy housework 30 mins+' ; A30MA06: Number of days heavy manual 30 mins+

HSE_files$hse2016 <- HSE_files$hse2016[, .(
  seriala, wt_int, wt_nurse, wt_blood, psu, cluster, 
  #age, # revise age16g5
  sex, 
  qimd,
  bmival, 
  # cholval12,# revise cholval13 
  omsysval, diabtotr, cigst1, startsmk, 
  # endsmoke, # revise to endsmokg
  # numsmok, #  Can't find
  bp1,
  # smokyrs, # revise to smokyrsg
  cigdyal,  
  # origin, # revise to origin2 
  # hdlval12, # revise to hdlval13 
  bpmedd2, diabete2, sha, vegpor, 
  # frtpor,# revise porfrt
  # glyhbval, # revise to glyhbvala
  alcohol = d7unitwg * 8, totalwu = totalwu * 8 / 7, xpsmok, eqv5,
  topqual3, 
  # statina, # can't find
  # statin_px, # can't find
  # statin_tkn, # can't find
  diage,
  age35g,
  frtpor15,
  cholval13,
  hdlval13,
  glyhbvala,
  origin2, # grouped category
  smokyrsg, # measured as grouped
  endsmokg, # measured as grouped
  kidfail16s, # revised for kidfailgp
  cre,
  uk1,
  una1, 
  a30to06
  )]
HSE_files$hse2016[, `:=`(year = 16L,
                        wt_urine = NA,
                         statina = NA, # can't find
                         statin_px = NA, # can't find
                         statin_tkn = NA, # can't find
                         numsmok = NA,
                         kiddiag = NA,
                         famcvd = NA
                         )]
setnames(HSE_files$hse2016, c("cholval13", "hdlval13", "bpmedd2", "origin2", "frtpor15", "glyhbvala", "kidfail16s", "cre", "uk1", "una1", "seriala"),
         c("cholval1", "hdlval1", "bpmedd",  "origin", "frtpor",  "glyhbval",  "kidfailgp", "creatin", "potass", "sodium", "id"))
# what is below
#replace_from_table(HSE_files$hse2016, "origin", 1:18, c(1L, 1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 6L, 7L, 9L, 9L, 9L))
HSE_files$hse2016[, sha := as.integer(sha)]

# 2017
HSE_files$hse2017 <- HSE_files$hse2017[, .(
  seriala, wt_int, wt_nurse, wt_blood, 
  psu,cluster,
  # age, # revise to age16g5
  sex, qimd,
  bmival, 
  # cholval12,revise to cholval13
  omsysval,
  diabtotr, 
  cigst1, 
  startsmk, 
  # endsmoke, # endsmokg
  # numsmok, # can't find
  # smokyrs, # smokyrsg
  cigdyal,
  #a30to06, # not exist
  # sodiumval, # don't have
  # potass, #don't have
  # creatin, #don't have
  # wt_urine, # don't have
  # origin, # origin 2
  # hdlval12, # revised to hdlval13
  bpmedd2, bp1, diabete2, sha, 
  # glyhbval, # revised to glyhbvala
  alcohol = d7unitwg * 8, totalwu = totalwu * 8 / 7,
  xpsmok, eqv5, topqual3, 
  # statina, # can't find
  # statin_px, # can't find
  # statin_tkn, # can't find
  # diage, # can't find
  age35g,
  frtpor15,
  vegpor,
  cholval13,
  hdlval13,
  glyhbvala,
  origin2, # grouped category
  smokyrsg, # measured as grouped
  endsmokg, # measured as grouped
  iregdef
  )]

HSE_files$hse2017[, `:=`(year = 17L, a30to06 = NA, 
                         sodium = NA, potass = NA,
                         creatin = NA,  wt_urine = NA,
                         statina = NA, # can't find
                         statin_px = NA, # can't find
                         statin_tkn = NA, # can't find
                         numsmok = NA,
                         diage = NA,
                         kiddiag = NA,
                         famcvd = NA
                         )] # can only find diageg

setnames(HSE_files$hse2017, c("cholval13", "hdlval13", "bpmedd2",  "origin2", "frtpor15", "glyhbvala", "seriala"),
         c("cholval1", "hdlval1", "bpmedd", "origin", "frtpor", "glyhbval", "id"))

HSE_files$hse2017[, sha := as.integer(sha)]

#replace_from_table(HSE_files$hse2017, "origin", 1:18, c(1L, 1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 6L, 7L, 9L, 9L, 9L))

# 2018
HSE_files$hse2018 <- HSE_files$hse2018[, .(
  seriala, wt_int, wt_nurse, wt_blood, 
  # psu, # revise to psu_scr
  # cluster, # revise to cluster195
  # TODO: how many levels do we need
  # age, # revise to age16g5
  sex, 
  qimd,
  bmival, 
  # cholval1, # revise to cholval13
  omsysval,
  diabtotr, 
  cigst1, 
  startsmk,
  # endsmoke, # revise to endsmokg
  vegpor,
  # frtpor, # revise to porfrt
  # numsmok, # don't have
  # smokyrs, # revise to smokyrsg
  cigdyal, 
  # origin, # revise to origin 2
  # hdlval1, # revise to hdlval13
  # iregdef, # TODO: What is THIS
  bpmedd2, 
  bp1,
  diabete2, 
  sha,
  # glyhbval, # revise to glyhbvala
  alcohol = d7unitwg * 8, totalwu = totalwu * 8 / 7,
  xpsmok, 
  eqv5,
  topqual3, 
  # statina, # can't find
  # statin_px, # can't find
  # statin_tkn, # can't find
  # diage,# revised to diageg
  age35g,
  frtpor15,
  cholval13,
  hdlval13,
  glyhbvala,
  origin2, # grouped category
  smokyrsg, # measured as grouped
  endsmokg, # measured as grouped
  diageg,
  psu_scr,
  cluster195
  )]

HSE_files$hse2018[, `:=`(year = 18L, a30to06 = NA, sodium = NA, potass = NA,
                         creatin = NA, wt_urine = 0,
                         statina = NA, # can't find
                         statin_px = NA, # can't find
                         statin_tkn = NA,
                         numsmok = NA,
                         kiddiag = NA,
                         famcvd = NA
                         )]
setnames(HSE_files$hse2018, c("cholval13", "hdlval13", "bpmedd2", "origin2", "frtpor15", "glyhbvala", "diageg", "psu_scr", "cluster195", "seriala"),
         c("cholval1", "hdlval1", "bpmedd", "origin", "frtpor", "glyhbval",  "diage", "psu", "cluster", "id"))

HSE_files$hse2018[, sha := as.integer(sha)]
#replace_from_table(HSE_files$hse2018, "origin", 1:18, c(1L, 1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 6L, 7L, 9L, 9L, 9L))

# 2014 copied
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

HSE_files$hse2014 <- HSE_files$hse2014[, .(
  pserial, age90, sex, qimd, bmival, 
  #cholval12, 
  omsysval, cigst1, startsmk, 
  endsmoke, # 2015 endsmokeg
  numsmok, 
  smokyrs, # 2015 smokyrsg
  cigdyal,  
  origin3, hdlval12, # they other don't have
  bpmedd2, diabete2, sha,
  vegpor, frtpor, sodiumval, potass, creatin, wt_int,wt_nurse, wt_blood,
  wt_urine, psu, cluster, glyhbval, diabtotr, alcohol = d7unitwg * 8,
  totalwu = totalwu * 8 / 7, xpsmok, eqv5, topqual3, bp1, statina,
  statin_px, statin_tkn, diage)]

setnames(HSE_files$hse2014, c("age90", "origin3", "sodiumval", 
                              #"cholval12",
                              "hdlval12", "bpmedd2", "pserial"),
         c("age", "origin", "sodium", 
           #"cholval1", 
           "hdlval1", "bpmedd", "id"))
HSE_files$hse2014[, `:=`(year = 14L, a30to06 = NA, famcvd = NA, kiddiag = NA)] # year 2000 = year 0
HSE_files$hse2014[, sha := as.integer(sha)]
#replace_from_table(HSE_files$hse2014, "origin", c(1:2, 4:18), c(1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 7L, 6L, 9L, 9L, 9L))

# 2013
# code bp1 (diagnosed htn)
HSE_files$hse2013[, bp1 := docbp]
HSE_files$hse2013[docbp == -1L, bp1 := 2]
HSE_files$hse2013[sex == 2L & othbp == 2L, bp1 := 2]
HSE_files$hse2013[docbp == -9L | pregbp == -9L | othbp == -9L, bp1 := -9]
HSE_files$hse2013[docbp == -8L | pregbp == -8L | othbp == -8L, bp1 := -8]

HSE_files$hse2013 <- HSE_files$hse2013[, .(
  pserial, wt_int, wt_nurse, wt_blood, psu, cluster, age, sex, qimd,
  bmival, cholval12, omsysval, diabtotr, cigst1, startsmk, endsmoke, numsmok, bp1,
  smokyrs, cigdyal,  origin, hdlval12, bpmedd2, diabete2, sha, vegpor, frtpor,
  glyhbval, alcohol = d7unitwg * 8, totalwu = totalwu * 8 / 7, xpsmok, eqv5,
  topqual3, statina, statin_px, statin_tkn, diage)]
HSE_files$hse2013[, `:=`(year = 13L, a30to06 = NA, sodium = NA, potass = NA,
                         creatin = NA,  wt_urine = NA, famcvd = NA, kiddiag = NA)]
setnames(HSE_files$hse2013, c("cholval12", "hdlval12", "bpmedd2", "pserial"),
         c("cholval1", "hdlval1", "bpmedd", "id"))
#replace_from_table(HSE_files$hse2013, "origin", 1:18, c(1L, 1L, 1L, 1L, 9L, 9L, 9L, 9L, 2L, 3L, 4L, 8L, 5L, 6L, 7L, 9L, 9L, 9L))

HSE_files$hse2013[, sha := as.integer(sha)]
HSE_files$hse2013[, length(unique(psu)), by = cluster][, table(V1)] # stratas can have many psus
HSE_files$hse2013[, length(unique(cluster)), by = psu][, table(V1)] # but here one psu has 2 stratas
HSE_files$hse2013[, length(unique(cluster)), by = psu][V1 == 2, ] #2131015
HSE_files$hse2013[psu == 2131015 , cluster]
HSE_files$hse2013[cluster == 213235, psu]
HSE_files$hse2013[psu == 2131015 & cluster == 213235, psu := 2131340]

# clean-up 
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
HSE_ts[famcvd    < 0, famcvd    := NA] # TODO: not found
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
HSE_ts[, ethnicity  := origin] #TODO: revised from :HSE_ts[is.na(ethnicity), ethnicity   := origin]
HSE_ts[(startsmk   <= 0) | (startsmk   == 97), startsmk   := NA] # 97 = never smoked regularly
HSE_ts[smokyrsg <= 0 | smokyrsg == 97, smokyrsg   := NA]
HSE_ts[endsmokg <= 0 | endsmokg == 97, endsmokg   := NA]
HSE_ts[age35g < 6, age35g := NA] # not use age category below 6 (age below 13)

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

HSE_ts_13_14 = HSE_ts[year == 13L | year == 14L]
HSE_ts_15_18 = HSE_ts[year == 15| year == 16 | year == 17 | year == 18]

HSE_ts_13_14[, `:=`(age35g = NULL, smokyrsg = NULL, endsmokg = NULL)]

set.seed(1234)
# age recode from age category

HSE_ts_15_18[age35g == 6, age := sample(c(13:15), size = .N, replace = TRUE)] # TODO: check HSE
HSE_ts_15_18[age35g == 7, age := sample(c(16:19), size = .N, replace = TRUE)]
HSE_ts_15_18[age35g == 8, age := sample(c(20:24), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 9, age := sample(c(25:29), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 10, age := sample(c(30:34), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 11, age := sample(c(35:39), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 12, age := sample(c(40:44), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 13, age := sample(c(45:49), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 14, age := sample(c(50:54), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 15, age := sample(c(55:59), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 16, age := sample(c(60:64), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 17, age := sample(c(65:69), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 18, age := sample(c(70:74), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 19, age := sample(c(75:79), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 20, age := sample(c(80:84), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 21, age := sample(c(85:89), size = .N , replace = TRUE)]
HSE_ts_15_18[age35g == 22, age := sample(c(90:99), size = .N,  replace = TRUE)] #TODO: represent 100+?
HSE_ts_15_18[, age35g := NULL]


# TODO: startsmkg not in category 

# smokyrs recode from smokyrsg category
HSE_ts_15_18[smokyrsg == 1, smokyrs := sample(c(0:4), size = .N, replace = TRUE)] 
HSE_ts_15_18[smokyrsg == 2, smokyrs := sample(c(5:9), size = .N, replace = TRUE)]
HSE_ts_15_18[smokyrsg == 3, smokyrs := sample(c(10:14), size = .N, replace = TRUE)]
HSE_ts_15_18[smokyrsg == 4, smokyrs := sample(c(15:19), size = .N, replace = TRUE)]
HSE_ts_15_18[smokyrsg == 5, smokyrs := sample(c(20:29), size = .N, replace = TRUE)]
HSE_ts_15_18[smokyrsg == 6, smokyrs := sample(c(30:39), size = .N, replace = TRUE)] 
HSE_ts_15_18[smokyrsg == 7, smokyrs := sample(c(40:49), size = .N, replace = TRUE)]
HSE_ts_15_18[smokyrsg == 8 & age < 65, age := 65] # TODO: some age assgined below 65
HSE_ts_15_18[smokyrsg == 8, smokyrs := age - 15] 
HSE_ts_15_18[, smokyrsg := NULL]

# endsmoke recode from endsmokeg category
HSE_ts_15_18[endsmokg == 1, endsmoke := sample(c(0:4), size = .N, replace = TRUE)] 
HSE_ts_15_18[endsmokg == 2, endsmoke := sample(c(5:9), size = .N, replace = TRUE)]
HSE_ts_15_18[endsmokg == 3, endsmoke := sample(c(10:14), size = .N, replace = TRUE)]
HSE_ts_15_18[endsmokg == 4, endsmoke := sample(c(15:19), size = .N, replace = TRUE)]
HSE_ts_15_18[endsmokg == 5, endsmoke := sample(c(20:29), size = .N, replace = TRUE)]
HSE_ts_15_18[endsmokg == 6, endsmoke := sample(c(30:39), size = .N, replace = TRUE)] 
HSE_ts_15_18[endsmokg == 7, endsmoke := sample(c(40:49), size = .N, replace = TRUE)]
HSE_ts_15_18[endsmokg == 8 & age < 65, age := 65] # some sampled with age below 65
HSE_ts_15_18[endsmokg == 8, endsmoke := age - 15] 
HSE_ts_15_18[, endsmokg := NULL]

remove(HSE_ts)
HSE_ts = rbind(HSE_ts_13_14, HSE_ts_15_18)

HSE_ts[, dm_dur := age - diage]
HSE_ts[, diage  := NULL]

HSE_ts <- HSE_ts[between(age, 15, 99), ]
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
summary(HSE_ts)

# Standartise weights to the number of participants every year
HSE_ts[, wt_int   := sum(wt_int > 0, na.rm = TRUE)   * wt_int/sum(wt_int, na.rm = TRUE),     keyby = year]
HSE_ts[, wt_nurse := sum(wt_nurse > 0, na.rm = TRUE) * wt_nurse/sum(wt_nurse, na.rm = TRUE), keyby = year]
HSE_ts[, wt_blood := sum(wt_blood > 0, na.rm = TRUE) * wt_blood/sum(wt_blood, na.rm = TRUE), keyby = year]
# HSE_ts[, .(sum(wt_blood), .N), keyby = year]

# TODO:  column with issues:   frtpor
# #  frtpor: check how the frtpor is counted (spss code)
# 
# debug_col = HSE_ts[, .( frtpor
# )]
# ref_col = original[, .(frtpor
# )]

rm(HSE_files)
write_fst(HSE_ts, "./output/HSE_ts_13_18.fst")

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
# intake ove a week of observation

# age35g
# 6: 13-15
# 7: 16,19
# 8: 20,24
# 9: 25,29 
# 10: 30,34
# 11: 35,39
# 12: 40,44
# 13: 45,49
# 14: 50,54
# 15: 55,59
# 16: 60,64
# 17: 65,69
# 18: 70,74
# 19: 75,79
# 20: 80,84
# 21: 85,89
# 22: above 90

# STARTSMKG: (D) Age started smoking categorised
# 1 12 or under
# 2 13 years old
# 3 14 years old
# 4 15 years old
# 5 16 years old
# 6 17 years old
# 7 18 years old
# 8 19-20 years old
# 9 21-24 years old
# 10 25orover

# smokyrsg
# Refused   -9            
# Don't know   -8         
# Not applicable   -1   
# 0-4              1      
# 5-9   2        
# 10-14       3
# 15-19       4
# 20-29       5
# 30-39       6         
# 40-49       7           
# 50+  8
# Never smoked regularly  97 

# endsmokeg
# Refused   -9            
# Don't know   -8         
# Not applicable   -1   
# 0-4         1      
# 5-9         2        
# 10-14       3
# 15-19       4
# 20-29       5
# 30-39       6         
# 40-49       7           
# 50+         8


# NOTE method to prescribed medication changed since 2012. Since then they record
# Actual use rather than prescription for all categories of medicine.
