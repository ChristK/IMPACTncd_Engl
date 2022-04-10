# From NOMIS ONS Census 2011 population by ward age sex ethnicity
library(data.table)
library(CKutils)
library(fst)

"./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/"

dt1 <- fread("./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/864518308 white.csv", skip = 9, fill = TRUE)
dt1[, sex := "women"]
dt1[1:34764, sex := "men"]
dt1[, ethnicity := "white"]
# View(dt1[c(1L, 34755:34778, 69532:.N)])
dt1 <- dt1[!c(1L, 34755:34778, 69532:.N)]
dt1 <- melt(dt1, c(1, 2, 7, 8), variable.name = "agegrp", value.name = "population_size")
dt1[, population_size := as.numeric(population_size)]
dt1[, unique(counts(sex))]
dt1[, unique(counts(agegrp))]

dt2 <- fread("./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/906210981 indian.csv", skip = 9, fill = TRUE)
dt2[, sex := "women"]
dt2[1:34764, sex := "men"]
dt2[, ethnicity := "indian"]
# View(dt2[c(1L, 34755:34778, 69532:.N)])
dt2 <- dt2[!c(1L, 34755:34778, 69532:.N)]
dt2 <- melt(dt2, c(1, 2, 7, 8), variable.name = "agegrp", value.name = "population_size")
dt2[, population_size := as.numeric(population_size)]
dt2[, unique(counts(sex))]
dt2[, unique(counts(agegrp))]

dt3 <- fread("./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/906210981 pakistani.csv", skip = 9, fill = TRUE)
dt3[, sex := "women"]
dt3[1:34764, sex := "men"]
dt3[, ethnicity := "pakistani"]
# View(dt3[c(1L, 34755:34778, 69532:.N)])
dt3 <- dt3[!c(1L, 34755:34778, 69532:.N)]
dt3 <- melt(dt3, c(1, 2, 7, 8), variable.name = "agegrp", value.name = "population_size")
dt3[, population_size := as.numeric(population_size)]
dt3[, unique(counts(sex))]
dt3[, unique(counts(agegrp))]

dt4 <- fread("./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/1021314074 bangladeshi.csv", skip = 9, fill = TRUE)
dt4[, sex := "women"]
dt4[1:34764, sex := "men"]
dt4[, ethnicity := "bangladeshi"]
# View(dt4[c(1L, 34755:34778, 69532:.N)])
dt4 <- dt4[!c(1L, 34755:34778, 69532:.N)]
dt4 <- melt(dt4, c(1, 2, 7, 8), variable.name = "agegrp", value.name = "population_size")
dt4[, population_size := as.numeric(population_size)]
dt4[, unique(counts(sex))]
dt4[, unique(counts(agegrp))]

dt5 <- fread("./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/1114011716 other asian.csv", skip = 9, fill = TRUE)
dt5[, sex := "women"]
dt5[1:34764, sex := "men"]
dt5[, ethnicity := "other asian"]
# View(dt5[c(1L, 34755:34778, 69532:.N)])
dt5 <- dt5[!c(1L, 34755:34778, 69532:.N)]
dt5 <- melt(dt5, c(1, 2, 7, 8), variable.name = "agegrp", value.name = "population_size")
dt5[, population_size := as.numeric(population_size)]
dt5[, unique(counts(sex))]
dt5[, unique(counts(agegrp))]

dt6 <- fread("./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/1067219951 chinese.csv", skip = 9, fill = TRUE)
dt6[, sex := "women"]
dt6[1:34764, sex := "men"]
dt6[, ethnicity := "chinese"]
# View(dt6[c(1L, 34755:34778, 69532:.N)])
dt6 <- dt6[!c(1L, 34755:34778, 69532:.N)]
dt6 <- melt(dt6, c(1, 2, 7, 8), variable.name = "agegrp", value.name = "population_size")
dt6[, population_size := as.numeric(population_size)]
dt6[, unique(counts(sex))]
dt6[, unique(counts(agegrp))]

dt7 <- fread("./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/1143515776 black african.csv", skip = 9, fill = TRUE)
dt7[, sex := "women"]
dt7[1:34764, sex := "men"]
dt7[, ethnicity := "black african"]
# View(dt7[c(1L, 34755:34778, 69532:.N)])
dt7 <- dt7[!c(1L, 34755:34778, 69532:.N)]
dt7 <- melt(dt7, c(1, 2, 7, 8), variable.name = "agegrp", value.name = "population_size")
dt7[, population_size := as.numeric(population_size)]
dt7[, unique(counts(sex))]
dt7[, unique(counts(agegrp))]

dt8 <- fread("./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/118255723 black caribbean.csv", skip = 9, fill = TRUE)
dt8[, sex := "women"]
dt8[1:34764, sex := "men"]
dt8[, ethnicity := "black caribbean"]
# View(dt8[c(1L, 34755:34778, 69532:.N)])
dt8 <- dt8[!c(1L, 34755:34778, 69532:.N)]
dt8 <- melt(dt8, c(1, 2, 7, 8), variable.name = "agegrp", value.name = "population_size")
dt8[, population_size := as.numeric(population_size)]
dt8[, unique(counts(sex))]
dt8[, unique(counts(agegrp))]

dt <- rbind(dt1, dt2, dt3, dt4, dt5, dt6, dt7, dt8)
dt[, unique(ethnicity)]

dt9 <- fread("./inputs/pop_estimates_lsoa/unprocessed/ethnic_group_by_sex_by_age_lsoa/591013939 all.csv", skip = 9, fill = TRUE)
dt9[, sex := "women"]
dt9[1:34764, sex := "men"]
dt9[, ethnicity := "other"] # atm all, but adter subtraction in next step will be other
# View(dt9[c(1L, 34755:34778, 69532:.N)])
dt9 <- dt9[!c(1L, 34755:34778, 69532:.N)]
dt9 <- melt(dt9, c(1, 2, 7, 8), variable.name = "agegrp", value.name = "population_size")
dt9[, population_size := as.numeric(population_size)]
dt9[, unique(counts(sex))]
dt9[, unique(counts(agegrp))]

dt[, uniqueN(`2011 super output area - lower layer`)]
tt <- dt[, sum(population_size), by = .(agegrp, sex, mnemonic, `2011 super output area - lower layer`)]
dt9[, sum(population_size)]
dt9[, .N] == tt[, .N]
dt9[tt, on = c("agegrp", "sex", "mnemonic", "2011 super output area - lower layer"), population_size := population_size - i.V1]
dt9[, summary(population_size)]
dt <- rbind(dt, dt9)

dt[, unique(counts(sex))]
dt[, unique(counts(agegrp))]
dt[, unique(counts(ethnicity))]
dt[, unique(counts(mnemonic))]
dt[, unique(counts(`2011 super output area - lower layer`))]

dt[, agegrp := as.character(agegrp)]
dt[, sort(unique(agegrp))]

replace_from_table(dt, "agegrp", dt[, sort(unique(agegrp))],
                   c("00-24", "25-49", "50-64", "65+"))
dt[, sort(unique(agegrp))]

# discard Wales wards
dt <- dt[substr(mnemonic, 1, 1) == "E", ]

dt[, ethn_prop := population_size/sum(population_size), by = c("2011 super output area - lower layer", "mnemonic", "sex", "agegrp")]
# From above, NAs are produced  View(et[is.na(prop), ]) for "Manchester 059B", "Leeds 111B", "Leeds 112D", "Sheffield 073E",
# "Manchester 057B", and "Manchester 013F" because no population existed over 65 (div by 0)
# For these areas I assume all white
dt[is.na(ethn_prop) & ethnicity == "white", ethn_prop := 1]
dt[is.na(ethn_prop), ethn_prop := 0]

tt <- data.table(age = 0:100, agegrp = c(rep("00-24", 25), rep("25-49", 25), rep("50-64", 15), rep("65+", 36)))
dt <- tt[dt, on = "agegrp", allow.cartesian = TRUE]
dt <-
  dt[, .(
    age,
    sex = factor(sex),
    ethnicity = factor(
      ethnicity,
      levels = c(
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
    LSOA11CD = factor(mnemonic),
    ethn_prop
  )]
dt <- dcast(dt, age + sex + LSOA11CD ~ ethnicity, value.var = "ethn_prop")

setkey(dt, LSOA11CD, age, sex)

write_fst(dt, "./inputs/pop_estimates_lsoa/ethn2011_pct.fst", compress = 100)
