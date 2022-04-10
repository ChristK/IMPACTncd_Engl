# From https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareamidyearpopulationestimates (accessed 1/4/2022)
library(data.table)
library(CKutils)
library(fst)
library(readxl)
library(ungroup)

lsoa_pop <- list()

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE8DT2a-LSOA-syoa-unformatted-males-mid2002-to-mid2006.xls"
lsoa_pop$mid_2002_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2002",
    skip = 0
  ))
lsoa_pop$mid_2002_m[, `:=` (sex = "men", year = 2L)]
lsoa_pop$mid_2003_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2003",
    skip = 0
  ))
lsoa_pop$mid_2003_m[, `:=` (sex = "men", year = 3L)]
lsoa_pop$mid_2004_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2004",
    skip = 0
  ))
lsoa_pop$mid_2004_m[, `:=` (sex = "men", year = 4L)]
lsoa_pop$mid_2005_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2005",
    skip = 0
  ))
lsoa_pop$mid_2005_m[, `:=` (sex = "men", year = 5L)]
lsoa_pop$mid_2006_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2006",
    skip = 0
  ))
lsoa_pop$mid_2006_m[, `:=` (sex = "men", year = 6L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE8DT3a-LSOA-syoa-unformatted-females-mid2002-to-mid2006.xls"
lsoa_pop$mid_2002_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2002",
    skip = 0
  ))
lsoa_pop$mid_2002_f[, `:=` (sex = "women", year = 2L)]
lsoa_pop$mid_2003_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2003",
    skip = 0
  ))
lsoa_pop$mid_2003_f[, `:=` (sex = "women", year = 3L)]
lsoa_pop$mid_2004_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2004",
    skip = 0
  ))
lsoa_pop$mid_2004_f[, `:=` (sex = "women", year = 4L)]
lsoa_pop$mid_2005_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2005",
    skip = 0
  ))
lsoa_pop$mid_2005_f[, `:=` (sex = "women", year = 5L)]
lsoa_pop$mid_2006_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2006",
    skip = 0
  ))
lsoa_pop$mid_2006_f[, `:=` (sex = "women", year = 6L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE8DT2b-LSOA-syoa-unformatted-males-mid2007-to-mid2010.xls"
lsoa_pop$mid_2007_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2007",
    skip = 0
  ))
lsoa_pop$mid_2007_m[, `:=` (sex = "men", year = 7L)]
lsoa_pop$mid_2008_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2008",
    skip = 0
  ))
lsoa_pop$mid_2008_m[, `:=` (sex = "men", year = 8L)]
lsoa_pop$mid_2009_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2009",
    skip = 0
  ))
lsoa_pop$mid_2009_m[, `:=` (sex = "men", year = 9L)]
lsoa_pop$mid_2010_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2010",
    skip = 0
  ))
lsoa_pop$mid_2010_m[, `:=` (sex = "men", year = 10L)]
lsoa_pop$mid_2011_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2011",
    skip = 0
  ))
lsoa_pop$mid_2011_m[, `:=` (sex = "men", year = 11L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE8DT3b-LSOA-syoa-unformatted-females-mid2007-to-mid2010.xls"
lsoa_pop$mid_2007_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2007",
    skip = 0
  ))
lsoa_pop$mid_2007_f[, `:=` (sex = "women", year = 7L)]
lsoa_pop$mid_2008_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2008",
    skip = 0
  ))
lsoa_pop$mid_2008_f[, `:=` (sex = "women", year = 8L)]
lsoa_pop$mid_2009_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2009",
    skip = 0
  ))
lsoa_pop$mid_2009_f[, `:=` (sex = "women", year = 9L)]
lsoa_pop$mid_2010_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2010",
    skip = 0
  ))
lsoa_pop$mid_2010_f[, `:=` (sex = "women", year = 10L)]
lsoa_pop$mid_2011_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2011",
    skip = 0
  ))
lsoa_pop$mid_2011_f[, `:=` (sex = "women", year = 11L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE20DT1-mid-2012-lsoa-syoa-estimates-formatted.xls"
lsoa_pop$mid_2012_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2012 Males",
    skip = 3
  ))
lsoa_pop$mid_2012_m[, `:=` (sex = "men", year = 12L)]
lsoa_pop$mid_2012_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2012 Females",
    skip = 3
  ))
lsoa_pop$mid_2012_f[, `:=` (sex = "women", year = 12L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE20DT1-mid-2013-lsoa-syoa-estimates-formatted.xls"
lsoa_pop$mid_2013_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2013 Males",
    skip = 3
  ))
lsoa_pop$mid_2013_m[, `:=` (sex = "men", year = 13L)]
lsoa_pop$mid_2013_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2013 Females",
    skip = 3
  ))
lsoa_pop$mid_2013_f[, `:=` (sex = "women", year = 13L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE20DT1-mid-2014-lsoa-syoa-estimates-formatted.xls"
lsoa_pop$mid_2014_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2014 Males",
    skip = 3
  ))
lsoa_pop$mid_2014_m[, `:=` (sex = "men", year = 14L)]
lsoa_pop$mid_2014_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2014 Females",
    skip = 3
  ))
lsoa_pop$mid_2014_f[, `:=` (sex = "women", year = 14L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE20DT1-mid-2015-lsoa-syoa-estimates-formatted.xls"
lsoa_pop$mid_2015_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2015 Males",
    skip = 3
  ))
lsoa_pop$mid_2015_m[, `:=` (sex = "men", year = 15L)]
lsoa_pop$mid_2015_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2015 Females",
    skip = 3
  ))
lsoa_pop$mid_2015_f[, `:=` (sex = "women", year = 15L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE20DT1-mid-2016-lsoa-syoa-estimates-formatted.xls"
lsoa_pop$mid_2016_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2016 Males",
    skip = 3
  ))
lsoa_pop$mid_2016_m[, `:=` (sex = "men", year = 16L)]
lsoa_pop$mid_2016_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2016 Females",
    skip = 3
  ))
lsoa_pop$mid_2016_f[, `:=` (sex = "women", year = 16L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE20DT1-mid-2017-lsoa-syoa-estimates-formatted.XLS"
lsoa_pop$mid_2017_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2017 Males",
    skip = 3
  ))
lsoa_pop$mid_2017_m[, `:=` (sex = "men", year = 17L)]
lsoa_pop$mid_2017_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2017 Females",
    skip = 3
  ))
lsoa_pop$mid_2017_f[, `:=` (sex = "women", year = 17L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE21DT1a-mid-2018-on-2019-LA-lsoa-syoa-estimates-formatted.xlsx"
lsoa_pop$mid_2018_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2018 Males",
    skip = 3
  ))
lsoa_pop$mid_2018_m[, `:=` (sex = "men", year = 18L)]
lsoa_pop$mid_2018_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2018 Females",
    skip = 3
  ))
lsoa_pop$mid_2018_f[, `:=` (sex = "women", year = 18L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/SAPE22DT2-mid-2019-lsoa-syoa-estimates-unformatted.xlsx"
lsoa_pop$mid_2019_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2019 Males",
    skip = 3
  ))
lsoa_pop$mid_2019_m[, `:=` (sex = "men", year = 19L)]
lsoa_pop$mid_2019_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2019 Females",
    skip = 3
  ))
lsoa_pop$mid_2019_f[, `:=` (sex = "women", year = 19L)]

file <- "./inputs/pop_estimates_lsoa/unprocessed/pop_size/sape23dt2mid2020lsoasyoaestimatesunformatted.xlsx"
lsoa_pop$mid_2020_m <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2020 Males",
    skip = 3
  ))
lsoa_pop$mid_2020_m[, `:=` (sex = "men", year = 20L)]
lsoa_pop$mid_2020_f <-
  as.data.table(read_excel(
    file,
    sheet = "Mid-2020 Females",
    skip = 3
  ))
lsoa_pop$mid_2020_f[, `:=` (sex = "women", year = 20L)]

lu1 <- unique(lsoa_pop$mid_2011_f[, .(LSOA11CD, LAD11CD, LAD11NM)])[grep("^E01", LSOA11CD)]
lu2 <- unique(lsoa_pop$mid_2019_f[, .(`LSOA Code`, `LSOA Name`, `LA Code (2019 boundaries)`,
                               `LA name (2019 boundaries)`, `LA Code (2020 boundaries)`,
                               `LA name (2020 boundaries)`)])[grep("^E01", `LSOA Code`)]
setnames(lu2, "LSOA Code", "LSOA11CD")
lu3 <- unique(lsoa_pop$mid_2020_f[, .(`LSOA Code`, `LSOA Name`, `LA Code (2018 boundaries)`,
                                      `LA name (2018 boundaries)`, `LA Code (2021 boundaries)`,
                                      `LA name (2021 boundaries)`)])[grep("^E01", `LSOA Code`)]
setnames(lu3, "LSOA Code", "LSOA11CD")
absorb_dt(lu1, lu2)
absorb_dt(lu1, lu3)
setnames(lu1, gsub(" ", "_", names(lu1)))
# View(lu1[LAD11NM != `LA name (2019 boundaries)`, .(LAD11NM, `LA name (2019 boundaries)`)])

lapply(lsoa_pop, function(x) names(x)[1:8])
lapply(lsoa_pop, function(x) {
  setnames(x, names(x), gsub("^m", "", names(x)))
  setnames(x, names(x), gsub("^f", "", names(x)))
  setnames(x, names(x), gsub("^90.*", "90", names(x)))
  if ("LSOA Code" %in% names(x)) setnames(x, "LSOA Code", "LSOA11CD")
  if ("Area Codes" %in% names(x)) setnames(x, "Area Codes", "LSOA11CD")

  x[, c("all_ages", "All Ages", "LA Code (2018 boundaries)", "LA name (2018 boundaries)",
        "LA Code (2021 boundaries)", "LA name (2021 boundaries)",
        "LA Code (2019 boundaries)", "LA name (2019 boundaries)", "LA (2019 boundaries)",
        "LA Code (2020 boundaries)", "LA name (2020 boundaries)", "LSOA", "LSOA Name",
        "LAD11CD", "LAD11NM", "Area Names", "...3") := NULL]
})

lsoa_pop <- rbindlist(lsoa_pop)
lsoa_pop <- lsoa_pop[grep("^E01", LSOA11CD), ]
lsoa_pop[, uniqueN(LSOA11CD), keyby = .(year, sex)]
lsoa_pop[, .N, keyby = .(year, sex)]
lsoa_pop[, `:=` (
  sex = factor(sex),
  LSOA11CD = factor(LSOA11CD)
)]
setkey(lsoa_pop, year, LSOA11CD, sex)
setcolorder(lsoa_pop)

# Distribute ages above 90 using national estimates from a PCLM2D model
tt <- melt(lsoa_pop, 1:3, variable.name = "age", value.name = "pops")
tt[, age := as.integer(as.character(age))]
tt <- tt[, .(pops = sum(pops)), keyby = .(year, age, sex)]
tt <- dcast(tt, age + sex ~ year, value.var = "pops")
ttm <- tt[sex == "men", pclm2D(age, as.data.frame(.SD), nlast = 10)$fitted, .SDcols = !c("age", "sex")]
ttm <- as.data.table(ttm)
ttm <- ttm[90:100, lapply(.SD, function(x) x/sum(x))] # rows ages 90-100
ttm <- transpose(cbind("age" = 90:100, ttm), make.names = "age")
ttm <- cbind("year" = 2:20, "sex" = "men", ttm)

ttw <- tt[sex == "women", pclm2D(age, as.data.frame(.SD), nlast = 10)$fitted, .SDcols = !c("age", "sex")]
ttw <- as.data.table(ttw)
ttw <- ttw[90:100, lapply(.SD, function(x) x/sum(x))] # rows ages 90-100
ttw <- transpose(cbind("age" = 90:100, ttw), make.names = "age")
ttw <- cbind("year" = 2:20, "sex" = "women", ttw)
tt <- rbind(ttm, ttw)


setnames(lsoa_pop, "90", "90+")
absorb_dt(lsoa_pop, tt)
for (j in as.character(90:100)) {
  set(lsoa_pop, NULL, j, lsoa_pop$`90+` * lsoa_pop[[j]])
}

replace_from_table(
  lu1,
  "LAD11NM",
  from = c(
    "Basingstoke and Dean",
    "Bristol, City of",
    "Herefordshire, County of",
    "Ise of Anglesey",
    "Kingston upon Hull, City of"
  ),
  to = c(
    "Basingstoke and Deane",
    "Bristol",
    "Herefordshire",
    "Isle of Anglesey",
    "Kingston upon Hull"
  )
)
lsoa_pop[, "90+" := NULL]
setkey(lsoa_pop, year, LSOA11CD, sex)
write_fst(lsoa_pop, "./inputs/pop_estimates_lsoa/LSOA_mid_year_population_estimates.fst", 100)

nat <- lsoa_pop[, lapply(.SD, sum), keyby = .(year, sex), .SDcols = as.character(0:100)]
nat <- melt(nat, c("year", "sex"), variable.name = "age", value.name = "pops")
nat[, age := as.character(as.integer(age))]
write_fst(nat, "./inputs/pop_estimates_lsoa/national_mid_year_population_estimates.fst", 100)

# Convert from mid-year (1st-2nd July) to 1st of April (3 months)
foo <- function(x) {
  3/4 * x + 1/4 * shift(x)
}
lsoa_pop[, c(as.character(0:100), "90+") := lapply(.SD, foo),
         keyby = .(sex, LSOA11CD),
         .SDcols = c(as.character(0:100), "90+")]
lsoa_pop <- lsoa_pop[year > 2]
setkey(lsoa_pop, year, LSOA11CD, sex)
write_fst(lsoa_pop, "./inputs/pop_estimates_lsoa/LSOA_1st_April_population_estimates.fst", 100)



# merge with Townsend
townsend <- fread("./inputs/aux_data/2011 UK Townsend Deprivation Scores LSOA.csv")
townsend <- townsend[grep("^E01", GEO_CODE)]
lu1[townsend, on = "LSOA11CD == GEO_CODE",
    `:=` (tds = i.TDS,
          tds_quintile = i.quintile)]


# merge with IMD2015
imd <- fread("./inputs/aux_data/imd2015eng.csv")
lu1[imd, on = "LSOA11CD == LSOA code (2011)", `:=` (
  imd = `i.Index of Multiple Deprivation (IMD) Score`,
  dimd = `i.Index of Multiple Deprivation (IMD) Decile (where 1 is most deprived 10% of LSOAs)`, # 1 most deprived, opposite of HSE
  LAD13CD = `i.Local Authority District code (2013)`,
  LAD13NM = `i.Local Authority District name (2013)`
)]
replace_from_table(lu1, "dimd", 1:10, rep(1:5, each = 2), "qimd")

# merge with LAD16NM/WD16CD
tt <- fread("./inputs/aux_data/Lower_Layer_Super_Output_Area_2011_to_Ward_to_Local_Authority_District_December_2016_Lookup_in_England_and_Wales.csv")
tt <- tt[grep("^E01", LSOA11CD)]
lu1[tt, on = "LSOA11CD", `:=` (
  WD16CD = i.WD16CD,
  WD16NM = i.WD16NM,
  LAD16CD = i.LAD16CD,
  LAD16NM = i.LAD16NM
)]

# merge with MSOA2011 SHA etc
tt <- fread("./inputs/aux_data/Output_Area_to_Lower_Layer_Super_Output_Area_to_Middle_Layer_Super_Output_Area_to_Local_Authority_District_December_2017_Lookup_in_Great_Britain__Classification_Version_2.csv")
ttt <- fread("./inputs/aux_data/Output_Area_to_Primary_Care_Organisation_to_Strategic_Health_Authority_December_2011_Lookup_in_England_and_Wales.csv")

tt[ttt, on = "OA11CD", `:=` (
  SHA11CD = i.SHA11CD,
  SHA11NM = i.SHA11NM
)]
tt <- tt[grep("^E01", LSOA11CD)]
tt[, c("OA11CD", "OAC11CD", "OAC11NM", "FID", "SOAC11CD", "SOAC11NM") := NULL]
tt <- unique(tt)
lu1[tt, on = "LSOA11CD", `:=` (
  MSOA11CD = i.MSOA11CD,
  MSOA11NM = i.MSOA11NM,
  LAD17CD = i.LAD17CD,
  LAD17NM = i.LAD17NM,
  RGN11CD = i.RGN11CD,
  RGN11NM = i.RGN11NM,
  SHA11CD = i.SHA11CD,
  SHA11NM = i.SHA11NM
)]

# LAD17 to CCG17
tt <- fread("./inputs/aux_data/Lower_Layer_Super_Output_Area_2011_to_Clinical_Commissioning_Group_to_Local_Authority_District_April_2017_Lookup_in_England_Version_4.csv", select = c("CCG17CD", "CCG17CDH", "CCG17NM", "LAD17CD"))
tt <- tt[, .N, by = .(LAD17CD, CCG17CD, CCG17CDH, CCG17NM)] # For areas were i LAD matches to many CCGs N is the weight.
setkey(tt, N)
tt <- tt[, first(.SD), by = LAD17CD] # Assuming every LAD is represented by the CCG with highest number of LSOAs
lu1[tt, on = "LAD17CD", `:=` (
  CCG17CD = i.CCG17CD,
  CCG17CDH = i.CCG17CDH,
  CCG17NM = i.CCG17NM
)]


# The Strategic Health Authorities are coterminous with government office regions (RGN11),
# except that the large South East England region is divided into two:
# South Central and South East Coast

lu1[, `:=` (
  qimd = factor(
    qimd,
    levels = 1:5,
    labels = c("1 most deprived", "2", "3", "4", "5 least deprived")
  ),
  dimd = factor(
    dimd,
    levels = 1:10,
    labels = c("1 most deprived", 2:9, "10 least deprived")
  ),
  LSOA_Name = factor(LSOA_Name),
  LAD13CD  = factor(LAD13CD),
  LAD13NM  = factor(LAD13NM),
  LAD16CD  = factor(LAD16CD),
  LAD16NM  = factor(LAD16NM),
  LAD17CD  = factor(LAD17CD),
  LAD17NM  = factor(LAD17NM),
  LAD11CD  = factor(LAD11CD),
  LAD11NM  = factor(LAD11NM),
  WD16CD   = factor(WD16CD),
  WD16NM   = factor(WD16NM),
  MSOA11CD = factor(MSOA11CD),
  MSOA11NM = factor(MSOA11NM),
  RGN11CD  = factor(RGN11CD),
  RGN11NM  = factor(RGN11NM),
  LSOA11CD = factor(LSOA11CD),
  SHA11CD  = factor(SHA11CD),
  SHA11NM  = factor(SHA11NM, levels = c(
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
  )),
  CCG17CD = factor(CCG17CD),
  CCG17CDH = factor(CCG17CDH),
  CCG17NM = factor(CCG17NM),
  `LA_Code_(2018_boundaries)` = factor(`LA_Code_(2018_boundaries)`),
  `LA_name_(2018_boundaries)` = factor(`LA_name_(2018_boundaries)`),
  `LA_Code_(2019_boundaries)` = factor(`LA_Code_(2019_boundaries)`),
  `LA_name_(2019_boundaries)` = factor(`LA_name_(2019_boundaries)`),
  `LA_Code_(2020_boundaries)` = factor(`LA_Code_(2020_boundaries)`),
  `LA_name_(2020_boundaries)` = factor(`LA_name_(2020_boundaries)`),
  `LA_Code_(2021_boundaries)` = factor(`LA_Code_(2021_boundaries)`),
  `LA_name_(2021_boundaries)` = factor(`LA_name_(2021_boundaries)`)
)]
setcolorder(lu1, c("LSOA11CD", "imd", "dimd", "qimd", "tds", "tds_quintile",
                   "LAD13CD", "LAD13NM", "LAD16CD", "LAD16NM", "LAD17CD", "LAD17NM",
                   "WD16CD", "WD16NM",  "MSOA11CD", "MSOA11NM",  "RGN11CD",
                   "RGN11NM", "SHA11CD", "SHA11NM", "CCG17CD", "CCG17CDH", "CCG17NM",
                   "LA_Code_(2019_boundaries)", "LA_name_(2019_boundaries)",
                   "LA_Code_(2020_boundaries)", "LA_name_(2020_boundaries)",
                   "LA_Code_(2018_boundaries)", "LA_name_(2018_boundaries)",
                   "LA_Code_(2021_boundaries)", "LA_name_(2021_boundaries)"
))

write_fst(lu1, "./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst", 100)
