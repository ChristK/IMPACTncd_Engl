library(fst)
library(data.table)

indx_hlp <-
            read_fst("./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst",
                     as.data.table = TRUE)

pops <- read_fst("./inputs/pop_estimates_lsoa/LSOA_1st_April_population_estimates.fst", as.data.table = TRUE)

pops <-
  melt(
    pops,
    grep("^[0-9]", names(pops), value = TRUE, invert = TRUE),
    variable.name = "age",
    value.name = "population_size",
    variable.factor = FALSE
  )
pops[, age := as.integer(age)]
pops[, year := as.integer(year)]


# For the initial year (2013) find the distribution of population by dimd and age
pops13 <- pops[year == 13L]
pops13[indx_hlp, on = "LSOA11CD", `:=` (dimd = i.dimd, LAD17CD = i.LAD17CD)]
pops13 <- pops13[, .(population_size = sum(population_size)), keyby = .(LAD17CD, age, sex, dimd)]
pops13[, denom := sum(population_size), by = .(LAD17CD, age, sex)]
pops13[, dimdprop := population_size / denom]
pops13 <- dcast(pops13, LAD17CD + age + sex ~ dimd, value.var = "dimdprop")
setnafill(pops13, type = "c", fill = 0, cols = c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived"))

# check the trends in areas with the highest and lowest deprivation
pops13[between(age, 30, 99), sum(`10 least deprived`)]

tt <- pops13[between(age, 30, 99), sum(`10 least deprived`), keyby =  LAD17CD]
setkey(tt, V1)
rich <- tail(tt, 10L)$LAD17CD
poor <- head(tt, 10L)$LAD17CD

# assume that the growth in population by LAD is the same across all dimd in the LAD
proj <- read_fst("./inputs/pop_projections/lad17_proj.fst", as.data.table = TRUE)[age < 100, ]
proj[LAD17CD %in% rich & between(age, 30, 99), sum(pops), keyby = year][, plot(year, V1)]
proj[LAD17CD %in% poor & between(age, 30, 99), sum(pops), keyby = year][, plot(year, V1)]


tt <- pops13[age == 30L, ][, age := NULL]
CKutils::absorb_dt(proj, tt)
pops13[, age_in_2013 := age]
for (yr in sort(unique(proj$year))) {
  print(yr)
  pops13[, year := yr]
  pops13[, age := age_in_2013 + year - 13L]
  CKutils::absorb_dt(proj, pops13[between(age, 0, 99)], exclude_col = c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived"))
  proj[, age_in_2013 := NULL]
  # proj <- proj[pops13[between(age, 0, 99)], on = c("year", "age", "sex", "LAD17CD")]
}
proj[year == yr]
dcast(proj[, anyNA(`5`), keyby = .(year, age)], age~year, value.var = "V1")

View(proj[between(age, 20, 35) & between(year, 13, 23) & sex == "men" & LAD17CD == "E07000007"])


proj[, c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived") := lapply(.SD, function(x) x * pops), .SDcols = c("1 most deprived", "2", "3", "4", "5", "6", "7", "8", "9", "10 least deprived")]
proj[, c("LAD17CD", "pops", "LAD17NM") := NULL]
proj <- melt(proj, 1:3, variable.name = "dimd", value.name = "population_size", variable.factor = TRUE)
proj <- proj[year >= 13L, .(population_size = sum(population_size)), keyby = .(year, age, sex, dimd)]
proj[, year := year + 2000L]
CKutils::to_agegrp(proj, 5L, 99L)
proj <- proj[, .(population_size = as.integer(sum(population_size))), keyby = .(year, agegrp, sex, dimd)]

# Calibrate to national population estimates and projections.
tt <- read_fst("./inputs/pop_projections/national_proj.fst", as.data.table = TRUE)
ttt <- read_fst("./inputs/pop_estimates_lsoa/national_pop_est.fst", as.data.table = TRUE)
tt <- rbind(ttt, tt)
CKutils::to_agegrp(tt, 5L, 99L)
# Assume 95-99 is 95+ as agreed with Toby & co
tt[is.na(agegrp), agegrp := "95-99"]
tt <- tt[, .(target_pops = sum(pops)), keyby = .(year, sex, agegrp)]
tt[, year := year + 2000L]

proj[, clbr := sum(population_size), keyby = .(year, sex, agegrp)]
proj[tt, on = c("year", "sex", "agegrp"), population_size_calibrated := as.numeric(population_size * i.target_pops / clbr)]
proj[, sum(population_size_calibrated), keyby = .(year, sex, agegrp)]
tt[between(year, 2013, 2043)]
fwrite(proj, "./pop_proj_by_dimd.csv")

proj[agegrp == "30-34" & sex == "men" & dimd == "1 most deprived", plot(year, population_size_calibrated / 1000, col = 1, ylim = c(100, 300), type = "l", lwd = 2, main = "Population projections by deprivation decile", xlab = "Year", ylab = "Population size (thousands)")]
proj[agegrp == "30-34" & sex == "men" & dimd == "10 least deprived", lines(year, population_size_calibrated / 1000, col = 2, type = "l", lwd = 2)]
proj[agegrp == "30-34" & sex == "men" & dimd == "5", lines(year, population_size_calibrated / 1000, col = 3, type = "l", lwd = 2)]

proj[, .(population_size = as.integer(sum(population_size_calibrated))), keyby = .(year, dimd)][dimd == "1 most deprived", plot(year, population_size/1000, col = 1, ylim = c(100, 7000), type = "l", lwd = 2, main = "Population projections by deprivation decile", xlab = "Year", ylab = "Population size (thousands)")]
proj[, .(population_size = as.integer(sum(population_size_calibrated))), keyby = .(year, dimd)][dimd == "10 least deprived", lines(year, population_size/1000, col = 2, type = "l", lwd = 2)]
proj[, .(population_size = as.integer(sum(population_size_calibrated))), keyby = .(year, dimd)][dimd == "5", lines(year, population_size/1000, col = 3, type = "l", lwd = 2)]



# Validate against modelled population
library(readxl)
library(data.table)
bn <- read_excel("C:/Users/Chris Kypridemos/Downloads/pop_proj_by_dimd.xlsx", sheet = "pop size by year-agegrp-sex-imd", skip = 1L)
setDT(bn)
bn <- bn[scenario == "sc0", .(year, agegrp,   sex,  dimd, `pop_size_50.0%`)]
setnames(bn, "pop_size_50.0%", "popsref")

tt <- fread("C:/Users/Chris Kypridemos/Downloads/pop_proj_by_dimd.csv")[agegrp %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99")]

oldtt <- read_excel("C:/Users/Chris Kypridemos/Downloads/pop_proj_by_dimd.xlsx", sheet = "pop_proj_by_dimd", skip = 1L)
setDT(oldtt)
oldtt <- oldtt[agegrp %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99")]

CKutils::absorb_dt(tt, bn)
CKutils::absorb_dt(oldtt, bn)

tt[, check := round((popsref - population_size)/popsref, 2)]
oldtt[, check := round((popsref - population_size)/popsref, 2)]

tt[year == 2013 & agegrp == "30-34"]
tt[year == 2043 & agegrp == "30-34"]
oldtt[year == 2043 & agegrp == "30-34"]


tt[dimd == "1 most deprived" & agegrp == "30-34" & sex == "men"][, plot(year, population_size)]
tt[dimd == "10 least deprived" & agegrp == "30-34" & sex == "men"][, plot(year, popsref)]

tt[dimd == "1 most deprived" , sum(popsref, na.rm = T), keyby = .(year)][, plot(year, V1)]
oldtt[dimd == "1 most deprived" , sum(population_size, na.rm = T), keyby = .(year)][, lines(year, V1, col = "red")]
tt[dimd == "1 most deprived" , sum(population_size, na.rm = T), keyby = .(year)][, lines(year, V1, col = "blue")]

plot_all_ages <- function(dimd_lvl = as.character(c("1 most deprived", 2:9, "10 least deprived"))) {
  tt[dimd %in% dimd_lvl , sum(popsref, na.rm = T), keyby = .(year)][, plot(year, V1, main = paste0("DIMD: ", paste(dimd_lvl, collapse = ", "), ", ages:30-99"))]
  oldtt[dimd %in% dimd_lvl , sum(population_size, na.rm = T), keyby = .(year)][, lines(year, V1, col = "red")]
  tt[dimd %in% dimd_lvl , sum(population_size, na.rm = T), keyby = .(year)][, lines(year, V1, col = "blue")]
}
plot_all_ages()

plot_all_ages("1 most deprived")
plot_all_ages("2")
plot_all_ages("3")
plot_all_ages("4")

plot_all_ages("5")
plot_all_ages("6")
plot_all_ages("7")
plot_all_ages("8")
plot_all_ages("9")
plot_all_ages("10 least deprived")

tt[dimd == "10 least deprived" , sum(popsref, na.rm = T), keyby = .(year)][, plot(year, V1)]
oldtt[dimd == "10 least deprived" , sum(population_size, na.rm = T), keyby = .(year)][, lines(year, V1, col = "red")]
tt[dimd == "10 least deprived" , sum(population_size, na.rm = T), keyby = .(year)][, lines(year, V1, col = "blue")]

tt[dimd == "1 most deprived" , sum(population_size, na.rm = T), keyby = .(year)][, plot(year, V1)]
tt[dimd == "10 least deprived" , sum(population_size, na.rm = T), keyby = .(year)][, plot(year, V1)]