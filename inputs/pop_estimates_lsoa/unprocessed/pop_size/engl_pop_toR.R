library(data.table)
library(CKutils)
library(fst)
library(readxl)


pop <- as.data.table(read_excel("inputs/pop_estimates_lsoa/unprocessed/pop_size/ukpopulationestimates18382022.xlsx",
    sheet = "Table 11", skip = 3
))
setnames(pop, tolower(names(pop)))
pop <- pop[sex != "Persons"]
pop <- pop[age != "All Ages"]
pop[, sex := factor(sex, c("Males", "Females"), c("men", "women"))]
pop <- melt(pop, c("age", "sex"), variable.name = "year", value.name = "pops")
pop[, year := as.integer(gsub("mid-", "", year))]
tt <- pop[age == "90+"]
tt

# get the distribution of 90+ from the national estimates at https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/ageing/bulletins/estimatesoftheveryoldincludingcentenarians/2002to2022
# ages <- 90:105
# women2002 <- c(61700,50610,41560,32770,25570,18470,13100,9090,6330,4190,2600,1570,960,500,260,280)
# women2002 <- women2002 / sum(women2002)
# women2022 <- c(76940,66100,55720,44220,34590,26890,20300,14850,10540,7370,5050,3380,2140,810,440,570)
# women2022 <- women2022 / sum(women2022)
# men2002   <- c(23260,17300,13110,9300,6640,4430,2840,1830,1050,620,350,210,110,40,20,20)
# men2002 <- men2002 / sum(men2002)
# men2022 <- c(45650,37040,28890,21340,15600,11290,7780,5200,3310,2090,1260,750,430,150,70,70)
# men2022 <- men2022 / sum(men2022)

# approx(c(2002, 2022), c(women2002[1], women2022[1]), xout = c(2002:2022))$y

# ttt <- CJ(age = 90:105, year = 2002:2022, sex = c("men", "women"))
# for (i in 91:105) {
#   ttt[age == i & sex == "women", prop := approx(c(2002, 2022), c(women2002[i - 89], women2022[i - 89]), xout = c(2002:2022))$y]
#   ttt[age == i & sex == "men", prop := approx(c(2002, 2022), c(men2002[i - 89], men2022[i - 89]), xout = c(2002:2022))$y]
# }
# ttt[, tmp := sum(prop, na.rm = TRUE), by = .(year, sex)]
# ttt[is.na(prop), prop := 1 - tmp]
# ttt[, tmp := sum(prop), by = .(year, sex)]
# ttt[, unique(tmp)]
# ttt[, tmp := NULL]
# ttt[tt, on = c("sex", "year"), pops := i.pops * prop]
# ttt[, prop := NULL]

# TODO use https://www.ons.gov.uk/peoplepopulationandcommunity/birthsdeathsandmarriages/ageing/datasets/midyearpopulationestimatesoftheveryoldincludingcentenariansengland
ttm <- as.data.table(read_excel("inputs/pop_estimates_lsoa/unprocessed/pop_size/englandevo2022.xls",
    sheet = "England Males", skip = 3
))
ttm[, c( "90 & over \n[note 4]", "90-99", "100 & over") := NULL]
setnames(ttm, c("Year", "105 & over"), c("year", "105"))
ttm[, sex := "men"]
ttm <- melt(ttm, id.vars = c("year", "sex"), variable.name = "age", value.name = "pops")
ttf <- as.data.table(read_excel("inputs/pop_estimates_lsoa/unprocessed/pop_size/englandevo2022.xls",
    sheet = "England Females", skip = 3
))
ttf[, c("90 & over \n[note 4]", "90-99", "100 & over") := NULL]
setnames(ttf, c("Year", "105 & over"), c("year", "105"))
ttf[, sex := "women"]
ttf <- melt(ttf, id.vars = c("year", "sex"), variable.name = "age", value.name = "pops")
ttt <- rbind(ttm, ttf)
ttt[, prop := pops / sum(pops), keyby = .(year, sex)]
ttt[tt, on = c("sex", "year"), pops := i.pops * prop]
ttt[, prop := NULL]
ttt[, age := as.integer(as.character(age))]

pop <- pop[age != "90+"]
pop <- rbind(pop[between(year, 2002, 2022)], ttt)
pop[, age := as.integer(as.character(age))]
pop[, year := as.integer(year) - 2000L]
setkey(pop, year, sex, age)
pop[, pops := as.integer(round(pops))]
write_fst(pop, "./inputs/pop_estimates_lsoa/national_pop_est.fst", 100)
