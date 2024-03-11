library(readxl)
library(data.table)
library(fst)
library(ungroup)

# National ----
pop <- as.data.table(read_excel("inputs/pop_projections/unprocessed/enpppopendata2021.xlsx",
                                sheet = "Population"))
setnames(pop, tolower(names(pop)))
pop[, c("2021", "2022") := NULL]
pop[, sex := factor(sex, c("Males", "Females"), c("men", "women"))]
pop <- melt(pop, c("age", "sex"), variable.name = "year", value.name = "pops")
pop[age == "105 - 109", age := 105]
pop[age == "110 and over", age := 110]
pop[, `:=` (
  age = as.integer(age),
  year = as.integer(as.character(year)) - 2000L
)]
setkey(pop, year, sex, age)
write_fst(pop, "./inputs/pop_projections/national_proj.fst", 100)

f <- fread("./inputs/pop_projections/unprocessed/2018 SNPP Population females.csv", header = TRUE)
setnames(f, tolower(names(f)))
f[, component := NULL]
f[, sex := "women"]
f[age_group == "90 and over", age_group := 90]
f[, age := as.integer(age_group)] # All_ages is converted to NA. That's desirable!
f[, age_group := NULL]
f <- na.omit(f)

m <- fread("./inputs/pop_projections/unprocessed/2018 SNPP Population males.csv", header = TRUE)
setnames(m, tolower(names(m)))
m[, component := NULL]
m[, sex := "men"]
m[age_group == "90 and over", age_group := 90]
m[, age := as.integer(age_group)]
m[, age_group := NULL]
m <- na.omit(m)

pops <- rbind(m, f)
pops <- melt(pops, c("area_code", "area_name", "age", "sex"), variable.name = "year", value.name = "pops")
indx <- read_fst("./inputs/pop_estimates_lsoa/lsoa_to_locality_indx.fst", as.data.table = TRUE)

# Regions ----
popr <- pops[grepl("^E12", area_code)]
all.equal(sort(unique(popr$area_code)), levels(indx$RGN11CD))
setnames(popr, c("area_code", "area_name"), c("RGN11CD", "RGN11NM"))

# Distribute ages above 90 using national estimates from a PCLM2D model
tt <- dcast(popr, age + sex + RGN11CD + RGN11NM ~ year, value.var = "pops")
ttm <- tt[, as.data.table(pclm2D(age, as.data.frame(.SD), nlast = 10)$fitted),
          .SDcols = !c("age", "sex", "RGN11CD", "RGN11NM"), keyby = .(sex, RGN11CD, RGN11NM)]
ttm[, age := rep(0:99, 2 * uniqueN(RGN11CD))]
ttm <- ttm[age >= 90]
ttm <- melt(ttm, c("RGN11CD", "RGN11NM", "age", "sex"), variable.name = "year", value.name = "pops")
popr <- popr[age < 90]
popr <- rbind(popr, ttm)

# include pop estimates before 2020
file <- "./inputs/pop_estimates_lsoa/LSOA_mid_year_population_estimates.fst"
ol <- read_fst(file, as.data.table = TRUE)
ol[indx, on = "LSOA11CD", `:=` (RGN11CD = i.RGN11CD, RGN11NM = RGN11NM)]
ol <- ol[, lapply(.SD, sum), .SDcols = as.character(0:100),
         keyby = .(year, RGN11CD, RGN11NM, sex)]
ol <- melt(ol, c("year", "RGN11CD", "RGN11NM", "sex"), variable.name = "age", value.name = "pops")
popr[, year := as.integer(as.character(year)) - 2000L]
popr <- rbind(popr[year > 20L], ol)
popr[, RGN11NM := NULL]
popr[indx, on = "RGN11CD", `:=` (RGN11NM = i.RGN11NM)]
popr[, `:=` (RGN11CD = factor(RGN11CD),
             RGN11NM = factor(RGN11NM),
             sex = factor(sex),
             year = as.integer(as.character(year)),
             age = as.integer(as.character(age))
)]
setkey(popr, year, sex, age, RGN11CD)
setcolorder(popr)
all.equal(sort(unique(popr$RGN11NM)), sort(unique(indx$RGN11NM)))
write_fst(popr, "./inputs/pop_projections/rgn11_proj.fst", 100)



# LAs ----
popl <- pops[grepl(paste("^E0", 6:9, collapse = "|", sep = ""), area_code)]
all.equal(sort(unique(popl$area_code)), levels(indx$LAD17CD))
setnames(popl, c("area_code", "area_name"), c("LAD17CD", "LAD17NM"))

# Distribute ages above 90 using national estimates from a PCLM2D model
tt <- dcast(popl, age + sex + LAD17CD + LAD17NM ~ year, value.var = "pops")
ttm <- tt[, as.data.table(pclm2D(age, as.data.frame(.SD), nlast = 10)$fitted),
          .SDcols = !c("age", "sex", "LAD17CD", "LAD17NM"), keyby = .(sex, LAD17CD, LAD17NM)]
ttm[, age := rep(0:99, 2 * uniqueN(LAD17CD))]
ttm <- ttm[age >= 90]
ttm <- melt(ttm, c("LAD17CD", "LAD17NM", "age", "sex"), variable.name = "year", value.name = "pops")
popl <- popl[age < 90]
popl <- rbind(popl, ttm)

# include pop estimates before 2020
file <- "./inputs/pop_estimates_lsoa/LSOA_mid_year_population_estimates.fst"
ol <- read_fst(file, as.data.table = TRUE)
ol[indx, on = "LSOA11CD", `:=` (LAD17CD = i.LAD17CD, LAD17NM = LAD17NM)]
ol <- ol[, lapply(.SD, sum), .SDcols = as.character(0:100),
         keyby = .(year, LAD17CD, LAD17NM, sex)]
ol <- melt(ol, c("year", "LAD17CD", "LAD17NM", "sex"), variable.name = "age", value.name = "pops")
popl[, year := as.integer(as.character(year)) - 2000L]
popl <- rbind(popl[year > 20L], ol)
popl[, LAD17NM := NULL]
popl[indx, on = "LAD17CD", `:=` (LAD17NM = i.LAD17NM)]
popl[, `:=` (LAD17CD = factor(LAD17CD),
             LAD17NM = factor(LAD17NM),
             sex = factor(sex),
             year = as.integer(as.character(year)),
             age = as.integer(as.character(age))
)]
setkey(popl, year, sex, age, LAD17CD)
setcolorder(popl)
all.equal(sort(unique(popl$LAD17NM)), sort(unique(indx$LAD17NM)))
write_fst(popl, "./inputs/pop_projections/lad17_proj.fst", 100)
