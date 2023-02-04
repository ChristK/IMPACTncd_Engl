library(data.table)
library(yaml)
library(CKutils)
library(ggplot2)
library(ggthemes)
library(scales)

prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)
theme_set(new = theme_economist())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))

# TODO standardised age/agegroup in filenames

# Prvl not standardised----
simulationParameters <- read_yaml(base::normalizePath("./inputs/sim_design.yaml", mustWork=TRUE))
sSummariesSubDirPath <- file.path(simulationParameters$output_dir,"summaries/")
sTablesSubDirPath <- file.path(simulationParameters$output_dir,"tables/")
output_dir <- simulationParameters$output_dir

smmrs_prvl_not_standardised <- function(
    strata = list(
      "year",
      c("year", "sex"),
      c("year", "dimd"),
      c("year", "agegroup"),
      c("year", "agegroup", "sex"),
      c("year", "agegroup", "sex", "dimd")
    ),
output_dir = output_dir,
prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)) {
  strata <- lapply(strata, function(x)
    c("mc", "scenario", x))

  tt <-
    fread(file.path(output_dir, "summaries/prvl_scaled_up.csv.gz")
          )[, `:=` (year = year + 2000L,
                    dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

  lapply(strata, function(x) {
    d <-
      tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"),
         keyby = eval(x)][, lapply(.SD, function(x)
        x / popsize), keyby = x]
    d <- melt(d, id.vars = x)
    setkey(d, "variable")
    d <-
      d[, fquantile_byid(value, prbl, id = as.character(variable)),
        keyby = eval(setdiff(x, "mc"))]
    setnames(d, c(
      setdiff(x, "mc"),
      "disease",
      percent(prbl, prefix = "prvl_rate_")
    ))
    d <- d[disease != "popsize"]
    setkeyv(d, setdiff(x, "mc"))
    fwrite(d, file.path(
      output_dir,
      paste0(
        "tables/prevalence by ",
        paste(setdiff(x, c(
          "mc", "scenario"
        )), collapse = "-"),
        " (not standardised).csv"
      )
    ))
  })
}

smmrs_prvl_not_standardised(list("year",
                                 c("year", "sex"),
                                 c("year", "dimd"),
                                 c("year", "agegrp"),
                                 c("year", "agegrp", "sex"),
                                 c("year", "agegrp", "sex", "dimd")), output_dir)

tt <- fread(paste0(sSummariesSubDirPath,"prvl_scaled_up.csv.gz")
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

#outstrata <- c("mc", "year", "scenario")
#d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
#][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
#d <- melt(d, id.vars = outstrata)
#setkey(d, "variable")
#d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
#setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
#d <- d[disease != "popsize"]
#setkeyv(d, setdiff(outstrata, "mc"))
#fwrite(d, paste0(sTablesSubDirPath,"prevalence by year (not standardised).csv"))
#outstrata <- c("mc", "year", "sex", "scenario")
#d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
#][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
#d <- melt(d, id.vars = outstrata)
#setkey(d, "variable")
#d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
#setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
#d <- d[disease != "popsize"]
#setkeyv(d, setdiff(outstrata, "mc"))
#fwrite(d, paste0(sTablesSubDirPath,"prevalence by year-sex (not standardised).csv"))
#
#outstrata <- c("mc", "year", "dimd", "scenario")
#d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
#][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
#d <- melt(d, id.vars = outstrata)
#setkey(d, "variable")
#d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
#setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
#d <- d[disease != "popsize"]
#setkeyv(d, setdiff(outstrata, "mc"))
#fwrite(d, paste0(sTablesSubDirPath,"prevalence by year-dimd (not standardised).csv"))

# Pop estimates ----
outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
d[, disease := NULL]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"pop size by year-agegroup-sex-dimd.csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
d[, disease := NULL]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"pop size by year-agegroup-sex.csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
d[, disease := NULL]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"pop size by year-sex.csv"))

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
d[, disease := NULL]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"pop size by year.csv"))

# Prvl change with 2019 baseline ----

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d19 <- d[year == 2019][, year := NULL]
d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d19 <- d[year == 2019][, year := NULL]
d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-sex (not standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d19 <- d[year == 2019][, year := NULL]
d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-dimd (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d19 <- d[year == 2019][, year := NULL]
d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-agegroup-sex.csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d19 <- d[year == 2019][, year := NULL]
d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-agegroup-sex-dimd.csv"))


# Prvl standardised----

tt <- fread(paste0(sSummariesSubDirPath,"prvl_esp.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence by year (age-sex-dimd standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence by year-sex (age-dimd standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence by year-dimd (age-sex standardised).csv"))


# Prvl change with 2019 baseline age standardised----

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d19 <- d[year == 2019][, year := NULL]
d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year (age-sex-dimd standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d19 <- d[year == 2019][, year := NULL]
d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-sex (age-dimd standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d19 <- d[year == 2019][, year := NULL]
d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-dimd (age-sex standardised).csv"))


# Incd not standardised----

tt <- fread(paste0(sSummariesSubDirPath,"incd_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"incidence by year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"incidence by year-sex (not standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"incidence by year-dimd (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"incidence by year-agegroup-sex.csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"incidence by year-agegroup-sex-dimd.csv"))

# Incd standardised----

tt <- fread(paste0(sSummariesSubDirPath,"incd_esp.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"incidence by year (age-sex-dimd standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"incidence by year-sex (age-dimd standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"incidence by year-dimd (age-sex standardised).csv"))



# Mrtl not standardised ----
tt <- fread(paste0(sSummariesSubDirPath,"mrtl_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mortality by year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mortality by year-sex (not standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mortality by year-dimd (not standardised).csv"))


outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mortality by year-agegroup-sex.csv"))


outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mortality by year-agegroup-sex-dimd.csv"))


# Mrtl standardised ----
tt <- fread(paste0(sSummariesSubDirPath,"mrtl_esp.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mortality by year (age-sex-dimd standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mortality by year-sex (age-dimd standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mortality by year-dimd (age-sex standardised).csv"))



# Mean CMS not standardised ----
# tt <- fread(paste0(sSummariesSubDirPath,"cms_score_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

# outstrata <- c("mc", "year", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year (not standardised).csv"))

# outstrata <- c("mc", "year", "sex", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-sex (not standardised).csv"))

# outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-agegroup-sex (not standardised).csv"))

# outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-agegroup-sex-dimd (not standardised).csv"))

# # Mean CMS change with 2019 baseline not standardised ----
# tt <- fread(paste0(sSummariesSubDirPath,"cms_score_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

# outstrata <- c("mc", "year", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e19 <- e[year == 2019][, year := NULL]
# e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year (not standardised).csv"))

# outstrata <- c("mc", "year", "sex", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e19 <- e[year == 2019][, year := NULL]
# e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-sex (not standardised).csv"))

# outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e19 <- e[year == 2019][, year := NULL]
# e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-agegroup-sex (not standardised).csv"))

# outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e19 <- e[year == 2019][, year := NULL]
# e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-agegroup-sex-dimd (not standardised).csv"))

# # Mean CMS standardised ----
# tt <- fread(paste0(sSummariesSubDirPath,"cms_score_esp.csv.gz"))[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

# outstrata <- c("mc", "year", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year (age-sex-dimd standardised).csv"))

# outstrata <- c("mc", "year", "dimd", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-dimd (age-sex standardised).csv"))


# # Mean CMS change wih 2019 as baseline standardised ----
# tt <- fread(paste0(sSummariesSubDirPath,"cms_score_esp.csv.gz"))[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

# outstrata <- c("mc", "year", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e19 <- e[year == 2019][, year := NULL]
# e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year (age-sex-dimd standardised).csv"))

# outstrata <- c("mc", "year", "dimd", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e19 <- e[year == 2019][, year := NULL]
# e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-dimd (age-sex standardised).csv"))






# Mean CMS by age not standardised ----
# tt <- fread(paste0(sSummariesSubDirPath,"cms_score_by_age_scaled_up.csv.gz")
# )[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]


# outstrata <- c("mc", "year", "age", "sex", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-age-sex (not standardised).csv"))

# outstrata <- c("mc", "year", "age", "sex", "dimd", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-age-sex-dimd (not standardised).csv"))

# # Mean CMS by age change with 2019 baseline not standardised ----
# tt <- fread(paste0(sSummariesSubDirPath,"cms_score_by_age_scaled_up.csv.gz")
# )[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]



# outstrata <- c("mc", "year", "age", "sex", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e19 <- e[year == 2019][, year := NULL]
# e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-age-sex (not standardised).csv"))

# outstrata <- c("mc", "year", "age", "sex", "dimd", "scenario")
# e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
#         keyby = eval(outstrata)]
# e19 <- e[year == 2019][, year := NULL]
# e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
# e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
# setkeyv(e, setdiff(outstrata, "mc"))
# fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-age-sex-dimd (not standardised).csv"))


# Ftlt standardised  ----
# Necessary for validation.
# Note that the denominator is prevalence
tt <- fread(paste0(sSummariesSubDirPath,"/dis_mrtl_esp.csv.gz")
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c(
            "1 most deprived", as.character(2:9), "10 least deprived"
          )))]
t1 <- fread(paste0(sSummariesSubDirPath,"prvl_esp.csv.gz") # for denom
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived"))
)]
setnames(t1, "popsize", "nonmodelled_prvl")
absorb_dt(tt, t1)

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_deaths$|_prvl$"),
        keyby = eval(outstrata)]
nm <- names(d)
nm <- grep("_deaths$", nm, value = TRUE)
nm <- gsub("_deaths$", "", nm)
nm <- setdiff(nm, "alive")
for (i in nm) {
  set(d, NULL, paste0(i, "_ftlt"),
      d[[paste0(i, "_deaths")]] / d[[paste0(i, "_prvl")]])
}
nm <- names(d)
nm <- grep("_deaths$|_prvl$", nm, value = TRUE)
d[, (nm) := NULL]
setnafill(d, "const", 0, cols = grep("_ftlt$", names(d), value = TRUE))
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "ftlt_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"fatality by year (age-sex-dimd standardised).csv"))







