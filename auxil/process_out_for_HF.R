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

tbl_smmrs <- function(
    what = c("prvl", "prvl_change", "incd", "incd_change",
             "ftlt", "ftlt_change" "pop"),
    type = c("ons", "esp"),
    strata = list(
      "year",
      c("year", "sex"),
      c("year", "dimd"),
      c("year", "agegroup"),
      c("year", "agegroup", "sex"),
      c("year", "agegroup", "sex", "dimd")
    ),
output_dir = output_dir,
prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
baseline_year = 2019L # only used for prvl_change etc.
) {
  strata <- lapply(strata, function(x)
    c("mc", "scenario", x))

    # construct file path to read from summaries
    str0 <- c(
            "prvl" = "prvl",
            "prvl_change" = "prvl",
            "pop" = "prvl",
            "incd" = "incd",
            "incd_change" = "incd",
            "ftlt" = "/dis_mrtl",
            "ftlt_change" = "/dis_mrtl"
    )
    str1 <- c("ons" = "_scaled_up.csv.gz", "esp" = "_esp.csv.gz")
    fpth <- file.path(output_dir, "summaries", paste0(str0[[what]], str1[[type]]))

# other useful strings
    str2 <- c(
            "prvl" = "_prvl$|^popsize$",
            "prvl_change" = "_prvl$|^popsize$",
            "pop" = "^popsize$",
            "incd" = "_incd$|^popsize$",
            "incd_change" = "_incd$|^popsize$",
            "ftlt" = "_deaths$|_prvl$",
            "ftlt_change" = "_deaths$|_prvl$"
    ) # used in grep
    str3 <- c(
            "prvl" = "prvl_rate_",
            "prvl_change" = "prct_change_",
            "pop" = "pop_size_",
            "incd" = "incd_rate_",
            "incd_change" = "prct_change_",
            "ftlt" = "ftlt_rate_",
            "ftlt_change" = "ftlt_rate_"
    ) # used to col name output
    str4 <- c(
            "prvl" = "prevalence by ",
            "prvl_change" = "prevalence change by ",
            "pop" = "pop size by ",
            "incd" = "incidence by ",
            "incd_change" = "incidence change by ",
            "ftlt" = "fatality by ",
            "ftlt_change" = "fatality change by "
    )
    str5 <- c(
     "ons" = " (not standardised).csv",
     "esp" = paste0(" (", paste(setdiff(c("age", x), c("mc", "scenario", "year")),
      collapse = "-"), " standardised).csv")
     )
    str6 <- paste0(
        str4[[what]],
        paste(setdiff(x, c("mc", "scenario")), collapse = "-"),
        str5[[type]]
      ) # used for output file name/path


  tt <-
    fread(fpth)[, `:=` (year = year + 2000L,
                    dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

  if (grepl("^ftlt"), what) {
    fpth <- file.path(output_dir, "summaries",  paste0(str0[["prvl"]], str1[[type]]))

   t1 <- fread(fpth)[, `:=` (year = year + 2000L,
                    dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
setnames(t1, "popsize", "nonmodelled_prvl")
absorb_dt(tt, t1)
  }

  lapply(strata, function(x) {
    d <-
      tt[, lapply(.SD, sum), .SDcols = patterns(str2[[what]]),
         keyby = eval(x)][, lapply(.SD, function(x)
        x / popsize), keyby = x]
    d <- melt(d, id.vars = x)

    if (grepl("_change$", what)) { # when calculating change
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year"), "variable"), value := value/i.value]
    }

    setkey(d, "variable")
    d <-
      d[, fquantile_byid(value, prbl, id = as.character(variable)),
        keyby = eval(setdiff(x, "mc"))]
    setnames(d, c(
      setdiff(x, "mc"),
      "disease",
      percent(prbl, prefix = str3[[what]])
    ))
    if (what == "pop") {
      d[, disease := NULL]
    } else {
      d <- d[disease != "popsize"]
    }
    setkeyv(d, setdiff(x, "mc"))
    fwrite(d, file.path(
      output_dir, "tables", str6
    ))
  })
}

tbl_smmrs("prvl", "ons", list("year",
                                 c("year", "sex"),
                                 c("year", "dimd"),
                                 c("year", "agegrp"),
                                 c("year", "agegrp", "sex"),
                                 c("year", "agegrp", "sex", "dimd")), output_dir)
tbl_smmrs("prvl_change", "ons", list("year",
                                 c("year", "sex"),
                                 c("year", "dimd"),
                                 c("year", "agegrp"),
                                 c("year", "agegrp", "sex"),
                                 c("year", "agegrp", "sex", "dimd")), output_dir)
tbl_smmrs("incd", "ons", list("year",
                                 c("year", "sex"),
                                 c("year", "dimd"),
                                 c("year", "agegrp"),
                                 c("year", "agegrp", "sex"),
                                 c("year", "agegrp", "sex", "dimd")), output_dir)
tbl_smmrs("pop", "ons", list("year",
                                 c("year", "sex"),
                                 c("year", "dimd"),
                                 c("year", "agegrp"),
                                 c("year", "agegrp", "sex"),
                                 c("year", "agegrp", "sex", "dimd")), output_dir)
tbl_smmrs("prvl", "esp", list("year",
                                 c("year", "sex"),
                                 c("year", "dimd"),
                                 c("year", "sex", "dimd")), output_dir)
tbl_smmrs("prvl_change", "esp", list("year",
                                 c("year", "sex"),
                                 c("year", "dimd"),
                                 c("year", "sex", "dimd")), output_dir)
tbl_smmrs("incd", "esp", list("year",
                                 c("year", "sex"),
                                 c("year", "dimd"),
                                 c("year", "sex", "dimd")), output_dir)


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

# tt <- fread(paste0(sSummariesSubDirPath,"prvl_scaled_up.csv.gz")
# )[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]



# Pop estimates ----
# outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
# d[, disease := NULL]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"pop size by year-agegroup-sex-dimd.csv"))

# outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
# d[, disease := NULL]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"pop size by year-agegroup-sex.csv"))

# outstrata <- c("mc", "year", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
# d[, disease := NULL]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"pop size by year-sex.csv"))

# outstrata <- c("mc", "year", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
# d[, disease := NULL]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"pop size by year.csv"))

# Prvl change with 2019 baseline ----

# outstrata <- c("mc", "year", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# d19 <- d[year == 2019][, year := NULL]
# d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year (not standardised).csv"))

# outstrata <- c("mc", "year", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# d19 <- d[year == 2019][, year := NULL]
# d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-sex (not standardised).csv"))

# outstrata <- c("mc", "year", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# d19 <- d[year == 2019][, year := NULL]
# d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-dimd (not standardised).csv"))

# outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# d19 <- d[year == 2019][, year := NULL]
# d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-agegroup-sex.csv"))

# outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# d19 <- d[year == 2019][, year := NULL]
# d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-agegroup-sex-dimd.csv"))


# Prvl standardised----

# tt <- fread(paste0(sSummariesSubDirPath,"prvl_esp.csv.gz"))[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

# outstrata <- c("mc", "year", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence by year (age-sex-dimd standardised).csv"))

# outstrata <- c("mc", "year", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence by year-sex (age-dimd standardised).csv"))

# outstrata <- c("mc", "year", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence by year-dimd (age-sex standardised).csv"))


# Prvl change with 2019 baseline age standardised----

# outstrata <- c("mc", "year", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# d19 <- d[year == 2019][, year := NULL]
# d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year (age-sex-dimd standardised).csv"))

# outstrata <- c("mc", "year", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# d19 <- d[year == 2019][, year := NULL]
# d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-sex (age-dimd standardised).csv"))

# outstrata <- c("mc", "year", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# d19 <- d[year == 2019][, year := NULL]
# d[d19, on = c(setdiff(outstrata, "year"), "variable"), value := value/i.value]
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"prevalence change by year-dimd (age-sex standardised).csv"))



# Incd not standardised----

# tt <- fread(paste0(sSummariesSubDirPath,"incd_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

# outstrata <- c("mc", "year", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"incidence by year (not standardised).csv"))

# outstrata <- c("mc", "year", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"incidence by year-sex (not standardised).csv"))

# outstrata <- c("mc", "year", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"incidence by year-dimd (not standardised).csv"))

# outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"incidence by year-agegroup-sex.csv"))

# outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"incidence by year-agegroup-sex-dimd.csv"))

# Incd standardised----

# tt <- fread(paste0(sSummariesSubDirPath,"incd_esp.csv.gz"))[, `:=` (year = year + 2000L,
#           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

# outstrata <- c("mc", "year", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"incidence by year (age-sex-dimd standardised).csv"))

# outstrata <- c("mc", "year", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"incidence by year-sex (age-dimd standardised).csv"))

# outstrata <- c("mc", "year", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_incd$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"incidence by year-dimd (age-sex standardised).csv"))



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
tt <- fread(paste0(sSummariesSubDirPath,"cms_score_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-sex (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-agegroup-sex (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-agegroup-sex-dimd (not standardised).csv"))

# Mean CMS change with 2019 baseline not standardised ----
tt <- fread(paste0(sSummariesSubDirPath,"cms_score_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-sex (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-agegroup-sex (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-agegroup-sex-dimd (not standardised).csv"))

# Mean CMS standardised ----
tt <- fread(paste0(sSummariesSubDirPath,"cms_score_esp.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year (age-sex-dimd standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-dimd (age-sex standardised).csv"))


# Mean CMS change wih 2019 as baseline standardised ----
tt <- fread(paste0(sSummariesSubDirPath,"cms_score_esp.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year (age-sex-dimd standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-dimd (age-sex standardised).csv"))






# Mean CMS by age not standardised ----
tt <- fread(paste0(sSummariesSubDirPath,"cms_score_by_age_scaled_up.csv.gz")
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]


outstrata <- c("mc", "year", "age", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-age-sex (not standardised).csv"))

outstrata <- c("mc", "year", "age", "sex", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score by year-age-sex-dimd (not standardised).csv"))

# Mean CMS by age change with 2019 baseline not standardised ----
tt <- fread(paste0(sSummariesSubDirPath,"cms_score_by_age_scaled_up.csv.gz")
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]



outstrata <- c("mc", "year", "age", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-age-sex (not standardised).csv"))

outstrata <- c("mc", "year", "age", "sex", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, paste0(sTablesSubDirPath,"mean CMS score change by year-age-sex-dimd (not standardised).csv"))

# All-cause mortality by disease not standardised----
tt <- fread(paste0(sSummariesSubDirPath,"all_cause_mrtl_by_dis_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-sex (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-agegroup-sex (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-agegroup-sex-dimd (not standardised).csv"))

# All-cause mortality by disease not standardised pop denominator----
tt <- fread(paste0(sSummariesSubDirPath,"all_cause_mrtl_by_dis_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
pp <- fread(paste0(sSummariesSubDirPath,"prvl_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
cases <- pp[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = outstrata, value := value/popsize]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year popdenom (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
cases <- pp[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = outstrata, value := value/popsize]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-sex popdenom (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
cases <- pp[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = outstrata, value := value/popsize]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-agegroup-sex popdenom (not standardised).csv"))

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
cases <- pp[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = outstrata, value := value/popsize]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-agegroup-sex-dimd popdenom (not standardised).csv"))

rm(pp)

# All-cause mortality by disease standardised----
tt <- fread(paste0(sSummariesSubDirPath,"all_cause_mrtl_by_dis_esp.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year (age-sex-dimd standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-sex (age-dimd standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-dimd (age-sex standardised).csv"))

outstrata <- c("mc", "year", "sex", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-sex-dimd (age standardised).csv"))
rm(cases)


# Disease characteristics non standardised ----
tt <- fread(paste0(sSummariesSubDirPath,"dis_characteristics_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")),
          mean_cms_count_cms1st_cont = as.numeric(mean_cms_count_cms1st_cont))]
d1 <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|dimd|^cases_")]
d1 <- melt(d1, id.vars = c("mc", "year", "scenario", "sex", "dimd"))
d1 <- unique(d1, by = c("mc", "year", "scenario", "sex", "dimd", "variable"))
d1[, `:=` (disease = gsub("^cases_", "", variable), variable = NULL)]
tt <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|dimd|^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_")]
tt <- melt(tt, id.vars = c("mc", "year", "scenario", "sex", "dimd"))
tt[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
tt[d1, on = c("mc", "year", "scenario", "sex", "dimd", "disease"), cases := i.value]
rm(d1) # NOTE mean_age_incd contains NAs

outstrata <- c("mc", "year", "scenario")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")] # na.rm = TRUE for mean_age_incd
setkey(d, "variable")
d <- d[, fquantile_byid(V1, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "variable", percent(prbl, prefix = "value_")))
d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
d[grep("^mean_duration_", variable), type := "mean_duration"]
d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
d[, variable := NULL]
setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
setcolorder(d)
fwrite(d, paste0(sTablesSubDirPath,"disease characteristics by year (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")]
setkey(d, "variable")
d <- d[, fquantile_byid(V1, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "variable", percent(prbl, prefix = "value_")))
d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
d[grep("^mean_duration_", variable), type := "mean_duration"]
d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
d[, variable := NULL]
setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
setcolorder(d)
fwrite(d, paste0(sTablesSubDirPath,"disease characteristics by year-sex (not standardised).csv"))

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")]
setkey(d, "variable")
d <- d[, fquantile_byid(V1, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "variable", percent(prbl, prefix = "value_")))
d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
d[grep("^mean_duration_", variable), type := "mean_duration"]
d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
d[, variable := NULL]
setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
setcolorder(d)
fwrite(d, paste0(sTablesSubDirPath,"disease characteristics by year-dimd (not standardised).csv"))

outstrata <- c("mc", "year", "sex", "dimd", "scenario")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")]
setkey(d, "variable")
d <- d[, fquantile_byid(V1, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "variable", percent(prbl, prefix = "value_")))
d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
d[grep("^mean_duration_", variable), type := "mean_duration"]
d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
d[, variable := NULL]
setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
setcolorder(d)
fwrite(d, paste0(sTablesSubDirPath,"disease characteristics by year-sex-dimd (not standardised).csv"))

rm(d, tt)

# XPS ----
# xps_tab <- fread(file.path(simulationParameters$output_dir,"xps/xps20.csv.gz"))
#
# xps_tab <- xps_tab[sex != "All" & agegrp20 != "All" & qimd != "All"]
# xps_names <- grep("_curr_xps$", names(xps_tab), value = TRUE)
#
# outstrata <- c("mc", "year", "agegrp20", "sex", "qimd", "scenario")
# d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"exposures by year-agegroup-sex-qimd (not standardised).csv"))
#
# outstrata <- c("mc", "year", "scenario")
# d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"exposures by year (not standardised).csv"))
#
# outstrata <- c("mc", "year", "agegrp20", "scenario")
# d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"exposures by year-agegroup (not standardised).csv"))
#
# outstrata <- c("mc", "year", "agegrp20", "sex", "scenario")
# d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"exposures by year-agegroup-sex (not standardised).csv"))
#
# outstrata <- c("mc", "year", "qimd", "scenario")
# d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"exposures by year-qimd (not standardised).csv"))
#
#
#
# # European standardised population 2013 (esp) weights
# xps_tab <- fread(file.path(simulationParameters$output_dir,"xps/xps5.csv.gz"))
#
# xps_tab <- xps_tab[sex != "All" & agegrp5 != "All" & qimd != "All"]
# xps_names <- grep("_curr_xps$", names(xps_tab), value = TRUE)
#
# tt <- data.table(agegrp5 = agegrp_name(0, 99),
#                  wt_esp  = c(1000, 4000, 5500, 5500, 5500, 6000, 6000, 6500,
#                              7000, 7000, 7000, 7000, 6500, 6000, 5500, 5000,
#                              4000, 2500, 1500, 800, 200))
# esp <- CJ(agegrp5 = agegrp_name(0, 99),
#           sex = c("men", "women"),
#           qimd = factor(c("1 most deprived", as.character(2:4), "5 least deprived")))
#
# absorb_dt(esp, tt)
#
# # esp[, `:=` (qimd = factor(qimd,
# #                           levels = c("1 most deprived",
# #                                      "2", "3", "4",
# #                                      "5 least deprived")),
# #             agegrp20 = ifelse(agegroup %in% c("<1", "01-04", "05-09", "10-14", "15-19", "20-24", "25-29"), "<30",
# #                               ifelse(agegroup %in% c("30-34", "35-39", "40-44", "45-49"), "30-49",
# #                                      ifelse(agegroup %in% c("50-54" ,"55-59", "60-64" ,"65-69" ), "50-69",
# #                                             ifelse(agegroup %in% c("70-74", "75-79", "80-84", "85-89" ), "70-89",
# #                                                    "90+"   )))))]
# # esp <- esp[, .(wt_esp = sum(wt_esp)), keyby = .(sex, qimd, agegrp20)]
# absorb_dt(xps_tab, esp)
#
#
#
#
# outstrata <- c("mc", "year", "scenario")
# d <- xps_tab[, lapply(.SD, weighted.mean, wt_esp), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)
# ]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"exposures by year (age-sex-qimd standardised).csv"))
#
# outstrata <- c("mc", "year", "sex", "scenario")
# d <- xps_tab[, lapply(.SD, weighted.mean, wt_esp), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)
# ]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"exposures by year-sex (age-qimd standardised).csv"))
#
# outstrata <- c("mc", "year", "qimd", "scenario")
# d <- xps_tab[, lapply(.SD, weighted.mean, wt_esp), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)
# ]
# d <- melt(d, id.vars = outstrata)
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"exposures by year-qimd (age-sex standardised).csv"))
#
# outstrata <- c("mc", "year", "sex", "qimd", "scenario")
# d <- xps_tab[, lapply(.SD, weighted.mean, wt_esp), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)
# ]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"exposures by year-sex-qimd (age standardised).csv"))



# library(data.table)
# ans <- data.table()
# for (i in 1:100) {
# lc <- fread(paste0("/mnt/storage_fast/output/hf_real/lifecourse/", i,"_lifecourse.csv.gz"), key = c("pid", "year"))[scenario == "sc0"]
# nm <- grep("_prvl$", names(lc), value = TRUE)
# nm <- grep("^cms", nm, value = TRUE, invert = TRUE)
# for (j in nm) {
#   set(lc, NULL, j, as.integer(lc[[j]] >= 1L))
#   # set(lc, NULL, j, as.character(lc[[j]]))
#
#   lc[, (j) := fifelse(get(j) == 1L, j, "")]
# }
# lc[, first_dis := paste0(t2dm_prvl, ra_prvl, ckd_prvl, asthma_prvl, helo_prvl, alcpr_prvl,
#                          prostate_ca_prvl, ibs_prvl, htn_prvl, copd_prvl, t1dm_prvl,
#                          ctd_prvl, dementia_prvl, colorectal_ca_prvl, breast_ca_prvl,
#                          lung_ca_prvl, chd_prvl,  dm_prvl, ctdra_prvl,  other_ca_prvl,
#                          af_prvl, hf_prvl,  pain_prvl,  cancer_prvl,  stroke_prvl,
#                          epilepsy_prvl, andep_prvl, psychosis_prvl, constipation_prvl)]
# out <- lc[cmsmm0_prvl == 1, .N, by = first_dis][,  N := N/sum(N)][order(N, decreasing = TRUE), ][1:10]
# out[, mc := i]
# ans <- rbind(ans, out)
# }
#
# View(lc[, .(pid, year, cmsmm0_prvl, cms_count)])
# View(lc[cmsmm0_prvl == 1 & first_dis == "", .(pid, year, cmsmm0_prvl, cms_count, t2dm_prvl, ra_prvl, ckd_prvl, asthma_prvl, helo_prvl, alcpr_prvl,
#                                               prostate_ca_prvl, ibs_prvl, htn_prvl, copd_prvl, t1dm_prvl,
#                                               ctd_prvl, dementia_prvl, colorectal_ca_prvl, breast_ca_prvl,
#                                               lung_ca_prvl, chd_prvl,  dm_prvl, ctdra_prvl,  other_ca_prvl,
#                                               af_prvl, hf_prvl,  pain_prvl,  cancer_prvl,  stroke_prvl,
#                                               epilepsy_prvl, andep_prvl, psychosis_prvl, constipation_prvl)])
# View(lc[cmsmm0_prvl == 1 & first_dis == "", ])
