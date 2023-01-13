library(data.table)
library(CKutils)
library(ggplot2)
library(ggthemes)
library(scales)

prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)
output_dir <- "/mnt/storage_fast/output/hf_real"
theme_set(new = theme_economist())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))

# TODO standardised age/agegroup in filenames

# Prvl not standardised----

smmrs_prvl_not_standardised <- function(strata = list(
  "year",
  c("year", "sex"),
  c("year", "dimd"),
  c("year", "agegroup"),
  c("year", "agegroup", "sex"),
  c("year", "agegroup", "sex", "dimd")
),
output_dir = "/mnt/storage_fast/output/hf_real",
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

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

# outstrata <- c("mc", "year", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year (not standardised).csv")
#
# outstrata <- c("mc", "year", "sex", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-sex (not standardised).csv")
#
# outstrata <- c("mc", "year", "dimd", "scenario")
# d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
# ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
# d <- melt(d, id.vars = outstrata)
# setkey(d, "variable")
# d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
# d <- d[disease != "popsize"]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-dimd (not standardised).csv")

# Pop estimates ----
outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
d[, disease := NULL]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/pop size by year-agegroup-sex-dimd.csv")

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
d[, disease := NULL]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/pop size by year-agegroup-sex.csv")

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
d[, disease := NULL]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/pop size by year-sex.csv")

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "pop_size_")))
d[, disease := NULL]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/pop size by year.csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence change by year (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence change by year-sex (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence change by year-dimd (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence change by year-agegroup-sex.csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence change by year-agegroup-sex-dimd.csv")


# Prvl standardised----

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_esp.csv.gz"
)[, `:=` (year = year + 2000L,
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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year (age-sex-dimd standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-sex (age-dimd standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-dimd (age-sex standardised).csv")


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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence change by year (age-sex-dimd standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence change by year-sex (age-dimd standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence change by year-dimd (age-sex standardised).csv")


# Incd not standardised----

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year (not standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-sex (not standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-dimd (not standardised).csv")

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-agegroup-sex.csv")

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-agegroup-sex-dimd.csv")

# Incd standardised----

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year (age-sex-dimd standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-sex (age-dimd standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-dimd (age-sex standardised).csv")



# Mrtl not standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year (not standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-sex (not standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-dimd (not standardised).csv")


outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-agegroup-sex.csv")


outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-agegroup-sex-dimd.csv")


# Mrtl standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year (age-sex-dimd standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-sex (age-dimd standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-dimd (age-sex standardised).csv")



# Mean CMS not standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/cms_score_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score by year (not standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score by year-sex (not standardised).csv")

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score by year-agegroup-sex (not standardised).csv")

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score by year-agegroup-sex-dimd (not standardised).csv")

# Mean CMS change with 2019 baseline not standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/cms_score_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score change by year (not standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score change by year-sex (not standardised).csv")

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score change by year-agegroup-sex (not standardised).csv")

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score change by year-agegroup-sex-dimd (not standardised).csv")

# Mean CMS standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/cms_score_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score by year (age-sex-dimd standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score by year-dimd (age-sex standardised).csv")


# Mean CMS change wih 2019 as baseline standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/cms_score_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score change by year (age-sex-dimd standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score change by year-dimd (age-sex standardised).csv")






# Mean CMS by age not standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/cms_score_by_age_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]


outstrata <- c("mc", "year", "age", "sex", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score by year-age-sex (not standardised).csv")

outstrata <- c("mc", "year", "age", "sex", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score by year-age-sex-dimd (not standardised).csv")

# Mean CMS by age change with 2019 baseline not standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/cms_score_by_age_scaled_up.csv.gz"
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
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score change by year-age-sex (not standardised).csv")

outstrata <- c("mc", "year", "age", "sex", "dimd", "scenario")
e <- tt[, .("mean_cms_score" = weighted.mean(cms_score, popsize)),
        keyby = eval(outstrata)]
e19 <- e[year == 2019][, year := NULL]
e[e19, on = setdiff(outstrata, "year"), mean_cms_score := mean_cms_score/i.mean_cms_score]
e <- e[, as.list(fquantile(mean_cms_score, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mean_cms_score_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mean CMS score change by year-age-sex-dimd (not standardised).csv")

# All-cause mortality by disease not standardised----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/all_cause_mrtl_by_dis_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year-sex (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year-agegroup-sex (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year-agegroup-sex-dimd (not standardised).csv")

# All-cause mortality by disease not standardised pop denominator----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/all_cause_mrtl_by_dis_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
pp <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year popdenom (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year-sex popdenom (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year-agegroup-sex popdenom (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year-agegroup-sex-dimd popdenom (not standardised).csv")

rm(pp)

# All-cause mortality by disease standardised----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/all_cause_mrtl_by_dis_esp.csv.gz"
)[, `:=` (year = year + 2000L,
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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year (age-sex-dimd standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year-sex (age-dimd standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year-dimd (age-sex standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/all-cause mrtl by disease-year-sex-dimd (age standardised).csv")
rm(cases)


# Disease characteristics non standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/dis_characteristics_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/disease characteristics by year (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/disease characteristics by year-sex (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/disease characteristics by year-dimd (not standardised).csv")

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
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/disease characteristics by year-sex-dimd (not standardised).csv")

rm(d, tt)
# XPS ----

xps_tab <- fread("/mnt/storage_fast/output/hf_real/xps/xps.csv.gz")

xps_tab <- xps_tab[sex != "All" & agegrp20 != "All" & qimd != "All" &
                       sha == "All" & ethnicity == "All"][
                           , c("ethnicity", "sha") := NULL]


# European standardised population 2013 (esp) weights
tt <- data.table(agegroup = agegrp_name(0, 99),
                 wt_esp  = c(1000, 4000, 5500, 5500, 5500, 6000, 6000, 6500,
                             7000, 7000, 7000, 7000, 6500, 6000, 5500, 5000,
                             4000, 2500, 1500, 800, 200))
esp <- CJ(agegroup = agegrp_name(0, 99, ),
          sex = c("men", "women"),
          qimd = c("1 most deprived", as.character(2:4), "5 least deprived")
)
absorb_dt(esp, tt)

esp[, `:=` (qimd = factor(qimd,
                          levels = c("1 most deprived",
                                     "2", "3", "4",
                                     "5 least deprived")),
            agegrp20 = ifelse(agegroup %in% c("<1", "01-04", "05-09", "10-14", "15-19", "20-24", "25-29"), "<30",
                              ifelse(agegroup %in% c("30-34", "35-39", "40-44", "45-49"), "30-49",
                                     ifelse(agegroup %in% c("50-54" ,"55-59", "60-64" ,"65-69" ), "50-69",
                                            ifelse(agegroup %in% c("70-74", "75-79", "80-84", "85-89" ), "70-89",
                                                   "90+"   )))))]
esp <- esp[, .(wt_esp = sum(wt_esp)), keyby = .(sex, qimd, agegrp20)]
absorb_dt(xps_tab, esp)


xps_names <- grep("_curr_xps$", names(xps_tab), value = TRUE)


outstrata <- c("mc", "year", "scenario")
d <- xps_tab[, lapply(.SD, weighted.mean, wt_esp), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/exposures by year (age-sex-qimd standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
d <- xps_tab[, lapply(.SD, weighted.mean, wt_esp), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/exposures by year-sex (age-qimd standardised).csv")

outstrata <- c("mc", "year", "qimd", "scenario")
d <- xps_tab[, lapply(.SD, weighted.mean, wt_esp), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/exposures by year-qimd (age-sex standardised).csv")

outstrata <- c("mc", "year", "sex", "qimd", "scenario")
d <- xps_tab[, lapply(.SD, weighted.mean, wt_esp), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)
]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/exposures by year-sex-qimd (age standardised).csv")

outstrata <- c("mc", "year", "agegrp20", "sex", "qimd", "scenario")
d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/exposures by year-agegroup-sex-qimd (not standardised).csv")

outstrata <- c("mc", "year", "scenario")
d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/exposures by year (not standardised).csv")

outstrata <- c("mc", "year", "agegrp20", "scenario")
d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/exposures by year-agegroup (not standardised).csv")

outstrata <- c("mc", "year", "agegrp20", "sex", "scenario")
d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/exposures by year-agegroup-sex (not standardised).csv")

outstrata <- c("mc", "year", "qimd", "scenario")
d <- xps_tab[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
setkey(d, "variable")
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/exposures by year-qimd (not standardised).csv")


library(data.table)
ans <- data.table()
for (i in 1:100) {
lc <- fread(paste0("/mnt/storage_fast/output/hf_real/lifecourse/", i,"_lifecourse.csv.gz"), key = c("pid", "year"))[scenario == "sc0"]
nm <- grep("_prvl$", names(lc), value = TRUE)
nm <- grep("^cms", nm, value = TRUE, invert = TRUE)
for (j in nm) {
  set(lc, NULL, j, as.integer(lc[[j]] >= 1L))
  # set(lc, NULL, j, as.character(lc[[j]]))

  lc[, (j) := fifelse(get(j) == 1L, j, "")]
}
lc[, first_dis := paste0(t2dm_prvl, ra_prvl, ckd_prvl, asthma_prvl, helo_prvl, alcpr_prvl,
                         prostate_ca_prvl, ibs_prvl, htn_prvl, copd_prvl, t1dm_prvl,
                         ctd_prvl, dementia_prvl, colorectal_ca_prvl, breast_ca_prvl,
                         lung_ca_prvl, chd_prvl,  dm_prvl, ctdra_prvl,  other_ca_prvl,
                         af_prvl, hf_prvl,  pain_prvl,  cancer_prvl,  stroke_prvl,
                         epilepsy_prvl, andep_prvl, psychosis_prvl, constipation_prvl)]
out <- lc[cmsmm0_prvl == 1, .N, by = first_dis][,  N := N/sum(N)][order(N, decreasing = TRUE), ][1:10]
out[, mc := i]
ans <- rbind(ans, out)
}

View(lc[, .(pid, year, cmsmm0_prvl, cms_count)])
View(lc[cmsmm0_prvl == 1 & first_dis == "", .(pid, year, cmsmm0_prvl, cms_count, t2dm_prvl, ra_prvl, ckd_prvl, asthma_prvl, helo_prvl, alcpr_prvl,
                                              prostate_ca_prvl, ibs_prvl, htn_prvl, copd_prvl, t1dm_prvl,
                                              ctd_prvl, dementia_prvl, colorectal_ca_prvl, breast_ca_prvl,
                                              lung_ca_prvl, chd_prvl,  dm_prvl, ctdra_prvl,  other_ca_prvl,
                                              af_prvl, hf_prvl,  pain_prvl,  cancer_prvl,  stroke_prvl,
                                              epilepsy_prvl, andep_prvl, psychosis_prvl, constipation_prvl)])
View(lc[cmsmm0_prvl == 1 & first_dis == "", ])
