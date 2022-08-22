library(data.table)
library(CKutils)
library(ggplot2)
library(ggthemes)
library(scales)

prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)
theme_set(new = theme_economist())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))

# Year ----
outstrata <- c("mc", "year", "scenario")

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year (not standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year (age-sex-dimd standardised).csv")



tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year (not standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year (age-sex-dimd standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_out.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year (not standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year (age-sex-dimd standardised).csv")

# Sex ----

outstrata <- c("mc", "year", "sex", "scenario")

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-sex (not standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-sex (age-dimd standardised).csv")



tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-sex (not standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-sex (age-dimd standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_out.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-sex (not standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-sex (age-dimd standardised).csv")



# DIMD ----
outstrata <- c("mc", "year", "dimd", "scenario")

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-dimd (not standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-dimd (age-sex standardised).csv")



tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-dimd (not standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-dimd (age-sex standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_out.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-dimd (not standardised).csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-dimd (age-sex standardised).csv")


# Age / sex ----

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-agegroup-sex.csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-agegroup-sex.csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_out.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-agegroup-sex.csv")



# Age / Sex /DIMD ----

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-agegroup-sex-dimd.csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-agegroup-sex-dimd.csv")


tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_out.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
e <- tt[, lapply(.SD, sum), .SDcols = patterns("_mrtl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e[, popsize := NULL]
e <- e[, as.list(fquantile(all_cause_mrtl, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mrtl_rate_")))
setkeyv(e, setdiff(outstrata, "mc"))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/mortality by year-agegroup-sex-dimd.csv")
