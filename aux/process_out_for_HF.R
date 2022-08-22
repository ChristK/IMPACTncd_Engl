library(data.table)
library(CKutils)
library(ggplot2)
library(ggthemes)
library(scales)

prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)
theme_set(new = theme_economist())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))

# Prvl not standardised----

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year (not standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-sex (not standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-dimd (not standardised).csv")

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-agegroup-sex.csv")

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-agegroup-sex-dimd.csv")

# Prvl age standardised----

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year (age-sex-dimd standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-sex (age-dimd standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year-dimd (age-sex standardised).csv")

# Incd not standardised----

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_scaled_up.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year (not standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-sex (not standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-dimd (not standardised).csv")

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-agegroup-sex.csv")

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-agegroup-sex-dimd.csv")

# Incd age standardised----

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_esp.csv.gz"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year (age-sex-dimd standardised).csv")

outstrata <- c("mc", "year", "sex", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-sex (age-dimd standardised).csv")

outstrata <- c("mc", "year", "dimd", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "incd_rate_")))
d <- d[disease != "popsize"]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/incidence by year-dimd (age-sex standardised).csv")



# Mrtl not standardised ----
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_out.csv.gz"
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


# Mrtl age standardised ----
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



# XPS ----
# dt <- fread("/mnt/storage_fast/output/hf_real/xps/xps.csv.gz")

# Median CMS ----

tt <- list.files("/mnt/storage_fast/output/hf_real/lifecourse/", full.names = TRUE)
out <- lapply(tt[1:2], function(x) {
    f <- fread(x)[, `:=` (year = year + 2000L,
                          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
    out <- list()
    out$yr <- f[, .(mc, mean_CMS = weighted.mean(cms_score, wt_esp)), keyby = .(year, scenario)]
    out$dimd <- f[, .(mc, mean_CMS = weighted.mean(cms_score, wt_esp)), keyby = .(year, dimd, scenario)]
    out
})
mean_cms_by_year <- rbindlist(lapply(out, `[[`, 1))[, as.list(fquantile(mean_CMS, prbl)), keyby = .(year, scenario)]
setnames(mean_cms_by_year, c("year", "scenario", percent(prbl, prefix = "mean_cms_")))
setkeyv(mean_cms_by_year, c("year", "scenario"))
fwrite(mean_cms_by_year, "/mnt/storage_fast/output/hf_real/tables/mean cms by year (age-sex-dimd standardised).csv")

mean_cms_by_year_dimd <- rbindlist(lapply(out, `[[`, 2))[, as.list(fquantile(mean_CMS, prbl)), keyby = .(year, dimd, scenario)]
setnames(mean_cms_by_year_dimd, c("year", "dimd", "scenario", percent(prbl, prefix = "mean_cms_")))
setkeyv(mean_cms_by_year_dimd,c("year", "dimd", "scenario"))
fwrite(mean_cms_by_year_dimd, "/mnt/storage_fast/output/hf_real/tables/mean cms by year-dimd (age-sex standardised).csv")
