library(data.table)
library(CKutils)
library(ggplot2)
library(ggthemes)
library(scales)

prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)
theme_set(new = theme_economist())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))



tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
fwrite(d, "/mnt/storage_fast/output/hf_real/tables/prevalence by year.csv")

d[disease == "cmsmm0_prvl", disease := "CMS score > 0"]
d[disease == "CMS score > 0", `:=` (
  `prvl_rate_50.0%` = 1 - `prvl_rate_50.0%`,
  `prvl_rate_2.5%` = 1 - `prvl_rate_97.5%`,
  `prvl_rate_97.5%` = 1 - `prvl_rate_2.5%`,
  `prvl_rate_10.0%` = 1 - `prvl_rate_90.0%`,
  `prvl_rate_90.0%` = 1 - `prvl_rate_10.0%`
)]
d[disease == "cmsmm1_prvl", disease := "CMS score > 1"]
d[disease == "cmsmm1.5_prvl", disease := "CMS score > 1.5"]
d[disease == "cmsmm2_prvl", disease := "CMS score > 2"]
d <- d[between(year, 2014, 2040) & disease %in% c("CMS score > 1", "CMS score > 1.5", "CMS score > 2", "CMS score = 0")]
d[scenario == "sc0", scenario := "Base-case"]
d[scenario == "sc1", scenario := "10% BMI reduction"]

ggplot(d[scenario == "Base-case"], aes(x = year, y = `prvl_rate_50.0%`, ymin = `prvl_rate_2.5%`,
              ymax = `prvl_rate_97.5%`, colour = disease, fill = disease)) +
  geom_ribbon(alpha = 2/5, colour = NA) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Prevalence rate", labels = percent) +
  ggtitle(paste0("Prevalence of ", "multimorbidity")) +
  expand_limits(y = 0) +
  theme(legend.title = element_blank())
ggsave("~/pCloudDrive/pCloud Sync/MM prevalence.png", scale = 1.5, width = 16/2, height = 9/2)


ggplot(d[disease == "CMS score > 1.5"], aes(x = year, y = `prvl_rate_50.0%`, ymin = `prvl_rate_2.5%`,
                                 ymax = `prvl_rate_97.5%`, colour = scenario, fill = scenario)) +
  geom_ribbon(alpha = 2/5, colour = NA) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Prevalence rate", labels = percent) +
  ggtitle(paste0("Prevalence of ", "multimorbidity (CMS score > 1.5)")) +
  expand_limits(y = 0) +
  theme(legend.title = element_blank())
ggsave("~/pCloudDrive/pCloud Sync/MM scenario.png", scale = 1.5, width = 16/2, height = 9/2)

d0 <- d[disease == "CMS score > 1",
        .(disease = "CMS score [0,1]",
          year,
          `prvl_rate_50.0%` = 1 - `prvl_rate_50.0%`,
          `prvl_rate_2.5%`  = 1 - `prvl_rate_97.5%`,
          `prvl_rate_97.5%` = 1 - `prvl_rate_2.5%`,
          `prvl_rate_10.0%` = 1 - `prvl_rate_90.0%`,
          `prvl_rate_90.0%` = 1 - `prvl_rate_10.0%`)
]

d1 <- d[disease == "CMS score > 1"]
d1.5 <- d[disease == "CMS score > 1.5"]
d2 <- d[disease == "CMS score > 2"]
d1 <- d1[d1.5, on = "year",  .(disease = "CMS score (1,1.5]",
                               year,
                               `prvl_rate_50.0%` = `prvl_rate_50.0%` - `i.prvl_rate_50.0%`,
                               `prvl_rate_2.5%`  = `prvl_rate_97.5%` - `i.prvl_rate_97.5%`,
                               `prvl_rate_97.5%` = `prvl_rate_2.5%` - `i.prvl_rate_2.5%`,
                               `prvl_rate_10.0%` = `prvl_rate_90.0%` - `i.prvl_rate_90.0%`,
                               `prvl_rate_90.0%` = `prvl_rate_10.0%` - `i.prvl_rate_10.0%`)]
d1.5 <- d1.5[d2, on = "year",  .(disease = "CMS score (1.5,2]",
                                 year,
                                 `prvl_rate_50.0%` = `prvl_rate_50.0%` - `i.prvl_rate_50.0%`,
                                 `prvl_rate_2.5%`  = `prvl_rate_97.5%` - `i.prvl_rate_97.5%`,
                                 `prvl_rate_97.5%` = `prvl_rate_2.5%` - `i.prvl_rate_2.5%`,
                                 `prvl_rate_10.0%` = `prvl_rate_90.0%` - `i.prvl_rate_90.0%`,
                                 `prvl_rate_90.0%` = `prvl_rate_10.0%` - `i.prvl_rate_10.0%`)]
d2[, disease := "CMS score > 2"]
d <- rbind(d0, d1, d1.5, d2)
d[, disease := factor(disease, c("CMS score > 2", "CMS score (1.5,2]", "CMS score (1,1.5]", "CMS score [0,1]"))]


ggplot(d, aes(x = year, y = `prvl_rate_50.0%`, fill = disease)) +
  geom_area(alpha = 0.6 , size = 1, colour = "black") +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Proportion", labels = percent) +
  ggtitle(paste0("CMS distribution")) +
  expand_limits(y = 0) +
  theme(legend.title = element_blank())

ggsave("~/pCloudDrive/pCloud Sync/CMS distr.png", scale = 1.5, width = 16/2, height = 9/2)

tt_esp <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_esp.csv"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt_esp[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e <- melt(e, id.vars = outstrata)
e <- e[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_esp_")))
fwrite(e, "/mnt/storage_fast/output/hf_real/tables/age-standardised prevalence by year and dimd.csv")

e[disease == "cmsmm0_prvl", disease := "CMS score > 0"]
e[disease == "cmsmm1_prvl", disease := "CMS score > 1"]
e[disease == "cmsmm1.5_prvl", disease := "CMS score > 1.5"]
e[disease == "cmsmm2_prvl", disease := "CMS score > 2"]
e <- e[between(year, 2014, 2040) & disease %in% c("CMS score > 1", "CMS score > 1.5", "CMS score > 2")]
e[scenario == "sc0", scenario := "Base-case"]
e[scenario == "sc1", scenario := "10% BMI reduction"]

ggplot(e[scenario == "Base-case" & disease %in% c("CMS score > 1", "CMS score > 1.5", "CMS score > 2")],
       aes(x = year, y = `prvl_rate_esp_50.0%`, ymin = `prvl_rate_esp_2.5%`,
           ymax = `prvl_rate_esp_97.5%`, colour = dimd, fill = dimd)) +
  # geom_ribbon(alpha = 2/5, colour = NA) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Prevalence rate", labels = percent) +
  ggtitle(paste0("Prevalence of ", "multimorbidity", " by dimd (age-standardised)"))+
  expand_limits(y = 0) +
  facet_grid(.~disease)
ggsave("~/pCloudDrive/pCloud Sync/MM prevalence by dimd.png", scale = 1.5, width = 16/2, height = 9/2)

e[year %in% c(2020, 2040) & dimd %in% c("1 most deprived", "10 least deprived") & disease == "CMS score > 1.5" & scenario == "Base-case"]

ggplot(e[disease %in% c("CMS score > 1.5")],
       aes(x = year, y = `prvl_rate_esp_50.0%`, ymin = `prvl_rate_esp_2.5%`,
           ymax = `prvl_rate_esp_97.5%`, colour = dimd, fill = dimd)) +
  # geom_ribbon(alpha = 2/5, colour = NA) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Prevalence rate", labels = percent) +
  ggtitle(paste0("Prevalence of ", "multimorbidity", " by dimd (age-standardised, CMS score > 1.5)"))+
  expand_limits(y = 0) +
  facet_grid(.~ scenario)
ggsave("~/pCloudDrive/pCloud Sync/MM prevalence by dimd scenario.png", scale = 1.5, width = 16/2, height = 9/2)

e[year %in% c(2020, 2040) & dimd %in% c("1 most deprived", "10 least deprived") & disease == "CMS score > 1.5" & scenario == "10% BMI reduction"]


fl <- list.files("/mnt/storage_fast/output/hf_real/lifecourse",
                 "_lifecourse.csv$", full.names = TRUE)


out <- rbindlist(lapply(fl, fread))[scenario == "sc0", ]
out[, dimd := factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived"))]
setkey(out, mc, pid, year)
out[, table(cmsmm1.5_prvl)]
View(out[cmsmm1.5_prvl > 0, .(mc, pid, year, cmsmm1.5_prvl) ])
t1 <- out[cmsmm1.5_prvl == 1 & year %in% c(20, 40), .(mm_free = weighted.mean(age, wt_esp)), keyby = .(year, sex)] # mean age of multimorbidity
t2 <- out[all_cause_mrtl > 0 & year %in% c(20, 40), .(le = weighted.mean(age, wt_esp)), keyby = .(year, sex)] # mean life expectancy
t1 <- t1[t2]

t1 <- out[cmsmm1.5_prvl == 1 & year %in% c(20, 40), .(mm_free = weighted.mean(age, wt_esp)), keyby = .(year, dimd)] # mean age of multimorbidity
t2 <- out[all_cause_mrtl > 0 & year %in% c(20, 40), .(le = weighted.mean(age, wt_esp)), keyby = .(year, dimd)] # mean life expectancy
t1 <- t1[t2]
View(t1)


tt_esp <- fread("/mnt/storage_fast/output/hf_real/summaries/incd_esp.csv"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
outstrata <- c("mc", "year", "dimd", "scenario")
e <- tt_esp[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e <- melt(e, id.vars = outstrata)
e <- e[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_esp_")))
# fwrite(e, "/mnt/storage_fast/output/hf_real/tables/age-standardised prevalence by year and dimd.csv")

e[disease == "cmsmm0_prvl", disease := "CMS score > 0"]
e[disease == "cmsmm1_prvl", disease := "CMS score > 1"]
e[disease == "cmsmm1.5_prvl", disease := "CMS score > 1.5"]
e[disease == "cmsmm2_prvl", disease := "CMS score > 2"]
e <- e[between(year, 2014, 2040) & disease %in% c("CMS score > 1", "CMS score > 1.5", "CMS score > 2")]
e[scenario == "sc0", scenario := "Base-case"]
e[scenario == "sc1", scenario := "10% BMI reduction"]

ggplot(e[scenario == "Base-case" & disease %in% c("CMS score > 1", "CMS score > 1.5", "CMS score > 2")],
       aes(x = year, y = `prvl_rate_esp_50.0%`, ymin = `prvl_rate_esp_2.5%`,
           ymax = `prvl_rate_esp_97.5%`, colour = dimd, fill = dimd)) +
  # geom_ribbon(alpha = 2/5, colour = NA) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Prevalence rate", labels = percent) +
  ggtitle(paste0("Prevalence of ", "multimorbidity", " by dimd (age-standardised)"))+
  expand_limits(y = 0) +
  facet_grid(.~disease)
ggsave("~/pCloudDrive/pCloud Sync/MM prevalence by dimd.png", scale = 1.5, width = 16/2, height = 9/2)




# for Tobys HLE
tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))][scenario == "sc0"]

outstrata <- c("mc", "year", "agegrp", "sex")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
setkey(d, year, sex, agegrp)
d[year == 2020 & disease == "cmsmm1.5_prvl"]
fwrite(d[year %in% c(2020, 2040) & disease == "cmsmm1.5_prvl", .(year, sex, agegrp, 1-`prvl_rate_50.0%`)], "~/pCloudDrive/pCloud Sync/mm.csv")
fwrite(d[year %in% c(2020, 2040) & disease == "cmsmm0_prvl", .(year, sex, agegrp, 1-`prvl_rate_50.0%`)], "~/pCloudDrive/pCloud Sync/mm0.csv")

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_scaled_up.csv"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))][scenario == "sc0"]

outstrata <- c("mc", "year", "agegrp", "sex", "dimd")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
d <- melt(d, id.vars = outstrata)
d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
setkey(d, year, dimd, sex, agegrp)
d[year == 2020 & disease == "cmsmm1.5_prvl"]
fwrite(d[year %in% c(2020, 2040) & dimd %in% c("1 most deprived", as.character(5), "10 least deprived") & disease == "cmsmm1.5_prvl", .(year, dimd, sex, agegrp, 1-`prvl_rate_50.0%`)], "~/pCloudDrive/pCloud Sync/mm_dimd.csv")
fwrite(d[year %in% c(2020, 2040) & dimd %in% c("1 most deprived", as.character(5), "10 least deprived")& disease == "cmsmm0_prvl", .(year, dimd, sex, agegrp, 1-`prvl_rate_50.0%`)], "~/pCloudDrive/pCloud Sync/mm0_dimd.csv")

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_out.csv"
)[, `:=` (year = year + 2000L,
          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))][scenario == "sc0"]
outstrata <- c("mc", "year", "agegrp", "sex")
e <- tt[, lapply(.SD, sum), .SDcols = patterns("V1"), keyby = eval(outstrata)
][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
e <- melt(e, id.vars = outstrata)
e <- e[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(e, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "prvl_rate_")))
setkey(e, year, sex, agegrp)
d[year == 2020 & disease == "cmsmm1.5_prvl"]

IMPACTncd$export_summaries()
