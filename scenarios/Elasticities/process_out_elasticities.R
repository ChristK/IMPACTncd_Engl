library(data.table)
library(yaml)
library(CKutils)
library(ggplot2)
library(ggthemes)
library(scales)
library(stringr)
setDTthreads(5)
prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)
# theme_set(new = theme_economist())
# theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))

# Setting the plot settings
theme_set(new = theme_economist_white(gray_bg = FALSE))
theme_update(axis.text = element_text(size = 10),
             axis.title = element_text(size = 10 ),
             legend.position = "bottom",
             legend.text=element_text(size=10),

             #plot.title = element_text(hjust = 0.5),
)


simulationParameters <- read_yaml(base::normalizePath("./inputs/sim_design_elasticities.yaml", mustWork=TRUE))
output_dir <- simulationParameters$output_dir
# output_dir <- "/mnt/storage_fast/output/hf_real_elasticities_tmp"
sSummariesSubDirPath <- file.path(output_dir,"summaries/")
sTablesSubDirPath <- file.path(output_dir,"tables/")
sPlotsSubDirPath <- file.path(output_dir,"tables/tmpplots/")


changeyr <- 23L
yrs <- c(changeyr + 5L, changeyr + 10L, changeyr + 15L, changeyr + 20L)
yrs2 <- c(yrs, changeyr + 25L, changeyr + 30L)
type <- "_scaled_up.csv.gz"
strata <- c("mc","scenario" ,"year","agegrp","sex", "dimd")

#prvl unstandardised
fpth <- file.path(output_dir, "summaries",  paste0("prvl", type))
prvl <- fread(fpth, header = T, select = c(strata, "popsize", "cmsmm1.5_prvl"))

#Extracting scenarios by type
scenarios <- prvl[, unique(scenario)]
sc_1pc <- str_subset(scenarios, "1")
sc_10pc <- str_subset(scenarios, "2")
sc_opt <- str_subset(scenarios, "parf")
#sc_opt <- c(sc_opt[1:5], "NA", sc_opt[6])

sc_nm <- c("Alcohol", "BMI", "Fruit & veg", "Physical activity", "SBP", "Smoking", "Total cholesterol")

lu <- data.table(sc_1pc = sc_1pc, sc_10pc = sc_10pc, sc_opt =sc_opt, sc_nm = sc_nm)
lu_l <- melt(lu, id.vars = "sc_nm")

# By year
outstrata <- c("mc", "scenario", "year" )
d <- prvl[, .(sum(cmsmm1.5_prvl)/sum(popsize) * 100), keyby = eval(outstrata)]
setkey(d, mc, scenario)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "prvl_chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in prevalence by year (unstandardised).csv"))

d[year %in% yrs, .(`prvl_chnge50.0%`, `prvl_chnge2.5%`, `prvl_chnge97.5%`), keyby = .(scenario, year)]
d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & year >= 23],
       aes(x = year + 2000,
           y = `prvl_chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in prevalence compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Prvl_change_10pc.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)


ggplot(d[scenario %in% sc_1pc & year >= 23],
       aes(x = year + 2000,
           y = `prvl_chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in prevalence compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("1% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Prvl_change_1pc.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)




#By year and DIMD
outstrata <- c("mc", "scenario", "year" , "dimd")
d <- prvl[, .(sum(cmsmm1.5_prvl)/sum(popsize) * 100), keyby = eval(outstrata)]
setkey(d, mc, scenario)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "dimd"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "prvl_chnge")))
d[, dimd := factor(dimd,
                   levels = c("1 most deprived", "2", "3", "4", "5",
                              "6", "7", "8", "9", "10 least deprived"))]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in prevalence by year, dimd (unstandardised).csv"))

d[year %in% yrs, .(`prvl_chnge50.0%`, `prvl_chnge2.5%`, `prvl_chnge97.5%`), keyby = .(scenario, year, dimd)]

d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & dimd %in% c("1 most deprived", "10 least deprived") & year >= 23],
       aes(x = year + 2000,
           y = `prvl_chnge50.0%`,
           col = sc_nm)) +
  facet_grid( ~dimd) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in prevalence compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Prvl_change_10pc_dimd.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)



#incd unstandardised
fpth <- file.path(output_dir, "summaries",  paste0("incd", type))
incd <- fread(fpth, header = T, select = c(strata, "popsize", "cmsmm1.5_incd"))
incd[prvl, on = strata, cmsmm1.5_prvl := i.cmsmm1.5_prvl ]
incd[, popsize := popsize - cmsmm1.5_prvl + cmsmm1.5_incd]

# By year
outstrata <- c("mc", "scenario", "year" )
d <- incd[, .(sum(cmsmm1.5_incd)/sum(popsize) * 10000), keyby = eval(outstrata)]
setkey(d, mc, scenario)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "incd_chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in incidence by year (unstandardised).csv"))
d[year %in% yrs, .(`incd_chnge50.0%`, `incd_chnge2.5%`, `incd_chnge97.5%`), keyby = .(scenario, year)]

d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & year >= 23],
       aes(x = year + 2000,
           y = `incd_chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in incidence compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")

ggsave(
  filename = paste0(sPlotsSubDirPath, "Incd_change_10pc.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)

ggplot(d[scenario %in% sc_1pc & year >= 23],
       aes(x = year + 2000,
           y = `incd_chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in incidence compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("1% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Incd_change_1pc.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)






#By year and DIMD
outstrata <- c("mc", "scenario", "year" , "dimd")
d <- incd[, .(sum(cmsmm1.5_incd)/sum(popsize) * 10000), keyby = eval(outstrata)]
setkey(d, mc, scenario)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "dimd"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "incd_chnge")))
d[, dimd := factor(dimd,
                   levels = c("1 most deprived", "2", "3", "4", "5",
                              "6", "7", "8", "9", "10 least deprived"))]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in incidence by year, dimd (unstandardised).csv"))
d[year %in% yrs, .(`incd_chnge50.0%`, `incd_chnge2.5%`, `incd_chnge97.5%`), keyby = .(scenario, year, dimd)]

d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & dimd %in% c("1 most deprived", "10 least deprived") & year >= 23],
       aes(x = year + 2000,
           y = `incd_chnge50.0%`,
           col = sc_nm)) +
  facet_grid( ~dimd) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in incidence compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Incd_change_10pc_dimd.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)





#By year & age
outstrata <- c("mc", "year", "agegrp", "scenario")
d <- incd[, .(sum(cmsmm1.5_incd)/sum(popsize) * 10000), keyby = eval(outstrata)]
setkey(d, mc, scenario, agegrp)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "agegrp"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "incd_chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in incidence by year, agegroup (unstandardised).csv"))
d[year %in% yrs & scenario == "bmi2", .(`incd_chnge50.0%`, `incd_chnge2.5%`, `incd_chnge97.5%`), keyby = .(scenario, year, agegrp)]



#Allcause mortality
fpth <- file.path(output_dir, "summaries",  paste0("mrtl", type))
tt <- fread(fpth, header = T)
# By year
outstrata <- c("mc", "scenario", "year" )
d <- tt[, .(sum(all_cause_mrtl)/sum(popsize) * 10000), keyby = eval(outstrata)]
setkey(d, mc, scenario)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in all-cause mortality by year (unstandardised).csv"))
d[year %in% yrs, .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year)]

d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]


ggplot(d[scenario %in% sc_10pc & year >= 23],
       aes(x = year + 2000,
           y = `chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in all-cause mortality compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")

ggsave(
  filename = paste0(sPlotsSubDirPath, "Mtlt_change_10pc.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)

ggplot(d[scenario %in% sc_1pc & year >= 23],
       aes(x = year + 2000,
           y = `chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in all-cause mortality compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("1% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Mtlt_change_1pc.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)






#By year and DIMD
outstrata <- c("mc", "scenario", "year" , "dimd")
d <- tt[, .(sum(all_cause_mrtl)/sum(popsize) * 10000), keyby = eval(outstrata)]
setkey(d, mc, scenario)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "dimd"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
d[, dimd := factor(dimd,
                   levels = c("1 most deprived", "2", "3", "4", "5",
                              "6", "7", "8", "9", "10 least deprived"))]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in all-cause mortality by year, dimd (unstandardised).csv"))

d[year %in% yrs, .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year, dimd)]
d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & dimd %in% c("1 most deprived", "10 least deprived") & year >= 23],
       aes(x = year + 2000,
           y = `chnge50.0%`,
           col = sc_nm)) +
  facet_grid( ~dimd) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in all-cause mortality compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Mtlt_change_10pc_dimd.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)





#By year & age
outstrata <- c("mc", "year", "agegrp", "scenario")
d <- tt[, .(sum(all_cause_mrtl)/sum(popsize) * 10000), keyby = eval(outstrata)]
setkey(d, mc, scenario, agegrp)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "agegrp"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "mtlt_chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in all-cause mortality by year (unstandardised).csv"))
d[year %in% yrs & scenario == "bmi2", .(`mtlt_chnge50.0%`, `mtlt_chnge2.5%`, `mtlt_chnge97.5%`), keyby = .(scenario, year, agegrp)]




#Allcause mortality with CMS >1.5
fpth <- file.path(output_dir, "summaries",  paste0("all_cause_mrtl_by_dis", type))
tt <- fread(fpth, header = T)
tt[, dimd := factor(dimd,
                   levels = c("1 most deprived", "2", "3", "4", "5",
                              "6", "7", "8", "9", "10 least deprived"))]

# By year
outstrata <- c("mc", "year", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")

outstrata <- c("mc", "scenario", "year" ,"variable")
setkey(d, mc, scenario, variable)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "variable"), bl := i.value][, change := ((value / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in disease all-cause mortality by year (unstandardised).csv"))
d[year %in% yrs & variable == "cmsmm1.5", .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year)]


#By year and DIMD
outstrata <- c("mc", "year", "dimd","scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")

outstrata <- c("mc", "scenario", "year","variable", "dimd")
setkey(d, mc, scenario, variable)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "variable", "dimd"), bl := i.value][, change := ((value / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in disease all-cause mortality by year, dimd (unstandardised).csv"))
d[year %in% yrs & variable == "cmsmm1.5", .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year, variable, dimd)]


#By year & age
outstrata <- c("mc", "year", "agegrp", "scenario")
d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
        keyby = eval(outstrata)]
d <- melt(d, id.vars = outstrata)
cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
d[cases, on = c(outstrata, "variable"), value := value/i.value]
setkey(d, "variable")

outstrata <- c("mc", "scenario", "year","variable", "agegrp")
setkey(d, mc, scenario, variable)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "variable", "agegrp"), bl := i.value][, change := ((value / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in disease all-cause mortality by year, agegroup (unstandardised).csv"))
d[year %in% yrs & variable == "cmsmm1.5" & scenario == "bmi2", .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year, variable, agegrp)]








# Disease characteristics -
fpth <- file.path(output_dir, "summaries",  paste0("dis_characteristics", type))
tt <- fread(fpth, header = T)[
  ,mean_cms_count_cms1st_cont := as.numeric(mean_cms_count_cms1st_cont)]
d1 <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|dimd|^cases_")]
d1 <- melt(d1, id.vars = c("mc", "year", "scenario", "sex", "dimd"))
d1 <- unique(d1, by = c("mc", "year", "scenario", "sex", "dimd", "variable"))
d1[, `:=` (disease = gsub("^cases_", "", variable), variable = NULL)]
tt <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|dimd|^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_")]
tt <- melt(tt, id.vars = c("mc", "year", "scenario", "sex", "dimd"))
tt[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
tt[d1, on = c("mc", "year", "scenario", "sex", "dimd", "disease"), cases := i.value]
rm(d1) # NOTE mean_age_incd contains NAs

outstrata <- c("mc", "year", "scenario", "variable")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata)] # na.rm = TRUE for mean_age_incd
setkey(d, mc, scenario, "variable")
d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "variable"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in disease characteristics by year (unstandardised).csv"))
d[year %in% yrs2 & variable == "mean_duration_cmsmm1.5",
  .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year)]

d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & year >= 23 & variable == "mean_duration_cmsmm1.5"],
       aes(x = year + 2000,
           y = `chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in mean duration of CMS>1.5 compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "CMS1.5dur_change_10pc_dimd.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)


d[year %in% yrs2 & variable == "mean_age_1st_onset_cmsmm1.5",
  .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year)]
ggplot(d[scenario %in% sc_10pc & year >= 23 & variable == "mean_age_1st_onset_cmsmm1.5"],
       aes(x = year + 2000,
           y = `chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in age of CMS>1.5 onset compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "CMS1.5onsetage_change_10pc_dimd.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)




# By year and dimd
outstrata <- c("mc", "year", "scenario","dimd", "variable")
d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata)] # na.rm = TRUE for mean_age_incd
setkey(d, mc, scenario, dimd , variable)
d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "dimd","variable"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
d[, dimd := factor(dimd,
                   levels = c("1 most deprived", "2", "3", "4", "5",
                              "6", "7", "8", "9", "10 least deprived"))]
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in disease characteristics by year, dimd (unstandardised).csv"))
d[year %in% yrs & variable == "mean_duration_cmsmm1.5",
  .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year, dimd)]
d[year %in% yrs & variable == "mean_age_1st_onset_cmsmm1.5",
  .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year, dimd)]



# Mean cms score
# Mean CMS not standardised ----
tt <- fread(paste0(sSummariesSubDirPath,"cms_score_scaled_up.csv.gz"))[, `:=` (
         dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
#by year
outstrata <- c("mc", "year", "scenario")
d <- tt[,  weighted.mean(cms_score, popsize),
         keyby = eval(outstrata)]
setkey(d, mc, scenario, year)
d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in mean CMS score by year (not standardised).csv"))
d[year %in% yrs, .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year)]
d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & year >= 23 ],
       aes(x = year + 2000,
           y = `chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in mean CMS score compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Meancms_change_10pc_dimd.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)


#by year & dimd
outstrata <- c("mc", "year", "scenario", "dimd")
d <- tt[,  weighted.mean(cms_score, popsize),
        keyby = eval(outstrata)]
setkey(d, mc, scenario, year)
d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "dimd"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in mean CMS score by year, dimd (not standardised).csv"))
d[year %in% yrs, .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year, dimd)]




# Mean cms count
tt <- fread(paste0(sSummariesSubDirPath,"cms_count_scaled_up.csv.gz"))[, `:=` (
  dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
#by year
outstrata <- c("mc", "year", "scenario")
d <- tt[,  weighted.mean(cms_count, popsize),
        keyby = eval(outstrata)]
setkey(d, mc, scenario, year)
d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in mean CMS count by year (not standardised).csv"))
d[year %in% yrs, .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year)]

d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & year >= 23 ],
       aes(x = year + 2000,
           y = `chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in mean CMS count compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Meancms_count_10pc_dimd.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)


#by year & dimd
outstrata <- c("mc", "year", "scenario", "dimd")
d <- tt[,  weighted.mean(cms_count, popsize),
        keyby = eval(outstrata)]
setkey(d, mc, scenario, year)
d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year", "dimd"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in mean CMS count by year, dimd (not standardised).csv"))
d[year %in% yrs, .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year, dimd)]



#No disease unstandardised
fpth <- file.path(output_dir, "summaries",  paste0("prvl", type))
prvl <- fread(fpth, header = T, select = c(strata, "popsize", "cmsmm0_prvl", "cmsmm1.5_prvl"))
#cmsmm0_prvl is the prev of people with CMS > 0
# By year
outstrata <- c("mc", "scenario", "year" )
d <- prvl[, .((1-sum(cmsmm0_prvl)/sum(popsize)) * 100), keyby = eval(outstrata)]
setkey(d, mc, scenario)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "prvl_chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in healthy prevalence by year (unstandardised).csv"))

d[year %in% yrs, .(`prvl_chnge50.0%`, `prvl_chnge2.5%`, `prvl_chnge97.5%`), keyby = .(scenario, year)]

d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & year >= 23 ],
       aes(x = year + 2000,
           y = `prvl_chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in healthy prevalence compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Healthyprvlcnge_10pc_dimd.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)


#Minor illness
# By year
outstrata <- c("mc", "scenario", "year" )
d <- prvl[, .((sum(cmsmm0_prvl - cmsmm1.5_prvl)/sum(popsize)) * 100), keyby = eval(outstrata)]
setkey(d, mc, scenario)

d_sc0 <- d[scenario == "sc0"]
d <- d[scenario != "sc0"]
d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]

d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "prvl_chnge")))
setkeyv(d, setdiff(outstrata, "mc"))
fwrite(d, paste0(sTablesSubDirPath,"scenario change in 0<CMS<1.5 prevalence by year (unstandardised).csv"))

d[year %in% yrs, .(`prvl_chnge50.0%`, `prvl_chnge2.5%`, `prvl_chnge97.5%`), keyby = .(scenario, year)]

d[lu_l, on = c("scenario" = "value"), sc_nm := i.sc_nm]

ggplot(d[scenario %in% sc_10pc & year >= 23 ],
       aes(x = year + 2000,
           y = `prvl_chnge50.0%`,
           col = sc_nm)) +
  geom_line() +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "Relative change in minor illness prevalence compared to baseline") +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  #scale_color_brewer(palette = "Set3") +
  theme(#legend.direction = "vertical",
    legend.title=element_blank(),
    #legend.box = "horizontal", legend.text = element_text(size = 12),
    legend.position = "bottom") +
  ggtitle("10% improvement in risk factors")
ggsave(
  filename = paste0(sPlotsSubDirPath, "Minorill_prvlcnge_10pc_dimd.png"),
  width = 12, height = 8, scale = 0.75, dpi = 600)



#
#
#
# #PARFS - for now these are the old ones
#
#
# # simulationParameters <- read_yaml(base::normalizePath("./auxil/sim_design_parf.yaml", mustWork=TRUE))
# # sSummariesSubDirPath <- file.path(simulationParameters$output_dir,"summaries/")
# # sTablesSubDirPath <- file.path(simulationParameters$output_dir,"tables/")
# # output_dir <- simulationParameters$output_dir
# output_dir <- "/mnt/storage_fast/output/hf_real_parf"
#
# changeyr <- 23L
# yrs <- c(changeyr + 5L, changeyr + 10L, changeyr + 15L, changeyr + 20L)
# type <- "_scaled_up.csv.gz"
# strata <- c("mc","scenario" ,"year","agegrp","sex", "dimd")
#
# #prvl unstandardised
# fpth <- file.path(output_dir, "summaries",  paste0("prvl", type))
# prvl <- fread(fpth, header = T, select = c(strata, "popsize", "cmsmm1.5_prvl"))[scenario %in% c("bmi", "sc0")]
#
# # By year
# outstrata <- c("mc", "scenario", "year" )
# d <- prvl[, .(sum(cmsmm1.5_prvl)/sum(popsize) * 100), keyby = eval(outstrata)]
# setkey(d, mc, scenario)
#
# d_sc0 <- d[scenario == "sc0"]
# d <- d[scenario != "sc0"]
# d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]
#
# d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "prvl_chnge")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"parf change in prevalence by year (unstandardised).csv"))
#
# d[year %in% yrs, .(`prvl_chnge50.0%`, `prvl_chnge2.5%`, `prvl_chnge97.5%`), keyby = .(scenario, year)]
#
#
# #By year and DIMD
# outstrata <- c("mc", "scenario", "year" , "dimd")
# d <- prvl[, .(sum(cmsmm1.5_prvl)/sum(popsize) * 100), keyby = eval(outstrata)]
# setkey(d, mc, scenario)
#
# d_sc0 <- d[scenario == "sc0"]
# d <- d[scenario != "sc0"]
# d[d_sc0, on = c("mc", "year", "dimd"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]
#
# d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "prvl_chnge")))
# d[, dimd := factor(dimd,
#                    levels = c("1 most deprived", "2", "3", "4", "5",
#                               "6", "7", "8", "9", "10 least deprived"))]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"parf change in prevalence by year, dimd (unstandardised).csv"))
#
# d[year %in% yrs, .(`prvl_chnge50.0%`, `prvl_chnge2.5%`, `prvl_chnge97.5%`), keyby = .(scenario, year, dimd)]
#
#
#
# #incd unstandardised
# fpth <- file.path(output_dir, "summaries",  paste0("incd", type))
# incd <- fread(fpth, header = T, select = c(strata, "popsize", "cmsmm1.5_prvl"))[scenario %in% c("bmi", "sc0")]
# setnames(incd, "cmsmm1.5_prvl", "cmsmm1.5_incd")
# incd[prvl, on = strata, cmsmm1.5_prvl := i.cmsmm1.5_prvl ]
# incd[, popsize := popsize - cmsmm1.5_prvl + cmsmm1.5_incd]
#
# # By year
# outstrata <- c("mc", "scenario", "year" )
# d <- incd[, .(sum(cmsmm1.5_incd)/sum(popsize) * 10000), keyby = eval(outstrata)]
# setkey(d, mc, scenario)
#
# d_sc0 <- d[scenario == "sc0"]
# d <- d[scenario != "sc0"]
# d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]
#
# d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "incd_chnge")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"parf change in incidence by year (unstandardised).csv"))
# d[year %in% yrs, .(`incd_chnge50.0%`, `incd_chnge2.5%`, `incd_chnge97.5%`), keyby = .(scenario, year)]
#
#
# #By year and DIMD
# outstrata <- c("mc", "scenario", "year" , "dimd")
# d <- incd[, .(sum(cmsmm1.5_incd)/sum(popsize) * 10000), keyby = eval(outstrata)]
# setkey(d, mc, scenario)
#
# d_sc0 <- d[scenario == "sc0"]
# d <- d[scenario != "sc0"]
# d[d_sc0, on = c("mc", "year", "dimd"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]
#
# d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "incd_chnge")))
# d[, dimd := factor(dimd,
#                    levels = c("1 most deprived", "2", "3", "4", "5",
#                               "6", "7", "8", "9", "10 least deprived"))]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"parf change in incidence by year, dimd (unstandardised).csv"))
# d[year %in% yrs, .(`incd_chnge50.0%`, `incd_chnge2.5%`, `incd_chnge97.5%`), keyby = .(scenario, year, dimd)]
#
#
#
# #Allcause mortality
# fpth <- file.path(output_dir, "summaries",  paste0("mrtl", type))
# tt <- fread(fpth, header = T)[scenario %in% c("bmi", "sc0")]
# # By year
# outstrata <- c("mc", "scenario", "year" )
# d <- tt[, .(sum(all_cause_mrtl)/sum(popsize) * 10000), keyby = eval(outstrata)]
# setkey(d, mc, scenario)
#
# d_sc0 <- d[scenario == "sc0"]
# d <- d[scenario != "sc0"]
# d[d_sc0, on = c("mc", "year"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]
#
# d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"parf change in all-cause mortality by year (unstandardised).csv"))
# d[year %in% yrs, .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year)]
#
#
# #By year and DIMD
# outstrata <- c("mc", "scenario", "year" , "dimd")
# d <- tt[, .(sum(all_cause_mrtl)/sum(popsize) * 10000), keyby = eval(outstrata)]
# setkey(d, mc, scenario)
#
# d_sc0 <- d[scenario == "sc0"]
# d <- d[scenario != "sc0"]
# d[d_sc0, on = c("mc", "year", "dimd"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]
#
# d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
# d[, dimd := factor(dimd,
#                    levels = c("1 most deprived", "2", "3", "4", "5",
#                               "6", "7", "8", "9", "10 least deprived"))]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"parf change in all-cause mortality by year, dimd (unstandardised).csv"))
#
# d[year %in% yrs, .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year, dimd)]
#
#
#

#PARF disease characteristics
# fpth <- file.path(output_dir, "summaries",  paste0("dis_characteristics", type))
# tt <- fread(fpth, header = T)[
#   ,mean_cms_count_cms1st_cont := as.numeric(mean_cms_count_cms1st_cont)]
# d1 <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|dimd|^cases_")]
# d1 <- melt(d1, id.vars = c("mc", "year", "scenario", "sex", "dimd"))
# d1 <- unique(d1, by = c("mc", "year", "scenario", "sex", "dimd", "variable"))
# d1[, `:=` (disease = gsub("^cases_", "", variable), variable = NULL)]
# tt <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|dimd|^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_")]
# tt <- melt(tt, id.vars = c("mc", "year", "scenario", "sex", "dimd"))
# tt[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
# tt[d1, on = c("mc", "year", "scenario", "sex", "dimd", "disease"), cases := i.value]
# rm(d1) # NOTE mean_age_incd contains NAs
#
# outstrata <- c("mc", "year", "scenario", "variable")
# d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata)] # na.rm = TRUE for mean_age_incd
# setkey(d, mc, scenario, "variable")
# d_sc0 <- d[scenario == "sc0"]
# d <- d[scenario != "sc0"]
# d[d_sc0, on = c("mc", "year", "variable"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]
#
# d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"scenario change in disease characteristics by year (unstandardised).csv"))
# d[year %in% yrs2 & variable == "mean_duration_cmsmm1.5",
#   .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year)]
#
# # By year and dimd
# outstrata <- c("mc", "year", "scenario","dimd", "variable")
# d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata)] # na.rm = TRUE for mean_age_incd
# setkey(d, mc, scenario, dimd , variable)
# d_sc0 <- d[scenario == "sc0"]
# d <- d[scenario != "sc0"]
# d[d_sc0, on = c("mc", "year", "dimd","variable"), bl := i.V1][, change := ((V1 / bl) - 1) * 100]
#
# d <- d[, as.list(fquantile(change, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
# setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "chnge")))
# d[, dimd := factor(dimd,
#                    levels = c("1 most deprived", "2", "3", "4", "5",
#                               "6", "7", "8", "9", "10 least deprived"))]
# setkeyv(d, setdiff(outstrata, "mc"))
# fwrite(d, paste0(sTablesSubDirPath,"scenario change in disease characteristics by year, dimd (unstandardised).csv"))
# d[year %in% yrs & variable == "mean_duration_cmsmm1.5",
#   .(`chnge50.0%`, `chnge2.5%`, `chnge97.5%`), keyby = .(scenario, year, dimd)]
#
#
#
#









# Some plots - need to make a function to do by disease
library(ggplot2)
tt <- fread(paste0(sTablesSubDirPath,"disease characteristics by year (not standardised).csv"))

tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to remove

tt[, scenario := factor(scenario,
                        levels = c("sc0", "bmi1", "bmi2"),
                        labels = c("Baseline", "1% BMI reduction", "10% BMI reduction"))]
dt <- tt[type == "mean_age_1st_onset" & disease %in% c("cmsmm1.5", "t2dm") & year != 23]
dt[, disease := factor(disease,
                        levels = c("t2dm", "cmsmm1.5"),
                        labels = c("Type 2 Diabetes", "CMS>1.5"))]

#Age of onset
ggplot(dt,
       aes(x = year,
  y = `value_50.0%`,
  col = scenario,
  linetype = disease )) +
  geom_smooth(linewidth = 0.75, se = F) +
  scale_x_continuous(name = "Year since change") +
  scale_y_continuous(name = "Mean age of onset", limits = c(60,71) ) +
  # expand_limits(y = 60) +
  # theme(legend.position="bottom", legend.text=element_text(size=12)) +
  scale_color_brewer(palette = "Set1") +
  theme(#legend.direction = "vertical",
        legend.title=element_blank(),
        #legend.box = "horizontal", legend.text = element_text(size = 12),
        legend.position = "right")
  ggsave(
    filename = paste0(sTablesSubDirPath, "Mean_age_onset.png"),
    width = 8, height = 8, scale = 0.75, dpi = 600)


  #Mean CMS score
  tt <- fread(paste0(sTablesSubDirPath,"mean CMS score by year (not standardised).csv"))
  tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to remove
  tt[, scenario := factor(scenario,
                          levels = c("sc0", "bmi1", "bmi2"),
                          labels = c("Baseline", "1% BMI reduction", "10% BMI reduction"))]

  ggplot(tt,
         aes(x = year,
             y = `mean_cms_score_50.0%`,
             col = scenario )) +
    geom_smooth(linewidth = 0.5, se = F) +
    scale_x_continuous(name = "Year since change") +
    scale_y_continuous(name = "Mean CMS score" ) +
    expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "Mean_CMS_score.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)


  #Mean CMS score
  tt <- fread(paste0(sTablesSubDirPath,"mean CMS score by year-age (not standardised).csv"))
  tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to remove
  tt[, scenario := factor(scenario,
                          levels = c("sc0", "bmi1", "bmi2"),
                          labels = c("Baseline", "1% BMI reduction", "10% BMI reduction"))]
  dt <- tt[age %in% c(40,60,80)]
  ggplot(dt,
         aes(x = year,
             y = `mean_cms_score_50.0%`,
             col = scenario )) +
    geom_smooth(linewidth = 0.5, se = F) +
    scale_x_continuous(name = "Year since change") +
    scale_y_continuous(name = "Mean CMS score" ) +
    facet_grid(cols = vars(age)) +
    expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "Mean_CMS_score_age.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)



  #Mean CMS count
  tt <- fread(paste0(sTablesSubDirPath,"mean CMS count by year (not standardised).csv"))
  tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to remove
  tt[, scenario := factor(scenario,
                          levels = c("sc0", "bmi1", "bmi2"),
                          labels = c("Baseline", "1% BMI reduction", "10% BMI reduction"))]

  ggplot(tt,
         aes(x = year,
             y = `mean_cms_count_50.0%`,
             col = scenario )) +
    geom_smooth(linewidth = 0.5, se = F) +
    scale_x_continuous(name = "Year since change") +
    scale_y_continuous(name = "Mean CMS count" ) +
    expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "Mean_CMS_count.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)
  allage <- copy(tt)
  allage[, agegrp := "all ages"]

  #Mean CMS score
  tt <- fread(paste0(sTablesSubDirPath,"mean CMS count by year-agegrp (not standardised).csv"))
  tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to remove
  tt[, scenario := factor(scenario,
                          levels = c("sc0", "bmi1", "bmi2"),
                          labels = c("Baseline", "1% BMI reduction", "10% BMI reduction"))]
  dt <- tt[agegrp %in% c("40-44","60-64","80-84")]
  dt <- rbind(dt, allage)
  dt[ ,agegrp := factor(agegrp,
                       levels = c("all ages", "40-44", "60-64", "80-84"))]

  ggplot(dt,
         aes(x = year,
             y = `mean_cms_count_50.0%`,
             col = scenario )) +
    geom_smooth(linewidth = 0.5, se = F) +
    scale_x_continuous(name = "Year since change") +
    scale_y_continuous(name = "Mean CMS count" ) +
    facet_grid(cols = vars(agegrp)) +
    expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "Mean_CMS_count_agegrp.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)



  #Incidence by agegroup
  tt <- fread(paste0(sTablesSubDirPath,"incidence by year-agegrp (not standardised).csv"))[
    disease == "cmsmm1.5_incd"]
  tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to remove
  tt[, scenario := factor(scenario,
                          levels = c("sc0", "bmi1", "bmi2"),
                          labels = c("Baseline", "1% BMI reduction", "10% BMI reduction"))]

  ggplot(tt[year != 2023],
         aes(x = year,
             y = `incd_rate_50.0%` * 10000,
             col = scenario )) +
    geom_smooth(linewidth = 0.5, se = F) +
   # geom_line(linewidth = 0.5) +
    facet_wrap("agegrp", scales = "free") +
    theme(strip.text = element_text(size = 9)) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Incidence of CMS>1.5 per 10,000\n unstandardised" ) +
    #expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "Incd_age_unstand.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)


  #Mortlaity by agegroup
  tt <- fread(paste0(sTablesSubDirPath,"mortality by year-agegrp (not standardised).csv"))
  tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to remove
  tt[, scenario := factor(scenario,
                          levels = c("sc0", "bmi1", "bmi2"),
                          labels = c("Baseline", "1% BMI reduction", "10% BMI reduction"))]

  ggplot(tt[year != 2023],
         aes(x = year,
             y = `mrtl_rate_50.0%` * 10000,
             col = scenario )) +
    geom_smooth(linewidth = 0.5, se = F) +
    # geom_line(linewidth = 0.5) +
    facet_wrap("agegrp", scales = "free") +
    theme(strip.text = element_text(size = 9)) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Mortality per 10,000\n unstandardised" ) +
    #expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "Mtlt_age_unstand.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)


  #Mortlaity by agegroup w CMS>1.5
  tt <- fread(paste0(sTablesSubDirPath,"all-cause mrtl by disease-year-agegrp (not standardised).csv"))[
    disease == "cmsmm1.5"]
  tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to remove
  tt[, scenario := factor(scenario,
                          levels = c("sc0", "bmi1", "bmi2"),
                          labels = c("Baseline", "1% BMI reduction", "10% BMI reduction"))]

  ggplot(tt[year != 2023],
         aes(x = year,
             y = `all_cause_mrtl_by_disease_rate_50.0%` * 10000,
             col = scenario )) +
    geom_smooth(linewidth = 0.5, se = F) +
    # geom_line(linewidth = 0.5) +
    facet_wrap("agegrp", scales = "free") +
    theme(strip.text = element_text(size = 9)) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Mortality per 10,000\n unstandardised" ) +
    #expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "CMS1.5_mtlt_age_unstand.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)


  #Change in prevalence of CMS==0, 0<CMS<=1.5, CMS > 1.5
  tt <- rbind(fread(paste0(sTablesSubDirPath,"scenario change in 0<CMS<1.5 prevalence by year (unstandardised).csv"))[, cms := "cms>0 & <=1.5"],
              fread(paste0(sTablesSubDirPath,"scenario change in prevalence by year (unstandardised).csv"))[, cms := "cms>1.5"],
              fread(paste0(sTablesSubDirPath,"scenario change in healthy prevalence by year (unstandardised).csv"))[, cms := "cms=0"])
  tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to remove
  ggplot(tt[year != 2023 & scenario == "bmi2"],
         aes(x = year + 2000,
             y = `prvl_chnge50.0%`,
             col = cms )) +
    geom_smooth(linewidth = 0.5, se = F) +
    # geom_line(linewidth = 0.5) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "% change compared to baseline" , limits = c(-5, 5)) +
    #expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "prvl_change.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)



  # Exposures
  tt <- fread(paste0(output_dir,"/xps/xps20.csv.gz"))[,
              .(mc, scenario, year, sex, agegrp20, qimd, bmi_curr_xps)]

  tt <- tt[scenario == "sc0" | scenario %like% "bmi"] #Need to change this

  tt[, scenario := factor(scenario,
                          levels = c("sc0", "bmi1", "bmi2", "bmi_parf"),
                          labels = c("Baseline", "1% BMI reduction", "10% BMI reduction", "Optimal BMI"))]

  outstrata <- c("mc", "year", "scenario")
  d <- tt[sex == "All" & agegrp20 == "All" & qimd == "All"] # This should depend on outstrata

  d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
  setkey(d, "mc", "year", "scenario")
  d <- d[, as.list(fquantile(bmi_curr_xps, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
  setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "xps_mean_")))
  setkeyv(d, setdiff(outstrata, "mc"))

  ggplot(d[year != 2023],
         aes(x = year,
             y = `xps_mean_50.0%`,
             col = scenario )) +
    geom_smooth(linewidth = 0.5, se = F) +
    # geom_line(linewidth = 0.5) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Mean BMI", limit = c(10, 30 )) +
    #expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "mean_bmi_xps.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)

  outstrata <- c("mc", "year", "scenario", "qimd")
  d <- tt[sex == "All" & agegrp20 == "All" & qimd != "All"] # This should depend on outstrata

  d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
  setkey(d, "mc", "year", "scenario")
  d <- d[, as.list(fquantile(bmi_curr_xps, prbl)), keyby = eval(setdiff(outstrata, "mc"))]
  setnames(d, c(setdiff(outstrata, "mc"), percent(prbl, prefix = "xps_mean_")))
  setkeyv(d, setdiff(outstrata, "mc"))

  ggplot(d[year != 2023],
         aes(x = year,
             y = `xps_mean_50.0%`,
             col = scenario )) +
    geom_smooth(linewidth = 0.5, se = F) +
    facet_grid(cols = vars(qimd)) +
    # geom_line(linewidth = 0.5) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Mean BMI", limit = c(10, 30 )) +
    #expand_limits(y = 0) +
    # theme(legend.position="bottom", legend.text=element_text(size=12)) +
    scale_color_brewer(palette = "Set1") +
    theme(
      legend.title=element_blank())
  ggsave(
    filename = paste0(sTablesSubDirPath, "mean_bmi_xps_qimd.png"),
    width = 12, height = 8, scale = 0.75, dpi = 600)



#Standardised prevalence - makes v little difference
  # tt <- fread(paste0(sSummariesSubDirPath,"prvl_esp.csv.gz"))[, `:=` (year = year + 2000L,
  #           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
  #
  # outstrata <- c("mc", "year", "scenario")
  # d <- tt[, lapply(.SD, sum), .SDcols = patterns("_prvl$|^popsize$"), keyby = eval(outstrata)
  # ][, lapply(.SD, function(x) x/popsize), keyby = outstrata]
  # d <- melt(d, id.vars = outstrata)
  # setkey(d, "variable")
  # d_sc0 <- d[scenario == "sc0"]
  # d <- d[scenario != "sc0"]
  # d[d_sc0, on = c("mc", "year", "variable"), bl := i.value][, change := ((value / bl) - 1) * 100]
  # d <- d[, fquantile_byid(change, prbl, id = as.character(variable)), keyby = eval(setdiff(c(outstrata), "mc"))]
  # setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "change")))
  # d <- d[disease != "popsize"]
  # setkeyv(d, setdiff(outstrata, "mc"))
  # fwrite(d, paste0(sTablesSubDirPath,"prevalence by year (age-sex-dimd standardised).csv"))
