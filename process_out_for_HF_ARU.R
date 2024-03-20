library(data.table)
library(CKutils)
library(scales)
library(yaml)


prbl <- c(0.5, 0.025, 0.975, 0.1, 0.9) # user input
## These parameters will come from sim_design.yaml and design class
simulationParameters <- read_yaml(base::normalizePath("./inputs/sim_design.yaml", mustWork = TRUE))
sSummariesSubDirPath <- file.path(simulationParameters$output_dir, "summaries/")
sTablesSubDirPath <- file.path(simulationParameters$output_dir, "tables/")
output_dir <- simulationParameters$output_dir

all_cause_mrtl_by_dis <- function(outstrata, tt, prbl, sTablesSubDirPath, file_name,
                                 pp, pop_denom) {
  d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
          keyby = eval(outstrata)]
  d <- melt(d, id.vars = outstrata)
  if(!pop_denom) {
    cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
  }
  d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
  if(pop_denom) {
    cases <- pp[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
    d[cases, on = outstrata, value := value/popsize]
  } else {
    d[cases, on = c(outstrata, "variable"), value := value/i.value]
  }
  setkey(d, "variable")
  d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
  setnames(d, c(setdiff(outstrata, "mc"), "disease", percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
  setkeyv(d, setdiff(outstrata, "mc"))
  fwrite(d, paste0(sTablesSubDirPath, file_name))
}

dis_chrs <- function(tt, outstrata, prbl, sTablesSubDirPath, file_name) {
  d1 <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|dimd|^cases_")]
  d1 <- melt(d1, id.vars = c("mc", "year", "scenario", "sex", "dimd"))
  d1 <- unique(d1, by = c("mc", "year", "scenario", "sex", "dimd", "variable"))
  d1[, `:=` (disease = gsub("^cases_", "", variable), variable = NULL)]
  tt <- tt[, .SD, .SDcols = patterns("mc|scenario|year|sex|dimd|^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_")]
  tt <- melt(tt, id.vars = c("mc", "year", "scenario", "sex", "dimd"))
  tt[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
  tt[d1, on = c("mc", "year", "scenario", "sex", "dimd", "disease"), cases := i.value]
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
  fwrite(d, paste0(sTablesSubDirPath, file_name))
  return(d)
}


# All-cause mortality by disease not standardised----
tt <- fread(paste0(sSummariesSubDirPath,"/all_cause_mrtl_by_dis_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
                                                                                            dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
pp <- fread(paste0(sSummariesSubDirPath,"/prvl_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
                                                                           dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]

outstrata <- c("mc", "year", "scenario")
file_name <-  "/all-cause mrtl by disease-year (not standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = F)

outstrata <- c("mc", "year", "sex", "scenario")
file_name <- "/all-cause mrtl by disease-year-sex (not standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = F)

outstrata <- c("mc", "year", "agegrp", "scenario")
file_name <- "/all-cause mrtl by disease-year-agegrp (not standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = F)

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
file_name <- "/all-cause mrtl by disease-year-agegroup-sex (not standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = F)

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
file_name <- "/all-cause mrtl by disease-year-agegroup-sex-dimd (not standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = F)


# All-cause mortality by disease not standardised pop denominator----

outstrata <- c("mc", "year", "scenario")
file_name <- "/all-cause mrtl by disease-year popdenom (not standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = T)

outstrata <- c("mc", "year", "sex", "scenario")
file_name <- "/all-cause mrtl by disease-year-sex popdenom (not standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = T)

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")
file_name <- "/all-cause mrtl by disease-year-agegroup-sex popdenom (not standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = T)

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")
file_name <- "/all-cause mrtl by disease-year-agegroup-sex-dimd popdenom (not standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = T)

rm(pp)

# All-cause mortality by disease standardised----

tt <- fread(paste0(sSummariesSubDirPath,"/all_cause_mrtl_by_dis_esp.csv.gz"))[, `:=` (year = year + 2000L,
                                                                                      dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
outstrata <- c("mc", "year", "scenario")
file_name <- "/all-cause mrtl by disease-year (age-sex-dimd standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = F)

outstrata <- c("mc", "year", "sex", "scenario")
file_name <- "/all-cause mrtl by disease-year-sex (age-dimd standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = F)

outstrata <- c("mc", "year", "dimd", "scenario")
file_name <- "/all-cause mrtl by disease-year-dimd (age-sex standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = F)

outstrata <- c("mc", "year", "sex", "dimd", "scenario")
file_name <- "/all-cause mrtl by disease-year-sex-dimd (age standardised).csv"
# all_cause_mrtl_by_dis(outstrata, tt, prbl, sTablesSubDirPath, file_name,
#                                   pp, pop_denom = F)


# Disease characteristics non standardised ----

tt <- fread(paste0(sSummariesSubDirPath,"/dis_characteristics_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
                                                                                          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")),
                                                                                         mean_cms_count_cms1st_cont = as.numeric(mean_cms_count_cms1st_cont))]
outstrata <- c("mc", "year", "scenario")
file_name <- "/disease characteristics by year (not standardised).csv"
# dis_chrs(tt, outstrata, prbl, sTablesSubDirPath, file_name)

outstrata <- c("mc", "year", "sex", "scenario")
file_name <- "/disease characteristics by year-sex (not standardised).csv"
# dis_chrs(tt, outstrata, prbl, sTablesSubDirPath, file_name)

outstrata <- c("mc", "year", "dimd", "scenario")
file_name <- "/disease characteristics by year-dimd (not standardised).csv"
# dis_chrs(tt, outstrata, prbl, sTablesSubDirPath, file_name)

outstrata <- c("mc", "year", "sex", "dimd", "scenario")
file_name <- "/disease characteristics by year-sex-dimd (not standardised).csv"
# dis_chrs(tt, outstrata, prbl, sTablesSubDirPath, file_name)

