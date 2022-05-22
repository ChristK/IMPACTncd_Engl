# CPRD Validation charts

library(data.table) # for fast data manipulation
library(fst)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(stringr)

source(paste0("/mnt/", Sys.info()[["user"]], "/UoL/CPRD2021/epi_models/scripts/aux_fn.R"))


# folder for plots
output_dir_incd <-
  function(x = character(0))
    paste0("./validation_internal/plots_incd/", x)

output_dir_prvl <-
  function(x = character(0))
    paste0("./validation_internal/plots_prvl/", x)

output_dir_ftlt <-
  function(x = character(0))
    paste0("./validation_internal/plots_ftlt/", x)


#Setting the plot settings
theme_set(new = theme_few())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))


fl <- list.files("/mnt/storage_fast/output/hf_real/lifecourse",
                 "_lifecourse.csv$", full.names = TRUE)

out <- rbindlist(lapply(fl, fread))[year >= 13L & age >= 30]
# out <- out[mc == 1] # REMOVE!!!
setkey(out, mc, pid, year)

strata <- c("year", "age", "sex", "dimd", "sha", "ethnicity")

out[, 1e5 * weighted.mean(t2dm_prvl == 1, wt), keyby = .(mc)]


dl_model <- c("t2dm_prvl", "af_prvl", "chd_prvl", "dementia_prvl",
              "lung_ca_prvl", "colorect_ca_prvl", "breast_ca_prvl", "stroke_prvl",
              "obesity_prvl", "htn_prvl", "copd_prvl", "asthma_prvl",
              "ckd_prvl", "prostate_ca_prvl", "andep_prvl", "other_ca_prvl", "hf_prvl")
dl_model_new <- str_replace(dl_model, "_prvl", "")

setnames(out, dl_model, dl_model_new)

dl_cprd <- c("Anxiety_Depression", "Asthma", "Atrial Fibrillation", "CHD",
             "COPD" , "Chronic Kidney Disease", "Dementia", "Heart failure",
             "Hypertension" ,  "Obesity", "Other cancers" ,
             "Primary Malignancy_Breast", "Primary Malignancy_Colorectal",
             "Primary Malignancy_Lung", "Primary Malignancy_Prostate",
             "Stroke", "Type 2 Diabetes Mellitus")
dl_cprd_new <- c("andep", "asthma", "af", "chd", "copd", "ckd", "dementia", "hf",
             "htn", "obesity", "other_ca", "breast_ca", "colorect_ca",
             "lung_ca", "prostate_ca", "stroke", "t2dm")

dl_mort <- c(15, 9, 5, 2, 7, 10, 8, 17, 6, NA, 16, 14, 12, 11, 13, 3, 4)

dl_plot <- c("Anxiety & Depression", "Asthma", "Atrial Fibrillation", "CHD",
             "COPD" , "Chronic Kidney Disease", "Dementia", "Heart failure",
             "Hypertension" ,  "Obesity", "Other cancers" ,
             "Breast cancer", "Colorectal cancer",
             "Lung cancer", "Prostate cancer",
             "Stroke", "Type 2 Diabetes")

dl_lookup <- data.table(cprd = dl_cprd, sim = dl_cprd_new, mortality = dl_mort, plot = dl_plot)


#Incidence
incd_dt <- harmonise(read_fst(input_path("panel_short_inc.fst"),
                         as.data.table = TRUE)[gender != "I"])[between(age, 30, 100) &
                                                                 year < 2020 , .SD, .SDcols = c(strata, dl_cprd)]

setnames(incd_dt, dl_cprd, dl_cprd_new )

for (i in dl_cprd_new) {
  plot_disnm <- dl_lookup[sim == i, plot]

  if (i == "breast_ca") {
  setnames(out, paste(i), "disease")

  subtab1 <- out[sex == "women" & disease < 2 , weighted.mean(disease== 1, wt), keyby = .(year, mc)][
    ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

  setnames(out, "disease", paste(i) )

  setnames(incd_dt, paste(i), "disease")
  subtab2 <- incd_dt[sex == "women" & disease != 2, sum(disease == 1)/.N, keyby = year ]
  setnames(incd_dt, "disease", paste(i) )

  } else {
      if(i == "prostate_ca"){
        setnames(out, paste(i), "disease")

        subtab1 <- out[sex == "men" & disease < 2 , weighted.mean(disease== 1, wt), keyby = .(year, mc)][
          ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

        setnames(out, "disease", paste(i) )

        setnames(incd_dt, paste(i), "disease")
        subtab2 <- incd_dt[sex == "men" & disease != 2, sum(disease == 1)/.N, keyby = year ]
        setnames(incd_dt, "disease", paste(i) )
        } else{

    setnames(out, paste(i), "disease")

    subtab1 <- out[disease < 2 , weighted.mean(disease == 1, wt), keyby = .(year, mc)][
      ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

    setnames(out, "disease", paste(i) )

    setnames(incd_dt, paste(i), "disease")
    subtab2 <- incd_dt[disease != 2, sum(disease == 1)/.N, keyby = year ]
    setnames(incd_dt, "disease", paste(i) )}
    }

  subtab1 <- rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])

  ggplot(subtab1, aes(x = year, y = V1 * 1e5, col = Type)) +
    geom_point() +
    geom_line() +
    geom_smooth(linetype = "dotdash" ) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Incidence per 100,000 persons") +
    ggtitle(paste0("Incidence of ", plot_disnm, " per 100,000 persons"))+
    expand_limits(y = 0)
  ggsave(filename = output_dir_incd(paste0(i,".png")), scale = 1.5)

  rm(subtab1, subtab2)
}

rm(incd_dt)
gc()

#Prevelance
prvl_dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                             as.data.table = TRUE)[gender != "I"])[between(age, 20, 100) &
                                                                     year < 2020 , .SD, .SDcols = c(strata, dl_cprd)]

setnames(prvl_dt, dl_cprd, dl_cprd_new )

for (i in dl_cprd_new){

  plot_disnm <- dl_lookup[sim == i, plot]

  if(i == "breast_ca"){
    setnames(out, paste(i), "disease")

    subtab1 <- out[sex == "women", weighted.mean(disease != 0, wt), keyby = .(year, mc)][
      ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

    setnames(out, "disease", paste(i) )

    setnames(prvl_dt, paste(i), "disease")
    subtab2 <- prvl_dt[sex == "women", sum(disease != 0, na.rm = TRUE)/.N, keyby = year ]
    setnames(prvl_dt, "disease", paste(i) )
  } else {
    if(i == "prostate_ca"){
      setnames(out, paste(i), "disease")

      subtab1 <- out[sex == "men", weighted.mean(disease != 0, wt), keyby = .(year, mc)][
        ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

      setnames(out, "disease", paste(i) )

      setnames(prvl_dt, paste(i), "disease")
      subtab2 <- prvl_dt[sex == "men", sum(disease != 0, na.rm = TRUE)/.N, keyby = year ]
      setnames(prvl_dt, "disease", paste(i) )
    } else{
    setnames(out, paste(i), "disease")
  subtab1 <- out[, weighted.mean(disease != 0, wt), keyby = .(year, mc)][
    ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

  setnames(out, "disease", paste(i) )

  setnames(prvl_dt, paste(i), "disease")
  subtab2 <- prvl_dt[, sum(disease != 0, na.rm = TRUE)/.N, keyby = year ]
  setnames(prvl_dt, "disease", paste(i) )
  }}
  subtab1 <- rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])



  ggplot(subtab1, aes(x = year, y = V1 * 100, col = Type)) +
    geom_point() +
    geom_line() +
    geom_smooth(linetype = "dotdash" ) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Prevalence (%)") +
    ggtitle(paste0("Prevalence of ", plot_disnm)) +
    expand_limits(y = 0) +
    scale_colour_brewer(type = "qual")
  ggsave(filename = output_dir_prvl(paste0(i,".png")), scale = 1.5)

  rm(subtab1, subtab2)
}

rm(prvl_dt)

# Overall mortality
strata_mort <- c("death_cause", "diedinyr", "diedinCPRD")
mort_dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                              as.data.table = TRUE)[gender != "I"])[between(age, 20, 100) &
                                                                      year < 2020 , .SD, .SDcols = c(strata, dl_cprd, strata_mort)]

subtab1 <- out[, weighted.mean(all_cause_mrtl != 0, wt), keyby = .(year, mc)][
  ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

subtab2 <- mort_dt[diedinCPRD == 1 | is.na(diedinCPRD),
                  sum(diedinyr == 1, na.rm = TRUE)/.N, keyby = year ]

subtab1 <- rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])

ggplot(subtab1, aes(x = year, y = V1 * 100000, col = Type)) +
  geom_point() +
  geom_line() +
  geom_smooth(linetype = "dotdash" ) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "All-cause mortality (per 100,000 persons)",
                     limits = c(0,2000)) +
  ggtitle("All-cause mortality (per 100,000 persons)") +
  scale_colour_brewer(type = "qual")
ggsave(filename = output_dir_ftlt("all_cause_mort.png"), scale = 1.5)

rm(subtab1, subtab2)

# Cause-specific case-fatality

for (i in dl_cprd_new[dl_cprd_new != "obesity"]){
  cprd_disnm <- dl_lookup[sim == i, cprd]
  cod_code <- dl_lookup[sim == i, mortality]
  plot_disnm <- dl_lookup[sim == i, plot]


  if(i == "breast_ca"){
    subtab1 <- out[sex == "women" & get(i) != 0L, weighted.mean(all_cause_mrtl == cod_code, wt), keyby = .(year, mc)][
      ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

    subtab2 <- mort_dt[sex == "women" & get(cprd_disnm) != 0L,
                       sum(death_cause == cprd_disnm, na.rm = TRUE)/.N, keyby = year]
  } else {
    if(i == "prostate_ca"){

      subtab1 <- out[sex == "men"& get(i) != 0L, weighted.mean(all_cause_mrtl == cod_code, wt), keyby = .(year, mc)][
        ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]


      subtab2 <- mort_dt[sex == "men" & get(cprd_disnm) != 0L,
                         sum(death_cause == cprd_disnm, na.rm = TRUE)/.N, keyby = year]

      } else{
      subtab1 <- out[ get(i) != 0L, weighted.mean(all_cause_mrtl == cod_code, wt), keyby = .(year, mc)][
        ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]


      subtab2 <- mort_dt[get(cprd_disnm) != 0L,
                         sum(death_cause == cprd_disnm, na.rm = TRUE)/.N, keyby = year]

      }}
  subtab1 <- rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])



  ggplot(subtab1, aes(x = year, y = V1 * 100, col = Type)) +
    geom_point() +
    geom_line() +
    #geom_smooth(linetype = "dotdash") +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Case-fatality (per 100,000 persons)") +
    ggtitle(paste0("Case fatality of ", plot_disnm)) +
    expand_limits(y = 0) +
    scale_colour_brewer(type = "qual")
  ggsave(filename = output_dir_ftlt(paste0("cf_", i,".png")), scale = 1.5)

  rm(subtab1, subtab2)
}

rm(mort_dt)

tt <- fread("/mnt/storage_fast/output/hf_real/xps/xps.csv")
tt <- tt[agegrp20 == "All" & qimd == "All" &
           ethnicity == "All" & sha == "All",]
tt[year == 2014, lapply(.SD, mean), .SDcols = is.numeric, keyby = mc]
tt[, median(ets_curr_xps), keyby = year]
