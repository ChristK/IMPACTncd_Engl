# CPRD Validation charts

library(data.table) # for fast data manipulation
library(fst)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(stringr)
library(CKutils)
theme_set(new = theme_economist())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))

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

# European standardised population 2013 (esp) weights
tt <- data.table(agegrp = agegrp_name(0, 99),
                  wt_esp  = c(1000, 4000, 5500, 5500, 5500, 6000, 6000, 6500,
                              7000, 7000, 7000, 7000, 6500, 6000, 5500, 5000,
                              4000, 2500, 1500, 800, 200))
esp <- CJ(agegrp = agegrp_name(0, 99),
          sex = c("men", "women"),
          dimd = c("1 most deprived", as.character(2:9), "10 least deprived")
          )
absorb_dt(esp, tt)



fl <- list.files("/mnt/storage_fast/output/hf_real/lifecourse",
                 "_lifecourse.csv.gz$", full.names = TRUE)


out <- rbindlist(lapply(fl, fread))
out[, dimd := factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived"))]

# out[, wt_esp := wt_esp * unique(wt_esp) / sum(wt_esp),
#     keyby = .(year, agegrp, sex, dimd)]
setkey(out, mc, pid, year)

strata <- c("year", "age", "sex", "dimd", "sha", "ethnicity")


dl_model <- c("t2dm_prvl", "af_prvl", "chd_prvl", "dementia_prvl",
              "lung_ca_prvl", "colorect_ca_prvl", "breast_ca_prvl", "stroke_prvl",
              "obesity_prvl", "htn_prvl", "copd_prvl", "asthma_prvl",
              "ckd45_prvl", "prostate_ca_prvl", "andep_prvl", "other_ca_prvl", "hf_prvl")
dl_model_new <- str_replace(dl_model, "_prvl", "")

setnames(out, dl_model, dl_model_new)

dl_cprd <- c("Anxiety_Depression", "Asthma", "Atrial Fibrillation", "CHD",
             "COPD" , "ckd45", "Dementia", "Heart failure",
             "Hypertension" ,  "Obesity", "Other cancers",
             "Primary Malignancy_Breast", "Primary Malignancy_Colorectal",
             "Primary Malignancy_Lung", "Primary Malignancy_Prostate",
             "Stroke", "Type 2 Diabetes Mellitus")
dl_cprd_new <- c("andep", "asthma", "af", "chd", "copd", "ckd45", "dementia", "hf",
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

# weights <- out[year < 20, sum(wt), keyby = .(age, sex, dimd, mc)][, .(wt = mean(V1)), keyby = .(age, sex, dimd)]

#Incidence plots ----
incd_dt <- harmonise(read_fst(input_path("panel_short_inc.fst"),
                         as.data.table = TRUE)[gender != "I"])[between(age, 30, 99) &
                                                                 year < 2020 , .SD, .SDcols = c(strata, dl_cprd)]
incd_dt[, dimd := factor(dimd, 10:1, levels(out$dimd))]
to_agegrp(incd_dt, 5, 99)
absorb_dt(incd_dt, esp)
incd_dt[, wt_esp := wt_esp * unique(wt_esp) / sum(wt_esp), keyby = .(year, agegrp, sex, dimd)]

setnames(incd_dt, dl_cprd, dl_cprd_new)

for (i in dl_cprd_new) {
  plot_disnm <- dl_lookup[sim == i, plot]

  if (i == "breast_ca") {
    setnames(out, paste(i), "disease")

    subtab1 <-
      out[sex == "women" &
            disease < 2 , weighted.mean(disease == 1, wt_esp), keyby = .(year, mc)
          ][,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

    setnames(out, "disease", paste(i))

    setnames(incd_dt, paste(i), "disease")
    subtab2 <-
      incd_dt[sex == "women" &
                disease < 2, weighted.mean(disease == 1, wt_esp), keyby = year]
    setnames(incd_dt, "disease", paste(i))

  } else if (i == "prostate_ca") {
    setnames(out, i, "disease")

    subtab1 <- out[sex == "men" & disease < 2 ,
                   weighted.mean(disease == 1, wt_esp), keyby = .(year, mc)
                   ][, quantile(V1, 0.5), keyby = year][, year := year + 2000]

    setnames(out, "disease", i)

    setnames(incd_dt, i, "disease")
    subtab2 <-
      incd_dt[sex == "men" &
                disease < 2, weighted.mean(disease == 1, wt_esp),
              keyby = year]
    setnames(incd_dt, "disease", i)
  } else {
    setnames(out, i, "disease")

    subtab1 <-
      out[disease < 2 , weighted.mean(disease == 1, wt_esp), keyby = .(year, mc)
          ][,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

    setnames(out, "disease", paste(i))

    setnames(incd_dt, paste(i), "disease")
    subtab2 <-
      incd_dt[disease < 2, weighted.mean(disease == 1, wt_esp), keyby = year]
    setnames(incd_dt, "disease", paste(i))
  }

  subtab1 <- rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])

  ggplot(subtab1, aes(x = year, y = V1 * 1e5, col = Type)) +
    geom_point() +
    # geom_line() +
    geom_smooth(linetype = "dotdash", se = FALSE) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Incidence per 100,000 persons") +
    ggtitle(paste0("Incidence of ", plot_disnm, " per 100,000 persons"))+
    expand_limits(y = 0)
  ggsave(filename = output_dir_incd(paste0(i,".png")), scale = 1.5)

  rm(subtab1, subtab2)
}

rm(incd_dt)
gc()

# Prevalence plots ----
prvl_dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                             as.data.table = TRUE)[gender != "I"])[between(age, 30, 99) &
                                                                     year < 2020 , .SD, .SDcols = c(strata, dl_cprd)]
prvl_dt[, dimd := factor(dimd, 10:1, levels(out$dimd))]
to_agegrp(prvl_dt, 5, 99)
absorb_dt(prvl_dt, esp)
prvl_dt[, wt_esp := wt_esp * unique(wt_esp) / sum(wt_esp), keyby = .(year, agegrp, sex, dimd)]

setnames(prvl_dt, dl_cprd, dl_cprd_new )

for (i in dl_cprd_new) {
  plot_disnm <- dl_lookup[sim == i, plot]

  if (i == "breast_ca") {
    setnames(out, i, "disease")

    subtab1 <-
      out[sex == "women", weighted.mean(disease != 0, wt_esp), keyby = .(year, mc)
          ][,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

    setnames(out, "disease", paste(i))

    setnames(prvl_dt, i, "disease")
    subtab2 <-
      prvl_dt[sex == "women", weighted.mean(disease != 0, wt_esp), keyby = year]
    setnames(prvl_dt, "disease", paste(i))
  } else if (i == "prostate_ca") {
    setnames(out, paste(i), "disease")

    subtab1 <-
      out[sex == "men", weighted.mean(disease != 0, wt_esp), keyby = .(year, mc)
          ][,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

    setnames(out, "disease", paste(i))

    setnames(prvl_dt, paste(i), "disease")
    subtab2 <-
      prvl_dt[sex == "men", weighted.mean(disease != 0, wt_esp), keyby = year]
    setnames(prvl_dt, "disease", paste(i))
  } else {
    setnames(out, paste(i), "disease")
    subtab1 <-
      out[, weighted.mean(disease != 0, wt_esp), keyby = .(year, mc)
          ][,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

    setnames(out, "disease", paste(i))

    setnames(prvl_dt, paste(i), "disease")
    subtab2 <-
      prvl_dt[, weighted.mean(disease > 0, wt_esp, na.rm = TRUE), keyby = year]
    setnames(prvl_dt, "disease", paste(i))
  }
  subtab1 <-
    rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])



  ggplot(subtab1, aes(x = year, y = V1 * 100, col = Type)) +
    geom_point() +
    # geom_line() +
    geom_smooth(linetype = "dotdash", se = FALSE) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "Prevalence (%)") +
    ggtitle(paste0("Prevalence of ", plot_disnm)) +
    expand_limits(y = 0) +
    scale_colour_brewer(type = "qual")
  ggsave(filename = output_dir_prvl(paste0(i, ".png")), scale = 1.5)

  rm(subtab1, subtab2)
}

rm(prvl_dt)

# Overall mortality plots ----
strata_mort <- c("death_cause", "diedinyr", "diedinCPRD")
mort_dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                              as.data.table = TRUE)[gender != "I"])[between(age, 30, 99) &
                                                                      year < 2020 , .SD, .SDcols = c(strata, dl_cprd, strata_mort)]
mort_dt[, dimd := factor(dimd, 10:1, levels(out$dimd))]
to_agegrp(mort_dt, 5, 99)
absorb_dt(mort_dt, esp)
mort_dt[, wt_esp := wt_esp * unique(wt_esp) / sum(wt_esp), keyby = .(year, agegrp, sex, dimd)]
setnafill(mort_dt, "c", 0L, cols = c("diedinCPRD", "diedinyr"))
subtab1 <- out[, weighted.mean(all_cause_mrtl > 0, wt_esp), keyby = .(year, mc)
               ][,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

subtab2 <- mort_dt[, # TODO should be diedinyr == 1 but looks strange
                   weighted.mean(diedinCPRD == 1, wt_esp), keyby = year]

subtab1 <- rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])

ggplot(subtab1, aes(x = year, y = V1 * 100000, col = Type)) +
  geom_point() +
  # geom_line() +
  geom_smooth(linetype = "dotdash", se = FALSE) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = "All-cause mortality (per 100,000 persons)",
                     limits = c(0,NA)) +
  ggtitle("All-cause mortality (per 100,000 persons)") +
  scale_colour_brewer(type = "qual")
ggsave(filename = output_dir_ftlt("all_cause_mort.png"), scale = 1.5, width = 16, height = 9)

rm(subtab1, subtab2)

# Cause-specific case-fatality plots ----
mort_dt[is.na(death_cause), death_cause := "alive"]
for (i in dl_cprd_new[dl_cprd_new != "obesity"]){
  cprd_disnm <- dl_lookup[sim == i, cprd]
  cod_code <- dl_lookup[sim == i, mortality]
  plot_disnm <- dl_lookup[sim == i, plot]


  if(i == "breast_ca"){
    subtab1 <- out[sex == "women" & get(i) > 0L, weighted.mean(all_cause_mrtl == cod_code, wt_esp), keyby = .(year, mc)][
      ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]

    subtab2 <- mort_dt[sex == "women" & get(cprd_disnm) > 0L,
                       weighted.mean(death_cause == cprd_disnm, wt_esp), keyby = year]
  } else {
    if(i == "prostate_ca"){

      subtab1 <- out[sex == "men"& get(i) > 0L, weighted.mean(all_cause_mrtl == cod_code, wt_esp), keyby = .(year, mc)][
        ,  quantile(V1, 0.5), keyby = year][, year := year + 2000]


      subtab2 <- mort_dt[sex == "men" & get(cprd_disnm) > 0L,
                         weighted.mean(death_cause == cprd_disnm, wt_esp), keyby = year]

      } else{
        subtab1 <-
          out[get(i) > 0L, weighted.mean(all_cause_mrtl == cod_code, wt_esp),
              keyby = .(year, mc)][,  quantile(V1, 0.5), keyby = year
                                   ][, year := year + 2000]


      subtab2 <- mort_dt[get(cprd_disnm) > 0L,
                         weighted.mean(death_cause == cprd_disnm, wt_esp),
                         keyby = year]

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

