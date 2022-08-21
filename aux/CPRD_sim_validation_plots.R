# CPRD Validation charts
library(data.table) # for fast data manipulation
library(fst)
library(ggplot2)
library(ggthemes)
library(scales)
library(RColorBrewer)
library(viridis)
# library(stringr)
library(CKutils)
library(future)

if (exists("IMPACTncd"))
  plan(multicore, workers = IMPACTncd$design$sim_prm$clusternumber)
if (exists("design"))
  plan(multicore, workers = design$sim_prm$clusternumber)

prvl_plot <- TRUE
incd_plot <- TRUE
mtlt_plot <- TRUE
ftlt_plot <- TRUE
out_pth <- "/mnt/storage_fast/output/hf_real/plots/validation/"
if (!dir.exists(out_pth)) dir.create(out_pth, recursive = TRUE)


#Setting the plot settings
theme_set(new = theme_economist())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))


# European standardised population 2013 (esp) weights
tt <- data.table(agegroup = agegrp_name(0, 99),
                 wt_esp  = c(1000, 4000, 5500, 5500, 5500, 6000, 6000, 6500,
                             7000, 7000, 7000, 7000, 6500, 6000, 5500, 5000,
                             4000, 2500, 1500, 800, 200))
esp <- CJ(agegroup = agegrp_name(0, 99, ),
          sex = c("men", "women"),
          dimd = c("1 most deprived", as.character(2:9), "10 least deprived")
)
absorb_dt(esp, tt)


# dl_plot <- c("Anxiety & Depression", "Asthma", "Atrial Fibrillation", "CHD",
#              "COPD" , "Chronic Kidney Disease", "Dementia", "Heart failure",
#              "Hypertension" ,  "Obesity", "Other cancers" ,
#              "Breast cancer", "Colorectal cancer",
#              "Lung cancer", "Prostate cancer",
#              "Stroke", "Type 2 Diabetes", "CMS > 0", "CMS > 1", "CMS > 1.5",
#              "CMS > 2")
#
# dl_lookup <- data.table(sim = dl_cprd_new, mortality = dl_mort, plot = dl_plot)


#strata <- c( "year", "agegroup"  , "sex", "dimd" )

cprdtab <- fread("./secure_data/summarytable.csv")[!agegroup %in% c("20-24", "25-29")]
# setnames(cprdtab, gsub("^colorect_", "colorectal_", names(cprdtab)))
# setnames(cprdtab, gsub("^psychos_", "psychosis_", names(cprdtab)))
# setnames(cprdtab, gsub("^constip_", "constipation_", names(cprdtab)))

absorb_dt(cprdtab, esp)
# cprdtab[, wt_esp := wt_esp * unique(wt_esp) / (wt_esp * popsize),
#         keyby = .(year, agegroup, sex, dimd)]
# cprdtab[, popsize_wtd := popsize * wt_esp]
cprdtab[, popsize_wtd := wt_esp]
cprdtab[, wt_esp := wt_esp / popsize]

tt <- fread("/mnt/storage_fast/output/hf_real/summaries/prvl_esp.csv.gz",
                  nrows = 0)
dl_model <- grep("_prvl$", copy(names(tt)), value = TRUE)
dl_model_new <- gsub("_prvl$", "", dl_model)



# Prvl ----
if (prvl_plot) {
  future({
    out_prvl <-
      fread("/mnt/storage_fast/output/hf_real/summaries/prvl_esp.csv.gz")[scenario == "sc0"]
    setnames(out_prvl, dl_model, dl_model_new)
    setkey(out_prvl, mc, year)
    for (i in dl_model_new) {
      # plot_disnm <- dl_lookup[sim == i, plot]

      if (i %in% c("dm", "ctdra", "cancer", "cms1st_cont"))
        next()
      setnames(out_prvl, i, "disease")
      setnames(cprdtab, paste0(i, "_prvl_N"), "disease")

      if (i == "breast_ca") {
        subtab1 <-
          out_prvl[sex == "women" , sum(disease) / sum(popsize), keyby = .(year, mc)][,  .(prvl = quantile(V1, 0.5)), keyby = year][, year := year + 2000]


        subtab2 <-
          cprdtab[sex == "women" , .(prvl = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = year]

      } else if (i == "prostate_ca") {
        subtab1 <-
          out_prvl[sex == "men" , sum(disease) / sum(popsize), keyby = .(year, mc)][,  .(prvl = quantile(V1, 0.5)), keyby = year][, year := year + 2000]
        subtab2 <-
          cprdtab[sex == "men", .(prvl = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = year]
      } else {
        subtab1 <-
          out_prvl[, sum(disease) / sum(popsize), keyby = .(year, mc)
                   ][,  .(prvl = quantile(V1, 0.5)), keyby = year
                     ][, year := year + 2000]
        subtab2 <-
          cprdtab[, .(prvl = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = year]

      }

      setnames(out_prvl, "disease", i)
      setnames(cprdtab, "disease", paste0(i, "_prvl_N"))


      subtab1 <-
        rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])

      ggplot(subtab1, aes(x = year, y = prvl, col = Type)) +
        geom_point() +
        geom_line() +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Prevalence", labels = scales::percent) +
        ggtitle(paste0("Prevalence of ", i, "\n(age-sex standardised)")) +
        expand_limits(y = 0) #+
      # scale_colour_brewer(type = "qual")

      ggsave(filename = paste0(out_pth, "prvl_", i, ".png"),
             scale = 3)
      # bigsumtab <- rbind(bigsumtab, subtab1[, condition := i])
    }
    rm(subtab1, subtab2)
  })
}


# Incd -----
if (incd_plot) {
  future({
    out_incd <-
      fread("/mnt/storage_fast/output/hf_real/summaries/incd_esp.csv.gz")[scenario == "sc0"]
    setnames(out_incd, dl_model, dl_model_new)
    setkey(out_incd, mc, year)

    # bigsumtab <- data.table() # Why do we need this?

    for (i in dl_model_new) {
      if (i %in% c("dm", "ctdra", "cancer", "cms1st_cont"))
        next()

      # plot_disnm <- dl_lookup[sim == i, plot]
      # tmp <- out_incd[, .SD, .SDcols = c( "agegrp", "mc", "year", "sex", "dimd", i)]
      setnames(out_incd, i, "disease")
      setnames(cprdtab, c(paste0(i, "_incd_N")), c("disease"))

      if (i == "breast_ca") {
        subtab1 <-
          out_incd[sex == "women" , sum(disease) / sum(popsize), keyby = .(year, mc)][,  .(incd = quantile(V1, 0.5)), keyby = year][, year := year + 2000]

        subtab2 <-
          cprdtab[sex == "women" , .(incd = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = year]

      } else if (i == "prostate_ca") {
        subtab1 <-
          out_incd[sex == "men" , sum(disease) / sum(popsize), keyby = .(year, mc)][,  .(incd = quantile(V1, 0.5)), keyby = year][, year := year + 2000]
        subtab2 <-
          cprdtab[sex == "men", .(incd = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = year]
      } else {
        subtab1 <-
          out_incd[, sum(disease) / sum(popsize), keyby = .(year, mc)][,  .(incd = quantile(V1, 0.5)), keyby = year][, year := year + 2000]
        subtab2 <-
          cprdtab[, .(incd = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = year]

      }

      setnames(out_incd, "disease", paste(i))
      setnames(cprdtab, c("disease"), c(paste0(i, "_incd_N")))


      subtab1 <-
        rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])

      ggplot(subtab1, aes(
        x = year,
        y = incd * 100000,
        col = Type
      )) +
        geom_point() +
        geom_line() +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ",
          i,
          " per 100,000 persons\n(age-sex standardised"
        )) +
        expand_limits(y = 0)
      ggsave(filename = paste0(out_pth, "incd_", i, ".png"),
             scale = 3)
      # bigsumtab <- rbind(bigsumtab, subtab1[,condition := i] )
    }
    rm(subtab1, subtab2)
  })
}


# Overall mrtl ----
if (mtlt_plot) {
  future({
    out_mrtl <-
      fread("/mnt/storage_fast/output/hf_real/summaries/mrtl_esp.csv.gz")[scenario == "sc0"]
    setkey(out_mrtl, mc, year)

    subtab1 <-
      out_mrtl[, .(sum(all_cause_mrtl) / sum(popsize)), keyby = .(year, mc)][, .(all_ftlt = quantile(V1, 0.5)), keyby = year][, year := year + 2000]

    subtab2 <-
      cprdtab[, .(all_ftlt = sum(all_ftlt_N * wt_esp) / sum(popsize_wtd)), keyby = year]

    subtab1 <-
      rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])

    ggplot(subtab1, aes(
      x = year,
      y = all_ftlt * 100000,
      col = Type
    )) +
      geom_point() +
      geom_line() +
      geom_smooth(linetype = "dotdash", se = FALSE) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "All-cause mortality (per 100,000 persons, age-sex standardised)",
                         limits = c(0, 2000)) +
      ggtitle("All-cause mortality (per 100,000 persons)") +
      scale_colour_brewer(type = "qual")
    ggsave(filename = paste0(out_pth, "all_cause_mrtl.png"),
           scale = 3)

    rm(subtab1, subtab2)
  })
}


# Cause-specific mrtl ----
if (ftlt_plot) {
  future({
    tt <- grep("_ftlt_N$", names(cprdtab), value = TRUE)
    # Ensure all_ftlt is always the first
    tt <- tt[order(match(tt, "all_ftlt_N"))]
    cprdtab[, nonmodelled_ftlt_N := Reduce(`-`, .SD), .SDcols = tt]

    out_ftlt <-
      fread("/mnt/storage_fast/output/hf_real/summaries/dis_mrtl_esp.csv.gz")[scenario == "sc0"]
    setkey(out_ftlt, mc, year)

    dsnam <- names(out_ftlt)[!names(out_ftlt) %in% c("agegrp",
                                                     "mc",
                                                     "scenario",
                                                     "year",
                                                     "sex",
                                                     "dimd",
                                                     "ethnicity",
                                                     "sha",
                                                     "popsize")]

    for (i in dsnam) {
      # if (i %in% c("dm", "ctdra", "cancer", "cms1st_cont"))
      #   next()

      setnames(out_ftlt, i, "disease")
      setnames(cprdtab, c(paste0(i, "_ftlt_N")), c("disease"))

      if (i == "breast_ca") {
        subtab1 <-
          out_ftlt[sex == "women" , sum(disease) / sum(popsize), keyby = .(year, mc)][,  .(ftlt = quantile(V1, 0.5)), keyby = year][, year := year + 2000]

        subtab2 <-
          cprdtab[sex == "women" , .(ftlt = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = year]

      } else if (i == "prostate_ca") {
        subtab1 <-
          out_ftlt[sex == "men" , sum(disease) / sum(popsize), keyby = .(year, mc)][,  .(ftlt = quantile(V1, 0.5)), keyby = year][, year := year + 2000]
        subtab2 <-
          cprdtab[sex == "men", .(ftlt = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = year]
      } else {
        subtab1 <-
          out_ftlt[, sum(disease) / sum(popsize), keyby = .(year, mc)][,  .(ftlt = quantile(V1, 0.5)), keyby = year][, year := year + 2000]
        subtab2 <-
          cprdtab[, .(ftlt = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = year]

      }

      setnames(out_ftlt, "disease", paste(i))
      setnames(cprdtab, c("disease"), c(paste0(i, "_ftlt_N")))


      subtab1 <-
        rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])

      ggplot(subtab1, aes(
        x = year,
        y = ftlt * 100000,
        col = Type
      )) +
        geom_point() +
        geom_line() +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0(
            "Cause-spesific mortality of ",
            i,
            " per 100,000 persons\n(age-sex standardised"
          )
        ) +
        expand_limits(y = 0)
      ggsave(filename = paste0(out_pth, "ftlt_", i, ".png"),
             scale = 3)
    }
    rm(subtab1, subtab2)
  })
}
