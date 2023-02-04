# CPRD Validation charts


# 2.	File path for the ONS mortality data on line 740: this will need updating to wherever you save it

prvl_plot <- TRUE
incd_plot <- TRUE
mtlt_plot <- TRUE
ftlt_plot <- TRUE

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

if (exists("IMPACTncd")) {
  plan(multicore, workers = IMPACTncd$design$sim_prm$clusternumber)
  mdl_rslts_pth <- paste0(IMPACTncd$design$sim_prm$output_dir, "/summaries/")
  out_pth <- paste0(IMPACTncd$design$sim_prm$output_dir, "/plots/validation/")

} else if (exists("design")) {
  plan(multicore, workers = design$sim_prm$clusternumber)
  mdl_rslts_pth <- paste0(design$sim_prm$output_dir, "/summaries/")
  out_pth <- paste0(design$sim_prm$output_dir, "/plots/validation/")

} else stop()

if (!dir.exists(out_pth)) dir.create(out_pth, recursive = TRUE)


# Setting the plot settings
theme_set(new = theme_economist())
theme_update(axis.text.x = element_text(size = 8), plot.title = element_text(hjust = 0.5))


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

esp90 <- copy(esp)
esp90[agegroup %in% c("90-94", "95-99"), agegroup := "90+"]
esp90 <- esp90[, .(wt_esp = sum(wt_esp)), by = .(agegroup, sex, dimd)]

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

cprdtab <- fread("./secure_data/summarytable.csv")[!agegroup %in% c("20-24", "25-29")][, wt_esp := NULL]
# setnames(cprdtab, gsub("^colorect_", "colorectal_", names(cprdtab)))
# setnames(cprdtab, gsub("^psychos_", "psychosis_", names(cprdtab)))
# setnames(cprdtab, gsub("^constip_", "constipation_", names(cprdtab)))

absorb_dt(cprdtab, esp)
# cprdtab[, wt_esp := wt_esp * unique(wt_esp) / (wt_esp * popsize),
#         keyby = .(year, agegroup, sex, dimd)]
# cprdtab[, popsize_wtd := popsize * wt_esp]
cprdtab[, popsize_wtd := wt_esp]
cprdtab[, wt_esp := wt_esp / popsize]

agetab <- data.table(agegroup = agegrp_name(30, 99, grp_width = 5),
                     agegroup10 = c("30-39", "30-39", "40-49", "40-49", "50-59",
                                    "50-59", "60-69", "60-69", "70-79", "70-79",
                                    "80-89", "80-89", "90-99", "90-99"))
agetab[, `:=` (agegroup = factor(agegroup),
               agegroup10 = factor(agegroup10))]
cprdtab <- cprdtab[agetab, on = 'agegroup' ]
cprdtab[, dimd := factor(dimd,
                          levels = c("1 most deprived", "2", "3", "4", "5",
                                     "6", "7", "8", "9", "10 least deprived"))]

tt <- fread(paste0(mdl_rslts_pth, "prvl_esp.csv.gz"), nrows = 0)
dl_model <- grep("_prvl$", names(tt), value = TRUE)
dl_model_new <- gsub("_prvl$", "", dl_model)
dsnm <- dl_model_new[!dl_model_new %in% c("dm", "ctdra", "cancer", "cms1st_cont")]



# General summarising function
prep <- function(i, data) {
  if (i == "breast_ca") {
  subtab1 <-
    data[sex == "women" ,
         sum(disease) / sum(popsize),
         keyby = strata_mc][
           ,  .(outcome = quantile(V1, 0.5)), keyby = outstrata][
             , year := year + 2000]

  subtab2 <-
    cprdtab[sex == "women" ,
            .(outcome = sum(disease * wt_esp) / sum(popsize_wtd)),
            keyby = outstrata]

} else if (i == "prostate_ca") {
  subtab1 <-
    data[sex == "men" ,
         sum(disease) / sum(popsize),
         keyby = strata_mc][
           ,  .(outcome = quantile(V1, 0.5)), keyby = outstrata][
             , year := year + 2000]

  subtab2 <-
    cprdtab[sex == "men" ,
            .(outcome = sum(disease * wt_esp) / sum(popsize_wtd)),
            keyby = outstrata]
} else {
  subtab1 <-
    data[, sum(disease) / sum(popsize), keyby = strata_mc
    ][,  .(outcome = quantile(V1, 0.5)), keyby = outstrata
    ][, year := year + 2000]
  subtab2 <-
    cprdtab[, .(outcome = sum(disease * wt_esp) / sum(popsize_wtd)), keyby = outstrata]
  }
  subtab1 <-
    rbind(subtab1[, Type := "Simulated"], subtab2[, Type := "CPRD"])
}


#For prev
prep_prev <- function(i, strata_mc, outstrata,
                      out_prvl = out_prvl,
                      cprdtab = cprdtab){

  setnames(out_prvl, i, "disease")
  setnames(cprdtab, paste0(i, "_prvl_N"), "disease")

  subtab1 <- prep(i, out_prvl)

  setnames(out_prvl, "disease", i)
  setnames(cprdtab, "disease", paste0(i, "_prvl_N"))
  subtab1
}

#For incd
prep_incd <- function(i, strata_mc, outstrata){

  setnames(out_incd, i, "disease")
  setnames(cprdtab, paste0(i, "_incd_N"), "disease")

  subtab1 <- prep(i, out_incd)

  setnames(out_incd, "disease", i)
  setnames(cprdtab, "disease", paste0(i, "_incd_N"))
  subtab1
}

#For cf
prep_cftlt <- function(i, strata_mc, outstrata){
  setnames(out_ftlt, i, "disease")
  setnames(cprdtab, c(paste0(i, "_ftlt_N")), c("disease"))

  subtab1 <- prep(i, out_ftlt)

   setnames(out_ftlt, "disease", i)
  setnames(cprdtab, "disease", paste0(i, "_ftlt_N"))
  subtab1
}

#For all-cause mortatlity
prep_mrtl <- function(strata_mc, outstrata){
  subtab1 <-
    out_mrtl[, .(sum(all_cause_mrtl) / sum(popsize)), keyby = strata_mc][
      , .(all_ftlt = quantile(V1, 0.5)), keyby = outstrata][, year := year + 2000]

  setkey(cprdtab, year, sex, dimd)
  subtab2 <-
    cprdtab[, .(all_ftlt = sum(all_ftlt_N * wt_esp) / sum(popsize_wtd)), keyby = outstrata]

  subtab3 <-
    onsmrtl[, .(all_ftlt = sum(deaths) / sum(popsize)), keyby = outstrata]


  rbind(subtab1[, Type := "Simulated"],
        subtab2[, Type := "CPRD"],
        subtab3[, Type := "ONS"])
}




# Prvl ----
if (prvl_plot) {
    out_prvl <- fread(paste0(mdl_rslts_pth,"prvl_esp.csv.gz"))[scenario == "sc0"]
    setnames(out_prvl, dl_model, dl_model_new)

    out_prvl[agetab, on = c("agegrp" = "agegroup"), agegroup10 := i.agegroup10]

    out_prvl[, dimd := factor(dimd,
                             levels = c("1 most deprived", "2", "3", "4", "5",
                                            "6", "7", "8", "9", "10 least deprived"))]



    # Overall
    future({
      strata_mc <- c("year", "mc")
      outstrata <- c("year")
      out_pth_prvl <- paste0(out_pth, "prvl/", paste(outstrata, collapse = "_"), "/")
      if (!dir.exists(out_pth_prvl)) dir.create(out_pth_prvl, recursive = TRUE)

      setkey(out_prvl, mc, year)
      for (i in dsnm) {
        subtab1 <- prep_prev(i, strata_mc, outstrata)

        ggplot(subtab1, aes(x = year, y = outcome, col = Type)) +
          geom_point(size = 0.5) +
          geom_line() +
          geom_smooth(linetype = "dotdash", se = FALSE) +
          scale_x_continuous(name = "Year") +
          scale_y_continuous(name = "Prevalence", labels = scales::percent) +
          ggtitle(paste0("Prevalence of ", i, "\n(age-sex standardised)")) +
          expand_limits(y = 0) #+
        # scale_colour_brewer(type = "qual")

        ggsave(filename = paste0(out_pth_prvl, "prvl_", i, ".png"),
               scale = 3)
      }
    })

    # by dimd
    future({
      strata_mc <- c("year", "mc", "dimd")
      outstrata <- c("year", "dimd")
      out_pth_prvl <- paste0(out_pth, "prvl/", paste(outstrata, collapse = "_"), "/")
      if (!dir.exists(out_pth_prvl)) dir.create(out_pth_prvl, recursive = TRUE)

      for (i in dsnm) {
        subtab1 <- prep_prev(i, strata_mc, outstrata)

        ggplot(subtab1, aes(x = year, y = outcome, col = Type)) +
          geom_point(size = 0.5) +
          geom_line() +
          facet_wrap( ~ dimd, scales = "free") +
          geom_smooth(linetype = "dotdash", se = FALSE) +
          scale_x_continuous(name = "Year") +
          scale_y_continuous(name = "Prevalence", labels = scales::percent) +
          ggtitle(paste0("Prevalence of ", i, "\n(age-sex standardised)")) +
          expand_limits(y = 0)  +
          theme(
            strip.text = element_text(size = 6),
            axis.text.x = element_text(
              angle = 90,
              vjust = 0.5,
              hjust = 1
            ),
            axis.text.y = element_text(size = 8)
          ) #+
        # scale_colour_brewer(type = "qual")

        ggsave(filename = paste0(out_pth_prvl, "prvl_dimd_", i, ".png"),
               scale = 3)
      }
    })


    # by sex
    future({
      strata_mc <- c("year", "mc", "sex")
      outstrata <- c("year", "sex")
      out_pth_prvl <- paste0(out_pth, "prvl/", paste(outstrata, collapse = "_"), "/")
      if (!dir.exists(out_pth_prvl)) dir.create(out_pth_prvl, recursive = TRUE)

      for (i in dsnm) {
        subtab1 <- prep_prev(i, strata_mc, outstrata)

        ggplot(subtab1, aes(
          x = year,
          y = outcome,
          col = sex,
          linetype = Type
        )) +
          geom_point(size = 0.5) +
          geom_line() +
          geom_smooth(linetype = "dotdash", se = FALSE) +
          scale_x_continuous(name = "Year") +
          scale_y_continuous(name = "Prevalence", labels = scales::percent) +
          ggtitle(paste0("Prevalence of ", i, "\n(age standardised)")) +
          expand_limits(y = 0) #+
        # scale_colour_brewer(type = "qual")

        ggsave(filename = paste0(out_pth_prvl, "prvl_sex_", i, ".png"),
               scale = 3)
      }
    })

    # by agegroup
    future({
      strata_mc <- c("year", "mc", "agegroup10")
      outstrata <- c("year", "agegroup10")
      out_pth_prvl <- paste0(out_pth, "prvl/", paste(outstrata, collapse = "_"), "/")
      if (!dir.exists(out_pth_prvl)) dir.create(out_pth_prvl, recursive = TRUE)

      for (i in dsnm) {
        subtab1 <- prep_prev(i, strata_mc, outstrata)

        ggplot(subtab1, aes(x = year, y = outcome, col = Type)) +
          geom_point(size = 0.5) +
          geom_line() +
          facet_wrap( ~ agegroup10, scales = "free") +
          geom_smooth(linetype = "dotdash", se = FALSE) +
          scale_x_continuous(name = "Year") +
          scale_y_continuous(name = "Prevalence", labels = scales::percent) +
          ggtitle(paste0("Prevalence of ", i, "\n(sex standardised)")) +
          expand_limits(y = 0)  +
          theme(
            strip.text = element_text(size = 6),
            axis.text.x = element_text(
              angle = 90,
              vjust = 0.5,
              hjust = 1
            ),
            axis.text.y = element_text(size = 8)
          )#+
        # scale_colour_brewer(type = "qual")

        ggsave(
          filename = paste0(out_pth_prvl, "prvl_agegroup_", i, ".png"),
          scale = 3
        )
      }
    })

    # by sex & imd
    future({
      strata_mc <- c("year", "mc", "sex", "dimd")
      outstrata <- c("year", "sex", "dimd")
      out_pth_prvl <- paste0(out_pth, "prvl/", paste(outstrata, collapse = "_"), "/")
      if (!dir.exists(out_pth_prvl)) dir.create(out_pth_prvl, recursive = TRUE)

      for (i in dsnm) {
        subtab1 <- prep_prev(i, strata_mc, outstrata)

        ggplot(subtab1, aes(x = year, y = outcome, col = Type)) +
          geom_point(size = 0.5) +
          geom_line() +
          facet_grid(dimd ~ sex, scales = "free") +
          geom_smooth(linetype = "dotdash", se = FALSE) +
          scale_x_continuous(name = "Year") +
          scale_y_continuous(name = "Prevalence", labels = scales::percent) +
          ggtitle(paste0("Prevalence of ", i, "\n(age standardised)")) +
          expand_limits(y = 0)  +
          theme(
            strip.text = element_text(size = 6),
            axis.text.x = element_text(
              angle = 90,
              vjust = 0.5,
              hjust = 1
            ),
            axis.text.y = element_text(size = 8)
          )#+
        # scale_colour_brewer(type = "qual")

        ggsave(
          filename = paste0(out_pth_prvl, "prvl_sex_dimd_", i, ".png"),
          scale = 3
        )
      }
    })

    # by sex & agegroup
    future({
      strata_mc <- c("year", "mc", "sex", "agegroup10")
      outstrata <- c("year", "sex", "agegroup10")
      out_pth_prvl <- paste0(out_pth, "prvl/", paste(outstrata, collapse = "_"), "/")
      if (!dir.exists(out_pth_prvl)) dir.create(out_pth_prvl, recursive = TRUE)

      for (i in dsnm) {
        subtab1 <- prep_prev(i, strata_mc, outstrata)

        ggplot(subtab1, aes(x = year, y = outcome, col = Type)) +
          geom_point(size = 0.5) +
          geom_line() +
          facet_grid(agegroup10 ~ sex, scales = "free") +
          geom_smooth(linetype = "dotdash", se = FALSE) +
          scale_x_continuous(name = "Year") +
          scale_y_continuous(name = "Prevalence", labels = scales::percent) +
          ggtitle(paste0("Prevalence of ", i)) +
          expand_limits(y = 0) +
          theme(
            strip.text = element_text(size = 6),
            axis.text.x = element_text(
              angle = 90,
              vjust = 0.5,
              hjust = 1
            ),
            axis.text.y = element_text(size = 8)
          ) #+
        # scale_colour_brewer(type = "qual")

        ggsave(
          filename = paste0(out_pth_prvl, "prvl_sex_agegroup_", i, ".png"),
          scale = 3
        )
      }
    })

    # by dimd & agegroup
    future({
      strata_mc <- c("year", "mc", "dimd", "agegroup10")
      outstrata <- c("year", "dimd", "agegroup10")
      out_pth_prvl <- paste0(out_pth, "prvl/", paste(outstrata, collapse = "_"), "/")
      if (!dir.exists(out_pth_prvl)) dir.create(out_pth_prvl, recursive = TRUE)

      for (i in dsnm) {
        subtab1 <- prep_prev(i, strata_mc, outstrata)

        ggplot(subtab1, aes(x = year, y = outcome, col  = Type)) +
          #geom_point(size = 0.5) + #removing because it makes it v busy
          geom_line(size = 0.5) +
          facet_grid(agegroup10 ~ dimd, scales = "free") +
          geom_smooth(linetype = "dotdash", se = FALSE) +
          scale_x_continuous(name = "Year") +
          scale_y_continuous(name = "Prevalence", labels = scales::percent) +
          ggtitle(paste0("Prevalence of ", i, "\n(sex standardised)")) +
          expand_limits(y = 0) +
          theme(
            strip.text = element_text(size = 8),
            axis.text.x = element_text(
              angle = 90,
              vjust = 0.5,
              hjust = 1
            ),
            axis.text.y = element_text(size = 10)
          )
        #+
        # scale_colour_brewer(type = "qual")

        ggsave(
          filename = paste0(out_pth_prvl, "prvl_dimd_agegroup_", i, ".png"),
          scale = 3
        )
      }
    })

    # by sex & dimd & agegroup
    future({
      strata_mc <- c("year", "mc", "sex", "dimd", "agegroup10")
      outstrata <- c("year", "sex", "dimd", "agegroup10")
      out_pth_prvl <- paste0(out_pth, "prvl/", paste(outstrata, collapse = "_"), "/")
      if (!dir.exists(out_pth_prvl)) dir.create(out_pth_prvl, recursive = TRUE)

      for (i in dsnm) {
        subtab1 <- prep_prev(i, strata_mc, outstrata)

        if (i != "breast_ca") {
          ggplot(subtab1[sex == "men"], aes(x = year, y = outcome, col = Type)) +
            #geom_point(size = 0.5) + #removing because it makes it v busy
            geom_line(size = 0.5) +
            facet_grid(agegroup10 ~ dimd , scales = "free") +
            geom_smooth(linetype = "dotdash", se = FALSE) +
            scale_x_continuous(name = "Year") +
            scale_y_continuous(name = "Prevalence", labels = scales::percent) +
            ggtitle(paste0("Prevalence of ", i, " in men")) +
            expand_limits(y = 0) +
            theme(
              strip.text = element_text(size = 8),
              axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 1
              ),
              axis.text.y = element_text(size = 8)
            )#+
          # scale_colour_brewer(type = "qual")

          ggsave(
            filename = paste0(out_pth_prvl, "prvl_dimd_agegroup_men_", i, ".png"),
            scale = 3
          )
        }

        if (i != "prostate_ca") {
          ggplot(subtab1[sex == "women"], aes(x = year, y = outcome, col = Type)) +
            #geom_point(size = 0.5) + #removing because it makes it v busy
            geom_line(size = 0.5) +
            facet_grid(agegroup10 ~ dimd , scales = "free") +
            geom_smooth(linetype = "dotdash", se = FALSE) +
            scale_x_continuous(name = "Year") +
            scale_y_continuous(name = "Prevalence", labels = scales::percent) +
            ggtitle(paste0("Prevalence of ", i, " in women")) +
            expand_limits(y = 0) +
            theme(
              strip.text = element_text(size = 8),
              axis.text.x = element_text(
                angle = 90,
                vjust = 0.5,
                hjust = 1
              ),
              axis.text.y = element_text(size = 8)
            ) #+
          # scale_colour_brewer(type = "qual")

          ggsave(
            filename = paste0(
              out_pth_prvl,
              "prvl_sex_dimd_agegroup_women_",
              i,
              ".png"
            ),
            scale = 3
          )
        }
      }
    })

    rm(subtab1)
 # })
}


# Incd -----
if (incd_plot) {
 # future({
    out_incd <-
      fread(paste0(mdl_rslts_pth,"incd_esp.csv.gz"))[scenario == "sc0"]
    setnames(out_incd, dl_model, dl_model_new)
    setkey(out_incd, mc, year)

    out_incd <- out_incd[agetab, on = c("agegrp" = "agegroup")]

    out_incd[, dimd := factor(dimd,
                              levels = c("1 most deprived", "2", "3", "4", "5",
                                         "6", "7", "8", "9", "10 least deprived"))]

    #Overall
    future({
    strata_mc <- c("year", "mc")
    outstrata <- c("year")
    out_pth_incd <- paste0(out_pth, "incd/", paste(outstrata, collapse = "_"), "/")
    if (!dir.exists(out_pth_incd)) dir.create(out_pth_incd, recursive = TRUE)

    setkey(out_incd, mc, year)
    for (i in dsnm) {

      subtab1 <- prep_incd(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000,
                          col = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ",
          i,
          " per 100,000 persons\n(age-sex standardised)"
        )) +
        expand_limits(y = 0)
      ggsave(filename = paste0(out_pth_incd, "incd_", i, ".png"),
             scale = 3)
    }
})

    #By sex
    future({
    strata_mc <- c("year", "mc", "sex")
    outstrata <- c("year", "sex")
    out_pth_incd <- paste0(out_pth, "incd/", paste(outstrata, collapse = "_"), "/")
    if (!dir.exists(out_pth_incd)) dir.create(out_pth_incd, recursive = TRUE)

    setkey(out_incd, mc, year, sex)
    for (i in dsnm) {

      subtab1 <- prep_incd(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000,
                          col = sex, linetype = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ",
          i,
          " per 100,000 persons\n(age standardised)"
        )) +
        expand_limits(y = 0)
      ggsave(filename = paste0(out_pth_incd, "incd_sex_", i, ".png"),
             scale = 3)
    }
})
    #By imd
    future({
    strata_mc <- c("year", "mc", "dimd")
    outstrata <- c("year", "dimd")
    out_pth_incd <- paste0(out_pth, "incd/", paste(outstrata, collapse = "_"), "/")
    if (!dir.exists(out_pth_incd)) dir.create(out_pth_incd, recursive = TRUE)

    setkey(out_incd, mc, year, dimd)
    for (i in dsnm) {

      subtab1 <- prep_incd(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000,
                          col = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        facet_wrap(~ dimd, scales = "free") +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ", i, " per 100,000 persons\n(age-sex standardised)"
        )) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 6),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))
      ggsave(filename = paste0(out_pth_incd, "incd_dimd_", i, ".png"),
             scale = 3)
    }
    })

    #By agegroup
    future({
    strata_mc <- c("year", "mc", "agegroup10")
    outstrata <- c("year", "agegroup10")
    out_pth_incd <- paste0(out_pth, "incd/", paste(outstrata, collapse = "_"), "/")
    if (!dir.exists(out_pth_incd)) dir.create(out_pth_incd, recursive = TRUE)

    setkey(out_incd, mc, year, agegroup10)
    for (i in dsnm) {

      subtab1 <- prep_incd(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000,
                          col = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        facet_wrap(~ agegroup10, scales = "free") +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ", i, " per 100,000 persons\n(sex standardised)"
        )) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 6),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))
      ggsave(filename = paste0(out_pth_incd, "incd_agegroup_", i, ".png"),
             scale = 3)
    }
})
    #By sex & imd
    future({
    strata_mc <- c("year", "mc", "sex", "dimd")
    outstrata <- c("year",  "sex", "dimd")
    out_pth_incd <- paste0(out_pth, "incd/", paste(outstrata, collapse = "_"), "/")
    if (!dir.exists(out_pth_incd)) dir.create(out_pth_incd, recursive = TRUE)

    setkey(out_incd, mc, year, sex, dimd)
    for (i in dsnm) {

      subtab1 <- prep_incd(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000,
                          col = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        facet_grid(dimd ~ sex, scales = "free") +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ", i, " per 100,000 persons\n(age standardised)"
        )) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 8),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))

      ggsave(filename = paste0(out_pth_incd, "incd_sex_dimd_", i, ".png"),
             scale = 3)
    }
})
    #By sex & agegroup
    future({
    strata_mc <- c("year", "mc", "sex", "agegroup10")
    outstrata <- c("year",  "sex", "agegroup10")
    out_pth_incd <- paste0(out_pth, "incd/", paste(outstrata, collapse = "_"), "/")
    if (!dir.exists(out_pth_incd)) dir.create(out_pth_incd, recursive = TRUE)

    setkey(out_incd, mc, year, sex, agegroup10)
    for (i in dsnm) {

      subtab1 <- prep_incd(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000,
                          col = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        facet_grid(agegroup10 ~ sex, scales = "free") +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ", i, " per 100,000 persons"
        )) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 8),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))

      ggsave(filename = paste0(out_pth_incd, "incd_sex_agegroup_", i, ".png"),
             scale = 3)
    }
})
    #By dimd & agegroup
    future({
    strata_mc <- c("year", "mc", "dimd", "agegroup10")
    outstrata <- c("year",  "dimd", "agegroup10")
    out_pth_incd <- paste0(out_pth, "incd/", paste(outstrata, collapse = "_"), "/")
    if (!dir.exists(out_pth_incd)) dir.create(out_pth_incd, recursive = TRUE)

    setkey(out_incd, mc, year, dimd, agegroup10)
    for (i in dsnm) {

      subtab1 <- prep_incd(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000,
                          col = Type)) +
        #geom_point(size = 0.5) +
        geom_line(size = 0.5) +
        facet_grid(dimd ~ agegroup10 , scales = "free") +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ", i, " per 100,000 persons\n(sex standardised)"
        )) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 8),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))

      ggsave(filename = paste0(out_pth_incd, "incd_dimd_agegroup_", i, ".png"),
             scale = 3)
    }
})
    #By sex, dimd & agegroup
    future({
    strata_mc <- c("year", "mc", "sex", "dimd", "agegroup10")
    outstrata <- c("year", "sex",  "dimd", "agegroup10")
    out_pth_incd <- paste0(out_pth, "incd/", paste(outstrata, collapse = "_"), "/")
    if (!dir.exists(out_pth_incd)) dir.create(out_pth_incd, recursive = TRUE)

    setkey(out_incd, mc, year, sex, dimd, agegroup10)
    for (i in dsnm) {

      subtab1 <- prep_incd(i, strata_mc, outstrata)

      if(i != "breast_ca"){
      ggplot(subtab1[sex == "men"], aes(x = year, y = outcome * 100000,
                                        col = Type)) +
        #geom_point(size = 0.5) +
        geom_line(size = 0.5) +
        facet_grid(dimd ~ agegroup10 , scales = "free") +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ", i, " in men\nper 100,000 persons"
        )) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 8),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))

      ggsave(filename = paste0(out_pth_incd, "incd_dimd_agegroup_men_", i, ".png"),
             scale = 3)
      }

      if(i != "prostate_ca"){
      ggplot(subtab1[sex == "women"], aes(x = year, y = outcome * 100000,
                          col = Type)) +
        #geom_point(size = 0.5) +
        geom_line(size = 0.5) +
        facet_grid(dimd ~ agegroup10 , scales = "free") +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Incidence per 100,000 persons") +
        ggtitle(paste0(
          "Incidence of ", i, " in women\nper 100,000 persons"
        )) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 8),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))

      ggsave(filename = paste0(out_pth_incd, "incd_dimd_agegroup_women_", i, ".png"),
             scale = 3)
      }
    }
    })
    rm(subtab1)
#  })
}


# Overall mrtl ----
if (mtlt_plot) {
    out_pth_all_mtlt <- paste0(out_pth, "all_mtlt/")
    if (!dir.exists(out_pth_all_mtlt)) dir.create(out_pth_all_mtlt, recursive = TRUE)

    out_mrtl <-
      fread(paste0(mdl_rslts_pth,"mrtl_esp.csv.gz"))[scenario == "sc0"]
    out_mrtl[, dimd := factor(dimd,
                              levels = c("1 most deprived", "2", "3", "4", "5",
                                         "6", "7", "8", "9", "10 least deprived"))]
    out_mrtl <- out_mrtl[agetab, on = c("agegrp" = "agegroup")]

    setkey(out_mrtl, mc, year)

    # onsmrtl <- fread("./inputs/mortality/mortality_dimd_2001-19.csv")[age > 29]
    # la <- data.table(age = 30:90, agegroup5 = agegrp_name(min_age = 30, 90, grp_width = 5,
    #                                                        match_input = TRUE, match_input_max_age = 90))
    # absorb_dt( onsmrtl, la)
    # onsmrtl <-  onsmrtl[, .(deaths = sum(deaths), pop_size = sum(pop_size)), by = .(year, agegroup5, sex, dimd)]
    # onsmrtl[esp90, on = c("agegroup5" = "agegroup", "sex", "dimd"), wt_esp := i.wt_esp]
    # onsmrtl[, deaths_esp := wt_esp * deaths/pop_size]
    # onsmrtl[agetab, on = c("agegroup5" = "agegroup"), agegroup10 := i.agegroup10]
    # onsmrtl[agegroup5 == "90+", agegroup10 := "90-99"]
    # onsmrtl <-  onsmrtl[, .(deaths = sum(deaths_esp), popsize = sum(wt_esp)), by = .(year, agegroup10, sex, dimd)]

    onsmrtl <- fread("./inputs/mortality/lt.csv")[between(year, 2002, 2043)]
    onsmrtl[esp, on = c("agegrp" = "agegroup", "sex", "dimd"), wt_esp := i.wt_esp]
    onsmrtl[, deaths_esp := wt_esp * mx_total]
    onsmrtl[agetab, on = c("agegrp" = "agegroup"), agegroup10 := i.agegroup10]
    onsmrtl <-  onsmrtl[!is.na(agegroup10), .(deaths = sum(deaths_esp), popsize = sum(wt_esp)), by = .(year, agegroup10, sex, dimd)]

    onsmrtl[, dimd := factor(dimd,
                             levels = c("1 most deprived", "2", "3", "4", "5",
                                        "6", "7", "8", "9", "10 least deprived"))]




    #Overall
    strata_mc <- c("year", "mc")
    outstrata <- c("year")
    setkey(out_mrtl, mc, year)

    subtab1 <- prep_mrtl(strata_mc, outstrata )


    ggplot(subtab1, aes(
      x = year,
      y = all_ftlt * 100000,
      col = Type
    )) +
      geom_point(size = 0.5) +
      geom_line() +
      geom_smooth(linetype = "dotdash", se = FALSE) +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "All-cause mortality (per 100,000 persons, age-sex standardised)",
                         limits = c(0, 2000)) +
      ggtitle("All-cause mortality (per 100,000 persons)") +
      scale_colour_brewer(type = "qual")
    ggsave(filename = paste0(out_pth_all_mtlt, "all_cause_mrtl.png"),
           scale = 3)

    rm(subtab1)





    #By imd

    strata_mc <- c("year", "mc", "dimd")
    outstrata <- c("year", "dimd")
    setkey(out_mrtl, mc, year, dimd)

    subtab1 <- prep_mrtl(strata_mc, outstrata )

    ggplot(subtab1[year != 2013], aes(
      x = year,
      y = all_ftlt * 100000,
      col = Type,
    )) +
      geom_point(size = 0.5) +
      geom_line() +
      geom_smooth( se = FALSE) +
      facet_wrap(~dimd, scales = "free") +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "All-cause mortality (per 100,000 persons, age-sex standardised)") +
      scale_colour_brewer(type = "qual") +
      ggtitle("All-cause mortality (per 100,000 persons)") +
      theme(strip.text = element_text(size = 6),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.text.y = element_text(size = 8))
    ggsave(filename = paste0(out_pth_all_mtlt, "all_cause_mrtl_dimd.png"),
           scale = 3)


    #By imd & 10-year age-group

    strata_mc <- c("year", "mc", "dimd", "agegroup10")
    outstrata <- c("year", "dimd", "agegroup10")
    setkey(out_mrtl, mc, year, dimd, agegroup10)

    subtab1 <- prep_mrtl(strata_mc, outstrata )

    ggplot(subtab1[year != 2013], aes(
      x = year,
      y = all_ftlt * 100000,
      col =  Type
    )) +
      #geom_point(size = 0.5) +
      geom_line(size = 0.5) +
      geom_smooth(linetype = "dotdash", se = FALSE) +
      facet_grid(agegroup10~ dimd, scales = "free") +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "All-cause mortality (per 100,000 persons)") +
      ggtitle("All-cause mortality\n(per 100,000 persons, sex-standardised)") +
      scale_colour_brewer(type = "qual") +
      labs(col = "Type") +
      theme(strip.text = element_text(size = 8),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.y = element_text(size = 10))
    ggsave(filename = paste0(out_pth_all_mtlt, "all_cause_mrtl_dimd_age10.png"),
           scale = 3)

    #By sex & 10-year age-group
    strata_mc <- c("year", "mc", "sex", "agegroup10")
    outstrata <- c("year", "sex", "agegroup10")
    setkey(out_mrtl, mc, year, sex, agegroup10)

    subtab1 <- prep_mrtl(strata_mc, outstrata )


    ggplot(subtab1[year != 2013], aes(
      x = year,
      y = all_ftlt * 100000,
      col = Type,
    )) +
      geom_point(size = 0.5) +
      geom_line() +
      geom_smooth(linetype = "dotdash", se = FALSE) +
      facet_grid(agegroup10~ sex, scales = "free") +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "All-cause mortality (per 100,000 persons)") +
      ggtitle("All-cause mortality (per 100,000 persons)") +
      scale_colour_brewer(type = "qual") +
      labs(col = "Type")+
      theme(strip.text = element_text(size = 8),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.y = element_text(size = 10))
    ggsave(filename = paste0(out_pth_all_mtlt, "all_cause_mrtl_sex_age10.png"),
           scale = 3)


    #By sex & dimd
    strata_mc <- c("year", "mc", "sex", "dimd")
    outstrata <- c("year", "sex", "dimd")
    setkey(out_mrtl, mc, year, sex, dimd)

    subtab1 <- prep_mrtl(strata_mc, outstrata )



     ggplot(subtab1[year != 2013], aes(
      x = year,
      y = all_ftlt * 100000,
      col =  Type
    )) +
      geom_point(size = 0.5) +
      geom_line() +
      geom_smooth(linetype = "dotdash", se = FALSE) +
      facet_grid(dimd ~ sex, scales = "free") +
      scale_x_continuous(name = "Year") +
      scale_y_continuous(name = "All-cause mortality (per 100,000 persons)") +
       ggtitle("All-cause mortality (per 100,000 persons, age standardised)") +
      labs(col = "Type") +
       scale_colour_brewer(type = "qual") +
       theme(strip.text = element_text(size = 6),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.text.y = element_text(size = 10))

    ggsave(filename = paste0(out_pth_all_mtlt, "all_cause_mrtl_sex_dimd.png"),
           scale = 3)


  #By sex & dimd & agegroup10
  strata_mc <- c("year", "mc", "sex", "dimd", "agegroup10")
  outstrata <- c("year", "sex", "dimd", "agegroup10")
  setkey(out_mrtl, mc, year, sex, dimd, agegroup10)

  subtab1 <- prep_mrtl(strata_mc, outstrata )


  #For men
  ggplot(subtab1[year != 2013 & sex == "men" ], aes(
    x = year,
    y = all_ftlt * 100000,
    col =  Type
  )) +
    geom_line() +
    geom_smooth(linetype = "dotdash", se = FALSE) +
    facet_grid(agegroup10 ~ dimd, scales = "free") +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "All-cause mortality in men\n(per 100,000 persons)") +
    ggtitle("All-cause mortality (per 100,000 persons)") +
    labs(col = "Type") +
    scale_colour_brewer(type = "qual") +
    theme(strip.text = element_text(size = 6),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 8))

  ggsave(filename = paste0(out_pth_all_mtlt, "all_cause_mrtl_agegrp_dimd_men.png"),
         scale = 3)

  #For women
  ggplot(subtab1[year != 2013 & sex == "women" ], aes(
    x = year,
    y = all_ftlt * 100000,
    col =  Type
  )) +
    geom_line() +
    geom_smooth(linetype = "dotdash", se = FALSE) +
    facet_grid(agegroup10 ~ dimd, scales = "free") +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = "All-cause mortality in women\n(per 100,000 persons)") +
    ggtitle("All-cause mortality (per 100,000 persons)") +
    labs(col = "Type") +
    scale_colour_brewer(type = "qual") +
    theme(strip.text = element_text(size = 6),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.text.y = element_text(size = 8))

  ggsave(filename = paste0(out_pth_all_mtlt, "all_cause_mrtl_agegrp_dimd_women.png"),
         scale = 3)
#})


}


# Cause-specific mrtl ----
if (ftlt_plot) {
 # future({
    out_pth_ftlt <- paste0(out_pth, "ftlt/")
    if (!dir.exists(out_pth_ftlt)) dir.create(out_pth_ftlt, recursive = TRUE)

    tt <- grep("_ftlt_N$", names(cprdtab), value = TRUE)
    # Ensure all_ftlt is always the first
    tt <- tt[order(match(tt, "all_ftlt_N"))]
    cprdtab[, nonmodelled_ftlt_N := Reduce(`-`, .SD), .SDcols = tt]

    out_ftlt <-
      fread(paste0(mdl_rslts_pth,"/dis_mrtl_esp.csv.gz"))[scenario == "sc0"]

    out_ftlt[, dimd := factor(dimd,
                             levels = c("1 most deprived", "2", "3", "4", "5",
                                        "6", "7", "8", "9", "10 least deprived"))]

    setkey(out_ftlt, mc, year)
    out_ftlt <- out_ftlt[agetab, on = c("agegrp" = "agegroup")]

    dsnm_cf <- names(out_ftlt)[!names(out_ftlt) %in% c("agegrp",
                                                     "mc",
                                                     "scenario",
                                                     "year",
                                                     "sex",
                                                     "dimd",
                                                     "ethnicity",
                                                     "sha",
                                                     "popsize",
                                                     "agegroup10")]

    #Overall
    strata_mc <- c("year", "mc")
    outstrata <- c("year")
    setkey(out_ftlt, mc, year)
    for (i in dsnm_cf) {

      subtab1 <- prep_cftlt(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000, col = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0(
            "Cause-specific mortality of ",
            i,
            " per 100,000 persons\n(age-sex standardised)"
          )
        ) +
        expand_limits(y = 0)
      ggsave(filename = paste0(out_pth_ftlt, "ftlt_", i, ".png"),
             scale = 3)
    }


    #by sex
    strata_mc <- c("year", "mc", "sex")
    outstrata <- c("year", "sex")
    setkey(out_ftlt, mc, year, sex)
    for (i in dsnm_cf) {

      subtab1 <- prep_cftlt(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000, col = sex, linetype = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0(
            "Cause-specific mortality of ",
            i,
            " per 100,000 persons\n(age standardised)"
          )
        ) +
        expand_limits(y = 0)
      ggsave(filename = paste0(out_pth_ftlt, "ftlt_sex_", i, ".png"),
             scale = 3)
    }

    #by imd
    strata_mc <- c("year", "mc", "dimd")
    outstrata <- c("year", "dimd")
    setkey(out_ftlt, mc, year, dimd)
    for (i in dsnm_cf) {

      subtab1 <- prep_cftlt(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000, col = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        facet_wrap(~ dimd, scales = "free") +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0(
            "Cause-specific mortality of ",
            i,
            " per 100,000 persons\n(age-sex standardised)"
          )
        ) +
        expand_limits(y = 0)  +
        theme(strip.text = element_text(size = 6),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))
      ggsave(filename = paste0(out_pth_ftlt, "ftlt_dimd_", i, ".png"),
             scale = 3)
    }

    #by agegroup
    strata_mc <- c("year", "mc", "agegroup10")
    outstrata <- c("year", "agegroup10")
    setkey(out_ftlt, mc, year, agegroup10)
    for (i in dsnm_cf) {

      subtab1 <- prep_cftlt(i, strata_mc, outstrata)

      ggplot(subtab1, aes(x = year, y = outcome * 100000, col = Type)) +
        geom_point(size = 0.5) +
        geom_line() +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        facet_wrap(~ agegroup10, scales = "free") +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0(
            "Cause-specific mortality of ",
            i,
            " per 100,000 persons\n(sex standardised)"
          )
        ) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 6),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))
      ggsave(filename = paste0(out_pth_ftlt, "ftlt_age_", i, ".png"),
             scale = 3)
    }



    #By sex & dimd
    strata_mc <- c("year", "mc", "sex", "dimd")
    outstrata <- c("year", "sex", "dimd")
    setkey(out_ftlt, mc, year, sex, dimd)
    for (i in dsnm_cf) {

      subtab1 <- prep_cftlt(i, strata_mc, outstrata)

      ggplot(subtab1, aes(
        x = year,
        y = outcome * 100000,
        col = Type
      )) +
      #  geom_point(size = 0.5) +
        geom_line(size = 0.5) +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        facet_grid( dimd ~ sex, scales = "free") +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0("Cause-specific mortality of ", i,
                 " per 100,000 persons\n(age standardised)")) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 6),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))
      ggsave(filename = paste0(out_pth_ftlt, "ftlt_sex_dimd_", i, ".png"),
             scale = 3)
    }


    #By sex & agegroup10
    strata_mc <- c("year", "mc", "sex", "agegroup10")
    outstrata <- c("year", "sex", "agegroup10")
    setkey(out_ftlt, mc, year, sex, agegroup10)
    for (i in dsnm_cf) {

      subtab1 <- prep_cftlt(i, strata_mc, outstrata)

      ggplot(subtab1, aes(
        x = year,
        y = outcome * 100000,
        col = Type
      )) +
        #  geom_point(size = 0.5) +
        geom_line(size = 0.5) +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        facet_grid(agegroup10 ~ sex, scales = "free") +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0("Cause-specific mortality of ", i,
                 " per 100,000 persons)")) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 6),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))
      ggsave(filename = paste0(out_pth_ftlt, "ftlt_sex_age_", i, ".png"),
             scale = 3)
    }

    #By dimd & agegroup10
    strata_mc <- c("year", "mc", "dimd", "agegroup10")
    outstrata <- c("year", "dimd", "agegroup10")
    setkey(out_ftlt, mc, year, dimd, agegroup10)
    for (i in dsnm_cf) {

      subtab1 <- prep_cftlt(i, strata_mc, outstrata)

      ggplot(subtab1, aes(
        x = year,
        y = outcome * 100000,
        col = Type
      )) +
        #  geom_point(size = 0.5) +
        geom_line(size = 0.5) +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        facet_grid(dimd ~ agegroup10, scales = "free") +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0("Cause-specific mortality of ", i,
                 " per 100,000 persons\n(Sex standardised)")) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 6),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))
      ggsave(filename = paste0(out_pth_ftlt, "ftlt_dimd_age_", i, ".png"),
             scale = 3)
    }

    #By sex, dimd & agegroup10
    strata_mc <- c("year", "mc", "sex", "dimd", "agegroup10")
    outstrata <- c("year", "sex","dimd", "agegroup10")
    setkey(out_ftlt, mc, year, sex, dimd, agegroup10)
    for (i in dsnm_cf) {

      subtab1 <- prep_cftlt(i, strata_mc, outstrata)

      if(i != "breast_ca"){
      ggplot(subtab1[sex == "men"], aes(
        x = year,
        y = outcome * 100000,
        col = Type
      )) +
        #  geom_point(size = 0.5) +
        geom_line(size = 0.5) +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        facet_grid(dimd ~ agegroup10, scales = "free") +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0("Cause-specific mortality of ", i,
                 " in men\nper 100,000 persons")) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 6),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))
      ggsave(filename = paste0(out_pth_ftlt, "ftlt_dimd_age_men_", i, ".png"),
             scale = 3)
      }

      if(i != "prostate_ca"){
      ggplot(subtab1[sex == "women"], aes(
        x = year,
        y = outcome * 100000,
        col = Type
      )) +
        #  geom_point(size = 0.5) +
        geom_line(size = 0.5) +
        geom_smooth(linetype = "dotdash", se = FALSE) +
        facet_grid(dimd ~ agegroup10, scales = "free") +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = "Mortality per 100,000 persons") +
        ggtitle(
          paste0("Cause-specific mortality of ", i,
                 " in women\nper 100,000 persons")) +
        expand_limits(y = 0) +
        theme(strip.text = element_text(size = 6),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              axis.text.y = element_text(size = 8))
      ggsave(filename = paste0(out_pth_ftlt, "ftlt_dimd_age_women_", i, ".png"),
             scale = 3)
      }
    }

    rm(subtab1)
 # })
}
