# CPRD Validation charts

library(data.table) # for fast data manipulation
library(fst)
library(ggplot2)
library(ggthemes)
library(scales)
library(RColorBrewer)
library(viridis)
library(CKutils)
library(doParallel)
library(foreach)

if (exists("IMPACTncd")) {
  workers = IMPACTncd$design$sim_prm$clusternumber
  dsnm <- names(IMPACTncd$design$sim_prm$diseases)
  mdl_rslts_pth <- paste0(IMPACTncd$design$sim_prm$output_dir, "/tables/")
  out_pth <- paste0(IMPACTncd$design$sim_prm$output_dir, "/plots/validation/")

} else if (exists("design")) {
  workers = design$sim_prm$clusternumber
  dsnm <- names(design$sim_prm$diseases)
  mdl_rslts_pth <- paste0(design$sim_prm$output_dir, "/tables/")
  out_pth <- paste0(design$sim_prm$output_dir, "/plots/validation/")

} else if (file.exists("./inputs/sim_design.yaml")) {
  design <- IMPACTncdEngl::Design$new("./inputs/sim_design.yaml")
  workers = design$sim_prm$clusternumber
  dsnm <- names(design$sim_prm$diseases)
  mdl_rslts_pth <- paste0(design$sim_prm$output_dir, "/tables/")
  out_pth <- paste0(design$sim_prm$output_dir, "/plots/validation/")

} else stop()

dsnm <- dsnm[!dsnm %in% c("dm", "ctdra", "cancer", "cms1st_cont", "nonmodelled")]
if (!dir.exists(out_pth)) dir.create(out_pth, recursive = TRUE)
# mdl_rslts_pth <- "past_results/Model_results_13-01-23_new_alcpr_def/tables/"
# out_pth <- "past_results/Model_results_13-01-23_new_alcpr_def//plots/validation/"

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

# esp90 <- copy(esp)
# esp90[agegroup %in% c("90-94", "95-99"), agegroup := "90+"]
# esp90 <- esp90[, .(wt_esp = sum(wt_esp)), by = .(agegroup, sex, dimd)]

# dl_plot <- c("Anxiety & Depression", "Asthma", "Atrial Fibrillation", "CHD",
#              "COPD" , "Chronic Kidney Disease", "Dementia", "Heart failure",
#              "Hypertension" ,  "Obesity", "Other cancers" ,
#              "Breast cancer", "Colorectal cancer",
#              "Lung cancer", "Prostate cancer",
#              "Stroke", "Type 2 Diabetes", "CMS > 0", "CMS > 1", "CMS > 1.5",
#              "CMS > 2")
#
# dl_lookup <- data.table(sim = dl_cprd_new, mortality = dl_mort, plot = dl_plot)


#strata <- c( "year", "agegroup" ,"sex", "dimd" )

cprdtab <- fread("./secure_data/summarytable.csv")[!agegroup %in% c("20-24", "25-29")][, wt_esp := NULL]
absorb_dt(cprdtab, esp)
# cprdtab[, wt_esp := wt_esp * unique(wt_esp) / (wt_esp * popsize),
#         keyby = .(year, agegroup, sex, dimd)]
# cprdtab[, popsize_wtd := popsize * wt_esp]
cprdtab[, popsize_wtd := wt_esp]
cprdtab[, wt_esp := wt_esp / popsize]

ll <- grep("^alcpr", names(cprdtab),  value = TRUE)
grep("^alcpr_rev", ll, invert = TRUE, value = TRUE)
cprdtab[, (grep("^alcpr_rev", ll, invert = TRUE, value = TRUE)) := NULL]
setnames(cprdtab, grep("^alcpr_rev", ll, value = TRUE), grep("^alcpr_rev", ll, invert = TRUE, value = TRUE))
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
tt <- grep("_ftlt_N$", names(cprdtab), value = TRUE)
# Ensure all_ftlt is always the first
tt <- tt[order(match(tt, "all_ftlt_N"))]
cprdtab[, nonmodelled_ftlt_N := Reduce(`-`, .SD), .SDcols = tt]

# outstrata <- "year"
# suffix <- "incd"
validation_plot <- function(outstrata, suffix) {
  if (suffix == "prvl") {
    term0 <- "Prevalence"
    term1 <- "prevalence by "
  } else if (suffix == "incd") {
    term0 <- "Incidence"
    term1 <- "incidence by "
  } else if (suffix == "mrtl") {
    term0 <- "Mortality"
    term1 <- "mortality by "
  } else if (suffix == "ftlt") {
    term0 <- "Fatality"
    term1 <- "fatality by "
    dsnm <- c(dsnm, "nonmodelled")
  } else stop()

  term2 <- paste(outstrata, collapse = "_")
  if (length(outstrata) == 1L && outstrata == "year") {
    term3 <- " (age-sex-dimd standardised).csv"
    term4 <- "(age-sex-dimd standardised)"
  }
  Y    <- paste0(suffix, "_rate_50.0%")
  Ymin <- paste0(suffix, "_rate_2.5%")
  Ymax <- paste0(suffix, "_rate_97.5%")

  tt <- fread(paste0(mdl_rslts_pth, term1, term2, term3))[scenario == "sc0"]
  tt[, `:=` (Type = "Modelled", scenario = NULL)]
  tt[, c(paste0(suffix, "_rate_10.0%"), paste0(suffix, "_rate_90.0%")) := NULL]

  out_pth_plot <- paste0(out_pth, paste0(suffix, "/"), term2, "/")
  if (!dir.exists(out_pth_plot)) dir.create(out_pth_plot, recursive = TRUE)

  # Plots for
  if (suffix == "mrtl") {
    onsmrtl <- fread("./inputs/mortality/lt.csv")[between(year, 2002, 2043)]
    onsmrtl[, dimd := factor(dimd,
                             levels = c("1 most deprived", "2", "3", "4", "5",
                                        "6", "7", "8", "9", "10 least deprived"))]
    onsmrtl[esp, on = c("agegrp" = "agegroup", "sex", "dimd"), wt_esp := i.wt_esp]
    onsmrtl[, deaths_esp := wt_esp * mx_total]
    onsmrtl[agetab, on = c("agegrp" = "agegroup"), agegrp := i.agegroup10]
    onsmrtl <-  onsmrtl[agegrp %in% agetab$agegroup10,
                        .(mx = sum(deaths_esp)/sum(wt_esp)),
                        by = outstrata]
    onsmrtl[, c(paste0(suffix, "_rate_2.5%"), paste0(suffix, "_rate_97.5%"),
                "Type") := .(
                  mx, mx, "Observed/Forecasted")]
    setnames(onsmrtl, "mx", paste0(suffix, "_rate_50.0%"))

    t1 <- rbind(tt, onsmrtl)

    # TODO avoid replication of plotting code
      ggplot(t1, aes(
        x = year,
        y = .data[[Y]],
        ymin = .data[[Ymin]],
        ymax = .data[[Ymax]],
        col = Type, group = Type
      )) +
        geom_point(size = 0.5) +
        geom_ribbon(fill = NA) +
        geom_line(linetype = "dotdash") +
        scale_x_continuous(name = "Year") +
        scale_y_continuous(name = term0, labels = scales::percent) +
        ggtitle(paste0("All-cause mortality",  "\n", term4)) +
        expand_limits(y = 0) #+
      # scale_colour_brewer(type = "qual")


    ggsave(
      filename = paste0(out_pth_plot, "All-cause_", suffix, ".png"),
      scale = 1.5
    )
  } else {

    if (Sys.info()["sysname"] == "Windows") {
      cl <-
        makeCluster(workers) # used for clustering. Windows compatible
      registerDoParallel(cl)
    } else {
      registerDoParallel(workers) # used for forking. Only Linux/OSX compatible
    }

    out <- foreach(
      i = dsnm,
      .inorder = FALSE,
      .options.multicore = list(preschedule = FALSE),
      .verbose = TRUE, # self$design$sim_prm$logs,
      .packages = c(
        "R6",
        "CKutils",
        "IMPACTncdEngl",
        "data.table",
        "ggplot2"
      ),
      .export = NULL,
      .noexport = NULL # c("time_mark")
    ) %dopar% {
      dp <- paste0(i, "_", suffix)
      dp2 <- paste0(i, "_", suffix, "_N") # NOTE need to use e suffix!!!


      t1 <- tt[disease == dp, .SD,
               .SDcols = c(outstrata, "Type", Y, Ymin, Ymax)]
      if (suffix == "ftlt" && i != "nonmodelled") {
        t2 <- cprdtab[, .("____y" = sum(get(dp2) * wt_esp) / sum(get(paste0(i, "_prvl_N")) * wt_esp)), keyby = outstrata]
      } else {
        t2 <- cprdtab[, .("____y" = sum(get(dp2) * wt_esp) / sum(popsize_wtd)), keyby = outstrata]
      }
      t2[, c(Ymin, Ymax, "Type") := .(`____y`, `____y`, "Observed")]
      setnames(t2, "____y", Y)
      t1 <- rbind(t1, t2)

      ggplot(t1, aes(
          x = year,
          y = .data[[Y]],
          ymin = .data[[Ymin]],
          ymax = .data[[Ymax]],
          col = Type, group = Type
        )) +
          geom_point(size = 0.5) +
          geom_ribbon(fill = NA) +
          geom_line(linetype = "dotdash") +
          scale_x_continuous(name = "Year") +
          scale_y_continuous(name = term0, labels = scales::percent) +
          ggtitle(paste0(term0, " of ", i, "\n", term4)) +
          expand_limits(y = 0) #+
        # scale_colour_brewer(type = "qual")

      ggsave(
        filename = paste0(out_pth_plot, suffix, "_", i, ".png"),
        scale = 1.5
      )

      NULL
    }
  }
  if (exists("cl")) stopCluster(cl)
}

validation_plot("year", "prvl")
validation_plot("year", "incd")
validation_plot("year", "ftlt")
validation_plot("year", "mrtl")

