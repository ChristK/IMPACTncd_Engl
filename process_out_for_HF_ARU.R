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

tbl_smmrs <- function(
    what = c(
      "prvl", "prvl_change", "incd", "incd_change",
      "ftlt", "ftlt_change", "mrtl", "mrtl_change",
      "cms_score", "cms_score_change", "cms_score_age",
      "cms_score_age_change", "cms_count", "cms_count_change",
      "pop"
    ),
    type = c("ons", "esp"),
    strata,
    output_dir = output_dir,
    prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
    baseline_year = 2019L, # only used for prvl_change etc.
    two_agegrps = FALSE # if TRUE, agegrp is 30-69 and 70-99
) {
  strata <- lapply(strata, function(x) {
    c("mc", "scenario", x)
  })

  # construct file path to read from summaries
  str0 <- c(
    "prvl" = "prvl",
    "prvl_change" = "prvl",
    "incd" = "incd",
    "incd_change" = "incd",
    "ftlt" = "/dis_mrtl",
    "ftlt_change" = "/dis_mrtl",
    "mrtl" = "mrtl",
    "mrtl_change" = "mrtl",
    "cms_score" = "cms_score",
    "cms_score_change" = "cms_score",
    "cms_score_age" = "cms_score_by_age",
    "cms_score_age_change" = "cms_score_by_age",
    "cms_count" = "cms_count",
    "cms_count_change" = "cms_count",
    "pop" = "prvl"
  )
  str1 <- c("ons" = "_scaled_up.csv.gz", "esp" = "_esp.csv.gz")
  fpth <- file.path(output_dir, "summaries", paste0(str0[[what]], str1[[type]]))
  if (!file.exists(fpth)) {
    message(fpth, " doesn't exist")
    return(NULL)
  }
  # other useful strings
  str2 <- c(
    "prvl" = "_prvl$|^popsize$",
    "prvl_change" = "_prvl$|^popsize$",
    "incd" = "_incd$|^popsize$",
    "incd_change" = "_incd$|^popsize$",
    "ftlt" = "_deaths$|_prvl$",
    "ftlt_change" = "_deaths$|_prvl$",
    "mrtl" = "_mrtl$|^popsize$",
    "mrtl_change" = "_mrtl$|^popsize$",
    "cms_score" = "cms_score",
    "cms_score_change" = "cms_score",
    "cms_score_age" = "cms_score",
    "cms_score_age_change" = "cms_score",
    "cms_count" = "cms_count",
    "cms_count_change" = "cms_count",
    "pop" = "^popsize$"
  ) # used in grep
  str3 <- c(
    "prvl" = "prvl_rate_",
    "prvl_change" = "prct_change_",
    "incd" = "incd_rate_",
    "incd_change" = "prct_change_",
    "ftlt" = "ftlt_rate_",
    "ftlt_change" = "ftlt_rate_",
    "mrtl" = "mrtl_rate_",
    "mrtl_change" = "mrtl_change_",
    "cms_score" = "mean_cms_score_",
    "cms_score_change" = "mean_cms_score_",
    "cms_score_age" = "mean_cms_score_",
    "cms_score_age_change" = "mean_cms_score_",
    "cms_count" = "mean_cms_count_",
    "cms_count_change" = "mean_cms_count_",
    "pop" = "pop_size_"
  ) # used to col name output
  str4 <- c(
    "prvl" = "prevalence by ",
    "prvl_change" = "prevalence change by ",
    "incd" = "incidence by ",
    "incd_change" = "incidence change by ",
    "ftlt" = "fatality by ",
    "ftlt_change" = "fatality change by ",
    "mrtl" = "mortality by ",
    "mrtl_change" = "mortality change by ",
    "cms_score" = "mean CMS score by ",
    "cms_score_change" = "mean CMS score change by ",
    "cms_score_age" = "mean CMS score by ",
    "cms_score_age_change" = "mean CMS score change by ",
    "cms_count" = "mean CMS count by ",
    "cms_count_change" = "mean CMS count change by ",
    "pop" = "pop size by "
  )



  tt <-
    fread(fpth)[, `:=`(
      year = year + 2000L,
      dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived"))
    )]

  if (two_agegrps) {
    sTablesSubDirPath <- file.path(simulationParameters$output_dir, "tables2agegrps/")
    if (!dir.exists(sTablesSubDirPath)) dir.create(sTablesSubDirPath)
    if ("agegrp" %in% names(tt)) {
      tt[agegrp %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69"), agegrp := "30-69"]
      tt[agegrp %in% c("70-74", "75-79", "80-84", "85-89", "90-94", "95-99"), agegrp := "70-99"]
    }
  }

  # For ftlt I need prvl for the denominator
  if (grepl("^ftlt", what)) {
    fpth <- file.path(output_dir, "summaries", paste0(str0[["prvl"]], str1[[type]]))
    if (!file.exists(fpth)) stop(fpth, " doesn't exist")

    t1 <- fread(fpth)[, `:=`(
      year = year + 2000L,
      dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived"))
    )]
    setnames(t1, "popsize", "nonmodelled_prvl")
    if (two_agegrps && "agegrp" %in% names(t1)) {
      t1[agegrp %in% c("30-34", "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69"), agegrp := "30-69"]
      t1[agegrp %in% c("70-74", "75-79", "80-84", "85-89", "90-94", "95-99"), agegrp := "70-99"]
    }
    absorb_dt(tt, t1)
    tt <- tt[nonmodelled_prvl > 0] # This is the denom. Cannot be 0 and it is meaningless anyway
  }

  lapply(strata, function(x) {
    if (grepl("^cms_", what)) {
      d <- tt[, .("value" = weighted.mean(get(str2[[what]]), popsize)),
              keyby = eval(x)
      ]

      if (grepl("_change$", what)) { # when calculating change
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year")), value := value / i.value]
      }
      d <- d[, as.list(fquantile(value, prbl)), keyby = eval(setdiff(x, "mc"))]
      setnames(d, c(setdiff(x, "mc"), percent(prbl, prefix = str3[[what]])))
    } else { # if not cms...
      d <- tt[, lapply(.SD, sum),
              .SDcols = patterns(str2[[what]]),
              keyby = x
      ]


      if (grepl("^ftlt", what)) {
        nm <- names(d)
        nm <- grep("_deaths$", nm, value = TRUE)
        nm <- gsub("_deaths$", "", nm)
        nm <- setdiff(nm, "alive")
        for (i in nm) {
          set(
            d, NULL, paste0(i, "_ftlt"),
            d[[paste0(i, "_deaths")]] / d[[paste0(i, "_prvl")]]
          )
        }

        nm <- names(d)
        nm <- grep("_deaths$|_prvl$", nm, value = TRUE)
        d[, (nm) := NULL]
        setnafill(d, "const", 0, cols = grep("_ftlt$", names(d), value = TRUE))
      } else if (what != "pop") { # if not ftlt related and not pop
        d <- d[, lapply(.SD, function(y) {
          y / popsize
        }), keyby = x]
      }

      d <- melt(d, id.vars = x)

      if (grepl("_change$", what)) { # when calculating change
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year"), "variable"), value := value / i.value]
      }

      setkey(d, "variable")
      d <-
        d[, fquantile_byid(value, prbl, id = as.character(variable), rounding = what == "pop"),
          keyby = eval(setdiff(x, "mc"))
        ]
      setnames(d, c(
        setdiff(x, "mc"),
        "disease",
        percent(prbl, prefix = str3[[what]])
      ))
      if (what == "pop") {
        d[, disease := NULL]
      } else {
        d <- d[disease != "popsize"]
      }
    }
    setkeyv(d, setdiff(x, "mc"))
    str5 <- c(
      "ons" = " (not standardised).csv",
      "esp" = paste0(" (", paste(setdiff(c("mc", "scenario", "year", "age", "sex", "dimd"), x),
                                 collapse = "-"
      ), " standardised).csv")
    )

    str6 <- paste0(
      str4[[what]],
      paste(setdiff(x, c("mc", "scenario")), collapse = "-"),
      str5[[type]]
    ) # used for output file name/path

    fwrite(d, file.path(
      output_dir,
      ifelse(two_agegrps, "tables2agegrps", "tables"), str6
    ))
  })
}

xps_summ <- function(outstrata, prbl, path,
                     sTablesSubDirPath, what = "xps",
                     type = c("ons", "esp")) {
  if("ons" %in% type){
    xps_tab <- fread(file.path(path, "xps/xps20.csv.gz"))
    if (all(c("agegrp20", "sex", "qimd") %in% outstrata)) {
      d <- xps_tab[sex != "All" & agegrp20 != "All" & qimd != "All"]
      file_name <- "/exposures by year-agegroup-sex-qimd (not standardised).csv"
    } else {
      if (all(c("qimd", "sex") %in% outstrata)) {
        d <- xps_tab[sex != "All" & agegrp20 == "All" & qimd != "All"]
        file_name <- "/exposures by year-sex-qimd (not standardised).csv"
      } else if (all(c("qimd", "agegrp20") %in% outstrata)) {
        d <- xps_tab[sex == "All" & agegrp20 != "All" & qimd != "All"]
        file_name <- "/exposures by year-agegroup-qimd (not standardised).csv"
      } else if (all(c("agegrp20", "sex") %in% outstrata)) {
        d <- xps_tab[sex != "All" & agegrp20 != "All" & qimd == "All"]
        file_name <- "/exposures by year-agegroup-sex (not standardised).csv"
      } else if ("agegrp20" %in% outstrata & !any(c("sex", "qimd") %in% outstrata)) {
        d <- xps_tab[sex == "All" & agegrp20 != "All" & qimd == "All"]
        file_name <- "/exposures by year-agegroup (not standardised).csv"
      } else if ("sex" %in% outstrata & !any(c("agegrp20", "qimd") %in% outstrata)) {
        d <- xps_tab[sex != "All" & agegrp20 == "All" & qimd == "All"]
        file_name <- "/exposures by year-sex (not standardised).csv"
      } else if ("qimd" %in% outstrata & !any(c("agegrp20", "sex") %in% outstrata)) {
        d <- xps_tab[sex == "All" & agegrp20 == "All" & qimd != "All"]
        file_name <- "/exposures by year-qimd (not standardised).csv"
      } else if(!all(c("sex", "qimd", "agegrp20") %in% outstrata)) {
        d <- xps_tab[sex == "All" & agegrp20 == "All" & qimd == "All"]
        file_name <- "/exposures by year (not standardised).csv"
      }
    }
  }
  if("esp" %in% type) {
    xps_tab <- fread(file.path(simulationParameters$output_dir, "xps/xps_esp.csv.gz"))
    if (all(c("sex", "qimd") %in% outstrata)) {
      d <- xps_tab[sex != "All" & qimd != "All"]
      file_name <- "/exposures by year-sex-qimd (age standardised).csv"
    } else {
      if ("sex" %in% outstrata & !any(c("qimd") %in% outstrata)) {
        d <- xps_tab[sex != "All" & qimd == "All"]
        file_name <- "/exposures by year-sex (age-qimd standardised).csv"
      } else if ("qimd" %in% outstrata & !any(c("sex") %in% outstrata)) {
        d <- xps_tab[sex == "All" & qimd != "All"]
        file_name <- "/exposures by year-qimd (age-sex standardised).csv"
      } else if(!all(c("sex", "qimd") %in% outstrata)) {
        d <- xps_tab[sex == "All" & qimd == "All"]
        file_name <- "/exposures by year (age-sex-qimd standardised).csv"
      }
    }
  }
  xps_names <- grep("_curr_xps$", names(xps_tab), value = TRUE)
  d <- d[, lapply(.SD, mean), .SDcols = patterns("_curr_xps$"), keyby = eval(outstrata)]
  d <- melt(d, id.vars = outstrata)
  setkey(d, "variable")
  d <- d[, fquantile_byid(value, prbl, id = as.character(variable)), keyby = eval(setdiff(outstrata, "mc"))]
  setnames(d, c(setdiff(outstrata, "mc"), "exposure", percent(prbl, prefix = "xps_mean_")))
  setkeyv(d, setdiff(outstrata, "mc"))
  fwrite(d, paste0(sTablesSubDirPath, file_name))
}

allcause_mrtl_by_dis <- function(outstrata, prbl, sTablesSubDirPath,
                                  pop_denom = F, what = "allcause_mrtl_by_dis",
                                  type = c("ons", "esp")) {
  if("ons" %in% type) {
    tt <- fread(paste0(sSummariesSubDirPath, "/all_cause_mrtl_by_dis_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
                                                                                                dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
    if (all(c("mc", "year", "scenario") %in% outstrata) &
        !(any(c("sex", "agegrp", "dimd") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year (not standardised).csv"
    } else if (all(c("mc", "year", "sex", "scenario") %in% outstrata)&
               !(any(c("agegrp", "dimd") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-sex (not standardised).csv"
    } else if (all(c("mc", "year", "agegrp", "scenario") %in% outstrata)&
               !(any(c("sex", "dimd") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-agegrp (not standardised).csv"
    } else if (all(c("mc", "year", "dimd", "scenario") %in% outstrata)&
               !(any(c("sex", "agegrp") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-dimd (not standardised).csv"
    } else if (all( c("mc", "year", "agegrp", "sex", "scenario") %in% outstrata)&
               !(any(c("dimd") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-agegroup-sex (not standardised).csv"
    } else if (all( c("mc", "year", "agegrp", "dimd", "scenario") %in% outstrata)&
               !(any(c("sex") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-agegroup-dimd (not standardised).csv"
    } else if (all( c("mc", "year", "sex", "dimd", "scenario") %in% outstrata)&
               !(any(c("agegrp") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-sex-dimd (not standardised).csv"
    } else if (all(c("mc", "year", "agegrp", "sex", "dimd", "scenario") %in% outstrata)) {
      file_name <- "/all-cause mrtl by disease-year-agegroup-sex-dimd (not standardised).csv"
    }
    if(pop_denom & "ons" %in% type) {
      pp <- fread(paste0(sSummariesSubDirPath, "/prvl_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
                                                                                 dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
      file_name <- gsub("(not standardised).csv", "popdenom (not standardised).csv", file_name)
    }
  } else if("esp" %in% type) {
    tt <- fread(paste0(sSummariesSubDirPath,"/all_cause_mrtl_by_dis_esp.csv.gz"))[, `:=` (year = year + 2000L,
                                                                                          dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")))]
    if (all(c("mc", "year", "scenario") %in% outstrata) &
        !(any(c("sex", "agegrp", "dimd") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year (age-sex-dimd standardised).csv"
    } else if (all(c("mc", "year", "sex", "scenario") %in% outstrata)&
               !(any(c("agegrp", "dimd") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-sex (age-dimd standardised).csv"
    } else if (all(c("mc", "year", "agegrp", "scenario") %in% outstrata)&
               !(any(c("sex", "dimd") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-agegrp (sex-dimd standardised).csv"
    } else if (all(c("mc", "year", "dimd", "scenario") %in% outstrata)&
               !(any(c("sex", "agegrp") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-dimd (age-sex standardised).csv"
    } else if (all( c("mc", "year", "agegrp", "sex", "scenario") %in% outstrata)&
               !(any(c("dimd") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-agegroup-sex (dimd standardised).csv"
    } else if (all( c("mc", "year", "agegrp", "dimd", "scenario") %in% outstrata)&
               !(any(c("sex") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-agegroup-dimd (sex standardised).csv"
    } else if (all( c("mc", "year", "sex", "dimd", "scenario") %in% outstrata)&
               !(any(c("agegrp") %in% outstrata))) {
      file_name <- "/all-cause mrtl by disease-year-sex-dimd (age standardised).csv"
    }
    # } else if (all(c("mc", "year", "agegrp", "sex", "dimd", "scenario") %in% outstrata)) {
    #   file_name <- "/all-cause mrtl by disease-year-agegroup-sex-dimd (all standardised).csv"
    # }
  }
  d <- tt[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"),
          keyby = eval(outstrata)]
  d <- melt(d, id.vars = outstrata)
  if(!pop_denom) {
    cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
  }
  d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
  if(pop_denom & "ons" %in% type) {
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

dis_chrs <- function(outstrata, prbl, sTablesSubDirPath,
                     what = "dis_char", type = c("ons")) {
  if("ons" %in% type) {
    tt <- fread(paste0(sSummariesSubDirPath,"/dis_characteristics_scaled_up.csv.gz"))[, `:=` (year = year + 2000L,
                                                                                              dimd = factor(dimd, c("1 most deprived", as.character(2:9), "10 least deprived")),
                                                                                              mean_cms_count_cms1st_cont = as.numeric(mean_cms_count_cms1st_cont))]
    if (all(c("mc", "year", "scenario") %in% outstrata) &
        !(any(c("sex", "dimd") %in% outstrata))) {
      file_name <- "/disease characteristics by year (not standardised).csv"
    } else if (all(c("mc", "year", "sex", "scenario") %in% outstrata)&
               !(any(c("dimd") %in% outstrata))) {
      file_name <- "/disease characteristics by year-sex (not standardised).csv"
    } else if (all(c("mc", "year", "dimd", "scenario") %in% outstrata)&
               !(any(c("sex") %in% outstrata))) {
      file_name <- "/disease characteristics by year-dimd (not standardised).csv"
    } else if (all( c("mc", "year", "sex", "dimd", "scenario") %in% outstrata)&
               !(any(c("agegrp") %in% outstrata))) {
      file_name <- "/disease characteristics by year-sex-dimd (not standardised).csv"
    }
  }
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
}

# All-cause mortality by disease not standardised ----
outstrata <- c("mc", "year", "scenario")
allcause_mrtl_by_dis(outstrata, prbl, sTablesSubDirPath,
                                  pop_denom = F, what = "allcause_mrtl_by_dis",
                                  type = c("ons"))

outstrata <- c("mc", "year", "sex", "scenario")

outstrata <- c("mc", "year", "agegrp", "scenario")

outstrata <- c("mc", "year", "agegrp", "sex", "scenario")

outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")


# All-cause mortality by disease not standardised pop denominator----

outstrata <- c("mc", "year", "scenario")
allcause_mrtl_by_dis(outstrata, prbl, sTablesSubDirPath,
                            pop_denom = T, what = "allcause_mrtl_by_dis",
                            type = c("ons"))

outstrata <- c("mc", "year", "sex", "scenario")


outstrata <- c("mc", "year", "agegrp", "sex", "scenario")


outstrata <- c("mc", "year", "agegrp", "sex", "dimd", "scenario")


# All-cause mortality by disease standardised ----

outstrata <- c("mc", "year", "scenario")
allcause_mrtl_by_dis(outstrata, prbl, sTablesSubDirPath,
                            pop_denom = F, what = "allcause_mrtl_by_dis",
                            type = c("esp"))

outstrata <- c("mc", "year", "sex", "scenario")

outstrata <- c("mc", "year", "dimd", "scenario")

outstrata <- c("mc", "year", "sex", "dimd", "scenario")


# Disease characteristics non standardised ----

outstrata <- c("mc", "year", "scenario")
dis_chrs(outstrata, prbl, sTablesSubDirPath,
         what = "dis_char", type = c("ons"))
outstrata <- c("mc", "year", "sex", "scenario")

outstrata <- c("mc", "year", "dimd", "scenario")

outstrata <- c("mc", "year", "sex", "dimd", "scenario")

# XPS ----
outstrata <- c("mc", "year", "scenario")
xps_summ(outstrata, prbl, path,
         sTablesSubDirPath, what = "xps",
         type = c("ons"))

outstrata <- c("mc", "year", "agegrp20", "scenario")

outstrata <- c("mc", "year", "agegrp20", "sex", "scenario")

outstrata <- c("mc", "year", "qimd", "scenario")

# XPS Standardised ----
outstrata <- c("mc", "year", "sex", "qimd", "scenario")
xps_summ(outstrata, prbl, path,
         sTablesSubDirPath, what = "xps",
         type = c("esp"))

outstrata <- c("mc", "year", "scenario")

outstrata <- c("mc", "year", "sex", "scenario")

outstrata <- c("mc", "year", "qimd", "scenario")


