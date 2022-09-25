library(data.table) # for fast data manipulation
library(fst) # Fast way to save and load data.tables
library(gamlss)
library(qs)
library(foreach)
library(doParallel)
registerDoParallel(5L)
setDTthreads(10)


disnm <- "nonmodelled" # disease name
overwrite_ftlt <- FALSE
overwrite_dpnd <- TRUE

# disease list
# "Anxiety_Depression"            "Asthma"
# "Atrial Fibrillation"           "CHD"
# "COPD"                          "Chronic Kidney Disease"
# "Dementia"                      "Heart failure"
# "Hypertension"                  "Obesity"
# "Other cancers"                 "Primary Malignancy_Breast"
# "Primary Malignancy_Colorectal" "Primary Malignancy_Lung"
# "Primary Malignancy_Prostate"   "Stroke"
# "Type 2 Diabetes Mellitus"

strata <- c("year", "age", "sex", "dimd", "sha", "ethnicity")
strata_ftlt <- c("year", "age", "sex", "dimd", "sha", "ethnicity")
strata_dpnd <- c("age", "sex", "dimd")

source(paste0("/mnt/", Sys.info()[["user"]], "/UoL/CPRD2021/epi_models/scripts/aux_fn.R"))


# Case Fatality 1st year ====
if (overwrite_ftlt ||
    # !file.exists(output_path(paste0(disnm, "_ftlt1.qs"))) ||
    !file.exists(output_path(paste0(disnm, "_ftlt2.qs")))) {


  # 2+ years case fatality
  dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                           as.data.table = TRUE)[gender != "I"]
  )[between(age, 20, 100)  &
      year < 2020, .SD, .SDcols = c(strata_ftlt, "death_cause", "diedinyr", "diedinCPRD")
  ][, .(ftlt2 = sum(diedinyr == 1L & (is.na(death_cause) | death_cause %in% c("Obesity", "Hypertension")), na.rm = TRUE),
        n = .N), keyby = strata_ftlt]
  dt[, no_ftlt2 := n - ftlt2]
  dt[, year := year - 2000]
  y <- cbind(dt$ftlt2, dt$no_ftlt2)
  dt[, c("ftlt2", "no_ftlt2", "n") := NULL]
  mod_max <- gamlss(
    y ~ (log(year) + pb(age) + pcat(sex) + dimd + pcat(sha) + pcat(ethnicity)) ^ 2,
    family = BI(),
    data = dt,
    method = mixed(20, 100)
  )


  validate_plots(dt, y, mod_max, "_ftlt2", disnm, strata_ftlt)
  qsave(mod_max, output_path(paste0(disnm, "_ftlt2.qs")), "archive")
  print(paste0(disnm, "_ftlt2 model saved!"))

  trms <- all.vars(formula(mod_max))[-1] # -1 excludes dependent var
  newdata2 <-
    CJ(
      age = 20:100,
      year = 3:100,
      sex = levels(dt$sex),
      dimd = levels(dt$dimd),
      ethnicity = levels(dt$ethnicity),
      sha = levels(dt$sha)
    )
  newdata2 <- split(newdata2, by = "dimd")
  newdata2 <-
    # assignment necessary! Copies of data.tables are happening
    lapply(newdata2, function(x)
      x[, "mu2" := predictAll(mod_max, .SD, data = dt), .SDcols = trms])
  newdata2 <- rbindlist(newdata2)


  newdata2[, dimd := factor(dimd, as.character(1:10))]
  setkeyv(newdata2, strata_ftlt)
  write_fst(newdata2, output_path(paste0(disnm, "_ftlt.fst")), 100L)
  print(paste0(disnm, "_ftlt model saved!"))
  rm(dt, mod_max, newdata2, trms)
}



# Dependencies ====
if (overwrite_dpnd ||
    !file.exists(output_path(paste0(disnm, "_dpnd.qs")))) {


  dpnd <- c("Alcohol problems", "Atrial Fibrillation" ,"CHD", "COPD"  ,
            "Chronic Kidney Disease" , "Connective tissue disorder", "Dementia" ,
            "Diabetes excl Type 2", "Epilepsy", "Hearing loss", "Heart failure",
            "IBS" , "Other cancers" , "Primary Malignancy_Breast",   "Primary Malignancy_Colorectal"  , "Primary Malignancy_Lung" ,
            "Primary Malignancy_Prostate"  , "Rheumatoid Arthritis"   ,         "Stroke" ,
            "Type 2 Diabetes Mellitus"  ,      "Hypertension_camb"       ,        "Psychosis_camb"            , "Asthma - spell"       ,
            "Anxiety_Depression - spell"    ,  "Pain"   ,     "Constipation")

  dpnd2 <- c("alcpr", "af" ,"chd", "copd"  ,
            "ckd" , "ctd", "dementia" ,
            "t1dm", "epilepsy", "helo", "hf",
            "ibs" , "other_ca" , "breast_ca",   "colorectal_ca"  , "lung_ca" ,
            "prostate_ca"  , "ra" , "stroke" ,
            "t2dm"  ,  "htn"   ,  "psychosis"   , "asthma" ,
            "andep" ,  "pain"  , "constipation")


  dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                           as.data.table = TRUE)[gender != "I"])[
                             between(age, 20, 99) &
                               year < 2020 &
                               (is.na(death_cause) | death_cause %in% c("Obesity", "Hypertension")) &
                               (is.na(diedinCPRD ) | diedinCPRD == 1), # Need to think about this ,
                             .SD,
                             .SDcols = c(#"death_cause",
                                         "diedinyr",
                                         #"diedinCPRD",
                                         dpnd,
                                         strata_dpnd)]



  CKutils::to_agegrp(dt, 5, 99, agegrp_colname = "agegroup")
  strata_dpnd <- gsub("^age$", "agegroup", strata_dpnd)

  for (j in dpnd) { #Need to make the exposures binary
    set(dt, NULL, j, fifelse(is.na(dt[[j]]), 0L, dt[[j]])) #setting those without spells to 0
    set(dt, NULL, j, fifelse(dt[[j]] > 0L, 1L, 0L)) #want incidence & prev
  }
  rm(j)

  #joining ra & ctd together
  dt[`Rheumatoid Arthritis` == 1L, `Connective tissue disorder` := 1L ][,
    `Rheumatoid Arthritis` := NULL]
  dpnd <- dpnd[!dpnd %in%  "Rheumatoid Arthritis" ]
  dpnd2 <- dpnd2[!dpnd2 %in%  "ra" ]

  adjusted <- CJ(sex = dt$sex,
                 agegroup = dt$agegroup,
                 dimd = dt$dimd,
                 unique = TRUE)

  disnm2 <- "diedinyr" # this is actually the outcome
  setnames(dt, dpnd, dpnd2) #getting rid of spaces so the collapse fn works


  results <- data.table()
  results <- foreach(j = 1:length(dpnd2), .combine = rbind) %dopar% {
    #Alternating through exposure condtiions
    tmptab <- copy(adjusted) #otherwise it just reoverwrites everything
    exps <- c(dpnd2[j], dpnd2[-j])
    frm <- as.formula(paste0(disnm2, "~", paste(exps, collapse = "+"), "+",
                             paste(strata_dpnd, collapse = "+")))

    tmp1 <- glm(frm, data = dt, family = binomial())
    summary(tmp1)
    qsave(tmp1, output_path(paste0("dependencies/mortality/",disnm, "fit", dpnd[j],".qs")), nthreads = 10)
    print(paste0(disnm, j, dpnd[j], "_dpndRR model fit saved!", Sys.time()))


    #Output the RR for all covariate options
    for(i in 1:nrow(adjusted)){
      output <- RRfn(tmp1, data = dt, fixcov = adjusted[i])
      tmptab[i, `:=` (rr = output$RR, lci_rr = output$LCI_RR, uci_rr = output$UCI_RR, var = output$delta.var, dpnd_on = dpnd2[j])]
    }
    qsave(tmptab, output_path(paste0("dependencies/mortality/",disnm, "_dpndRR_", dpnd[j],".qs")), nthreads = 10)
    print(paste0(disnm, j, dpnd[j], "_dpndRR table saved!", Sys.time()))

    tmptab
    }


  #Swapping the dimds round for results purposes
  dimdlabs <- c("1 most deprived" ,  "2" ,   "3","4","5" ,  "6"  , "7" ,  "8" , "9" , "10 least deprived")
  results[, dimd := factor(dimd,
                           levels = 10:1,
                           labels = dimdlabs)]

  results[dpnd_on == "breast_ca" & sex == "men", c("rr", "lci_rr" , "uci_rr") := 1][, var := 0]
  results[dpnd_on == "prostate_ca" & sex == "women", c("rr", "lci_rr" , "uci_rr") := 1][, var := 0]


  qsave(results, output_path(paste0("dependencies/",disnm, "_dpndRR.qs")), nthreads = 10)
  print(paste0(disnm, "_dpndRR table saved!"))

  #Only want to save the .csvy files if some RRs are statistically sig
  results[, statsig := ifelse( #add a flag
    (lci_rr < 1 & uci_rr <1) | (lci_rr > 1 & uci_rr >1),
    1, 0)]
  keep <- results[, max(statsig), by = dpnd_on][V1 == 1, dpnd_on]

  results <- results[dpnd_on %in% keep , .(sex, agegroup, dimd, rr = round(rr, digits = 2), ci_rr = round(uci_rr, digits = 2), dpnd_on )]
  dpnd2 <- dpnd2[dpnd2 %in% keep]

  type <- "_prvl"

  for(j in 1:length(dpnd2)){
    write_xps_tmplte_file(results, j, dpnd2, disnm, type, output_path(paste0("dependencies/",dpnd2[j],"~",disnm,".csvy")))
  }
  rm(j)

}
