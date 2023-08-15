#Age-sex standardised exposure trends

library(data.table)
library(CKutils)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(scales)
library(cowplot)

bl_pth <- "/mnt/storage_fast/output/hf_real/"
parf_pth <- "/mnt/storage_fast/output/hf_real_parf/"


check_pop_str <- FALSE

if(check_pop_str){
  #Comparing popsize/dist of various scenarios and baseline
  strata <- c("agegrp", "mc", "scenario", "year" ,"sex" ,  "dimd")
  pop_tab <- fread(paste0(parf_pth,"summaries/incd_scaled_up.csv.gz"),
                   select = c("agegrp", "mc", "scenario", "year" ,"sex" ,  "dimd", "popsize"))
  pop_tab <- pop_tab[, .(popsize = sum(popsize)),
                     by = .(agegrp , scenario, year,mc )][,
                                                          .(popsize = quantile(popsize, 0.5)), by = .(agegrp, scenario, year)]
  pop_tab[scenario == "alc", popsize := popsize/2]  #for some reason alc has double the number of people WHY???
  
  library(ggplot2)
  ggplot(data = pop_tab[year == 23], aes( x = scenario , y = popsize, fill = agegrp)) +
    geom_col(position = "fill")
  ggplot(data = pop_tab[year == 23], aes( x = scenario , y = popsize, fill = agegrp)) +
    geom_col(position = "stack")
  
  pop_tab_w <- dcast(pop_tab, formula =agegrp + year  ~ scenario, value.var = "popsize")
  pop_tab_w
}

#Scenarios with risk factors set to optimal levels
optimal_tab <- fread(paste0(parf_pth,"summaries/incd_scaled_up.csv.gz"))
dl <- grep("_prvl$", copy(names(optimal_tab)), value = TRUE)
dl_new <- gsub("_prvl$", "_incd", dl)
dl_gen <- gsub("_prvl$", "", dl)
setnames(optimal_tab, dl, dl_new)

#Prev tab for calculating incidence correctly
prev_tab <- fread(paste0(parf_pth,"summaries/prvl_scaled_up.csv.gz"))

strata <- c("agegrp", "mc", "scenario", "year" ,"sex" ,  "dimd")

optimal_tab <- optimal_tab[prev_tab[, -c("popsize")],
                           on = strata]

#Calculating incidence
incd_tab <- CJ(year = optimal_tab$year,
               mc = optimal_tab$mc,
               scenario = optimal_tab$scenario,
               unique = TRUE)
strata2 <- names(incd_tab)
for(i in dl_gen){
  tmp <- optimal_tab[, .SD, .SDcols = c(strata, "popsize",
                                        paste0(i, "_incd"),
                                        paste0(i, "_prvl")
  )]
  setnames(tmp, c(paste0(i, "_incd"), paste0(i, "_prvl") ),
           c("incd", "prvl"))
  tmp[, atrisk := popsize - prvl + incd]
  tmp <- tmp[, .(sum(incd)/sum(atrisk)), by = strata2]
  setnames(tmp, "V1", paste0(i))
  incd_tab <- incd_tab[tmp, on = strata2]
}


#Calculating PARFs
baseline <- incd_tab[year == 23 & scenario == "sc0", -c("scenario")]
incd_tab_sc <- incd_tab[year == 23 & scenario != "sc0"]
setnames(baseline, dl_gen, paste0(dl_gen, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))

parf_tab <- CJ(year = 23,
               #mc = 1:100,
               scenario = incd_tab_sc$scenario,
               unique = TRUE)

strata3 <- c("mc","year", "scenario")

for(i in dl_gen){
  tmp <- incd_tab_sc[, .SD,
                     .SDcols = c(strata3, i, paste0(i, "_bl"))]
  setnames(tmp, c(i, paste0(i, "_bl")), c("sc", "bl"))
  tmp <- tmp[, .(parf = (bl - sc)/bl), by = strata3][, .(parf = quantile(parf, 0.5, na.rm = TRUE )), keyby = .(year, scenario)]
  setnames(tmp, "parf", i)
  parf_tab <- parf_tab[tmp, on = .(year, scenario)]
}
parf_tab[, year := year + 2000]
cols <- names(parf_tab)[3:37]
parf_tab[,(cols) := round(.SD,4), .SDcols=cols]

# tmp_path <- #because I can't write to the actual directory
#   function(x = character(0))
#     paste0("/mnt/", Sys.info()[["user"]],
#            "/UoL/CPRD2021/alhead_cprd/HFdemandmodeldata/tmpfiles/", x)
#
# fwrite(parf_tab, tmp_path("/PARFs in 2023.csv"))
fwrite(parf_tab, paste0(parf_pth, "postprocess/PARFs in 2023.csv"))



# PARFs by disease at max year of median lag time
#
sum_tab <- CJ(year = c(17, 18, 22),
              mc = 1:200,
              scenario = c(parf_tab$scenario),
              unique = TRUE)

#because going to apply the correction, need to include all risk factors for now
risk_fa <- incd_tab[scenario != "sc0", unique(scenario)]
risk_fa_tab <- data.table() #make a look-up table for the risk factors actually modelled

#Breast Ca  - 9 years
disnm <- "breast_ca"
risk_fa_lst <- c("pa", "alc", "bmi", "smk") #smk & ets combined
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))

#Colorectal Ca  - 9 years
disnm <- "colorectal_ca"
risk_fa_lst <- c("pa", "alc", "bmi", "smk")
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))

#Lung Ca  - 9 years
disnm <- "lung_ca"
risk_fa_lst <- c("frvg", "smk") #smk & ets combined; fruit & veg combined
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario %in% c(risk_fa), .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))

#Prostate Ca  - 9 years
disnm <- "prostate_ca"
risk_fa_lst <- c( "smk")
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario %in% c(risk_fa), .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))


#Dementia  - 9 years
disnm <- "dementia"
risk_fa_lst <- c( "smk", "bmi")
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year ==  (13 + y) & scenario %in% c(risk_fa), .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))

#CHD - 4 years
disnm <- "chd"
risk_fa_lst <- risk_fa[risk_fa != "all"]
y <- 4
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario != "sc0", .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))

#Stroke - 4 years
disnm <- "stroke"
risk_fa_lst <- risk_fa[risk_fa != "all"]
y <- 4
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario != "sc0", .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))

#AF - 4 years
disnm <- "af"
risk_fa_lst <- c("smk", "bmi", "alc", "sbp")
y <- 4
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario %in% risk_fa, .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))


#T2dm - 5 years
disnm <- "t2dm"
risk_fa_lst <- c("smk", "bmi", "alc", "pa", "frvg") #smk & ets combined, fr & veg combined
y <- 5
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario %in% risk_fa, .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))

#Asthma - 5 years
disnm <- "asthma"
risk_fa_lst <- c("smk", "bmi")
y <- 5
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario %in% risk_fa, .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))


#CKD - 5 years
disnm <- "ckd"
risk_fa_lst <- c( "sbp", "bmi")
y <- 5
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario %in% risk_fa, .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))

#COPD - 5 years
disnm <- "copd"
risk_fa_lst <- c("smk") #smk & ets combined
y <- 5
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)

baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario %in% risk_fa, .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))



parf_tab_spec <- CJ(year = c(17,18,22),
                    # mc = parf_tab$mc,
                    scenario = sum_tab$scenario,
                    unique = TRUE)

strata_in <- names(parf_tab_spec)
strata_out <- "year"
prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)

dl_gen_modeled <- dl_gen[dl_gen %in% names(sum_tab)]

parf_tab_sum <- data.table()
parf_tab_sum_corrected <- data.table()
for(i in dl_gen_modeled){
  tmp <- sum_tab[, .SD,
                 .SDcols = c(strata_in, i, paste0(i, "_bl"), "mc")]
  setnames(tmp, c(i, paste0(i, "_bl")), c("sc", "bl"))
  tmp <- tmp[!is.na(sc)][, parf := (bl - sc)/bl][, c("sc", "bl") := NULL]
  setkey(tmp, scenario)
  tmp <- tmp[, fquantile_byid(parf, prbl, id = as.character(scenario)), keyby = eval(setdiff(strata_out, "mc"))]
  setnames(tmp, c(setdiff(strata_out, "mc"), "scenario", percent(prbl, prefix = "parf_")))
  tmp[, disease := paste0(i)]
  parf_tab_sum <- rbind(parf_tab_sum, tmp)
  
  
  #making wide so easier to do the correction
  risk_fa_w <- tmp[scenario != "all", unique(scenario)]
  tmp_w <- dcast(tmp, disease + year  ~ scenario, value.var = "parf_50.0%")
  tmp_w[, sumparf := Reduce("+", lapply(.SD, function(x) x)), .SDcols = risk_fa_w]
  tmp_w[, correction := all/sumparf]
  for (j in risk_fa_w) {
    set(tmp_w, NULL, j, tmp_w[[j]] * tmp_w$correction)
  }
  tmp_w[, c("sumparf", "correction") := NULL]
  parf_tab_sum_corrected <- rbind(parf_tab_sum_corrected, tmp_w )
  
}
parf_tab_sum[, year := year + 2000]
cols <- names(parf_tab_sum)[3:7]
parf_tab_sum[,(cols) := round(.SD,4), .SDcols=cols]

parf_tab_sum_corrected[, year := year + 2000]
cols <- names(parf_tab_sum_corrected)[3:8]
parf_tab_sum_corrected[,(cols) := round(.SD,4), .SDcols=cols]


##### APPLY THESE TO INCIDENCE
#Incidence rates
dl_gen_modeled_prvl <- paste0(dl_gen_modeled, "_prvl")

# tt_incd_rt <- fread("/mnt/storage_fast/output/hf_real 25-09-22/tables/incidence by year (age-sex-dimd standardised).csv")[
#   year %in% c(2017,2018,2022) & scenario == "sc0" & disease %in% dl_gen_modeled_prvl][
#     , .(year, disease, incd_p100k = `incd_rate_50.0%` * 100000)]
tt_incd_rt <- fread(paste0(bl_pth,"tables/incidence by year (not standardised).csv"))[
  year %in% c(2017,2018,2022) & scenario == "sc0" & disease %in% dl_gen_modeled_prvl][
    , .(year, disease, incd_p100k = `incd_rate_50.0%` * 100000)]
tt_incd_rt[, disease := gsub("_prvl", "", disease)]

parf_tab_sum_corrected_shrt <- data.table()
for(i in dl_gen_modeled){
  risk_fa_lst <- risk_fa_tab[disease == i, unique(risk_fa)]
  tmp <- parf_tab_sum_corrected[ disease == i, .SD, .SDcols = c("disease", "year", "all", risk_fa_lst)]
  tmp[, other := all - Reduce("+", lapply(.SD, function(x) x)), .SDcols = risk_fa_lst]
  
  
  tmp <- tt_incd_rt[tmp, on = c("year", "disease")]
  for (j in c(risk_fa_lst, "all", "other")) {
    set(tmp, NULL, j, tmp[[j]] * tmp$incd_p100k)
  }
  tmp[, not_attrib := incd_p100k - all + other][, c("all", "other") := NULL]
  parf_tab_sum_corrected_shrt <- rbind(parf_tab_sum_corrected_shrt, tmp, fill = T)
}
setnafill(parf_tab_sum_corrected_shrt, type = "const", fill = 0L, cols = risk_fa[!risk_fa %in% "all"])
setorder(parf_tab_sum_corrected_shrt, -incd_p100k        )
dis_ord <-parf_tab_sum_corrected_shrt[, disease]

parf_tab_sum_corrected_l <- melt(parf_tab_sum_corrected_shrt,
                                 id.vars = c("year", "disease", "incd_p100k"),
                                 measure.vars = c( "alc", "bmi", "frvg", "pa",
                                                   "smk", "not_attrib",
                                                   "sbp", "tchol"),
                                 variable.name = "risk_fa")

parf_tab_sum_corrected_l[, disease := factor(disease, levels = dis_ord)]
parf_tab_sum_corrected_l[, risk_fa := factor(risk_fa,
                                             levels = c("not_attrib", "smk", "bmi", "sbp", "tchol", "alc", "pa", "frvg"),
                                             labels = c("Not attributable",  "Smoking", "BMI", "SBP", "Total cholesterol",
                                                        "Alcohol", "Physical activity", "Fruit & veg"))]
parf_tab_sum_corrected_l

disnm <- c("af", "asthma", "breast_ca", "chd", "ckd", "colorectal_ca", "copd",
           "dementia", "lung_ca", "prostate_ca" , "stroke" , "t2dm" )

disnm2 <- c("Atrial Fibrillation", "Asthma", "Primary Malignancy_Breast", "CHD",
            "Chronic Kidney Disease",  "Primary Malignancy_Colorectal" , "COPD",
            "Dementia", "Primary Malignancy_Lung", "Primary Malignancy_Prostate",
            "Stroke" ,    "Type 2 Diabetes Mellitus"   )

disnm_long <- c("Atrial\nFibrillation", "Asthma", "Breast\nCancer", "CHD", "CKD",
                "Colorectal\nCancer", "COPD", "Dementia", "Lung\nCancer",
                "Prostate\nCancer" , "Stroke" , "T2\nDiabetes" )
disnm_tab <- data.table(disnm= disnm, disnm2 = disnm2, disnm_long = disnm_long )
disnm_tab[, disnm := factor(disnm,
                            levels = c(dis_ord))]
setorder(disnm_tab, disnm)
disnm_tab[, `:=` (disnm2 = factor(disnm2,
                                  levels = unique(disnm_tab$disnm2)),
                  disnm_long = factor(disnm_long,
                                      levels = unique(disnm_tab$disnm_long)))]
parf_tab_sum_corrected_l <- parf_tab_sum_corrected_l[disnm_tab, on = c("disease" = "disnm")]

#Setting the plot settings

theme_set(new = theme_few())
theme_update(plot.title = element_text(hjust = 0.5))
ggcust <- function(...){
  ggplot(...) +
    scale_fill_brewer(type = "qual") +
    scale_y_continuous(name = "Cases", labels = comma_format())
}

ggplot(parf_tab_sum_corrected_l, aes(x = disnm_long, y = value, fill = risk_fa)) +
  geom_col(position = "stack") +
  scale_y_continuous(name = "Incidence proportion (per 100,000)", labels = comma_format()) +
  scale_x_discrete(name = "Condition", guide = guide_axis(angle = 45 ) ) +
  labs(fill = "Risk factor") +
  ggtitle(paste0("Attibutable cases by risk factor\n(incidence proportion, not standardised)"))
ggsave2(filename = paste0(parf_pth, "postprocess/PARFcases_incd_all_crude.png"), scale = 0.8)
fwrite(parf_tab_sum_corrected_l,paste0(parf_pth, "postprocess/PARFcases_incd_all_crude.csv"))






### By DIMD

#Calculating incidence
incd_tab <- CJ(year = optimal_tab$year,
               dimd = optimal_tab$dimd,
               mc = optimal_tab$mc,
               scenario = optimal_tab$scenario,
               unique = TRUE)
strata2 <- names(incd_tab)
for(i in dl_gen){
  tmp <- optimal_tab[, .SD, .SDcols = c(strata, "popsize",
                                        paste0(i, "_incd"),
                                        paste0(i, "_prvl")
  )]
  setnames(tmp, c(paste0(i, "_incd"), paste0(i, "_prvl") ),
           c("incd", "prvl"))
  tmp[, atrisk := popsize - prvl + incd]
  tmp <- tmp[, .(sum(incd)/sum(atrisk)), by = strata2]
  setnames(tmp, "V1", paste0(i))
  incd_tab <- incd_tab[tmp, on = strata2]
}


#Calculating PARFs in 2023
# baseline <- incd_tab[year == 23 & scenario == "sc0", -c("scenario")]
# incd_tab_sc <- incd_tab[year == 23 & scenario != "sc0"]
# setnames(baseline, dl_gen, paste0(dl_gen, "_bl"))
# absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
#
# parf_tab <- CJ(year = 23,
#                #mc = 1:100,
#                scenario = incd_tab_sc$scenario,
#                unique = TRUE)
#
# strata3 <- names(parf_tab)
#
# for(i in dl_gen){
#   tmp <- incd_tab_sc[, .SD,
#                      .SDcols = c(strata3, i, paste0(i, "_bl"))]
#   setnames(tmp, c(i, paste0(i, "_bl")), c("sc", "bl"))
#   tmp <- tmp[, .(parf = (bl - sc)/bl), by = strata3][, .(parf = quantile(parf, 0.5, na.rm = TRUE )), keyby = .(year, scenario)]
#   setnames(tmp, "parf", i)
#   parf_tab <- parf_tab[tmp, on = .(year, scenario)]
# }
# parf_tab[, year := year + 2000]
# cols <- names(parf_tab)[3:37]
# parf_tab[,(cols) := round(.SD,4), .SDcols=cols]
#
# # tmp_path <- #because I can't write to the actual directory
# #   function(x = character(0))
# #     paste0("/mnt/", Sys.info()[["user"]],
# #            "/UoL/CPRD2021/alhead_cprd/HFdemandmodeldata/tmpfiles/", x)
# #
# # fwrite(parf_tab, tmp_path("/PARFs in 2023.csv"))
# fwrite(parf_tab, paste0(parf_pth, "PARFs in 2023.csv"))



# PARFs by disease at max year of median lag time
#
sum_tab <- CJ(year = c(17, 18, 22),
              mc = 1:200,
              dimd = c(optimal_tab$dimd),
              scenario = c(optimal_tab$scenario),
              unique = TRUE)

#because going to apply the correction, need to include all risk factors for now
# risk_fa <- incd_tab[scenario != "sc0", unique(scenario)]
# risk_fa_tab <- data.table() #make a look-up table for the risk factors actually modelled

strata_sc <- names(copy(sum_tab)) #have to copy it otherwise it does weird stuff with the assign
strata_bl <- strata_sc[strata_sc != "scenario"]


#Breast Ca  - 9 years
disnm <- "breast_ca"
risk_fa_lst <- c("pa", "alc", "bmi", "smk") #smk & ets combined
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)

#Colorectal Ca  - 9 years
disnm <- "colorectal_ca"
risk_fa_lst <- c("pa", "alc", "bmi", "smk")
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)

#Lung Ca  - 9 years
disnm <- "lung_ca"
risk_fa_lst <- c("frvg", "smk") #smk & ets combined; fruit & veg combined
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)

#Prostate Ca  - 9 years
disnm <- "prostate_ca"
risk_fa_lst <- c( "smk")
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)


#Dementia  - 9 years
disnm <- "dementia"
risk_fa_lst <- c( "smk", "bmi")
y <- 9
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)

#CHD - 4 years
disnm <- "chd"
risk_fa_lst <- risk_fa[risk_fa != "all"]
y <- 4
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)

#Stroke - 4 years
disnm <- "stroke"
risk_fa_lst <- risk_fa[risk_fa != "all"]
y <- 4
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)

#AF - 4 years
disnm <- "af"
risk_fa_lst <- c("smk", "bmi", "alc", "sbp")
y <- 4
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)


#T2dm - 5 years
disnm <- "t2dm"
risk_fa_lst <- c("smk", "bmi", "alc", "pa", "frvg") #smk & ets combined, fr & veg combined
y <- 5
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)

#Asthma - 5 years
disnm <- "asthma"
risk_fa_lst <- c("smk", "bmi")
y <- 5
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == (13 + y) & scenario == "sc0", .SD, .SDcols = c("year", "mc", disnm)]
incd_tab_sc <- incd_tab[year == (13 + y) & scenario %in% risk_fa, .SD, .SDcols = c("year", "mc", "scenario", disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = c("year", "mc"))
absorb_dt(sum_tab,incd_tab_sc, on = c("year", "mc", "scenario"))


#CKD - 5 years
disnm <- "ckd"
risk_fa_lst <- c( "sbp", "bmi")
y <- 5
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)
baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)

#COPD - 5 years
disnm <- "copd"
risk_fa_lst <- c("smk") #smk & ets combined
y <- 5
tmp <- CJ(disease = disnm,
          year = 13 + y,
          risk_fa = risk_fa_lst,
          unique = TRUE)
risk_fa_tab <- rbind(risk_fa_tab, tmp)

baseline <- incd_tab[year == 13 + y & scenario == "sc0", .SD, .SDcols = c(strata_bl, disnm)]
incd_tab_sc <- incd_tab[year == 13 + y & scenario %in% c(risk_fa), .SD, .SDcols = c(strata_sc, disnm)]
setnames(baseline, disnm, paste0(disnm, "_bl"))
absorb_dt(incd_tab_sc, baseline, on = strata_bl)
absorb_dt(sum_tab,incd_tab_sc, on = strata_sc)

parf_tab_spec <- CJ(year = c(17,18,22),
                    # mc = parf_tab$mc,
                    dimd = sum_tab$dimd,
                    scenario = sum_tab$scenario,
                    unique = TRUE)

strata_in <- names(parf_tab_spec)
strata_out <- c("year", "dimd")
prbl = c(0.5, 0.025, 0.975, 0.1, 0.9)

dl_gen_modeled <- dl_gen[dl_gen %in% names(sum_tab)]

parf_tab_sum <- data.table()
parf_tab_sum_corrected <- data.table()
for(i in dl_gen_modeled){
  tmp <- sum_tab[, .SD,
                 .SDcols = c(strata_in, i, paste0(i, "_bl"), "mc")]
  setnames(tmp, c(i, paste0(i, "_bl")), c("sc", "bl"))
  tmp <- tmp[!is.na(sc)][, parf := (bl - sc)/bl][, c("sc", "bl") := NULL]
  setkey(tmp, scenario)
  tmp <- tmp[, fquantile_byid(parf, prbl, id = as.character(scenario)), keyby = eval(setdiff(strata_out, "mc"))]
  setnames(tmp, c(setdiff(strata_out, "mc"), "scenario", percent(prbl, prefix = "parf_")))
  tmp[, disease := paste0(i)]
  parf_tab_sum <- rbind(parf_tab_sum, tmp)
  
  
  #making wide so easier to do the correction
  risk_fa_w <- tmp[scenario != "all", unique(scenario)]
  tmp_w <- dcast(tmp, disease + year +dimd ~ scenario, value.var = "parf_50.0%")
  tmp_w[, sumparf := Reduce("+", lapply(.SD, function(x) x)), .SDcols = risk_fa_w]
  tmp_w[, correction := all/sumparf]
  for (j in risk_fa) {
    set(tmp_w, NULL, j, tmp_w[[j]] * tmp_w$correction)
  }
  tmp_w[, c("sumparf", "correction") := NULL]
  parf_tab_sum_corrected <- rbind(parf_tab_sum_corrected, tmp_w )
  
}
parf_tab_sum[, year := year + 2000]
cols <- names(parf_tab_sum)[4:8]
parf_tab_sum[,(cols) := round(.SD,4), .SDcols=cols]

parf_tab_sum_corrected[, year := year + 2000]
cols <- names(parf_tab_sum_corrected)[4:9]
parf_tab_sum_corrected[,(cols) := round(.SD,4), .SDcols=cols]


##### APPLY THESE TO INCIDENCE
#Incidence rates
dl_gen_modeled_prvl <- paste0(dl_gen_modeled, "_prvl")
tt_incd_rt <- fread(paste0(bl_pth,"tables/incidence by year-dimd (not standardised).csv"))[
  year %in% c(2017,2018,2022) & scenario == "sc0" & disease %in% dl_gen_modeled_prvl][
    , .(year, disease, dimd, incd_p100k = `incd_rate_50.0%` * 100000)]
tt_incd_rt[, disease := gsub("_prvl", "", disease)]


parf_tab_sum_corrected[, dimd := factor(dimd,
                                        levels = unique(tt_incd_rt$dimd))]
parf_tab_sum_corrected_shrt <- data.table()
for(i in dl_gen_modeled){
  risk_fa_lst <- risk_fa_tab[disease == i, unique(risk_fa)]
  tmp <- parf_tab_sum_corrected[ disease == i, .SD, .SDcols = c("disease", "year", "dimd","all", risk_fa_lst)]
  tmp[, other := all - Reduce("+", lapply(.SD, function(x) x)), .SDcols = risk_fa_lst]
  
  
  tmp <- tt_incd_rt[tmp, on = c("year", "disease", "dimd")]
  for (j in c(risk_fa_lst, "all", "other")) {
    set(tmp, NULL, j, tmp[[j]] * tmp$incd_p100k)
  }
  tmp[, not_attrib := incd_p100k - all + other][, c("all", "other") := NULL]
  parf_tab_sum_corrected_shrt <- rbind(parf_tab_sum_corrected_shrt, tmp, fill = T)
}
setnafill(parf_tab_sum_corrected_shrt, type = "const", fill = 0L, cols = risk_fa[risk_fa != "all"])
setorder(parf_tab_sum_corrected_shrt, -incd_p100k        )
dis_ord <-parf_tab_sum_corrected_shrt[, unique(disease)]

parf_tab_sum_corrected_l <- melt(parf_tab_sum_corrected_shrt,
                                 id.vars = c("year", "disease", "dimd", "incd_p100k"),
                                 measure.vars = c( "alc", "bmi", "frvg", "pa",
                                                   "smk", "not_attrib",
                                                   "sbp", "tchol"),
                                 variable.name = "risk_fa")

parf_tab_sum_corrected_l[, disease := factor(disease, levels = dis_ord)]
parf_tab_sum_corrected_l[, risk_fa := factor(risk_fa,
                                             levels = c("not_attrib", "smk", "bmi", "sbp", "tchol", "alc", "pa", "frvg"),
                                             labels = c("Not attributable",  "Smoking", "BMI", "SBP", "Total cholesterol",
                                                        "Alcohol", "Physical activity", "Fruit & veg"))]
parf_tab_sum_corrected_l

disnm <- c("af", "asthma", "breast_ca", "chd", "ckd", "colorectal_ca", "copd",
           "dementia", "lung_ca", "prostate_ca" , "stroke" , "t2dm" )

disnm2 <- c("Atrial Fibrillation", "Asthma", "Primary Malignancy_Breast", "CHD",
            "Chronic Kidney Disease",  "Primary Malignancy_Colorectal" , "COPD",
            "Dementia", "Primary Malignancy_Lung", "Primary Malignancy_Prostate",
            "Stroke" ,    "Type 2 Diabetes Mellitus"   )

disnm_long <- c("Atrial\nFibrillation", "Asthma", "Breast\nCancer", "CHD", "CKD",
                "Colorectal\nCancer", "COPD", "Dementia", "Lung\nCancer",
                "Prostate\nCancer" , "Stroke" , "T2\nDiabetes" )
disnm_tab <- data.table(disnm= disnm, disnm2 = disnm2, disnm_long = disnm_long )
disnm_tab[, disnm := factor(disnm,
                            levels = c(dis_ord))]
setorder(disnm_tab, disnm)
disnm_tab[, `:=` (disnm2 = factor(disnm2,
                                  levels = unique(disnm_tab$disnm2)),
                  disnm_long = factor(disnm_long,
                                      levels = unique(disnm_tab$disnm_long)))]
parf_tab_sum_corrected_l <- parf_tab_sum_corrected_l[disnm_tab, on = c("disease" = "disnm")]
parf_tab_sum_corrected_l[, dimd := factor(dimd,
                                          levels = unique(tt_incd_rt$dimd),
                                          labels = c("1 most\ndeprived", "2", "3", "4", "5",
                                                     "6", "7", "8", "9", "10 least\ndeprived"))]
#Setting the plot settings

theme_set(new = theme_few())
theme_update(plot.title = element_text(hjust = 0.5),
             axis.text.x = element_text(size = 9))
ggcust <- function(...){
  ggplot(...) +
    scale_fill_brewer(type = "qual") +
    scale_y_continuous(name = "Cases", labels = comma_format())
  
}

ggplot(parf_tab_sum_corrected_l, aes(x = dimd, y = value, fill = risk_fa)) +
  geom_col(position = "stack") +
  facet_wrap(~disnm_long, scales = "free") +
  scale_y_continuous(name = "Incidence proportion (per 100,000)", labels = comma_format()) +
  scale_x_discrete(name = "Decile of IMD", guide = guide_axis(angle = 45 ) ) +
  labs(fill = "Risk factor") +
  ggtitle(paste0("Attibutable cases by risk factor and dimd\n(incidence proportion, not standardised)"))
ggsave2(filename = paste0(parf_pth, "postprocess/PARFcases_incd_imd_crude.png"), scale = 2)
fwrite(parf_tab_sum_corrected_l,paste0(parf_pth, "postprocess/PARFcases_incd_imd_crude.csv"))








#For mortality

#Scenarios with risk factors set to optimal levels
optimal_tab <- fread(paste0(parf_pth,"summaries/mrtl_scaled_up.csv.gz"))

strata <- c("agegrp", "mc", "scenario", "year" ,"sex" ,  "dimd")
strata2 <- c("mc", "scenario", "year" )

#Calculating mort

mort_tab <- optimal_tab[, .(mrtl = sum(all_cause_mrtl)/sum(popsize )), by = strata2]

#Calculating PARFs

baseline <- mort_tab[year == 23 & scenario == "sc0", -c("scenario")]
baseline_sc <- mort_tab[year == 23 & scenario != "sc0"]
setnames(baseline, "mrtl", "mrtl_bl")
absorb_dt(baseline_sc, baseline, on = c("year", "mc"))


strata3 <- c("year", "scenario")

parf_tab <- baseline_sc[, .(parf = (mrtl_bl - mrtl)/mrtl_bl), by = strata3][
  , .(parf = quantile(parf, 0.5, na.rm = TRUE )), keyby = .(year, scenario)]

parf_tab[, parf:= round(parf,4)]

# tmp_path <- #because I can't write to the actual directory
#   function(x = character(0))
#     paste0("/mnt/", Sys.info()[["user"]],
#            "/UoL/CPRD2021/alhead_cprd/HFdemandmodeldata/tmpfiles/", x)
#
# fwrite(parf_tab, tmp_path("/PARFs in 2023.csv"))
fwrite(parf_tab, paste0(parf_pth, "postprocess/Mortality PARFs in 2023.csv"))

