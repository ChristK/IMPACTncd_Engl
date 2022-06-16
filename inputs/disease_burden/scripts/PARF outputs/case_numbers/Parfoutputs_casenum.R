# PARF outputs as case numbers 

library(data.table) # for fast data manipulation
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(scales)
library(cowplot)
library(fst)


# folder for plots 
output_dir <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2021/epi_models/scripts/PARF outputs/case_numbers/", x)


#Setting the plot settings
theme_set(new = theme_few())
theme_update(plot.title = element_text(hjust = 0.5))
ggcust <- function(...){
  ggplot(...) +
    scale_fill_brewer(type = "qual") +
    scale_y_continuous(name = "Cases", labels = comma_format())    
}

#parf <- fread("/mnt/storage_fast/output/hf_real_parf/parf/parf_final_corrected.csv") 
parf <- fread("/mnt/storage_fast/output/hf_real_parf_24May/parf/parf_final_corrected.csv") 


strata <- c("year", "age", "sex", "dimd", "sha", "ethnicity")
source(paste0("/mnt/", Sys.info()[["user"]], "/UoL/CPRD2021/epi_models/scripts/aux_fn.R"))

#Names for graphs 
disnm <- c("af", "asthma", "breast_ca", "chd", "ckd", "colorect_ca", "copd",
           "dementia", "lung_ca", "prostate_ca" , "stroke" , "t2dm" )

disnm2 <- c("Atrial Fibrillation", "Asthma", "Primary Malignancy_Breast", "CHD",
            "Chronic Kidney Disease",  "Primary Malignancy_Colorectal" , "COPD",
            "Dementia", "Primary Malignancy_Lung", "Primary Malignancy_Prostate",
            "Stroke" ,    "Type 2 Diabetes Mellitus"   )

disnm_long <- c("Atrial\nFibrillation", "Asthma", "Breast\nCancer", "CHD", "CKD", 
                "Colorectal\nCancer", "COPD", "Dementia", "Lung\nCancer", 
                "Prostate\nCancer" , "Stroke" , "T2\nDiabetes" )

popest <- fread("/mnt/alhead/UoL/CPRD2021/epi_models/scripts/PARF outputs/case_numbers/pop_estimates.csv")
popest <- popest[year == 2013]



parfsum <- parf[, weighted.mean(parf, pop_size), keyby = .(disease,  age,sex,dimd, mc)][
  , .(parf = quantile(V1,0.5)), keyby = .(disease,  age,sex,dimd)]
parfsum


# Prevalence

dt_prev <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                         as.data.table = TRUE)[
                           gender != "I" & year == 2013])[
                             between(age, 20, 100), .SD, .SDcols = c(disnm2, strata)]

dt_prev <- na.omit(dt_prev )
setnames(dt_prev, disnm2, disnm)



dimdlabs <- c("1 most deprived" ,  "2" ,   "3","4","5" ,  "6"  , "7" ,  "8" , "9" , "10 least deprived")
dt_prev[, dimd := factor(dimd,
                         levels = 10:1,
                         labels = dimdlabs)]


cprdpop <- dt_prev[, .N, keyby = .(age, sex,dimd)]
cprdpop[popest, on = c("age", "sex","dimd"), ONSpop_size := i.pop_size]
cprdpop[, wt := ONSpop_size/N]


prev_tab <- data.table()
for (i in disnm){
subtab <- dt_prev[, .(cases = sum(get(i) != 0L)), 
                  keyby = .(age, sex, dimd) ][
                    , disease := paste0(i)]
prev_tab <- rbind(prev_tab,subtab)
prev_tab  
}

prev_tab[parfsum, on = c("disease","age", "sex", "dimd"), parf := i.parf]
prev_tab <- prev_tab[!is.na(parf), ]
prev_tab[, `:=` (attrib = cases * parf,
                 nonattrib = cases * (1-parf))]


prev_tab[cprdpop, 
            on = c("age", "sex","dimd"), 
            wt := i.wt]
#prev_tab[, wt := wt*cases/sum(wt, na.rm = TRUE), by = disease]

#checking that the weights add up to the number of people in each year. 
#prev_tab[, .(.N, sum(wt)), by = disease]




prev_tab_sum <- prev_tab[, .(cases = sum(cases*wt),
                             attrib  = sum(attrib*wt),
                             nonattrib = sum(nonattrib*wt)),
                         keyby = .(disease,dimd)]
prev_tab_sum <- melt(prev_tab_sum, 
                     measure.vars = c("attrib", "nonattrib"), 
                     value.name = "casenum",
                     variable.name = "type")
#prev_tab_sum[, cases := NULL]
prev_tab_sum[, type:= factor(type,
                             levels = c("nonattrib", "attrib"),
                             labels = c("Not preventable", "Preventable"))]

prev_tab_sum[, disease := factor(disease,
                                 levels = c("prostate_ca", "lung_ca", "colorect_ca",  "dementia", "breast_ca",  "af" , "copd",   
                                            "stroke"  ,    "asthma"   ,"chd"      ,    "t2dm"      ,   "ckd" ),
                                 labels = c("prostate_ca", "lung_ca", "colorect_ca",  "dementia", "breast_ca",  "af" , "copd",   
                                            "stroke"  ,    "asthma"   ,"chd"      ,    "t2dm"      ,   "ckd" ))]


labnames <- c("Prostate\nCancer", "Lung\nCancer", "Colorectal\nCancer",  "Dementia", "Breast\nCancer",  "Atrial\nFibrillation" , "COPD",   
           "Stroke"  ,   "Asthma"  ,"CHD"      ,    "T2\nDiabetes"      ,   "CKD" )

ggcust(prev_tab_sum[, .(cases = sum(casenum)), keyby = .(disease, type)], 
       aes(x = disease, y = cases, fill =type)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = "Condition", guide = guide_axis(angle = 45 ), labels = labnames ) + 
  scale_fill_brewer(type = "qual", palette =  "Paired") +
    ggtitle(paste0("Preventable cases")) + 
  theme(legend.title = element_blank())

ggsave2(filename = output_dir(paste0("PARFcases_prev_all.png")), scale = 0.8)


#Can't get the labeller function to work so being lazy
labnames2 <- c("Prostate Cancer", "Lung Cancer", "Colorectal Cancer",  "Dementia", "Breast Cancer",  "Atrial Fibrillation" , "COPD",   
              "Stroke"  ,   "Asthma"  ,"CHD"      ,    "T2 Diabetes"      ,   "CKD" )

prev_tab_sum[, disease2 := disease]
prev_tab_sum[, disease2 := factor(disease2,
                                 levels = c("prostate_ca", "lung_ca", "colorect_ca",  "dementia", "breast_ca",  "af" , "copd",   
                                            "stroke"  ,    "asthma"   ,"chd"      ,    "t2dm"      ,   "ckd" ),
                                 labels = labnames2)]


ggcust(prev_tab_sum[, .(cases = sum(casenum)), keyby = .(disease2, type, dimd)], 
       aes(x = dimd, y = cases, fill =type)) +
  facet_wrap( vars(disease2)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = "IMD decile", guide = guide_axis(angle = 45 ) ) +
  scale_fill_brewer(type = "qual", palette =  "Paired") +
  ggtitle(paste0("Preventable cases by condition and IMD")) + 
  theme(legend.title = element_blank())
ggsave2(filename = output_dir(paste0("PARFcases_prev_imd_all.png")), scale = 0.8)



#Stacked with risk factor attributions 
#Renaming risk factors for graphs 
riskfactors <- c("overall", "physical activity", "alcohol",  "BMI",  "fruit",  
                 "SBP", "total cholesterol",  "veg" , "smoking" , "ETS")
setnames(parf, 
         c("parf", "parf_pa", "parf_alc",  "parf_bmi",  "parf_fru",  
           "parf_sbp", "parf_tch",  "parf_veg" , "parf_smo" , "parf_ets"),
         riskfactors)


prev_tab2 <- prev_tab[, .(disease,age, sex,dimd, cases, wt)]

parfbyrisk <- data.table()
for (i in riskfactors[riskfactors != "overall"]){
parfdata <- parf[, weighted.mean(get(i)* parf_correction, pop_size), keyby = .(disease,  age,sex,dimd, mc)][
  , .(parf = quantile(V1,0.5, na.rm = T)), keyby = .(disease,  age,sex,dimd)][, riskfactor := paste0(i)]
subtab <- copy(prev_tab2)
subtab[parfdata, on = c("disease","age", "sex", "dimd"), parf := i.parf]
subtab <- subtab[!is.na(parf), ]
subtab[, attrib := cases * parf][, riskfactor := paste0(i)]
parfbyrisk <- rbind(parfbyrisk, subtab)
}




parfbyrisksum <- parfbyrisk[, .( attrib  = sum(attrib*wt)),
                         keyby = .(disease,riskfactor, dimd)]
parfbyrisksum <- rbind(parfbyrisksum, prev_tab_sum[type == "Not preventable"  , .(disease, dimd, attrib = casenum, riskfactor = "Not preventable")])
parfbyrisksum[, riskfactor:= factor(riskfactor,
                             levels = c("Not preventable", "alcohol",  "BMI",  "ETS", "fruit",  "physical activity",
                                          "SBP", "smoking" , "total cholesterol",  "veg"   ),
                             labels = c("Not preventable", "Alcohol",  "BMI",  "ETS", "Fruit",  "Physical activity",
                                        "SBP", "Smoking" , "Total cholesterol",  "Veg"  ))]

parfbyrisksum[, disease := factor(disease,
                                 levels = c("prostate_ca", "lung_ca", "colorect_ca",  "dementia", "breast_ca",  "af" , "copd",   
                                            "stroke"  ,    "asthma"   ,"chd"      ,    "t2dm"      ,   "ckd" ),
                                 labels = labnames2)]


ggcust(parfbyrisksum[, .(cases = sum(attrib)), keyby = .(disease, riskfactor)], 
       aes(x = disease, y = cases, fill =riskfactor)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = "Condition", guide = guide_axis(angle = 45 ), labels = labnames ) +
  scale_fill_brewer(type = "qual", palette = "Paired") + 
  labs(fill = "Preventable\nrisk factor") +
  ggtitle(paste0("Preventable cases by risk factor")) 
ggsave2(filename = output_dir(paste0("PARFcases_prev_riskfactor.png")), scale = 0.8)


ggcust(parfbyrisksum[, .(cases = sum(attrib)), keyby = .(disease, riskfactor, dimd)], 
       aes(x = dimd, y = cases, fill =riskfactor)) +
  geom_col(position = "stack") +
  facet_wrap( vars(disease)) +
  scale_x_discrete(name = "IMD decile", guide = guide_axis(angle = 45 ) ) +
  scale_fill_brewer(type = "qual", palette = "Paired") + 
  labs(fill = "Preventable\nrisk factor") +
  ggtitle(paste0("Preventable cases by risk factor and IMD")) 
ggsave2(filename = output_dir(paste0("PARFcases_prev_riskfactor_imd.png")), scale = 0.8)






## Incidence


dt_inc <- harmonise(read_fst(input_path("panel_short_inc.fst"),
                              as.data.table = TRUE)[
                                gender != "I" & year == 2013])[
                                  between(age, 20, 100), .SD, .SDcols = c(disnm2, strata)]

dt_inc <- na.omit(dt_inc )
setnames(dt_inc, disnm2, disnm)


dimdlabs <- c("1 most deprived" ,  "2" ,   "3","4","5" ,  "6"  , "7" ,  "8" , "9" , "10 least deprived")
dt_inc[, dimd := factor(dimd,
                         levels = 10:1,
                         labels = dimdlabs)]

inc_tab <- data.table()
for (i in disnm){
  subtab <- dt_inc[, .(cases = sum(get(i) == 1L)), 
                    keyby = .(age, sex, dimd) ][
                      , disease := paste0(i)]
  inc_tab <- rbind(inc_tab,subtab)
  inc_tab  
}

inc_tab[parfsum, on = c("disease","age", "sex", "dimd"), parf := i.parf]
inc_tab <- inc_tab[!is.na(parf), ]
inc_tab[, `:=` (attrib = cases * parf,
                 nonattrib = cases * (1-parf))]

#Are these the right weights? These are weights of the whole cprd population (i.e. the same I used for prevalence)
inc_tab[cprdpop, 
         on = c("age", "sex","dimd"), 
         wt := i.wt]



inc_tab_sum <- inc_tab[, .(cases = sum(cases*wt),
                             attrib  = sum(attrib*wt),
                             nonattrib = sum(nonattrib*wt)),
                         keyby = .(disease,dimd)]
inc_tab_sum <- melt(inc_tab_sum, 
                     measure.vars = c("attrib", "nonattrib"), 
                     value.name = "casenum",
                     variable.name = "type")
inc_tab_sum[, type:= factor(type,
                             levels = c("nonattrib", "attrib"),
                             labels = c("Not preventable", "Preventable"))]

inc_tab_sum[, disease := factor(disease,
                                 levels = c("prostate_ca", "colorect_ca", "lung_ca",  "breast_ca", "dementia", "copd", "af" ,    
                                            "asthma"   , "stroke"  ,    "chd"      ,    "t2dm"      ,   "ckd" ),
                                 labels = c("prostate_ca", "colorect_ca", "lung_ca",  "breast_ca", "dementia",  "copd", "af" ,  
                                            "asthma" , "stroke"      ,"chd"      ,    "t2dm"      ,   "ckd" ))]


labnames <- c("Prostate\nCancer", "Colorectal\nCancer", "Lung\nCancer",   "Breast\nCancer", "Dementia",  "COPD",   "Atrial\nFibrillation" ,
              "Asthma"  , "Stroke"    ,"CHD"      ,    "T2\nDiabetes"      ,   "CKD" )

ggcust(inc_tab_sum[, .(cases = sum(casenum)), keyby = .(disease, type)], 
       aes(x = disease, y = cases, fill =type)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = "Condition", guide = guide_axis(angle = 45 ), labels = labnames ) + 
  scale_fill_brewer(type = "qual", palette =  "Paired") +
  ggtitle(paste0("Preventable cases (incident cases)")) + 
  theme(legend.title = element_blank())

ggsave2(filename = output_dir(paste0("PARFcases_inc_all.png")), scale = 0.8)


#Can't get the labeller function to work so being lazy
labnames2 <- c("Prostate Cancer",  "Colorectal Cancer", "Lung Cancer","Breast Cancer",  "Dementia", "COPD",    "Atrial Fibrillation" , 
               "Asthma"  , "Stroke"  ,   "CHD"      ,    "T2 Diabetes"      ,   "CKD" )

inc_tab_sum[, disease2 := disease]
inc_tab_sum[, disease2 := factor(disease2,
                                  levels = c("prostate_ca", "colorect_ca", "lung_ca",  "breast_ca", "dementia", "copd", "af" ,    
                                             "asthma"   , "stroke"  ,    "chd"      ,    "t2dm"      ,   "ckd" ),
                                  labels = labnames2)]


ggcust(inc_tab_sum[, .(cases = sum(casenum)), keyby = .(disease2, type, dimd)], 
       aes(x = dimd, y = cases, fill =type)) +
  facet_wrap( vars(disease2)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = "IMD decile", guide = guide_axis(angle = 45 ) ) +
  scale_fill_brewer(type = "qual", palette =  "Paired") +
  ggtitle(paste0("Preventable cases by condition and IMD\n(incident cases)")) + 
  theme(legend.title = element_blank())
ggsave2(filename = output_dir(paste0("PARFcases_inc_imd_all.png")), scale = 0.8)



#Stacked with risk factor attributions 
#Renaming risk factors for graphs 
riskfactors <- c("overall", "physical activity", "alcohol",  "BMI",  "fruit",  
                 "SBP", "total cholesterol",  "veg" , "smoking" , "ETS")
#setnames(parf, 
#         c("parf", "parf_pa", "parf_alc",  "parf_bmi",  "parf_fru",  
#           "parf_sbp", "parf_tch",  "parf_veg" , "parf_smo" , "parf_ets"),
#         riskfactors)


inc_tab2 <- inc_tab[, .(disease,age, sex,dimd, cases, wt)]

parfbyrisk_inc <- data.table()
for (i in riskfactors[riskfactors != "overall"]){
  parfdata <- parf[, weighted.mean(get(i)* parf_correction, pop_size), keyby = .(disease,  age,sex,dimd, mc)][
    , .(parf = quantile(V1,0.5, na.rm = T)), keyby = .(disease,  age,sex,dimd)][, riskfactor := paste0(i)]
  subtab <- copy(inc_tab2)
  subtab[parfdata, on = c("disease","age", "sex", "dimd"), parf := i.parf]
  subtab <- subtab[!is.na(parf), ]
  subtab[, attrib := cases * parf][, riskfactor := paste0(i)]
  parfbyrisk_inc <- rbind(parfbyrisk_inc, subtab)
}




parfbyrisksum_inc <- parfbyrisk_inc[, .( attrib  = sum(attrib*wt)),
                            keyby = .(disease,riskfactor, dimd)]
parfbyrisksum_inc <- rbind(parfbyrisksum_inc, inc_tab_sum[type == "Not preventable"  , .(disease, dimd, attrib = casenum, riskfactor = "Not preventable")])
parfbyrisksum_inc[, riskfactor:= factor(riskfactor,
                                    levels = c("Not preventable", "alcohol",  "BMI",  "ETS", "fruit",  "physical activity",
                                               "SBP", "smoking" , "total cholesterol",  "veg"   ),
                                    labels = c("Not preventable", "Alcohol",  "BMI",  "ETS", "Fruit",  "Physical activity",
                                               "SBP", "Smoking" , "Total cholesterol",  "Veg"  ))]

parfbyrisksum_inc[, disease := factor(disease,
                                  levels =  c("prostate_ca", "colorect_ca", "lung_ca",  "breast_ca", "dementia", "copd", "af" ,    
                                              "asthma"   , "stroke"  ,    "chd"      ,    "t2dm"      ,   "ckd" ),
                                  labels = labnames2)]


ggcust(parfbyrisksum_inc[, .(cases = sum(attrib)), keyby = .(disease, riskfactor)], 
       aes(x = disease, y = cases, fill =riskfactor)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = "Condition", guide = guide_axis(angle = 45 ), labels = labnames ) +
  scale_fill_brewer(type = "qual", palette = "Paired") + 
  labs(fill = "Preventable\nrisk factor") +
  ggtitle(paste0("Preventable cases by risk factor\n(incident cases)")) 
ggsave2(filename = output_dir(paste0("PARFcases_inc_riskfactor.png")), scale = 0.8)

ggcust(parfbyrisksum_inc[, .(cases = sum(attrib)), keyby = .(disease, riskfactor, dimd)], 
       aes(x = dimd, y = cases, fill =riskfactor)) +
  geom_col(position = "stack") +
  facet_wrap( vars(disease)) +
  scale_x_discrete(name = "IMD decile", guide = guide_axis(angle = 45 ) ) +
  scale_fill_brewer(type = "qual", palette = "Paired") + 
  labs(fill = "Preventable\nrisk factor") +
  ggtitle(paste0("Preventable cases by risk factor and IMD\n(incident cases)")) 
ggsave2(filename = output_dir(paste0("PARFcases_inc_riskfactor_imd.png")), scale = 0.8)












# Age standardised incident PARFs - THESE ARE COMPLETELY WORNG 
#Downloaded the 2013 european pop from here https://www.opendata.nhs.scot/en/dataset/standard-populations/resource/29ce4cda-a831-40f4-af24-636196e05c1a 
standpop <- fread("/mnt/alhead/UoL/CPRD2021/epi_models/scripts/PARF outputs/case_numbers/european_standard_population_by_sex.csv")

agegroup5yfn <- function(x) { 
  lev <- c("20-24 years", "25-29 years", "30-34 years",
           "35-39 years", "40-44 years", "45-49 years", "50-54 years", "55-59 years",
           "60-64 years", "65-69 years", "70-74 years", "75-79 years", "80-84 years",
           "85-89 years", "90plus years")
  x[, agegrp5 := fcase(
    (age) >= 90L,           factor(lev[15], levels = lev),
    between(age, 85, 89), factor(lev[14], levels = lev),
    between(age, 80, 84), factor(lev[13], levels = lev),
    between(age, 75, 79), factor(lev[12], levels = lev),
    between(age, 70, 74), factor(lev[11], levels = lev),
    between(age, 65, 69), factor(lev[10], levels = lev),
    between(age, 60, 64), factor(lev[9], levels = lev),
    between(age, 55, 59), factor(lev[8], levels = lev),
    between(age, 50, 54), factor(lev[7], levels = lev),
    between(age, 45, 49), factor(lev[6], levels = lev),
    between(age, 40, 44), factor(lev[5], levels = lev),
    between(age, 35, 39), factor(lev[4], levels = lev),
    between(age, 30, 34), factor(lev[3], levels = lev),
    between(age, 25, 29), factor(lev[2], levels = lev),
    between(age, 20, 24), factor(lev[1], levels = lev)
  )
  ]
}
agegroup5yfn(cprdpop)
cprdpop_stand <- cprdpop[, .(N = sum(N)), keyby = .(dimd, sex, agegrp5)]
cprdpop_stand[standpop, on = c("sex", "agegrp5"), standpop := i.EuropeanStandardPopulation ]
cprdpop_stand[, wt := standpop/N]
cprdpop_stand[, wt := wt*sum(standpop)/sum(wt, na.rm = TRUE), by = .(sex, dimd)]
cprdpop_stand[, sum(wt), keyby = .(sex, dimd)] #this adds up to the standard pop, which is 78500 (not 100,000) because missing some groups 
cprdpop_stand[, wt2:= wt/78500] #making weights add up to1 

## Incidence
dt_inc <- harmonise(read_fst(input_path("panel_short_inc.fst"),
                             as.data.table = TRUE)[
                               gender != "I" & year == 2013])[
                                 between(age, 20, 100), .SD, .SDcols = c(disnm2, strata)]

dt_inc <- na.omit(dt_inc )
setnames(dt_inc, disnm2, disnm)


dimdlabs <- c("1 most deprived" ,  "2" ,   "3","4","5" ,  "6"  , "7" ,  "8" , "9" , "10 least deprived")
dt_inc[, dimd := factor(dimd,
                        levels = 10:1,
                        labels = dimdlabs)]

agegroup5yfn(dt_inc)
dt_inc[cprdpop_stand, on = c("sex", "dimd", "agegrp5"), wt := i.wt2]

inc_tab <- data.table()
for (i in disnm){
  subtab <- dt_inc[, .(cases = sum(ifelse(get(i) == 1L, wt, 0)),
                       atrisk = sum(ifelse(get(i) != 2L, wt, 0))), 
                   keyby = .(agegrp5, sex, dimd) ][
                     , disease := paste0(i)]
  inc_tab <- rbind(inc_tab,subtab)
  inc_tab  
}

inc_tab[sex == "men" & disease == "breast_ca", atrisk := 0]
inc_tab[sex == "women" & disease == "prostate_ca", atrisk := 0]




inc_tab[, incrate := cases/atrisk * 100000]
inc_tab[, .(incrate = cases/atrisk * 100000), keyby = .(disease,sex, dimd)]
#I think these are standardised incident rates


#BUT, how do I apply the parfs to these? 
parf <- fread("/mnt/storage_fast/output/hf_real_parf_24May/parf/parf_final_corrected.csv") 

agegroup5yfn(parf)        
parf[cprdpop_stand, on = c("sex", "dimd", "agegrp5"), wt := i.wt2]


parfsum <- parf[, weighted.mean(parf, wt), keyby = .(disease,  agegrp5,sex,dimd, mc)][
  , .(parf = quantile(V1,0.5)), keyby = .(disease,  agegrp5,sex,dimd)]
parfsum

        


inc_tab[parfsum, on = c("disease","agegrp5", "sex", "dimd"), parf := i.parf]
inc_tab <- inc_tab[!is.na(parf), ]
inc_tab[, `:=` (attrib = cases       * parf,
                nonattrib = cases       * (1-parf))]


inc_tab_sum <- inc_tab[, .(totcases = sum(cases),
                           atrisk  = sum(atrisk),
                           attribcases  = sum(attrib),
                           nonattribcases = sum(nonattrib)),
                       keyby = .(disease,dimd)]
inc_tab_sum <- melt(inc_tab_sum, 
                    measure.vars = c("attribcases", "nonattribcases"), 
                    value.name = "cases",
                    variable.name = "type")
inc_tab_sum[, type:= factor(type,
                            levels = c("nonattribcases", "attribcases"),
                            labels = c("Not preventable", "Preventable"))]

inc_tab_sum[, disease := factor(disease,
                                levels = c("prostate_ca", "colorect_ca", "lung_ca",  "asthma"   ,  
                                "copd", "dementia",  "breast_ca", "af" ,    "stroke"  ,       "t2dm"      ,"chd"      ,    "ckd" ),
                                labels = c("prostate_ca", "colorect_ca", "lung_ca",  "asthma"   ,  
                                           "copd", "dementia",  "breast_ca", "af" ,    "stroke"  ,       "t2dm"      ,"chd"      ,    "ckd" ))]


labnames <- c("Prostate\nCancer", "Colorectal\nCancer", "Lung\nCancer",  "Asthma"  , 
              "COPD", "Dementia",   "Breast\nCancer",    "Atrial\nFibrillation" ,
             "Stroke"    ,    "T2\nDiabetes"      , "CHD"      ,  "CKD" )

ggcust(inc_tab_sum[, .(rate = sum(cases)/sum(atrisk)*100000), keyby = .(disease, type)], 
       aes(x = disease, y = rate, fill =type)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = "Condition", guide = guide_axis(angle = 45 ), labels = labnames ) + 
  scale_y_continuous(name = "Incidence proportion (per 100,000)", labels = comma_format()) +
  scale_fill_brewer(type = "qual", palette =  "Paired") +
  ggtitle(paste0("Preventable cases (incidence proportion, standardised)")) + 
  theme(legend.title = element_blank())

ggsave2(filename = output_dir(paste0("PARFcases_inc_stand_all.png")), scale = 0.8)


#Can't get the labeller function to work so being lazy
labnames2 <- c("Prostate Cancer",  "Colorectal Cancer", "Lung Cancer", "Asthma"  , "COPD",   "Dementia", "Breast Cancer",  "Atrial Fibrillation" , 
               "Stroke"  ,   "T2 Diabetes"      ,  "CHD"      ,     "CKD" )

inc_tab_sum[, disease2 := disease]
inc_tab_sum[, disease2 := factor(disease2,
                                 levels = c("prostate_ca", "colorect_ca", "lung_ca",  "asthma"   ,  
                                            "copd", "dementia",   "breast_ca", "af" ,    "stroke"  ,       "t2dm"      ,"chd"      ,    "ckd"),
                                 labels = labnames2)]


ggcust(inc_tab_sum[, .(rate =sum(cases)/sum(atrisk)*100000), keyby = .(disease2, type, dimd)], 
       aes(x = dimd, y = rate, fill =type)) +
  facet_wrap( vars(disease2)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = "IMD decile", guide = guide_axis(angle = 45 ) ) +
  scale_y_continuous(name = "Incidence proportion (per 100,000)", labels = comma_format()) +
  scale_fill_brewer(type = "qual", palette =  "Paired") +
  ggtitle(paste0("Preventable cases by condition and IMD\n(incidence proportion, standardised)")) + 
  theme(legend.title = element_blank())
ggsave2(filename = output_dir(paste0("PARFcases_inc_stand_imd_all.png")), scale = 0.8)



#Stacked with risk factor attributions 
#Renaming risk factors for graphs 
riskfactors <- c("overall", "physical activity", "alcohol",  "BMI",  "fruit",  
                 "SBP", "total cholesterol",  "veg" , "smoking" , "ETS")
#setnames(parf, 
#         c("parf", "parf_pa", "parf_alc",  "parf_bmi",  "parf_fru",  
#           "parf_sbp", "parf_tch",  "parf_veg" , "parf_smo" , "parf_ets"),
#         riskfactors)


inc_tab2 <- inc_tab[, .(disease,agegrp5, sex,dimd, cases, atrisk)]

parf <- fread("/mnt/storage_fast/output/hf_real_parf_24May/parf/parf_final_corrected.csv") 
agegroup5yfn(parf)
setnames(parf, 
         c("parf", "parf_pa", "parf_alc",  "parf_bmi",  "parf_fru",  
           "parf_sbp", "parf_tch",  "parf_veg" , "parf_smo" , "parf_ets"),
         riskfactors)
parf[cprdpop_stand, on = c("sex", "dimd", "agegrp5"), wt := i.wt2]



parfbyrisk_inc <- data.table()
for (i in riskfactors[riskfactors != "overall"]){
  parfdata <- parf[, weighted.mean(get(i)* parf_correction, wt), keyby = .(disease,  agegrp5,sex,dimd, mc)][
    , .(parf = quantile(V1,0.5, na.rm = T)), keyby = .(disease,  agegrp5,sex,dimd)][, riskfactor := paste0(i)]
  subtab <- copy(inc_tab2)
  subtab[parfdata, on = c("disease","agegrp5", "sex", "dimd"), parf := i.parf]
  subtab <- subtab[!is.na(parf), ]
  subtab[, attrib := cases * parf][, riskfactor := paste0(i)]
  parfbyrisk_inc <- rbind(parfbyrisk_inc, subtab)
}




parfbyrisksum_inc <- parfbyrisk_inc[, .( attrib  = sum(attrib), atrisk = sum(atrisk)),
                                    keyby = .(disease,riskfactor, dimd)]
parfbyrisksum_inc <- rbind(parfbyrisksum_inc, inc_tab_sum[type == "Not preventable"  , .(disease, dimd, attrib = cases,atrisk, riskfactor = "Not preventable")])
parfbyrisksum_inc[, riskfactor:= factor(riskfactor,
                                        levels = c("Not preventable", "alcohol",  "BMI",  "ETS", "fruit",  "physical activity",
                                                   "SBP", "smoking" , "total cholesterol",  "veg"   ),
                                        labels = c("Not preventable", "Alcohol",  "BMI",  "ETS", "Fruit",  "Physical activity",
                                                   "SBP", "Smoking" , "Total cholesterol",  "Veg"  ))]

parfbyrisksum_inc[, disease2 := factor(disease,
                                      levels =  c("prostate_ca", "colorect_ca", "lung_ca",   "asthma"   ,  
                                                  "copd", "dementia","breast_ca",  "af" ,    "stroke"  ,       "t2dm"      ,"chd"      ,    "ckd" ),
                                      labels = labnames2)]


ggcust(parfbyrisksum_inc[, .(rate = sum(attrib)/sum(atrisk)*100000), keyby = .(disease2, riskfactor)], 
       aes(x = disease2, y = rate, fill =riskfactor)) +
  geom_col(position = "stack") +
  scale_x_discrete(name = "Condition", guide = guide_axis(angle = 45 ), labels = labnames ) +
  scale_y_continuous(name = "Incidence proportion (per 100,000)", labels = comma_format()) +
  scale_fill_brewer(type = "qual", palette = "Paired") + 
  labs(fill = "Preventable\nrisk factor") +
  ggtitle(paste0("Preventable cases by risk factor\n(incidence proportion, standardised)")) 
ggsave2(filename = output_dir(paste0("PARFcases_inc_stand_riskfactor.png")), scale = 0.8)

ggcust(parfbyrisksum_inc[, .(rate = sum(attrib)/sum(atrisk)*100000), keyby = .(disease2, riskfactor, dimd)], 
       aes(x = dimd, y = rate, fill =riskfactor)) +
  geom_col(position = "stack") +
  facet_wrap( vars(disease2)) +
  scale_x_discrete(name = "IMD decile", guide = guide_axis(angle = 45 ) ) +
  scale_y_continuous(name = "Incidence proportion (per 100,000)", labels = comma_format()) +
  scale_fill_brewer(type = "qual", palette = "Paired") + 
  labs(fill = "Preventable\nrisk factor") +
  ggtitle(paste0("Preventable cases by risk factor and IMD\n(incidence proportion, standardised)")) 
ggsave2(filename = output_dir(paste0("PARFcases_inc_stand_riskfactor_imd.png")), scale = 0.8)






