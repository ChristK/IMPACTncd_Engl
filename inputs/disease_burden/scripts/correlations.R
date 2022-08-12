#Correlations between conditions 

library(data.table) # for fast data manipulation
library(fst) # Fast way to save and load data.tables
library(qs)

source(paste0("/mnt/", Sys.info()[["user"]], "/UoL/CPRD2021/epi_models/scripts/aux_fn.R"))


dpnd <- c("Alcohol problems", "Atrial Fibrillation" ,"CHD", "COPD"  ,      
          "Chronic Kidney Disease" , "Connective tissue disorder", "Dementia" ,                      
          "Diabetes excl Type 2", "Epilepsy", "Hearing loss", "Heart failure",  
          "IBS" , "Other cancers" , "Primary Malignancy_Breast",   "Primary Malignancy_Colorectal"  , "Primary Malignancy_Lung" ,       
          "Primary Malignancy_Prostate"  , "Rheumatoid Arthritis"   ,         "Stroke" ,                        
          "Type 2 Diabetes Mellitus"  ,      "Hypertension_camb"       ,        "Psychosis_camb"            , "Asthma - spell"       ,
          "Anxiety_Depression - spell"    ,  "Pain"   ,                         "Constipation") 


strata <- c("age", "dimd")

dt <- harmonise(read_fst(input_path("panel_short_prev.fst"),
                         as.data.table = TRUE)[gender != "I"])[
                           between(age, 20, 99) &
                             year == 2018,
                           .SD, 
                           .SDcols = c(dpnd, strata)]



for (j in dpnd) { #Need to make the exposures binary
  set(dt, NULL, j, fifelse(is.na(dt[[j]]), 0L, dt[[j]])) #setting those without spells to 0 
  set(dt, NULL, j, fifelse(dt[[j]] > 0L, 1L, 0L)) #want incidence & prev
}
rm(j)



library("Hmisc")
library("corrplot")

dt_mtrx <-as.matrix(dt[, 1: 26])
dt_rcorr <-rcorr(dt_mtrx, type = "spearman")
dt_rcorr_p <-round(dt_rcorr[["P"]], 3)
dt_rcorr_p
qsave(dt_rcorr, output_path("dependencies/corrmat.qs"), nthreads = 10)

corrplot(dt_rcorr$r, method = "circle")

imd <- unique(dt$dimd)
for (i in 1: length(imd)){
dtnew <- dt[age < 60 & dimd %in% imd[i]]
dt_mtrxnew <-as.matrix(dtnew[, 1: 26])
dt_rcorr <-rcorr(dt_mtrxnew, type = "spearman")
corrplot(dt_rcorr$r, method = "circle")
qsave(corrplot, output_path(paste0("dependencies/corrplot_under60_",imd[i],".qs")), nthreads = 10)
}



