# PARF outputs 

library(data.table) # for fast data manipulation
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(scales)
library(cowplot)



# folder for plots 
output_dir <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2021/epi_models/scripts/PARF outputs/", x)


parf <- fread("/mnt/storage_fast/output/hf_real_parf/parf/parf_final.csv") #Ignore hypertension


#20year age groups 
parf[, agegrp20 := cut_width(age, width = 20, boundary = 29)]
agegrps <- levels(parf$agegrp20)
parf[, agegrp20 := factor(agegrp20, 
                          levels = agegrps,
                          labels = c("30-49", "50-69", "70-89", "90+"))]
rm(agegrps)

#Changing sex names to have capitals for plots
parf[, sex := factor(sex, 
                     levels = c("men", "women"),
                     labels = c("Men", "Women"))]

#DIMD as ordered factor
parf[, dimd := factor(dimd, 
                      levels = c("1 most deprived", "2", "3", "4", "5",
                                 "6", "7", "8", "9", "10 least deprived"),
                      labels = c("1 most deprived", "2", "3", "4", "5",
                                 "6", "7", "8", "9", "10 least deprived")
                      )]
#Ethnicity as a factor
parf[, ethnicity := factor(ethnicity,
                           levels = c("bangladeshi", "black african", 
                                      "black caribbean", "chinese", "indian",
                                      "other", "other asian", "pakistani",
                                      "white"),
                           labels =  c("Bangladeshi", "Black african", 
                                       "Black caribbean", "Chinese", "Indian",
                                       "Other", "Other asian", "Pakistani",
                                       "White"))]

#Renaming risk factors for graphs 
riskfactors <- c("overall", "physical activity", "alcohol",  "BMI",  "fruit",  
                 "SBP", "total cholesterol",  "veg" , "smoking" , "ETS")
setnames(parf, 
         c("parf", "parf_pa", "parf_alc",  "parf_bmi",  "parf_fru",  
           "parf_sbp", "parf_tch",  "parf_veg" , "parf_smo" , "parf_ets"),
         riskfactors)


#Names for graphs 
disnm <- c("af", "asthma", "breast_ca", "chd", "ckd", "colorectal_ca", "copd",
           "dementia", "lung_ca", "prostate_ca" , "stroke" , "t2dm" )

disnm_long <- c("Atrial\nFibrillation", "Asthma", "Breast\nCancer", "CHD", "CKD", 
                "Colorectal\nCancer", "COPD", "Dementia", "Lung\nCancer", 
                "Prostate\nCancer" , "Stroke" , "T2\nDiabetes" )
imd_long <- c("1\nmost deprived", "2", "3", "4", "5",
              "6", "7", "8", "9", "10\nleast deprived")
imd_short <- c("1", "2", "3", "4", "5",
              "6", "7", "8", "9", "10")
ethnicity <-  c("Bangladeshi", "Black african", "Black caribbean", "Chinese", "Indian",
"Other", "Other asian", "Pakistani", "White")
region_short <- c("East Midlands", "East of England", "London", "North East",
                  "North West", "South Central", "South East Coast", "South West", 
                  "West Midlands", "Yorkshire + Humber")

  
#And for facet grids 
disnm_long2 <- c("Atrial Fibrillation", "Asthma", "Breast Cancer", "CHD", "CKD", 
                "Colorectal Cancer", "COPD", "Dementia", "Lung Cancer", 
                "Prostate Cancer" , "Stroke" , "T2 Diabetes" )
names(disnm_long2) <- disnm


#Setting the plot settings
theme_set(new = theme_few())
theme_update(plot.title = element_text(hjust = 0.5))
ggcust <- function(...){
  ggplot(...) +
    scale_fill_brewer(type = "qual") +
    scale_y_continuous(name = "PARF", labels = percent_format())    
}



#By sex
for (i in riskfactors) {
tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(sex, disease, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, sex)]
tabdata

ggcust(tabdata, 
  aes(x = disease, y = `50%`, fill = sex)) +
  geom_col(position = "dodge") +
    labs(fill = "Sex") +
  scale_x_discrete(name = "Condition", labels = disnm_long, guide = guide_axis(angle = 45 ) ) +
  ggtitle(paste0("Population Attributable Risk Fraction by Sex - ", i)) 
ggsave2(filename = output_dir(paste0("PARF_sex_all_",i,".png")), scale = 0.8)
}


#By IMD
for (i in riskfactors) {
  tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(dimd, disease, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, dimd)]

tabdata

ggcust(tabdata, 
       aes(x = dimd, y = `50%`, fill = dimd)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(disease), labeller = labeller(disease = disnm_long2)) + 
  scale_x_discrete(name = "IMD Deciles (1 most deprived)", labels = imd_short) +
  ggtitle(paste0("Population Attributable Risk Fraction by IMD deciles - ", i)) + 
  scale_fill_viridis(discrete = T) + 
  theme(legend.position="none")
ggsave2(filename = output_dir(paste0("PARF_imd_all_",i,".png")), scale = 0.8)
}

#By Ethnicity
for (i in riskfactors) {
  tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(ethnicity, disease, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, ethnicity)]

tabdata

ggcust(tabdata, 
       aes(x = ethnicity, y = `50%`, fill = ethnicity)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(disease), labeller = labeller(disease = disnm_long2)) + 
  scale_x_discrete(name = "Ethnicity", labels = ethnicity, 
                   guide = guide_axis(angle = 45 )) +     
  scale_fill_brewer(type = "qual", palette = "Set1") +
  theme(legend.position="none") +
  ggtitle(paste0("Population Attributable Risk Fraction by Ethnicity - ", i)) 

ggsave2(filename = output_dir(paste0("PARF_eth_all_",i,".png")), scale = 0.8)
}

#By Region
for (i in riskfactors) {
tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(sha, disease, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, sha)]
tabdata

ggcust(tabdata, 
       aes(x = sha, y = `50%`, fill = sha)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(disease), labeller = labeller(disease = disnm_long2)) + 
  scale_x_discrete(name = "Region", labels = region_short, 
                   guide = guide_axis(angle = 45 )) +
  ggtitle(paste0("Population Attributable Risk Fraction by Region - ",i )) +     
  scale_fill_brewer(type = "qual", palette = "Paired") +
  theme(legend.position="none")
ggsave2(filename = output_dir(paste0("PARF_region_all_",i,".png")), scale = 0.8)
}

#By Agegroup
for (i in riskfactors) {
tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(agegrp20, disease, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, agegrp20)]
tabdata

ggcust(tabdata, 
       aes(x = disease, y = `50%`, fill = agegrp20)) +
  geom_col(position = "dodge") +
  ggtitle(paste0("Population Attributable Risk Fraction\nby 20-year age-group - ", i )) +
  labs(fill = "20-year\nage-group") +
  scale_x_discrete(name = "Condition", labels = disnm_long,
                   guide = guide_axis(angle = 45 )) 
ggsave2(filename = output_dir(paste0("PARF_age_all_",i,".png")), scale = 0.8)
}

## Stratified by sex + xxx

# + IMD
for (i in riskfactors) {

tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(dimd, disease, sex, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, sex, dimd)]

tabdata
ggcust(tabdata, 
  aes(x = dimd, y = `50%`, fill = sex)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(disease), labeller = labeller(disease = disnm_long2)) + 
  labs(fill = "Sex") +
  ggtitle(paste0("Population Attributable Risk Fraction\nby Sex & IMD - ", i)) +
  scale_x_discrete(name = "Deciles of IMD (1 most deprived)", labels = imd_short ) 
ggsave2(filename = output_dir(paste0("PARF_sex_imd_all_",i,".png")), scale = 0.8)
}



# + ethnicity
for (i in riskfactors) {
  
tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(ethnicity, disease, sex, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, sex, ethnicity)]

tabdata

ggcust(tabdata, 
  aes(x = ethnicity, y = `50%`, fill = sex)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(disease), labeller = labeller(disease = disnm_long2)) + 
  labs(fill = "Sex") +
  ggtitle(paste0("Population Attributable Risk Fraction\nby Sex & Ethnicity - ", i)) +
  scale_x_discrete(name = "Ethnicity", labels = ethnicity,
                   guide = guide_axis(angle = 45 ) )
ggsave2(filename = output_dir(paste0("PARF_sex_eth_all_",i,".png")), scale = 0.8)

}

# + region
for (i in riskfactors) {
  
tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(sha, disease, sex, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, sex, sha)]

tabdata

ggcust(tabdata, 
       aes(x = sha, y = `50%`, fill = sex)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(disease), labeller = labeller(disease = disnm_long2)) + 
  labs(fill = "Sex") +
  ggtitle(paste0("Population Attributable Risk Fraction\nby Sex & Region - ", i)) +
  scale_x_discrete(name = "Region",
                   labels = region_short,
                   guide = guide_axis(angle = 45 ) )
ggsave2(filename = output_dir(paste0("PARF_sex_region_all_",i,".png")), scale = 0.8)

}

# + age group
for (i in riskfactors) {
  
tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(agegrp20, disease, sex, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, sex, agegrp20)]

tabdata

ggcust(tabdata, 
       aes(x = agegrp20, y = `50%`, fill = sex)) +
  geom_col(position = "dodge") +
  facet_wrap(vars(disease), labeller = labeller(disease = disnm_long2)) + 
  labs(fill = "Sex") +
  ggtitle(paste0("Population Attributable Risk Fraction\nby Sex & 20-year age-group - ", i)) +
  scale_x_discrete(name = "20-year age-group",
                   guide = guide_axis(angle = 45 ) )
ggsave2(filename = output_dir(paste0("PARF_sex_age_all_", i,".png")), scale = 0.8)
}




## By imd & ethnicity & sex
for (i in riskfactors) {
  
tabdata <- parf[, weighted.mean(get(i), pop_size), keyby = .(dimd, disease, ethnicity, sex, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(disease, ethnicity, dimd, sex)]
tabdata

#adding in long diseasenames for ease 
tabdata[, disease_long := factor(disease, levels = disnm, labels = disnm_long2)]

for (j in disnm_long2){
ggcust(tabdata[disease_long == j], 
       aes(x = dimd, y = `50%`, fill = sex)) + 
  geom_col(position = "dodge") + 
  facet_wrap(vars(ethnicity)) + 
  labs(fill = "Sex") +
  ggtitle(paste0("Population Attributable Risk Fraction for ", j, "\nby Sex, IMD & Ethnicity - ", i)) +
  scale_x_discrete(name = "Deciles of IMD (1 most deprived)", labels = imd_short ) 
  tmp <- tabdata[disease_long == j, unique(disease)]
  ggsave2(filename = output_dir(paste0("PARF_sex_imd_eth_", tmp ,"_",i, ".png")), scale = 0.8)
  rm(tmp)
  }
}



