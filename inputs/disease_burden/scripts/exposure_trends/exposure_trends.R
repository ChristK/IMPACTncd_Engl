# Trends in exposure

library(data.table) # for fast data manipulation
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(cowplot)


# folder for plots
output_dir <-
  function(x = character(0))
    paste0("/mnt/", Sys.info()[["user"]],
           "/UoL/CPRD2021/epi_models/scripts/exposure_trends/", x)


myDT <- fread("/mnt/storage_fast/output/hf_real/xps/xps20.csv.gz") #Ignore hypertension


#Changing sex names to have capitals for plots
myDT[, sex := factor(sex,
                     levels = c("All", "men", "women"),
                     labels = c("All", "Men", "Women"))]


#Agegrp20 as ordered factor
myDT[, agegrp20 := factor(agegrp20,
                      levels = c("All", "30-49", "50-69", "70-89","90+"  ),
                      labels = c("All", "30-49", "50-69", "70-89","90+"  )
)]

#DIMD as ordered factor
myDT[, qimd := factor(qimd,
                      levels = c("All","1 most deprived", "2", "3", "4", "5 least deprived"),
                      labels = c("All", "1 most deprived", "2", "3", "4", "5 least deprived")
)]

#Ethnicity as a factor
myDT[, ethnicity := factor(ethnicity,
                           levels = c("All","bangladeshi", "black african",
                                      "black caribbean", "chinese", "indian",
                                      "other", "other asian", "pakistani",
                                      "white"),
                           labels =  c("All","Bangladeshi", "Black african",
                                       "Black caribbean", "Chinese", "Indian",
                                       "Other", "Other asian", "Pakistani",
                                       "White"))]

#Region as a factor
myDT[, sha := factor(sha,
                           levels = c("All","East Midlands", "East of England", "London", "North East",
                                      "North West", "South Central", "South East Coast", "South West",
                                      "West Midlands", "Yorkshire and the Humber"),
                           labels =  c("All", "East Midlands", "East of England", "London", "North East",
                                                "North West", "South Central", "South East Coast", "South West",
                                                "West Midlands", "Yorkshire + Humber"))]



#Names for graphs
exp_nm <- c("active_days_curr_xps", "fruit_curr_xps",  "veg_curr_xps",  "alcohol_curr_xps",
            "bmi_curr_xps" , "sbp_curr_xps" , "tchol_curr_xps",  "statin_px_curr_xps"   ,
            "smok_never_curr_xps"  ,  "smok_active_curr_xps" )

exp_nm_title <- c("Mean days active", "Mean fruit consumption",  "Mean veg consumption",  "Mean alcohol consumption",
            "Mean BMI" , "Mean SBP" , "Mean total cholesterol",  "Proportion on statins"   ,
            "Proportion of never smokers"  ,  "Proportion of active smokers" )

exp_nm_y <- c("Days active", "Fruit consumption (g)",  "Veg consumption (g)",  "Alcohol consumption (g)",
                  "BMI" , "SBP" , "Total cholesterol",  "Proportion on statins"   ,
                  "Proportion of never smokers"  ,  "Proportion of active smokers" )

ylowlim <- c(0,0,0,0,18,50,0,0,0,0)

lookuptab <- data.table(exp = exp_nm,
                        title = exp_nm_title,
                        yaxis = exp_nm_y,
                        ylowlim = ylowlim
                        )

#Setting the plot settings
theme_set(new = theme_few())
theme_update(axis.text.x = element_text(size = 9), plot.title = element_text(hjust = 0.5))




# By sex
for (i in exp_nm) {
tabdata <- myDT[sex != "All" & agegrp20 == "All" & qimd == "All" & ethnicity == "All" & sha == "All", mean(get(i)), keyby = .(sex, year, mc)][
  , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(year, sex)]
tabdata

title <- lookuptab[exp == i, title]
yaxis <- lookuptab[exp == i, yaxis]
ylowlim <- lookuptab[exp == i, ylowlim]

ggplot(tabdata, aes(x = year, y = `50%`, col = sex)) +
  geom_smooth(se = F) +
  scale_x_continuous(name = "Year") +
  scale_y_continuous(name = yaxis) +
  expand_limits(y = ylowlim) +
  ggtitle(paste0(title, " by sex")) +
  scale_color_brewer(type = "qual") +
  labs(col = "Sex")

ggsave(filename = output_dir(paste0("exp_sex_",i,".png")))
}
rm(i)


# By age group
for (i in exp_nm) {
  tabdata <- myDT[sex == "All" & agegrp20 != "All" & qimd == "All" & ethnicity == "All" & sha == "All", mean(get(i)), keyby = .(agegrp20, year, mc)][
    , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(year, agegrp20)]
  tabdata

  title <- lookuptab[exp == i, title]
  yaxis <- lookuptab[exp == i, yaxis]
  ylowlim <- lookuptab[exp == i, ylowlim]


  ggplot(tabdata, aes(x = year, y = `50%`, col = agegrp20)) +
    geom_smooth(se = F) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = yaxis) +
    expand_limits(y = ylowlim) +
    ggtitle(paste0(title, " by age group")) +
    scale_color_brewer(type = "qual") +
    labs(col = "20-year\nage group")

  ggsave(filename = output_dir(paste0("exp_agegro_",i,".png")))
}
rm(i)

# By qimd
for (i in exp_nm) {
  tabdata <- myDT[sex == "All" & agegrp20 == "All" & qimd != "All" & ethnicity == "All" & sha == "All", mean(get(i)), keyby = .(qimd, year, mc)][
    , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(year, qimd)]
  tabdata

  title <- lookuptab[exp == i, title]
  yaxis <- lookuptab[exp == i, yaxis]
  ylowlim <- lookuptab[exp == i, ylowlim]


  ggplot(tabdata, aes(x = year, y = `50%`, col = qimd)) +
    geom_smooth(se = F) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = yaxis) +
    expand_limits(y = ylowlim) +
    ggtitle(paste0(title, " by IMD")) +
    scale_color_viridis(discrete = T) +
    labs(col = "IMD quintile")

  ggsave(filename = output_dir(paste0("exp_qimd_",i,".png")))
}
rm(i)

# By ethnicity
for (i in exp_nm) {
  tabdata <- myDT[sex == "All" & agegrp20 == "All" & qimd == "All" & ethnicity != "All" & sha == "All", mean(get(i)), keyby = .(ethnicity, year, mc)][
    , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(year, ethnicity)]
  tabdata

  title <- lookuptab[exp == i, title]
  yaxis <- lookuptab[exp == i, yaxis]
  ylowlim <- lookuptab[exp == i, ylowlim]


  ggplot(tabdata, aes(x = year, y = `50%`, col = ethnicity)) +
    geom_smooth(se = F) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = yaxis) +
    expand_limits(y = ylowlim) +
    ggtitle(paste0(title, " by ethnicity")) +
    scale_colour_brewer(type = "qual", palette = "Set1") +
    labs(col = "Ethnicity")

  ggsave(filename = output_dir(paste0("exp_ethn_",i,".png")))
}
rm(i)


# By sha
for (i in exp_nm) {
  tabdata <- myDT[sex == "All" & agegrp20 == "All" & qimd == "All" & ethnicity == "All" & sha != "All", mean(get(i)), keyby = .(sha, year, mc)][
    , as.list(quantile(V1, c(0.5, 0.025, 0.975))), keyby = .(year, sha)]
  tabdata

  title <- lookuptab[exp == i, title]
  yaxis <- lookuptab[exp == i, yaxis]
  ylowlim <- lookuptab[exp == i, ylowlim]


  ggplot(tabdata, aes(x = year, y = `50%`, col = sha)) +
    geom_smooth(se = F) +
    scale_x_continuous(name = "Year") +
    scale_y_continuous(name = yaxis) +
    expand_limits(y = ylowlim) +
    ggtitle(paste0(title, " by region")) +
    scale_colour_brewer(type = "qual", palette = "Set3") +
    labs(col = "Region")

  ggsave(filename = output_dir(paste0("exp_sha_",i,".png")))
}
rm(i)
