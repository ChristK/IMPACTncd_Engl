library(data.table)
library(fst)
tt <- fread("/home/ckyprid/My_Models/IMPACTncd_Engl/secure_data/HSE_data/hse11ai.tab")

# Calculate the impact of cardiovascular disease (cvdis) and diabetes on various health outcomes
tt[omsysval< 50, omsysval := NA]
tt[cholval1 < 1, cholval1 := NA]
tt[bmival < 10, bmival := NA]
tt[ag16g10 < 1, ag16g10 := NA]
tt[cvdis < 1, cvdis := NA]
tt[cvdis == 2, cvdis := 0]
tt[diabete2r < 1, diabete2r := NA]
tt[diabete2r == 2, diabete2r := 0]
summary(glm(omsysval ~ factor(ag16g10) + Sex + factor(qimd) + factor(cvdis), data = tt, family = gaussian(link = "identity"), weights = wt_nurse))
# -3.6 mmHg
summary(glm(cholval1 ~ factor(ag16g10) + Sex + factor(qimd) + factor(cvdis), data = tt, family = gaussian(link = "identity"), weights = wt_blood))
# -0.78 mmol/L = -30 mg/DL
summary(glm(bmival ~ factor(ag16g10) + Sex + factor(qimd) + factor(cvdis), data = tt, family = gaussian(link = "identity"), weights = wt_int))
# +0.9 kg/m^2

summary(glm(omsysval ~ factor(ag16g10) + Sex + factor(qimd) + factor(diabete2r), data = tt, family = gaussian(link = "identity"), weights = wt_nurse))
# not significant
summary(glm(cholval1 ~ factor(ag16g10) + Sex + factor(qimd) + factor(diabete2r), data = tt, family = gaussian(link = "identity"), weights = wt_blood))
# - 1.0 mmol/L = - 38.67 mg/DL
summary(glm(bmival ~ factor(ag16g10) + Sex + factor(qimd) + factor(diabete2r), data = tt, family = gaussian(link = "identity"), weights = wt_int))
# +3.7 kg/m^2
