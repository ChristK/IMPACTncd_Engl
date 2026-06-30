library(data.table) # for fast data manipulation
library(fst) # Fast way to save and load data.tables
library(gamlss)
library(qs)
library(foreach)
library(doParallel)

source(paste0("/mnt/", Sys.info()[["user"]], "/UoL/CPRD2021/epi_models/scripts/aux_fn.R"))

dt <- harmonise(read_fst(input_path("recovery_dur.fst"),
                         as.data.table = TRUE
)[gender != "I"])[between(age, 20, 100)
                  ][, dimd := factor(dimd,
                    levels = as.character(10:1),
                    labels = c("1 most deprived", 2:9, "10 least deprived"))]
dt[, ethnicity := NULL]
dt[, year := year(date) - 2000]
dt[, dur := as.integer(floor(dur/365) - 1)] # during the sim will add 1
dt[, hist(dur/365, 10)] # during the sim will add 1

# Asthma ----
dt[disease == "Asthma - spell", hist(dur)]
dt[disease == "Asthma - spell", summary(dur)]

# marg_distr <- gamlss::fitDist(
#   dt[disease == "Asthma - spell"]$dur,
#   log(nrow(dt)),
#   type = "counts",
#   try.gamlss = TRUE,
#   trace = TRUE
# )
# head(marg_distr$fits) # best distributions ZANBI

m1 <- gamlss(
  dur ~ log(year) + log(age) + sex + dimd,
  ~ log(year) + log(age) + sex + dimd,
  ~ log(year) + log(age) + sex + dimd,
  family = "ZANBI",
  data = dt[disease == "Asthma - spell"],
  method = mixed(20, 100)
)

m <- stepGAIC(m1, parallel = "multicore", ncpus = 16L)


cc <- list(
  "distr" = "ZANBI",
  "note" = "remember to add 1 to the final result",
  "mu" = list(
    "intercept" = m$mu.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$mu.coefficients["log(age)"])) 0 else m$mu.coefficients["log(age)"],
    "log(year)" = if (is.na(m$mu.coefficients["log(year)"])) 0 else m$mu.coefficients["log(year)"],
    "sex" = if (is.na(m$mu.coefficients["sexwomen"])) rep(0, 2) else c(0, m$mu.coefficients["sexwomen"]),
    "dimd" = if (length(m$mu.coefficients[grep("dimd", names(m$mu.coefficients))]) == 0) rep(0, 10) else c(0, m$mu.coefficients[grep("dimd", names(m$mu.coefficients))])
  ),
  "sigma" = list(
    "intercept" = m$sigma.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$sigma.coefficients["log(age)"])) 0 else m$sigma.coefficients["log(age)"],
    "log(year)" = if (is.na(m$sigma.coefficients["log(year)"])) 0 else m$sigma.coefficients["log(year)"],
    "sex" = if (is.na(m$sigma.coefficients["sexwomen"])) rep(0, 2) else c(0, m$sigma.coefficients["sexwomen"]),
    "dimd" = if (length(m$sigma.coefficients[grep("dimd", names(m$sigma.coefficients))]) == 0) rep(0, 10) else c(0, m$sigma.coefficients[grep("dimd", names(m$sigma.coefficients))])
  ),
  "nu" = list(
    "intercept" = m$nu.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$nu.coefficients["log(age)"])) 0 else m$nu.coefficients["log(age)"],
    "log(year)" = if (is.na(m$nu.coefficients["log(year)"])) 0 else m$nu.coefficients["log(year)"],
    "sex" = if (is.na(m$nu.coefficients["sexwomen"])) rep(0, 2) else c(0, m$nu.coefficients["sexwomen"]),
    "dimd" = if (length(m$nu.coefficients[grep("dimd", names(m$nu.coefficients))]) == 0) rep(0, 10) else c(0, m$nu.coefficients[grep("dimd", names(m$nu.coefficients))])
  )
)

yaml::write_yaml(cc, "~/My_Models/IMPACTncd_Engl/inputs/disease_burden/asthma_dur_forward2.yaml")
# yaml::write_yaml(cc, output_path("asthma_dur_forward2.yaml"))

# Andep ----
dt[disease == "Anxiety_Depression - spell", summary(dur)]
m1 <- gamlss(
  dur ~ log(year) + log(age) + sex + dimd,
  ~ log(year) + log(age) + sex + dimd,
  ~ log(year) + log(age) + sex + dimd,
  family = "ZANBI",
  data = dt[disease == "Anxiety_Depression - spell"],
  method = mixed(20, 100)
)

m <- stepGAIC(m1, parallel = "multicore", ncpus = 16L)



cc <- list(
  "distr" = "ZANBI",
  "note" = "remember to add 1 to the final result",
  "mu" = list(
    "intercept" = m$mu.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$mu.coefficients["log(age)"])) 0 else m$mu.coefficients["log(age)"],
    "log(year)" = if (is.na(m$mu.coefficients["log(year)"])) 0 else m$mu.coefficients["log(year)"],
    "sex" = if (is.na(m$mu.coefficients["sexwomen"])) rep(0, 2) else c(0, m$mu.coefficients["sexwomen"]),
    "dimd" = if (length(m$mu.coefficients[grep("dimd", names(m$mu.coefficients))]) == 0) rep(0, 10) else c(0, m$mu.coefficients[grep("dimd", names(m$mu.coefficients))])
  ),
  "sigma" = list(
    "intercept" = m$sigma.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$sigma.coefficients["log(age)"])) 0 else m$sigma.coefficients["log(age)"],
    "log(year)" = if (is.na(m$sigma.coefficients["log(year)"])) 0 else m$sigma.coefficients["log(year)"],
    "sex" = if (is.na(m$sigma.coefficients["sexwomen"])) rep(0, 2) else c(0, m$sigma.coefficients["sexwomen"]),
    "dimd" = if (length(m$sigma.coefficients[grep("dimd", names(m$sigma.coefficients))]) == 0) rep(0, 10) else c(0, m$sigma.coefficients[grep("dimd", names(m$sigma.coefficients))])
  ),
  "nu" = list(
    "intercept" = m$nu.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$nu.coefficients["log(age)"])) 0 else m$nu.coefficients["log(age)"],
    "log(year)" = if (is.na(m$nu.coefficients["log(year)"])) 0 else m$nu.coefficients["log(year)"],
    "sex" = if (is.na(m$nu.coefficients["sexwomen"])) rep(0, 2) else c(0, m$nu.coefficients["sexwomen"]),
    "dimd" = if (length(m$nu.coefficients[grep("dimd", names(m$nu.coefficients))]) == 0) rep(0, 10) else c(0, m$nu.coefficients[grep("dimd", names(m$nu.coefficients))])
  )
)

yaml::write_yaml(cc, "~/My_Models/IMPACTncd_Engl/inputs/disease_burden/andep_dur_forward2.yaml")
# yaml::write_yaml(cc, output_path("andep_dur_forward2.yaml"))


# Pain ----
dt[disease == "Pain", summary(dur)]
m1 <- gamlss(
  dur ~ log(year) + log(age) + sex + dimd,
  ~ log(year) + log(age) + sex + dimd,
  ~ log(year) + log(age) + sex + dimd,
  family = "ZANBI",
  data = dt[disease == "Pain"],
  method = mixed(20, 100)
)

m <- stepGAIC(m1, parallel = "multicore", ncpus = 16L)


cc <- list(
  "distr" = "ZANBI",
  "note" = "remember to add 1 to the final result",
  "mu" = list(
    "intercept" = m$mu.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$mu.coefficients["log(age)"])) 0 else m$mu.coefficients["log(age)"],
    "log(year)" = if (is.na(m$mu.coefficients["log(year)"])) 0 else m$mu.coefficients["log(year)"],
    "sex" = if (is.na(m$mu.coefficients["sexwomen"])) rep(0, 2) else c(0, m$mu.coefficients["sexwomen"]),
    "dimd" = if (length(m$mu.coefficients[grep("dimd", names(m$mu.coefficients))]) == 0) rep(0, 10) else c(0, m$mu.coefficients[grep("dimd", names(m$mu.coefficients))])
  ),
  "sigma" = list(
    "intercept" = m$sigma.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$sigma.coefficients["log(age)"])) 0 else m$sigma.coefficients["log(age)"],
    "log(year)" = if (is.na(m$sigma.coefficients["log(year)"])) 0 else m$sigma.coefficients["log(year)"],
    "sex" = if (is.na(m$sigma.coefficients["sexwomen"])) rep(0, 2) else c(0, m$sigma.coefficients["sexwomen"]),
    "dimd" = if (length(m$sigma.coefficients[grep("dimd", names(m$sigma.coefficients))]) == 0) rep(0, 10) else c(0, m$sigma.coefficients[grep("dimd", names(m$sigma.coefficients))])
  ),
  "nu" = list(
    "intercept" = m$nu.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$nu.coefficients["log(age)"])) 0 else m$nu.coefficients["log(age)"],
    "log(year)" = if (is.na(m$nu.coefficients["log(year)"])) 0 else m$nu.coefficients["log(year)"],
    "sex" = if (is.na(m$nu.coefficients["sexwomen"])) rep(0, 2) else c(0, m$nu.coefficients["sexwomen"]),
    "dimd" = if (length(m$nu.coefficients[grep("dimd", names(m$nu.coefficients))]) == 0) rep(0, 10) else c(0, m$nu.coefficients[grep("dimd", names(m$nu.coefficients))])
  )
)

yaml::write_yaml(cc, "~/My_Models/IMPACTncd_Engl/inputs/disease_burden/pain_dur_forward2.yaml")
# yaml::write_yaml(cc, output_path("pain_dur_forward2.yaml"))

# Constipation ----
dt[disease == "Constipation", summary(dur)]
m1 <- gamlss(
  dur ~ log(year) + log(age) + sex + dimd,
  ~ log(year) + log(age) + sex + dimd,
  ~ log(year) + log(age) + sex + dimd,
  family = "ZANBI",
  data = dt[disease == "Constipation"],
  method = mixed(20, 100)
)

m <- stepGAIC(m1, parallel = "multicore", ncpus = 16L)



cc <- list(
  "distr" = "ZANBI",
  "note" = "remember to add 1 to the final result",
  "mu" = list(
    "intercept" = m$mu.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$mu.coefficients["log(age)"])) 0 else m$mu.coefficients["log(age)"],
    "log(year)" = if (is.na(m$mu.coefficients["log(year)"])) 0 else m$mu.coefficients["log(year)"],
    "sex" = if (is.na(m$mu.coefficients["sexwomen"])) rep(0, 2) else c(0, m$mu.coefficients["sexwomen"]),
    "dimd" = if (length(m$mu.coefficients[grep("dimd", names(m$mu.coefficients))]) == 0) rep(0, 10) else c(0, m$mu.coefficients[grep("dimd", names(m$mu.coefficients))])
  ),
  "sigma" = list(
    "intercept" = m$sigma.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$sigma.coefficients["log(age)"])) 0 else m$sigma.coefficients["log(age)"],
    "log(year)" = if (is.na(m$sigma.coefficients["log(year)"])) 0 else m$sigma.coefficients["log(year)"],
    "sex" = if (is.na(m$sigma.coefficients["sexwomen"])) rep(0, 2) else c(0, m$sigma.coefficients["sexwomen"]),
    "dimd" = if (length(m$sigma.coefficients[grep("dimd", names(m$sigma.coefficients))]) == 0) rep(0, 10) else c(0, m$sigma.coefficients[grep("dimd", names(m$sigma.coefficients))])
  ),
  "nu" = list(
    "intercept" = m$nu.coefficients["(Intercept)"],
    "log(age)" = if (is.na(m$nu.coefficients["log(age)"])) 0 else m$nu.coefficients["log(age)"],
    "log(year)" = if (is.na(m$nu.coefficients["log(year)"])) 0 else m$nu.coefficients["log(year)"],
    "sex" = if (is.na(m$nu.coefficients["sexwomen"])) rep(0, 2) else c(0, m$nu.coefficients["sexwomen"]),
    "dimd" = if (length(m$nu.coefficients[grep("dimd", names(m$nu.coefficients))]) == 0) rep(0, 10) else c(0, m$nu.coefficients[grep("dimd", names(m$nu.coefficients))])
  )
)

yaml::write_yaml(cc, "~/My_Models/IMPACTncd_Engl/inputs/disease_burden/constipation_dur_forward2.yaml")
# yaml::write_yaml(cc, output_path("constipation_dur_forward2.yaml"))
