## IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncd_Engl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncd_Engl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.


setwd("/home/ckyprid/My Models/IMPACTncd_Engl/")
if (!require(CKutils)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_bitbucket("ChristK/CKutils")
}
dependencies(c("fst", "gamlss", "gamlss.add", "reldist", "data.table"))

HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- na.omit(HSE_ts[wt_int > 0 & smok_status %in% c("2", "3") & age > 15, .(
  smok_dur_ex, smok_status, year, age, agegrp, sex, qimd, ethnicity, sha, wt_int)]
)
dt[, smok_dur_ex := round(smok_dur_ex)]
dt[, age := scale(age, 50, 16.8)]
set.seed(49)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
diagnostics <- FALSE


distr_nam <- "DPO"
con1 <- gamlss.control(c.crit = 1e-1) # increase for faster exploratory analysis.


setMKLthreads(12L)
mod_min <- gamlss(
  smok_dur_ex ~ 1,
  family = distr_nam,
  weights = dt_trn$wt_int,
  data = dt_trn,
  method = mixed(5, 100),
  control = con1
)

setMKLthreads(1L)
smok_dur_ex_model <- stepTGDAll.A(
  mod_min,
  scope = list(
    lower =  ~ 1,
    upper =  ~ (log(year) + age + sex + smok_status + pcat(qimd) +
        pcat(ethnicity) + pcat(sha)
    ) ^ 2
  ),
  sigma.scope = list(
    lower =  ~ 1,
    upper =  ~ (log(year) + age + sex + smok_status + pcat(qimd) +
                  pcat(ethnicity) + pcat(sha)
    ) ^ 2
  ),
  # nu.scope = list(
  #   lower =  ~ 1,
  #   upper =  ~ (log(year) + age + sex + smok_status + pcat(qimd) +
  #                 pcat(ethnicity) + pcat(sha)
  #   ) ^ 2
  # ),
  # tau.scope = list(
  #   lower =  ~ 1,
  #   upper =  ~ ga( ~ s(log(year), age, by = sex)) + (
  #     log(year) + pb(age) + sex + pcat(qimd) +
  #       pcat(ethnicity) + pcat(sha)
  #   ) ^ 2
  # ),
  newdata = dt_crv,
  parallel = "multicore",
  ncpus = 16L
)


con2 <- gamlss.control(c.crit = 1e-3)
smok_dur_ex_model <- update(smok_dur_ex_model, control = con2,  method = mixed(5, 200))

smok_dur_ex_model$data <- dt_trn

saveRDS(smok_dur_ex_model, "./secure_data/lifecourse_models/smok_dur_ex_model.rds")
print("Model saved.")

if (diagnostics) {
  smok_dur_ex_model <- readRDS("./secure_data/lifecourse_models/smok_dur_ex_model.rds")

  wp(smok_dur_ex_model)
  wp(smok_dur_ex_model, xvar = age)
  plot(smok_dur_ex_model)

  zz <- validate_gamlss(dt, smok_dur_ex_model, 10, smok_dur_ex_model$data)[smok_dur_ex < 100]
  zz[, weight := wt_int/sum(wt_int), by = type]
  reldist_diagnostics(zz[type == "Observed", smok_dur_ex],
                      zz[type == "Modelled", smok_dur_ex],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = expression(bold(Alcohol~(g/d))),
                      discrete = TRUE)

  dependencies("ggplot2")
  zz[, weight := wt_int/sum(wt_int), by = .(type, agegrp)]
  ggplot(zz, aes(smok_dur_ex, colour = type, weight = weight, linetype = type)) +
    # geom_freqpoly() +
    stat_ecdf() +
    facet_wrap(.~agegrp, nrow = 3) + ggtitle("Age group")

  zz[, weight := wt_int/sum(wt_int), by = .(type, year)]
  ggplot(zz, aes(smok_dur_ex, colour = type, weight = weight)) +
    geom_freqpoly() +
    facet_wrap(.~year, nrow = 3) + ggtitle("Year")

  zz[, weight := wt_int/sum(wt_int), by = .(type, qimd)]
  ggplot(zz, aes(smok_dur_ex, colour = type, weight = weight)) +
    # geom_freqpoly() +
    stat_ecdf() +
    facet_wrap(.~qimd, nrow = 3) + ggtitle("QIMD")

  zz[, weight := wt_int/sum(wt_int), by = .(type, sha)]
  ggplot(zz, aes(smok_dur_ex, colour = type, weight = weight)) +
    # geom_freqpoly() +
    stat_ecdf() +
    facet_wrap(.~sha, nrow = 3) + ggtitle("SHA")

  zz[, weight := wt_int/sum(wt_int), by = .(type, ethnicity)]
  ggplot(zz, aes(smok_dur_ex, colour = type, weight = weight)) +
    geom_freqpoly() +
    facet_wrap(.~ethnicity, nrow = 3) + ggtitle("Ethnicity")

  zz[, weight := wt_int/sum(wt_int), by = .(type, year, agegrp)]
  ggplot(zz, aes(smok_dur_ex, colour = type, weight = weight)) +
    # geom_freqpoly() +
    stat_ecdf() +
    facet_grid(year~agegrp) + ggtitle("Year ~ Age group")

  zz[, weight := wt_int/sum(wt_int), by = .(type, year, sha)]
  ggplot(zz, aes(smok_dur_ex, colour = type, weight = weight)) +
    # geom_freqpoly() +
    stat_ecdf() +
    facet_grid(year~sha) + ggtitle("Year ~ SHA")
}



