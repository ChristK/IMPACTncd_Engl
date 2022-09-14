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


setwd("/home/ckyprid/My_Models/IMPACTncd_Engl/")
univariate_analysis <- FALSE
diagnostics         <- FALSE
plots               <- FALSE
seed                <- 43L



if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}
if (!require(IMPACTncdEngl)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/",
                         force = TRUE,
                         upgrade = "never")
  }
dependencies(c("Rcpp", "dqrng", "fst", "qs", "logisticRR", "reldist", "future", "future.apply", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multiprocess)

if (file.exists("./secure_data/HSE_ts.fst")) {
  HSE_ts <- read_fst("./secure_data/HSE_ts.fst", as.data.table = TRUE)
} else {
  source("./secure_data/preprocess_HSE.R", local = TRUE)
}

distr_nam <- "BI"
BI(mu.link = "log") # To estimate RR instead of OR
dt <- na.omit(HSE_ts[wt_nurse > 0 & between(age, 20, 90),
                     .(sbp, year, age, agegrp10, sex, qimd, ethnicity, sha, bp1, wt_nurse)])
# dt[, age := scale(age, 52.1, 17.1)]
to_agegrp(dt,  5L, 85L, "age", "agegrp5" , to_factor = TRUE)
set.seed(seed)

dt[, .(dgnhtn_median = wtd.mean(as.integer(bp1) - 1, weight = wt_nurse)), keyby = .(round(sbp))
   ][, scatter.smooth(round, dgnhtn_median)]


dt[, htn := as.integer(as.character(bp1))]
dt[, lnsbp := log(sbp)]
dt[, invsbp := 1/sbp]
dt[, csbp := sqrt(1/sbp)]

m1 <- logisticRR(htn ~ sbp + agegrp5 + sex +  qimd, basecov = 140, data = dt, boot = F)
m2 <- logisticRR(htn ~ lnsbp + agegrp5 + sex +  qimd, basecov = log(140), data = dt, boot = F)
m3 <- logisticRR(htn ~ invsbp + agegrp5 + sex +  qimd, basecov = 1/140, data = dt, boot = F)
m4 <- logisticRR(htn ~ csbp + agegrp5 + sex +  qimd, basecov = sqrt(1/140), data = dt, boot = F)

AIC(m1$fit, m2$fit, m3$fit, m4$fit) # m3 is better but RR weird

m1$RR
quantile(m1$boot.rr, 0.975)
m1$delta.var
plot(m1$RR ^ ((100:200) - 140))
lines(m2$RR ^ (log(100:200) - log(140)), col = "red")
lines(m3$RR ^ (1/(100:200) - 1/(140)), col = "blue")


## print out conditioning
cn <- CJ(agegrp5 = dt$agegrp5, sex = dt$sex, qimd = dt$qimd, unique = TRUE)
results <- list()
for(i in 1:nrow(cn)){
  results[[i]] <- logisticRR(htn ~ lnsbp + agegrp5 + sex +  qimd, basecov = log(140),
                             data = dt, fixcov = cn[i,], boot = FALSE)

}
cn[, rr := sapply(results, `[[`, "RR")]
cn[, var := sapply(results, `[[`, "delta.var")]
cn[, ci_rr := exp(log(rr) + qnorm(0.975) * var)]
cn[, var := NULL]
setnames(cn, "agegrp5", "agegroup")
View(cn)

# ---
#   xps_name: sbp
# outcome: htn
# lag: 0
# distribution: lognormal
# source: HSE own estimation
# notes: ''
# apply_rr_extra_fn: >
#   function(sp) {
#     if (!inherits(sp, 'SynthPop')) stop('Argument sp needs to be a SynthPop object.')
#     ideal_xps_lvl <- log(140)
#     sp$pop[, sbp_rr := sbp_rr^(log(sbp) - ideal_xps_lvl)]
#   }
# ---
ll <- agegrp_name(min_age = 0, 19, grp_width = 5,
            match_input = F, match_input_max_age = 19)
cn2 <- CJ(agegroup = ll, sex = dt$sex, qimd = dt$qimd, unique = TRUE, rr = 1, ci_rr = 1)
cn <- rbind(cn2, cn)
setkey(cn, agegroup, sex, qimd)
fwrite(cn, "./inputs/rr/sbp~htn.csvy")
