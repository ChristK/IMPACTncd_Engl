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


setwd("/home/ckypridMy_ModelsIMPACTncd_Engl/")
if (!require(CKutils)) {
  if (!require(remotes))
    install.packages("remotes")
  remotes::install_bitbucket("ChristK/CKutils")
}
dependencies(c("fst", "gamlss", "reldist", "data.table"))

HSE_ts <- read_fst("./preparatory_work/HSE_ts.fst", as.data.table = TRUE)
dt <- HSE_ts[wt_int > 0, .(
  wt_int, smok_init_age, smok_status,
  agegrp, year, age, sex, qimd, ethnicity, sha, smok_dur, smok_quit_yrs)]
dt[, wt_int := wt_int*10000/sum(wt_int), by = year] # make pop of each hse =10000
dt[smok_status == 4 & between(age - smok_init_age, 0, 3),
   `:=`(year = year - (age - smok_init_age))]
dt[smok_status == 4 & between(age - smok_init_age, 0, 3),
   `:=`(age = smok_init_age, smok_incid = 1L)]
dt[smok_status %in% c("2", "3") & between(smok_quit_yrs + smok_dur, 0, 3),
   `:=`(
     age  = age  - smok_quit_yrs - smok_dur,
     year = year - smok_quit_yrs - smok_dur,
     smok_incid = 1L
   )]
dt[smok_status == "1", smok_incid := 0L]

dt <- na.omit(dt[between(year, 3, 11) & age > 17, .(
  wt_int, smok_incid, agegrp, year, age, sex, qimd, ethnicity, sha)])
# dt[, age := scale(age, 50, 16.8)]

set.seed(445)
lns <- sample(nrow(dt), nrow(dt) * 0.8)
dt_trn   <- dt[lns] # train dataset
dt_crv   <- dt[!lns]  # cross-validation dataset
diagnostics <- FALSE

distr_nam <- "BI" # from preparatory analysis
con1 <- gamlss.control(c.crit = 1e-3) # increase for faster exploratory analysis.

setMKLthreads(12L)
mod_min <- gamlss(
  smok_incid ~ 1,
  family = distr_nam,
  weights = dt_trn$wt_int,
  data = dt_trn,
  method = mixed(5, 100),
  control = con1
)

getTGD <- function (object, newdata = NULL, ...) {
  if (!is.gamlss(object))
    stop("not a gamlss object")
  fname <- as.gamlss.family(object$family[[1]])
  dfun <- paste("d", fname$family[[1]], sep = "")
  pfun <- paste("p", fname$family[[1]], sep = "")
  lpar <- length(fname$parameters)
  if (is.null(newdata))
    stop("no newdata is set in VGD")
  nfitted <- predictAll(object, newdata = newdata, ...)
  if (is.null(nfitted$y))
    stop("the response variables is missing in the newdata")
  if (fname$family[1] %in% .gamlss.bi.list) {
    if (NCOL(nfitted$y) == 1) {
      y1 <- nfitted$y
      bd <- rep(max(y1), length(y1))
    } else {
      bd <- nfitted$y[, 1] + nfitted$y[, 2]
      y1 <- nfitted$y[, 1]
    }
  } else {
    y1 <- nfitted$y
  }
  if (lpar == 1) {
    if (fname$family[[1]] %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfitted$mu, bd = bd,
                   log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, bd = bd)
    } else {
      devi <- call(dfun, x = y1, mu = nfitted$mu, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu)
    }
  } else if (lpar == 2) {
    if (fname$family[[1]] %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$nmu, sigma = nfitted$sigma,
                   , bd = bd)
    } else {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma)
    }
  } else if (lpar == 3) {
    if (fname$family[[1]] %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, bd = bd, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, bd = bd)
    } else {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu)
    }
  } else {
    if (fname$family[[1]] %in% .gamlss.bi.list) {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, tau = nfitted$tau, bd = bd,
                   log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, tau = nfitted$tau, bd = bd)
    } else {
      devi <- call(dfun, x = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, tau = nfitted$tau, log = TRUE)
      ures <- call(pfun, q = y1, mu = nfitted$mu, sigma = nfitted$sigma,
                   nu = nfitted$nu, tau = nfitted$tau)
    }
  }
  Vresid <- qNO(eval(ures))
  dev <- -2 * sum(eval(devi))
  out <- list()
  out <- list(TGD = dev, predictError = dev/dim(newdata)[1],
              resid = Vresid)
  class(out) <- "gamlssTGD"
  out
}
assignInNamespace("getTGD", getTGD, "gamlss") # Fixes a bug in gamlss

setMKLthreads(1L)
smok_incid_model <- stepTGDAll.A(
  mod_min,
  scope = list(
    lower =  ~ 1,
    upper =  ~ (log(year) + pcat(qimd) + log(age) + sex +
        pcat(sha) + pcat(ethnicity)
    ) ^ 2
  ),
  newdata = dt_crv,
  parallel = "multicore",
  ncpus = 16L
)

smok_incid_model$data <- dt_trn

saveRDS(smok_incid_model, "./secure_data/lifecourse_models/smok_incid_model.rds")
print("Model saved.")

if (diagnostics) {
  smok_incid_model <- readRDS("./secure_data/lifecourse_models/smok_incid_model.rds")

  wp(smok_incid_model)
  wp(smok_incid_model, xvar = age)
  plot(smok_incid_model)

  zz <- validate_gamlss(dt_crv, smok_incid_model, 10, dt_trn)
  zz[, weight := wt_int/sum(wt_int), by = type]
  reldist_diagnostics(zz[type == "Observed", smok_incid],
                      zz[type == "Modelled", smok_incid],
                      zz[type == "Observed", weight],
                      zz[type == "Modelled", weight],
                      main = expression(bold(smok_incid)),
                      discrete = TRUE)
  dependencies("ggplot2")

  zz[, weight := wt_int/sum(wt_int), by = .(type, agegrp)]
  ggplot(zz, aes(smok_incid, colour = type, weight = weight, linetype = type)) +
    stat_ecdf() +
    facet_wrap(.~agegrp, nrow = 3) + ggtitle("Age group")

  zz[, weight := wt_int/sum(wt_int), by = .(type, year)]
  ggplot(zz, aes(smok_incid, colour = type, weight = weight)) +
    stat_ecdf() +
    facet_wrap(.~year, nrow = 3) + ggtitle("Year")

  zz[, weight := wt_int/sum(wt_int), by = .(type, qimd)]
  ggplot(zz, aes(smok_incid, colour = type, weight = weight)) +
    stat_ecdf() + ylim(c(0.95 ,1)) +
    facet_wrap(.~qimd, nrow = 1) + ggtitle("QIMD")

  zz[, weight := wt_int/sum(wt_int), by = .(type, sha)]
  ggplot(zz, aes(smok_incid, colour = type, weight = weight)) +
    stat_ecdf() +
    facet_wrap(.~sha, nrow = 3) + ggtitle("SHA")

  zz[, weight := wt_int/sum(wt_int), by = .(type, ethnicity)]
  ggplot(zz, aes(smok_incid, colour = type, weight = weight)) +
    stat_ecdf() +
    facet_wrap(.~ethnicity, nrow = 3) + ggtitle("Ethnicity")

  zz[, weight := wt_int/sum(wt_int), by = .(type, year, agegrp)]
  ggplot(zz, aes(smok_incid, colour = type, weight = weight)) +
    stat_ecdf() +
    facet_grid(year~agegrp) + ggtitle("Year ~ Age group")

  zz[, weight := wt_int/sum(wt_int), by = .(type, year, sha)]
  ggplot(zz, aes(smok_incid, colour = type, weight = weight)) +
    stat_ecdf() +
    facet_grid(year~sha) + ggtitle("Year ~ SHA")

  zz[, weight := wt_int/sum(wt_int), by = type]
  tt <- zz[, prop_if(smok_incid == 0), by = .(type, year)]
  ggplot(tt, aes(year, V1, colour = type)) +
    geom_line() + ylim(0, 1)

  tt <- zz[, prop_if(smok_incid == 0), by = .(type, year, agegrp)]
  ggplot(tt, aes(year, V1, colour = type)) +
    geom_line() + ylim(0, 1) +
    facet_grid(agegrp~.) + ggtitle("Year ~ Age group")

    tt <- zz[, prop_if(smok_incid == 0), by = .(type, year, qimd, agegrp)]
  ggplot(tt, aes(year, V1, colour = type)) +
    geom_line() + ylim(0, 1) +
    facet_grid(qimd~agegrp) + ggtitle("Year ~ qimd")

  zz[, .N, by = .(type, year, qimd, agegrp)][, hist(N)]

}



