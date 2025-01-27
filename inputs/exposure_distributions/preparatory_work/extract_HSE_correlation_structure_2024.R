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



n_cpu <- 20L
impute <- TRUE
iterations <- 100L


if (!require(CKutils)) {
  if (!require(remotes)) install.packages("remotes")
  remotes::install_github("ChristK/CKutils")
}

# if (!require(IMPACTncdEnglmisc)) {
#   if (!require(remotes))
#     install.packages("remotes")
#   remotes::install_local("./Rpackage/IMPACTncd_Engl_model_pkg/")
#   library(IMPACTncdEnglmisc)
# }

dependencies(c("Rcpp", "dqrng", "fst", "qs", "gamlss.dist", "future",
               "future.apply", "ggplot2", "data.table"))
options(future.fork.enable = TRUE) # enable fork in Rstudio
plan(multisession, workers = n_cpu)

  HSE_ts <- read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)[between(age, 20L, 90L)]


cm <- vector("list", iterations)
cm <-
  future_lapply(1:iterations, function(x) {
  # Extract percentiles & impute missing values
  HSE_ts <-
    read_fst("./secure_data/HSE_ts_03_19.fst", as.data.table = TRUE)[between(age, 20L, 90L)]
 # HSE_ts[, smok_dur_curr := as.integer(round(age - smok_init_age))] this is tricky
  limit_year <- HSE_ts[, .(min = min(year), max = max(year))]
  HSE_ts[is.na(smok_status), c("smok_cig_ex", "smok_cig_curr",
                               "smok_init_age", "smok_dur_curr") := NA]
  setindexv(HSE_ts, c("year", "age", "sex", "sha", "qimd", "ethnicity_grp", 
                      "urban_rural", "smok_status"))


  dqRNGkind("pcg64") # dqRNGkind("Xoroshiro128+") ~10% faster
  SEED <- 387705L # sample(1e7, 1)
  dqset.seed(SEED, x) # Ensure that seed hasn't been used elsewhere in the model
  set.seed(SEED + x)

  # impute ethnicity
  tbl <-
    read_fst("./inputs/exposure_distributions/ethnicity_grp_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  HSE_ts[, rank_var := dqrunif(.N)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(ethnicity_grp),
         ethnicity_grp := levels(ethnicity_grp)[(rank_var > white) +
                                          (rank_var > black) +
                                          (rank_var > asian) +
                                          (rank_var > mixed) +
                                          (rank_var > other) + 1]]
  HSE_ts[, c(col_nam, "rank_var") := NULL]

  # impute education
  tbl <-
    read_fst("./inputs/exposure_distributions/education_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  HSE_ts[, rank_var := dqrunif(.N)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(education),
         education := levels(education)[(rank_var > ed1) +
                                          (rank_var > ed2) +
                                          (rank_var > ed3) +
                                          (rank_var > ed4) +
                                          (rank_var > ed5) +
                                          (rank_var > ed6) + 1]]
  HSE_ts[, c(col_nam, "rank_var") := NULL]

  HSE_ts[, education_l := as.integer(education)]
  tbl <-
    melt(
      tbl,
      id.vars = 1:sum(sapply(tbl, class) != "numeric"),
      variable.name = "education_l",
      value.name = "quant_max"
    )
  tbl[, education_l := as.integer(gsub("^ed", "", education_l, perl = TRUE))]
  absorb_dt(HSE_ts, tbl)
  HSE_ts[education_l == max(education_l), quant_max := 1]
  tbl[, education_l := education_l + 1L]
  setnames(tbl, "quant_max", "quant_min")
  absorb_dt(HSE_ts, tbl)
  HSE_ts[education_l == min(education_l), quant_min := 0]
  HSE_ts[, education_r := runif(.N, quant_min, quant_max)]
  HSE_ts[, c("education_l", "quant_min", "quant_max") := NULL]


  # impute income
  tbl <-
    read_fst("./inputs/exposure_distributions/income_table.fst",
             as.data.table = TRUE)
  HSE_ts[, rank_var := dqrunif(.N)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(income),
         income := levels(income)[(rank_var > inc1) +
                                    (rank_var > inc2) +
                                    (rank_var > inc3) +
                                    (rank_var > inc4) + 1]]
  HSE_ts[, c(col_nam, "rank_var") := NULL]

  HSE_ts[, income_l := as.integer(income)]
  tbl <-
    melt(
      tbl,
      id.vars = 1:sum(sapply(tbl, class) != "numeric"),
      variable.name = "income_l",
      value.name = "quant_max"
    )
  tbl[, income_l := as.integer(gsub("^inc", "", income_l, perl = TRUE))]
  absorb_dt(HSE_ts, tbl)
  HSE_ts[income_l == max(income_l), quant_max := 1]
  tbl[, income_l := income_l + 1L]
  setnames(tbl, "quant_max", "quant_min")
  absorb_dt(HSE_ts, tbl)
  HSE_ts[income_l == min(income_l), quant_min := 0]
  HSE_ts[, income_r := runif(.N, quant_min, quant_max)]
  HSE_ts[, c("income_l", "quant_min", "quant_max") := NULL]

  # impute pa [polr]
  tbl <-
    read_fst("./inputs/exposure_distributions/active_days_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  HSE_ts[, rank_var := dqrunif(.N)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(active_days),
         active_days := (rank_var > pa1) +
           (rank_var > pa2) +
           (rank_var > pa3) +
           (rank_var > pa4) +
           (rank_var > pa5) +
           (rank_var > pa6) + 1]
  HSE_ts[, c(col_nam, "rank_var") := NULL]

  HSE_ts[, active_days_l := as.integer(active_days)]
  tbl <-
    melt(
      tbl,
      id.vars = 1:sum(sapply(tbl, class) != "numeric"),
      variable.name = "active_days_l",
      value.name = "quant_max"
    )
  tbl[, active_days_l := as.integer(gsub("^pa", "", active_days_l, perl = TRUE))]
  absorb_dt(HSE_ts, tbl)
  HSE_ts[active_days_l == max(active_days_l), quant_max := 1]
  tbl[, active_days_l := active_days_l + 1L]
  setnames(tbl, "quant_max", "quant_min")
  absorb_dt(HSE_ts, tbl)
  HSE_ts[active_days_l == min(active_days_l), quant_min := 0]
  HSE_ts[, active_days_r := runif(.N, quant_min, quant_max)]
  HSE_ts[, c("active_days_l", "quant_min", "quant_max") := NULL]

  # impute fruit [ZINBI]
  tbl <-
    read_fst("./inputs/exposure_distributions/frtpor_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(frtpor), frtpor := qZINBI(dqrunif(.N), mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[, frtpor_r := pZINBI(frtpor, mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]

  # impute veg
  tbl <-
    read_fst("./inputs/exposure_distributions/vegpor_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(vegpor), vegpor := my_qDEL(dqrunif(.N), mu, sigma, nu)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[, vegpor_r := my_pDEL(vegpor, mu, sigma, nu)]
  HSE_ts[, (col_nam) := NULL]

  # impute smoking status [MN4]
  tbl <-
    read_fst("./inputs/exposure_distributions/smok_status_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(smok_status), smok_status := my_qMN4(dqrunif(.N), mu, sigma, nu)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[, smok_status_r := pMN4(as.integer(smok_status), mu, sigma, nu)]
  HSE_ts[, (col_nam) := NULL]

  # impute quit years
  tbl <-
    read_fst("./inputs/exposure_distributions/smok_quit_yrs_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(smok_quit_yrs) & smok_status %in% c("2", "3"),
         smok_quit_yrs := my_qDPO(dqrunif(.N, 0, 0.99), mu, sigma)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[smok_status %in% c("2", "3"),
         smok_quit_yrs_r := my_pDPO(smok_quit_yrs, mu, sigma)]
  HSE_ts[, (col_nam) := NULL]

  # impute smok duration for ex
  tbl <-
    read_fst("./inputs/exposure_distributions/smok_dur_ex_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  tbl[, smok_status := factor(smok_status, levels = 1:4)]
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(smok_dur_ex) & smok_status %in% c("2", "3"),
         smok_dur_ex := my_qDPO(dqrunif(.N, 0, 0.99), mu, sigma)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[smok_status %in% c("2", "3"),
         smok_dur_ex_r := my_pDPO(smok_dur_ex, mu, sigma)]
  HSE_ts[, (col_nam) := NULL]

  # impute smok duration for curr
  tbl <-
    read_fst("./inputs/exposure_distributions/smok_dur_curr_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(smok_dur_curr) & smok_status %in% c("4"),
         smok_dur_curr := as.integer(round(qNBI(dqrunif(.N, 0, 0.88), mu, sigma)))]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[smok_status %in% c("4"),
         smok_dur_curr_r := pNBI(smok_dur_curr, mu, sigma)]
  HSE_ts[, (col_nam) := NULL]

  # impute smok cig ex
  tbl <-
    read_fst("./inputs/exposure_distributions/smok_cig_ex_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(smok_cig_ex) & smok_status %in% c("3"),
         smok_cig_ex := my_qZABNB(dqrunif(.N), mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[smok_status %in% c("3"), # for "2" is NaN
         smok_cig_ex_r := pZABNB(smok_cig_ex, mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]

  # impute smok cig curr
  tbl <-
    read_fst("./inputs/exposure_distributions/smok_cig_curr_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(smok_cig_curr) & smok_status %in% c("4"),
         smok_cig_curr := qZINBI(dqrunif(.N), mu, sigma, nu)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[smok_status %in% c("4"),
         smok_cig_curr_r := pZINBI(smok_cig_curr, mu, sigma, nu)]
  HSE_ts[, (col_nam) := NULL]

  # impute ETS [BI]
  tbl <-
    read_fst("./inputs/exposure_distributions/ets_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(ets), ets := factor(qBI(dqrunif(.N), 1, mu))]

  HSE_ts[, ets_r := pBI(as.integer(as.character(ets)), 1, mu)]
  HSE_ts[, (col_nam) := NULL]


  # impute alcohol (totalwu) [ZINBI]
  tbl <-
    read_fst("./inputs/exposure_distributions/alcohol_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(totalwu), totalwu := qZINBI(dqrunif(.N), mu, sigma, nu)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[, totalwu_r := pZINBI(totalwu, mu, sigma, nu)]
  HSE_ts[, (col_nam) := NULL]

  # impute bmi [GB2]
  tbl <-
    read_fst("./inputs/exposure_distributions/bmi_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(bmi), bmi := qGB2(dqrunif(.N), mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[, bmi_r := pGB2(bmi, mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]

  # impute sbp [GB2]
  tbl <-
    read_fst("./inputs/exposure_distributions/sbp_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(sbp), sbp := qGB2(dqrunif(.N), mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[, sbp_r := pGB2(sbp, mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]

  # impute bp medication
  HSE_ts[, `:=` (sbp_acc = sbp,
              sbp = round(clamp(sbp, 110, 200), -1))]
  tbl <-
    read_fst("./inputs/exposure_distributions/bp_med_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <- setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(bp_med), bp_med := factor(qbinom(dqrunif(.N), 1, mu))]

  HSE_ts[, bp_med_r := pbinom(as.integer(as.character(bp_med)), 1, mu)]
  HSE_ts[, (col_nam) := NULL]

  HSE_ts[, `:=` (sbp = sbp_acc,
              sbp_acc = NULL)]

  # impute tchol [ST5]
  tbl <-
    read_fst("./inputs/exposure_distributions/tchol_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(tchol), tchol := qST5(dqrunif(.N), mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]

  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[, tchol_r := pST5(tchol, mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]


  # HDL to Tchol
  HSE_ts[, hdl_to_tchol := hdl / tchol]
  HSE_ts[hdl_to_tchol >= 1, hdl_to_tchol := NA]
  tbl <-
    read_fst("./inputs/exposure_distributions/hdl_to_tchol_table.fst",
             as.data.table = TRUE)
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(hdl_to_tchol), hdl_to_tchol := qGB1(dqrunif(.N), mu, sigma, nu, tau)]

  HSE_ts[, hdl_to_tchol_r := pGB1(hdl_to_tchol, mu, sigma, nu, tau)]
  HSE_ts[, (col_nam) := NULL]


  # impute AF_diagn
  tbl <-
    read_fst("./inputs/exposure_distributions/af_dgn_table.fst",
             as.data.table = TRUE)
  col_nam <- setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(af), af := factor(qbinom(dqrunif(.N), 1, mu))]

  HSE_ts[, af_r := pbinom(as.integer(as.character(af)), 1, mu)]
  HSE_ts[, (col_nam) := NULL]

  # impute famcvd
  tbl <-
    read_fst("./inputs/exposure_distributions/famcvd_table.fst",
             as.data.table = TRUE)
  col_nam <- setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(famcvd), famcvd := factor(qbinom(dqrunif(.N), 1, mu))]

  HSE_ts[, famcvd_r := pbinom(as.integer(as.character(famcvd)), 1, mu)]
  HSE_ts[, (col_nam) := NULL]

  # impute t2dm
  HSE_ts[, `:=` (bmi_acc = bmi,
                 bmi = round(clamp(bmi, 18, 50)))]
  tbl <-
    read_fst("./inputs/exposure_distributions/dm_table.fst",
             as.data.table = TRUE)
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(dm), dm := factor(qbinom(dqrunif(.N), 1, mu))]

  HSE_ts[, dm_r := pbinom(as.integer(as.character(dm)), 1, mu)]
  HSE_ts[, (col_nam) := NULL]


  # impute t2dm dgn
  HSE_ts[, `:=` (bmi = round(clamp(bmi, 18, 50), -1))]
  tbl <-
    read_fst("./inputs/exposure_distributions/dm_dgn_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(dm_dgn) & dm == "1", dm_dgn := factor(qbinom(dqrunif(.N), 1, mu))]

  HSE_ts[, dm_dgn_r := pbinom(as.integer(as.character(dm_dgn)), 1, mu)]
  HSE_ts[, (col_nam) := NULL]
  HSE_ts[, `:=` (bmi = bmi_acc,
                 bmi_acc = NULL)]

  # statin
  HSE_ts[, `:=` (tchol_acc = tchol,
                 tchol = round(clamp(tchol, 2, 12), 0))]
  tbl <-
    read_fst("./inputs/exposure_distributions/statin_px_table.fst",
             as.data.table = TRUE)[between(year, limit_year$min, limit_year$max)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(statin_px), statin_px := qbinom(dqrunif(.N), 1, mu)]

  HSE_ts[, statin_px_r := pbinom(statin_px, 1, mu)]
  HSE_ts[, (col_nam) := NULL]
  HSE_ts[, `:=` (tchol = tchol_acc,
              tchol_acc = NULL)]

  # impute CKD
  tbl <-
    read_fst("./inputs/exposure_distributions/ckd_table.fst",
             as.data.table = TRUE)
  HSE_ts[, rank_var := dqrunif(.N)]
  col_nam <-
    setdiff(names(tbl), intersect(names(HSE_ts), names(tbl)))
  absorb_dt(HSE_ts, tbl)
  HSE_ts[is.na(ckd),
         ckd := levels(ckd)[(rank_var > ckd0) +
                            (rank_var > ckd1) +
                            (rank_var > ckd2) +
                            (rank_var > ckd3) + 1]]
  HSE_ts[, c(col_nam, "rank_var") := NULL]

  HSE_ts[, ckd_l := as.integer(ckd)]
  tbl <-
    melt(
      tbl,
      id.vars = 1:sum(sapply(tbl, class) != "numeric"),
      variable.name = "ckd_l",
      value.name = "quant_max"
    )
  tbl[, ckd_l := as.integer(gsub("^ckd", "", ckd_l, perl = TRUE))]
  tbl[, ckd_l := ckd_l + 1L]
  absorb_dt(HSE_ts, tbl)
  HSE_ts[ckd_l == max(ckd_l), quant_max := 1]
  tbl[, ckd_l := ckd_l + 1L]
  setnames(tbl, "quant_max", "quant_min")
  absorb_dt(HSE_ts, tbl)
  HSE_ts[ckd_l == min(ckd_l), quant_min := 0]
  HSE_ts[, ckd_r := runif(.N, quant_min, quant_max)]
  HSE_ts[, c("ckd_l", "quant_min", "quant_max") := NULL]

  data.table(HSE_ts[, cor(.SD, use = "p", method = "s"),
                      .SDcols = patterns("_r$")], keep.rownames = TRUE)[, mc := x]
}
)

exposure_corr <- rbindlist(cm, idcol = FALSE)
write_fst(exposure_corr, "./inputs/exposure_distributions/exposure_corr_mc.fst", 100L)

exposure_corr[, mc := NULL]
cm_sd <- exposure_corr[, lapply(.SD, sd), by = rn]
write_fst(cm_sd, "./inputs/exposure_distributions/exposure_corr_sd.fst", 100L)
cm_sd[, lapply(.SD, function(x) round(max(x, na.rm = TRUE), 3)), .SDcols = -"rn"]


# SD very small so imputation has little impact and just a few iterations are OK
# Moreover simply using the mean for synth pop generation is adequate
# Note that NAs are expected and have no impact on the generation of RNs

cm_mean <- exposure_corr[, lapply(.SD, mean), by = rn]
tt <- which(sapply(cm_mean, typeof, USE.NAMES = FALSE) == "double")
setnafill(cm_mean, "c", 0, cols = tt)
if (min(eigen(as.matrix(cm_mean,
  rownames = "rn"))$values) >= 0) print("matrix is semi-positive definite")
write_fst(cm_mean, "./inputs/exposure_distributions/exposure_corr_mean.fst", 100L)
cm_mean[, lapply(.SD, function(x) round(max(x[x < 1], na.rm = TRUE), 3)), .SDcols = -"rn"]


# plot results (marginals)
# from http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
cm_mean <- read_fst(normalizePath("./inputs/exposure_distributions/exposure_corr_mean.fst"), as.data.table = TRUE)

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
corrplot::corrplot(as.matrix(cm_mean, rownames = TRUE), method="color", col=col(200),
  type = "full", order = "a",
  #addCoef.col = "black", # Add coefficient of correlation
  tl.col = "black", tl.srt = 45, tl.cex = 0.7, #Text label color and rotation
  # hide correlation coefficient on the principal diagonal
  diag = FALSE, title = "HSE"
)





