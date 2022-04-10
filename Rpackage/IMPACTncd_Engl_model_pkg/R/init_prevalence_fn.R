## IMPACTncdEngl is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz. This work has been
## funded by NIHR  HTA Project: 16/165/01 - IMPACTncdEngl: Health Outcomes
## Research Simulation Environment.  The views expressed are those of the
## authors and not necessarily those of the NHS, the NIHR or the Department of
## Health.
##
## Copyright (C) 2018-2020 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngl is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.



#' @export
init_prevalence <- function(mc, dt, design_, timing = TRUE) {
  if (timing) ptm <- proc.time()

  # Generate chd prevalence and duration (BI, NBI) ----
  simulate_init_prvl(mc, "chd", design_, dt)

  tbl <-
    read_fst("./inputs/exposure_distributions/chd_duration_table.fst",
      as.data.table = TRUE)
  tbl[, chd_prvl := 1L]
  col_nam <- setdiff(names(tbl), names(dt))
  absorb_dt(dt, tbl)
  dt[chd_prvl == 1L,
     chd_prvl := 2L + qNBI(dqrunif(.N), clamp(mu - 2, 0, Inf), sigma)] # +2 to avoid 0 and bias prevalence
  dt[chd_prvl > age, chd_prvl := age]
  dt[, (col_nam) := NULL]
  dt[, chd_prvl := carry_backward(chd_prvl, pid_mrk)] # Not strictly necessary
  # dt[year == 13, table(chd_prvl)]
  # View(dt[pid %in% c(1346, 1347), .(pid, pid_mrk, year, chd_prvl)])
  # dt[, chd_prvl:=NULL]

  # Generate stroke prevalence and duration (BI, PIG) ----
  # simulate_init_prvl(mc, "stroke", design_, dt)
  #
  # tbl <-
  #   read_fst("./inputs/exposure_distributions/stroke_duration_table.fst",
  #     as.data.table = TRUE)
  # tbl[, stroke_prvl := 1L]
  # col_nam <- setdiff(names(tbl), names(dt))
  # absorb_dt(dt, tbl)
  # dt[stroke_prvl == 1L,
  #    stroke_prvl := 2L + qPIG(dqrunif(.N), clamp(mu - 2, 0, Inf), sigma)] # +2 to avoid 0 and bias prevalence
  # dt[stroke_prvl > age, stroke_prvl := age]
  # dt[, (col_nam) := NULL]
  # dt[, stroke_prvl := carry_backward(stroke_prvl, pid_mrk)] # Not strictly necessary

  # Generate poststroke dementia prevalence and duration ----
  # tbl <- RR$stroke_dementia$get_rr(mc, design_, drop = TRUE)
  # setnames(tbl, "stroke_rr", "poststroke_dementia_prb")
  # absorb_dt(dt, tbl)
  # setnafill(dt, "c", 0, cols = "poststroke_dementia_prb")
  # dt[stroke_prvl > 1L, poststroke_dementia_prvl := as.integer(dqrunif(.N) < poststroke_dementia_prb)]
  # dt[poststroke_dementia_prvl == 1L, poststroke_dementia_prvl := stroke_prvl - 1L]
  # setnafill(dt, "c", 0L, cols = "poststroke_dementia_prvl")
  # dt[, poststroke_dementia_prb := NULL]

  # Generate copd prevalence and duration (BI, GEOM) ----
  # simulate_init_prvl(mc, "copd", design_, dt)
  #
  # tbl <-
  #   read_fst("./inputs/exposure_distributions/copd_duration_table.fst",
  #     as.data.table = TRUE)
  # tbl[, copd_prvl := 1L]
  # col_nam <- setdiff(names(tbl), names(dt))
  # absorb_dt(dt, tbl)
  # dt[copd_prvl == 1L,
  #    copd_prvl := 2L + qGEOM(dqrunif(.N), clamp(mu - 2, 0, Inf), sigma)] # +2 to avoid 0 and bias prevalence
  # dt[copd_prvl > age, copd_prvl := age]
  # dt[, (col_nam) := NULL]
  # dt[, copd_prvl := carry_backward(copd_prvl, pid_mrk)] # Not strictly necessary

  # Generate breast_cancer prevalence and duration (BI, PO) ----
  # simulate_init_prvl(mc, "breast_ca", design_, dt)
  #
  # # duration
  # tbl <- get_disease_epi_mc(mc, "breast_ca", "d", "m", FALSE)
  # tbl[, year := design_$sim_prm$init_year] # given init year is 13 and close enough to 11 which I have epi for
  # col_nam <- setdiff(names(tbl), names(dt))
  # absorb_dt(dt, tbl)
  # dt[breast_ca_prvl == 1L, breast_ca_prvl := 2L + rpois(.N, clamp(duration - 2, 0, 100))]
  # dt[breast_ca_prvl > design_$sim_prm$cancer_cure,
  #    breast_ca_prvl := dqsample(2:design_$sim_prm$cancer_cure, .N, TRUE)]
  # dt[breast_ca_prvl > age, breast_ca_prvl := age]
  # dt[, (col_nam) := NULL]
  # dt[, breast_ca_prvl := carry_backward(breast_ca_prvl, pid_mrk)] # Not strictly necessary

  # Generate colon_cancer prevalence and duration (BI, PO) ----
  # simulate_init_prvl(mc, "colon_ca", design_, dt)
  #
  # # duration
  # tbl <- get_disease_epi_mc(mc, "colon_ca", "d", "m", FALSE)
  # tbl[, year := design_$sim_prm$init_year] # given init year is 13 and close enough to 11 which I have epi for
  # col_nam <- setdiff(names(tbl), names(dt))
  # absorb_dt(dt, tbl)
  # dt[colon_ca_prvl == 1L, colon_ca_prvl := 2L + rpois(.N, clamp(duration - 2, 0, 100))]
  # dt[colon_ca_prvl > design_$sim_prm$cancer_cure,
  #    colon_ca_prvl := dqsample(2:design_$sim_prm$cancer_cure, .N, TRUE)]
  # dt[colon_ca_prvl > age, colon_ca_prvl := age]
  # dt[, (col_nam) := NULL]
  # dt[, colon_ca_prvl := carry_backward(colon_ca_prvl, pid_mrk)] # Not strictly necessary


  # Generate lung_cancer prevalence and duration (BI, PO) ----
  # simulate_init_prvl(mc, "lung_ca", design_, dt)
  #
  # # duration
  # tbl <- get_disease_epi_mc(mc, "lung_ca", "d", "m", FALSE)
  # tbl[, year := design_$sim_prm$init_year] # given init year is 13 and close enough to 11 which I have epi for
  # col_nam <- setdiff(names(tbl), names(dt))
  # absorb_dt(dt, tbl)
  # dt[lung_ca_prvl == 1L, lung_ca_prvl := 2L + rpois(.N, clamp(duration - 2, 0, 100))]
  # dt[lung_ca_prvl > design_$sim_prm$cancer_cure,
  #    lung_ca_prvl := dqsample(2:design_$sim_prm$cancer_cure, .N, TRUE)]
  # dt[lung_ca_prvl > age, lung_ca_prvl := age]
  # dt[, (col_nam) := NULL]
  # dt[, lung_ca_prvl := carry_backward(lung_ca_prvl, pid_mrk)] # Not strictly necessary
  #
  if (timing) message(proc.time() - ptm, units = "sec")
}
