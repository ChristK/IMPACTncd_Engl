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
most_cost_effective <- function(dt) {
  names(dt[year == max(year), .(nmb_cml, mc, friendly_name)
  ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
  ][, friendly_name[which.max(nmb_cml)], by = mc
  ][, first(sort(table(V1), decreasing = TRUE))])
}

#' @export
rank_cost_effective <- function(dt) {
  names(dt[year == max(year), .(nmb_cml, mc, friendly_name)][,
    sum_dt(.SD, c("mc", "friendly_name"), character(0))][,
      friendly_name[which.max(nmb_cml)], by = mc][,
        sort(table(V1), decreasing = TRUE)])
}

#' @export
most_effective <- function(dt) {
  names(dt[year == max(year),
    .(eq5d_cml, mc, friendly_name)
  ][, friendly_name[which.max(eq5d_cml)], by = mc
  ][, first(sort(table(V1), decreasing = TRUE))])
}

#' @export
# based on absolute equity
most_equitable <- function(dt) {
  tt <- dt[year == max(year), .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc, friendly_name, qimd)
  ][, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name"))
  calc_sei(tt, c("mc", "friendly_name"))
  names(tt[, .(nmb_cml = sum(nmb_cml), sei = first(sei)), by = c("mc", "friendly_name")][, friendly_name[which.max(sei)], by = mc][, head(sort(table(V1), decreasing = TRUE), 1L)])
}

#' @export
most_equitable_rel <- function(dt) {
  tt <-
    dt[year == max(year), .(net_utility_cml,
      nmb_cml, eq5d_cml, pops_cml, mc, friendly_name, qimd)
      ][, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name"))
  calc_sei(tt, c("mc", "friendly_name"))
  calc_rei(tt, c("mc", "friendly_name"))
  names(tt[, .(
    nmb_cml = sum(nmb_cml),
    sei = first(sei),
    rei = first(rei)
  ), by = c("mc", "friendly_name")
    ][, .(abs = friendly_name[which.max(sei)],
    rel = friendly_name[which.max(rei)]),
    by = mc
      ][rel == abs, head(sort(table(rel), decreasing = TRUE), 1L)])
}

#' @export
most_equitable_abs <- function(dt) {
  tt <- dt[year == max(year), .(net_utility_cml, nmb_cml, eq5d_cml, pops_cml, mc, friendly_name, qimd)
  ][, sum_dt(.SD, c("mc", "friendly_name", "qimd"), character(0))]
  calc_rigit_scores(tt, c("mc", "friendly_name"))
  calc_sei(tt, c("mc", "friendly_name"))
  names(tt[, .(nmb_cml = sum(nmb_cml), sei = first(sei)), by = c("mc", "friendly_name")][, friendly_name[which.max(sei)], by = mc][, head(sort(table(V1), decreasing = TRUE), 1L)])
}

#' @export
incr_cost_effect_ratio <- function() " '£X, fill me!' "

tot_net_monet_benef <- function(dt) {
  dt[year == max(year),
    .(
      nmb_cml = sum(nmb_cml)
    ), by = .(mc, friendly_name)
  ][, mean(nmb_cml), by = friendly_name][, signif(first(sort(V1, decreasing = TRUE)), 2)]
}

#' @export
reduce_rel_index_ineq <- function() " 'X, fill me!' "   # TODO

#' @export
increase_abs_index_ineq <- function() " 'X, fill me!' " # TODO

#' @export
social_care_cost_sav <- function(dt) {
  dt[year == max(year), signif(median(net_socialcare_cost_cml), 2), keyby = friendly_name
  ]
}

#' @export
prod_benef <- function(dt) {
  dt[year == max(year), signif(median(net_productivity_cost_cml), 2), keyby = friendly_name
  ]
}

#' @export
inform_care_cost_sav <- function(dt) {
  dt[year == max(year), signif(median(net_informal_care_cost_cml), 2), keyby = friendly_name
    ]
}

#' @export
benefit_cost_ratio_cml <-
  function(dt,
    perspective = c("Societal perspective",
      "Health and social care perspective",
      "Healthcare perspective"),
    wtp = input$out_wtp_box) {
    # perspective input$health_econ_perspective_checkbox
    # wtp input$out_wtp_box
    if (perspective == "Societal perspective") {
      out <- dt[, .(bcr_cml = (
        -(
          net_healthcare_cost_cml + net_socialcare_cost_cml + net_informal_care_cost_cml - net_productivity_cost_cml
        ) + (net_utility_cml * wtp)
      ) / net_policy_cost_cml)]
    } else if (perspective == "Health and social care perspective") {
      out <- dt[, .(bcr_cml = (
        -(
          net_healthcare_cost_cml + net_socialcare_cost_cml
        ) + (net_utility_cml * wtp)
      ) / net_policy_cost_cml)]
    } else if (perspective == "Healthcare perspective") {
      out <- dt[, .(bcr_cml = (
        -(
          net_healthcare_cost_cml
        ) + (net_utility_cml * wtp)
      ) / net_policy_cost_cml)]
    }
    invisible(out)
  }

#' @export
most_benefit_cost_ratio <- function(dt,
  perspective = c("Societal perspective",
    "Health and social care perspective",
    "Healthcare perspective"),
  wtp = input$out_wtp_box,
  how_many = 1L) {

  names(dt[year == max(year), .(net_healthcare_cost_cml, net_socialcare_cost_cml,
    net_informal_care_cost_cml, net_productivity_cost_cml,
    net_utility_cml, net_policy_cost_cml,
    mc, friendly_name)
  ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
  ][, bcr_cml := benefit_cost_ratio_cml(.SD, perspective = perspective, wtp = wtp)
  ][, friendly_name[which.max(bcr_cml)], by = mc
  ][, head(sort(table(V1), decreasing = TRUE), how_many)])
}

#' @export
benefit_cost_ratio_value <- function(dt,
  perspective = c("Societal perspective",
    "Health and social care perspective",
    "Healthcare perspective"),
  wtp = input$out_wtp_box,
  friendly_nam) {

  dt[year == max(year), .(net_healthcare_cost_cml, net_socialcare_cost_cml,
    net_informal_care_cost_cml, net_productivity_cost_cml,
    net_utility_cml, net_policy_cost_cml,
    mc, friendly_name)
  ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
  ][, bcr_cml := benefit_cost_ratio_cml(.SD, perspective = perspective, wtp = wtp)
  ][, median(bcr_cml), by = friendly_name
  ][, setorder(na.omit(.SD), -V1)][friendly_name %in% friendly_nam, V1]
}

#' @export
rank_benefit_cost_ratio <- function(dt,
  perspective = c("Societal perspective",
    "Health and social care perspective",
    "Healthcare perspective"),
  wtp = input$out_wtp_box) {

  names(dt[year == max(year), .(net_healthcare_cost_cml, net_socialcare_cost_cml,
    net_informal_care_cost_cml, net_productivity_cost_cml,
    net_utility_cml, net_policy_cost_cml,
    mc, friendly_name)
  ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
  ][, bcr_cml := benefit_cost_ratio_cml(.SD, perspective = perspective, wtp = wtp)
  ]
    [, friendly_name[which.max(bcr_cml)], by = mc][, sort(table(V1), decreasing = TRUE)])
}

#' @export
icer_cml <-
  function(dt,
    perspective = c("Societal perspective",
      "Health and social care perspective",
      "Healthcare perspective")) {
    # perspective input$health_econ_perspective_checkbox
    # wtp input$out_wtp_box
    if (perspective == "Societal perspective") {
      out <- dt[, .(icer_cml = round(net_policy_cost_cml + net_healthcare_cost_cml +
          net_socialcare_cost_cml + net_informal_care_cost_cml -
          net_productivity_cost_cml) / net_utility_cml)]
    } else if (perspective == "Health and social care perspective") {
      out <- dt[, .(icer_cml = round(net_policy_cost_cml + net_healthcare_cost_cml +
          net_socialcare_cost_cml) / net_utility_cml)]
    } else if (perspective == "Healthcare perspective") {
      out <- dt[, .(icer_cml = round(net_policy_cost_cml + net_healthcare_cost_cml) / net_utility_cml)]
    }
    invisible(out)
  }

#' @export
scn_icer_cml <-
  function(dt,
    friendly_name_,
    perspective = c("Societal perspective",
      "Health and social care perspective",
      "Healthcare perspective")) {
    dt[year == max(year) & friendly_name == friendly_name_,
      .(
        net_healthcare_cost_cml,
        net_socialcare_cost_cml,
        net_informal_care_cost_cml,
        net_productivity_cost_cml,
        net_utility_cml,
        net_policy_cost_cml,
        mc
      )][, sum_dt(.SD, c("mc"), character(0))
      ][, icer_cml := icer_cml(.SD, perspective = perspective)
      ][, mean(icer_cml, na.rm = TRUE)]
  }

#' @export
net_monetary_benefit_cml <-
  function(dt,
    perspective = c("Societal perspective",
      "Health and social care perspective",
      "Healthcare perspective"),
    wtp = input$out_wtp_box) {
    # perspective input$health_econ_perspective_checkbox <- "Societal perspective"
    # wtp input$out_wtp_box
    if (perspective == "Societal perspective") {
      out <-dt[, .(nmb_cml = (-(net_healthcare_cost_cml + net_socialcare_cost_cml + net_informal_care_cost_cml - net_productivity_cost_cml) + (net_utility_cml * wtp)) - net_policy_cost_cml)]
    } else if (perspective == "Health and social care perspective") {
      out <- dt[, .(nmb_cml = (-(net_healthcare_cost_cml + net_socialcare_cost_cml) + (net_utility_cml * wtp)) - net_policy_cost_cml)]
    } else if (perspective == "Healthcare perspective") {
      out <- dt[, .(nmb_cml = (-(net_healthcare_cost_cml) + (net_utility_cml * wtp)) - net_policy_cost_cml)]
    }
    invisible(out)
  }

#' @export
most_net_monetary_benefit_value <- function(dt,
  perspective = c("Societal perspective",
    "Health and social care perspective",
    "Healthcare perspective"),
  wtp = input$out_wtp_box,
  how_many = 1L) {

  names(dt[year == max(year), .(net_healthcare_cost_cml, net_socialcare_cost_cml,
    net_informal_care_cost_cml, net_productivity_cost_cml,
    net_utility_cml, net_policy_cost_cml,
    mc, friendly_name)
  ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
  ][, nmb_cml := net_monetary_benefit_cml(.SD, perspective = perspective, wtp = wtp)
  ][, friendly_name[which.max(nmb_cml)], by = mc
  ][, head(sort(table(V1), decreasing = TRUE), how_many)])
}

#' @export
net_monetary_benefit_value <- function(dt,
  perspective = c("Societal perspective",
    "Health and social care perspective",
    "Healthcare perspective"),
  wtp = input$out_wtp_box,
  friendly_nam) {

  dt[year == max(year), .(net_healthcare_cost_cml, net_socialcare_cost_cml,
    net_informal_care_cost_cml, net_productivity_cost_cml,
    net_utility_cml, net_policy_cost_cml,
    mc, friendly_name)
  ][, sum_dt(.SD, c("mc", "friendly_name"), character(0))
  ][, nmb_cml := net_monetary_benefit_cml(.SD, perspective = perspective, wtp = wtp)
  ][, median(nmb_cml), by = friendly_name
  ][, setorder(na.omit(.SD), -V1)][friendly_name %in% friendly_nam, V1]
}


#' @export
sum_dt <- function(dt, by = c("year", "friendly_name", "mc"), cols_to_exclude = c("scenario", "agegrp", "sex" , "qimd", "ethnicity")) {
  cols_to_exclude <- intersect(cols_to_exclude, names(dt))

  dt[, lapply(.SD, function(x) {
    if (is.numeric(x)) sum(x) else x
  }), keyby = by,
    .SDcols = -cols_to_exclude]
}

#' @export
median_dt <- function(dt, by = c("year", "friendly_name"), cols_to_exclude = c("scenario", "agegrp", "sex" , "qimd", "ethnicity", "mc"), digits = 0) {
  cols_to_exclude <- intersect(cols_to_exclude, names(dt))

  dt[, lapply(.SD, function(x) {
    if (is.numeric(x)) round(median(x), digits) else x
  }), keyby = by,
    .SDcols = -cols_to_exclude]
}

#' @export
mean_dt <- function(dt, by = c("year", "friendly_name"), cols_to_exclude = c("scenario", "agegrp", "sex" , "qimd", "ethnicity", "mc"), digits = 0) {
  cols_to_exclude <- intersect(cols_to_exclude, names(dt))

  dt[, lapply(.SD, function(x) {
    if (is.numeric(x)) round(mean(x), digits) else x
  }), keyby = by,
    .SDcols = -cols_to_exclude]
}

#' @export
# rigit scores
calc_rigit_scores <- function(dt, by = c("mc", "year", "friendly_name")) {
  keydt <- key(dt)
  keydt <- intersect(keydt, names(dt))
  setkey(dt, mc, qimd, friendly_name)
  dt[, qimd2 := cumsum(pops_cml)/sum(pops_cml),
    keyby = by]
  dt[, qimd2 := (qimd2 + c(0, qimd2[seq_len(.N - 1)]))/2,
    keyby = by]
  # dt[, qimd2_median := median(qimd2, na.rm = T), keyby = .(year, scenario, qimd)]
  if (length(keydt)) setkeyv(dt, keydt)
  dt
}

#' @export
# Absolute equity slope index
calc_sei <- function(dt, by = c("mc", "year", "friendly_name")) {
  dt[, sei := -lm(net_utility_cml~qimd2)$coefficients["qimd2"], by = by]
  # -lm(...)coefficients["qimd2"] because qimd 1 == most deprived
}

#' @export
# Relative equity slope index
calc_rei <- function(dt, by = c("mc", "year", "friendly_name")) {
  dt[, proport_utility_cml := net_utility_cml/eq5d_cml]
  dt[, rei := -lm(proport_utility_cml~qimd2)$coefficients["qimd2"], by = by]
}
