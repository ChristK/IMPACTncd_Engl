## IMPACTncdEngland is an implementation of the IMPACTncd framework, developed by Chris
## Kypridemos with contributions from Peter Crowther (Melandra Ltd), Maria
## Guzman-Castillo, Amandine Robert, and Piotr Bandosz.
##
## Copyright (C) 2018-2026 University of Liverpool, Chris Kypridemos
##
## IMPACTncdEngland is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version. This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details. You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/> or write
## to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
## Boston, MA 02110-1301 USA.

# -----------------------------------------------------------------------------
# Simulation class calibration methods
# This file adds calibration methods to the Simulation class using $set()
# -----------------------------------------------------------------------------


#' @description
#' Calibrate incidence and case fatality rates.
#'
#' Reconstructs large input files, runs a Monte Carlo set, then iteratively
#' adjusts disease-specific calibration factors (`*_incd_clbr_fctr`,
#' `*_ftlt_clbr_fctr`) so simulated incidence and mortality match observed
#' targets. Writes results to `./simulation/calibration_prms.csv`.
#'
#' @param mc Integer vector of Monte Carlo iterations to run for calibration.
#' @param replace Logical. If `TRUE`, overwrite the existing
#'   `calibration_prms.csv`; otherwise extend it.
#' @return The `Simulation` object, invisibly.
Simulation$set("public", "calibrate_incd_ftlt", function(mc, replace = FALSE) {
  # recombine the chunks of large files
  self$reconstruct_large_files()

  export_xps <- self$design$sim_prm$export_xps # save the original value to be restored later
  self$design$sim_prm$export_xps <- FALSE # turn off export_xps to speed up the calibration
  private$create_empty_calibration_prms_file(replace = replace)
  clbr <- fread(
    "./simulation/calibration_prms.csv",
    colClasses = list(
      numeric = c(
        "chd_incd_clbr_fctr",
        "stroke_incd_clbr_fctr",
        "chd_ftlt_clbr_fctr",
        "stroke_ftlt_clbr_fctr",
        "nonmodelled_ftlt_clbr_fctr"
      )
    )
  )

  memedian <- function(x) {
    out <- median(x)
    if (out == 0L) {
      # For rare events, consider alternative estimators
      nonzero_vals <- x[x > 0]
      if (length(nonzero_vals) > 0) {
        prob_occurrence <- length(nonzero_vals) / length(x)
        median_nonzero <- median(nonzero_vals)
        out <- prob_occurrence * median_nonzero
      } else {
        out <- 0
      }
    }
    out
  }

  if (replace) {
    age_start <- self$design$sim_prm$ageL
  } else {
    # if replace == FALSE
    # if all ages exist skip calibration
    if (
      dim(clbr[
        chd_incd_clbr_fctr == 1 |
          stroke_incd_clbr_fctr == 1 |
          chd_ftlt_clbr_fctr == 1 |
          stroke_ftlt_clbr_fctr == 1 |
          nonmodelled_ftlt_clbr_fctr == 1
      ])[1] ==
        0
    ) {
      if (self$design$sim_prm$logs) {
        message("All ages have been calibrated. Skipping calibration.")
      }
      return(invisible(self))
    }
    age_start <- clbr[
      chd_incd_clbr_fctr == 1 |
        stroke_incd_clbr_fctr == 1 |
        chd_ftlt_clbr_fctr == 1 |
        stroke_ftlt_clbr_fctr == 1 |
        nonmodelled_ftlt_clbr_fctr == 1,
      min(age)
    ]
    if (self$design$sim_prm$logs) {
      message("Starting calibration from age ", age_start, ".")
    }
  }

  # Run the simulation from min to max age
  for (age_ in age_start:self$design$sim_prm$ageH) {
    # Run the simulation and export summaries
    self$del_logs()$del_outputs()$run(
      mc,
      multicore = TRUE,
      "sc0"
    )$export_summaries(
      multicore = TRUE,
      type = c("incd", "prvl", "dis_mrtl"),
      single_year_of_age = TRUE
    )

    # Incidence calibration
    # load the uncalibrated results
    unclbr <- open_dataset(file.path(
      self$design$sim_prm$output_dir,
      "summaries",
      "incd_scaled_up"
    )) %>%
      filter(age == age_) %>%
      select(
        "year",
        "age",
        "sex",
        "mc",
        "popsize",
        "chd_incd",
        "stroke_incd"
      ) %>%
      collect()
    setDT(unclbr)

    unclbr <- unclbr[,
      .(
        age,
        sex,
        year,
        mc,
        chd_incd = chd_incd / popsize,
        stroke_incd = stroke_incd / popsize
      )
    ][,
      .(
        chd_incd = memedian(chd_incd),
        stroke_incd = memedian(stroke_incd)
      ),
      keyby = .(age, sex, year)
    ]

    # for CHD
    # fit a log-log linear model to the uncalibrated results and store the coefficients
    tt <- unclbr[
      chd_incd > 0,
      as.list(coef(lm(
        log(chd_incd) ~ log(year)
      ))),
      by = sex
    ]

    unclbr[
      tt,
      on = "sex",
      c("intercept_unclbr", "trend_unclbr") := .(
        `(Intercept)`,
        `log(year)`
      )
    ]
    rm(tt)

    # load benchmark
    benchmark <- read_fst(
      file.path("./inputs/disease_burden", "chd_incd.fst"),
      columns = c("age", "sex", "year", "mu"),
      as.data.table = TRUE
    )[age == age_, ]
    # fit a log-log linear model to the benchmark incidence and store the coefficients
    benchmark[
      year >= self$design$sim_prm$init_year_long,
      c("intercept_bnchmrk", "trend_bnchmrk") := as.list(coef(lm(
        log(mu) ~ log(year)
      ))),
      by = sex
    ]

    # calculate the calibration factors
    unclbr[
      benchmark[year == max(year)],
      chd_incd_clbr_fctr := exp(
        intercept_bnchmrk + trend_bnchmrk * log(year)
      ) /
        exp(intercept_unclbr + trend_unclbr * log(year)),
      on = c("age", "sex")
    ]
    unclbr[, c("intercept_unclbr", "trend_unclbr") := NULL]

    # Repeat for stroke
    tt <- unclbr[
      stroke_incd > 0,
      as.list(coef(lm(
        log(stroke_incd) ~ log(year)
      ))),
      by = sex
    ]

    unclbr[
      tt,
      on = "sex",
      c("intercept_unclbr", "trend_unclbr") := .(
        `(Intercept)`,
        `log(year)`
      )
    ]
    rm(tt)

    benchmark <- read_fst(
      file.path("./inputs/disease_burden", "stroke_incd.fst"),
      columns = c("age", "sex", "year", "mu"),
      as.data.table = TRUE
    )[age == age_, ]
    benchmark[
      year >= self$design$sim_prm$init_year_long,
      c(
        "intercept_bnchmrk",
        "trend_bnchmrk"
      ) := as.list(coef(lm(
        log(mu) ~ log(year)
      ))),
      by = sex
    ]

    unclbr[
      benchmark[year == max(year)],
      stroke_incd_clbr_fctr := exp(
        intercept_bnchmrk + trend_bnchmrk * log(year)
      ) /
        exp(intercept_unclbr + trend_unclbr * log(year)),
      on = c("age", "sex")
    ]

    unclbr[, `:=`(
      chd_prvl_correction = chd_incd * (chd_incd_clbr_fctr - 1),
      stroke_prvl_correction = stroke_incd * (stroke_incd_clbr_fctr - 1),
      chd_incd = NULL,
      stroke_incd = NULL,
      intercept_unclbr = NULL,
      trend_unclbr = NULL
    )]
    clbr[
      unclbr,
      on = c("year", "age", "sex"),
      `:=`(
        chd_incd_clbr_fctr = i.chd_incd_clbr_fctr,
        stroke_incd_clbr_fctr = i.stroke_incd_clbr_fctr
      )
    ]

    # Case fatality calibration
    # Because we do incd and case fatality correction in the same step, we
    # need to estimate the expected changes on prvl because of the incd
    # calibration, before we proceed with the case fatality calibration.
    # Note that the calibration factor (multiplier) is 1/prvl as we
    # currently have mortality rates in the ftlt files.
    prvl <- open_dataset(file.path(
      self$design$sim_prm$output_dir,
      "summaries",
      "prvl_scaled_up"
    )) %>%
      filter(age == age_) %>%
      select(
        "year",
        "age",
        "sex",
        "mc",
        "popsize",
        "chd_prvl",
        "stroke_prvl"
      ) %>%
      collect()
    setDT(prvl)

    prvl <- prvl[, .(
      chd_prvl = chd_prvl / popsize,
      stroke_prvl = stroke_prvl / popsize,
      popsize,
      age,
      sex,
      year,
      mc
    )][,
      .(
        chd_prvl = memedian(chd_prvl),
        stroke_prvl = memedian(stroke_prvl),
        popsize = memedian(popsize)
      ),
      keyby = .(age, sex, year)
    ]
    prvl[
      unclbr,
      on = c("year", "age", "sex"),
      `:=`(
        chd_prvl_correction = i.chd_prvl_correction,
        stroke_prvl_correction = i.stroke_prvl_correction
      )
    ]
    benchmark <- read_fst(
      file.path("./inputs/disease_burden", "chd_ftlt.fst"),
      columns = c("age", "sex", "year", "mu2"),
      as.data.table = TRUE
    )[age == age_, ]
    prvl[benchmark, on = c("age", "sex", "year"), chd_mrtl := mu2]
    benchmark <- read_fst(
      file.path("./inputs/disease_burden", "stroke_ftlt.fst"),
      columns = c("age", "sex", "year", "mu2"),
      as.data.table = TRUE
    )[age == age_, ]
    prvl[benchmark, on = c("age", "sex", "year"), stroke_mrtl := mu2]
    benchmark <- read_fst(
      file.path("./inputs/disease_burden", "nonmodelled_ftlt.fst"),
      columns = c("age", "sex", "year", "mu2"),
      as.data.table = TRUE
    )[age == age_, ]
    prvl[benchmark, on = c("age", "sex", "year"), nonmodelled_mrtl := mu2]

    prvl[, `:=`(
      chd_ftlt_clbr_fctr = 1 / (chd_prvl + chd_prvl_correction),
      stroke_ftlt_clbr_fctr = 1 / (stroke_prvl + stroke_prvl_correction),
      nonmodelled_ftlt_clbr_fctr = 1 / (1 - chd_mrtl - stroke_mrtl)
    )]

    # Fix the calibration factors for the ages that have been calibrated
    if (age_ > age_start) {
      mrtl <- open_dataset(file.path(
        self$design$sim_prm$output_dir,
        "summaries",
        "dis_mrtl_scaled_up"
      )) %>%
        filter(age == age_ - 1L) %>%
        select(
          "year",
          "age",
          "sex",
          "mc",
          "popsize",
          "chd_deaths",
          "stroke_deaths",
          "nonmodelled_deaths"
        ) %>%
        collect()
      setDT(mrtl)

      mrtl <- mrtl[, .(
        chd_mrtl = chd_deaths / popsize,
        stroke_mrtl = stroke_deaths / popsize,
        nonmodelled_mrtl = nonmodelled_deaths / popsize,
        popsize,
        age,
        sex,
        year,
        mc
      )][,
        .(
          chd_mrtl = memedian(chd_mrtl),
          stroke_mrtl = memedian(stroke_mrtl),
          nonmodelled_mrtl = memedian(nonmodelled_mrtl),
          popsize = memedian(popsize)
        ),
        keyby = .(age, sex, year)
      ]
      benchmark <- read_fst(
        file.path("./inputs/disease_burden", "chd_ftlt.fst"),
        columns = c("age", "sex", "year", "mu2"),
        as.data.table = TRUE
      )[age == age_ - 1L, ]
      mrtl[
        benchmark,
        on = c("age", "sex", "year"),
        chd_ftlt_clbr_fctr := mu2 / chd_mrtl
      ]
      benchmark <- read_fst(
        file.path("./inputs/disease_burden", "stroke_ftlt.fst"),
        columns = c("age", "sex", "year", "mu2"),
        as.data.table = TRUE
      )[age == age_ - 1L, ]
      mrtl[
        benchmark,
        on = c("age", "sex", "year"),
        stroke_ftlt_clbr_fctr := mu2 / stroke_mrtl
      ]
      benchmark <- read_fst(
        file.path("./inputs/disease_burden", "nonmodelled_ftlt.fst"),
        columns = c("age", "sex", "year", "mu2"),
        as.data.table = TRUE
      )[age == age_ - 1L, ]
      mrtl[
        benchmark,
        on = c("age", "sex", "year"),
        nonmodelled_ftlt_clbr_fctr := mu2 / nonmodelled_mrtl
      ]
      mrtl[chd_ftlt_clbr_fctr == Inf, chd_ftlt_clbr_fctr := 1] # to avoid Inf through division by 0
      mrtl[stroke_ftlt_clbr_fctr == Inf, stroke_ftlt_clbr_fctr := 1] # to avoid Inf through division by 0
      mrtl[
        nonmodelled_ftlt_clbr_fctr == Inf,
        nonmodelled_ftlt_clbr_fctr := 1
      ] # to avoid Inf through division by 0

      clbr[
        mrtl,
        on = c("year", "age", "sex"),
        `:=`(
          chd_ftlt_clbr_fctr = i.chd_ftlt_clbr_fctr * chd_ftlt_clbr_fctr,
          stroke_ftlt_clbr_fctr = i.stroke_ftlt_clbr_fctr *
            stroke_ftlt_clbr_fctr,
          nonmodelled_ftlt_clbr_fctr = i.nonmodelled_ftlt_clbr_fctr *
            nonmodelled_ftlt_clbr_fctr
        )
      ]
    }

    if (age_ == self$design$sim_prm$ageH) {
      # shortcut for age == 99 hopefully with tiny bias
      mrtl[, age := age + 1L]
      prvl[
        mrtl,
        on = c("year", "age", "sex"),
        `:=`(
          chd_ftlt_clbr_fctr = i.chd_ftlt_clbr_fctr * chd_ftlt_clbr_fctr,
          stroke_ftlt_clbr_fctr = i.stroke_ftlt_clbr_fctr *
            stroke_ftlt_clbr_fctr,
          nonmodelled_ftlt_clbr_fctr = i.nonmodelled_ftlt_clbr_fctr *
            nonmodelled_ftlt_clbr_fctr
        )
      ]
    }

    clbr[
      prvl,
      on = c("year", "age", "sex"),
      `:=`(
        chd_ftlt_clbr_fctr = i.chd_ftlt_clbr_fctr,
        stroke_ftlt_clbr_fctr = i.stroke_ftlt_clbr_fctr,
        nonmodelled_ftlt_clbr_fctr = i.nonmodelled_ftlt_clbr_fctr
      )
    ]

    fwrite(clbr, "./simulation/calibration_prms.csv") # NOTE this needs to be inside the loop so it influences the simulation during the loop over ages
  } # end loop over ages

  self$design$sim_prm$export_xps <- export_xps # restore the original value
  invisible(self)
})


# create_empty_calibration_prms_file ----
# if replace is FALSE then it creates a calibration parameters when it is
# missing, file filed with 1. If replace = TRUE it overwrites the existin
# file
# returns invisible(self)
# TODO Automate based on diseases in design.yaml
Simulation$set("private", "create_empty_calibration_prms_file", function(replace = FALSE) {
  if (replace || !file.exists("./simulation/calibration_prms.csv")) {
    clbr <- CJ(
      year = self$design$sim_prm$init_year_long:self$design$sim_prm$sim_horizon_max,
      age = self$design$sim_prm$ageL:self$design$sim_prm$ageH,
      sex = c("men", "women"),
      chd_incd_clbr_fctr = 1,
      stroke_incd_clbr_fctr = 1,
      chd_ftlt_clbr_fctr = 1,
      stroke_ftlt_clbr_fctr = 1,
      nonmodelled_ftlt_clbr_fctr = 1
    )
    fwrite(clbr, "./simulation/calibration_prms.csv")
  }
  invisible(self)
})
