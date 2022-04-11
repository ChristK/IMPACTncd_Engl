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


`:=` = function(...)
  NULL # due to NSE notes in R CMD check

.onUnload <- function(libpath) {
  library.dynam.unload("IMPACTncdEnglmisc", libpath)
}

# Ensures that when fwrite appends file colnames of file to be written, match
# those already in the file
#' @export
fwrite_safe <- function(x,
                        file = "",
                        append = TRUE,
                        ...) {
  if (append) {
    if (file.exists(file)) {
      col_names_disk <- names(fread(file, nrows = 0))
      col_names_file <- names(x)
      col_names <- outersect(col_names_disk, col_names_file)
      if (length(col_names) > 0)
        x[, (col_names) := NA]
      setcolorder(x, col_names_disk)
    }
  }
  fwrite(x, file, append, ...)
}


#' @export
inflate <- function(x, percentage_rate, year, baseline_year) {
  x * (1 + percentage_rate / 100) ^ (year - baseline_year)
}
# inflate(1000, 3, 2011:2020, 2013)

#' @export
deflate <- function(x, percentage_rate, year, baseline_year) {
  x * (1 - percentage_rate / 100) ^ (year - baseline_year)
}

# Necessary aux functions
# plots mean x by y stratified by z
#' @export
plot_cor <- function(x, y, z, wt, dt) {
  xlab  <- toupper(x)
  ylab  <- toupper(y)
  title <- paste0(ylab, "~", xlab, "|", toupper(z))

  dt[, .(y = weighted.mean(get(y), get(wt), na.rm = TRUE)),
     keyby = .(x = round(get(x)), z = get(z))][, qplot(
    x,
    y,
    data = na.omit(.SD),
    facets = as.formula("~z"),
    geom = "smooth",
    main = title,
    xlab = xlab,
    ylab = ylab
  )]
}

# helper func for gamlss::fitDistr
#' @export
distr_best_fit <-
  function(dt,
           var,
           wt,
           distr_family,
           distr_extra = NULL,
           pred = FALSE,
           seed = NULL,
           trace = TRUE) {
    if (pred) {
      print("Selection based on minimum prediction global deviance")
      if (!is.null(seed))
        set.seed(seed)
      lns <- sample(nrow(dt), round(nrow(dt) * 0.8))
      dt_trn   <- dt[lns,] # train dataset
      dt_crv   <- dt[!lns,]  # cross-validation dataset
      marg_distr <- gamlss::fitDistPred(
        dt_trn[[var]],
        type = distr_family,
        weights = dt_trn[[wt]],
        extra = distr_extra,
        try.gamlss = TRUE,
        trace = trace,
        newdata = dt_crv[[var]]
      )
    } else {
      print("Selection based on BIC")
      marg_distr <-
        gamlss::fitDist(
          dt[[var]],
          log(nrow(dt)),
          type = distr_family,
          weights = dt[[wt]],
          extra = distr_extra,
          try.gamlss = TRUE,
          trace = trace
        )
    }
    marg_distr
  }

#' @export
distr_validation <-
  function(marg_distr,
           dt,
           title,
           discrete = FALSE,
           smooth = 0.35) {
    params <- vector("list", length(marg_distr$parameters) + 1L)
    names(params) <- c(marg_distr$parameters, "p")
    for (i in marg_distr$parameters)
      params[[i]] <- get(i, marg_distr)

    params$p <- seq(0.001, 0.999, length.out = nrow(dt))
    distr_nam <- marg_distr$family[[1]]
    y <- do.call(paste0("q", distr_nam), params)
    y <- y[between(y, dt[, min(var)], dt[, max(var)])]
    y_wt <- rep(1 / length(y), length(y))
    # validation plots
    # see http://www.csss.washington.edu/files/working-papers/2002/wp27.pdf
    out <- reldist_diagnostics(
      dt$var,
      y,
      dt[, wt / sum(wt)],
      y_wt,
      main = title,
      discrete = discrete,
      smooth = smooth
    )
    out
  }

#' @export
centile_predictAll <-
  function(gamlss_model, orig_data, newdata, cent = 0.5) {
    stopifnot("gamlss" %in% class(gamlss_model))
    tt <- predictAll(gamlss_model, newdata, data = orig_data)
    tt$p <- cent
    distr_nam <- attr(tt, "family")[[1]]
    out <- do.call(paste0("q", distr_nam), tt)
    return(out)
  }

#' @export
mean_predictAll <- function(gamlss_model, orig_data, newdata) {
  stopifnot("gamlss" %in% class(gamlss_model))
  tt <- predictAll(gamlss_model, newdata, data = orig_data)
  setDT(tt)
  tt[, n := 1e5]
  setcolorder(tt, "n")
  distr_nam <- attr(tt, "family")[[1]]
  out <- apply(tt, 1, function(x) {
    mean(do.call(paste0("r", distr_nam), as.list(x)))
  })
  return(out)
}

#' @export
plot_synthpop_val <-
  function(dt,
           x,
           grp,
           wt,
           title,
           x_label,
           standardised_to_grp = c(FALSE, TRUE),
           print_to_screen = c(FALSE, TRUE)) {
    if (standardised_to_grp) {
      dt[, weight := get(wt) / sum(get(wt)), by = c("type", grp)]
    } else {
      dt[, weight := get(wt) / sum(get(wt)), by = c("type")]
    }

    if ("year" %in% grp)
      dt[, year := year + 2000L]
    x <- enquo(x)
    p <- ggplot(dt, aes(
      !!x,
      colour = type,
      weight = weight,
      linetype = type
    )) +
      geom_density() +
      facet_grid(grp) +
      xlab(x_label) +
      ggtitle(title)
    if (print_to_screen)
      print(p)
    suppressWarnings(
      cowplot::ggsave2(
        paste0(gsub(" ", "_", title), "_density.png"),
        p,
        width = 16,
        height = 9,
        units = "cm",
        scale = 2,
        dpi = 300,
        path = "./validation/synthpop_models"
      )
    )

    p <- ggplot(dt, aes(
      !!x,
      colour = type,
      weight = weight,
      linetype = type
    )) +
      stat_ecdf() +
      facet_grid(grp) +
      xlab(x_label) +
      ggtitle(title)
    if (print_to_screen)
      print(p)

    if ("year" %in% grp)
      dt[, year := year - 2000L]

    suppressWarnings(
      cowplot::ggsave2(
        paste0(gsub(" ", "_", title), "_cdf.png"),
        p,
        width = 16,
        height = 9,
        units = "cm",
        scale = 2,
        dpi = 300,
        path = "./validation/synthpop_models"
      )
    )
  }

#' @export
validate_gamlss_tbl <-
  function(dt,
           gamlss_tbl,
           mc = 10L,
           colname,
           distr_nam = distr_nam) {
    stopifnot(is.data.table(dt), is.data.table(gamlss_tbl), mc >= 1)
    nam_var <- intersect(names(dt), names(gamlss_tbl))
    nam_param <- setdiff(names(gamlss_tbl), nam_var)
    dt[, age := as.integer(age)]
    x <- copy(dt)
    x[, `:=`(type, "Observed")]
    z <- copy(dt)
    z[gamlss_tbl, (nam_param) := mget(nam_param), on = nam_var]
    if (z[is.na(mu), .N] > 0)
      stop("NAs produced in join")
    z[, `:=`(type, "Modelled")]
    z <- rbindlist(rep(list(z), mc))
    z[, p := dqrunif(.N, 0, 0.999)]
    z[, (colname) := do.call((distr_nam), .SD), .SDcols = c("p", nam_param)]
    z[, c("p", (nam_param)) := NULL]
    out <- rbind(x, z, use.names = TRUE, fill = TRUE)
  }

# dt <- dt[age >= 20L]
# gamlss_tbl <- copy(smok_incid_model_tbl)
# colname <- "smok_incid"
# z[, .SD, .SDcols = c("p", nam_param)]
# z[, do.call("qBI", .SD), .SDcols = c("p", nam_param)]

#' @export
shift_bypid <-
  function(x, lag, id, replace = NA) {
    if (typeof(x) == "integer") {
      return(shift_bypidInt(x, lag, replace, id))
    } else if (typeof(x) == "logical") {
      return(shift_bypidBool(x, lag, replace, id))
    } else if (typeof(x) == "double") {
      return(shift_bypidNum(x, lag, replace, id))
    } else
      stop("type of x not supported")
  }

#' @export
lung_ca_rr <- function(smok_status = 3,
                       num_cig = 3,
                       smok_dur = 3,
                       quit_yrs = 40,
                       age = 62,
                       education = 4L,
                       ethnicity = 1L,
                       bmi = 27,
                       copd = FALSE,
                       personal_history_of_cancer = FALSE,
                       family_history_of_lung_cancer = FALSE) {
  ethn_vec <- c(
    0,
    # "white"-0.5241286,
    # "indian"-0.5241286,
    # "pakistani"-0.5241286,
    # "bangladeshi"-0.5241286,
    # "other asian"
    0.3211605,
    # "black caribbean"
    0.3211605,
    # "black african"-0.5241286,
    # "chinese"
    0          # "other"

  )
  # Original education in PLCO
  # 1=less than high school completed; 2=high school graduate;
  # 3=post high school training; 4=some college; 5=college graduate;
  # 6=postgraduate or professional degree
  edu_vec <- c(5L, 5L, 2L, 2L, 1L, 1L, 1L)
  smok_vec <- c(0, 2.542472, 2.542472, 2.799727)

  # prb_a is the probability of lung cancer if not a smoker (assumes also no
  # COPD)
  a_nocopd <- -7.02198 + 0.079597 * (age - 62) -
    0.0879289 * (edu_vec[education] - 4) -
    0.028948 * (bmi - 27) + ethn_vec[ethnicity] +
    0.4845208 * personal_history_of_cancer +
    0.5856777 * family_history_of_lung_cancer

  # a_copd       <- a_nocopd + 0.3457265
  # needs special treatment as it depends on another disease
  # a            <- a_nocopd + 0.3457265 * copd
  exp_a_nocopd <- exp(a_nocopd)
  # exp_a_copd   <- exp(a_copd)
  # exp_a        <- exp(a)
  # TODO simplification as if not a smoker the likelihood of COPD is decreasing
  # but there is residual

  prb_a_nocopd <- exp_a_nocopd / (1 + exp_a_nocopd)
  # prb_a_copd <- exp_a_copd / (1 + exp_a_copd)
  # prb_a <- exp_a / (1 + exp_a)

  # prb_b is the probability of lung cancer for observed smoking status
  x <- fifelse(smok_status != 1L, {
    smok_vec[smok_status] +
      (-0.1815486 * (((num_cig / 100) ^ -1) - 4.021541613)) +
      0.0305566 * (smok_dur - 27)
  }, 0)

  x <-
    x - fifelse(smok_status %in% (2:3), 0.0321362 * (quit_yrs - 8.593417626), 0)

  b        <- a_nocopd + x + 0.3457265 * copd
  b_nocopd <- a_nocopd + x
  b_copd   <- a_nocopd + x + 0.3457265 # COPD coef

  exp_b_nocopd <- exp(b_nocopd)
  exp_b_copd   <- exp(b_copd)
  exp_b        <- exp(b)
  prb_b_nocopd <- exp_b_nocopd / (1 + exp_b_nocopd)
  prb_b_copd   <- exp_b_copd / (1 + exp_b_copd)
  prb_b        <- exp_b / (1 + exp_b)

  # Then RR is:
  rr <- list(
    "rr_for_parf"  = clamp(prb_b / prb_a_nocopd, 1, 20),
    # otherwise light ex smokers had RR < 1
    "rr_no_copd"   = clamp(prb_b_nocopd / prb_a_nocopd, 1, 20),
    "rr_with_copd" = clamp(prb_b_copd / prb_a_nocopd, 1, 20)
  )

  return(rr)
}

# get observed mortality from ONS data
#' @export
get_ons_mrtl <-
  function(disease,
           type = c("rate", "absolute"),
           agegrp_width = c(5L, 20L)) {
    if (type == "rate")
      prefix <- "Mx_"
    if (type == "absolute")
      prefix <- "deaths_"

    colnams <-
      c(
        "year",
        "sex",
        "qimd",
        paste0("agegrp", agegrp_width),
        paste0(prefix, disease),
        "pops"
      )
    if (agegrp_width == 5)
      tt <-
      read_fst(
        "./ONS_data/mortality_by_agegrp5_for_validation.fst",
        columns = colnams,
        as.data.table = TRUE
      )
    if (agegrp_width == 20)
      tt <-
      read_fst(
        "./ONS_data/mortality_by_agegrp20_for_validation.fst",
        columns = colnams,
        as.data.table = TRUE
      )
    tt
  }
# get_ons_mrtl("chd", "rate", 20)

# get observed cancer incidence from ONS data
#' @export
get_ons_incd <-
  function(disease,
           type = c("rate", "absolute"),
           agegrp_width = 20L) {
    if (type == "rate")
      prefix <- "rate_"
    if (type == "absolute")
      prefix <- "cases_"

    colnams <-
      c(
        "year",
        "sex",
        "qimd",
        paste0("agegrp", agegrp_width),
        paste0(prefix, disease),
        "pops"
      )
    if (agegrp_width == 5)
      tt <-
      read_fst(
        "./ONS_data/cancer_incd_by_agegrp5_for_validation.fst",
        columns = colnams,
        as.data.table = TRUE
      )
    if (agegrp_width == 20)
      tt <-
      read_fst(
        "./ONS_data/cancer_incd_by_agegrp20_for_validation.fst",
        columns = colnams,
        as.data.table = TRUE
      )
    tt
  }

# get disease epi parameters for mc
#' @export
get_disease_epi_mc <-
  function(mc,
           disease,
           epi_par = c("incidence", "prevalence", "fatality", "duration"),
           what = c("value", "median_value", "p"),
           # p is the percentile of the uncertainty
           stochastic = TRUE) {
    # stochastic = F returns median_value even if what == "value"
    stopifnot(between(mc, 1, 1e3))
    # TODO make 1e3 derived from disease_epi_indx.fst
    epi_par2 <- epi_par <- match.arg(epi_par)
    what <- match.arg(what)
    if (epi_par == "fatality") {
      epi_par <- paste0("case_", epi_par)
    }
    if (epi_par == "duration") {
      epi_par <-
        paste0(epi_par, "_years")
    } else {
      epi_par <- paste0(epi_par, "_rates")
    }
    if (!stochastic && what == "value") {
      what <- "median_value"
    }
    indx <-
      read_fst(
        "./disease_epidemiology/disease_epi_indx.fst",
        from = mc,
        to = mc,
        as.data.table = TRUE
      )

    if (epi_par2 != "duration") {
      colnam <-
        c("age",
          "sex",
          "qimd",
          paste0(what, "_", disease, "_", epi_par))
      out <-
        read_fst(
          "./disease_epidemiology/disease_epi_l.fst",
          columns = colnam,
          from = indx$from,
          to = indx$to,
          as.data.table = TRUE
        )
      setnames(out, tail(colnam, 1), epi_par2)

    } else {
      # if duration

      dur_colnam <- paste0(what, "_", disease, "_", epi_par)
      prev_colnam <-
        paste0(what, "_", disease, "_", "prevalence_rates")
      colnam <- c("age", "sex", "qimd", dur_colnam, prev_colnam)

      out <-
        read_fst(
          "./disease_epidemiology/disease_epi_l.fst",
          columns = colnam,
          from = indx$from,
          to = indx$to,
          as.data.table = TRUE
        )
      age_range <- out[, range(age)]
      setnames(out, prev_colnam, "prevalence_rates")
      setnames(out, dur_colnam, "duration")

      out <-
        out[prevalence_rates > 0] # To exclude men in breast ca. May have

      # get the prevalent cases. Each row is now a case
      # logic to get at least 10 cases per strata and prevent 0
      n <-
        out[, 10 / min(prevalence_rates)]
      out <- out[rep(seq_len(.N), prevalence_rates * n)]
      out[, `:=` (pid = .I - 1L, prevalence_rates = NULL)]
      # Now project the pids that each case will have copd
      out <- out[rep(seq_len(.N), duration)]

      out[, age := age + seq_len(.N) - 1L, by = pid]
      out[, disease_years := seq_len(.N), by = pid]
      out <-
        out[between(age, age_range[1], age_range[2]),
            .("duration" = mean(disease_years)), keyby = .(age, sex, qimd)]
    }

    invisible(out)
  }
# get_disease_epi_mc(1, "chd", "i", "v")

# get disease disease epi parameters for mc

# Extract causal pathways from RR list
#' @export
get_causal_paths <- function(RR) {
  exposures <- sapply(RR, `[[`, "name", USE.NAMES = FALSE)
  # tobacco_lung_ca RR is a function and handled separately
  exposures <- c(exposures, "tobacco")
  outcomes <- sapply(RR, `[[`, "outcome", USE.NAMES = FALSE)
  outcomes <- c(outcomes, "lung_ca")

  out <- list()
  out$dt <- data.table(exposures, outcomes, key = "exposures")
  out$by_outcome <- unstack(out$dt, exposures ~ outcomes)
  out$by_exps <- unstack(out$dt, outcomes ~ exposures)
  out$matrix <- as.matrix(table(exposures, outcomes))
  out
}

# causal_paths <- get_causal_paths(RR)
# causal_paths$by_outcome$chd

#' @export
get_lifetable_all <-
  function(mc, disease, design, type = c("qx", "mx")) {
    if (disease %in% c("allcause", "nonmodelled")) {
      disease2 <- "chd"
    } else {
      disease2 <- disease
    }
    prb <-
      get_disease_epi_mc(mc,
                         disease = disease2,
                         "fatality",
                         "p",
                         stochastic = design$stochastic)
    # NOTE fatality column is poorly named. It is a probability which is
    # correlated with incidence & prevalence
    colnam <- c("year", "age", "sex", "qimd", "disease",
                paste0(
                  type,
                  c(
                    "_total",
                    "_total_1",
                    "_total_99",
                    "_total_10",
                    "_total_20",
                    "_total_30",
                    "_total_40",
                    "_total_60",
                    "_total_70",
                    "_total_80",
                    "_total_90"
                  )
                ))
    disease_ <- disease
    indx <-
      read_fst("./inputs/exposure_distributions/mortality_projections_indx.fst",
               as.data.table = TRUE)[disease == disease_]

    lifetable_all <-
      read_fst(
        "./inputs/exposure_distributions/mortality_projections.fst",
        colnam,
        from = indx$from,
        to = indx$to,
        as.data.table = TRUE
      )[between(year,
                design$init_year,
                design$init_year + design$sim_horizon_max) &
          between(age, design$ageL, design$ageH)]
    if (design$stochastic) {
      absorb_dt(lifetable_all, prb)
      setnames(lifetable_all, colnam, gsub(paste0("^", type, "_"), "", colnam))

      lifetable_all[fatality < 0.1, qx_mc :=
                      qunif(normalise(c(fatality, 0, 0.1))[1],
                            total_1, total_10)]
      lifetable_all[between(fatality, 0.1, 0.2), qx_mc :=
                      qunif(normalise(c(fatality, 0.1, 0.2))[1],
                            total_10, total_20)]
      lifetable_all[between(fatality, 0.2, 0.3), qx_mc :=
                      qunif(normalise(c(fatality, 0.2, 0.3))[1],
                            total_20, total_30)]
      lifetable_all[between(fatality, 0.3, 0.4), qx_mc :=
                      qunif(normalise(c(fatality, 0.3, 0.4))[1],
                            total_30, total_40)]
      lifetable_all[between(fatality, 0.4, 0.5), qx_mc :=
                      qunif(normalise(c(fatality, 0.4, 0.5))[1],
                            total_40, total)]
      lifetable_all[between(fatality, 0.5, 0.6), qx_mc :=
                      qunif(normalise(c(fatality, 0.5, 0.6))[1],
                            total, total_60)]
      lifetable_all[between(fatality, 0.6, 0.7), qx_mc :=
                      qunif(normalise(c(fatality, 0.6, 0.7))[1],
                            total_60, total_70)]
      lifetable_all[between(fatality, 0.7, 0.8), qx_mc :=
                      qunif(normalise(c(fatality, 0.7, 0.8))[1],
                            total_70, total_80)]
      lifetable_all[between(fatality, 0.8, 0.9), qx_mc :=
                      qunif(normalise(c(fatality, 0.8, 0.9))[1],
                            total_80, total_90)]
      lifetable_all[between(fatality, 0.9, 1), qx_mc :=
                      qunif(normalise(c(fatality, 0.9, 1))[1],
                            total_90, total_99)]
      lifetable_all[fatality == 0.5, qx_mc := total]
      return(lifetable_all[, .(year, age, sex, qimd, qx_mc)])
    } else {
      setnames(lifetable_all, paste0(type, "_total"), "qx_mc")
      return(lifetable_all[, .(year, age, sex, qimd, qx_mc)])
    }
  }

# get_lifetable_all(4, "chd", design, "qx")[age == 60 & qimd == "3" & sex ==
# "men", plot(year, qx_mc)]


# Get population estimates adjusted for mortality
#' @export
generate_pop_adj_for_mrtl <-
  function(mc, dt, design, update_dt = FALSE) {
    tt <- get_lifetable_all(mc = mc, "allcause", design = design, "qx")
    orig_pops <-
      pops <-
      dt[, .(N = as.numeric(.N)), keyby = .(year, age, sex, qimd)]
    absorb_dt(pops, tt)

    out <- data.table()

    ttt <-
      pops[year == design$init_year &
             between(age, design$ageL, design$ageH)][, N_adj := as.numeric(N)]
    out <- rbind(out, ttt)
    for (i in seq_along(unique(pops$year))) {
      ttt1 <- # + 1 here correct, not + i
        copy(ttt)[, `:=`(year = year + 1L, age = age + 1L)]
      ttt <-
        pops[year == design$init_year + i &
               between(age, design$ageL, design$ageH)]
      ttt[ttt1, on = c("age", "sex", "qimd"), N := i.N_adj]
      ttt[, N_adj := N * (1 - qx_mc)]
      out <- rbind(out, ttt)
    }
    out[orig_pops, N := i.N, on = c("year", "age", "sex", "qimd")]
    if (update_dt) {
      dt[out, on = c("year", "age", "sex", "qimd"), pops_adj := N_adj]
      message("dt was updated.")
    }
    out[, qx_mc := NULL]
    out[]
  }
# generate_pop_adj_for_mrtl(mc, dt, design, TRUE)
#
# dt[is.na(pops_adj) & between(year, design$init_year, design$init_year +
# design$sim_horizon) & between(age, design$ageL, design$ageH), .N]


#' @export
generate_rns <- function(mc, dt, colnams) {
  dqRNGkind("pcg64") # dqRNGkind("Xoroshiro128+") ~10% faster
  SEED <- 4719349L # sample(1e7, 1)
  set.seed(SEED + mc)
  dqset.seed(SEED, mc)

  nrows <- nrow(dt)
  for (nam in colnams)
    set(dt, NULL, nam, dqrunif(nrows))
  invisible(dt)
}


# Given a correlation matrix (Pearson), produces a matrix of
# correlated uniforms

#' @export
generate_corr_unifs <- function(n, M) {
  # generate normals, check correlations
  # from http://comisef.wikidot.com/tutorial:correlateduniformvariates
  stopifnot(is.matrix(M))
  # Check that matrix is semi-positive definite
  # NOTE next line crashes frequently!! see
  # https://stat.ethz.ch/pipermail/r-help/2006-March/102703.html for a
  # workaround and https://stat.ethz.ch/pipermail/r-help/2006-March/102647.html
  # for some explanation
  # stopifnot(min(eigen(M, only.values = TRUE)$values) >= 0)

  M_original <- M


  # adjust correlations for uniforms
  for (i in seq_len(dim(M)[[1L]])) {
    for (j in seq_len(dim(M)[[2L]])) {
      if (i != j) {
        M[i, j] <- 2 * sin(pi * M[i, j] / 6)
        M[j, i] <- 2 * sin(pi * M[j, i] / 6)
      }
    }
  }

  X <- matrix(dqrnorm(n * dim(M)[[2]]), n)
  colnames(X) <- colnames(M)

  # induce correlation, check correlations
  Y <- pnorm(X %*% chol(M))

  # message(paste0("Mean square error is: ", signif(sum((cor(Y) - M_original) ^
  # 2), 3)))
  return(Y)
}


# Estimate health utility decreaments
#' @export
generate_eq5d_decr <- function(dt) {
  # From Sullivan et al. 2011
  # TODO add copd and cancers
  utility_pop_norms <-
    c(
      0.922,
      0.922,
      0.922,
      0.922,
      0.922,
      0.922,
      0.922,
      0.914,
      0.914,
      0.914,
      0.914,
      0.914,
      0.914,
      0.914,
      0.914,
      0.914,
      0.914,
      0.888,
      0.888,
      0.888,
      0.888,
      0.888,
      0.888,
      0.888,
      0.888,
      0.888,
      0.888,
      0.854,
      0.854,
      0.854,
      0.854,
      0.854,
      0.854,
      0.854,
      0.854,
      0.854,
      0.854,
      0.814,
      0.814,
      0.814,
      0.814,
      0.814,
      0.814,
      0.814,
      0.814,
      0.814,
      0.814,
      0.775,
      0.775,
      0.775,
      0.775,
      0.775,
      0.775,
      0.775,
      0.775,
      0.775,
      0.775,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706,
      0.706
    )
  utility_income    <-
    c(0.012572, 0.038783, 0.0396568, 0.0396568, 0.0408501)
  utility_education <-
    c(0.0060444,
      0.0056836,
      0.0028418,
      0.0028418,
      0.0028418,
      0.0056836,
      0)
  utility_ncc       <-
    c(0,
      0,
      -0.0528,
      -0.0415,
      -0.0203,
      0.0083,
      0.04087,
      0.06687,
      0.11589,
      0.13444,
      0.18361)

  new_ncc <-
    fclamp_int(
      dt$ncc + (dt$htn_prvl > 0L) + (dt$af_prvl > 0L) +
        (dt$t2dm_prvl > 0L) +  (dt$chd_prvl > 0L) +
        (dt$stroke_prvl > 0L) + (dt$poststroke_dementia_prvl > 0L) +
        (dt$copd_prvl > 0L) + (dt$lung_ca_prvl > 0L) +
        (dt$colon_ca_prvl > 0L) + (dt$breast_ca_prvl > 0L),
      0L,
      10L
    )


  out <-
    utility_pop_norms[dt$age - 17L] + utility_income[dt$income] +
    utility_education[dt$education] +
    0.0010046 * (dt$sex == "men") +
    utility_ncc[new_ncc + 1L] - 0.0460 * (dt$htn_prvl > 0L) -
    0.0384 * (dt$af_prvl > 0L) - 0.0714 * (dt$t2dm_prvl > 0L) -
    0.0679 * (dt$chd_prvl > 0L) - 0.0578 * (dt$stroke_prvl > 0L) -
    0.0957 * (dt$copd_prvl > 0L) - 0.1192427 * (dt$lung_ca_prvl > 0L) -
    0.0673908 * (dt$colon_ca_prvl > 0L) - 0.0194279 * (dt$breast_ca_prvl > 0L) -
    0.2165659 * (dt$poststroke_dementia_prvl > 0L)

  # half eq5d for year of death
  out <-
    out / ((dt$all_cause_mrtl > 0) + 1L) # out/1 for alive and out/2 for dead

  clamp(out, 0, 1, TRUE)
  dt[, eq5d := out]
}
# generate_eq5d_decr(output)
# output[, summary(eq5d)]

get_healthcare_costs <- function(mc) {
  if (mc < 1L | mc > 1000L) stop("mc need to be between 1 and 1000")
  out <- list()
  tt <- read_fst("./simulation/health_econ/healthcare_costs_indx.fst",
    from = mc, to = mc,
    as.data.table = TRUE)
  tt <- read_fst("./simulation/health_econ/healthcare_costs_l.fst",
    from = tt$from, to = tt$to,
    columns = c("disease", "years_since_diagnosis", "healthcare_cost"),
    as.data.table = TRUE)
  out$lung_ca_costs_young <- tt[disease == "lung_ca_18_64", c(0, healthcare_cost)] # (first element is for prevalence 0)
  out$lung_ca_costs_old <- tt[disease == "lung_ca_65+", c(0, healthcare_cost)]
  out$colon_ca_costs_young <- tt[disease == "colon_ca_18_64", c(0, healthcare_cost)]
  out$colon_ca_costs_old <- tt[disease == "colon_ca_65+", c(0, healthcare_cost)]
  out$breast_ca_costs_young <- tt[disease == "breast_ca_18_64", c(0, healthcare_cost)]
  out$breast_ca_costs_old <- tt[disease == "breast_ca_65+", c(0, healthcare_cost)]
  out$other <- tt[disease == "other", healthcare_cost]
  out$htn <- tt[disease == "hypertension", healthcare_cost]
  out$af <- tt[disease == "atrial_fibrillation", healthcare_cost]
  out$t2dm <- tt[disease == "t2dm", healthcare_cost]
  out$chd <- tt[disease == "chd", healthcare_cost]
  out$stroke_y1 <- tt[disease %in% c("stroke_acute_event", "stroke_first_year"),
    sum(healthcare_cost)]
  out$stroke_posty1 <- tt[disease == "stroke_year_2+", healthcare_cost]
  out$poststroke_dementia <- tt[disease == "dementia", healthcare_cost]
  out$copd <- tt[disease == "copd", healthcare_cost]
  out
}

generate_healthcare_costs <- function(dt, mc) {
  costs <- get_healthcare_costs(mc)

  out <- costs$other + costs$htn * (dt$htn_prvl > 0L) +
    costs$t2dm * (dt$t2dm_prvl > 0L) + costs$chd * (dt$chd_prvl > 0L) +
    costs$stroke_y1 * (dt$stroke_prvl == 1L) +
    costs$stroke_posty1 * (dt$stroke_prvl > 1L) +
    costs$copd * (dt$copd_prvl > 0L) + costs$af * (dt$af_prvl > 0L) +
    costs$lung_ca_costs_young[1L + dt$lung_ca_prvl] * (dt$age < 65) +
    costs$lung_ca_costs_old[1L + dt$lung_ca_prvl] * (dt$age >= 65) +
    costs$colon_ca_costs_young[1L + dt$colon_ca_prvl] * (dt$age < 65) +
    costs$colon_ca_costs_old[1L + dt$colon_ca_prvl] * (dt$age >= 65) +
    costs$breast_ca_costs_young[1L + dt$breast_ca_prvl] * (dt$age < 65) +
    costs$breast_ca_costs_old[1L + dt$breast_ca_prvl] * (dt$age >= 65) +
    costs$poststroke_dementia * (dt$poststroke_dementia_prvl > 0L)

  # half cost for year of death
  out <-
    out / ((dt$all_cause_mrtl > 0) + 1L) # out/1 for alive and out/2 for dead
  dt[, healthcare_cost := out]
}


get_socialcare_costs <- function(mc) {
  if (mc < 1L | mc > 1000L) stop("mc need to be between 1 and 1000")
  out <- list()
  tt <- read_fst("./simulation/health_econ/socialcare_costs_indx.fst",
    from = mc, to = mc,
    as.data.table = TRUE)
  out$socialcare_cost <- read_fst("./simulation/health_econ/socialcare_costs_l.fst",
    from = tt$from, to = tt$to,
    columns = c("socialcare_cost"),
    as.data.table = FALSE)$socialcare_cost
  tt <- read_fst("./simulation/health_econ/socialcare_costs_added_diseases_indx.fst",
    from = mc, to = mc,
    as.data.table = TRUE)
  tt <- read_fst("./simulation/health_econ/socialcare_costs_added_diseases_l.fst",
    from = tt$from, to = tt$to,
    columns = c("disease", "socialcare_cost"),
    as.data.table = TRUE)
  out$socialcare_cost_stroke <- tt[disease == "stroke", socialcare_cost]
  out$socialcare_cost_poststroke_dementia <- tt[disease == "poststroke_dementia", socialcare_cost]
  out
}

generate_socialcare_costs <- function(dt, mc) {
  costs <- get_socialcare_costs(mc)
  # for ages 18 to 100. Hence, age - 17L
  out <-
    costs$socialcare_cost[dt$age - 17L] + costs$socialcare_cost_stroke * (dt$stroke_prvl > 0L) +
    costs$socialcare_cost_poststroke_dementia * (dt$poststroke_dementia_prvl > 0L)

  # half cost for year of death
  out <-
    out / ((dt$all_cause_mrtl > 0) + 1L) # out/1 for alive and out/2 for dead
  dt[, socialcare_cost := out]
  dt
}

get_productivity_costs <- function(mc) {
  if (mc < 1L | mc > 1000L) stop("mc need to be between 1 and 1000")
  tt <- read_fst("./simulation/health_econ/productivity_costs_indx.fst",
    from = mc, to = mc,
    as.data.table = TRUE)
  read_fst("./simulation/health_econ/productivity_costs_l.fst",
    from = tt$from, to = tt$to,
    columns = c("age", "sex", "eq5d_r", "productivity_cost"),
    as.data.table = TRUE)
}

generate_productivity_costs <- function(dt, mc) {
  tt <- get_productivity_costs(mc)
  tt[, eq5d_r := as.integer(100L * eq5d_r)]
  dt[, eq5d_r := as.integer(100L * round(eq5d / 0.05) * 0.05)]
  absorb_dt(dt, tt)
  # half cost for year of death
  dt[all_cause_mrtl > 0, productivity_cost := productivity_cost / 2]
  dt[, eq5d_r := NULL]
  dt
}


get_informal_care_costs <- function(mc) {
  if (mc < 1L | mc > 1000L) stop("mc need to be between 1 and 1000")
  tt <- read_fst("./simulation/health_econ/informal_care_costs_indx.fst",
    from = mc, to = mc,
    as.data.table = TRUE)
  read_fst("./simulation/health_econ/informal_care_costs_l.fst",
    from = tt$from, to = tt$to,
    columns = c("age", "sex", "eq5d_r", "informal_care_cost"),
    as.data.table = TRUE)
}

generate_informal_care_costs <- function(dt, mc) {
  tt <- get_informal_care_costs(mc)
  tt[, eq5d_r := as.integer(100L * eq5d_r)]
  dt[, eq5d_r := as.integer(100 * round(eq5d / 0.05) * 0.05)]
  absorb_dt(dt, tt)
  # half cost for year of death
  dt[all_cause_mrtl > 0, informal_care_cost := informal_care_cost / 2]
  dt[, eq5d_r := NULL]
  dt
}

#' @export
generate_health_econ <- function(dt, mc) {
  generate_eq5d_decr(dt)
  generate_healthcare_costs(dt, mc)
  generate_socialcare_costs(dt, mc)
  generate_productivity_costs(dt, mc)
  generate_informal_care_costs(dt, mc)
  invisible(dt)
}

#' @export
set_eligible <- function(scenario_parms, dt, hlp, env = parent.frame()) {
  colnam <- "eligible_sc"
  if (scenario_parms$sc_ens_is && colnam %in% names(dt)) {
    env$hlp$previous_elig <- clamp(hlp$previous_elig + dt$eligible_sc)
  }

  set(dt, NULL, colnam, 0L)
  if (!scenario_parms$sc_eligib_noone) {
    if (!scenario_parms$sc_ens_parallel_is) {
      # if not parallel ensemble (following works OK with serial)
      dt[between(year,
                 scenario_parms$sc_init_year - 2000L,
                 scenario_parms$sc_last_year - 2000L) &
           between(age,
                   scenario_parms$sc_eligib_age[[1]],
                   scenario_parms$sc_eligib_age[[2]]) &
           dead == FALSE &
           # statin_px_curr_xps == 0L &
           # af_dgn == 0L &
           htn_dgn < fifelse(scenario_parms$sc_eligib_htn, Inf, 1L) &
           t2dm_dgn < fifelse(scenario_parms$sc_eligib_diab, Inf, 1L) &
           # ckd_prvl_curr_xps <= 3L &
           # ra_prvl == 0L &
           chd_dgn == 0L &
           stroke_dgn == 0L,
         (colnam) := 1L]
    } else { # if parallel ensemble
      dt[between(year,
                 scenario_parms$sc_init_year - 2000L,
                 scenario_parms$sc_last_year - 2000L) &
           between(age,
                   scenario_parms$sc_eligib_age[[1]],
                   scenario_parms$sc_eligib_age[[2]]) &
           dead == FALSE &
           # statin_px_curr_xps == 0L &
           # af_dgn == 0L &
           htn_dgn < fifelse(scenario_parms$sc_eligib_htn, Inf, 1) &
           t2dm_dgn < fifelse(scenario_parms$sc_eligib_diab, Inf, 1L) &
           # ckd_prvl_curr_xps <= 3L &
           # ra_prvl == 0L &
           chd_dgn == 0L &
           stroke_dgn == 0L &
           pid %in% hlp$sc_alloc[[scenario_parms$sc_name]],
         (colnam) := 1L]
    }
  }
  invisible(dt)
  # dt[year == 20 &
  #      between(age, scenario_parms$sc_eligib_age[[1]],
  # scenario_parms$sc_eligib_age[[2]]), prop_if(eligible_sc1 == 1)]
}
# set_eligible("sc2", POP, parameters_dt) POP[between(year, 18, 35) &
# between(age, 40, 74) & dead == FALSE, prop_if(eligible_sc1 == 1)] Around 70%
# of the population should be eligible. From
# https://www.healthcheck.nhs.uk/commissioners-and-providers/data/total-eligible-population/
# POP[between(year, 18, 35) & between(age, 40, 74) & dead == FALSE, 1 -
# prop_if(eligible_sc1 == 1), keyby = age]

#' @export
set_invitees <- function(scenario_parms, dt, hlp, env = parent.frame()) {
  colnam <- "invitees_sc"
  colnam_cost <- "invitation_cost_sc"
  elig_colnam <- "eligible_sc"
  if (scenario_parms$sc_ens_is && colnam %in% names(dt)) {
      env$hlp$previous_invitees <- clamp(hlp$previous_invitees + dt$invitees_sc)
  }
  set(dt, NULL, colnam, 0L)

  # TODO find better approach

  if (scenario_parms$sc_invit_detailed) {
    tt1 <- data.table(
      year = (scenario_parms$sc_init_year:scenario_parms$sc_last_year) - 2000L,
      qimd = "1 most deprived",
      mu = scenario_parms$sc_invit_qimd1,
      freq = scenario_parms$sc_eligib_freq
    )
    # for serial ensembles take into account previous scenarios
    if (scenario_parms$sc_ens_serial_is) {
      tt1 <- rbind(hlp$invit_tbl1, tt1)
      env$hlp$invit_tbl1 <- copy(tt1)
    }
    elig <- vector("numeric", nrow(tt1))
    for (i in seq_len(nrow(tt1))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt1[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt1[i, mu])
      if (i > tt1[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt1[i - 1, mu] + tt1[i - tt1[i, freq], mu])
    }
    tt1[, mu := clamp(mu / elig)]
    tt1[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt2 <-
      data.table(
        year = (scenario_parms$sc_init_year:scenario_parms$sc_last_year) - 2000L,
        qimd = "2",
        mu = scenario_parms$sc_invit_qimd2,
        freq = scenario_parms$sc_eligib_freq
      )
    if (scenario_parms$sc_ens_serial_is) {
      tt2 <- rbind(hlp$invit_tbl2, tt2)
      env$hlp$invit_tbl2 <- copy(tt2)
    }
    elig <- vector("numeric", nrow(tt2))
    for (i in seq_len(nrow(tt2))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt2[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt2[i, mu])
      if (i > tt2[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt2[i - 1, mu] + tt2[i - tt2[i, freq], mu])
    }
    tt2[, mu := clamp(mu / elig)]
    tt2[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt3 <-
      data.table(
        year = (scenario_parms$sc_init_year:scenario_parms$sc_last_year) - 2000L,
        qimd = "3",
        mu = scenario_parms$sc_invit_qimd3,
        freq = scenario_parms$sc_eligib_freq
      )
    if (scenario_parms$sc_ens_serial_is) {
      tt3 <- rbind(hlp$invit_tbl3, tt3)
      env$hlp$invit_tbl3 <- copy(tt3)
    }
    elig <- vector("numeric", nrow(tt3))
    for (i in seq_len(nrow(tt3))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt3[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt3[i, mu])
      if (i > tt3[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt3[i - 1, mu] + tt3[i - tt3[i, freq], mu])
    }
    tt3[, mu := clamp(mu / elig)]
    tt3[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt4 <-
      data.table(
        year = (scenario_parms$sc_init_year:scenario_parms$sc_last_year) - 2000L,
        qimd = "4",
        mu = scenario_parms$sc_invit_qimd4,
        freq = scenario_parms$sc_eligib_freq
      )
    if (scenario_parms$sc_ens_serial_is) {
      tt4 <- rbind(hlp$invit_tbl4, tt4)
      env$hlp$invit_tbl4 <- copy(tt4)
    }
    elig <- vector("numeric", nrow(tt4))
    for (i in seq_len(nrow(tt4))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt4[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt4[i, mu])
      if (i > tt4[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt4[i - 1, mu] + tt4[i - tt4[i, freq], mu])
    }
    tt4[, mu := clamp(mu / elig)]
    tt4[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt5 <-
      data.table(
        year = (scenario_parms$sc_init_year:scenario_parms$sc_last_year) - 2000L,
        qimd = "5 least deprived",
        mu = scenario_parms$sc_invit_qimd5,
        freq = scenario_parms$sc_eligib_freq
      )
    if (scenario_parms$sc_ens_serial_is) {
      tt5 <- rbind(hlp$invit_tbl5, tt5)
      env$hlp$invit_tbl5 <- copy(tt5)
    }
    elig <- vector("numeric", nrow(tt5))
    for (i in seq_len(nrow(tt5))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt5[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt5[i, mu])
      if (i > tt5[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt5[i - 1, mu] + tt5[i - tt5[i, freq], mu])
    }
    tt5[, mu := clamp(mu / elig)]
    tt5[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25

    tt <- rbind(tt1, tt2, tt3, tt4, tt5)
    tt <- tt[between(year, scenario_parms$sc_init_year - 2000L, scenario_parms$sc_last_year - 2000L)]

    ttcost <-
      data.table(
        V1 = 1L,
        V2 = c("1 most deprived", "2", "3", "4", "5 least deprived"),
        V3 = c(
          scenario_parms$sc_invit_qimd1_cost,
          scenario_parms$sc_invit_qimd2_cost,
          scenario_parms$sc_invit_qimd3_cost,
          scenario_parms$sc_invit_qimd4_cost,
          scenario_parms$sc_invit_qimd5_cost
        )
      )
    setnames(ttcost, c(colnam, "qimd", colnam_cost))

  } else { # if not detailed invites (input not by qimd)
    # The method below is robust in changes of mu and/or freq to define
    # probability of being invited given that eligible population is reduced the
    # more you invite for a check
    tt <- data.table(
      year = (scenario_parms$sc_init_year:scenario_parms$sc_last_year) - 2000L,
      mu = scenario_parms$sc_invit_qimdall,
      freq = scenario_parms$sc_eligib_freq
    )
    # for serial ensembles take into account previous scenarios
    if (scenario_parms$sc_ens_serial_is) {
      tt <- rbind(hlp$invit_tbl, tt)
      env$hlp$invit_tbl <- copy(tt)
    }

    elig <- vector("numeric", nrow(tt))
    for (i in seq_len(nrow(tt))) {
      if (i == 1L)
        elig[i] <- 1 - 0
      if (between(i, 2, tt[i, freq]))
        elig[i] <- clamp(elig[i - 1L] - tt[i, mu])
      if (i > tt[i, freq])
        elig[i] <-
          clamp(elig[i - 1L] - tt[i - 1, mu] + tt[i - tt[i, freq], mu])
    }
    tt[, mu := clamp(mu / elig)]
    tt[is.na(mu), mu := 0] # i.e. for freq 5 and prb 0.25
    tt <- tt[between(year, scenario_parms$sc_init_year - 2000L, scenario_parms$sc_last_year - 2000L)]

    ttcost <- data.table(1L, scenario_parms$sc_invit_qimdall_cost)
    setnames(ttcost, c(colnam, colnam_cost))
  }

  absorb_dt(dt, tt)
  setnafill(dt, "c", 0, cols = c("mu", "freq"))
  dt[, (colnam) :=
      identify_invitees(eligible_sc, hlp$previous_invitees, mu, freq, pid_mrk)]
  dt[, c("mu", "freq") := NULL]
  absorb_dt(dt, ttcost, on = setdiff(names(ttcost), colnam_cost))
  setnafill(dt, "c", 0, cols = colnam_cost)
  invisible(dt)
  # dt[ eligible_sc1 == 1L &
  #      between(age, scenario_parms$sc_eligib_age[[1]], scenario_parms$sc_eligib_age[[2]]),
  #      prop_if(invitees_sc1== 1), keyby = year]

}
# set_invitees("sc2", POP, parameters_dt)

#' @export
set_attendees <- function(scenario_parms, dt, scenario_nam, parameters_dt,
  design, hlp, env = parent.frame()) {
  colnam       <- "attendees_sc"
  colnam_cost  <- "attendees_cost_sc"
  invit_colnam <- "invitees_sc"
  if (scenario_parms$sc_ens_is && colnam %in% names(dt)) {
    env$hlp$previous_attendees <- clamp(hlp$previous_attendees + dt$attendees_sc)
  }
  set(dt, NULL, colnam_cost, 0)


  if (scenario_parms$sc_uptake_detailed) {
    agegrp <- fromGUI_uptake_table_agegrps(scenario_nam = scenario_nam,
                                           parameters_dt = parameters_dt)
    absorb_dt(dt, agegrp)
    dt[, "Qrisk2_cat" := Qrisk2(.SD,
                                scenario_parms$sc_qrisk_ignore_bmi,
                                scenario_parms$sc_qrisk_ignore_sbp,
                                scenario_parms$sc_qrisk_ignore_tchol)$Qrisk2_cat]
    if (scenario_parms$sc_uptake_structural0s) { # if structural 0s
      absorb_dt(dt, scenario_parms$sc_uptake)
      setnafill(dt, "c", 0, cols = "uptake_wt")
      dt[, c("Qrisk2_cat", "agegrp10") := NULL]
      tt <- sort(dt[invitees_sc == 1L & uptake_wt > 0,
                    sample_int_expj(.N, as.integer(round(scenario_parms$sc_uptake_all * .N)),
                                    uptake_wt)])
      tt <-
        dt[invitees_sc == 1L &
             uptake_wt > 0, .(year, pid)][tt, ][, (colnam) := 1L]
      absorb_dt(dt, tt, on = c("pid", "year"))
    } else { # if no structural 0s
      # Abuse of the rule of 3
      absorb_dt(dt, scenario_parms$sc_uptake[uptake_wt == 0,
                                uptake_wt := 0.5 * 3 /
                                  dt[year == scenario_parms$sc_init_year - 2000L &
                                       invitees_sc == 1L,
                                     sum(wt) *
                                       design$sim_prm$n_synthpop_aggregation]])
      setnafill(dt, "c", 0, cols = "uptake_wt")
      dt[, c("Qrisk2_cat", "agegrp10") := NULL]

      tt <- sort(dt[invitees_sc == 1L,
                    sample_int_expj(.N, as.integer(round(scenario_parms$sc_uptake_all * .N)),
                                    uptake_wt)])
      tt <-
        dt[invitees_sc == 1L, .(year, pid)][tt, ][, (colnam) := 1L]
      absorb_dt(dt, tt, on = c("pid", "year"))
    }
  } else {
    set(dt, NULL, "uptake_wt", scenario_parms$sc_uptake_all)
    dt[invitees_sc == 1L, (colnam) := rbinom(.N, 1, uptake_wt)]
  }

  setnafill(dt, "c", 0, cols = colnam)
  dt[attendees_sc == 1L, (colnam_cost) := scenario_parms$sc_uptake_all_cost]
  dt[, uptake_wt := NULL]
  invisible(dt)

  # dt[ invitees_sc1 == 1L & between(age, scenario_parms$sc_eligib_age[[1]],
  # scenario_parms$sc_eligib_age[[2]]), prop_if(attendees_sc1== 1), keyby = year]
}

#' @export
set_px <- function(scenario_parms, dt, mc, design_) {
  dt[, "Qrisk2_cat" := Qrisk2(.SD, FALSE, FALSE, FALSE)$Qrisk2_cat]
  atte_colnam <- "attendees_sc"

  # for statins
  colnam     <- "statin_px_sc"
  colnam_bio <- "tchol_sc"

  if (scenario_parms$sc_px_detailed) {
    absorb_dt(dt, scenario_parms$sc_px_statins_wt)

    # Below assumes people on statin_px_curr_xps but undertreated will titrate
    # statin treatment.

    # Adjusted prb where the denominator changed from all attendees to those
    # eligible for statins
    tt <- dt[attendees_sc == 1L,
             clamp(scenario_parms$sc_px_statins / prop_if(tchol_curr_xps >= 5))]


    tt <- sort(dt[attendees_sc == 1L & tchol_curr_xps >= 5,
                  sample_int_expj(.N, as.integer(round(tt * .N)),
                                  px_statins_wt)])
    # rows that will have statin px
    tt <-
      dt[attendees_sc == 1L & tchol_curr_xps >= 5,
         .(year, pid)][tt, ][, (colnam) := 1L]
    absorb_dt(dt, tt, on = c("pid", "year"))
  } else {
    tt <- dt[attendees_sc == 1L,
             clamp(scenario_parms$sc_px_statins / prop_if(tchol_curr_xps >= 5 &
                                               Qrisk2_cat != "low"))]

    set(dt, NULL, "px_statins_wt", tt)
    dt[attendees_sc == 1L &
         tchol_curr_xps >= 5 & Qrisk2_cat != "low",
       (colnam) := rbinom(.N, 1L, px_statins_wt)]
  }
  setnafill(dt, "c", 0, cols = colnam)
  dt[, (colnam) := hc_effect(statin_px_sc, 0.9749866, pid_mrk)]

  # 0.9749866 comes from the following study in Wales.
  # King W, Lacey A, White J, Farewell D, Dunstan F, Fone D. Socioeconomic
  # inequality in medication persistence in primary and secondary prevention of
  # coronary heart disease – A population-wide electronic cohort study. PLOS ONE
  # 2018;13:e0194081.
  # from 33228 individuals px a statin in Wales for primary prevention, 5378 had
  # discontinued it within 7 years without socioeconomic gradient statin <-
  # c(0.96, (33228 - 5378)/33228) year <- c(0.1, 7) m1 <- glm(statin~ -1 + year,
  # family = gaussian(link = "log")) exp(m1$coefficients) x <- predict(m1,
  # newdata = data.table(year = 0:7), type = "re") shift(x, -1)/x plot(0:70,
  # predict(m1, newdata = data.table(year = 0:70), type = "re"), ylim = c(0, 1))
  # Every year 0.9749866 of those taking statin continue next year


  # estimate tchol change atorvastatin effect from Law MR, et al. Quantifying
  # effect of statins on low density lipoprotein cholesterol, ischaemic heart
  # disease, and stroke: systematic review and meta-analysis. BMJ 2003;326:1423.
  # table 2. 43% (0.3958 - 0.46875) reduction of ldl. to convert to tc, tc/ldl =
  # 0.27/0.36 from Edwards JE, et al. Statins in hypercholesterolaemia: A
  # dose-specific meta-analysis of lipid changes in randomised, double blind
  # trials. BMC Family Practice 2003;4:18.

  atorv_eff <- RR$statins_tchol$get_rr(mc, design_, drop = TRUE)

  # adherence <- rpert(1e6, 0.5, 0.9, 1, 8)
  # proportion of prescribed dose taken
  # or to avoid dependency for rpert
  # adherence <- rBE(1e6, 0.9, 0.2) # proportion of prescribed dose taken

  if (!colnam_bio %in% names(dt)) dt[, (colnam_bio) := tchol_curr_xps]
  dt[statin_px_sc == 1L &
       # statin_px_curr_xps == 0L & # Assume GP titrates treatment
       year >= scenario_parms$sc_init_year - 2000L,
     (colnam_bio) := tchol_curr_xps * (1 - atorv_eff * statin_adherence)]

  # for bpmed
  colnam     <- "bpmed_px_sc"
  colnam_bio <- "sbp_sc"

  if (scenario_parms$sc_px_detailed) {
    absorb_dt(dt, scenario_parms$sc_px_antihtn_wt)

    # Adjusted prb where the denominator changed from all attendees to those
    # eligible for bpmed
    tt <- dt[attendees_sc == 1L,
             clamp(scenario_parms$sc_px_antihtn / prop_if(sbp_curr_xps >= 135))]

    tt <- sort(dt[attendees_sc == 1L & sbp_curr_xps >= 135,
                  sample_int_expj(.N, as.integer(round(tt * .N)),
                                  px_antihtn_wt)])
    tt <-
      dt[attendees_sc == 1L & sbp_curr_xps >= 135,
         .(year, pid)][tt, ][, (colnam) := 1L]
    absorb_dt(dt, tt, on = c("pid", "year"))
  } else {
    tt <- dt[attendees_sc == 1L,
             clamp(scenario_parms$sc_px_antihtn / prop_if(sbp_curr_xps >= 135))]

    set(dt, NULL, "px_antihtn_wt", tt)
    dt[attendees_sc == 1L & sbp_curr_xps >= 135,
       (colnam) := rbinom(.N, 1L, px_antihtn_wt)]
  }
  setnafill(dt, "c", 0, cols = colnam)
  dt[, (colnam) := hc_effect(bpmed_px_sc, 0.9749866, pid_mrk)]
  # assume same prb as statins

  # Estimate sbp change
  if (!colnam_bio %in% names(dt)) dt[, (colnam_bio) := sbp_curr_xps]
  dt[bpmed_px_sc == 1L &
       # bpmed_curr_xps == 0L & # Assume GP titrates treatment
       year >= scenario_parms$sc_init_year - 2000L,
     (colnam_bio) := sbp_curr_xps - clamp((sbp_curr_xps - 135) *
         bpmed_adherence, 0, 1e3)]
  # Assume that antihtn medication can potentially achieve sbp 135 for all. Not
  # 110 to account for residual risk

  dt[, c("Qrisk2_cat", "px_statins_wt", "px_antihtn_wt") := NULL]

  invisible(dt)
}

#' @export
set_lifestyle <-
  function(scenario_parms, dt, design) {
    atte_colnam <- "attendees_sc"

    # PA ----
    colnam      <- "active_days_sc"
    colnam_cost <- "active_days_cost_sc"
    if (!"hc_eff_pa" %in% names(dt))      set(dt, NULL, "hc_eff_pa", 0L)
    if (!colnam %in% names(dt))        dt[, (colnam) := active_days_curr_xps]
    if (!colnam_cost %in% names(dt))   set(dt, NULL, colnam_cost, 0)
    dt[attendees_sc == 1L, hc_eff_pa := rbinom(.N, 1, scenario_parms$sc_ls_papct)]
    dt[hc_eff_pa == 1L & year >= scenario_parms$sc_init_year - 2000L,
       (colnam_cost) := scenario_parms$sc_ls_pa_cost_ind]
    # Cost only the year of referral
    dt[,
      hc_eff_pa := hc_effect(hc_eff_pa, (1 - scenario_parms$sc_ls_attrition), pid_mrk)]
    dt[hc_eff_pa == 1L, (colnam) :=
         as.integer(round(clamp(active_days_sc + scenario_parms$sc_ls_papincr, 0, 7)))]

    # Weight management ----
    colnam      <- "bmi_sc"
    colnam_cost <- "bmi_cost_sc"
    if (!"hc_eff_wm" %in% names(dt)) set(dt, NULL, "hc_eff_wm", 0L)
    if (!colnam %in% names(dt))      dt[, (colnam) := bmi_curr_xps]
    if (!colnam_cost %in% names(dt)) set(dt, NULL, colnam_cost, 0)
    dt[attendees_sc == 1L &
         bmi_curr_xps > 30, hc_eff_wm := rbinom(.N, 1, scenario_parms$sc_ls_wghtpct)]
    dt[hc_eff_wm == 1L  & year >= scenario_parms$sc_init_year - 2000L,
       (colnam_cost) := scenario_parms$sc_ls_wghtloss_cost_ind]
    # Cost only the year of referral
    dt[,
      hc_eff_wm := hc_effect(hc_eff_wm,
        (1 - scenario_parms$sc_ls_attrition), pid_mrk)]
    dt[hc_eff_wm == 1L, (colnam) := bmi_sc * (1 - scenario_parms$sc_ls_wghtreduc)]

    # Alcohol ----
    colnam      <- "alcohol_sc"
    colnam_cost <- "alcohol_cost_sc"
    if (!"hc_eff_al" %in% names(dt)) set(dt, NULL, "hc_eff_al", 0L)
    if (!colnam %in% names(dt))      dt[, (colnam) := alcohol_curr_xps]
    if (!colnam_cost %in% names(dt)) set(dt, NULL, colnam_cost, 0)
    dt[attendees_sc == 1L &
         alcohol_curr_xps >= 16, hc_eff_al := rbinom(.N, 1, scenario_parms$sc_ls_alcoholpct)]
    dt[hc_eff_al == 1L & year >= scenario_parms$sc_init_year - 2000L,
       (colnam_cost) := scenario_parms$sc_ls_alcoholreduc_cost_ind]
    # Cost only the year of referral
    dt[,
       hc_eff_al := hc_effect(hc_eff_al, (1 - scenario_parms$sc_ls_attrition), pid_mrk)]
    dt[hc_eff_al == 1L, (colnam) :=
         as.integer(round(alcohol_sc * (1 - scenario_parms$sc_ls_alcoholreduc)))]


    # Smoking cessation ----
    colnam_status   <- "smok_status_sc"
    colnam_quit_yrs <- "smok_quit_yrs_sc"
    colnam_dur      <- "smok_dur_sc"
    colnam_cig      <- "smok_cig_sc"
    colnam_cost     <- "smoking_cost_sc"
    if (!"hc_eff_sm" %in% names(dt))     set(dt, NULL, "hc_eff_sm", 0L)
    if (!colnam_status %in% names(dt))   dt[, (colnam_status) := smok_status_curr_xps]
    if (!colnam_quit_yrs %in% names(dt)) dt[, (colnam_quit_yrs) := smok_quit_yrs_curr_xps]
    if (!colnam_dur %in% names(dt))      dt[, (colnam_dur) := smok_dur_curr_xps]
    if (!colnam_cig %in% names(dt))      dt[, (colnam_cig) := smok_cig_curr_xps]
    if (!colnam_cost %in% names(dt)) set(dt, NULL, colnam_cost, 0)
    dt[attendees_sc == 1L & smok_status_curr_xps == "4",
       hc_eff_sm := rbinom(.N, 1, scenario_parms$sc_ls_smkcess)]
    dt[hc_eff_sm == 1L  & year >= scenario_parms$sc_init_year - 2000L,
       (colnam_cost) := scenario_parms$sc_ls_smkcess_cost_ind]
    # Cost only the year of referral

    # Handle smok_relapse probabilities
    tbl <-
      read_fst("./inputs/exposure_distributions/smok_relapse_table.fst",
               as.data.table = TRUE)
    tbl <- dcast(tbl, sex + qimd ~ smok_quit_yrs, value.var = "pr")
    nam <- tbl[, paste0(sex, " ", qimd)]
    tbl <- as.matrix(tbl[, mget(paste0(1:15))], rownames = nam)

    dt[,
       (c(colnam_status, colnam_quit_yrs, colnam_dur)) :=
         simsmok_cessation(
           smok_status_sc,
           smok_quit_yrs_sc,
           smok_dur_sc,
           sex,
           qimd,
           pid_mrk,
           hc_eff_sm,
           dqrunif(.N),
           tbl,
           design$sim_prm$smoking_relapse_limit
         )]

    dt[, smok_status_sc := factor(smok_status_sc)]
    # needed for QRisk and QDrisk
    dt[, smoke_cat_sc := 0L]
    dt[smok_status_sc == "3", smoke_cat_sc := 1L]
    dt[smok_status_sc == "4", smoke_cat_sc := 3L]
    dt[smok_status_sc == "4" & smok_cig_sc < 10L, smoke_cat_sc := 2L]
    dt[smok_status_sc == "4" & smok_cig_sc > 19L, smoke_cat_sc := 4L]

    invisible(dt)
  }

#' @export
set_structural <-
  function(scenario_parms, dt, design) {
    if (any(scenario_parms[grepl("^sc_str_", names(scenario_parms))] != 0)) {
      row_sel <-
        dt[between(year,
          scenario_parms$sc_init_year - 2000L,
          scenario_parms$sc_last_year - 2000L) &
            dead == FALSE, which = TRUE]
    }

    # smoking ----
    if (scenario_parms$sc_str_smk_change != 0) {
      if (!"smok_status_sc" %in% names(dt))
        set(dt, NULL, "smok_status_sc", dt$smok_status_curr_xps)
      if (!"smok_quit_yrs_sc" %in% names(dt))
        set(dt, NULL, "smok_quit_yrs_sc", dt$smok_quit_yrs_curr_xps)
      if (!"smok_dur_sc" %in% names(dt))
        set(dt, NULL, "smok_dur_sc", dt$smok_dur_curr_xps)
      if (!"smok_cig_sc" %in% names(dt))
        set(dt, NULL, "smok_cig_sc", dt$smok_cig_curr_xps)


      if (scenario_parms$sc_str_smk_change < 0) {

        dt[between(year, scenario_parms$sc_init_year - 2000L,
                   scenario_parms$sc_last_year - 2000L) &
             dead == FALSE & smok_status_sc == "4",
           hc_eff := rbinom(.N, 1, -scenario_parms$sc_str_smk_change)]

        dt[, (c("smok_status_sc", "smok_quit_yrs_sc", "smok_dur_sc", "smok_cig_sc")) :=
             simsmok_policy_impact_decr(
               smok_status_sc,
               smok_quit_yrs_sc,
               smok_dur_sc,
               smok_cig_sc,
               pid_mrk,
               hc_eff
               )]
      }

      if (scenario_parms$sc_str_smk_change > 0) {
        # calculate policy effect with those quit smoking recently be more
        # likely to relapse
        tt <- dt[between(year, scenario_parms$sc_init_year - 2000L,
                   scenario_parms$sc_last_year - 2000L) &
             dead == FALSE, .("ex"   = sum(smok_status_sc == "3"),
                              "curr" = sum(smok_status_sc == "4")), keyby = year]
        tt[, impacted := round(curr * scenario_parms$sc_str_smk_change)]

        # Make change to add up every year (for Vincy's SCC abstract)
        # tt[, impacted := round(curr * scenario_parms$sc_str_smk_change *
        #     (year - min(year) + 1L))]

        dt[tt, `:=`(impacted = i.impacted,
                    ex = i.ex), on = "year"]
        dt[between(year, scenario_parms$sc_init_year - 2000L,
                   scenario_parms$sc_last_year - 2000L) &
             dead == FALSE & smok_status_sc == "3",
           rid := 1:.N, by = year]
        dt[, hc_eff := 0L]
        tt <- dt[between(year,
                         scenario_parms$sc_init_year - 2000L,
                         scenario_parms$sc_last_year - 2000L) &
                   dead == FALSE & smok_status_sc == "3",
                 .(rid = sample_int_expj(first(ex), first(impacted),
                                         (smok_quit_yrs_sc + 1) ^
                                           -1)),
                 keyby = year]
        dt[tt, hc_eff := 1L, on = .(year, rid)]
        dt[, c("impacted", "ex", "rid") := NULL]

        dt[, (c("smok_status_sc", "smok_quit_yrs_sc", "smok_dur_sc")) :=
             simsmok_policy_impact_incr(
               smok_status_sc,
               smok_quit_yrs_sc,
               smok_dur_sc,
               pid_mrk,
               hc_eff
             )]
      }

      dt[, smok_status_sc := factor(smok_status_sc)]
      # needed for QRisk and QDrisk
      dt[, smoke_cat_sc := 0L]
      dt[smok_status_sc == "3", smoke_cat_sc := 1L]
      dt[smok_status_sc == "4", smoke_cat_sc := 3L]
      dt[smok_status_sc == "4" & smok_cig_sc < 10L, smoke_cat_sc := 2L]
      dt[smok_status_sc == "4" & smok_cig_sc > 19L, smoke_cat_sc := 4L]

      dt[, hc_eff := NULL]
}

    # fv ----
    if (scenario_parms$sc_str_fv_change != 0) {
      if (!"fruit_sc" %in% names(dt))
        set(dt, NULL, "fruit_sc", dt$fruit_curr_xps)
      if (!"veg_sc" %in% names(dt))
        set(dt, NULL, "veg_sc", dt$veg_curr_xps)
      dt[row_sel,
         `:=`(
           fruit_sc = as.integer(round(fruit_sc * (
             1 + scenario_parms$sc_str_fv_change
           ))),
           veg_sc = as.integer(round(veg_sc * (
             1 + scenario_parms$sc_str_fv_change
           ))))]
    }

    # alcohol ----
    if (scenario_parms$sc_str_alcohol_change != 0) {
      if (!"alcohol_sc" %in% names(dt))
        set(dt, NULL, "alcohol_sc", dt$alcohol_curr_xps)
      dt[row_sel,
         `:=`(
           alcohol_sc = as.integer(round(alcohol_sc * (
             1 + scenario_parms$sc_str_alcohol_change
           )))
           )]
    }

    # active_days ----
    if (scenario_parms$sc_str_pa_change != 0) {
      if (!"active_days_sc" %in% names(dt))
        set(dt, NULL, "active_days_sc", dt$active_days_curr_xps)
      dt[row_sel,
         `:=`(active_days_sc = active_days_sc + scenario_parms$sc_str_pa_change
         )]
    }

    # bmi ----
    if (scenario_parms$sc_str_bmi_change != 0) {
      if (!"bmi_sc" %in% names(dt))
        set(dt, NULL, "bmi_sc", dt$bmi_curr_xps)
      dt[row_sel,
         `:=`(
           bmi_sc = bmi_sc * (1 + scenario_parms$sc_str_bmi_change)
         )]
    }

    # sbp ----
    if (scenario_parms$sc_str_sbp_change != 0) {
      if (!"sbp_sc" %in% names(dt))
        set(dt, NULL, "sbp_sc", dt$sbp_curr_xps)
      dt[row_sel,
         `:=`(
           sbp_sc = sbp_sc * (1 + scenario_parms$sc_str_sbp_change)
         )]
    }

    # tchol ----
    if (scenario_parms$sc_str_tchol_change != 0) {
      if (!"tchol_sc" %in% names(dt))
        set(dt, NULL, "tchol_sc", dt$tchol_curr_xps)
      dt[row_sel,
         `:=`(
           tchol_sc = tchol_sc * (1 + scenario_parms$sc_str_tchol_change)
         )]
    }

    invisible(dt)
  }

#' @export
set_social <- function(scenario_parms, dt, design) {
  # bypass if irrelevant
  if (all(
    scenario_parms$sc_soc_qimd1_change == 1L,
    scenario_parms$sc_soc_qimd2_change == 2L,
    scenario_parms$sc_soc_qimd3_change == 3L,
    scenario_parms$sc_soc_qimd4_change == 4L,
    scenario_parms$sc_soc_qimd5_change == 5L
  )) {
    return(invisible(dt))
  } else {
    # if any relevant scenario input

    # Manipulate qimd per user input ----
    set(dt, NULL, "qimd_sc", dt$qimd) # create new scenario qimd
    l <- levels(dt$qimd)
    if (scenario_parms$sc_soc_qimd1_change != 1L)
      dt[qimd == l[1], qimd_sc := l[scenario_parms$sc_soc_qimd1_change]]
    if (scenario_parms$sc_soc_qimd2_change != 2L)
      dt[qimd == l[2], qimd_sc := l[scenario_parms$sc_soc_qimd2_change]]
    if (scenario_parms$sc_soc_qimd3_change != 3L)
      dt[qimd == l[3], qimd_sc := l[scenario_parms$sc_soc_qimd3_change]]
    if (scenario_parms$sc_soc_qimd4_change != 4L)
      dt[qimd == l[4], qimd_sc := l[scenario_parms$sc_soc_qimd4_change]]
    if (scenario_parms$sc_soc_qimd5_change != 5L)
      dt[qimd == l[5], qimd_sc := l[scenario_parms$sc_soc_qimd5_change]]



    row_sel <- # Indices of eligible rows
      dt[between(year,
        scenario_parms$sc_init_year - 2000L,
        scenario_parms$sc_last_year - 2000L) &
          # dead == FALSE & # NOTE this is appropriate
          qimd != qimd_sc, which = TRUE]


    # smoking ----
    # Assumes that from the smoking initiation/cessation/relapse probabilities
    # change, not smoking prevalence. Smoking intensity also changes.
    # pid_mrk needs to be recalculated for row_sel
    if ("smok" %in% scenario_parms$sc_soc_qimd_rf_change) {
      if (!"smok_status_sc" %in% names(dt))
        set(dt, NULL, "smok_status_sc", dt$smok_status_curr_xps)
      if (!"smok_quit_yrs_sc" %in% names(dt))
        set(dt, NULL, "smok_quit_yrs_sc", dt$smok_quit_yrs_curr_xps)
      if (!"smok_dur_sc" %in% names(dt))
        set(dt, NULL, "smok_dur_sc", dt$smok_dur_curr_xps)
      if (!"smok_cig_sc" %in% names(dt))
        set(dt, NULL, "smok_cig_sc", dt$smok_cig_curr_xps)

      dt[row_sel, pid_mrk_sc := mk_new_simulant_markers(pid)]

      # Assign smok_incid probabilities
      lutbl <-
        read_fst("./inputs/exposure_distributions/smok_incid_table.fst",
          as.data.table = TRUE)
      setnames(lutbl, c("qimd", "mu"), c("qimd_sc", "prb_smok_incid_sc"))
      lookup_dt(dt, lutbl)


      # Assign smok_cessation probabilities
      lutbl <-
        read_fst("./inputs/exposure_distributions/smok_cess_table.fst",
          as.data.table = TRUE)
      setnames(lutbl, c("qimd", "mu"), c("qimd_sc", "prb_smok_cess_sc"))
      lookup_dt(dt, lutbl)

      # Handle smok_relapse probabilities
      # No need to use qimd_sc here. It happens at the simsmok_sc side
      tbl <-
        read_fst("./inputs/exposure_distributions/smok_relapse_table.fst",
          as.data.table = TRUE)
      tbl <-
        dcast(tbl, sex + qimd ~ smok_quit_yrs, value.var = "pr")
      nam <- tbl[, paste0(sex, " ", qimd)]
      tbl <-
        as.matrix(tbl[, mget(paste0(1:15))], rownames = nam)

      simsmok_sc(dt, tbl, design$sim_prm$smoking_relapse_limit, row_sel)

      dt[, c("prb_smok_incid_sc", "prb_smok_cess_sc") := NULL]

      # smok intensity
      lutbl <-
        read_fst("./inputs/exposure_distributions/smok_cig_curr_table.fst",
          as.data.table = TRUE)
      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))

      dt[row_sel, mrk := TRUE]

      dt[(mrk) & smok_status_sc == "4",
        smok_cig_sc := qZINBI(rankstat_smok_cig_curr, mu, sigma, nu)]

      lutbl <-
        read_fst("./inputs/exposure_distributions/smok_cig_ex_table.fst",
          as.data.table = TRUE)
      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[(pid_mrk_sc) & # no need for mrk as superseded by pid_mrk_sc
          smok_status_sc == "3",
        smok_cig_sc := my_qZABNB(rankstat_smok_cig_ex,
          mu,
          sigma,
          nu,
          tau,
          n_cpu = design$sim_prm$n_cpu)]

      simsmok_cig_sc(dt, row_sel) # carry forward smok_cig if smok_status == 3
      dt[smok_cig_sc == 0L & smok_status_sc != "1", smok_cig_sc := 1L]
      dt[, mrk := NULL]


      dt[, smok_status_sc := factor(smok_status_sc)]
      dt[, smoke_cat_sc := 0L]
      dt[smok_status_sc == "3", smoke_cat_sc := 1L]
      dt[smok_status_sc == "4", smoke_cat_sc := 3L]
      dt[smok_status_sc == "4" & smok_cig_sc < 10L, smoke_cat_sc := 2L]
      dt[smok_status_sc == "4" & smok_cig_sc > 19L, smoke_cat_sc := 4L]
    }

    # ets ----
    if ("ets" %in% scenario_parms$sc_soc_qimd_rf_change) {
      if (!"ets_sc" %in% names(dt)) set(dt, NULL, "ets_sc", dt$ets_curr_xps)

      lutbl <-
        read_fst("./inputs/exposure_distributions/ets_table.fst", as.data.table = TRUE)
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, rank := pbinom(ets_curr_xps, 1, mu)]

      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, ets_sc := as.integer(rank < mu)]
    }

    # fv ----
    if ("fv" %in% scenario_parms$sc_soc_qimd_rf_change) {
      if (!"fruit_sc" %in% names(dt))
        set(dt, NULL, "fruit_sc", dt$fruit_curr_xps)
      if (!"veg_sc" %in% names(dt))
        set(dt, NULL, "veg_sc", dt$veg_curr_xps)


      lutbl <-
        read_fst("./inputs/exposure_distributions/frtpor_table.fst",
          as.data.table = TRUE)
      # is_valid_lookup_tbl(lutbl, c("year", "age", "sex", "sha", "qimd", "ethnicity"))
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, rank :=
          my_pZISICHEL(fruit_curr_xps / 80,
            mu,
            sigma,
            nu,
            tau,
            n_cpu = design$sim_prm$n_cpu)]
      # rn not uniformly distributed because it is discrete distr. That's expected
      # and without consequences.

      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, fruit_sc :=
          my_qZISICHEL(rank,
            mu, sigma, nu, tau, n_cpu = design$sim_prm$n_cpu) * 80L]  # g/d

      lutbl <-
        read_fst("./inputs/exposure_distributions/vegpor_table.fst",
          as.data.table = TRUE)
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, rank :=
          my_pDEL(veg_curr_xps / 80,
            mu, sigma, nu, n_cpu = design$sim_prm$n_cpu)]

      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, veg_sc :=
          my_qDEL(rank,
            mu, sigma, nu, n_cpu = design$sim_prm$n_cpu) * 80L]  # g/d
    }

    # alcohol ----
    if ("alc" %in% scenario_parms$sc_soc_qimd_rf_change) {
      if (!"alcohol_sc" %in% names(dt))
        set(dt, NULL, "alcohol_sc", dt$alcohol_curr_xps)

      lutbl <-
        read_fst("./inputs/exposure_distributions/alcohol_table.fst",
          as.data.table = TRUE)
      # is_valid_lookup_tbl(lutbl, c("year", "age", "sex", "sha", "qimd", "ethnicity"))
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, rank := pZINBI(alcohol_curr_xps, mu, sigma, nu)]

      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, alcohol_sc := qZINBI(rank, mu, sigma, nu)]
    }

    # active_days ----
    if ("pa" %in% scenario_parms$sc_soc_qimd_rf_change) {
      if (!"active_days_sc" %in% names(dt))
        set(dt, NULL, "active_days_sc", dt$active_days_curr_xps)

      dt[, c("mu", "sigma", "nu", "tau") := NULL]

      lutbl <-
        read_fst("./inputs/exposure_distributions/active_days_table.fst",
          as.data.table = TRUE)
      lookup_dt(dt, lutbl)
      dt[row_sel, rank := fcase(
        active_days_curr_xps == 0,
        pa0 - 1e-5,
        # -1e-5 for safety
        active_days_curr_xps == 1,
        pa1 - 1e-5,
        active_days_curr_xps == 2,
        pa2 - 1e-5,
        active_days_curr_xps == 3,
        pa3 - 1e-5,
        active_days_curr_xps == 4,
        pa4 - 1e-5,
        active_days_curr_xps == 5,
        pa5 - 1e-5,
        active_days_curr_xps == 6,
        pa6 - 1e-5,
        active_days_curr_xps == 7,
        1
      )]

      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = (paste0("pa", 0:6)))


      dt[row_sel, active_days_sc := (rank > pa0) + (rank > pa1) + (rank > pa2) +
          (rank > pa3) + (rank > pa4) + (rank > pa5) + (rank > pa6)]

      dt[, (paste0("pa", 0:6)) := NULL]
    }

    # bmi ----
    if ("bmi" %in% scenario_parms$sc_soc_qimd_rf_change) {
      if (!"bmi_sc" %in% names(dt))
        set(dt, NULL, "bmi_sc", dt$bmi_curr_xps)

      lutbl <-
        read_fst("./inputs/exposure_distributions/bmi_table.fst", as.data.table = TRUE)
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, rank :=
          my_pBCPEo(bmi_curr_xps, mu, sigma, nu, tau,
            n_cpu = design$sim_prm$n_cpu)]

      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, bmi_sc :=
          my_qBCPEo(rank, mu, sigma, nu, tau, n_cpu = design$sim_prm$n_cpu)]
    }

    # sbp ----
    if ("sbp" %in% scenario_parms$sc_soc_qimd_rf_change) {
      if (!"sbp_sc" %in% names(dt))
        set(dt, NULL, "sbp_sc", dt$sbp_curr_xps)

      lutbl <-
        read_fst("./inputs/exposure_distributions/sbp_table.fst", as.data.table = TRUE)
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, rank :=
          my_pBCPEo(sbp_curr_xps, mu, sigma, nu, tau,
            n_cpu = design$sim_prm$n_cpu)]

      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, sbp_sc :=
          my_qBCPEo(rank, mu, sigma, nu, tau, n_cpu = design$sim_prm$n_cpu)]
    }

    # tchol ----
    if ("tchol" %in% scenario_parms$sc_soc_qimd_rf_change) {
      if (!"tchol_sc" %in% names(dt))
        set(dt, NULL, "tchol_sc", dt$sbp_curr_xps)

      lutbl <-
        read_fst("./inputs/exposure_distributions/tchol_table.fst", as.data.table = TRUE)
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, rank :=
          my_pBCT(tchol_curr_xps, mu, sigma, nu, tau,
            n_cpu = design$sim_prm$n_cpu)]

      setnames(lutbl, "qimd", "qimd_sc")
      lookup_dt(dt, lutbl, exclude_col = c("mu", "sigma", "nu", "tau"))
      dt[row_sel, tchol_sc :=
          my_qBCT(rank, mu, sigma, nu, tau, n_cpu = design$sim_prm$n_cpu)]

      dt[, c("mu", "sigma", "nu", "tau", "rank") := NULL]
    }

    # case fatality ----

    if (any(
      scenario_parms$sc_soc_qimd_disease_fatality_change,
      scenario_parms$sc_soc_qimd_nonmodelled_fatality_change
    ) &&
        length(scenario_parms$sc_soc_qimd_rf_change) > 0L) {

      # only affect case fatalities for diseases which are linked to the selected exposures
      tt <-
        data.table(
          exposures = c(
            "tobacco",
            "packyears",
            "ets",
            "fv",
            "fruit",
            "alcohol",
            "pa",
            "sbp",
            "bmi",
            "tchol"
          ),
          short_exps = c(
            "tob",
            "tob",
            "ets",
            "fv",
            "fv",
            "alc",
            "pa",
            "sbp",
            "bmi",
            "tchol"
          )
        )
      tt <- tt[short_exps %in% scenario_parms$sc_soc_qimd_rf_change]

      causal_paths <- get_causal_paths(RR)
      affected_diseases <-
        causal_paths$dt[exposures %in% tt$exposures, unique(outcomes)]

      if (scenario_parms$sc_soc_qimd_disease_fatality_change) {
        nam <- affected_diseases[affected_diseases != "nonmodelled"]
        nam <- paste0("prb_", nam, "_mrtl")

        for (i in nam) {
          if (grepl("chd|stroke", i)) {
            lag <- design$lags_mc$cvd_lag
          } else if (grepl("copd", i)) {
            lag <- design$lags_mc$copd_lag
          } else if (grepl("_ca_mrtl$", i)) {
            lag <- design$lags_mc$cancer_lag
          }

          dt[, qimd_sc_lagged := # for lagged case fatalities
              fifelse(year >= (scenario_parms$sc_init_year - 2000L + lag),
                qimd_sc,
                qimd)]

          nc <- paste0(i, "_sc")
          setnames(dt, i, "mod____") # to avoid using get() due to performance issues
          tt <- dt[, min(mod____), keyby = .(age, qimd, sex, year)]
          setnames(dt, "mod____", i)

          lutbl <- tt[, {
            # because some combinations are missing
            l <- lapply(.SD, unique)
            setDT(expand.grid(l))
          }, .SDcols = c("age", "qimd", "sex", "year")]
          absorb_dt(lutbl, tt)

          setnames(lutbl, c("qimd", "V1"), c("qimd_sc_lagged", nc))

          lookup_dt(dt, lutbl, exclude_col = nc) # row_sel not appropriate here
          dt[, ("qimd_sc_lagged") := NULL]
        } # for loop
      } # if diseases fatalities in scenario

      if (scenario_parms$sc_soc_qimd_nonmodelled_fatality_change &&
          "nonmodelled" %in% affected_diseases) {
        i <- "p0_nonmodelled"
        lag <- design$lags_mc$nonmodelled_lag

        dt[, qimd_sc_lagged := # for lagged case fatalities
            fifelse(year >= scenario_parms$sc_init_year - 2000L + lag,
              qimd_sc,
              qimd)]

        nc <- paste0(i, "_sc")
        setnames(dt, i, "mod____") # to avoid using get() due to performance issues
        tt <- dt[, min(mod____), keyby = .(age, qimd, sex, year)]
        setnames(dt, "mod____", i)

        lutbl <- tt[, {
          # because some combinations are missing
          l <- lapply(.SD, unique)
          setDT(expand.grid(l))
        }, .SDcols = c("age", "qimd", "sex", "year")]
        absorb_dt(lutbl, tt)

        setnames(lutbl, c("qimd", "V1"), c("qimd_sc_lagged", nc))

        lookup_dt(dt, lutbl, exclude_col = nc) # row_sel not appropriate here
        dt[, ("qimd_sc_lagged") := NULL]
      }
    }

    return(invisible(dt))
  }
}



#' @export
run_scenario <-
  function(scenario_nam, # This is true_scenario names
           mc, # need to be mc_aggr
           dt,
           parameters_dt,
           design,
           output,
           timing = c(TRUE, FALSE)) {
    if (timing[[1]])
      ptm <- proc.time()

    basic_sc_nam <- parameters_dt[true_scenario == scenario_nam, unique(true_scenario), keyby = scenario]$scenario

    scenario_parms <- lapply(basic_sc_nam, fromGUI_scenario_parms, parameters_dt)
    names(scenario_parms) <- basic_sc_nam
    # Sort scenarios in chronological order (important for serial ensembles)
    basic_sc_nam <- names(sort(unlist(lapply(scenario_parms, `[[`, "sc_init_year"))))
    scenario_parms <- scenario_parms[basic_sc_nam]
    hlp <- list() # aux object to pass information between scenarios
    hlp$previous_invitees <- hlp$previous_elig <- hlp$previous_attendees <- rep(0L, nrow(dt$pop))

    # logic for parallel ensemble
    tt <- unlist(lapply(scenario_parms, `[[`, "sc_ens_parallel_prc"))
    if (length(tt) > 0) { # if parallel ensemble
      unique_pid <- dt$pop[, unique(pid)] # all pid

      if (sum(tt) == 1) {
        ttt <- sample(names(tt), length(unique_pid), TRUE, tt)
        hlp$sc_alloc <- lapply(names(tt), function(x) unique_pid[ttt==x])
        names(hlp$sc_alloc) <- names(tt)
      } else if (sum(tt) < 1) {
        tt <- c(tt, 1-sum(tt))
        names(tt) <- c(head(names(tt), -1), "excluded_")
        ttt <- sample(names(tt), length(unique_pid), TRUE, tt)
        hlp$sc_alloc <- lapply(names(tt), function(x) unique_pid[ttt==x])
        names(hlp$sc_alloc) <- names(tt)
      } else { # if > 1
        hlp$sc_alloc <- lapply(names(tt),
               function(x, tt) {
                 unique_pid[as.logical(rbinom(length(unique_pid), 1, tt[x]))]
        }, tt
        )
        names(hlp$sc_alloc) <- names(tt)
      }
      rm(unique_pid)
    }


    for (sc in basic_sc_nam) {
      set_eligible(scenario_parms[[sc]], dt$pop, hlp)
      set_invitees(scenario_parms[[sc]], dt$pop, hlp)
      set_attendees(scenario_parms[[sc]], dt$pop, scenario_nam, parameters_dt, design, hlp)
      set_px(scenario_parms[[sc]], dt$pop, mc, design) # slow
      set_lifestyle(scenario_parms[[sc]], dt$pop, design)
      set_structural(scenario_parms[[sc]], dt$pop, design)
      set_social(scenario_parms[[sc]], dt$pop, design)
    }

    dt$pop[, eligible_sc  := clamp(eligible_sc + hlp$previous_elig)]
    dt$pop[, invitees_sc  := clamp(invitees_sc + hlp$previous_invitees)]
    dt$pop[, attendees_sc := clamp(attendees_sc + hlp$previous_attendees)]
    dt$pop[, c("hc_eff_pa", "hc_eff_wm", "hc_eff_al", "hc_eff_sm") := NULL]
    # TODO I can calculate the effect of xps change to disease prb for
    # efficiency No need to recalculate disease probability for everyone only
    # apply disease impact on attendees (works only with kismet == TRUE)

    af_model(                 scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    htn_model(                scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    t2dm_model(               scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    chd_model(                scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    stroke_model(             scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    poststroke_dementia_model(scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    copd_model(               scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    lung_ca_model(            scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    colon_ca_model(           scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    breast_ca_model(          scenario_nam, mc, dt$pop, design, timing = timing[[2]])
    nonmodelled_model(        scenario_nam, mc, dt$pop, design, timing = timing[[2]])

    # if (!"smok_status_sc" %in% names(dt$pop))
    #   set(dt$pop, NULL, "smok_status_sc", dt$pop$smok_status_curr_xps)
    # if (!"smok_quit_yrs_sc" %in% names(dt$pop))
    #   set(dt$pop, NULL, "smok_quit_yrs_sc", dt$pop$smok_quit_yrs_curr_xps)
    # if (!"smok_dur_sc" %in% names(dt$pop))
    #   set(dt$pop, NULL, "smok_dur_sc", dt$pop$smok_dur_curr_xps)
    # if (!"smok_cig_sc" %in% names(dt$pop))
    #   set(dt$pop, NULL, "smok_cig_sc", dt$pop$smok_cig_curr_xps)

    output <- gen_output(scenario_nam, design$sim_prm, design$lags_mc, dt$pop, output)

    dt$pop[, (grep("_sc$", names(dt$pop), value = TRUE)) := NULL]

    if (timing[[1]])
      print(proc.time() - ptm)
    invisible(output)
  }


# Function for timing log
#' @export
time_mark <- function(x, file_nam = output_dir("times.txt")) {
  sink(
    file = file_nam,
    append = TRUE,
    type = "output",
    split = FALSE
  )
  cat(paste0(x, " at: ", Sys.time(), "\n"))
  sink()
}

#' @export
export_xps <- function(mc_,
                       dt,
                       write_to_disk = TRUE,
                       filenam = "val_xps_output.csv",
                       reweighted_to_hse = FALSE) {
  to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
  setnames(dt, "t2dm_prvl_curr_xps", "t2dm_prvl_original")
  if ("t2dm_prvl" %in% names(dt)) {
    dt[, t2dm_prvl_curr_xps := fifelse(t2dm_prvl == 0L, 0L, 1L)]
  } else {
    dt[, t2dm_prvl_curr_xps := fifelse(t2dm_prvl_original == 0L, 0L, 1L)]
  }
  dt[, smok_never_curr_xps := fifelse(smok_status_curr_xps == "1", 1L, 0L)]
  dt[, smok_active_curr_xps := fifelse(smok_status_curr_xps == "4", 1L, 0L)]
  if (reweighted_to_hse) {
    wt <- read_fst("./synthpop/hse_sociodemographics.fst",
                   as.data.table = TRUE)
    absorb_dt(dt, wt)
  } else {
    dt[, hse_wt := 1]
  }

  xps <- grep("_curr_xps$", names(dt), value = TRUE)
  xps <- xps[-which(xps == "smok_status_curr_xps")]
  out_xps <- groupingsets(
    dt,
    j = lapply(.SD, weighted.mean, hse_wt),
    by = c("year", "sex", "agegrp20", "qimd"),
    .SDcols = xps,
    sets = list(
      c("year", "sex", "agegrp20", "qimd"),
      c("year", "sex"),
      c("year", "agegrp20"),
      c("year", "qimd")
    )
  )[, `:=` (year = year + 2000L, mc = mc_)]
  for (j in seq_len(ncol(out_xps)))
    set(out_xps, which(is.na(out_xps[[j]])), j, "All")
  dt[, c(
    "agegrp20",
    "smok_never_curr_xps",
    "smok_active_curr_xps",
    "t2dm_prvl_curr_xps",
    "hse_wt"
  ) := NULL]
  setnames(dt, "t2dm_prvl_original", "t2dm_prvl_curr_xps")

  setkey(out_xps, year)
  if (write_to_disk)
    fwrite_safe(out_xps, output_dir(filenam))
  invisible(out_xps)
}

#' @export
export_mrtl <-
  function(mc_,
           dt,
           write_to_disk = TRUE,
           filenam = "val_mrtl_output.csv") {
    to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
    out <- groupingsets(
      dt,
      j = .N,
      by = c("year", "sex", "agegrp20", "qimd", "all_cause_mrtl"),
      sets = list(
        c("year", "sex", "agegrp20", "qimd", "all_cause_mrtl"),
        c("year", "sex", "all_cause_mrtl"),
        c("year", "agegrp20", "all_cause_mrtl"),
        c("year", "qimd", "all_cause_mrtl")
      )
    )[, `:=` (year = year + 2000L, mc = mc_)]

    out <-
      dcast(out, year + agegrp20 + sex + qimd + mc ~
              all_cause_mrtl, value.var = "N")
    setnafill(out, "c", 0L, cols = paste0(0:7))
    out[, pops := do_cols_dt(out, paste0(0:7), "+")]
    for (j in seq_len(ncol(out)))
      set(out, which(is.na(out[[j]])), j, "All")
    out[, `0` := NULL]
    setnames(out,
             paste0(1:7),
             paste0(
               "deaths_",
               c(
                 "nonmodelled",
                 "chd",
                 "stroke",
                 "copd",
                 "lung_ca",
                 "colon_ca",
                 "breast_ca"
               )
             ))

    dt[, c("agegrp20") := NULL]
    setkey(out, year)
    if (write_to_disk)
      fwrite_safe(out, output_dir(filenam))
    invisible(out)
  }

#' @export
export_incd <-
  function(mc_,
           dt,
           write_to_disk = TRUE,
           filenam = "val_incd_output.csv") {
    to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
    xps <- grep("_ca_prvl$", names(dt), value = TRUE)

    out <- groupingsets(
      dt,
      j = lapply(.SD, function(x)
        sum(x == 1L)),
      by = c("year", "sex", "agegrp20", "qimd"),
      .SDcols = xps,
      sets = list(
        c("year", "sex", "agegrp20", "qimd"),
        c("year", "sex"),
        c("year", "agegrp20"),
        c("year", "qimd")
      )
    )[, `:=` (year = year + 2000L, mc = mc_)]
    for (j in seq_len(ncol(out)))
      set(out, which(is.na(out[[j]])), j, "All")

    pop <- groupingsets(
      dt,
      j = .N,
      by = c("year", "sex", "agegrp20", "qimd"),
      sets = list(
        c("year", "sex", "agegrp20", "qimd"),
        c("year", "sex"),
        c("year", "agegrp20"),
        c("year", "qimd")
      )
    )[, `:=` (year = year + 2000L)]
    for (j in seq_len(ncol(pop)))
      set(pop, which(is.na(pop[[j]])), j, "All")
    setnames(pop, "N", "pops")
    absorb_dt(out, pop)
    xps2 <-  gsub("_prvl$", "", xps)
    xps2 <- paste0("cases_", xps2)
    setnames(out, xps, xps2)

    dt[, c("agegrp20") := NULL]
    setkey(out, year)
    if (write_to_disk)
      fwrite_safe(out, output_dir(filenam))
    invisible(out)
  }


#' @export
export_all_incd <-
  function(mc_,
           dt,
           write_to_disk = TRUE,
           filenam = "val_all_incd_output.csv") {
    to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
    xps <- grep("_prvl$", names(dt), value = TRUE)

    out <- groupingsets(
      dt,
      j = lapply(.SD, function(x)
        sum(x == 1L)),
      by = c("year", "sex", "agegrp20", "qimd"),
      .SDcols = xps,
      sets = list(
        c("year", "sex", "agegrp20", "qimd"),
        c("year", "sex"),
        c("year", "agegrp20"),
        c("year", "qimd")
      )
    )[, `:=` (year = year + 2000L, mc = mc_)]
    for (j in seq_len(ncol(out)))
      set(out, which(is.na(out[[j]])), j, "All")

    pop <- groupingsets(
      dt,
      j = .N,
      by = c("year", "sex", "agegrp20", "qimd"),
      sets = list(
        c("year", "sex", "agegrp20", "qimd"),
        c("year", "sex"),
        c("year", "agegrp20"),
        c("year", "qimd")
      )
    )[, `:=` (year = year + 2000L)]
    for (j in seq_len(ncol(pop)))
      set(pop, which(is.na(pop[[j]])), j, "All")
    setnames(pop, "N", "pops")
    absorb_dt(out, pop)
    xps2 <-  gsub("_prvl$", "", xps)
    xps2 <- paste0("cases_", xps2)
    setnames(out, xps, xps2)

    dt[, c("agegrp20") := NULL]
    setkey(out, year)
    if (write_to_disk)
      fwrite_safe(out, output_dir(filenam))
    invisible(out)
  }

#' @export
export_all_prvl <-
  function(mc_,
           dt,
           write_to_disk = TRUE,
           filenam = "val_all_prvl_output.csv") {
    to_agegrp(dt, 20L, 89L, "age", "agegrp20", to_factor = TRUE)
    xps <- grep("_prvl$", names(dt), value = TRUE)

    out <- groupingsets(
      dt,
      j = lapply(.SD, function(x)
        sum(x > 0L)),
      by = c("year", "sex", "agegrp20", "qimd"),
      .SDcols = xps,
      sets = list(
        c("year", "sex", "agegrp20", "qimd"),
        c("year", "sex"),
        c("year", "agegrp20"),
        c("year", "qimd")
      )
    )[, `:=` (year = year + 2000L, mc = mc_)]
    for (j in seq_len(ncol(out)))
      set(out, which(is.na(out[[j]])), j, "All")

    pop <- groupingsets(
      dt,
      j = .N,
      by = c("year", "sex", "agegrp20", "qimd"),
      sets = list(
        c("year", "sex", "agegrp20", "qimd"),
        c("year", "sex"),
        c("year", "agegrp20"),
        c("year", "qimd")
      )
    )[, `:=` (year = year + 2000L)]
    for (j in seq_len(ncol(pop)))
      set(pop, which(is.na(pop[[j]])), j, "All")
    setnames(pop, "N", "pops")
    absorb_dt(out, pop)
    xps2 <-  gsub("_prvl$", "", xps)
    xps2 <- paste0("cases_", xps2)
    setnames(out, xps, xps2)

    dt[, c("agegrp20") := NULL]
    setkey(out, year)
    if (write_to_disk)
      fwrite_safe(out, output_dir(filenam))
    invisible(out)
  }

#' @export
simulate_init_prvl <-
  function(mc,
           disease,
           design_ = design_$sim_prm,
           dt = POP) {
    tbl <- get_disease_epi_mc(mc, disease, "p", "v", design_$sim_prm$stochastic)
    incd <-
      get_disease_epi_mc(mc, disease, "i", "v", design_$sim_prm$stochastic)
    absorb_dt(tbl, incd)
    tbl[, prevalence := prevalence - incidence]
    tbl[, year := design_$sim_prm$init_year]
    # given init year is 13 and close enough to 11 which I have epi for
    col_nam <- setdiff(names(tbl), names(dt))
    absorb_dt(dt, tbl)
    nam <- paste0(disease, "_prvl")
    dt[year == design_$sim_prm$init_year, (nam) := as.integer(dqrunif(.N) < prevalence)]
    # faster than rbinom
    dt[year > design_$sim_prm$init_year &
         age == design_$sim_prm$ageL, (nam) := as.integer(dqrunif(.N) < prevalence)]
    # faster than rbinom
    setnafill(dt, "c", 0L, cols = nam)
    dt[, (col_nam) := NULL]
  }

#' @export
output_dir <- function(x = "") {
  file.path(design$sim_prm$output_dir, x)
}

#' @export
run_simulation <- function(parameters, design, final = FALSE) {
  # NOTE final = F for the exploratory runs
  # output_dir <<- function(x = "") {
  #   file.path(design$sim_prm$output_dir, x)
  # }

  on.exit(file.remove(output_dir("intermediate_out.csv")))

  parameters_dt <- fromGUI_to_dt(parameters)

  if (design$sim_prm$logs)
    time_mark("start parallelisation")

  # Parallelisation ----
  foreach(
    mc_iter = 1:(ifelse(final,
                     design$sim_prm$iteration_n_final,
                     design$sim_prm$iteration_n) *
                   design$sim_prm$n_synthpop_aggregation),
    .inorder = FALSE,
    .verbose = TRUE,
    .packages = c(
      "R6",
      "data.table",
      "IMPACTncdEnglmisc",
      "gamlss.dist", # For distr in prevalence.R
      "dqrng",
      "qs",
      "fst",
      "wrswoR",
      "CKutils"
    )
  ) %dopar% {
    setDTthreads(design$sim_prm$n_cpus)
    threads_fst(design$sim_prm$n_cpus)

    if (design$sim_prm$logs) {
      if (!dir.exists(output_dir("logs/"))) {
        dir.create(normalizePath(output_dir("logs/")), FALSE, TRUE)
      }
      sink(
        file = normalizePath(output_dir(
          paste0("logs/log", mc_iter, ".txt")
        ), mustWork = FALSE),
        append = TRUE,
        type = "output",
        split = FALSE
      )
    }
    mc_aggr <- as.integer(ceiling(mc_iter / design$sim_prm$n_synthpop_aggregation))
    design$get_lags(mc_aggr)

    POP <- SynthPop$new(mc_iter, design)
    if (!parameters$national_qimd_checkbox) {
      setnames(POP$pop, c("qimd", "lqimd"), c("nqimd", "qimd"))
    }

    # Run scenarios
    if (design$sim_prm$logs)
      print("scenario outputs")
    output_chunk <- list()
    output_chunk <- sapply(
      sort(fromGUI_scenario_names(parameters_dt)$true_scenario),
      # parameters_dt[grepl("^sc[1-9]", scenario), sort(unique(scenario))],
      run_scenario,
      mc_aggr,
      POP,
      parameters_dt,
      design,
      output_chunk,
      USE.NAMES = FALSE
    )


    rm(POP)
    # invisible(gc(verbose = FALSE, full = TRUE))

    invisible(lapply(output_chunk, setDT))
    invisible(lapply(output_chunk, function(x) {
      oldnam <- grep("_sc.*$", names(x), value = TRUE)
      newnam <- gsub("_sc.*$", "", oldnam)
      setnames(x, oldnam, newnam)
    }))
    output_chunk <- rbindlist(output_chunk, idcol = "scenario")

    setkey(output_chunk, scenario, pid, year) # for identify_longdead
    output_chunk <-
      output_chunk[output_chunk$year >= design$sim_prm$init_year_fromGUI &
                     between(output_chunk$age,
                             design$sim_prm$ageL,
                             design$sim_prm$ageH) &
                     !identify_longdead(all_cause_mrtl, pid_mrk), ]
    output_chunk[, pid_mrk  := mk_new_simulant_markers(pid)]
    output_chunk[, scenario := factor(scenario)]
    generate_health_econ(output_chunk, mc_aggr)
    output_chunk[, c("ncc", "pid", "income", "education") := NULL]

    # gen incd
    for (nam in grep("_prvl$", names(output_chunk), value = TRUE)) {
      newnam <- gsub("_prvl$", "_incd", nam)
      set(output_chunk,
          NULL,
          newnam,
          fifelse(output_chunk[[nam]] == 1L, 1L, 0L))
    }

    if (design$sim_prm$logs)
      print("transform prvl/dgn/mrtl to be summed")
    invisible(output_chunk[, lapply(.SD, fclamp_int, 0L, 1L, TRUE),
                           .SDcols = patterns("_prvl$|_dgn$|_mrtl$")])
    if ("lqimd" %in% names(output_chunk)) {
      output_chunk[, ("lqimd") := NULL]
    } else {
      output_chunk[, ("nqimd") := NULL]
    }

    output_chunk[, year := year + 2000L]
    to_agegrp(output_chunk, 20L, 89L, "age", "agegrp", TRUE, 30L) # TODO use case_when

    # Scale-up to ONS population projections
    # output_chunk[, sum(wt), keyby = .(year, scenario)]
    for (nam in grep("_prvl$|_dgn$|_incd|_mrtl$|_cost$|^eq5d",
                     names(output_chunk),
                     value = TRUE)) {
      set(output_chunk, NULL, nam, output_chunk[[nam]] * output_chunk$wt)
    } # TODO loop over numeric indices rather than col names

    # summarise by strata
    output_chunk[, c("age", "pid_mrk") := NULL]
    output_chunk <-
      output_chunk[, lapply(.SD, sum),
                   keyby = eval(design$sim_prm$strata_for_output)]
    setnames(output_chunk, "wt", "pops")
    output_chunk[, mc := mc_aggr]


    fwrite_safe(output_chunk, output_dir("intermediate_out.csv"))

    # write_fst(output_chunk, output_dir(
    #   paste0("intermediate_out_",
    #     mc_iter,
    #     ".fst"
    #   )
    # ), 100)

    rm(output_chunk)
    if (design$sim_prm$logs)
      sink()
    NULL
  }

  # End of parallelisation ----

  if (design$sim_prm$logs)
    time_mark("End of parallelisation")

  while (sink.number() > 0L)
    sink()

  output <- fread(output_dir("intermediate_out.csv"), stringsAsFactors = TRUE)

  tt <- fromGUI_scenario_names(parameters_dt)
  tt <- tt[, lapply(.SD, factor), .SDcols = is.character]
  setnames(tt, "true_scenario", "scenario")
  absorb_dt(output, tt)

  strata <- c("mc", design$sim_prm$strata_for_output, "friendly_name")
  output <- output[, lapply(.SD, sum), keyby = eval(strata)]

  # spread overhead policy costs to individuals
  tt <-
    parameters_dt[scenario != "global", .(scenario = unique(scenario)),
      keyby = c("true_scenario", "friendly_name")]

  for (sc_nam in tt$scenario) {
    l <- fromGUI_scenario_parms(sc_nam, parameters_dt)
    if (!l$sc_ens_is) {
      for (nam in grep("_ovrhd$", names(l), value = TRUE)) {
        newnam <- gsub("^sc_ls_|_cost_ovrhd$", "", nam)
        newnam <- paste0(newnam, "_ovrhd_cost")
        output[scenario == tt[scenario == sc_nam]$true_scenario,
          (newnam) := pops * l[[nam]] / sum(pops), by = year]
      }
    }
    if (l$sc_ens_serial_is) {
      for (nam in grep("_ovrhd$", names(l), value = TRUE)) {
        newnam <- gsub("^sc_ls_|_cost_ovrhd$", "", nam)
        newnam <- paste0(newnam, "_ovrhd_cost")
        output[scenario == tt[scenario == sc_nam]$true_scenario &
            between(year, l$sc_init_year, l$sc_last_year),
          (newnam) := pops * l[[nam]] / sum(pops), by = year]
      }
    }
    if (l$sc_ens_parallel_is) {
      for (nam in grep("_ovrhd$", names(l), value = TRUE)) {
        newnam <- gsub("^sc_ls_|_cost_ovrhd$", "", nam)
        newnam <- paste0(newnam, "_ovrhd_cost")
        # TODO better logic on how to spread overheads for parallel scenarios
        # currently it assumes parallel scenarios have identical overheads and
        # those of the last scenario ovewrite the previous overheads.
        output[scenario == tt[scenario == sc_nam]$true_scenario &
            between(year, l$sc_init_year, l$sc_last_year),
          (newnam) := pops * l[[nam]] / sum(pops), by = year]
      }
    }
  }

  # calculate useful cost indices
  output[, policy_cost := invitation_cost + attendees_cost +
           active_days_cost + bmi_cost + alcohol_cost + smoking_cost +
           alcoholreduc_ovrhd_cost + smkcess_ovrhd_cost +
           wghtloss_ovrhd_cost + pa_ovrhd_cost]

  setkey(output, year, friendly_name)

  # calculate net cost/effect
  tt <-
    output[scenario == fromGUI_baseline_scenario(parameters_dt)]

  strata <- setdiff(strata, c("scenario", "friendly_name"))

  # ensure baseline has all combinations of strata observed in other scenarios.
  # Otherwise NAs occur from the subtraction
  baseline <- CJ(mc = output$mc, year = output$year, agegrp = output$agegrp,
    sex = output$sex, qimd = output$qimd, ethnicity = output$ethnicity,
    unique = TRUE)
  absorb_dt(baseline, tt)
  baseline[, c("scenario") := NULL]
  setnafill(baseline, "const", 0, cols = names(baseline)[sapply(baseline, is.numeric)])
  rm(tt)

  output[baseline, on = strata, `:=` (
    net_utility = eq5d - i.eq5d,
    net_policy_cost = policy_cost - i.policy_cost,
    net_healthcare_cost = healthcare_cost - i.healthcare_cost,
    net_socialcare_cost = socialcare_cost - i.socialcare_cost,
    net_informal_care_cost = informal_care_cost - i.informal_care_cost,
    net_productivity_cost = productivity_cost - i.productivity_cost,
    # (higher is better)
    cpp_cvd = i.cvd_incd - cvd_incd,
    cpp_chd = i.chd_incd - chd_incd,
    cpp_stroke = i.stroke_incd - stroke_incd,
    cpp_poststroke_dementia = i.poststroke_dementia_incd - poststroke_dementia_incd,
    cpp_copd = i.copd_incd - copd_incd,
    cpp_lung_ca = i.lung_ca_incd - lung_ca_incd,
    cpp_colon_ca = i.colon_ca_incd - colon_ca_incd,
    cpp_breast_ca = i.breast_ca_incd - breast_ca_incd,
    cpp_htn = i.htn_incd - htn_incd,
    cpp_af = i.af_incd - af_incd,
    cpp_t2dm = i.t2dm_incd - t2dm_incd,

    cypp_cvd = i.cvd_prvl - cvd_prvl,
    cypp_chd = i.chd_prvl - chd_prvl,
    cypp_stroke = i.stroke_prvl - stroke_prvl,
    cypp_poststroke_dementia = i.poststroke_dementia_prvl - poststroke_dementia_prvl,
    cypp_copd = i.copd_prvl - copd_prvl,
    cypp_lung_ca = i.lung_ca_prvl - lung_ca_prvl,
    cypp_colon_ca = i.colon_ca_prvl - colon_ca_prvl,
    cypp_breast_ca = i.breast_ca_prvl - breast_ca_prvl,
    cypp_htn = i.htn_prvl - htn_prvl,
    cypp_af = i.af_prvl - af_prvl,
    cypp_t2dm = i.t2dm_prvl - t2dm_prvl,

    dpp_nonmodelled = i.nonmodelled_mrtl - nonmodelled_mrtl,
    dpp_chd = i.chd_mrtl - chd_mrtl,
    dpp_stroke = i.stroke_mrtl - stroke_mrtl,
    dpp_copd = i.copd_mrtl - copd_mrtl,
    dpp_lung_ca = i.lung_ca_mrtl - lung_ca_mrtl,
    dpp_colon_ca = i.colon_ca_mrtl - colon_ca_mrtl,
    dpp_breast_ca = i.breast_ca_mrtl - breast_ca_mrtl,
    dpp_all_cause = i.all_cause_mrtl - all_cause_mrtl
  )]

  output[, `:=` (
    total_hcp_cost = net_policy_cost + net_healthcare_cost,
    # healthcare perspective
    total_hscp_cost = net_policy_cost + net_healthcare_cost +  net_socialcare_cost,
    societal_cost = net_policy_cost + net_healthcare_cost +  net_socialcare_cost +
      net_informal_care_cost - net_productivity_cost
  )]

  output[, mc := factor(mc)] # helpful later when I calculate cumsums
  write_fst(output, output_dir("results.fst"), 100)
  return(output)
}

# Function is borrowed from the shinythemes:::allThemes(), because it is not
# exported there
allthemes <- function() {
  themes <- dir(system.file("shinythemes/css", package = "shinythemes"),
                "*.min.css")
  sub(".min.css", "", themes)
}

# from https://stackoverflow.com/questions/47827337/prodution-ready-themeselector-for-shiny
#' @export
mythemeSelector <- function() {
  div(
    div(
      # tags$style(type='text/css', ".shiny-input-container { height: 20px; line-height: 0.86; font-size: 12px; margin: 0px; padding: 0px; border: 0px;}"),
      #
      # tags$style(type='text/css', ".form-control { height: 20px; line-height: 1.1; font-size: 12px; margin-left: 20px; padding: 0px; border: 0px;}"),
      selectInput("shinytheme_selector", "Please select another theme",
                  c("default", allthemes()),
                  selectize = FALSE, width = "100%"
      )
    ),
    tags$script(
      "$('#shinytheme_selector')
        .on('change', function(el) {
        var allThemes = $(this).find('option').map(function() {
        if ($(this).val() === 'default')
        return 'bootstrap';
        else
        return $(this).val();
        });
        // Find the current theme
        var curTheme = el.target.value;
        if (curTheme === 'default') {
        curTheme = 'bootstrap';
        curThemePath = 'shared/bootstrap/css/bootstrap.min.css';
        } else {
        curThemePath = 'shinythemes/css/' + curTheme + '.min.css';
        }
        // Find the <link> element with that has the bootstrap.css
        var $link = $('link').filter(function() {
        var theme = $(this).attr('href');
        theme = theme.replace(/^.*\\//, '').replace(/(\\.min)?\\.css$/, '');
        return $.inArray(theme, allThemes) !== -1;
        });
        // Set it to the correct path
        $link.attr('href', curThemePath);
        });"
    )
  )
}

#' Get the Recent System Load
#'
#' @return A named numeric vector with five non-negativeelements `1min`,
#'   `5min`, and `15min` average CPU utilisation, used RAM, and available RAM.
#' The first values represent estimates of the CPU load during the last
#' minute, the last five minutes, and the last fifteen minutes \[1\]. An idle
#' system have values close to zero, and a heavily loaded system have values
#' near one`. If they are unknown, missing values are returned. Only works for
#' Linux systems.
#'
#' @details
#' This function works only Unix-like system with \file{/proc/loadavg}. It is
#' heavily based on parallely::cpuLoad
#' (\url{https://github.com/HenrikBengtsson/parallelly})
#'
#' @references
#' 1. Linux Load Averages: Solving the Mystery,
#'    Brendan Gregg's Blog, 2017-08-08,
#'    \url{http://www.brendangregg.com/blog/2017-08-08/linux-load-averages.html}
#'
#' @keywords internal
#' @export
sysLoad <- function() {
  if (file.exists("/proc/loadavg")) {
    res <- readLines("/proc/loadavg", n = 1L)
    res <- strsplit(res, split=" ", fixed = TRUE)[[1]]
    res <- as.numeric(res[1:3])
    res <- signif(100 * res/parallel::detectCores(), 2)
    x <- system2('free', args = '-m', stdout = TRUE) # only for linux
    x <- strsplit(x[2], " +")[[1]][3:4]
    x <- round(as.numeric(x)/1024, 2)
    res <- as.numeric(c(res, x))
  } else {
    res <- rep(NA_real_, times = 5L)
  }
  names(res) <- c("1minAvgCPU(%)", "5minAvgCPU(%)", "15minAvgCPU(%)", "UsedRAM(Gb)", "FreeRAM(Gb)")
  res
}

