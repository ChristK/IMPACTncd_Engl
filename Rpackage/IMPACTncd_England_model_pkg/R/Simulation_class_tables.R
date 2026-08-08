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
# Simulation class table export methods
# This file adds table export methods to the Simulation class using $set()
# -----------------------------------------------------------------------------


# safe_fquantile_byid ----
# Safe wrapper for CKutils::fquantile_byid that handles empty vectors.
# The underlying C++ function crashes with a segfault when passed empty vectors.
# This wrapper returns an empty data.table with the expected structure instead.
safe_fquantile_byid <- function(x, q, id, rounding = FALSE) {
  # Handle empty vector case - return empty data.table with expected structure
  if (length(x) == 0L) {
    result <- data.table::data.table(id = character(0))
    for (qi in q) {
      result[, (as.character(qi)) := numeric(0)]
    }
    return(result)
  }

  # Call the original function for non-empty vectors
  fquantile_byid(x, q, id, rounding)
}


# Discounting helpers ---------------------------------------------------------
# Costs and QALYs are reported at several discount levels side by side: every
# affected table (qalys, net_qalys, costs, net_costs and the CEA tables) gains a
# `discount` column and carries one block of rows per level, so the undiscounted
# (0%) figures and the discounted ones live in the same file.
#
# The equity slope-index tables are deliberately NOT expanded: their metric set
# mixes discountable QALYs with plain event counts (CPP / CYPP / DPP), which
# have no present value, so one `discount` column could not describe the file.


# fmt_discount_pct ----
# Rate -> the label used in the `discount` column ("0%", "3.5%", "1.5%").
# Trailing zeros are dropped so the label round-trips the rate that was asked
# for.
fmt_discount_pct <- function(rate) {
  paste0(format(rate, trim = TRUE, scientific = FALSE, drop0trailing = TRUE), "%")
}


# build_discount_levels ----
# Normalises the user-supplied rates into one row per discount level.
# `qaly_discount_rate` and `cost_discount_rate` are paired element-wise (a
# length-1 rate is recycled against a longer one), so the default `c(0, 3.5)`
# for both gives two levels: undiscounted, and NICE's 3.5% on QALYs and costs
# alike. Differential rates are supported by passing equal-length vectors, e.g.
# `qaly_discount_rate = c(0, 1.5)` with `cost_discount_rate = c(0, 3.5)`.
build_discount_levels <- function(qaly_discount_rate, cost_discount_rate) {
  q <- as.numeric(qaly_discount_rate)
  cst <- as.numeric(cost_discount_rate)
  if (length(q) == 0L) q <- 0
  if (length(cst) == 0L) cst <- 0

  n <- max(length(q), length(cst))
  # The contract is EQUAL lengths, or one length-1 rate recycled against the
  # other -- exactly what the error below states. A divisibility test
  # (`n %% length(q) != 0L`) is NOT that contract: it also admits 4-vs-2 and
  # 6-vs-3, where rep_len() then invents pairings the caller never asked for.
  # `build_discount_levels(c(0, 1.5, 3.5, 5), c(0, 3.5))` recycled the costs to
  # c(0, 3.5, 0, 3.5) and published levels "QALYs 3.5%/costs 0%" and
  # "QALYs 5%/costs 3.5%" -- discounting health but not money, which is not a
  # coherent economic analysis. It also made rejection non-monotonic in the
  # input: 4-vs-2 was accepted while 3-vs-2 errored.
  if (length(q) != length(cst) && length(q) != 1L && length(cst) != 1L) {
    stop("qaly_discount_rate (length ", length(q), ") and cost_discount_rate ",
         "(length ", length(cst), ") must be the same length, or one of them ",
         "must be a single rate to recycle.", call. = FALSE)
  }
  q <- rep_len(q, n)
  cst <- rep_len(cst, n)

  if (!all(is.finite(q)) || !all(is.finite(cst)) || any(c(q, cst) <= -100)) {
    stop("Discount rates must be finite percentages greater than -100.",
         call. = FALSE)
  }

  lvl <- unique(data.table(qaly = q, cost = cst), by = c("qaly", "cost"))
  lvl[, `:=`(qaly_label = fmt_discount_pct(qaly),
             cost_label = fmt_discount_pct(cost))]
  # The CEA tables mix both rates, so their label spells the pair out whenever
  # the two differ.
  lvl[, label := fifelse(qaly == cost, qaly_label,
                         paste0("QALYs ", qaly_label, "/costs ", cost_label))]
  lvl[]
}


# promote_short_year ----
# Year ARGUMENTS accept a two-digit shorthand, because the summaries store short
# years (19 rather than 2019) and callers naturally pass what they see. A value
# <= 100 is promoted by 2000; anything else is returned untouched, so a real
# year, a NULL or a non-numeric passes through unchanged.
#
# This exists as a helper because the promotion has to happen for EVERY year
# argument or not at all. `discount_from_year` was missed when multi-level
# discounting landed, and the failure is silent rather than loud: an unpromoted
# 19 makes the discount exponent `year - 19` ~ 2000, 1/1.035^2000 underflows to
# 0, and every discounted figure is published as zero with no error.
promote_short_year <- function(year) {
  if (is.numeric(year) && length(year) == 1L && !is.na(year) && year <= 100) {
    return(year + 2000L)
  }
  year
}


# discount_factor ----
# Present-value factor: PV = FV * discount_factor(year, rate, from_year), i.e.
# 1 / (1 + rate/100)^max(0, year - from_year). Years at or before `from_year`
# are undiscounted, and a 0% rate leaves every year unchanged.
#
# `year` must be in FULL format (2019, not 19) - as it is by the time
# tbl_smmrs_core() and export_cea_tables() see it, both of which promote it.
discount_factor <- function(year, rate, from_year) {
  1 / (1 + rate / 100)^pmax(0, year - from_year)
}


# expand_discount_levels ----
# Replicates a per-year long table once per discount level, scaling `value_cols`
# by that level's present-value factor and tagging the copy with a `discount`
# label column. `rate_for` selects which of the paired rates applies: QALY
# tables discount at the QALY rate and cost tables at the cost rate, so levels
# sharing a rate collapse to a single block of rows.
#
# Scaling here happens on the per-year values, i.e. before any cumsum(), so the
# cumulative columns accumulate present values rather than discounting a total.
#
# CALLERS: add "discount" to the id.vars of the melt() that follows. Those melts
# name id.vars but not measure.vars, so a column left out of id.vars is absorbed
# as a measure variable - the same trap the `.cmp__` note in tbl_smmrs_core
# warns about.
expand_discount_levels <- function(d, value_cols, levels, from_year,
                                   rate_for = c("qaly", "cost")) {
  rate_for <- match.arg(rate_for)
  if (!"year" %in% names(d)) {
    stop("Discounting needs a `year` column; every stratum must include \"year\".",
         call. = FALSE)
  }

  rates <- levels[[rate_for]]
  labels <- levels[[paste0(rate_for, "_label")]]
  keep <- !duplicated(rates)
  rates <- rates[keep]
  labels <- labels[keep]

  rbindlist(lapply(seq_along(rates), function(i) {
    di <- copy(d)
    fctr <- discount_factor(di[["year"]], rates[i], from_year)
    for (cl in value_cols) set(di, NULL, cl, di[[cl]] * fctr)
    set(di, NULL, "discount", labels[i])
    di
  }))
}


# calc_equity_slope_indices ----
# Absolute and relative equity slope indices (analogues of the Slope Index of
# Inequality, SII, and Relative Index of Inequality, RII) for a benefit
# distributed across ordered socioeconomic groups (DIMD deciles).
#
# Statistical framework: each socioeconomic group is assigned a
# population-weighted midpoint cumulative rank (a *ridit* score) in [0, 1],
# ordered from LEAST deprived (r -> 0) to MOST deprived (r -> 1), and the
# absolute equity index is the population-WEIGHTED OLS slope of the PER-CAPITA
# benefit on that score -- exactly the SII. The grouped-data weighted
# least-squares estimator this implements is Kakwani, Wagstaff & van Doorslaer
# 1997 (J Econometrics 77(1):87-103, Eqs 3-4 and 6): "the grouped nature of the
# data calls for Weighted Least Squares". The ridit/structured-regression setup
# follows Moreno-Betancur et al. 2015 (Epidemiology 26(4):518-527), and applying
# the index to a cumulative *modelled benefit* rather than a health level
# follows Cookson et al. 2026 (PharmacoEconomics 44(3):317-328, Eq. 5).
#
# The closed-form weighted slope beta = Sxy / Sxx is identical to
# `coef(lm(y ~ ridit, weights = N))["ridit"]` but is computed by pure
# data.table aggregation so it can be evaluated over millions of
# (mc, scenario, year, ...) groups without a per-group lm() call.
#
# Why it is legitimate to fit this to a BENEFIT (a difference between two
# scenarios) rather than to a health level: beta is exactly linear in y, with
# coefficients c_k = N_k (r_k - xbar) / Sxx that depend only on the group sizes
# and sum to zero. So with the ridit and the weights held FIXED across arms --
# which the comparator reference population guarantees -- SII(y_post) -
# SII(y_pre) = SII(y_post - y_pre) identically, and a benefit that is equal per
# head in every group scores exactly zero.
#
# Sign convention: all supported metrics (CPP, CYPP, DPP, net QALYs) are
# "goods", so a POSITIVE index means the most deprived gain more (pro-poor,
# inequality-reducing); a negative index is pro-rich (inequality-widening).
#
# @param d  A data.table with (at least) columns `rank` (integer deprivation
#   rank, 1 = MOST deprived .. K = LEAST deprived, i.e. `as.integer(dimd)`),
#   `B` (cumulative benefit for that group) and `N` (reference population for
#   that group, used both as the regression weight and the ridit denominator),
#   plus the grouping columns named in `by`.
# @param by  Character vector of columns identifying one gradient per group
#   (must include the Monte-Carlo id `mc`; DIMD is NOT in `by` -- it is the
#   gradient axis consumed into the index).
# @return A data.table with one row per `by` group and columns:
#   AEI_total (the SII scaled to the whole reference population, SII * n: the
#     HYPOTHETICAL total benefit difference if every person received the fitted
#     most-deprived per-capita benefit instead of the fitted least-deprived one.
#     It is NOT the benefit difference between the extreme deciles themselves,
#     which is roughly an order of magnitude smaller);
#   AEI_per100k (the SII per 100,000 people alive in the reporting year);
#   total_benefit (the total undiscounted benefit the gradient was fit to --
#     reported because the absolute indices are translation-invariant and so
#     cannot distinguish levelling up from levelling down);
#   REI_rel (SII / population-weighted mean benefit; NA unless that mean is
#     POSITIVE);
#   RII_ratio (fitted benefit at the most-deprived end / least-deprived end; NA
#     unless the fitted benefit is positive at BOTH ends);
#   fit_R2 (population-weighted R^2 of the straight-line fit -- how much of the
#     between-group variation the gradient actually describes).
calc_equity_slope_indices <- function(d, by) {
  d <- copy(d)
  # Order least -> most deprived within each group (rank descending) so the
  # cumulative population share runs from the bottom to the top of the ridit.
  setorderv(d, c(by, "rank"), order = c(rep(1L, length(by)), -1L))
  d[, .totN := sum(N), by = by]
  d[, .csum := cumsum(N), by = by]
  d[, ridit := (.csum - N / 2) / .totN]

  agg <- d[, .(
    Sw   = sum(N),
    Swx  = sum(N * ridit),
    Swy  = sum(B),
    Swxx = sum(N * ridit * ridit),
    Swxy = sum(ridit * B),
    # sum(N * y^2) with y = B/N, for the weighted total sum of squares. Groups
    # of zero population carry zero weight, so they contribute nothing rather
    # than dividing by zero: this function has no `N > 0` precondition of its
    # own (export_equity_tables() enforces one, the tinytest callers do not).
    Swyy = sum(fifelse(N > 0, B * B / N, 0))
  ), by = by]

  agg[, xbar := Swx / Sw]
  agg[, ybar := Swy / Sw]                       # population-weighted mean benefit
  agg[, Sxx := Swxx - Sw * xbar * xbar]
  agg[, Sxy := Swxy - Sw * xbar * ybar]
  agg[, Syy := Swyy - Sw * ybar * ybar]
  agg[, beta := fifelse(Sxx > 0, Sxy / Sxx, NA_real_)]
  agg[, fit0 := ybar - beta * xbar]             # fitted per-capita benefit at r = 0 (least deprived)
  agg[, fit1 := ybar + beta * (1 - xbar)]       # fitted per-capita benefit at r = 1 (most deprived)

  agg[, AEI_total     := beta * Sw]
  agg[, AEI_per100k   := beta * 1e5]
  agg[, total_benefit := Swy]

  # Both relative indices need a POSITIVE fitted benefit, not merely a non-zero
  # one. beta / ybar with ybar < 0 flips sign, and fit1 / fit0 with BOTH ends
  # negative inverts the ordering -- so a net-harm scenario in which the most
  # deprived lose least (pro-poor under the sign convention above) would be
  # published as pro-rich. The earlier `ybar != 0` and
  # `sign(fit0) == sign(fit1)` guards both pass in exactly that case. Every
  # metric here is a signed contrast (a disinvestment scenario is legitimate),
  # so this is reachable, not hypothetical. A relative index has no meaning on a
  # variable that is not positive (O'Donnell et al. 2008, World Bank ch. 8: "if
  # the mean of the variable were 0, the concentration index would not be
  # defined"), so return NA rather than an inverted number. The ABSOLUTE indices
  # are unaffected -- they stay correct and signed throughout.
  agg[, REI_rel   := fifelse(!is.na(beta) & ybar > 0, beta / ybar, NA_real_)]
  agg[, RII_ratio := fifelse(!is.na(beta) & fit0 > 0 & fit1 > 0,
                             fit1 / fit0, NA_real_)]

  # Population-weighted R^2 of the straight line. The whole index family assumes
  # the per-capita benefit is linear in the ridit, and non-linearity is the
  # SII's principal known vulnerability (Sergeant & Firth 2006), so publish how
  # well the line actually fits. Unlike the other outputs this is NOT invariant
  # under the decile -> quintile collapse.
  agg[, fit_R2 := fifelse(!is.na(beta) & Sxx > 0 & Syy > 0,
                          (Sxy * Sxy) / (Sxx * Syy), NA_real_)]

  agg[, .SD, .SDcols = c(by, "AEI_total", "AEI_per100k", "total_benefit",
                         "REI_rel", "RII_ratio", "fit_R2")]
}


# Deprivation axes the model reports on, with their canonical factor levels
# (repeated across the package, e.g. Simulation_class.R, Disease_class.R,
# ExposureEffect_class.R). Order matters: as.integer() on these gives the
# deprivation rank 1 = most deprived .. K = least deprived.
.deprivation_vars <- c("dimd", "qimd")
.deprivation_levels <- list(
  dimd = c("1 most deprived", as.character(2:9), "10 least deprived"),
  qimd = c("1 most deprived", as.character(2:4), "5 least deprived")
)

# Standard decile -> quintile collapse (deciles 1,2 -> quintile 1, and so on),
# matching `Simulation_class.R` run_sim() and `Disease_class.R`.
.dimd_to_qimd <- c(1L, 1L, 2L, 2L, 3L, 3L, 4L, 4L, 5L, 5L)


# deprivation_rank ----
# Deprivation rank (1 = most deprived .. K = least deprived) for a `dimd`/`qimd`
# column. factor() silently returns NA for labels outside `levels`, and callers
# use the rank only to order or collapse, so an unexpected level set would yield
# a plausible-but-wrong result rather than any error (decile and quintile labels
# even share four of five labels). Fail loudly instead.
deprivation_rank <- function(x, var) {
  levs <- .deprivation_levels[[var]]
  r <- as.integer(factor(as.character(x), levels = levs))
  if (anyNA(r)) {
    stop(sprintf(
      paste0("unexpected `%s` level(s) in the summaries: %s. Expected exactly: ",
             "%s. The deprivation rank drives the ordering and the decile -> ",
             "quintile collapse, so an unrecognised level would silently ",
             "produce a wrong result."),
      var,
      paste(sQuote(unique(as.character(x)[is.na(r)])), collapse = ", "),
      paste(sQuote(levs), collapse = ", ")
    ), call. = FALSE)
  }
  r
}


# umbrella_disease_names ----
# Diseases whose incidence is declared `type: 0` in the design YAML, i.e. that
# are defined purely as the union of other modelled diseases. In the stock
# England design these are `dm` (= t1dm + t2dm), `ctdra` (= ctd + ra) and
# `cancer` (= the five site-specific cancers).
#
# They are legitimate output rows in their own right, but they must never enter
# a sum that already contains their components: `all_diseases_sum` is a sum of
# event counts that deliberately carries comorbidity multiplicity (two DIFFERENT
# diseases in one person count twice), which is a documented property -- but
# counting one person's lung cancer once as `lung_ca` and again as `cancer` is
# the same event twice, which is not. On the stock design the three umbrellas
# are ~9.5% of the summed events, and because they concentrate in three disease
# families with their own deprivation gradients they distort the SLOPE, not just
# the level.
#
# The rule is read off the design rather than hard-coded, so it stays correct
# for `sim_design_elast.yaml` and for user-supplied YAMLs. A disease counts as
# an umbrella only when at least one of its components is actually present in
# the summary being summed -- otherwise it is the only representative of its
# family and excluding it would lose those events entirely.
#
# @param diseases  `self$design$sim_prm$diseases` (a list of disease configs).
#   NULL or empty gives character(0), so callers that mock the design out (see
#   inst/tinytest/test_export_equity_tables.R) behave exactly as before.
# @param present  Character vector of disease names present in the data.
umbrella_disease_names <- function(diseases, present) {
  if (!is.list(diseases) || length(diseases) == 0L) return(character(0))
  out <- vapply(diseases, function(x) {
    ty  <- x[["meta"]][["incidence"]][["type"]]
    dep <- x[["meta"]][["incidence"]][["influenced_by_disease_name"]]
    if (length(ty) == 1L && !is.na(ty) && ty == 0L &&
        length(dep) > 0L && any(as.character(dep) %chin% present)) {
      as.character(x[["name"]])
    } else {
      NA_character_
    }
  }, character(1))
  unname(out[!is.na(out)])
}


# add_qimd_from_dimd ----
# Adds a `qimd` (IMD quintile) column derived from `dimd` (IMD decile), in
# place, so ANY table family can be stratified by quintiles without `qimd`
# having been named in `strata_for_output`.
#
# This is a relabelling only -- the actual collapsing is then done by whatever
# group-by the caller already performs, which is what makes it exact: summing
# two deciles' counts (and popsize) reproduces precisely what aggregating by
# quintile at summary time would have given, and rates built as
# summed-numerator / summed-denominator follow automatically. The converse is
# impossible: deciles cannot be recovered from quintiles.
#
# No-op (returning FALSE) when `dt` is NULL, already has `qimd`, or has no
# `dimd` to derive from. Adding the column is harmless when unused: every
# aggregation in this file selects value columns by anchored pattern, so `qimd`
# is only ever carried into a result when it is named in the strata.
add_qimd_from_dimd <- function(dt) {
  if (is.null(dt) || "qimd" %in% names(dt) || !"dimd" %in% names(dt)) {
    return(invisible(FALSE))
  }
  lv <- .deprivation_levels$qimd
  dt[, qimd := factor(lv[.dimd_to_qimd][deprivation_rank(dimd, "dimd")],
                      levels = lv)]
  invisible(TRUE)
}


# =============================================================================
# Two-age-group collapse (`two_agegrps`)
# -----------------------------------------------------------------------------
# `export_tables(two_agegrps = TRUE)` reports `agegrp` as TWO coarse groups
# instead of the design's 5-year bands, writing to `tables2agegrps/`. The split
# is configurable via `two_agegrps_split_age` (default 65L, the historical
# 30-64 / 65-99).
#
# This is a RELABELLING only -- exactly like add_qimd_from_dimd() -- so the
# actual collapsing is done by whatever group-by the caller already performs.
# That is what makes it exact for counts: summing the constituent bands
# reproduces what banding at summary time would have given, and rates built as
# summed-numerator / summed-denominator follow. (The equity family is the one
# exception: it FITS A MODEL, so see export_equity_tables().)
#
# The band bounds are read from the DATA, never hard-coded to 30..99: `ageL` /
# `ageH` are design parameters, and a run reporting 20-89 must not be
# relabelled "30-64" / "65-99".
# =============================================================================

# .agegrp_bounds ----
# Numeric bounds of age-band labels ("30-34" -> lo 30, hi 34; an open-ended
# "90+" -> lo 90, hi NA). Fails loudly on anything unparseable rather than
# returning NA and dropping those rows into neither group -- same reasoning as
# deprivation_rank().
.agegrp_bounds <- function(labs, col = "agegrp") {
  labs <- as.character(labs)
  open <- grepl("^[0-9]+\\+$", labs)
  lo <- suppressWarnings(as.integer(sub("^([0-9]+).*$", "\\1", labs)))
  hi <- suppressWarnings(as.integer(sub("^[0-9]+-([0-9]+)$", "\\1", labs)))
  hi[open] <- NA_integer_
  bad <- is.na(lo) | (is.na(hi) & !open)
  if (any(bad)) {
    stop(sprintf(paste0(
      "unparseable `%s` level(s): %s. Expected labels like '30-34', or an ",
      "open-ended '95+'. The two-age-group collapse splits on these bounds, so ",
      "an unrecognised label would silently land in the wrong group."),
      col, paste(sQuote(unique(labs[bad])), collapse = ", ")), call. = FALSE)
  }
  list(lo = lo, hi = hi, open = open)
}


# two_agegrp_map ----
# Map every age-band label to one of two coarse groups, split at `split_age` --
# which becomes the FIRST year of the UPPER group. Returns a named character
# vector (label -> coarse label) carrying the ordered coarse levels in the
# "levels_ordered" attribute.
two_agegrp_map <- function(labs, split_age = 65L, col = "agegrp") {
  labs <- unique(as.character(labs[!is.na(labs)]))
  b <- .agegrp_bounds(labs, col)
  # The split must fall BETWEEN bands, never inside one. Test CONTAINMENT
  # rather than "is it a band start?": a small locality can easily have an
  # empty 5-year band, and a run whose data happen not to contain the band
  # beginning exactly at `split_age` is still perfectly splittable there.
  cut <- b$lo < split_age & (is.na(b$hi) | split_age <= b$hi)
  if (any(cut)) {
    stop(sprintf(paste0(
      "`two_agegrps_split_age` = %s falls INSIDE the `%s` band(s) %s, which the ",
      "collapse would have to cut in half. Choose an age that starts a band -- ",
      "those present start at: %s."),
      split_age, col, paste(sQuote(labs[cut]), collapse = ", "),
      paste(sort(unique(b$lo)), collapse = ", ")), call. = FALSE)
  }
  lower <- b$lo < split_age
  # A split outside the range gives ONE group, not two. That is always a
  # mistake under an argument called `two_agegrps`, and it would also make the
  # empty group's label undefined (max() over nothing).
  if (!any(lower) || all(lower)) {
    stop(sprintf(paste0(
      "`two_agegrps_split_age` = %s leaves one of the two age groups EMPTY: ",
      "every `%s` band present falls %s it. The bands run %s-%s."),
      split_age, col, if (all(lower)) "below" else "at or above", min(b$lo),
      if (any(b$open)) paste0(max(b$lo), "+") else max(b$hi)), call. = FALSE)
  }
  lo_lab <- paste0(min(b$lo), "-", split_age - 1L)
  hi_lab <- if (any(b$open[!lower])) {
    paste0(split_age, "+")
  } else {
    paste0(split_age, "-", max(b$hi[!lower]))
  }
  res <- fifelse(lower, lo_lab, hi_lab)
  names(res) <- labs
  attr(res, "levels_ordered") <- c(lo_lab, hi_lab)
  res
}


# collapse_two_agegrps ----
# Relabel `col` IN PLACE into the two coarse groups. A no-op returning FALSE
# when `dt` is NULL or has no such column, so callers need no guard of their
# own -- `dis_characteristics` carries no `agegrp` at all, and the xps family
# uses `agegrp20` (20-year bands) which this deliberately does not touch.
#
# Idempotent: re-running it on already-collapsed labels reproduces them, which
# matters because export_main_tables() applies it per METRIC to copies of one
# shared base table.
collapse_two_agegrps <- function(dt, split_age = 65L, col = "agegrp") {
  if (is.null(dt) || !col %chin% names(dt)) return(invisible(FALSE))
  v <- dt[[col]]
  m <- two_agegrp_map(if (is.factor(v)) levels(v) else v, split_age, col)
  new <- unname(m[as.character(v)])
  if (is.factor(v)) {
    new <- factor(new, levels = intersect(attr(m, "levels_ordered"), new))
  }
  data.table::set(dt, j = col, value = new)
  invisible(TRUE)
}


# .promote_integer_measures ----
# Promote the INTEGER measure columns of an aggregated table to double, leaving
# the key columns alone.
#
# Not cosmetic, and the hazard is narrower than it looks. Every consumer below
# divides one summed column by another and writes the quotient back, but the
# form of the assignment decides whether that is safe:
#
#   d[, value := value / denom]        replaces the WHOLE column ("plonk"), so
#                                      it is promoted to double. Harmless.
#   d[cases, on = ..., value := ...]   a JOIN-UPDATE writes into a SUBSET of
#                                      the existing column and must coerce to
#                                      that column's type -- an integer column
#                                      TRUNCATES, turning a case-fatality rate
#                                      of 0.034 into 0.
#
# Every division site in tbl_smmrs_core() and export_all_cause_mrtl_tables() is
# the second form. data.table warns, but a whole table of zeros is
# publishable-looking output and the warning is one of thousands in a run. The
# stock summaries store these counts as double so nothing truncates today --
# but that is a property of whoever wrote the parquet, not a guarantee of this
# code.
.promote_integer_measures <- function(dt, keys) {
  meas <- setdiff(names(dt), keys)
  if (length(meas) == 0L) return(invisible(FALSE))
  int <- meas[vapply(meas, function(n) is.integer(dt[[n]]), logical(1))]
  if (length(int) == 0L) return(invisible(FALSE))
  dt[, (int) := lapply(.SD, as.numeric), .SDcols = int]
  invisible(TRUE)
}


# validate_two_agegrps_split_age ----
# STRUCTURAL, data-free check, run in the PARENT next to
# validate_comparator_scenario(): a bad argument must not surface hours later
# from inside a forked worker. The band-alignment check needs the data and so
# lives in two_agegrp_map(), with export_tables() probing the summaries up front
# to raise it in the parent too.
validate_two_agegrps_split_age <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || x <= 0) {
    stop("`two_agegrps_split_age` must be a single positive whole number (the ",
         "first year of the upper age group), not ",
         paste(deparse(x), collapse = " "), ".", call. = FALSE)
  }
  invisible(as.integer(x))
}


# Columns the equity tables manage themselves, so they cannot also be named as
# an output stratum: `mc` is the Monte-Carlo axis the indices are quantiled
# over, and `scenario` is always present because every index is a contrast
# against *its* comparator.
.equity_reserved_vars <- c("mc", "scenario")


# =============================================================================
# Per-scenario comparator map
# -----------------------------------------------------------------------------
# `comparator_scenario` is either
#   * a bare scenario name -- the historical scalar form, in which every OTHER
#     scenario is contrasted against it; or
#   * a NAMED character vector mapping intervention -> comparator, with at most
#     one UNNAMED element acting as the default for every scenario not named.
#
# There is always EXACTLY ONE comparator per scenario. That is the load-bearing
# invariant: it makes `comparator` a FUNCTION of `scenario`, which is why it can
# be attached as a cosmetic column *after* quantiling rather than entering the
# key vector `x` -- and therefore never touching a group-by, one of the many
# positional setnames() calls, or a filename. It also keeps `scenario` a unique
# row key, so downstream code keyed on (scenario, year, type) is unaffected.
#
# Chains (prs -> qrisk -> sc0) are legal; cycles are not.
# =============================================================================

# Columns the comparator machinery owns; a summary must not already carry them.
.comparator_reserved_cols <- c("comparator", ".cmp__")


# .comparator_is_map ----
# TRUE iff the caller asked for a map. Deliberately SYNTACTIC -- never derived
# from the data -- so the output schema cannot vary between summary families,
# between localities or between reruns, and a bare "sc0" stays byte-identical
# to every table published to date.
.comparator_is_map <- function(cs) {
  !is.null(names(cs)) && any(nzchar(names(cs)), na.rm = TRUE)
}


# .resolve_comparators ----
# Concretise the argument against the scenarios a given summary ACTUALLY
# carries (families can legitimately differ). Returns a named character vector,
# names = intervention, values = comparator, with an attribute "dropped"
# listing requested pairs this summary cannot build.
#
# The scalar case is the identity: .resolve_comparators("sc0", scns) is
# setNames(rep("sc0", n), setdiff(scns, "sc0")) -- precisely the set the old
# `scenario != comparator_scenario` predicate selected.
.resolve_comparators <- function(comparator_scenario, scenarios) {
  cs <- comparator_scenario
  scenarios <- as.character(scenarios)
  nm <- names(cs)
  if (is.null(nm)) nm <- rep("", length(cs))
  nm[is.na(nm)] <- ""
  default  <- unname(cs[!nzchar(nm)])      # length 0 or 1 (validated up front)
  explicit <- cs[nzchar(nm)]

  # Requested pairs whose INTERVENTION this summary does not carry.
  dropped <- character(0)
  gone_i <- setdiff(names(explicit), scenarios)
  if (length(gone_i)) {
    dropped <- paste0(gone_i, " -> ", unname(explicit[gone_i]))
  }

  keys <- intersect(scenarios, names(explicit))   # keeps `scenarios` order
  out  <- stats::setNames(unname(explicit[keys]), keys)
  if (length(default)) {
    rest <- setdiff(scenarios, c(names(out), default))
    out  <- c(out, stats::setNames(rep(default, length(rest)), rest))
  }

  # ... and pairs whose COMPARATOR it does not carry. Dropping them is what
  # stops the cpp/cypp/dpp join leaving a raw LEVEL in a column published as
  # "cases prevented" (see the anti-join guard at that site).
  keep <- out %chin% scenarios
  if (any(!keep)) {
    dropped <- c(dropped, paste0(names(out)[!keep], " -> ", unname(out)[!keep]))
  }
  out <- out[keep]
  attr(out, "dropped") <- dropped          # set AFTER subsetting: `[` drops it
  out
}


# .comparator_cycles ----
# The map is a functional graph (out-degree 1), so a cycle is found by walking
# forward from each node. Chains terminate at a scenario that is not itself an
# intervention (sc0), which is the normal shape.
.comparator_cycles <- function(m) {
  bad <- character(0)
  for (s in names(m)) {
    seen <- s
    cur <- s
    repeat {
      nxt <- unname(m[cur])
      if (is.na(nxt)) break
      if (nxt %chin% seen) { bad <- c(bad, s); break }
      seen <- c(seen, nxt)
      cur <- nxt
      if (!cur %chin% names(m)) break
    }
  }
  unique(bad)
}


# validate_comparator_scenario ----
# STRUCTURAL, data-free rules. Called from export_tables() next to
# build_strata_config(), i.e. in the PARENT process before any task is
# dispatched -- the same slot and the same reason as validate_equity_strata().
validate_comparator_scenario <- function(cs, arg = "comparator_scenario") {
  if (!is.character(cs) || length(cs) == 0L) {
    stop("`", arg, "` must be a character vector of length >= 1; got ",
         class(cs)[1L], " of length ", length(cs),
         ". Use a bare scenario name, or a named vector mapping intervention ",
         "-> comparator, e.g. c(\"sc0\", sc2 = \"sc1\").", call. = FALSE)
  }
  if (anyNA(cs) || !all(nzchar(cs))) {
    stop("`", arg, "` must not contain NA or empty scenario names.",
         call. = FALSE)
  }
  nm <- names(cs)
  if (is.null(nm)) nm <- rep("", length(cs))
  if (anyNA(nm)) {
    stop("`", arg, "` has an NA name. An NA name is a typo, not an unnamed ",
         "default element.", call. = FALSE)
  }
  n_default <- sum(!nzchar(nm))
  if (n_default > 1L) {
    stop("`", arg, "` may have AT MOST ONE unnamed element (the default ",
         "comparator); got ", n_default, ": ",
         paste(sQuote(unname(cs[!nzchar(nm)])), collapse = ", "),
         ". Name each intervention explicitly, e.g. c(sc1 = \"sc0\"). Note a ",
         "two-element unnamed vector is NOT read positionally as ",
         "(intervention, comparator): position would decide the SIGN of every ",
         "reported benefit.", call. = FALSE)
  }
  named <- nm[nzchar(nm)]
  if (anyDuplicated(named)) {
    dup <- unique(named[duplicated(named)])
    stop("`", arg, "` names intervention scenario(s) ",
         paste(sQuote(dup), collapse = ", "),
         " more than once (comparators: ",
         paste(sQuote(unname(cs[nm %chin% dup])), collapse = ", "),
         "). Each scenario has exactly one comparator.", call. = FALSE)
  }
  self <- nzchar(nm) & nm == unname(cs)
  if (any(self)) {
    stop("`", arg, "`: a scenario cannot be its own comparator (",
         paste(sQuote(nm[self]), collapse = ", "),
         "). The contrast would be identically zero and would be published as ",
         "a table of zeros rather than failing.", call. = FALSE)
  }
  invisible(TRUE)
}


# validate_equity_strata ----
# Each entry of the equity strata list is one output table. Any variable carried
# by the summaries (`strata_for_output` in the sim_design YAML) can be an output
# stratum -- the gradient is then fit *within* each stratum -- except that
# `dimd`/`qimd` are not stratifications at all: they select the deprivation axis
# the gradient is fit over. Whether a stratum variable actually exists can only
# be judged against the summaries, so that check lives in export_equity_tables();
# this validates the structural rules that hold regardless of the data.
validate_equity_strata <- function(strata) {
  if (!is.list(strata)) {
    stop("`strata$equity` must be a list of character vectors.", call. = FALSE)
  }
  for (i in seq_along(strata)) {
    s <- strata[[i]]
    lbl <- sprintf("`strata$equity[[%d]]`", i)
    if (!is.character(s) || length(s) == 0L) {
      stop(lbl, " must be a non-empty character vector.", call. = FALSE)
    }
    if (anyDuplicated(s)) {
      stop(lbl, " has duplicated entries: ",
           paste(unique(s[duplicated(s)]), collapse = ", "), ".", call. = FALSE)
    }
    bad <- intersect(s, .equity_reserved_vars)
    if (length(bad) > 0L) {
      stop(lbl, " names reserved column(s): ", paste(bad, collapse = ", "),
           ". `mc` is the Monte-Carlo axis the indices are quantiled over and",
           " `scenario` is always included, so neither can be an output",
           " stratum.", call. = FALSE)
    }
    if (!"year" %in% s) {
      stop(lbl, " must include \"year\". The benefit is accumulated over years",
           " and each output row reports the benefit cumulated to that year,",
           " against that year's reference population -- so `year` defines the",
           " rows rather than merely subsetting them.", call. = FALSE)
    }
  }
  invisible(TRUE)
}


# build_equity_plans ----
# Expands an equity strata list into one plan per output table. `dimd`/`qimd`
# in an entry select the gradient axis; an entry naming both yields one table
# per axis. An entry naming neither falls back to `dimd`, which keeps the
# historical defaults -- and their filenames -- unchanged.
build_equity_plans <- function(strata) {
  unlist(lapply(strata, function(s) {
    out_vars <- setdiff(s, .deprivation_vars)   # keeps the caller's order
    grads <- intersect(.deprivation_vars, s)    # canonical dimd, qimd order
    implicit <- length(grads) == 0L
    if (implicit) grads <- "dimd"
    lapply(grads, function(g) list(
      out_vars = out_vars,
      gradient = g,
      # The gradient token is echoed in the filename only when the caller wrote
      # it, so `list("year", c("year", "sex"))` still produces exactly the
      # filenames it always has. `agegrp` is spelled `agegroup` in the suffix,
      # matching every other table (see make_strata_configs()).
      suffix = gsub(
        "agegrp", "agegroup",
        paste(if (implicit) out_vars else c(out_vars, g), collapse = "-")
      )
    ))
  }), recursive = FALSE)
}


# export_tables ----
# Exports summary tables from simulation summaries.
# See main class documentation in Simulation_class.R for details.
#
# @param strata A named list specifying stratification levels for different table types.
#   If NULL (default), uses built-in defaults matching process_out_for_NotinghamLA.R.
#   The list can contain:
#   - ons: List of character vectors for non-standardised main tables
#   - esp: List of character vectors for standardised (ESP) main tables
#   - mrtl_ons: List of character vectors for non-standardised all-cause mortality tables
#   - mrtl_esp: List of character vectors for standardised all-cause mortality tables
#   - disease_char: List of character vectors for disease characteristics tables
#   - xps_ons: List of character vectors for non-standardised exposure tables
#   - xps_esp: List of character vectors for standardised exposure tables
#   - equity: List of character vectors for equity slope-index tables. `dimd`
#     and `qimd` mean something different here -- see below -- and the list is
#     validated.
#
#   Valid stratification variables:
#   - year: Simulation year (always required)
#   - sex: Sex (men/women)
#   - agegrp: Age groups (5-year bands)
#   - dimd: IMD deciles (10 levels: "1 most deprived" to "10 least deprived")
#   - agegrp20: 20-year age groups (xps tables only)
#   - qimd: IMD quintiles (5 levels)
#
#   `qimd` works in EVERY table family and does NOT need to be in
#   `strata_for_output`: whenever the summaries carry `dimd` it is derived by
#   the standard decile-pair collapse (see add_qimd_from_dimd()), which is
#   exact. The converse does not hold -- `dimd` must really be in the
#   summaries, and is unavailable in the xps tables, which are written
#   quintile-based by export_xps().
#
#   For `equity`, any variable the summaries carry can be an output stratum
#   (the gradient is then fit within it), and `year` is required. The
#   exception is `dimd`/`qimd`, which do NOT stratify -- they select the
#   deprivation gradient the index is fit over. Naming both in one entry writes
#   one table per gradient; naming neither falls back to `dimd`.
#
# @examples
# # Use default strata (includes dimd and qimd stratification)
# sim$export_tables()
#
# # Custom strata - minimal outputs
# sim$export_tables(strata = list(
#   ons = list("year", c("year", "sex")),
#   esp = list("year")
# ))
#
# # Full stratification with dimd
# sim$export_tables(strata = list(
#   ons = list("year", c("year", "sex"), c("year", "dimd"),
#              c("year", "agegrp"), c("year", "agegrp", "sex"),
#              c("year", "agegrp", "sex", "dimd")),
#   esp = list("year", c("year", "sex"), c("year", "dimd"),
#              c("year", "sex", "dimd"))
# ))
#
# # Every table family by IMD quintiles. `qimd` need not be in
# # strata_for_output -- it is derived from dimd whenever the summaries carry it
# sim$export_tables(strata = list(
#   ons          = list("year", c("year", "qimd"), c("year", "agegrp", "sex", "qimd")),
#   esp          = list("year", c("year", "qimd")),
#   mrtl_ons     = list("year", c("year", "qimd")),
#   mrtl_esp     = list("year", c("year", "qimd")),
#   disease_char = list("year", c("year", "qimd")),
#   equity       = list(c("year", "qimd"))
# ))
#
# # Equity indices over IMD quintiles instead of deciles, plus both gradients
# # by year, and a gradient fit within each age group
# sim$export_tables(strata = list(
#   equity = list(c("year", "dimd", "qimd"), c("year", "sex", "qimd"),
#                 c("year", "agegrp"))
# ))
# # -> equity cpp slope index by year-dimd (not standardised).csv
# #    equity cpp slope index by year-qimd (not standardised).csv
# #    equity cpp slope index by year-sex-qimd (not standardised).csv
# #    equity cpp slope index by year-agegroup (not standardised).csv
#' @description
#' Export summary tables for policy analysis.
#'
#' Builds main, all-cause mortality, disease-characteristics, exposure,
#' (optionally) cost-effectiveness and (optionally) equity slope-index tables
#' from the per-summary files produced by `$export_summaries()`. Outputs are
#' written to `output_dir/tables/` (or `output_dir/tables2agegrps/` when
#' `two_agegrps = TRUE`).
#'
#' When `equity = TRUE`, equity slope-index tables are written as
#' `equity <metric> slope index by <strata> (not standardised).csv` for each of
#' CPP, CYPP, DPP and net QALYs. These summarise how each cumulative benefit is
#' distributed across deprivation groups -- DIMD deciles by default, or IMD
#' quintiles (`qimd`) if asked for via `strata$equity` -- using absolute and
#' relative analogues of the Slope Index of Inequality and Relative Index of
#' Inequality (see the `equity` argument).
#'
#' When `cea = TRUE`, cost-effectiveness tables are written as
#' `cost-effectiveness by <strata> (<perspective>-<scale>) (not standardised).csv`,
#' containing the cumulative discounted incremental QALYs (`dQALYs_cuml`) and
#' costs (`dCosts_cuml`), the incremental cost-effectiveness ratio (`ICER`), and
#' the net monetary benefit (`NMB_at_wtp_<threshold>`) at each willingness-to-pay
#' threshold, for each QALY scale (EQ5D5L) and from three cost perspectives:
#' \itemize{
#'   \item \emph{societal} -- `total_cost`, which already nets economic output
#'     out of healthcare + social care + informal care;
#'   \item \emph{healthcare} -- `healthcare_cost`, direct treatment costs only;
#'   \item \emph{healthcare_socialcare} -- `healthcare_cost` +
#'     `socialcare_cost`. This is NICE's reference case, "NHS and Personal
#'     Social Services": formal social care is included, informal care and lost
#'     productivity are not.
#' }
#' User-defined `*_costs` columns are always added to the societal perspective
#' and, if named in `custom_costs_in_healthcare`, to the healthcare and
#' healthcare_socialcare perspectives as well. A perspective whose required
#' built-in columns are absent from the costs summary is skipped rather than
#' silently falling back to a narrower one.
#'
#' @param baseline_year_for_change_outputs Integer. Reference year used for
#'   computing change-from-baseline columns and the start year for the
#'   cumulative cost-effectiveness calculations. Two-digit values (e.g. `19`)
#'   are auto-promoted to four-digit (`2019`).
#' @param prbl Numeric vector of probability levels for output quantiles
#'   (median plus uncertainty bounds). Default
#'   `c(0.5, 0.025, 0.975, 0.1, 0.9)`.
#' @param comparator_scenario Character. The scenario(s) each intervention is
#'   contrasted against. Two forms are accepted.
#'
#'   **A bare scenario name** (the default, `"sc0"`) is the historical
#'   behaviour: every *other* scenario is contrasted against it. Output is
#'   unchanged in every respect -- no extra column, no filename change.
#'
#'   **A named character vector** maps intervention -> comparator, with at most
#'   one *unnamed* element acting as the default for every scenario not named:
#'
#'   ```r
#'   export_tables(comparator_scenario = c(
#'     "sc0",                                  # default: other arms vs sc0
#'     mult_base_prs    = "mult_base_qrisk",   # the within-variant increment
#'     mult_lowzero_prs = "mult_lowzero_qrisk"
#'   ))
#'   ```
#'
#'   A scenario may be both an intervention and somebody else's comparator --
#'   above, `mult_base_qrisk` is contrasted against `sc0` *and* serves as
#'   `mult_base_prs`'s reference. Chains are legal; cycles are rejected.
#'   Omitting the unnamed element builds only the pairs you list, and the
#'   scenarios left out are named in a `warning()`.
#'
#'   There is exactly **one comparator per scenario**, so `scenario` remains a
#'   unique row key and no row multiplication occurs. In map mode the contrast
#'   tables (`cpp`, `cypp`, `dpp`, `net_qalys`, `net_costs`, the
#'   cost-effectiveness tables and the equity slope-index tables) carry an extra
#'   `comparator` column immediately after `scenario`, and a `comparators.csv`
#'   manifest of the resolved map is written to the tables directory. Sign
#'   conventions are unchanged: `cpp`/`cypp`/`dpp` are comparator minus
#'   intervention; `net_qalys`/`net_costs` and the CEA increments are
#'   intervention minus comparator.
#'
#'   Note this cannot be reproduced after the fact by post-processing tables
#'   built against a single comparator: `(a - c) - (b - c) = a - b` holds per
#'   Monte-Carlo draw, but the published tables are already quantiled across
#'   draws, and neither a median nor an ICER commutes with differencing.
#'
#'   Scenario names are validated against the summaries up front, so a typo
#'   fails immediately rather than hours into the export.
#' @param two_agegrps Logical. If `TRUE`, `agegrp` is reported as TWO coarse
#'   groups instead of the design's 5-year bands, and the tables are written to
#'   `tables2agegrps/` instead of `tables/`. The split age is
#'   `two_agegrps_split_age`.
#'
#'   It is a pure **relabelling** applied before each family's own group-by, so
#'   for counts and for rates built as summed-numerator / summed-denominator it
#'   is exactly equivalent to having banded at summary time. The one exception
#'   is `equity`, which fits a *model*: the index is invariant only under an
#'   exactly linear gradient, so a two-group and a 5-year-band fit are different
#'   estimands (see `export_equity_tables()`).
#'
#'   Applies to every family whose summaries carry `agegrp`: the main tables,
#'   all-cause mortality by disease, the CEA tables and the equity tables.
#'   `disease characteristics` has no `agegrp` column, and the exposure (`xps`)
#'   tables use `agegrp20` (20-year bands, plus `90+` and `All`) which is a
#'   different variable and is deliberately left untouched.
#' @param two_agegrps_split_age Integer, default `65L`. The first year of the
#'   UPPER age group when `two_agegrps = TRUE`, i.e. the default reproduces the
#'   historical `30-64` / `65-99` split. Ignored when `two_agegrps = FALSE`.
#'
#'   It must fall *between* two age bands rather than inside one: a split of
#'   `62` would have to cut `60-64` in half and is refused. It need not start a
#'   band that is actually present -- a locality whose data happen to have an
#'   empty band at the split is still splittable there. A split that would leave
#'   one of the two groups empty is also refused, since that is not two groups.
#'   Both the structural check and the band-alignment check run in the parent
#'   process before any task is dispatched, so a typo fails in seconds rather
#'   than from inside a forked worker minutes later.
#'
#'   The group labels are derived from the data: `30-64` / `65-99` for the stock
#'   England design, but a run reporting ages 20-89 gets `20-64` / `65-89`, and
#'   an open-ended top band gives `65+`.
#' @param strata Optional named list overriding the default stratification
#'   configuration. See examples in the file header for shape; passed
#'   through `private$build_strata_config()`.
#'
#'   Any element may be stratified by `"qimd"` (IMD quintiles) instead of
#'   `"dimd"` (deciles), and `qimd` does **not** need to be listed in
#'   `strata_for_output` in the sim_design YAML: whenever the loaded summaries
#'   carry `dimd`, a `qimd` column is derived from it by the standard
#'   decile-pair collapse (deciles 1-2 -> quintile 1, and so on). That is exact
#'   rather than approximate, because every aggregation here groups first and
#'   then sums, so summing two deciles reproduces precisely what aggregating by
#'   quintile at summary time would have given -- and rates, which are formed as
#'   summed-numerator / summed-denominator after the group-by, follow
#'   automatically. The converse is impossible: deciles cannot be recovered from
#'   quintiles, so a `dimd` stratum requires `dimd` in the summaries, and is
#'   unavailable in the `xps_*` elements because `$export_xps()` writes those
#'   datasets quintile-based already.
#'
#'   In the `equity` element only, `"dimd"`/`"qimd"` select the deprivation
#'   gradient rather than a stratification, and `"year"` is required -- see the
#'   `equity` argument.
#' @param multicore Logical. If `TRUE`, runs table-building tasks in
#'   parallel with single-threaded workers; otherwise runs sequentially.
#' @param cea Logical. If `TRUE` (default), also build the cost-effectiveness
#'   (ICER / NMB) tables from the `qalys` and `costs` summaries.
#'   On a run with no intervention arm (only the comparator scenario present)
#'   there is nothing to contrast, so no CEA table is written -- but the worker
#'   still READS `qalys` and `costs` before it can discover that, since the
#'   scenario list is a property of the data. Pass `cea = FALSE` explicitly on a
#'   baseline-only run to skip the read and free one of the parallel task slots.
#' @param wtp Numeric vector of willingness-to-pay thresholds (currency per
#'   QALY) at which the net monetary benefit is computed. Default
#'   `c(25000, 35000)` (GBP/QALY). NICE's standard reference range is
#'   £20,000-£30,000 per QALY; pass `wtp = c(20000, 30000)` for that. Note the
#'   thresholds appear in the output **column names**
#'   (`NMB_at_wtp_25000`, ...), so changing them changes the CEA table schema.
#'
#'   `wtp` and the discount levels are **crossed**: every threshold is evaluated
#'   at every discount level, so a CEA table holds `length(wtp)` x (number of
#'   discount levels) NMB rows per stratum cell. (Contrast the two discount-rate
#'   vectors, which are *paired* with each other - see `qaly_discount_rate`.)
#'   Filter on both `discount` and `type` to isolate one series.
#'
#'   A threshold is a ratio of present values (currency per QALY), so it is
#'   *not* itself discounted: it multiplies the already-discounted
#'   `dQALYs_cuml`, making `NMB = wtp * dQALYs_cuml - dCosts_cuml` a present
#'   value at that level's rates. `ICER` does not depend on `wtp` at all
#'   (it is `dCosts_cuml / dQALYs_cuml`), so it appears once per discount level
#'   rather than once per threshold, and its rows are absent wherever
#'   `dQALYs_cuml` is 0 and the ratio is undefined.
#' @param qaly_discount_rate,cost_discount_rate Numeric *vectors*. Annual
#'   discount rates (percent) applied to QALYs and costs. Default `c(0, 3.5)`
#'   each, i.e. undiscounted and UK NICE's 3.5%. These arguments are the single
#'   source of truth for discounting; there is no `discounting` block in
#'   `sim_design.yaml`.
#'
#'   The two vectors are paired **element-wise** into discount *levels* (a
#'   length-1 rate is recycled against a longer one), so the defaults give two
#'   levels: undiscounted, and 3.5% on QALYs and costs alike. Differential rates
#'   are set by passing equal-length vectors, e.g.
#'   `qaly_discount_rate = c(0, 1.5)` with `cost_discount_rate = c(0, 3.5)`.
#'
#'   Every level appears in the `qalys`, `net_qalys`, `costs`, `net_costs` **and**
#'   cost-effectiveness tables - not just the cost-effectiveness ones. The file
#'   set is unchanged: each of those tables gains a `discount` column and carries
#'   one block of rows per level, labelled `"0%"`, `"3.5%"`, or
#'   `"QALYs 1.5%/costs 3.5%"` when a level's two rates differ. The `0%` level
#'   reproduces the undiscounted figures, so nothing is lost by adding
#'   discounted ones. Filter on `discount` to pick a level.
#'
#'   The prevalence, incidence, mortality and exposure tables carry no
#'   `discount` column: they are rates and counts, not monetary or QALY flows,
#'   so a present value is not defined for them. The equity slope-index tables
#'   are likewise not expanded - their metric set mixes QALYs with plain event
#'   counts (CPP / CYPP / DPP), which a single `discount` column could not
#'   describe.
#' @param discount_from_year Integer. First year from which present values are
#'   discounted (`PV = FV / (1 + rate/100)^max(0, year - base)`), applied to each
#'   year's flow *before* the cumulative columns are accumulated, so the `_cuml`
#'   columns are sums of present values. When `NULL` (default) it defaults to
#'   `baseline_year_for_change_outputs`. Years at or before it are undiscounted.
#'   Accepts the same two-digit shorthand as `baseline_year_for_change_outputs`:
#'   a value `<= 100` is promoted by 2000, so `19` means 2019.
#' @param custom_costs_in_healthcare Character vector of user-defined `*_costs`
#'   column names to add to the healthcare **and** healthcare_socialcare
#'   perspectives (in addition to the built-in cost columns those perspectives
#'   already carry). `NULL`/`FALSE` (default) adds none; `TRUE` adds all
#'   user-defined cost columns. User-defined cost columns are always included in
#'   the societal perspective regardless of this argument. The two narrower
#'   perspectives share this setting deliberately: a programme cost belongs to
#'   the payer in both, so letting them diverge would make the pair
#'   incomparable.
#' @param equity Logical. If `TRUE` (default), also build the equity
#'   slope-index tables (absolute and relative analogues of the Slope Index of
#'   Inequality / Relative Index of Inequality) for the cumulative CPP, CYPP,
#'   DPP and net-QALYs benefits, distributed across deprivation groups.
#'   Every index is a benefit-vs-comparator contrast, so as with `cea` a
#'   baseline-only run writes no equity table -- yet still reads the `incd`,
#'   `prvl`, `mrtl` and `qalys` summaries in turn before finding that out. Pass
#'   `equity = FALSE` explicitly when there is no intervention arm.
#'   Written as `equity <metric> slope index by <strata> (not standardised).csv`
#'   with a `type` column giving:
#'   * `AEI_total` -- the SII scaled to the whole reference population
#'     (`SII * n`): the *hypothetical* total benefit difference if every person
#'     received the fitted most-deprived per-capita benefit instead of the fitted
#'     least-deprived one. It is **not** the benefit difference between the
#'     most- and least-deprived deciles themselves, which is roughly an order of
#'     magnitude smaller. Because it scales with `n` it is not comparable across
#'     localities, strata or years -- use `AEI_per100k` or `REI_rel` for that.
#'   * `AEI_per100k` -- the SII per 100,000 people *alive in the reporting year*.
#'     The regressand is a benefit cumulated from `baseline_year`, so this is a
#'     stock per head, not a rate per 100,000 person-years.
#'   * `total_benefit` -- the total undiscounted benefit the gradient was fit to.
#'     Reported because the absolute indices are translation-invariant and so
#'     cannot distinguish levelling up from levelling down.
#'   * `REI_rel` -- the SII divided by the population-weighted mean benefit. `NA`
#'     unless that mean is positive, and unstable when it is near zero.
#'   * `RII_ratio` -- the ratio of fitted benefit at the most- to the
#'     least-deprived end. `NA` unless the fitted benefit is positive at *both*
#'     ends, which is common rather than exceptional.
#'   * `fit_R2` -- the population-weighted R-squared of the straight-line fit.
#'
#'   An `n_mc` column records how many Monte-Carlo draws each published quantile
#'   rests on; the finite-value filter is applied per `type`, so the relative
#'   indices can be summarised over a strict subset of the draws the absolute
#'   ones use. By the pro-poor sign convention a positive index means the most
#'   deprived gain more (inequality-reducing). The gradient is fit per
#'   Monte-Carlo iteration and quantiled across iterations. The grouped-data
#'   weighted least-squares estimator is that of Kakwani, Wagstaff & van
#'   Doorslaer (1997); the ridit/structured-regression setup follows
#'   Moreno-Betancur et al. (2015); and applying the index to a cumulative
#'   modelled benefit follows Cookson et al. (2026). Note `RII_ratio` is the
#'   Kunst-Mackenbach ratio of fitted extremes, which Moreno-Betancur et al.
#'   call a "previous definition" and criticise -- their own RII is `exp(beta)`
#'   from a log-linear model, which cannot be fit to a signed quantity such as
#'   net QALYs -- so they are cited here for the framework, not for that ratio.
#'   `REI_rel` is a Pamuk/Wagstaff-type relative SII, not an RII of any kind.
#'
#'   A gradient requires at least two deprivation groups **present in the
#'   summaries**; where only one is (e.g. a single-decile locality) the index is
#'   undefined and those rows are omitted, reported when `logs` are enabled. An
#'   intervention *targeted* at one decile is a different case: the untouched
#'   deciles contribute zero benefit, so all groups are present and a fully
#'   extrapolated index is produced. Read those with care and see
#'   `vignette("understanding_model_outputs")`.
#'
#'   **Output strata (`strata$equity`).** Each entry of the list is one output
#'   table, and any variable the summaries carry (`strata_for_output` in the
#'   sim_design YAML: `agegrp`, `sex`, `ethnicity`, `sha`, ...) can be an output
#'   stratum -- the gradient is then fit *within* each stratum, exactly as
#'   filtering the summaries to that stratum and fitting without it would. Every
#'   entry must name `"year"`, because the benefit is accumulated over years and
#'   each row reports the benefit cumulated to that year against that year's
#'   reference population, so `year` defines the rows rather than subsetting
#'   them. `"mc"` and `"scenario"` are reserved. Variables the summaries do not
#'   carry raise a `warning()` and that table is skipped.
#'
#'   **Choosing the gradient axis (`dimd` vs `qimd`).** The gradient is fit over
#'   DIMD deciles by default. Name `"dimd"` and/or `"qimd"` in an entry of
#'   `strata$equity` to choose explicitly; unlike every other variable these do
#'   *not* stratify, because the deprivation axis is consumed into the index. An
#'   entry naming both writes one table per axis, and an entry naming neither
#'   falls back to `dimd` so the built-in defaults keep producing exactly the
#'   files they always have. The chosen axis is echoed in a `gradient` column in
#'   every file, and (when the caller wrote it) in the filename suffix:
#'   `list(c("year", "sex", "qimd"))` writes
#'   `equity cpp slope index by year-sex-qimd (not standardised).csv`. As in
#'   every other table family, `agegrp` is spelled `agegroup` in the suffix.
#'   Quintiles give a coarser, lower-variance gradient with less power to detect
#'   departures from linearity; deciles are usually preferable when available.
#'   `qimd` does **not** need to be in `strata_for_output` -- if the summaries
#'   carry `dimd` it is derived by the standard decile-pair collapse. The derived
#'   quintile *counts* are exact (event counts and populations are additive, and
#'   the collapse preserves the population-weighted ridit midpoint exactly), but
#'   the resulting *index* is invariant only when the per-capita benefit is
#'   exactly linear in the ridit -- under curvature a weighted fit to five points
#'   is a different estimand from a fit to ten. Treat the decile and quintile
#'   tables as two analyses of the same data, not two views of one number. The
#'   converse does not hold: a `dimd` gradient requires `dimd` in the summaries. If the
#'   requested axis is unavailable, or the summaries carry no deprivation column
#'   at all, a `warning()` is raised (not a `logs`-gated message) naming the
#'   tables that were skipped.
#' @param equity_ridit_reference Character, one of `"comparator"` (default) or
#'   `"scenario"`. Chooses the population used to build the ridit (deprivation)
#'   ranks and the population weights for the equity slope-index regression --
#'   the reference against which each decile's share of the population is
#'   measured. `"comparator"` uses the comparator scenario's decile populations,
#'   so the ridit scores and regression weights are *identical* across the
#'   scenarios being compared; `"scenario"` uses each intervention scenario's
#'   own decile populations, reflecting the socioeconomic composition that
#'   scenario actually produces.
#'
#'   They are deliberately **not** pinned across years: `N` is the comparator's
#'   population *of the reporting year*, so the ridit scores and weights are
#'   re-derived each year (the ordinal decile ordering is of course fixed).
#'   Pinning `N` at `baseline_year` was considered and rejected -- it would
#'   divide a benefit generated by the year-t population by the year-0
#'   population, so a benefit that is exactly equal per head in year t would
#'   register a spurious non-zero index purely from differential population
#'   growth between deciles.
#'
#'   **Which to choose (Renard et al. 2019).** Renard et al. showed that the SII
#'   and RII can move purely because the socioeconomic *composition* of the
#'   population shifts, even when group-specific rates are unchanged -- so a
#'   changing denominator population can masquerade as a change in inequality.
#'   For the usual question here -- *how does an intervention redistribute
#'   benefits across a fixed deprivation structure?* -- keep the default
#'   `"comparator"`: a common, fixed reference makes between-scenario and
#'   between-year comparisons clean and attributes differences to the benefits,
#'   not to composition. Choose `"scenario"` only when the intervention itself is
#'   expected to change the deprivation composition (e.g. large differential
#'   survival shifting who is alive in each decile) and you specifically want the
#'   index to reflect each scenario's own realised population; interpret
#'   between-scenario differences with the Renard caveat in mind.
#'
#'   Formally, `"scenario"` violates the precondition Cookson et al. (2026) state
#'   for their impact SII -- that the comparator and intervention arms have the
#'   same group sizes -- so treat it as a **sensitivity analysis or diagnostic**
#'   rather than a co-equal estimator: if it moves the index materially, that is
#'   evidence the intervention is reshaping the population, not evidence about
#'   how the benefit was distributed. Under `"scenario"` each arm gets its own
#'   ridit axis and weights, so cross-scenario comparisons are then made on
#'   different axes. In this model DIMD deciles are structurally near-fixed --
#'   in the ONS projection input the decile shares of England aged 30-99 move at
#'   most 0.18 percentage points between 2019 and 2043 -- so the two options
#'   usually differ only slightly.
#' @details When `logs: yes` in the design YAML, console output from this
#'   method (including `foreach`'s `.verbose` bookkeeping) is appended to
#'   `<output_dir>/logs/console.txt` and restored on exit. See
#'   `Simulation$run()` for the rationale and the caveats.
#' @return The `Simulation` object, invisibly.
Simulation$set("public", "export_tables", function(
    baseline_year_for_change_outputs = 2019L,
    prbl = c(0.5, 0.025, 0.975, 0.1, 0.9),
    comparator_scenario = "sc0",
    two_agegrps = FALSE,
    two_agegrps_split_age = 65L,
    strata = NULL,
    multicore = TRUE,
    cea = TRUE,
    wtp = c(25000, 35000),
    qaly_discount_rate = c(0, 3.5),
    cost_discount_rate = c(0, 3.5),
    discount_from_year = NULL,
    custom_costs_in_healthcare = NULL,
    equity = TRUE,
    equity_ridit_reference = c("comparator", "scenario")
) {
  # Console output -> logs/console.txt while this phase runs (no-op unless
  # sim_prm$logs). See private$start_console_log() for the rationale.
  private$start_console_log("export_tables()")
  on.exit(private$stop_console_log(), add = TRUE)

  equity_ridit_reference <- match.arg(equity_ridit_reference)

  # Structural (data-free) rules for `comparator_scenario`, checked HERE in the
  # parent -- same slot and same reason as validate_equity_strata(): a bad
  # argument must not surface hours later from inside a parallel worker.
  validate_comparator_scenario(comparator_scenario)

  # Same slot, same reason, for the two-age-group split. The structural check is
  # free; the band-alignment one reads a single column and runs ONLY under
  # two_agegrps, so the default path is unaffected.
  two_agegrps_split_age <- validate_two_agegrps_split_age(two_agegrps_split_age)
  if (two_agegrps) private$check_two_agegrps_split_age(two_agegrps_split_age)

  # Ensure baseline year is in full format (e.g. 2019, not 19)
  # Data is converted to full year format in export_main_tables()
  baseline_year_for_change_outputs <-
    promote_short_year(baseline_year_for_change_outputs)

  # Discounting is controlled solely by the arguments below (there is
  # intentionally no `discounting` block in sim_design.yaml). It reaches the
  # QALY / cost / net-QALY / net-cost tables AND the cost-effectiveness tables:
  # each reports every requested discount level side by side, tagged in a
  # `discount` column. Prevalence / incidence / mortality / exposure tables are
  # rates and counts rather than flows, so they carry no discount column; the
  # equity tables are left alone because they mix QALYs with event counts.
  # `discount_from_year` takes the SAME two-digit shorthand as
  # `baseline_year_for_change_outputs` above, and must be promoted the same way.
  # Without this, `discount_from_year = 19` is treated as the year 19 AD: the
  # exponent becomes `year - 19` ~ 2000, so 1/1.035^2000 underflows to 0 and
  # every discounted figure is silently published as zero. There is no error and
  # no warning -- a table of zeros looks like a result -- which is exactly the
  # failure mode `.promote_integer_measures()` guards against elsewhere.
  discount_from_year <- if (is.null(discount_from_year)) {
    baseline_year_for_change_outputs
  } else {
    promote_short_year(discount_from_year)
  }

  # One row per discount level, shared by the main (qalys/costs) and the CEA
  # tables so both report the same set of levels.
  discount_levels <- build_discount_levels(qaly_discount_rate, cost_discount_rate)

  # Thread control for parallel execution
  if (multicore) {
    arrow::set_cpu_count(1L)
    data.table::setDTthreads(threads = 1L, restore_after_fork = NULL)
    fst::threads_fst(nr_of_threads = 1L, reset_after_fork = NULL)
  } else {
    arrow::set_cpu_count(self$design$sim_prm$clusternumber_export)
    data.table::setDTthreads(
      threads = self$design$sim_prm$clusternumber_export,
      restore_after_fork = NULL
    )
    fst::threads_fst(
      nr_of_threads = self$design$sim_prm$clusternumber_export,
      reset_after_fork = NULL
    )
  }

  # Build strata configuration (merge user-provided with defaults)
  strata_cfg <- private$build_strata_config(strata, two_agegrps)

  tables_subdir <- if (two_agegrps) "tables2agegrps" else "tables"
  tables_dir <- private$output_dir(tables_subdir)
  private$create_new_folder(tables_dir)

  # Data-dependent checks for a comparator MAP (unknown scenario names, cycles,
  # orphans, reserved-column collisions) plus the `comparators.csv` manifest.
  # Map mode only: a bare `comparator_scenario` reaches none of this, so no
  # call that works today can start failing and no new file appears in an
  # existing output tree.
  if (.comparator_is_map(comparator_scenario)) {
    private$check_comparator_map(comparator_scenario, tables_dir)
  } else {
    # A manifest left by an earlier map run would describe contrasts these
    # tables were NOT built with. Tables are rewritten in place, so it would
    # otherwise survive and mislead. Removing our own stale artifact is part of
    # writing this directory; nothing else is touched.
    stale <- file.path(tables_dir, "comparators.csv")
    if (file.exists(stale)) {
      unlink(stale)
      # Say so: with cea = FALSE / equity = FALSE this call does NOT rewrite the
      # CEA and equity tables, so the manifest it removes still described THOSE
      # correctly even though it no longer describes the main contrast tables
      # this call has just rewritten against a single comparator.
      if (self$design$sim_prm$logs) {
        message("Removed a stale comparators.csv left by an earlier ",
                "comparator-map run; the tables rewritten now use the single ",
                "comparator '", comparator_scenario[[1L]], "'.")
      }
    }
  }

  # Build task list for parallel execution
  tasks <- list(
    list(
      id = 1L,
      type = "main",
      prbl = prbl,
      baseline_year = baseline_year_for_change_outputs,
      output_dir = private$output_dir(),
      tables_dir = tables_dir,
      comparator_scenario = comparator_scenario,
      two_agegrps = two_agegrps,
      two_agegrps_split_age = two_agegrps_split_age,
      discount_levels = discount_levels,
      discount_from_year = discount_from_year,
      strata_ons = strata_cfg$ons,
      strata_esp = strata_cfg$esp
    ),
    list(
      id = 2L,
      type = "all_cause_mrtl",
      prbl = prbl,
      summaries_dir = private$output_dir("summaries"),
      tables_dir = tables_dir,
      two_agegrps = two_agegrps,
      two_agegrps_split_age = two_agegrps_split_age,
      strata_ons = strata_cfg$mrtl_ons,
      strata_esp = strata_cfg$mrtl_esp
    ),
    list(
      id = 3L,
      type = "disease_char",
      prbl = prbl,
      summaries_dir = private$output_dir("summaries"),
      tables_dir = tables_dir,
      strata = strata_cfg$disease_char
    ),
    list(
      id = 4L,
      type = "xps",
      prbl = prbl,
      output_dir = private$output_dir(),
      tables_dir = tables_dir,
      strata_ons = strata_cfg$xps_ons,
      strata_esp = strata_cfg$xps_esp
    )
  )

  # Cost-effectiveness (ICER / NMB) tables, built from qalys + costs summaries.
  if (cea) {
    tasks[[length(tasks) + 1L]] <- list(
      id = 5L,
      type = "cea",
      prbl = prbl,
      summaries_dir = private$output_dir("summaries"),
      tables_dir = tables_dir,
      comparator_scenario = comparator_scenario,
      baseline_year = baseline_year_for_change_outputs,
      wtp = wtp,
      discount_levels = discount_levels,
      discount_from_year = discount_from_year,
      custom_costs_in_healthcare = custom_costs_in_healthcare,
      two_agegrps = two_agegrps,
      two_agegrps_split_age = two_agegrps_split_age,
      strata = strata_cfg$ons
    )
  }

  # Equity slope-index (SII / RII analogue) tables, built from the prvl / incd /
  # mrtl / qalys summaries. The deprivation gradient is fit across DIMD deciles.
  if (equity) {
    tasks[[length(tasks) + 1L]] <- list(
      id = 6L,
      type = "equity",
      prbl = prbl,
      summaries_dir = private$output_dir("summaries"),
      tables_dir = tables_dir,
      comparator_scenario = comparator_scenario,
      baseline_year = baseline_year_for_change_outputs,
      ridit_reference = equity_ridit_reference,
      strata = strata_cfg$equity,
      two_agegrps = two_agegrps,
      two_agegrps_split_age = two_agegrps_split_age
    )
  }

  if (multicore) {
    if (self$design$sim_prm$logs) {
      private$time_mark("Start exporting tables (parallel)")
    }

    n_cores <- min(length(tasks), self$design$sim_prm$clusternumber_export)

    if (.Platform$OS.type == "windows") {
      cl <- parallelly::makeClusterPSOCK(
        n_cores,
        dryrun = FALSE,
        quiet = !self$design$sim_prm$logs,
        rscript_startup = quote(local({
          library(CKutils)
          library(IMPACTncdEngland)
          library(R6)
          library(data.table)
          library(scales)
        })),
        rscript_args = c("--no-init-file", "--no-site-file", "--no-environ"),
        setup_strategy = "parallel"
      )
      on.exit(parallel::stopCluster(cl), add = TRUE)

      parallel::parLapplyLB(
        cl = cl,
        X = tasks,
        fun = function(task) {
          private$export_tables_hlpr(task, implicit_parallelism = FALSE)
          NULL
        }
      )
    } else {
      # Linux/macOS: forking
      doParallel::registerDoParallel(n_cores)
      foreach::foreach(
        task = tasks,
        .inorder = FALSE,
        .packages = c("R6", "CKutils", "IMPACTncdEngland", "data.table", "scales"),
        .verbose = self$design$sim_prm$logs
      ) %dopar% {
        private$export_tables_hlpr(task, implicit_parallelism = FALSE)
        NULL
      }
    }

    if (self$design$sim_prm$logs) {
      private$time_mark("End exporting tables (parallel)")
    }
  } else {
    # Sequential execution
    lapply(tasks, function(task) {
      private$export_tables_hlpr(task, implicit_parallelism = TRUE)
    })
  }

  invisible(self)
})


# export_tables_hlpr ----
# Helper function for parallel table export. Dispatches to the appropriate
# export function based on task type.
Simulation$set("private", "export_tables_hlpr", function(task, implicit_parallelism) {
  # Thread control
  if (implicit_parallelism) {
    arrow::set_cpu_count(self$design$sim_prm$clusternumber_export)
    data.table::setDTthreads(
      threads = self$design$sim_prm$clusternumber_export,
      restore_after_fork = NULL
    )
    fst::threads_fst(
      nr_of_threads = self$design$sim_prm$clusternumber_export,
      reset_after_fork = NULL
    )
  } else {
    arrow::set_cpu_count(1L)
    data.table::setDTthreads(threads = 1L, restore_after_fork = NULL)
    fst::threads_fst(nr_of_threads = 1L, reset_after_fork = NULL)
  }

  # Dispatch to appropriate export function
  switch(task$type,
    "main" = private$export_main_tables(
      prbl = task$prbl,
      baseline_year = task$baseline_year,
      output_dir = task$output_dir,
      tables_dir = task$tables_dir,
      comparator_scenario = task$comparator_scenario,
      two_agegrps = task$two_agegrps,
      two_agegrps_split_age = task$two_agegrps_split_age,
      discount_levels = task$discount_levels,
      discount_from_year = task$discount_from_year,
      strata_ons = task$strata_ons,
      strata_esp = task$strata_esp
    ),
    "all_cause_mrtl" = private$export_all_cause_mrtl_tables(
      prbl = task$prbl,
      summaries_dir = task$summaries_dir,
      tables_dir = task$tables_dir,
      two_agegrps = task$two_agegrps,
      two_agegrps_split_age = task$two_agegrps_split_age,
      strata_ons = task$strata_ons,
      strata_esp = task$strata_esp
    ),
    "disease_char" = private$export_disease_characteristics_tables(
      prbl = task$prbl,
      summaries_dir = task$summaries_dir,
      tables_dir = task$tables_dir,
      strata = task$strata
    ),
    "xps" = private$export_xps_tables(
      prbl = task$prbl,
      output_dir = task$output_dir,
      tables_dir = task$tables_dir,
      strata_ons = task$strata_ons,
      strata_esp = task$strata_esp
    ),
    "cea" = private$export_cea_tables(
      prbl = task$prbl,
      summaries_dir = task$summaries_dir,
      tables_dir = task$tables_dir,
      comparator_scenario = task$comparator_scenario,
      baseline_year = task$baseline_year,
      wtp = task$wtp,
      discount_levels = task$discount_levels,
      discount_from_year = task$discount_from_year,
      custom_costs_in_healthcare = task$custom_costs_in_healthcare,
      two_agegrps = task$two_agegrps,
      two_agegrps_split_age = task$two_agegrps_split_age,
      strata = task$strata
    ),
    "equity" = private$export_equity_tables(
      prbl = task$prbl,
      summaries_dir = task$summaries_dir,
      tables_dir = task$tables_dir,
      comparator_scenario = task$comparator_scenario,
      baseline_year = task$baseline_year,
      ridit_reference = task$ridit_reference,
      strata = task$strata,
      two_agegrps = task$two_agegrps,
      two_agegrps_split_age = task$two_agegrps_split_age
    )
  )

  gc(verbose = FALSE)
  invisible(NULL)
})


# probe_summary_scenarios ----
# The scenario names and column names present in the summaries that carry a
# comparator contrast. A projected one-column arrow scan, materialised via
# as.data.frame(): this is a pure aggregate, never a join, so it is outside the
# arrow-collect join-corruption hazard documented in CLAUDE.md -- and
# materialising means the result is a plain character vector regardless.
# Costs ~1 s per summary, and runs in MAP MODE ONLY, so the scalar path (every
# table published to date) pays nothing.
Simulation$set("private", "probe_summary_scenarios", function() {
  srcs <- c("incd", "prvl", "mrtl", "qalys", "costs")
  scns <- character(0)
  cols <- character(0)
  for (s in srcs) {
    p <- private$output_dir(paste0("summaries/", s, "_scaled_up"))
    if (!dir.exists(p)) next
    ds <- tryCatch(arrow::open_dataset(p), error = function(e) NULL)
    if (is.null(ds)) next
    cols <- union(cols, names(ds))
    if (!"scenario" %chin% names(ds)) next
    tb <- arrow::Scanner$create(ds, projection = "scenario")$ToTable()
    scns <- union(scns, as.character(as.data.frame(tb)$scenario))
  }
  list(scenarios = sort(scns), columns = cols)
})


# check_two_agegrps_split_age ----
# Up-front, DATA-DEPENDENT validation of `two_agegrps_split_age`, run in the
# PARENT for the same reason as check_comparator_map(): a split age that does
# not start an age band is a typo, and the user should learn that in a second
# rather than from a forked worker several minutes in. TWO_AGEGRPS MODE ONLY,
# so the default path pays nothing.
#
# Reads one summary's `agegrp` column via a projected arrow scan (an aggregate,
# never a join -- see probe_summary_scenarios). Silently does nothing if no
# summary carries `agegrp`: there is then nothing for the collapse to act on,
# and collapse_two_agegrps() is a no-op everywhere.
Simulation$set("private", "check_two_agegrps_split_age", function(split_age) {
  for (s in c("prvl", "incd", "mrtl", "all_cause_mrtl_by_dis", "qalys", "costs")) {
    p <- private$output_dir(paste0("summaries/", s, "_scaled_up"))
    if (!dir.exists(p)) next
    ds <- tryCatch(arrow::open_dataset(p), error = function(e) NULL)
    if (is.null(ds) || !"agegrp" %chin% names(ds)) next
    tb <- arrow::Scanner$create(ds, projection = "agegrp")$ToTable()
    labs <- unique(as.character(as.data.frame(tb)$agegrp))
    two_agegrp_map(labs, split_age, "agegrp")   # errors if it splits no band
    return(invisible(TRUE))
  }
  invisible(FALSE)
})


# check_comparator_map ----
# Up-front, DATA-DEPENDENT validation of a comparator map. Hard errors are
# raised in the PARENT, before any task is dispatched, so nothing new can throw
# inside a forked/PSOCK worker where a message may not surface reliably.
# MAP MODE ONLY -- the scalar path keeps its historical warn-and-skip behaviour.
Simulation$set("private", "check_comparator_map", function(cs, tables_dir) {
  probe <- private$probe_summary_scenarios()
  avail <- probe$scenarios

  clash <- intersect(.comparator_reserved_cols, probe$columns)
  if (length(clash) > 0L) {
    stop("`comparator_scenario`: the summaries already carry column(s) ",
         paste(sQuote(clash), collapse = ", "),
         ", which the comparator map needs. Rename them in ",
         "`strata_for_output` and re-run $export_summaries().", call. = FALSE)
  }

  if (length(avail) == 0L) {
    warning("`comparator_scenario`: no summaries found to validate the map ",
            "against; the per-family checks remain the backstop.",
            call. = FALSE, immediate. = TRUE)
    return(invisible(NULL))
  }

  nm <- names(cs)
  if (is.null(nm)) nm <- rep("", length(cs))
  unknown <- setdiff(unique(c(nm[nzchar(nm)], unname(cs))), avail)
  if (length(unknown) > 0L) {
    # A NAMED vector is a map, so a name that is not a scenario is an error --
    # including the easy slip of naming a single element for readability
    # (c(baseline = "sc0")), which used to be ignored and is now significant.
    hint <- if (length(cs) == 1L && nzchar(nm[1L])) {
      paste0(" Note a NAMED element is read as <intervention> = <comparator>; ",
             "for the historical behaviour (every other scenario against it) ",
             "pass the bare name, i.e. \"", unname(cs[1L]), "\".")
    } else {
      ""
    }
    stop("`comparator_scenario`: scenario(s) ",
         paste(sQuote(unknown), collapse = ", "),
         " are not present in any summary. Available: ",
         paste(sQuote(avail), collapse = ", "), ".", hint, call. = FALSE)
  }

  cmap <- .resolve_comparators(cs, avail)

  cyc <- .comparator_cycles(cmap)
  if (length(cyc) > 0L) {
    stop("`comparator_scenario`: comparator cycle through ",
         paste(sQuote(cyc), collapse = ", "),
         ". Every contrast must terminate; chains such as ",
         "prs -> qrisk -> sc0 are fine.", call. = FALSE)
  }

  # The design's one real footgun: with no unnamed default, unlisted scenarios
  # silently get no contrast at all. Warn here, in the parent, naming them --
  # so it lands at the TOP of logs/console.txt rather than buried.
  orphans <- setdiff(setdiff(avail, names(cmap)), unname(cmap))
  if (length(orphans) > 0L) {
    warning("`comparator_scenario`: scenario(s) ",
            paste(sQuote(orphans), collapse = ", "),
            " have no comparator and will produce NO comparison rows. Add an ",
            "unnamed default element (e.g. c(\"sc0\", ...)) to contrast them ",
            "too.", call. = FALSE, immediate. = TRUE)
  }

  # Machine-readable record of what was requested, written BEFORE any task so
  # it exists even if a family later fails. Map mode only, so it cannot appear
  # in an existing published output tree.
  fwrite(data.table(scenario = names(cmap), comparator = unname(cmap)),
         file.path(tables_dir, "comparators.csv"))
  invisible(cmap)
})


# build_strata_config ----
# Builds the strata configuration by merging user-provided strata with defaults.
# Defaults match the stratification from process_out_for_NotinghamLA.R
Simulation$set("private", "build_strata_config", function(user_strata, two_agegrps = FALSE) {
  # Default strata configurations (from process_out_for_NotinghamLA.R)
  if (two_agegrps) {
    # For two_agegrps mode, only include agegrp-based strata
    defaults <- list(
      ons = list(
        c("year", "agegrp"),
        c("year", "agegrp", "sex"),
        c("year", "agegrp", "sex", "dimd")
      ),
      esp = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      mrtl_ons = list(
        c("year", "agegrp"),
        c("year", "agegrp", "sex"),
        c("year", "agegrp", "sex", "dimd")
      ),
      mrtl_esp = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      disease_char = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      xps_ons = list(
        "year",
        c("year", "agegrp20"),
        c("year", "sex"),
        c("year", "qimd"),
        c("year", "agegrp20", "sex"),
        c("year", "agegrp20", "sex", "qimd")
      ),
      xps_esp = list(
        "year",
        c("year", "sex"),
        c("year", "qimd"),
        c("year", "sex", "qimd")
      ),
      # Equity slope-index output strata. DIMD deciles are always the gradient
      # axis (consumed into the index), so it never appears here.
      equity = list(
        "year",
        c("year", "sex")
      )
    )
  } else {
    # Standard strata configurations
    defaults <- list(
      ons = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "agegrp"),
        c("year", "agegrp", "sex"),
        c("year", "agegrp", "sex", "dimd")
      ),
      esp = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      mrtl_ons = list(
        "year",
        c("year", "sex"),
        c("year", "agegrp"),
        c("year", "agegrp", "sex"),
        c("year", "agegrp", "sex", "dimd")
      ),
      mrtl_esp = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      disease_char = list(
        "year",
        c("year", "sex"),
        c("year", "dimd"),
        c("year", "sex", "dimd")
      ),
      xps_ons = list(
        "year",
        c("year", "agegrp20"),
        c("year", "sex"),
        c("year", "qimd"),
        c("year", "agegrp20", "sex"),
        c("year", "agegrp20", "sex", "qimd")
      ),
      xps_esp = list(
        "year",
        c("year", "sex"),
        c("year", "qimd"),
        c("year", "sex", "qimd")
      ),
      # Equity slope-index output strata. DIMD deciles are always the gradient
      # axis (consumed into the index), so it never appears here.
      equity = list(
        "year",
        c("year", "sex")
      )
    )
  }

  # Merge user-provided strata with defaults (user overrides defaults)
  result <- defaults
  if (!is.null(user_strata)) {
    for (name in names(user_strata)) {
      if (name %in% names(defaults)) {
        result[[name]] <- user_strata[[name]]
      }
    }
  }

  # Equity strata accept a much smaller vocabulary than the rest, and silently
  # ignoring an unsupported variable produced a mislabelled duplicate table.
  # Validate here so it fails fast, before any table-building work starts.
  validate_equity_strata(result$equity)

  return(result)
})


# .tbl_col_prefixes ----
# Prefix applied to each quantile column of an output table, keyed by metric
# (`what`). Lifted out of tbl_smmrs_core() so the naming convention is
# testable -- it is pure data, and the bug it guards against is invisible at
# the call site.
#
# Convention: the level metric takes `<family>_rate_` (or `_mean_`), and BOTH
# of its change variants take a prefix distinct from the level and from each
# other, always containing "change" so that a downstream
# `grep("change", names(dt))` is a safe way to select change columns. Every
# RELATIVE prefix is additionally unique across families.
#
# Two entries have been corrected, both of which published a column whose name
# asserted something the values did not support:
#
#  * `ftlt_change_relative` mapped to `ftlt_rate_` -- the same prefix as the
#    `ftlt` LEVEL table. It published a ratio (0.72 = 72% of the baseline year)
#    under a name claiming to be a rate, and, containing neither "change" nor
#    "rel", hid those six tables from any name-pattern selection.
#
#  * `prvl_change_relative` / `incd_change_relative` both mapped to
#    `prct_change_relative_`. The relative change is computed as
#    `value / i.value` (see the `_change_relative$` branch of tbl_smmrs_core),
#    so the quantity is a RATIO to the baseline year -- 1.0 at baseline, 0.9 for
#    a 10% reduction. It is not a percentage and never was: a reader taking
#    `prct_` at face value misreads a 10% fall as a 0.9% fall, an order-of-
#    magnitude error that looks entirely plausible. The shared prefix also left
#    the column unable to say which family it came from.
#
# `abs_change_` is deliberately NOT renamed alongside them: an absolute change
# is exactly what that name claims, so it is merely generic (shared by `prvl`
# and `incd`), not wrong. Renaming it would break consumers for no correctness
# gain. Hence the uniqueness assertion in the tests covers relative only.
# Tests: test_tbl_col_prefixes.R.
.tbl_col_prefixes <- function() {
  c(
    "prvl" = "prvl_rate_", "prvl_change_relative" = "prvl_change_relative_", "prvl_change_absolute" = "abs_change_",
    "incd" = "incd_rate_", "incd_change_relative" = "incd_change_relative_", "incd_change_absolute" = "abs_change_",
    "ftlt" = "ftlt_rate_", "ftlt_change_relative" = "ftlt_change_relative_", "ftlt_change_absolute" = "ftlt_abs_change_",
    "mrtl" = "mrtl_rate_", "mrtl_change_relative" = "mrtl_change_relative_", "mrtl_change_absolute" = "mrtl_abs_change_",
    "dis_mrtl" = "disease_mrtl_rate_", "dis_mrtl_change_relative" = "disease_mrtl_change_relative_", "dis_mrtl_change_absolute" = "disease_mrtl_abs_change_",
    "qalys" = "qalys_", "net_qalys" = "net_qalys_",
    "costs" = "costs_", "net_costs" = "net_costs_",
    "contd" = "contd_mean_", "contd_change_relative" = "contd_change_relative_", "contd_change_absolute" = "contd_abs_change_",
    "cypp" = "cypp_", "cpp" = "cpp_", "dpp" = "dpp_",
    "pop" = "pop_size_"
  )
}


# tbl_smmrs_core ----
# Core table summary logic adapted from process_out_Bradford.R tbl_smmrs()
# Handles aggregation, rate calculation, and quantile computation
Simulation$set("private", "tbl_smmrs_core", function(
    tt,                    # data.table with summary data
    what,                  # metric type
    population,            # "ons" or "esp"
    strata,                # strata list (e.g., list("year", c("year", "sex")))
    prbl,                  # quantile probabilities
    baseline_year,         # for _change calculations
    comparator_scenario,   # for comparison metrics
    comparison_starting_year,
    tables_dir,            # output directory
    two_agegrps = FALSE,
    discount_levels = NULL, # one row per level (build_discount_levels)
    discount_from_year = NULL # first year at which discounting bites
) {
  # The money/QALY tables (qalys, net_qalys, costs, net_costs) report EVERY
  # discount level in `discount_levels` side by side, tagged in a `discount`
  # column; the 0% level is the undiscounted figure. Rate/count tables
  # (prevalence, incidence, mortality, ...) have no present value and are left
  # alone. `year` is already in full format here (export_main_tables promotes
  # it), so it is on the same scale as `discount_from_year`.
  if (is.null(discount_levels)) discount_levels <- build_discount_levels(0, 0)
  if (is.null(discount_from_year)) discount_from_year <- baseline_year
  # String mappings for file paths and column patterns (from process_out_Bradford.R)
  str0 <- c(
    "prvl" = "prvl", "prvl_change_relative" = "prvl", "prvl_change_absolute" = "prvl",
    "incd" = "incd", "incd_change_relative" = "incd", "incd_change_absolute" = "incd",
    "ftlt" = "dis_mrtl", "ftlt_change_relative" = "dis_mrtl", "ftlt_change_absolute" = "dis_mrtl",
    "mrtl" = "mrtl", "mrtl_change_relative" = "mrtl", "mrtl_change_absolute" = "mrtl",
    "dis_mrtl" = "dis_mrtl", "dis_mrtl_change_relative" = "dis_mrtl", "dis_mrtl_change_absolute" = "dis_mrtl",
    "qalys" = "qalys", "net_qalys" = "qalys",
    "costs" = "costs", "net_costs" = "costs",
    "cypp" = "prvl", "cpp" = "incd", "dpp" = "mrtl",
    "pop" = "prvl"
  )

  # Column patterns for grep
  str2 <- c(
    "prvl" = "_prvl$|^popsize$", "prvl_change_relative" = "_prvl$|^popsize$", "prvl_change_absolute" = "_prvl$|^popsize$",
    "incd" = "_incd$|^popsize$", "incd_change_relative" = "_incd$|^popsize$", "incd_change_absolute" = "_incd$|^popsize$",
    "ftlt" = "_deaths$|_prvl$", "ftlt_change_relative" = "_deaths$|_prvl$", "ftlt_change_absolute" = "_deaths$|_prvl$",
    "mrtl" = "_mrtl$|^popsize$", "mrtl_change_relative" = "_mrtl$|^popsize$", "mrtl_change_absolute" = "_mrtl$|^popsize$",
    "dis_mrtl" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
    "dis_mrtl_change_relative" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
    "dis_mrtl_change_absolute" = "^nonmodelled_deaths$|^chd_deaths$|^stroke_deaths$|^popsize$",
    "qalys" = "^EQ5D5L$", "net_qalys" = "^EQ5D5L$",
    "costs" = "_cost$|_costs$|^economic_output$", "net_costs" = "_cost$|_costs$|^economic_output$",
    "contd" = "_contd$", "contd_change_relative" = "_contd$", "contd_change_absolute" = "_contd$",
    "cypp" = "_prvl$", "cpp" = "_incd$", "dpp" = "_mrtl$",
    "pop" = "^popsize$"
  )

  # Output column name prefixes (see .tbl_col_prefixes above)
  str3 <- .tbl_col_prefixes()

  # Output file descriptions
  str4 <- c(
    "prvl" = "prevalence by ", "prvl_change_relative" = "prevalence relative change by ", "prvl_change_absolute" = "prevalence absolute change by ",
    "incd" = "incidence by ", "incd_change_relative" = "incidence relative change by ", "incd_change_absolute" = "incidence absolute change by ",
    "ftlt" = "case fatality by ", "ftlt_change_relative" = "case fatality relative change by ", "ftlt_change_absolute" = "case fatality absolute change by ",
    "mrtl" = "all-cause mortality by ", "mrtl_change_relative" = "all-cause mortality relative change by ", "mrtl_change_absolute" = "all-cause mortality absolute change by ",
    "dis_mrtl" = "disease-specific mortality by ",
    "dis_mrtl_change_relative" = "disease-specific mortality relative change by ",
    "dis_mrtl_change_absolute" = "disease-specific mortality absolute change by ",
    "qalys" = "QALYs by ", "net_qalys" = "net QALYs by ",
    "costs" = "costs by ", "net_costs" = "net costs by ",
    "contd" = "continuous outcome by ",
    "contd_change_relative" = "continuous outcome relative change by ",
    "contd_change_absolute" = "continuous outcome absolute change by ",
    "cypp" = "case-years prevented or postponed by ",
    "cpp" = "cases prevented or postponed by ",
    "dpp" = "deaths prevented or postponed by ",
    "pop" = "pop size by "
  )

  # Derived, not passed as a formal: it can never desynchronise from the
  # argument it describes, and the tinytest harnesses that re-parent and call
  # these private methods directly keep working with no edit.
  comparator_is_map <- .comparator_is_map(comparator_scenario)

  # Per-scenario comparator map. A bare `comparator_scenario` resolves to the
  # historical {every other scenario -> it}. Resolved from `tt` once: every `d`
  # below is an aggregation of `tt`, so it carries the same scenarios.
  cmap <- .resolve_comparators(comparator_scenario, unique(tt$scenario))

  # Split a table into intervention rows and comparator rows. These are NO
  # LONGER complementary: under a map a scenario can be BOTH (an intervention
  # against sc0 AND the comparator for another arm), so `cmp` must be taken
  # from the FULL table rather than from the complement of `d`.
  split_cmp <- function(d, x, extra_on) {
    cmp <- d[scenario %chin% unique(unname(cmap)) &
               year >= comparison_starting_year]
    dd  <- d[scenario %chin% names(cmap) & year >= comparison_starting_year]
    # as.character() is MANDATORY: indexing a named vector with a FACTOR uses
    # the LEVEL CODES and mis-maps silently.
    dd[, .cmp__ := unname(cmap[as.character(scenario)])]
    list(d = dd, cmp = cmp,
         on = c(".cmp__" = "scenario", setdiff(x, "scenario"), extra_on))
  }

  # Add mc and scenario to strata
  strata <- lapply(strata, function(x) c("mc", "scenario", x))

  # Process each strata combination
  lapply(strata, function(x) {
    if (grepl("^qalys$", what)) {
      # QALYs processing (EQ5D5L only - HUI3 not implemented)
      d <- tt[, .("EQ5D5L" = sum(EQ5D5L)), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "scale", value.name = "QALYs")
      # One block of rows per discount level, tagged in a `discount` column.
      # Discounting the annual flows before cumsum() makes the cumulative
      # column a sum of present values.
      d <- expand_discount_levels(d, "QALYs", discount_levels,
                                  discount_from_year, "qaly")
      # `discount` MUST join id.vars below for the same reason `.cmp__` must be
      # dropped in the net_qalys branch: that melt names id.vars but not
      # measure.vars, so a column outside id.vars is absorbed as a measure.
      x <- c(x, "scale", "discount")
      setkeyv(d, c(setdiff(x, "year"), "year"))
      d[, cumulative := cumsum(QALYs), keyby = setdiff(x, "year")]
      d <- melt(d, id.vars = x, variable.name = "type")
      d[, type := fifelse(type == "cumulative", "QALYs_cuml", "QALYs")]
      setkey(d, "type", "scale")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
             keyby = eval(setdiff(x, "mc"))]
      setnames(d, c(setdiff(x, "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(x, "mc")))
      setcolorder(d, setdiff(x, "mc"))

    } else if (grepl("^net_qalys$", what)) {
      # Net QALYs (intervention - baseline, EQ5D5L only - HUI3 not implemented)
      d <- tt[, .("EQ5D5L" = sum(EQ5D5L)), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "scale", value.name = "QALYs")
      p <- split_cmp(d, x, "scale")
      d <- p$d[p$cmp, on = p$on, net_QALYs := QALYs - i.QALYs]
      # `.cmp__` MUST go before the melt() below: that melt names id.vars but
      # not measure.vars, so a stray helper column is absorbed as a measure
      # variable, doubling the rows and coercing `value` to character.
      d[, c("QALYs", ".cmp__") := NULL]
      # Discounting is linear in the flows, so discounting the difference is
      # identical to differencing the two discounted arms (both share `year`).
      # Runs AFTER `.cmp__` is dropped so the helper column is not replicated.
      d <- expand_discount_levels(d, "net_QALYs", discount_levels,
                                  discount_from_year, "qaly")
      x <- c(x, "scale", "discount")
      setkeyv(d, c(setdiff(x, "year"), "year"))
      d[, cumulative := cumsum(net_QALYs), keyby = setdiff(x, "year")]
      d <- melt(d, id.vars = x, variable.name = "type")
      d[type == "cumulative", type := "net_QALYs_cuml"]
      setkey(d, "type", "scale")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = FALSE),
             keyby = eval(setdiff(x, "mc"))]
      setnames(d, c(setdiff(x, "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(x, "mc")))
      setcolorder(d, setdiff(x, "mc"))

    } else if (grepl("^costs$", what)) {
      # Costs processing
      d <- tt[, lapply(.SD, sum), .SDcols = patterns("_cost$|_costs$|^economic_output$"), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "costs_type", value.name = "costs")
      d <- expand_discount_levels(d, "costs", discount_levels,
                                  discount_from_year, "cost")
      x <- c(x, "costs_type", "discount")
      setkeyv(d, c(setdiff(x, "year"), "year"))
      d[, cumulative := cumsum(costs), keyby = setdiff(x, "year")]
      d <- melt(d, id.vars = x, variable.name = "type")
      d[type == "cumulative", type := "costs_cuml"]
      setkey(d, "type", "costs_type")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = TRUE),
             keyby = eval(setdiff(x, "mc"))]
      setnames(d, c(setdiff(x, "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(x, "mc")))
      setcolorder(d, setdiff(x, "mc"))

    } else if (grepl("^net_costs$", what)) {
      # Net costs (intervention - baseline)
      d <- tt[, lapply(.SD, sum), .SDcols = patterns("_cost$|_costs$|^economic_output$"), keyby = eval(x)]
      d <- melt(d, id.vars = x, variable.name = "costs_type", value.name = "value")
      p <- split_cmp(d, x, "costs_type")
      d <- p$d[p$cmp, on = p$on, net_costs := value - i.value]
      d[, c("value", ".cmp__") := NULL]
      # After `.cmp__` is dropped, so the helper column is not replicated.
      d <- expand_discount_levels(d, "net_costs", discount_levels,
                                  discount_from_year, "cost")
      x <- c(x, "costs_type", "discount")
      setkeyv(d, c(setdiff(x, "year"), "year"))
      d[, cumulative := cumsum(net_costs), keyby = setdiff(x, "year")]
      d <- melt(d, id.vars = x, variable.name = "type")
      d[type == "cumulative", type := "net_costs_cuml"]
      setkey(d, "type", "costs_type")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(type), rounding = TRUE),
             keyby = eval(setdiff(x, "mc"))]
      setnames(d, c(setdiff(x, "mc"), "type", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, c("type", setdiff(x, "mc")))
      setcolorder(d, setdiff(x, "mc"))

    } else if (grepl("^contd", what)) {
      # User-defined continuous outcomes (*_contd columns). The contd summaries
      # already hold a population-weighted mean per stratum, so collapsing to
      # coarser strata re-weights by popsize (NOT a sum, and no division by
      # popsize).
      contd_cols <- grep("_contd$", names(tt), value = TRUE)
      d <- tt[, lapply(.SD, function(v) weighted.mean(v, popsize, na.rm = TRUE)),
              .SDcols = contd_cols, keyby = eval(x)]
      d <- melt(d, id.vars = x)

      if (grepl("_change_relative$", what)) {
        # Relative change from baseline year
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year"), "variable"), value := value / i.value]
      }
      if (grepl("_change_absolute$", what)) {
        # Absolute change from baseline year
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year"), "variable"), value := value - i.value]
      }

      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable), rounding = FALSE),
             keyby = eval(setdiff(x, "mc"))]
      setnames(d, c(setdiff(x, "mc"), "outcome", scales::percent(prbl, prefix = str3[[what]])))
      setkeyv(d, setdiff(x, "mc"))
      setcolorder(d, setdiff(x, "mc"))

    } else {
      # All other metrics (prvl, incd, mrtl, dis_mrtl, ftlt, pop, cypp, cpp, dpp)
      d <- tt[, lapply(.SD, sum), .SDcols = patterns(str2[[what]]), keyby = x]

      # Convert integer columns to numeric. Same guard, same reason, as the one
      # export_all_cause_mrtl_tables() applies -- shared so the two cannot
      # drift apart.
      .promote_integer_measures(d, x)

      if (grepl("^ftlt", what)) {
        # Case fatality: deaths / prevalence
        nm <- names(d)
        nm <- grep("_deaths$", nm, value = TRUE)
        nm <- gsub("_deaths$", "", nm)
        nm <- setdiff(nm, "alive")
        for (i in nm) {
          set(d, NULL, paste0(i, "_ftlt"),
              d[[paste0(i, "_deaths")]] / d[[paste0(i, "_prvl")]])
        }
        nm <- names(d)
        nm <- grep("_deaths$|_prvl$", nm, value = TRUE)
        d[, (nm) := NULL]
        setnafill(d, "const", 0, cols = grep("_ftlt$", names(d), value = TRUE))
      } else if (!what %in% c("pop", "cypp", "cpp", "dpp")) {
        # Calculate rates (divide by popsize)
        d <- d[, lapply(.SD, function(y) y / popsize), keyby = x]
      }

      d <- melt(d, id.vars = x)

      if (grepl("_change_relative$", what)) {
        # Relative change from baseline year
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year"), "variable"), value := value / i.value]
      }

      if (grepl("_change_absolute$", what)) {
        # Absolute change from baseline year
        d19 <- d[year == baseline_year][, year := NULL]
        d[d19, on = c(setdiff(x, "year"), "variable"), value := value - i.value]
      }

      if (grepl("^cypp$|^cpp$|^dpp$", what)) {
        # Comparison metrics: baseline - intervention
        p <- split_cmp(d, x, "variable")
        d <- p$d
        # This is the ONLY contrast site whose update TARGET is `value` itself.
        # data.table's x[i, :=] touches only MATCHED rows, so an intervention
        # cell with no comparator counterpart would silently keep its raw LEVEL
        # -- a prevalence or incidence count published, with the WRONG SIGN, in
        # a file called "cases prevented or postponed".
        #
        # The fix is to treat the missing cell as ZERO, because that is what
        # absence MEANS here: `d` is a sum aggregation keyed by the strata, so a
        # combination with no row is a combination with no events, not one with
        # an unknown count. The contrast is therefore 0 - value = -value: the
        # intervention arm has cases the comparator arm does not, i.e. cases
        # CAUSED. That is a real quantity and must not be erased -- NA-ing these
        # rows would drop the draw from that cell's quantiles for this year and,
        # via the cumsum() below, every later year of the same series, turning a
        # genuine harm into an apparently certain null.
        #
        # Cell coverage differs between arms only in SPARSE runs -- small `n`,
        # few Monte-Carlo draws, fine strata -- where differential survival
        # leaves a (year, agegrp, sex, dimd) cell occupied in one arm and empty
        # in the other. Measured as exactly zero across 131M rows of a
        # 100-iteration production run, so this changes no production output.
        miss_idx <- d[!p$cmp, on = p$on, which = TRUE]
        d[p$cmp, on = p$on, value := i.value - value]
        if (length(miss_idx) > 0L) {
          set(d, miss_idx, "value", -d$value[miss_idx])
          if (self$design$sim_prm$logs) {
            message("    ", what, ": ", length(miss_idx), " of ", nrow(d),
                    " row(s) have no comparator cell; treated as a comparator ",
                    "count of zero (contrast = -intervention)")
          }
        }
        set(d, j = ".cmp__", value = NULL)   # MUST precede the melt() below
        d[, variable := gsub(paste0("_", str0[[what]]), "", variable)]
        setkeyv(d, c(x[x != "year"], "variable", "year"))
        d[, cumulative := cumsum(value), keyby = c(setdiff(x, "year"), "variable")]
        d <- melt(d, id.vars = c(x, "variable"), variable.name = "type")
        d[, type := fifelse(type == "cumulative", paste0(what, "_cuml"), what)]
        setkey(d, "type", "variable")
        d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable),
                                rounding = (what %in% c("pop", "cypp", "cpp", "dpp"))),
               keyby = eval(setdiff(c(x, "type"), "mc"))]
        x <- c(x, "type")
        setnames(d, c(setdiff(x, "mc"), "disease", scales::percent(prbl, prefix = str3[[what]])))
      } else {
        setkey(d, "variable")
        d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable),
                                rounding = what == "pop"),
               keyby = eval(setdiff(x, "mc"))]
        setnames(d, c(setdiff(x, "mc"), "disease", scales::percent(prbl, prefix = str3[[what]])))
      }

      if (what == "pop") {
        d[, disease := NULL]
      } else {
        if ("popsize" %in% d$disease) d <- d[disease != "popsize"]
      }
      setkeyv(d, setdiff(x, "mc"))
      setcolorder(d, setdiff(x, "mc"))
    }

    # Build output filename
    str5 <- c(
      "ons" = " (not standardised).csv",
      "esp" = paste0(" (", paste(setdiff(c("mc", "scenario", "year", "age", "sex"), x),
                                 collapse = "-"), " standardised).csv")
    )
    # Discount levels live in the `discount` column rather than in the file
    # name, so the file set is the same whether one or several are requested.
    str6 <- paste0(
      str4[[what]],
      paste(setdiff(x, c("mc", "scenario", "comparator", "type", "scale",
                         "costs_type", "discount")), collapse = "-"),
      str5[[population]]
    )

    # `comparator` is a self-describing column, exactly like `gradient` in the
    # equity tables, and is written ONLY in map mode so the scalar output stays
    # byte-identical. Attached HERE, after quantiling, because one comparator
    # per scenario makes it a FUNCTION of `scenario` -- so it never enters `x`
    # and therefore never touches a group-by, one of the positional setnames()
    # calls, or a filename.
    if (comparator_is_map &&
        what %chin% c("cypp", "cpp", "dpp", "net_qalys", "net_costs")) {
      d[, comparator := unname(cmap[as.character(scenario)])]
      nms <- setdiff(names(d), "comparator")
      setcolorder(d, append(nms, "comparator", after = match("scenario", nms)))
    }

    fwrite(d, file.path(tables_dir, str6))
  })

  invisible(NULL)
})


# export_main_tables ----
# Generate main summary tables (prevalence, incidence, mortality, etc.)
# Memory-optimized: processes by source dataset to minimize simultaneous memory
Simulation$set("private", "export_main_tables", function(
    prbl,
    baseline_year,
    output_dir,
    tables_dir,
    comparator_scenario = "sc0",
    two_agegrps = FALSE,
    two_agegrps_split_age = 65L,
    discount_levels = NULL,
    discount_from_year = NULL,
    strata_ons = NULL,
    strata_esp = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating main summary tables...")
  }

  # Undiscounted-only fallback, so the method still works if called directly.
  if (is.null(discount_levels)) discount_levels <- build_discount_levels(0, 0)
  if (is.null(discount_from_year)) discount_from_year <- baseline_year

  # String mappings for source datasets
  str0 <- c(
    "prvl" = "prvl", "prvl_change_relative" = "prvl", "prvl_change_absolute" = "prvl",
    "incd" = "incd", "incd_change_relative" = "incd", "incd_change_absolute" = "incd",
    "ftlt" = "dis_mrtl", "ftlt_change_relative" = "dis_mrtl", "ftlt_change_absolute" = "dis_mrtl",
    "mrtl" = "mrtl", "mrtl_change_relative" = "mrtl", "mrtl_change_absolute" = "mrtl",
    "dis_mrtl" = "dis_mrtl", "dis_mrtl_change_relative" = "dis_mrtl", "dis_mrtl_change_absolute" = "dis_mrtl",
    "qalys" = "qalys", "net_qalys" = "qalys",
    "costs" = "costs", "net_costs" = "costs",
    "cypp" = "prvl", "cpp" = "incd", "dpp" = "mrtl",
    "pop" = "prvl"
  )
  str1 <- c("ons" = "scaled_up", "esp" = "esp")

  # Group metrics by source dataset for efficient memory usage
  source_to_metrics <- list(
    prvl = c("prvl", "prvl_change_relative", "prvl_change_absolute", "cypp", "pop"),
    incd = c("incd", "incd_change_relative", "incd_change_absolute", "cpp"),
    dis_mrtl = c("ftlt", "ftlt_change_relative", "ftlt_change_absolute", "dis_mrtl", "dis_mrtl_change_relative", "dis_mrtl_change_absolute"),
    mrtl = c("mrtl", "mrtl_change_relative", "mrtl_change_absolute", "dpp"),
    qalys = c("qalys", "net_qalys"),
    costs = c("costs", "net_costs"),
    # User-defined continuous outcomes (*_contd columns). Skipped silently when
    # the contd summary was not produced (read_summary_dataset returns NULL).
    contd = c("contd", "contd_change_relative", "contd_change_absolute")
  )

  # Process each source dataset group
  for (source_name in names(source_to_metrics)) {
    metrics_for_source <- source_to_metrics[[source_name]]

    # Process both populations for this source
    for (pop_name in c("ons", "esp")) {
      pop_key <- str1[[pop_name]]

      # Load dataset once for this source/population combination
      tt_base <- private$read_summary_dataset(source_name, pop_key)
      if (is.null(tt_base)) next
      # Convert year from short format (19) to full format (2019)
      tt_base[, year := year + 2000L]
      # Allow `qimd` strata even when the summaries only carry `dimd`. Every
      # aggregation below selects value columns by anchored pattern, so the
      # extra column is inert unless it is named in the strata.
      add_qimd_from_dimd(tt_base)

      # Also load prvl for ftlt metrics (case fatality denominator)
      prvl_for_ftlt <- NULL
      if (source_name == "dis_mrtl") {
        prvl_for_ftlt <- private$read_summary_dataset("prvl", pop_key)
        if (!is.null(prvl_for_ftlt)) {
          prvl_for_ftlt[, year := year + 2000L]
          add_qimd_from_dimd(prvl_for_ftlt)
        }
      }

      # Process each metric that uses this source
      for (what in metrics_for_source) {
        # Skip pop for esp
        if (what == "pop" && pop_name == "esp") next
        if (grepl("_age", what) && pop_name == "esp") next

        # Use configurable strata (already filtered by two_agegrps in build_strata_config)
        if (pop_name == "ons") {
          strata <- strata_ons
        } else {
          strata <- strata_esp
        }

        if (self$design$sim_prm$logs) {
          message(paste0("  ", what, "-", pop_name))
        }

        # Check if comparison metrics can be computed. This guard runs BEFORE
        # the copy() below, deliberately: on a single-scenario run EVERY
        # comparison metric is skipped, and copying `tt_base` first only to
        # rm() it again allocated a full duplicate of the summary for nothing --
        # ten times per call (five metrics x two populations), 0.6-2.2 GB each
        # on a production run. Reading the scenario list off `tt_base` is
        # equivalent: `tt` was a verbatim copy and nothing had touched it yet.
        comparison_metrics <- c("cypp", "cpp", "dpp", "net_qalys", "net_costs")
        if (what %in% comparison_metrics) {
          available_scenarios <- unique(tt_base$scenario)
          if (.comparator_is_map(comparator_scenario)) {
            cm <- .resolve_comparators(comparator_scenario, available_scenarios)
            drops <- attr(cm, "dropped")
            if (length(drops) > 0L) {
              warning(sprintf(paste0(
                "main tables: comparator pair(s) %s cannot be built from this ",
                "summary (found: %s); those contrasts were skipped for %s."),
                paste(sQuote(drops), collapse = ", "),
                paste(available_scenarios, collapse = ", "), what),
                call. = FALSE, immediate. = TRUE)
            }
            if (length(cm) == 0L) {
              next
            }
          } else {
            # A comparator absent from this summary leaves every contrast row
            # unmatched, so the file comes out empty. That used to be silent --
            # and before the join was corrected it was worse than silent, since
            # the unmatched rows were written out with blank values for BOTH
            # arms. Say it plainly; the equity family already does.
            if (!comparator_scenario %chin% available_scenarios) {
              warning(sprintf(paste0(
                "main tables: comparator scenario '%s' is not in this summary ",
                "(found: %s), so '%s' has no rows to contrast and its file will ",
                "be empty."),
                comparator_scenario, paste(available_scenarios, collapse = ", "),
                what), call. = FALSE, immediate. = TRUE)
            }
            ## ---- ORIGINAL, VERBATIM, bar the dropped `rm(tt)` -----------
            non_comparator_scenarios <- setdiff(available_scenarios, comparator_scenario)
            if (length(non_comparator_scenarios) == 0) {
              if (self$design$sim_prm$logs) {
                message("    Skipping ", what, " - no intervention scenarios (only '",
                        comparator_scenario, "' found)")
              }
              next
            }
            ## -------------------------------------------------------------
          }
        }

        # Working copy, taken only once the metric is known to be buildable.
        # It MUST be a copy: the two_agegrps relabel below, the absorb_dt() in
        # the ftlt branch and tbl_smmrs_core() all mutate `tt` BY REFERENCE,
        # while `tt_base` is reused by every later metric in this source group.
        tt <- copy(tt_base)

        # Handle two_agegrps transformation
        if (two_agegrps) collapse_two_agegrps(tt, two_agegrps_split_age)

        # For case fatality, add prevalence denominator. The denominator MUST be
        # collapsed on the same split, or absorb_dt()'s natural join on `agegrp`
        # would match nothing and every case-fatality row would be dropped.
        if (grepl("^ftlt", what) && !is.null(prvl_for_ftlt)) {
          t1 <- copy(prvl_for_ftlt)
          setnames(t1, "popsize", "nonmodelled_prvl")
          if (two_agegrps) collapse_two_agegrps(t1, two_agegrps_split_age)
          absorb_dt(tt, t1)
          tt <- tt[nonmodelled_prvl > 0]
          rm(t1)
        }

        # Generate tables. The qalys / net_qalys / costs / net_costs metrics are
        # expanded to one block of rows per discount level, tagged in a
        # `discount` column; the other metrics here are rates and counts and are
        # left untouched. (Discounting stopped being CEA-only when multi-level
        # discounting landed -- `discount_levels` is passed straight below.)
        private$tbl_smmrs_core(
          tt = tt,
          what = what,
          population = pop_name,
          strata = strata,
          prbl = prbl,
          baseline_year = baseline_year,
          comparator_scenario = comparator_scenario,
          comparison_starting_year = baseline_year,
          tables_dir = tables_dir,
          two_agegrps = two_agegrps,
          discount_levels = discount_levels,
          discount_from_year = discount_from_year
        )
        rm(tt)
      }

      # Cleanup after processing this source/population
      rm(tt_base)
      if (!is.null(prvl_for_ftlt)) rm(prvl_for_ftlt)
    }

    # Garbage collect after each source group
    gc(verbose = FALSE)
  }

  invisible(NULL)
})


# export_all_cause_mrtl_tables ----
# Generate all-cause mortality by disease tables
# Memory-optimized: loads datasets once and reuses them
Simulation$set("private", "export_all_cause_mrtl_tables", function(
    prbl,
    summaries_dir,
    tables_dir,
    two_agegrps = FALSE,
    two_agegrps_split_age = 65L,
    strata_ons = NULL,
    strata_esp = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating all-cause mortality by disease tables...")
  }

  # Helper function to convert user strata to internal format with mc and scenario
  make_strata_configs <- function(strata_list, standardised = FALSE) {
    lapply(strata_list, function(s) {
      outstrata <- c("mc", s, "scenario")
      suffix <- paste(s, collapse = "-")
      # Map agegrp to agegroup in suffix for backwards compatibility
      suffix <- gsub("agegrp", "agegroup", suffix)
      if (standardised) {
        # Determine what was standardised by (variables NOT in strata).
        # `qimd` retains the deprivation axis just as `dimd` does (coarser),
        # so a qimd stratum must not be reported as standardised over dimd.
        possible_vars <- c("age", "sex", "dimd")
        s_axes <- if ("qimd" %in% s) c(s, "dimd") else s
        standardised_vars <- setdiff(possible_vars, s_axes)
        std_suffix <- paste(standardised_vars, collapse = "-")
        list(strata = outstrata, suffix = suffix, std = std_suffix)
      } else {
        list(strata = outstrata, suffix = suffix)
      }
    })
  }

  # Load datasets once (avoid duplicate reads)
  tt_scaled <- private$read_summary_dataset("all_cause_mrtl_by_dis", "scaled_up")
  pp_scaled <- private$read_summary_dataset("prvl", "scaled_up")
  tt_esp <- private$read_summary_dataset("all_cause_mrtl_by_dis", "esp")

  # Convert year from short format (19) to full format (2019)
  if (!is.null(tt_scaled)) tt_scaled[, year := year + 2000L]
  if (!is.null(pp_scaled)) pp_scaled[, year := year + 2000L]
  if (!is.null(tt_esp)) tt_esp[, year := year + 2000L]
  # Allow `qimd` strata even when the summaries only carry `dimd`
  add_qimd_from_dimd(tt_scaled)
  add_qimd_from_dimd(pp_scaled)
  add_qimd_from_dimd(tt_esp)
  # Match export_main_tables(): under two_agegrps, `agegrp` means the coarse
  # split. `pp_scaled` is the population denominator and is joined to
  # `tt_scaled` on the strata, so it MUST be collapsed on the same split.
  if (two_agegrps) {
    collapse_two_agegrps(tt_scaled, two_agegrps_split_age)
    collapse_two_agegrps(pp_scaled, two_agegrps_split_age)
    collapse_two_agegrps(tt_esp, two_agegrps_split_age)
  }

  # ---- Non-standardised with disease denominator ----
  if (!is.null(tt_scaled)) {
    strata_configs <- make_strata_configs(strata_ons, standardised = FALSE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- tt_scaled[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"), keyby = eval(outstrata)]
      # Before the melt, so `value` is born numeric and the division below is
      # not truncated in place. See .promote_integer_measures().
      .promote_integer_measures(d, outstrata)
      d <- melt(d, id.vars = outstrata)
      cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
      d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
      d[cases, on = c(outstrata, "variable"), value := value / i.value]
      rm(cases)  # Cleanup intermediate
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "disease",
                    scales::percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("all-cause mortality given disease-", cfg$suffix, " (not standardised).csv")))
      rm(d)
    }
  }

  # ---- Non-standardised with population denominator ----
  if (!is.null(tt_scaled) && !is.null(pp_scaled)) {
    strata_configs <- make_strata_configs(strata_ons, standardised = FALSE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      cases <- pp_scaled[, lapply(.SD, sum), .SDcols = patterns("^popsize$"), keyby = eval(outstrata)]
      d <- tt_scaled[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"), keyby = eval(outstrata)]
      .promote_integer_measures(cases, outstrata)
      .promote_integer_measures(d, outstrata)
      d <- melt(d, id.vars = outstrata)
      d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
      d[cases, on = outstrata, value := value / popsize]
      rm(cases)  # Cleanup intermediate
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "disease",
                    scales::percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("all-cause mortality given disease-", cfg$suffix, " popdenom (not standardised).csv")))
      rm(d)
    }
  }

  # Cleanup scaled_up datasets before ESP processing
  rm(tt_scaled, pp_scaled)
  gc(verbose = FALSE)

  # ---- Standardised (ESP) ----
  if (!is.null(tt_esp)) {
    strata_configs <- make_strata_configs(strata_esp, standardised = TRUE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- tt_esp[, lapply(.SD, sum), .SDcols = patterns("^deaths_|^cases_"), keyby = eval(outstrata)]
      .promote_integer_measures(d, outstrata)
      d <- melt(d, id.vars = outstrata)
      cases <- d[grep("^cases_", variable)][, variable := gsub("^cases_", "", variable)]
      d <- d[grep("^deaths_", variable)][, variable := gsub("^deaths_", "", variable)]
      d[cases, on = c(outstrata, "variable"), value := value / i.value]
      rm(cases)  # Cleanup intermediate
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "disease",
                    scales::percent(prbl, prefix = "all_cause_mrtl_by_disease_rate_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("all-cause mortality given disease-", cfg$suffix,
                                 " (", cfg$std, " standardised).csv")))
      rm(d)
    }
  }

  rm(tt_esp)
  invisible(NULL)
})


# export_disease_characteristics_tables ----
# Generate disease characteristics tables (duration, age metrics, CMS)
Simulation$set("private", "export_disease_characteristics_tables", function(
    prbl,
    summaries_dir,
    tables_dir,
    strata = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating disease characteristics tables...")
  }

  tt <- private$read_summary_dataset("dis_characteristics", "scaled_up")
  if (is.null(tt)) return(invisible(NULL))

  # Convert year from short format (19) to full format (2019)
  tt[, year := year + 2000L]
  # Allow `qimd` strata even when the summaries only carry `dimd`
  add_qimd_from_dimd(tt)

  # Type conversion if needed
  if ("mean_cms_count_cms1st_cont" %in% names(tt)) {
    tt[, mean_cms_count_cms1st_cont := as.numeric(mean_cms_count_cms1st_cont)]
  }

  # Derive id variables dynamically from requested strata, keeping only
  # columns that actually exist in the data (e.g. agegrp is excluded from

  # disease characteristics summaries via strata_noagegrp)
  all_strata_vars <- unique(unlist(strata))
  id_vars <- intersect(unique(c("mc", "scenario", all_strata_vars)), names(tt))
  id_pattern <- paste(id_vars, collapse = "|")

  # Extract case counts for weighting
  d1 <- tt[, .SD, .SDcols = patterns(paste0(id_pattern, "|^cases_"))]
  d1 <- melt(d1, id.vars = id_vars)
  d1 <- unique(d1, by = c(id_vars, "variable"))
  d1[, disease := gsub("^cases_", "", variable)]
  d1[, variable := NULL]

  # Extract characteristics columns
  char_patterns <- paste0(id_pattern, "|^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_")
  tt <- tt[, .SD, .SDcols = patterns(char_patterns)]

  if ("mean_cms_count_cmsmm1" %in% names(tt)) {
    tt[, mean_cms_count_cmsmm1 := as.double(mean_cms_count_cmsmm1)]
  }

  tt <- melt(tt, id.vars = id_vars)
  tt[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
  tt[d1, on = c(id_vars, "disease"), cases := i.value]

  # Helper function to convert user strata to internal format
  make_strata_configs <- function(strata_list) {
    lapply(strata_list, function(s) {
      outstrata <- c("mc", s, "scenario")
      suffix <- paste(s, collapse = "-")
      list(strata = outstrata, suffix = suffix)
    })
  }

  # Process each strata
  strata_configs <- make_strata_configs(strata)

  for (cfg in strata_configs) {
    outstrata <- cfg$strata
    d <- tt[, weighted.mean(value, cases, na.rm = TRUE), keyby = c(outstrata, "variable")]
    setkey(d, "variable")
    d <- d[, safe_fquantile_byid(V1, prbl, id = as.character(variable)),
           keyby = eval(setdiff(outstrata, "mc"))]
    setnames(d, c(setdiff(outstrata, "mc"), "variable", scales::percent(prbl, prefix = "value_")))

    # Parse variable name to extract disease and type
    d[, disease := gsub("^mean_duration_|^mean_age_incd_|^mean_age_1st_onset_|^mean_age_prvl_|^mean_cms_score_|^mean_cms_count_", "", variable)]
    d[grep("^mean_duration_", variable), type := "mean_duration"]
    d[grep("^mean_age_incd_", variable), type := "mean_age_incd"]
    d[grep("^mean_age_1st_onset_", variable), type := "mean_age_1st_onset"]
    d[grep("^mean_age_prvl_", variable), type := "mean_age_prvl"]
    d[grep("^mean_cms_score_", variable), type := "mean_cms_score"]
    d[grep("^mean_cms_count_", variable), type := "mean_cms_count"]
    d[, variable := NULL]
    setkeyv(d, c(setdiff(outstrata, "mc"), "disease", "type"))
    setcolorder(d)

    fwrite(d, file.path(tables_dir,
                        paste0("disease characteristics by ", cfg$suffix, " (not standardised).csv")))
    rm(d)  # Cleanup after each output
  }

  # Cleanup
  rm(tt, d1)
  invisible(NULL)
})


# export_xps_tables ----
# Generate exposure summary tables
# Memory-optimized: adds cleanup between sections
Simulation$set("private", "export_xps_tables", function(
    prbl,
    output_dir,
    tables_dir,
    strata_ons = NULL,
    strata_esp = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating exposure tables...")
  }

  # Helper function to build filter expression from strata
  # Variables in strata -> filter on != "All"
  # Variables not in strata -> filter on == "All"
  make_xps_strata_configs <- function(strata_list, filterable_vars, standardised = FALSE) {
    lapply(strata_list, function(s) {
      outstrata <- c("mc", s, "scenario")
      suffix <- paste(setdiff(s, "year"), collapse = "-")
      if (suffix == "") suffix <- "year" else suffix <- paste0("year-", suffix)
      # Map agegrp20 to agegroup in suffix for backwards compatibility
      suffix <- gsub("agegrp20", "agegroup", suffix)

      # Build filter expression: vars in strata -> != "All", vars not in strata -> == "All"
      filter_parts <- character(0)
      for (v in filterable_vars) {
        if (v %in% s) {
          filter_parts <- c(filter_parts, paste0(v, " != 'All'"))
        } else {
          filter_parts <- c(filter_parts, paste0(v, " == 'All'"))
        }
      }
      filter_str <- paste(filter_parts, collapse = " & ")
      filter_expr <- if (length(filter_parts) > 0) parse(text = filter_str)[[1]] else quote(TRUE)

      if (standardised) {
        # Determine what was standardised by: all filterable vars not in the strata
        # Convert agegrp20 to age for standardisation naming
        s_for_std <- gsub("agegrp20", "age", s)
        all_std_vars <- gsub("agegrp20", "age", filterable_vars)
        standardised_vars <- setdiff(all_std_vars, s_for_std)
        std_suffix <- paste(standardised_vars, collapse = "-")
        list(strata = outstrata, suffix = suffix, filter_expr = filter_expr, std = std_suffix)
      } else {
        list(strata = outstrata, suffix = suffix, filter_expr = filter_expr)
      }
    })
  }

  # Helper to detect filterable vars: columns with "All" values from groupingsets
  detect_filterable_vars <- function(dt) {
    candidates <- setdiff(names(dt), c("mc", "scenario", "year"))
    candidates[vapply(candidates, function(v) {
      is.character(dt[[v]]) && "All" %in% dt[[v]]
    }, logical(1))]
  }

  # ---- Non-standardised (xps20) ----
  xps_path <- file.path(output_dir, "xps", "xps20")
  if (dir.exists(xps_path)) {
    xps_tab <- CKutils::read_parquet_dt(xps_path)
    # Convert year from short format (19) to full format (2019)
    xps_tab[, year := year + 2000L]

    # Identify exposure columns: _curr_xps columns plus any extra from exposures_for_output
    xps_cols <- grep("_curr_xps$", names(xps_tab), value = TRUE)
    extra_in_data <- intersect(self$design$sim_prm$exposures_for_output, names(xps_tab))
    extra_in_data <- extra_in_data[
      vapply(extra_in_data, function(col) is.numeric(xps_tab[[col]]), logical(1))
    ]
    xps_cols <- unique(c(xps_cols, extra_in_data))

    # Detect which columns have "All" marginals from groupingsets
    filt_vars <- detect_filterable_vars(xps_tab)
    strata_configs <- make_xps_strata_configs(strata_ons, filt_vars, standardised = FALSE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- xps_tab[eval(cfg$filter_expr)]
      if (nrow(d) == 0) next
      d <- d[, lapply(.SD, mean), .SDcols = xps_cols, keyby = eval(outstrata)]
      d <- melt(d, id.vars = outstrata)
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "exposure", scales::percent(prbl, prefix = "xps_mean_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("exposures by ", cfg$suffix, " (not standardised).csv")))
      rm(d)  # Cleanup after each output
    }
    rm(xps_tab)  # Cleanup before next section
    gc(verbose = FALSE)
  }

  # ---- Standardised (xps5 / ESP) ----
  xps_path <- file.path(output_dir, "xps", "xps5")
  if (dir.exists(xps_path)) {
    xps_tab <- CKutils::read_parquet_dt(xps_path)
    # Convert year from short format (19) to full format (2019)
    xps_tab[, year := year + 2000L]

    # Identify exposure columns: _curr_xps columns plus any extra from exposures_for_output
    xps_cols <- grep("_curr_xps$", names(xps_tab), value = TRUE)
    extra_in_data <- intersect(self$design$sim_prm$exposures_for_output, names(xps_tab))
    extra_in_data <- extra_in_data[
      vapply(extra_in_data, function(col) is.numeric(xps_tab[[col]]), logical(1))
    ]
    xps_cols <- unique(c(xps_cols, extra_in_data))

    # Detect which columns have "All" marginals from groupingsets
    filt_vars <- detect_filterable_vars(xps_tab)
    strata_configs <- make_xps_strata_configs(strata_esp, filt_vars, standardised = TRUE)

    for (cfg in strata_configs) {
      outstrata <- cfg$strata
      d <- xps_tab[eval(cfg$filter_expr)]
      if (nrow(d) == 0) next
      d <- d[, lapply(.SD, mean), .SDcols = xps_cols, keyby = eval(outstrata)]
      d <- melt(d, id.vars = outstrata)
      setkey(d, "variable")
      d <- d[, safe_fquantile_byid(value, prbl, id = as.character(variable)),
             keyby = eval(setdiff(outstrata, "mc"))]
      setnames(d, c(setdiff(outstrata, "mc"), "exposure", scales::percent(prbl, prefix = "xps_mean_")))
      setkeyv(d, setdiff(outstrata, "mc"))
      fwrite(d, file.path(tables_dir,
                          paste0("exposures by ", cfg$suffix, " (", cfg$std, " standardised).csv")))
      rm(d)  # Cleanup after each output
    }
    rm(xps_tab)  # Final cleanup
  }

  invisible(NULL)
})


# export_cea_tables ----
# Generate cost-effectiveness (ICER / NMB) tables from the qalys and costs
# summaries. For each stratum, perspective and QALY scale (EQ5D5L; HUI3 if
# present) it computes, per Monte-Carlo iteration, the cumulative discounted
# incremental QALYs and costs versus the comparator scenario, then the ICER and
# the net monetary benefit (NMB) at each willingness-to-pay threshold, and
# quantiles those across iterations.
#
# Cost perspectives (one file each, named in the filename):
#   societal              = total_cost (+ ALL user *_costs columns).
#                           total_cost already nets economic_output out of
#                           healthcare + socialcare + informalcare costs.
#   healthcare            = healthcare_cost (+ the user *_costs columns named in
#                           custom_costs_in_healthcare) -- direct treatment
#                           costs only.
#   healthcare_socialcare = healthcare_cost + socialcare_cost (+ the same user
#                           columns). This is NICE's reference case, "NHS and
#                           Personal Social Services": formal social care is in,
#                           informal care and lost productivity are not (they
#                           belong to the societal perspective).
#
# A perspective whose REQUIRED built-in columns are absent is skipped, not
# silently degraded -- otherwise a run without socialcare_cost would publish
# healthcare_socialcare tables identical to the healthcare ones.
#
# Discounting: PV = FV / (1 + rate/100)^max(0, year - discount_from_year),
# with separate rates for QALYs and costs. Every discount level in
# `discount_levels` is reported, tagged in a `discount` column. Actual
# (scaled_up) population only.
#
# Output layout. Discount levels and WTP thresholds are CROSSED, and both end up
# on the row axis (the tables are long), via different mechanisms:
#
#   discount -> rows, via rbindlist(): one block per level, in a `discount`
#               column. Applied BEFORE the cumsum, so *_cuml are sums of
#               present values.
#   wtp      -> columns, via the set() loop below, then rows once melt() folds
#               the metric columns into `type`.
#
# Because the wtp loop runs on the already-discount-expanded table, each
# NMB_at_wtp_* column is computed for every level with no nested loop - the
# cross product falls out of the column write. So a stratum cell holds
# length(wtp) x nrow(discount_levels) NMB rows, plus dCosts_cuml, dQALYs_cuml
# and ICER once per level:
#
#   rows = scenarios x strata x nrow(discount_levels) x (3 + length(wtp))
#
# minus the non-finite draws dropped before quantiling - in practice the ICER
# rows where dQALYs_cuml == 0.
#
# A WTP threshold is a ratio of present values (currency per QALY) and so is NOT
# discounted; it multiplies the already-discounted dQALYs_cuml. Under
# differential rates one NMB figure therefore mixes both (QALYs at one rate,
# costs at the other), which is why such a level's label spells both out.
Simulation$set("private", "export_cea_tables", function(
    prbl,
    summaries_dir,
    tables_dir,
    comparator_scenario = "sc0",
    baseline_year = 2019L,
    wtp = c(25000, 35000),
    discount_levels = NULL,
    discount_from_year = NULL,
    custom_costs_in_healthcare = NULL,
    two_agegrps = FALSE,
    two_agegrps_split_age = 65L,
    strata = NULL
) {
  if (self$design$sim_prm$logs) {
    message("Generating cost-effectiveness (ICER/NMB) tables...")
  }

  if (is.null(discount_from_year)) discount_from_year <- baseline_year
  if (is.null(discount_levels)) discount_levels <- build_discount_levels(0, 0)

  qalys <- private$read_summary_dataset("qalys", "scaled_up")
  costs <- private$read_summary_dataset("costs", "scaled_up")
  if (is.null(qalys) || is.null(costs)) {
    if (self$design$sim_prm$logs) {
      message("  qalys or costs summary missing; skipping CEA tables")
    }
    return(invisible(NULL))
  }

  # Summaries store year in short format (e.g. 19); promote to full (2019) so
  # baseline-year filtering and discounting use the same scale as the args.
  qalys[, year := year + 2000L]
  costs[, year := year + 2000L]
  # Allow `qimd` strata even when the summaries only carry `dimd`
  add_qimd_from_dimd(qalys)
  add_qimd_from_dimd(costs)
  # Match export_main_tables(): the CEA strata default to the `ons` list, which
  # may carry `agegrp`. Both frames are aggregated on the same strata and then
  # joined, so both must be collapsed on the same split.
  if (two_agegrps) {
    collapse_two_agegrps(qalys, two_agegrps_split_age)
    collapse_two_agegrps(costs, two_agegrps_split_age)
  }

  # Derived, not passed as a formal (see tbl_smmrs_core for why).
  comparator_is_map <- .comparator_is_map(comparator_scenario)

  # Need at least one complete (intervention, comparator) pair
  if (comparator_is_map) {
    if (length(.resolve_comparators(comparator_scenario,
                                    unique(qalys$scenario))) == 0L) {
      warning("CEA tables: no complete (intervention, comparator) pair in the ",
              "`qalys` summary (found: ",
              paste(unique(qalys$scenario), collapse = ", "),
              "); skipping cost-effectiveness tables.",
              call. = FALSE, immediate. = TRUE)
      return(invisible(NULL))
    }
  } else {
    ## ---- ORIGINAL, VERBATIM ------------------------------------------------
    non_comparator <- setdiff(unique(qalys$scenario), comparator_scenario)
    if (length(non_comparator) == 0L) {
      if (self$design$sim_prm$logs) {
        message("  no intervention scenarios (only '", comparator_scenario,
                "' found); skipping CEA tables")
      }
      return(invisible(NULL))
    }
    ## -------------------------------------------------------------------------
  }

  # QALY scales actually present in the summary (England: EQ5D5L only, but
  # HUI3 is honoured if a future summary provides it).
  scales_avail <- intersect(c("EQ5D5L", "HUI3"), names(qalys))
  if (length(scales_avail) == 0L) {
    if (self$design$sim_prm$logs) {
      message("  no EQ5D5L/HUI3 columns in qalys summary; skipping CEA tables")
    }
    return(invisible(NULL))
  }

  # User-defined cost columns use the plural `_costs$` suffix; the built-in cost
  # columns use the singular `_cost` suffix, so this grep matches only the
  # user-defined columns.
  custom_cost_cols <- grep("_costs$", names(costs), value = TRUE)

  # Resolve which custom cost columns to add to the healthcare perspective.
  # `custom_costs_in_healthcare` accepts a character vector of (custom) cost
  # column names to include there, in addition to the always-present
  # healthcare_cost. For convenience a logical is also honoured:
  # NULL/FALSE -> none (default), TRUE -> all user-defined custom cost columns.
  if (is.null(custom_costs_in_healthcare) ||
      isFALSE(custom_costs_in_healthcare)) {
    healthcare_custom_cols <- character(0)
  } else if (isTRUE(custom_costs_in_healthcare)) {
    healthcare_custom_cols <- custom_cost_cols
  } else {
    requested <- as.character(custom_costs_in_healthcare)
    healthcare_custom_cols <- intersect(requested, custom_cost_cols)
    unknown <- setdiff(requested, custom_cost_cols)
    if (length(unknown) > 0L && self$design$sim_prm$logs) {
      message("  custom_costs_in_healthcare: ignoring name(s) not matching a ",
              "user-defined cost column: ", paste(unknown, collapse = ", "))
    }
  }

  # Cost perspectives, each with the BUILT-IN columns it requires. The required
  # set is checked per perspective rather than by the old
  # "is total_cost or healthcare_cost present" test, because a perspective whose
  # distinguishing column is missing must be SKIPPED, not silently degraded:
  # without the check, a run lacking `socialcare_cost` would publish
  # `healthcare_socialcare` tables byte-identical to the `healthcare` ones under
  # a name claiming to include social care.
  perspective_cfg <- list(
    societal = list(
      cols     = c("total_cost", custom_cost_cols),
      required = "total_cost"),
    healthcare = list(
      cols     = c("healthcare_cost", healthcare_custom_cols),
      required = "healthcare_cost"),
    # NICE's reference case is NHS and Personal Social Services, i.e. direct
    # treatment costs plus formal social care -- but NOT informal care or lost
    # productivity, which sit in the societal perspective only.
    healthcare_socialcare = list(
      cols     = c("healthcare_cost", "socialcare_cost", healthcare_custom_cols),
      required = c("healthcare_cost", "socialcare_cost"))
  )

  # WTP -> NMB column-name labels (plain integer form, e.g. 20000)
  wtp_labels <- paste0(
    "NMB_at_wtp_",
    vapply(wtp, function(w) format(w, scientific = FALSE, trim = TRUE,
                                   big.mark = ""), character(1))
  )

  for (s in strata) {
    x <- c("mc", "scenario", s) # s always contains "year"

    for (persp in names(perspective_cfg)) {
      missing_req <- setdiff(perspective_cfg[[persp]]$required, names(costs))
      if (length(missing_req) > 0L) {
        if (self$design$sim_prm$logs) {
          message("  ", persp, ": required cost column(s) ",
                  paste(missing_req, collapse = ", "), " missing; skipping")
        }
        next
      }
      pcols <- intersect(perspective_cfg[[persp]]$cols, names(costs))

      # Aggregate costs for this perspective, once per stratum. Discounting is
      # deferred until after the incremental difference, so a single aggregate
      # serves every discount level.
      cc <- copy(costs)
      cc[, .cost := Reduce(`+`, .SD), .SDcols = pcols]
      cc <- cc[, .(C = sum(.cost)), keyby = eval(x)]

      for (scale in scales_avail) {
        # Aggregate QALYs for this scale
        qq <- qalys[, .(Q = sum(get(scale))), keyby = eval(x)]

        d <- merge(qq, cc, by = x, all = TRUE)
        d[is.na(Q), Q := 0][is.na(C), C := 0]

        # Incremental vs comparator (intervention - comparator), from baseline.
        # Resolved from `d` -- the MERGED qalys+costs table -- not from `qalys`
        # alone: `costs` can carry a scenario `qalys` lacks, which the old
        # `scenario != comparator_scenario` predicate kept (with Q = 0).
        # `cmp` KEEPS its `scenario` column: it is now the join key, and under a
        # map `d` and `cmp` overlap rather than partitioning.
        cmap <- .resolve_comparators(comparator_scenario, unique(d$scenario))
        cmp <- d[scenario %chin% unique(unname(cmap)) & year >= baseline_year]
        d   <- d[scenario %chin% names(cmap) & year >= baseline_year]
        if (nrow(d) == 0L || nrow(cmp) == 0L) next
        d[, .cmp__ := unname(cmap[as.character(scenario)])]
        d[cmp, on = c(".cmp__" = "scenario", setdiff(x, "scenario")),
          `:=`(dQ = Q - i.Q, dC = C - i.C)]
        set(d, j = ".cmp__", value = NULL)
        d <- d[!is.na(dQ) & !is.na(dC)]
        if (nrow(d) == 0L) next

        # One block of rows per discount level. Unlike the single-outcome
        # tables, a CEA row mixes both rates, so QALYs and costs are each
        # scaled by their own factor before being cumulated. Runs AFTER
        # `.cmp__` is dropped, so the helper column is not replicated.
        d <- rbindlist(lapply(seq_len(nrow(discount_levels)), function(i) {
          di <- copy(d)
          di[, `:=`(
            dQ = dQ * discount_factor(year, discount_levels$qaly[i],
                                      discount_from_year),
            dC = dC * discount_factor(year, discount_levels$cost[i],
                                      discount_from_year),
            discount = discount_levels$label[i]
          )]
          di
        }))
        # `x` itself stays untouched: it still names the filename suffix and
        # the map-mode `comparator` insert point below.
        xd <- c(x, "discount")

        # Cumulative over year within (mc, scenario, other strata, discount)
        setkeyv(d, c(setdiff(xd, "year"), "year"))
        d[, `:=`(dQALYs_cuml = cumsum(dQ), dCosts_cuml = cumsum(dC)),
          by = setdiff(xd, "year")]

        # ICER and NMB at each WTP. The loop runs over the already-expanded
        # table, so every threshold is computed at every discount level.
        d[, ICER := fifelse(dQALYs_cuml == 0, NA_real_,
                            dCosts_cuml / dQALYs_cuml)]
        for (i in seq_along(wtp)) {
          set(d, NULL, wtp_labels[i], wtp[i] * d$dQALYs_cuml - d$dCosts_cuml)
        }

        metric_cols <- c("dCosts_cuml", "dQALYs_cuml", "ICER", wtp_labels)
        # `discount` must be an id.var: this melt names measure.vars, so a
        # column left out of id.vars is silently DROPPED.
        dm <- melt(d, id.vars = xd, measure.vars = metric_cols,
                   variable.name = "type", value.name = "value")
        # Drop non-finite draws (e.g. ICER when dQALYs_cuml == 0) so the
        # quantile is taken over the finite Monte-Carlo iterations only.
        dm <- dm[is.finite(value)]
        if (nrow(dm) == 0L) next

        setkey(dm, "type")
        out <- dm[, safe_fquantile_byid(value, prbl, id = as.character(type),
                                        rounding = FALSE),
                  keyby = eval(setdiff(xd, "mc"))]
        setnames(out, c(setdiff(xd, "mc"), "type",
                        scales::percent(prbl, prefix = "value_")))
        setkeyv(out, c("type", setdiff(xd, "mc")))
        setcolorder(out, setdiff(xd, "mc"))

        if (comparator_is_map) {
          out[, comparator := unname(cmap[as.character(scenario)])]
          nms <- setdiff(names(out), "comparator")
          setcolorder(out, append(nms, "comparator",
                                  after = match("scenario", nms)))
        }
        suffix <- paste(setdiff(x, c("mc", "scenario", "comparator")),
                        collapse = "-")
        fwrite(out, file.path(
          tables_dir,
          paste0("cost-effectiveness by ", suffix,
                 " (", persp, "-", scale, ") (not standardised).csv")
        ))
        rm(d, dm, out, qq)
      }
      rm(cc)
    }
  }

  rm(qalys, costs)
  invisible(NULL)
})


# export_equity_tables ----
# Generate equity slope-index tables (absolute and relative analogues of the
# Slope Index of Inequality / Relative Index of Inequality) for the cumulative
# policy benefits CPP, CYPP, DPP and net QALYs, distributed across DIMD
# deprivation deciles. See calc_equity_slope_indices() for the statistical
# framework and sign convention.
#
# For each metric and each requested output stratum (e.g. "year" or
# c("year", "sex")) the deprivation gradient is fit PER Monte-Carlo iteration
# across the 10 DIMD deciles, and the indices are then quantiled across the
# Monte-Carlo iterations, mirroring every other table. Those percentiles are
# Monte-Carlo percentile intervals covering parameter and synthpop variation --
# NOT frequentist confidence intervals, and not a complete accounting of
# uncertainty: structural uncertainty (the disease model, the effect sizes, and
# the assumption that the gradient is linear in the ridit) sits outside them.
# There is deliberately no within-draw regression standard error; see the
# vignette for why one would double-count and would price misspecification as
# sampling error.
#
# Uses the actual (scaled_up) population only -- "not standardised". This is a
# choice of ESTIMAND, not a claim that age composition is irrelevant: the
# monitoring SII is fitted to age-standardised rates because its target is risk
# AT A GIVEN AGE, whereas the target here is the health gain actually accruing
# to actual people, in whole cases / deaths / QALYs, as in the distributional
# cost-effectiveness convention. (Direct standardisation would in any case
# rescale only the numerator, leaving the decile populations, ridit ranks and
# weights untouched.) For an age-confounding-free view, request an age-stratified
# gradient with strata$equity = list(c("year", "agegrp")); those are age-specific
# SIIs, not a decomposition, and do not sum to the whole-population index.
Simulation$set("private", "export_equity_tables", function(
    prbl,
    summaries_dir,
    tables_dir,
    comparator_scenario = "sc0",
    baseline_year = 2019L,
    ridit_reference = "comparator",
    strata = NULL,
    two_agegrps = FALSE,
    two_agegrps_split_age = 65L
) {
  ridit_reference <- match.arg(ridit_reference, c("comparator", "scenario"))
  # Derived, not passed as a formal (see tbl_smmrs_core for why).
  comparator_is_map <- .comparator_is_map(comparator_scenario)
  if (self$design$sim_prm$logs) {
    message("Generating equity slope-index tables (ridit reference: ",
            ridit_reference, ")...")
  }

  if (is.null(strata)) strata <- list("year", c("year", "sex"))
  validate_equity_strata(strata)
  # One plan == one output table: the output stratification, the deprivation
  # gradient axis it is fit over, and the filename suffix.
  plans <- build_equity_plans(strata)
  wanted_grads <- unique(vapply(plans, `[[`, character(1), "gradient"))

  # Metric -> source dataset, value-column pattern, suffix to strip for the
  # disease label, and benefit direction. "prevented" = comparator - intervention
  # (CPP/CYPP/DPP), "gained" = intervention - comparator (net QALYs).
  metric_cfg <- list(
    cpp       = list(source = "incd",  pattern = "_incd$",   strip = "_incd$",  dir = "prevented"),
    cypp      = list(source = "prvl",  pattern = "_prvl$",   strip = "_prvl$",  dir = "prevented"),
    dpp       = list(source = "mrtl",  pattern = "_mrtl$",   strip = "_mrtl$",  dir = "prevented"),
    net_qalys = list(source = "qalys", pattern = "^EQ5D5L$", strip = NA,        dir = "gained")
  )

  for (metric in names(metric_cfg)) {
    cfg <- metric_cfg[[metric]]

    tt <- private$read_summary_dataset(cfg$source, "scaled_up")
    if (is.null(tt)) {
      if (self$design$sim_prm$logs) {
        message("  ", metric, ": ", cfg$source, " summary missing; skipping")
      }
      next
    }
    # A missing deprivation column means NO equity table at all, so say so
    # unconditionally rather than under `logs` -- an otherwise silent no-op is
    # indistinguishable from the feature having run.
    if (!any(.deprivation_vars %in% names(tt))) {
      warning(sprintf(
        paste0("equity tables: the `%s` summary has no deprivation column ",
               "(looked for %s), so no equity %s tables were written. Add ",
               "`dimd` to `strata_for_output` in the sim_design YAML and re-run ",
               "$export_summaries()."),
        cfg$source, paste(.deprivation_vars, collapse = "/"), metric),
        call. = FALSE, immediate. = TRUE)
      rm(tt); next
    }
    valcols <- grep(cfg$pattern, names(tt), value = TRUE)
    if (length(valcols) == 0L) {
      if (self$design$sim_prm$logs) {
        message("  ", metric, ": no value columns matching ", cfg$pattern, "; skipping")
      }
      rm(tt); next
    }

    # Promote short year (19) to full (2019) to match baseline_year.
    tt[, year := year + 2000L]

    # Match export_main_tables(): under two_agegrps, `agegrp` means the coarse
    # split. Relevant now that agegrp can be an equity stratum.
    if (two_agegrps) collapse_two_agegrps(tt, two_agegrps_split_age)

    # The index is a benefit-vs-comparator contrast, so both the comparator and
    # at least one intervention scenario must be present.
    scns <- unique(tt$scenario)
    cmap <- .resolve_comparators(comparator_scenario, scns)
    if (comparator_is_map) {
      drops <- attr(cmap, "dropped")
      if (length(drops) > 0L) {
        warning(sprintf(paste0(
          "equity tables: comparator pair(s) %s cannot be built from the `%s` ",
          "summary (found: %s); those %s contrasts were skipped."),
          paste(sQuote(drops), collapse = ", "), cfg$source,
          paste(scns, collapse = ", "), metric),
          call. = FALSE, immediate. = TRUE)
      }
      if (length(cmap) == 0L) { rm(tt); next }
    } else {
      ## ---- ORIGINAL, VERBATIM (the warning text is pinned by
      ## ---- inst/tinytest/test_export_equity_tables.R) ---------------------
      if (!comparator_scenario %in% scns) {
        warning(sprintf(
          paste0("equity tables: comparator scenario '%s' is not in the `%s` ",
                 "summary (found: %s), so no equity %s tables were written."),
          comparator_scenario, cfg$source, paste(scns, collapse = ", "), metric),
          call. = FALSE, immediate. = TRUE)
        rm(tt); next
      }
      non_comparator <- setdiff(scns, comparator_scenario)
      if (length(non_comparator) == 0L) {
        if (self$design$sim_prm$logs) {
          message("  ", metric, ": no intervention scenarios (only '",
                  comparator_scenario, "'); skipping")
        }
        rm(tt); next
      }
      ## ---------------------------------------------------------------------
    }

    # Derive `qimd` on demand so a quintile gradient does not require `qimd` in
    # `strata_for_output` (see add_qimd_from_dimd()).
    if ("qimd" %in% wanted_grads && add_qimd_from_dimd(tt) &&
        self$design$sim_prm$logs) {
      message("  ", metric, ": derived qimd quintiles from dimd deciles")
    }

    for (p in plans) {                        # out_vars always contains "year"
      grad <- p$gradient
      if (!grad %in% names(tt)) {
        warning(sprintf(
          paste0("equity tables: gradient '%s' was requested but the `%s` ",
                 "summary has no `%s` column and it cannot be derived from the ",
                 "columns present (%s); skipping 'equity %s slope index by %s'."),
          grad, cfg$source, grad,
          paste(intersect(.deprivation_vars, names(tt)), collapse = ", "),
          metric, p$suffix),
          call. = FALSE, immediate. = TRUE)
        next
      }
      out_vars <- p$out_vars                    # always contains "year"
      missing_vars <- setdiff(out_vars, names(tt))
      if (length(missing_vars) > 0L) {
        warning(sprintf(
          paste0("equity tables: stratification variable(s) %s are not in the ",
                 "`%s` summary (available: %s); skipping 'equity %s slope index ",
                 "by %s'. Add them to `strata_for_output` in the sim_design YAML ",
                 "and re-run $export_summaries()."),
          paste(sQuote(missing_vars), collapse = ", "), cfg$source,
          paste(setdiff(names(tt), c(valcols, "popsize")), collapse = ", "),
          metric, p$suffix),
          call. = FALSE, immediate. = TRUE)
        next
      }
      keep <- c("mc", "scenario", out_vars, grad)

      # Collapse every summary variable NOT named in this stratum (agegrp, sex,
      # ...); sum the disease value columns and the deprivation-group
      # population. When `grad` is `qimd` this same aggregation is what
      # collapses the deciles.
      agg <- tt[, c(lapply(.SD, sum), .(popsize = sum(popsize))),
                .SDcols = valcols, by = keep]

      long <- melt(agg, id.vars = c(keep, "popsize"),
                   measure.vars = valcols, variable.name = "disease",
                   value.name = "val")
      long[, disease := as.character(disease)]
      if (!is.na(cfg$strip)) long[, disease := gsub(cfg$strip, "", disease)]

      # Benefit versus comparator, from the baseline year onwards.
      idcols <- setdiff(c(keep, "disease"), "scenario")   # mc, <out_vars>, <grad>, disease
      # `cmp` now KEEPS `scenario` -- it is the join key. It is taken from the
      # FULL `long`, so a scenario that is both an intervention and somebody
      # else's comparator appears in BOTH `d` and `cmp` (`[` copies, so
      # mutating one cannot corrupt the other).
      cmp <- long[scenario %chin% unique(unname(cmap)) & year >= baseline_year]
      setnames(cmp, c("val", "popsize"), c("val_cmp", "pop_cmp"))
      d <- long[scenario %chin% names(cmap) & year >= baseline_year]
      if (nrow(d) == 0L || nrow(cmp) == 0L) next
      d[, .cmp__ := unname(cmap[as.character(scenario)])]
      on_cmp <- c(".cmp__" = "scenario", idcols)

      if (cfg$dir == "prevented") {
        d[cmp, on = on_cmp, `:=`(B_year = i.val_cmp - val, N_cmp = i.pop_cmp)]
      } else {
        d[cmp, on = on_cmp, `:=`(B_year = val - i.val_cmp, N_cmp = i.pop_cmp)]
      }
      set(d, j = ".cmp__", value = NULL)
      # Reference population for the ridit ranking, the population weights of the
      # weighted regression and the per-capita denominator. "comparator" uses
      # the comparator scenario's population so the ridit scores and weights are
      # IDENTICAL across the scenarios being compared (recommended: it prevents
      # intervention-induced shifts in the socioeconomic composition from being
      # read as changes in inequality, per Renard et al. 2019, and it is what
      # makes SII(benefit) exactly equal the policy-induced change in the SII).
      #
      # They are deliberately NOT pinned across YEARS: `cmp` is joined on
      # `idcols`, which includes `year`, so N is the comparator population OF THE
      # REPORTING YEAR and the ridit is re-derived each year (the ordinal decile
      # ordering is of course fixed throughout). Pinning N at `baseline_year` was
      # considered and rejected -- it would divide a benefit generated by the
      # year-t population by the year-0 population, so a benefit that is exactly
      # equal per head in year t would score a spurious non-zero index purely
      # from differential population growth between deciles.
      #
      # "scenario" uses each intervention scenario's own population (the
      # composition it actually produces); `popsize` here is that population.
      # It violates the precondition Cookson et al. 2026 state for the impact
      # SII (equal group sizes across arms), so treat it as a sensitivity
      # analysis rather than a co-equal estimator.
      d[, N := if (ridit_reference == "comparator") N_cmp else popsize]
      d <- d[is.finite(B_year) & is.finite(N) & N > 0]
      if (nrow(d) == 0L) next

      # Cumulative benefit over year within
      # (mc, scenario, <out_vars except year>, <grad>, disease).
      by_cum <- c("mc", "scenario", setdiff(out_vars, "year"), grad, "disease")
      setkeyv(d, c(by_cum, "year"))
      d[, B := cumsum(B_year), by = by_cum]

      # All-diseases summed row for CPP/CYPP (plain event-count sum; carries
      # comorbidity multiplicity -- not unique persons). DPP is already
      # all-cause; net QALYs is already a total. Two families are held out of
      # the sum: the CMS multimorbidity indicators (cms1st_cont, cmsmm*, cmscs*),
      # which are derived metrics rather than incident diseases; and the design's
      # umbrella conditions (incidence `type: 0`, e.g. `cancer` alongside the
      # five site-specific cancers), which would count the same event twice --
      # see umbrella_disease_names(). Both are kept as their own rows.
      if (metric %in% c("cpp", "cypp")) {
        by_tot <- c("mc", "scenario", out_vars, grad)
        umbrella <- umbrella_disease_names(self$design$sim_prm$diseases,
                                           unique(d$disease))
        if (length(umbrella) > 0L && self$design$sim_prm$logs) {
          message("  ", metric, ": excluding umbrella condition(s) from ",
                  "all_diseases_sum (components are summed instead): ",
                  paste(sort(umbrella), collapse = ", "))
        }
        tot <- d[!grepl("^cms", disease) & !disease %chin% umbrella,
                 .(disease = "all_diseases_sum", B = sum(B),
                   B_year = sum(B_year), N = N[1L]), by = by_tot]
        d <- rbind(d, tot, fill = TRUE)
      }

      set(d, j = "rank", value = deprivation_rank(d[[grad]], grad))

      # A deprivation gradient needs at least two deprivation groups PRESENT IN
      # THE SUMMARIES. This counts distinct levels present, so it fires only for
      # a genuinely single-group extract (e.g. a single-decile locality), whose
      # slope index is undefined and is dropped at the finite-value filter.
      # NOTE it does NOT fire for an intervention TARGETED at one decile: the
      # untouched deciles contribute B_year = 0, which is finite, so all ten
      # groups survive the filter above and a fully extrapolated index is
      # produced. That is intended -- but see the vignette, which warns that the
      # indices are then extreme-end extrapolations from a line fit to a spike.
      # Report the omission so it is transparent rather than silent.
      if (self$design$sim_prm$logs) {
        cov_grp <- c("scenario", out_vars)
        ndec <- unique(d[, c(cov_grp, grad), with = FALSE])[, .(nd = .N),
                                                            by = cov_grp]
        insuff <- ndec[nd < 2L]
        if (nrow(insuff) > 0L) {
          message("  ", metric, " by ", p$suffix, ": omitting ",
                  nrow(insuff), " output group(s) with <2 ", grad,
                  " groups (undefined gradient); scenario(s): ",
                  paste(sort(unique(insuff$scenario)), collapse = ", "))
        }
      }

      by_grp <- c("mc", "scenario", out_vars, "disease")
      idx <- calc_equity_slope_indices(
        d[, .SD, .SDcols = c(by_grp, "rank", "B", "N")], by = by_grp)

      m <- melt(idx, id.vars = by_grp,
                measure.vars = c("AEI_total", "AEI_per100k", "total_benefit",
                                 "REI_rel", "RII_ratio", "fit_R2"),
                variable.name = "type", value.name = "value")
      # Quantile over the finite Monte-Carlo draws only (REI_rel / RII_ratio /
      # fit_R2 are NA on degenerate or non-positive gradients). The filter is
      # deliberately PER TYPE: dropping a whole draw because one derived ratio
      # is undefined would import that ratio's definedness selection into
      # AEI_total, the always-defined statistic the literature treats as primary
      # (King, Harper & Young 2012). The cost is that different types can rest
      # on different numbers of draws, so `n_mc` below records how many.
      m <- m[is.finite(value)]
      if (nrow(m) == 0L) next

      setkey(m, "type")
      out <- m[, safe_fquantile_byid(value, prbl, id = as.character(type),
                                     rounding = FALSE),
               keyby = eval(setdiff(by_grp, "mc"))]
      setnames(out, c(setdiff(by_grp, "mc"), "type",
                      scales::percent(prbl, prefix = "equity_")))
      setkeyv(out, c(setdiff(by_grp, "mc"), "type"))

      # How many Monte-Carlo draws each published quantile actually rests on.
      # This is NOT the number of MC iterations run: it counts draws where the
      # group existed AND that type was defined. A low n_mc on a relative index
      # means it has been truncated by the positivity guards and should not be
      # quoted on its own. It does not flag the other failure mode -- a
      # near-zero mean benefit yields a full-n_mc but wildly unstable REI_rel.
      nmc <- m[, .(n_mc = .N), keyby = c(setdiff(by_grp, "mc"), "type")]
      nmc[, type := as.character(type)]   # melt yields a factor; the CSV is character
      out[nmc, on = c(setdiff(by_grp, "mc"), "type"), n_mc := i.n_mc]
      # The filename suffix reads as "stratified by", but the gradient token in
      # it means "gradient over"; the `gradient` column makes each file
      # self-describing regardless of how it was named.
      out[, gradient := grad]
      if (comparator_is_map) {
        # Sits next to `gradient` -- the same design move for the same reason:
        # the file must be self-describing about which contrast each row is.
        out[, comparator := unname(cmap[as.character(scenario)])]
        setcolorder(out, c("gradient", "scenario", "comparator",
                           setdiff(by_grp, c("mc", "scenario")), "type"))
      } else {
        setcolorder(out, c("gradient", setdiff(by_grp, "mc"), "type"))
      }

      fwrite(out, file.path(
        tables_dir,
        paste0("equity ", metric, " slope index by ", p$suffix,
               " (not standardised).csv")
      ))
      rm(agg, long, cmp, d, idx, m, out)
    }

    rm(tt)
    gc(verbose = FALSE)
  }

  invisible(NULL)
})
