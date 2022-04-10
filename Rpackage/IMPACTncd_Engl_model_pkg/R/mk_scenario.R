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
mk_scenario_init <- function(scenario_name, design, lags_mc) {
  # scenario_suffix_for_pop <- paste0("_", scenario_name) # TODO get suffix from design
  scenario_suffix_for_pop <- paste0(scenario_name)
  list(
    "exposures"          = design$exposures,
    "scenarios"          = design$scenarios, # to be generated programmatically
    "scenario"           = scenario_name,
    "kismet"             = design$kismet, # If TRUE random numbers are the same for each scenario.
    "init_year"          = design$init_year,
    "pids"               = list("column_name" = "pid"),
    "years"              = list("column_name" = "year"),
    "rn_all_cause_mrtl"  = list("column_name" = paste0("rn_multi_mrtl", scenario_suffix_for_pop)),
    "strata_for_outputs" = c("pid", "year", "age", "sex", "qimd"),
    "diseases" = list(
      # One element per disease to be included in the model.
      # Output columns will be named with the disease name as a prefix, and the kind of result as the suffix
      "htn" = list(
        "incidence" = mk_incidence_type1("htn", scenario_suffix_for_pop, design$kismet),
        "diagnosis" = mk_diagnosis_type1("htn", scenario_suffix_for_pop, design$kismet)
      ),
      "af" = list(
        "incidence" = mk_incidence_type1("af", scenario_suffix_for_pop, design$kismet),
        "diagnosis" = mk_diagnosis_type1("af", scenario_suffix_for_pop, design$kismet)
      ),
      "t2dm" = list(
        "incidence" = mk_incidence_type3_specific_column("t2dm", scenario_suffix_for_pop, "cvd", lags_mc$cvd_lag, "t2dm_incd_cvd_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("t2dm", scenario_suffix_for_pop, design$kismet)
        # t2dm never kills directly, so there's no mortality
      ),
      "chd" = list(
        "incidence" = mk_incidence_type3_specific_column("chd", scenario_suffix_for_pop, "t2dm", lags_mc$cvd_lag, "chd_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("chd", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type1("chd", scenario_suffix_for_pop, 2L, design$kismet)
      ),
      "stroke" = list(
        "incidence" = mk_incidence_type3_specific_column("stroke", scenario_suffix_for_pop, "t2dm", lags_mc$cvd_lag, "stroke_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("stroke", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type1("stroke", scenario_suffix_for_pop, 3L, design$kismet)
      ),
      "poststroke_dementia" = list(
        "incidence" = mk_incidence_type5_specific_column("poststroke_dementia", scenario_suffix_for_pop, "stroke", 2L, "poststroke_dementia_incd_stroke_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("stroke", scenario_suffix_for_pop, design$kismet)
      ),
      "cvd" = list(
        "incidence" = mk_incidence_type4("cvd", scenario_suffix_for_pop, c("chd", "stroke"), design$kismet)
        # CVD is essentially a label for "CHD or stroke", so there's no independent diagnosis or mortality.
      ),
      "copd" = list(
        "incidence" = mk_incidence_type2("copd", FALSE, scenario_suffix_for_pop, design$kismet),
        "diagnosis" = mk_diagnosis_type1("copd", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type1("copd", scenario_suffix_for_pop, 4L, design$kismet)
      ),
      "lung_ca" = list(
        "incidence" = mk_incidence_type3_specific_column("lung_ca", scenario_suffix_for_pop,
                                                         "copd", lags_mc$plco_lag, "lung_ca_incd_copd_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("lung_ca", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type2("lung_ca", scenario_suffix_for_pop, 5L, design$cancer_cure, design$kismet)
      ),
      "colon_ca" = list(
        "incidence" = mk_incidence_type3_specific_column("colon_ca", scenario_suffix_for_pop,
                                                         "t2dm", lags_mc$cancer_lag, "colon_ca_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("colon_ca", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type2("colon_ca", scenario_suffix_for_pop, 6L, design$cancer_cure, design$kismet)
      ),
      "breast_ca" = list(
        "incidence" = mk_incidence_type3_specific_column("breast_ca", scenario_suffix_for_pop,
                                                         "t2dm", lags_mc$cancer_lag, "breast_ca_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("breast_ca", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type2("breast_ca", scenario_suffix_for_pop, 7L, design$cancer_cure, design$kismet)
      ),
      "nonmodelled" = list(
        # Mortality from non-modelled causes.
        # Everybody suffers from being alive until they're dead.
        # This approach means we can use the normal disease model for "everything else" rather than custom code.
        "incidence" = list("type" = "Universal"),
        "mortality" = mk_mortality_type3_specific_column("nonmodelled", scenario_suffix_for_pop, 1L, "t2dm", lags_mc$nonmodelled_lag, "nonmodelled_mrtl_t2dm_mltp", design$kismet)
      )
    )
  )
}

#' @export
mk_scenario <- function(scenario_name, design, lags_mc, POP) {
  scenario_suffix_for_pop <-
    fifelse(nzchar(scenario_name), "_sc", scenario_name)
  list(
    "exposures"          = design$exposures,
    "scenarios"          = design$scenarios, # to be generated programmatically
    "scenario"           = scenario_name,
    "kismet"             = design$kismet, # If TRUE random numbers are the same for each scenario.
    "init_year"          = design$init_year_fromGUI,
    "pids"               = list("column_name" = "pid"),
    "years"              = list("column_name" = "year"),
    "rn_all_cause_mrtl"  = list("column_name" =
                                  paste0("rn_multi_mrtl",
                                         fifelse(design$kismet, "", scenario_suffix_for_pop))),
    "strata_for_outputs" = c(design$cols_for_output,
                             grep(paste0("cost", scenario_suffix_for_pop),
                                  names(POP), value = TRUE)),
    "diseases" = list(
      # One element per disease to be included in the model.
      # Output columns will be named with the disease name as a prefix, and the kind of result as the suffix
      "htn" = list(
        "incidence" = mk_incidence_type1("htn", scenario_suffix_for_pop, design$kismet),
        "diagnosis" = mk_diagnosis_type1("htn", scenario_suffix_for_pop, design$kismet)
      ),
      "af" = list(
        "incidence" = mk_incidence_type1("af", scenario_suffix_for_pop, design$kismet),
        "diagnosis" = mk_diagnosis_type1("af", scenario_suffix_for_pop, design$kismet)
      ),
      "t2dm" = list(
        "incidence" = mk_incidence_type3_specific_column("t2dm", scenario_suffix_for_pop, "cvd", lags_mc$cvd_lag, "t2dm_incd_cvd_mltp_sc", design$kismet),
        "diagnosis" = mk_diagnosis_type1("t2dm", scenario_suffix_for_pop, design$kismet)
        # t2dm never kills directly, so there's no mortality
      ),
      "chd" = list(
        "incidence" = mk_incidence_type3_specific_column("chd", scenario_suffix_for_pop, "t2dm", lags_mc$cvd_lag, "chd_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("chd", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type1("chd", scenario_suffix_for_pop, 2L, design$kismet)
      ),
      "stroke" = list(
        "incidence" = mk_incidence_type3_specific_column("stroke", scenario_suffix_for_pop, "t2dm", lags_mc$cvd_lag, "stroke_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("stroke", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type1("stroke", scenario_suffix_for_pop, 3L, design$kismet)
      ),
      "poststroke_dementia" = list(
        "incidence" = mk_incidence_type5_specific_column("poststroke_dementia", scenario_suffix_for_pop, "stroke", 2L, "poststroke_dementia_incd_stroke_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("stroke", scenario_suffix_for_pop, design$kismet)
      ),
      "cvd" = list(
        "incidence" = mk_incidence_type4("cvd", scenario_suffix_for_pop, c("chd", "stroke"), design$kismet)
        # CVD is essentially a label for "CHD or stroke", so there's no independent diagnosis or mortality.
      ),
      "copd" = list(
        "incidence" = mk_incidence_type2("copd", FALSE, scenario_suffix_for_pop, design$kismet),
        "diagnosis" = mk_diagnosis_type1("copd", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type1("copd", scenario_suffix_for_pop, 4L, design$kismet)
      ),
      "lung_ca" = list(
        "incidence" = mk_incidence_type3_specific_column("lung_ca", scenario_suffix_for_pop,
                                                         "copd", lags_mc$plco_lag, "lung_ca_incd_copd_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("lung_ca", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type2("lung_ca", scenario_suffix_for_pop, 5L, design$cancer_cure, design$kismet)
      ),
      "colon_ca" = list(
        "incidence" = mk_incidence_type3_specific_column("colon_ca", scenario_suffix_for_pop,
                                                         "t2dm", lags_mc$cancer_lag, "colon_ca_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("colon_ca", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type2("colon_ca", scenario_suffix_for_pop, 6L, design$cancer_cure, design$kismet)
      ),
      "breast_ca" = list(
        "incidence" = mk_incidence_type3_specific_column("breast_ca", scenario_suffix_for_pop,
                                                         "t2dm", lags_mc$cancer_lag, "breast_ca_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("breast_ca", scenario_suffix_for_pop, design$kismet),
        "mortality" = mk_mortality_type2("breast_ca", scenario_suffix_for_pop, 7L, design$cancer_cure, design$kismet)
      ),
      "nonmodelled" = list(
        # Mortality from non-modelled causes.
        # Everybody suffers from being alive until they're dead.
        # This approach means we can use the normal disease model for "everything else" rather than custom code.
        "incidence" = list("type" = "Universal"),
        "mortality" = mk_mortality_type3_specific_column("nonmodelled", scenario_suffix_for_pop, 1L, "t2dm", lags_mc$nonmodelled_lag, "nonmodelled_mrtl_t2dm_mltp", design$kismet)
      )
    )
  )
}

#' @export
mk_scenario_mortality_calib <- function(scenario_name, design, lags_mc) {
  # scenario_suffix_for_pop <- paste0("_", scenario_name) # TODO get suffix from design
  scenario_suffix_for_pop <- paste0(scenario_name)
  list(
    "exposures"          = design$exposures,
    "scenarios"          = design$scenarios, # to be generated programmatically
    "scenario"           = scenario_name,
    "kismet"             = design$kismet, # If TRUE random numbers are the same for each scenario.
    "init_year"          = design$init_year,
    "pids"               = list("column_name" = "pid"),
    "years"              = list("column_name" = "year"),
    "rn_all_cause_mrtl"  = list("column_name" = paste0("rn_multi_mrtl", scenario_suffix_for_pop)),
    "strata_for_outputs" = c("pid", "year", "age", "sex", "qimd"),
    "diseases" = list(
      # One element per disease to be included in the model.
      # Output columns will be named with the disease name as a prefix, and the kind of result as the suffix
      "htn" = list(
        "incidence" = mk_incidence_type1("htn", scenario_suffix_for_pop, design$kismet),
        "diagnosis" = mk_diagnosis_type1("htn", scenario_suffix_for_pop, design$kismet)
      ),
      "af" = list(
        "incidence" = mk_incidence_type1("af", scenario_suffix_for_pop, design$kismet),
        "diagnosis" = mk_diagnosis_type1("af", scenario_suffix_for_pop, design$kismet)
      ),
      "t2dm" = list(
        "incidence" = mk_incidence_type3_specific_column("t2dm", scenario_suffix_for_pop, "cvd", lags_mc$cvd_lag, "t2dm_incd_cvd_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("t2dm", scenario_suffix_for_pop, design$kismet)
        # t2dm never kills directly, so there's no mortality
      ),
      "chd" = list(
        "incidence" = mk_incidence_type3_specific_column("chd", scenario_suffix_for_pop, "t2dm", lags_mc$cvd_lag, "chd_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("chd", scenario_suffix_for_pop, design$kismet)
      ),
      "stroke" = list(
        "incidence" = mk_incidence_type3_specific_column("stroke", scenario_suffix_for_pop, "t2dm", lags_mc$cvd_lag, "stroke_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("stroke", scenario_suffix_for_pop, design$kismet)
      ),
      "poststroke_dementia" = list(
        "incidence" = mk_incidence_type5_specific_column("poststroke_dementia", scenario_suffix_for_pop, "stroke", 2L, "poststroke_dementia_incd_stroke_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("stroke", scenario_suffix_for_pop, design$kismet)
      ),
      "cvd" = list(
        "incidence" = mk_incidence_type4("cvd", scenario_suffix_for_pop, c("chd", "stroke"), design$kismet)
        # CVD is essentially a label for "CHD or stroke", so there's no independent diagnosis or mortality.
      ),
      "copd" = list(
        "incidence" = mk_incidence_type2("copd", FALSE, scenario_suffix_for_pop, design$kismet),
        "diagnosis" = mk_diagnosis_type1("copd", scenario_suffix_for_pop, design$kismet)
      ),
      "lung_ca" = list(
        "incidence" = mk_incidence_type3_specific_column("lung_ca", scenario_suffix_for_pop,
                                                         "copd", lags_mc$plco_lag, "lung_ca_incd_copd_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("lung_ca", scenario_suffix_for_pop, design$kismet)
      ),
      "colon_ca" = list(
        "incidence" = mk_incidence_type3_specific_column("colon_ca", scenario_suffix_for_pop,
                                                         "t2dm", lags_mc$cancer_lag, "colon_ca_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("colon_ca", scenario_suffix_for_pop, design$kismet)
      ),
      "breast_ca" = list(
        "incidence" = mk_incidence_type3_specific_column("breast_ca", scenario_suffix_for_pop,
                                                         "t2dm", lags_mc$cancer_lag, "breast_ca_incd_t2dm_mltp", design$kismet),
        "diagnosis" = mk_diagnosis_type1("breast_ca", scenario_suffix_for_pop, design$kismet)
      )
    )
  )
}

#' @export
gen_output <- function(scenario_nam, design, lags_mc, dt,
                       output, calibration_run = FALSE) {
  if (calibration_run) {
    run_impactncd_simulation(
      mk_scenario_mortality_calib(scenario_nam, design, lags_mc),
      dt,
      output)
  } else {
    if (!nzchar(scenario_nam)) {
      message("baseline simulation")
      run_impactncd_simulation(
        mk_scenario_init(scenario_nam, design, lags_mc),
        dt,
        output)
    } else {
      message("scenario simulation")
      run_impactncd_simulation(
        mk_scenario(scenario_nam, design, lags_mc, dt),
        dt,
        output)
    }
  }
}




# The functions below have been developed by Peter Crowther from Melandra Ltd
# Convenience functions to generate definitions - these encode conventions on how columns will be named.
# TODO: Use these to generate shorthands for entire diseases.

#' @export
mk_incidence_type1 <- function(disease_name, scenario_suffix, kismet) {
  # kismet is not used. Just for consistency
  list(
    "type" = "Type1",
    "prevalence" = list("column_name" = paste0(disease_name, "_prvl", scenario_suffix)),
    # "prevalence" = list("column_name" = paste0(disease_name, "_prvl")),

    "incidence" = list("column_name" = paste0("prb_", disease_name, "_incd", scenario_suffix))
  )
}

#' @export
mk_incidence_type2 <- function(disease_name, can_recur, scenario_suffix, kismet) {
  list(
    "type" = "Type2",
    "prevalence" = list("column_name" = paste0(disease_name, "_prvl", scenario_suffix)),
    # "prevalence" = list("column_name" = paste0(disease_name, "_prvl")),

    "incidence" = list("column_name" = paste0("prb_", disease_name, "_incd", scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_incd",
                                       fifelse(kismet, "", scenario_suffix))),
    "can_recur" = can_recur
  )
}

#' @export
mk_incidence_type3 <- function(disease_name, scenario_suffix, influenced_by_disease_name, lag, kismet) {
  influenced_by <- list()
  influenced_by[[influenced_by_disease_name]] <-
    list(
      "multiplier" = list("column_name" = paste0(disease_name, "_", influenced_by_disease_name, "_incd_mltp")),
      "lag" = lag # Ensure from R side that never 0
    )
  list(
    "type" = "Type3",
    "prevalence" = list("column_name" = paste0(disease_name, "_prvl", scenario_suffix)),
    # "prevalence" = list("column_name" = paste0(disease_name, "_prvl")),

    "incidence" = list("column_name" = paste0("prb_", disease_name, "_incd", scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_incd",
                                       fifelse(kismet, "", scenario_suffix))),
    "influenced_by" = influenced_by
  )
}

#' @export
mk_incidence_type3_specific_column <- function(disease_name, scenario_suffix, influenced_by_disease_name, lag, column_name, kismet) {
  influenced_by <- list()
  influenced_by[[influenced_by_disease_name]] <-
    list(
      "multiplier" = list("column_name" = column_name),
      "lag" = lag # Ensure from R side that never 0
    )
  list(
    "type" = "Type3",
    "prevalence" = list("column_name" = paste0(disease_name, "_prvl", scenario_suffix)),
    # "prevalence" = list("column_name" = paste0(disease_name, "_prvl")),

    "incidence" = list("column_name" = paste0("prb_", disease_name, "_incd_no",
                                              influenced_by_disease_name, scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_incd",
                                       fifelse(kismet, "", scenario_suffix))),
    "influenced_by" = influenced_by
  )
}

#' @export
mk_incidence_type4 <- function(disease_name, scenario_suffix, aggregates, kismet) {
  list(
    "type" = "Type4",
    "aggregates" = aggregates
  )
}

# like type 3 but multiplier applies only for incident cases, not prevalent
# Lag needs to be 2L. For some reason it doesn't work with lag == 1L
# To be used in post-stroke dementia
#' @export
mk_incidence_type5 <- function(disease_name, scenario_suffix, influenced_by_disease_name, lag = 1L, kismet) {
  influenced_by <- list()
  influenced_by[[influenced_by_disease_name]] <-
    list(
      "multiplier" = list("column_name" = paste0(disease_name, "_", influenced_by_disease_name, "_incd_mltp")),
      "lag" = lag # Ensure from R side that never 0
    )
  list(
    "type" = "Type5",
    "prevalence" = list("column_name" = paste0(disease_name, "_prvl", scenario_suffix)),
    # "prevalence" = list("column_name" = paste0(disease_name, "_prvl")),

    "incidence" = list("column_name" = paste0("prb_", disease_name, "_incd", scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_incd",
                                       fifelse(kismet, "", scenario_suffix))),
    "influenced_by" = influenced_by
  )
}

#' @export
mk_incidence_type5_specific_column <- function(disease_name, scenario_suffix, influenced_by_disease_name, lag, column_name, kismet) {
  influenced_by <- list()
  influenced_by[[influenced_by_disease_name]] <-
    list(
      "multiplier" = list("column_name" = column_name),
      "lag" = lag # Ensure from R side that never 0
    )
  list(
    "type" = "Type5",
    "prevalence" = list("column_name" = paste0(disease_name, "_prvl", scenario_suffix)),
    # "prevalence" = list("column_name" = paste0(disease_name, "_prvl")),

    "incidence" = list("column_name" = paste0("prb_", disease_name, "_incd_no",
                                              influenced_by_disease_name, scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_incd",
                                       fifelse(kismet, "", scenario_suffix))),
    "influenced_by" = influenced_by
  )
}

#' @export
mk_diagnosis_type1 <- function(disease_name, scenario_suffix, kismet) {
  list(
    "type" = "Type1",
    "diagnosed" = list("column_name" = paste0(disease_name, "_dgn", scenario_suffix)),
    "probability" = list("column_name" = paste0("prb_", disease_name, "_dgn", scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_dgn",
                                       fifelse(kismet, "", scenario_suffix)))
  )
}



#' @export
mk_mortality_type1 <- function(disease_name, scenario_suffix, code, kismet) {
  list(
    "type" = "Type1",
    "probability" = list("column_name" = paste0("prb_", disease_name, "_mrtl", scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_mrtl", fifelse(kismet, "", scenario_suffix))),
    "code" = code
  )
}

#' @export
mk_mortality_type2 <- function(disease_name, scenario_suffix, code, cure, kismet) {
  list(
    "type" = "Type2",
    "probability" = list("column_name" = paste0("prb_", disease_name, "_mrtl", scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_mrtl", fifelse(kismet, "", scenario_suffix))),
    "cure" = cure,
    "code" = code
  )
}

#' @export
mk_mortality_type3 <- function(disease_name, scenario_suffix, code, influenced_by_disease_name, lag, kismet) {
  influenced_by <- list()
  influenced_by[[influenced_by_disease_name]] <-
    list(
      "multiplier" = list("column_name" = paste0(disease_name, "_", influenced_by_disease_name, "_mrtl_mltp")),
      "lag" = lag # Ensure from R side that never 0
    )
  list(
    "type" = "Type3",
    "probability" = list("column_name" = paste0("prb_", disease_name, "_mrtl", scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_mrtl", fifelse(kismet, "", scenario_suffix))),
    "influenced_by" = influenced_by,
    "code" = code
  )
}

#' @export
mk_mortality_type3_specific_column <- function(disease_name, scenario_suffix, code, influenced_by_disease_name, lag, column_name, kismet) {
  influenced_by <- list()
  influenced_by[[influenced_by_disease_name]] <-
    list(
      "multiplier" = list("column_name" = column_name),
      "lag" = lag # Ensure from R side that never 0
    )
  list(
    "type" = "Type3",
    "probability" = list("column_name" = paste0("prb_", disease_name, "_mrtl_no", influenced_by_disease_name, scenario_suffix)),
    "rn" = list("column_name" = paste0("rn_", disease_name, "_mrtl", fifelse(kismet, "", scenario_suffix))),
    "influenced_by" = influenced_by,
    "code" = code
  )
}

