# Graph Report - .  (2026-04-30)

## Corpus Check
- 210 files · ~449,015 words
- Verdict: corpus is large enough that graph structure adds value.

## Summary
- 725 nodes · 1198 edges · 44 communities detected
- Extraction: 86% EXTRACTED · 14% INFERRED · 0% AMBIGUOUS · INFERRED: 169 edges (avg confidence: 0.79)
- Token cost: 0 input · 0 output

## Community Hubs (Navigation)
- [[_COMMUNITY_Cardiometabolic Risk Fitting (T2DMHDLTChol)|Cardiometabolic Risk Fitting (T2DM/HDL/TChol)]]
- [[_COMMUNITY_CPRD Validation & HFSS Energy Models|CPRD Validation & HFSS Energy Models]]
- [[_COMMUNITY_Core Architecture & Bootstrap|Core Architecture & Bootstrap]]
- [[_COMMUNITY_Simulation Testing Scripts|Simulation Testing Scripts]]
- [[_COMMUNITY_HSE Data Preparation|HSE Data Preparation]]
- [[_COMMUNITY_Sociodemographic Distribution Fits|Sociodemographic Distribution Fits]]
- [[_COMMUNITY_Disease Epi Helpers (panel data)|Disease Epi Helpers (panel data)]]
- [[_COMMUNITY_Zenodo Asset Management|Zenodo Asset Management]]
- [[_COMMUNITY_Population Estimates & SynthPop Seeding|Population Estimates & SynthPop Seeding]]
- [[_COMMUNITY_Cancer Epidemiology Modeling|Cancer Epidemiology Modeling]]
- [[_COMMUNITY_C++ Simulation Engine (top-level)|C++ Simulation Engine (top-level)]]
- [[_COMMUNITY_SynthPop R6 Class Methods|SynthPop R6 Class Methods]]
- [[_COMMUNITY_Simulation Export & Summaries|Simulation Export & Summaries]]
- [[_COMMUNITY_C++ Disease Incidence Core|C++ Disease Incidence Core]]
- [[_COMMUNITY_Health Economics Generators|Health Economics Generators]]
- [[_COMMUNITY_C++ Disease Incidence Types|C++ Disease Incidence Types]]
- [[_COMMUNITY_Reproducibility Testing|Reproducibility Testing]]
- [[_COMMUNITY_Cancer History Inputs (PLCO)|Cancer History Inputs (PLCO)]]
- [[_COMMUNITY_Zenodo Upload Tests|Zenodo Upload Tests]]
- [[_COMMUNITY_VSCode Debug Env Capture|VSCode Debug Env Capture]]
- [[_COMMUNITY_Compression Helpers (deflateinflate)|Compression Helpers (deflate/inflate)]]
- [[_COMMUNITY_GAMLSS Prediction Helpers|GAMLSS Prediction Helpers]]
- [[_COMMUNITY_ONS IncidenceMortality Fetchers|ONS Incidence/Mortality Fetchers]]
- [[_COMMUNITY_Smoking Cessation Sim|Smoking Cessation Sim]]
- [[_COMMUNITY_Smoking Policy Impact Sim|Smoking Policy Impact Sim]]
- [[_COMMUNITY_Test Framework (tinytest)|Test Framework (tinytest)]]
- [[_COMMUNITY_SynthPop Validation Plots|SynthPop Validation Plots]]
- [[_COMMUNITY_Project Critical Rules|Project Critical Rules]]
- [[_COMMUNITY_R6 Hierarchy Doc|R6 Hierarchy Doc]]
- [[_COMMUNITY_Data Flow Diagram|Data Flow Diagram]]
- [[_COMMUNITY_Correlation Plot Helper|Correlation Plot Helper]]
- [[_COMMUNITY_Distribution Validation Helper|Distribution Validation Helper]]
- [[_COMMUNITY_RNG Generation Helper|RNG Generation Helper]]
- [[_COMMUNITY_Output Directory Helper|Output Directory Helper]]
- [[_COMMUNITY_Smoking Cigarette Sim|Smoking Cigarette Sim]]
- [[_COMMUNITY_Smoking Cigarette Scenario Sim|Smoking Cigarette Scenario Sim]]
- [[_COMMUNITY_Smoking Cessation Cigarette Sim|Smoking Cessation Cigarette Sim]]
- [[_COMMUNITY_Rcpp Auto-Generated Bindings|Rcpp Auto-Generated Bindings]]
- [[_COMMUNITY_C++ Fortify Header (2)|C++ Fortify Header (2)]]
- [[_COMMUNITY_Stratification (AgeSexDIMD)|Stratification (Age/Sex/DIMD)]]
- [[_COMMUNITY_Disease Burden Estimation|Disease Burden Estimation]]
- [[_COMMUNITY_Data Harmonisation|Data Harmonisation]]
- [[_COMMUNITY_Age-Sex-DIMD-SHA-Ethnicity Stratification|Age-Sex-DIMD-SHA-Ethnicity Stratification]]
- [[_COMMUNITY_GAMLSS Methodology|GAMLSS Methodology]]

## God Nodes (most connected - your core abstractions)
1. `ZenodoAssetManager (R6 class)` - 43 edges
2. `SynthPop (R6 class)` - 37 edges
3. `Simulation (R6 class, partial - summaries)` - 25 edges
4. `lsoa_pop_toR.R - LSOA mid-year population estimates ingestion` - 25 edges
5. `Auxiliary Functions (aux_fn.R)` - 21 edges
6. `HSE_ts.fst (Health Survey England time series)` - 18 edges
7. `R package: gamlss` - 15 edges
8. `CLAUDE.md project instructions` - 14 edges
9. `Design R6 class` - 14 edges
10. `SynthPop$initialize` - 14 edges

## Surprising Connections (you probably didn't know these)
- `root .vscode/debug.R env capture` --semantically_similar_to--> `package .vscode/debug.R env capture`  [INFERRED] [semantically similar]
  .vscode/debug.R → Rpackage/IMPACTncd_England_model_pkg/.vscode/debug.R
- `SynthPop$pop dataframe (Rcpp DataFrame dt)` --conceptually_related_to--> `SynthPop R6 class`  [INFERRED]
  auxil/IMPACTncd_sim.cpp → Rpackage/IMPACTncd_England_model_pkg/R/SynthPop_class.R
- `Disease.to_cpp() output (List l)` --conceptually_related_to--> `Disease R6 class`  [INFERRED]
  auxil/IMPACTncd_sim.cpp → Rpackage/IMPACTncd_England_model_pkg/R/Disease_class.R
- `Disease engine reproducibility (byte-identical runs)` --tests--> `simcpp() (C++ disease engine)`  [INFERRED]
  testing/test_reproducibility.R → Rpackage/IMPACTncd_England_model_pkg/src/IMPACTncd_sim.cpp
- `simulation.R runner script` --references--> `SynthPop synthetic population concept`  [INFERRED]
  simulation.R → CLAUDE.md

## Hyperedges (group relationships)
- **Core R6 class hierarchy** — design_class, disease_class, exposure_class, exposureeffect_class, simulation_class [EXTRACTED 1.00]
- **Initialization chain: global.R loads pkg, simulation.R instantiates Simulation from YAML** — global_init_script, simulation_runner_script, simulation_class, design_class, yaml_sim_design [EXTRACTED 1.00]
- **Inputs pipeline: YAML config + RR CSVY + exposure FST feed Design** — yaml_sim_design, rr_csvy_inputs, exposure_distributions_inputs, design_class, exposureeffect_class, exposure_class [EXTRACTED 1.00]
- **C++ simulation engines exposed to R** — rcppexports, rcpp_simcpp, rcpp_simcpp_year_based, rcpp_simsmok, simulation_class [EXTRACTED 1.00]
- **SynthPop generation/load lifecycle** — synthpop_class_initialize, synthpop_class_gen_synthpop, synthpop_class_gen_synthpop_filename, synthpop_class_gen_checksum, synthpop_class_write_synthpop, synthpop_class_delete_incomplete_synthpop [INFERRED 0.85]
- **Zenodo archive creation and upload pipeline** — zenodomanager_class_create_input_archives, zenodomanager_class_create_grouped_archive, zenodomanager_class_create_archive, zenodomanager_class_upload_archive, zenodomanager_class_create_new_record, zenodomanager_class_create_new_version, zenodomanager_class_publish_record [INFERRED 0.85]
- **Zenodo manifest compute/save/load/compare workflow** — zenodomanager_class_compute_manifest, zenodomanager_class_save_manifest, zenodomanager_class_load_manifest, zenodomanager_class_compare_manifests, zenodomanager_class_compute_file_hash, zenodomanager_class_check_source_changes [INFERRED 0.85]
- **Zenodo input sync (check/download/upload)** — zenodomanager_class_sync_inputs, zenodomanager_class_sync_check, zenodomanager_class_sync_download, zenodomanager_class_sync_upload [INFERRED 0.90]
- **Simulation summary exports (LE, prevalence, incidence, mortality, QALY, cost)** — simulation_class_export_summaries, simulation_class_export_summaries_hlpr, simulation_class_export_le_summaries, simulation_class_export_hle_summaries, simulation_class_export_prvl_summaries, simulation_class_export_incd_summaries, simulation_class_export_mrtl_summaries, simulation_class_export_qalys_summaries, simulation_class_export_costs_summaries, simulation_class_export_cms_summaries [INFERRED 0.90]
- **Simulation table exports (main, mortality, disease characteristics, exposures)** — simulation_class_export_tables, simulation_class_export_tables_hlpr, simulation_class_tbl_smmrs_core, simulation_class_build_strata_config, simulation_class_export_main_tables, simulation_class_export_all_cause_mrtl_tables, simulation_class_export_disease_characteristics_tables, simulation_class_export_xps_tables [INFERRED 0.90]
- **Smoking cessation/policy simulation engine** — simsmok_simsmok, simsmok_simsmok_sc, simsmok_simsmok_postcalibration, simsmok_simsmok_cessation, simsmok_simsmok_complete_cessation, simsmok_simsmok_policy_impact_incr, simsmok_simsmok_policy_impact_decr [INFERRED 0.90]
- **Core disease microsimulation engine (C++)** — impactncd_sim_simcpp, impactncd_sim_simcpp_year_based, impactncd_sim_get_simul_meta, impactncd_sim_get_disease_meta [INFERRED 0.95]
- **Disease epi scripts implement GAMLSS-based incidence/prevalence/fatality/duration estimation** — alcohol_misuse_epi, anxiety_depression_spell_epi, asthma_spell_epi, atrial_fibrillation_epi, chd_epi, copd_epi, chronic_kidney_disease_epi, connective_tissue_disorder_epi, constipation_epi, dementia_epi, diabetes_excl_type_2_epi, epilepsy_epi, hearing_loss_epi, heart_failure_epi, hypertension_epi, ibs_epi, nonmodelled_epi, obesity_epi, other_cancers_epi, pain_epi, primary_malignancy_breast_epi, method_gamlss_fitting, method_incd_prvl_ftlt_dur [INFERRED 0.85]
- **Common stratification strata: year, age, sex, dimd, sha, ethnicity** — alcohol_misuse_epi, anxiety_depression_spell_epi, asthma_spell_epi, atrial_fibrillation_epi, chd_epi, copd_epi, chronic_kidney_disease_epi, dementia_epi, heart_failure_epi, hypertension_epi, obesity_epi, primary_malignancy_breast_epi, method_age_sex_dimd_stratification [EXTRACTED 0.95]
- **All disease epi scripts consume CPRD 2021 panel_short data assets** — alcohol_misuse_epi, anxiety_depression_spell_epi, asthma_spell_epi, atrial_fibrillation_epi, chd_epi, copd_epi, chronic_kidney_disease_epi, connective_tissue_disorder_epi, constipation_epi, dementia_epi, diabetes_excl_type_2_epi, epilepsy_epi, hearing_loss_epi, heart_failure_epi, hypertension_epi, ibs_epi, nonmodelled_epi, obesity_epi, other_cancers_epi, pain_epi, primary_malignancy_breast_epi, data_source_cprd_2021 [EXTRACTED 0.95]
- **All disease epi scripts call harmonise() helper** — alcohol_misuse_epi, anxiety_depression_spell_epi, asthma_spell_epi, atrial_fibrillation_epi, chd_epi, copd_epi, chronic_kidney_disease_epi, connective_tissue_disorder_epi, constipation_epi, dementia_epi, diabetes_excl_type_2_epi, epilepsy_epi, hearing_loss_epi, heart_failure_epi, hypertension_epi, ibs_epi, nonmodelled_epi, obesity_epi, other_cancers_epi, pain_epi, primary_malignancy_breast_epi, method_harmonise, aux_fn_r [EXTRACTED 0.90]
- **GAMLSS-based disease epi modeling pipeline** — primary_malignancy_colorectal_epi_disease, primary_malignancy_lung_epi_disease, primary_malignancy_prostate_epi_disease, psychosis_epi_disease, rheumatoid_arthritis_epi_disease, stroke_epi_disease, type_2_diabetes_mellitus_epi_disease, concept_gamlss_methodology [INFERRED 0.85]
- **Shared age-sex-DIMD-SHA-ethnicity stratification** — primary_malignancy_colorectal_epi_disease, primary_malignancy_lung_epi_disease, psychosis_epi_disease, rheumatoid_arthritis_epi_disease, stroke_epi_disease, type_2_diabetes_mellitus_epi_disease, anxiety_depression_epi_disease, asthma_epi_disease, concept_age_sex_dimd_strata [INFERRED 0.90]
- **Scripts sourcing aux_fn.R helpers** — primary_malignancy_colorectal_epi_disease, primary_malignancy_lung_epi_disease, primary_malignancy_prostate_epi_disease, psychosis_epi_disease, rheumatoid_arthritis_epi_disease, stroke_epi_disease, type_2_diabetes_mellitus_epi_disease, correlations_disease_corr, recover_dur_recovery, anxiety_depression_epi_disease, asthma_epi_disease, cld_epi_disease, ckd45_epi_disease, chronic_kidney_disease_s3_epi_disease, parfoutputs_casenum_plots, aux_fn_helpers [EXTRACTED 1.00]
- **Scripts consuming CPRD 2021 panel data** — primary_malignancy_colorectal_epi_disease, primary_malignancy_lung_epi_disease, primary_malignancy_prostate_epi_disease, psychosis_epi_disease, rheumatoid_arthritis_epi_disease, stroke_epi_disease, type_2_diabetes_mellitus_epi_disease, correlations_disease_corr, anxiety_depression_epi_disease, asthma_epi_disease, cld_epi_disease, chronic_kidney_disease_s3_epi_disease, data_source_cprd2021 [EXTRACTED 1.00]
- **Cancer epi modeling family** — primary_malignancy_colorectal_epi_disease, primary_malignancy_lung_epi_disease, primary_malignancy_prostate_epi_disease [INFERRED 0.80]
- **CKD stage epi modeling family** — chronic_kidney_disease_s3_epi_disease, chronic_kidney_disease_s4_epi_disease, chronic_kidney_disease_s5_epi_disease, ckd45_epi_disease [INFERRED 0.85]
- **HFSS dietary intake fitting suite** — fit_hfss_energy_model, fit_hfss_food_weight_model, fit_hfss_salt_model [INFERRED 0.85]
- **Cardiometabolic exposure fits (BMI, alcohol, BP-med, AF, CKD, famCVD)** — fit_bmi_model, fit_bmi_model_2024, fit_alcohol_model, fit_alcohol_model_2024, fit_bpmed_model, fit_af_diag_model, fit_ckd_model, fit_famcvd_model [INFERRED 0.80]
- **Sociodemographic fits (education, ethnicity)** — fit_education_model, fit_education_model_2024, fit_education_model_for_imputation, fit_ethnicity_model, fit_ethnicitygrp_model_2024 [INFERRED 0.85]
- **Behavioral exposure fits (alcohol, ETS, fruit)** — fit_alcohol_model, fit_alcohol_model_2024, fit_ets_model, fit_ets_model_2024, fit_fruit_model, fit_fruit_model_2024 [INFERRED 0.80]
- **Scripts consuming HSE_ts data** — extract_hse_correlation_structure_2024, extract_hse_sociodemographics_for_validation, fit_af_diag_model, fit_alcohol_model, fit_alcohol_model_2024, fit_bmi_model, fit_bmi_model_2024, fit_bpmed_model, fit_ckd_model, fit_education_model, fit_education_model_2024, fit_education_model_for_imputation, fit_ethnicity_model, fit_ethnicitygrp_model_2024, fit_ets_model, fit_ets_model_2024, fit_famcvd_model, fit_fruit_model, fit_fruit_model_2024 [INFERRED 0.90]
- **GAMLSS distribution fitting scripts** — fit_hfss_energy_model, fit_hfss_food_weight_model, fit_hfss_salt_model, fit_af_diag_model, fit_alcohol_model, fit_alcohol_model_2024, fit_bmi_model, fit_bmi_model_2024, fit_bpmed_model, fit_ets_model, fit_ets_model_2024, fit_famcvd_model, fit_fruit_model, fit_fruit_model_2024 [INFERRED 0.85]
- **2024 refresh of fitting scripts (use HSE_ts_03_19.fst directly)** — extract_hse_correlation_structure_2024, fit_alcohol_model_2024, fit_bmi_model_2024, fit_education_model_2024, fit_ethnicitygrp_model_2024, fit_ets_model_2024, fit_fruit_model_2024 [INFERRED 0.85]
- **T2DM fitting workflow (prevalence + diagnosis + duration)** — fit_t2dm_model_script, fit_t2dm_dgn_model_script, fit_t2dm_dur_model_script, concept_t2dm [INFERRED 0.85]
- **HSE time series construction pipeline** — generate_hse_ts_script, generate_hse_ts_15_18_script, generate_hse_ts_3_19_script, datasource_hse [INFERRED 0.80]
- **Mortality projection and calibration pipeline** — preprocess_mrtl_script, mrtl_projections_script, mrtl_projections_coherent_script, prepare_calibr_file_script [INFERRED 0.80]
- **2024 exposure model refresh (HSE 2003-2019)** — fit_tchol_model_2024_script, fit_veg_model_2024_script, datasource_hse_ts_03_19 [INFERRED 0.75]

## Communities

### Community 0 - "Cardiometabolic Risk Fitting (T2DM/HDL/TChol)"
Cohesion: 0.06
Nodes (39): All-cause mortality projections, LSOA ethnicity-by-age-sex population, Mortality calibration, Type 2 Diabetes Mellitus (T2DM), T2DM diagnosis (dm_dgn), T2DM duration (dm_dur), Total cholesterol (tchol), Vegetable portions (vegpor) (+31 more)

### Community 1 - "CPRD Validation & HFSS Energy Models"
Cohesion: 0.04
Nodes (38): HFSS energy intake, HFSS food weight (g), HFSS salt intake, Smoking relapse probability, breast_ca, *_incd incidence outputs, *_mrtl mortality outputs, prostate_ca (+30 more)

### Community 2 - "Core Architecture & Bootstrap"
Cohesion: 0.05
Nodes (57): arrow package, CKutils package dependency, CLAUDE.md project instructions, Key YAML Parameters, EQ-5D-5L QALY calculation, Scenario Definition Pattern, config: clusternumber, config: diseases (+49 more)

### Community 3 - "Simulation Testing Scripts"
Cohesion: 0.05
Nodes (39): Export tables comparison test, Health Survey for England (HSE), Design$new(), Discounting logic for QALYs and costs, active_days_curr_xps (physical activity exposure), alcohol_curr_xps (alcohol exposure), bmi_curr_xps (BMI exposure), fruit_curr_xps (fruit exposure) (+31 more)

### Community 4 - "HSE Data Preparation"
Cohesion: 0.07
Nodes (34): Atrial fibrillation diagnosis, Alcohol exposure, BMI exposure, Blood pressure medication, CKD chronic kidney disease, Education, Ethnicity, ETS environmental tobacco smoke (+26 more)

### Community 5 - "Sociodemographic Distribution Fits"
Cohesion: 0.08
Nodes (30): distr_best_fit (gamlss helper), Physical activity days per week, Income (ordinal), Systolic blood pressure (mmHg), Relative risk of hypertension diagnosis given SBP, Smoking cessation event probability, Cigarettes per day (current smokers), Cigarettes per day (ex-smokers) (+22 more)

### Community 6 - "Disease Epi Helpers (panel data)"
Cohesion: 0.12
Nodes (37): Alcohol Misuse (epi), Anxiety/Depression spell (epi), Asthma spell (epi), Atrial Fibrillation (epi), aux_fn.R (disease_burden helpers), aux_fn.R helper functions, chd disease, Chronic Kidney Disease (epi) (+29 more)

### Community 7 - "Zenodo Asset Management"
Cohesion: 0.05
Nodes (44): package: httr2, package: zen4R, Vignette: Inputs manifest system, Vignette: Zenodo data management, ZenodoAssetManager$check_connection (private), ZenodoAssetManager$check_source_changes, ZenodoAssetManager$compare_manifests, ZenodoAssetManager$compare_with_inputs_manifest (+36 more)

### Community 8 - "Population Estimates & SynthPop Seeding"
Cohesion: 0.07
Nodes (44): CCG - Clinical Commissioning Group, DIMD - IMD Decile (1 most deprived to 10 least), IMD - Index of Multiple Deprivation, LAD - Local Authority District, LSOA - Lower Super Output Area (~1500 people), MSOA - Middle Layer Super Output Area, QIMD - IMD Quintile, RGN - Government Office Region (+36 more)

### Community 9 - "Cancer Epidemiology Modeling"
Cohesion: 0.1
Nodes (43): Anxiety/Depression Epi (OLD), Asthma Epi (OLD), Auxiliary Functions (aux_fn.R), CKD Stage 3 Epi (OLD), CKD Stage 4 Epi (OLD), CKD Stage 5 Epi (OLD), CKD Stages 4-5 Epi (OLD), Chronic Liver Disease Epi (OLD) (+35 more)

### Community 10 - "C++ Simulation Engine (top-level)"
Cohesion: 0.08
Nodes (35): DiseaseIncidenceType2(), DiseaseIncidenceType3(), EvalDiagnosis(), EvalDiseaseIncidence(), EvalMortality(), get_disease_meta(), get_simul_meta(), GetStackTrace() (+27 more)

### Community 11 - "SynthPop R6 Class Methods"
Cohesion: 0.06
Nodes (41): Design (R6 class, defined elsewhere), package: R6, SynthPop$check_integrity, SynthPop$count_synthpop, SynthPop$deep_clone (private), SynthPop$del_incomplete (private), SynthPop$delete_incomplete_synthpop, SynthPop$delete_synthpop (+33 more)

### Community 12 - "Simulation Export & Summaries"
Cohesion: 0.13
Nodes (27): Simulation$build_strata_config (private), Simulation$calc_costs (private), Simulation$calc_QALYs (private), Simulation$export_all_cause_mrtl_by_dis_summaries (private), Simulation$export_all_cause_mrtl_tables (private), Simulation$export_cms_summaries (private, CMS), Simulation$export_costs_summaries (private), Simulation$export_dis_char_summaries (private) (+19 more)

### Community 13 - "C++ Disease Incidence Core"
Cohesion: 0.11
Nodes (24): Disease diagnosis simulation, Disease.to_cpp() output (List l), Disease incidence simulation, Disease mortality simulation, DiseaseIncidenceType2, DiseaseIncidenceType3, disease_meta struct, EvalDiagnosis (+16 more)

### Community 14 - "Health Economics Generators"
Cohesion: 0.14
Nodes (6): generate_eq5d_decr(), generate_health_econ(), generate_healthcare_costs(), generate_informal_care_costs(), generate_productivity_costs(), generate_socialcare_costs()

### Community 15 - "C++ Disease Incidence Types"
Cohesion: 0.24
Nodes (13): DiseaseIncidenceType2(), DiseaseIncidenceType3(), EvalDiagnosis(), EvalDiseaseIncidence(), EvalMortality(), get_disease_meta(), get_dur_forward(), get_dur_forward_prvl() (+5 more)

### Community 16 - "Reproducibility Testing"
Cohesion: 0.67
Nodes (2): Disease engine reproducibility (byte-identical runs), simcpp() (C++ disease engine)

### Community 17 - "Cancer History Inputs (PLCO)"
Cohesion: 0.67
Nodes (2): History of cancer (PLCO input), Maddams et al. 2012 cancer prevalence projections

### Community 18 - "Zenodo Upload Tests"
Cohesion: 0.67
Nodes (2): keyring (R package), ZenodoAssetManager (R6 class)

### Community 19 - "VSCode Debug Env Capture"
Cohesion: 1.0
Nodes (2): package .vscode/debug.R env capture, root .vscode/debug.R env capture

### Community 20 - "Compression Helpers (deflate/inflate)"
Cohesion: 1.0
Nodes (2): deflate (helper), inflate (helper)

### Community 21 - "GAMLSS Prediction Helpers"
Cohesion: 1.0
Nodes (2): centile_predictAll (helper), mean_predictAll (helper)

### Community 22 - "ONS Incidence/Mortality Fetchers"
Cohesion: 1.0
Nodes (2): get_ons_incd (helper), get_ons_mrtl (helper)

### Community 23 - "Smoking Cessation Sim"
Cohesion: 1.0
Nodes (2): simsmok_cessation (Rcpp export), simsmok_complete_cessation (Rcpp export)

### Community 24 - "Smoking Policy Impact Sim"
Cohesion: 1.0
Nodes (2): simsmok_policy_impact_decr (Rcpp export), simsmok_policy_impact_incr (Rcpp export)

### Community 25 - "Test Framework (tinytest)"
Cohesion: 1.0
Nodes (2): package: tinytest, tinytest.R runner

### Community 26 - "SynthPop Validation Plots"
Cohesion: 1.0
Nodes (2): plot_synthpop_val (helper), package: cowplot

### Community 28 - "Project Critical Rules"
Cohesion: 1.0
Nodes (1): Critical Rules (no delete, always test)

### Community 29 - "R6 Hierarchy Doc"
Cohesion: 1.0
Nodes (1): R6 Class Hierarchy section

### Community 30 - "Data Flow Diagram"
Cohesion: 1.0
Nodes (1): Data Flow (Inputs to Outputs)

### Community 31 - "Correlation Plot Helper"
Cohesion: 1.0
Nodes (1): plot_cor (helper)

### Community 32 - "Distribution Validation Helper"
Cohesion: 1.0
Nodes (1): distr_validation (helper)

### Community 33 - "RNG Generation Helper"
Cohesion: 1.0
Nodes (1): generate_rns (RNG helper)

### Community 34 - "Output Directory Helper"
Cohesion: 1.0
Nodes (1): output_dir (helper)

### Community 35 - "Smoking Cigarette Sim"
Cohesion: 1.0
Nodes (1): simsmok_cig (Rcpp export)

### Community 36 - "Smoking Cigarette Scenario Sim"
Cohesion: 1.0
Nodes (1): simsmok_cig_sc (Rcpp export)

### Community 37 - "Smoking Cessation Cigarette Sim"
Cohesion: 1.0
Nodes (1): simsmok_complete_cessation_cig (Rcpp export)

### Community 38 - "Rcpp Auto-Generated Bindings"
Cohesion: 1.0
Nodes (1): RcppExports.cpp (auto-generated bindings)

### Community 39 - "C++ Fortify Header (2)"
Cohesion: 1.0
Nodes (1): undef_fortify.h (header)

### Community 40 - "Stratification (Age/Sex/DIMD)"
Cohesion: 1.0
Nodes (1): Age/sex/dimd/sha/ethnicity stratification

### Community 41 - "Disease Burden Estimation"
Cohesion: 1.0
Nodes (1): Incidence/Prevalence/Fatality/Duration estimation pipeline

### Community 42 - "Data Harmonisation"
Cohesion: 1.0
Nodes (1): harmonise() data preparation

### Community 43 - "Age-Sex-DIMD-SHA-Ethnicity Stratification"
Cohesion: 1.0
Nodes (1): Age-Sex-DIMD-SHA-Ethnicity Stratification

### Community 44 - "GAMLSS Methodology"
Cohesion: 1.0
Nodes (1): GAMLSS distribution fitting methodology

## Knowledge Gaps
- **245 isolated node(s):** `Critical Rules (no delete, always test)`, `R6 Class Hierarchy section`, `Data Flow (Inputs to Outputs)`, `TODO transition notes`, `root .vscode/debug.R env capture` (+240 more)
  These have ≤1 connection - possible missing edges or undocumented components.
- **Thin community `Reproducibility Testing`** (4 nodes): `mk_scenario_init2()`, `Disease engine reproducibility (byte-identical runs)`, `simcpp() (C++ disease engine)`, `test_reproducibility.R`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Cancer History Inputs (PLCO)`** (3 nodes): `History of cancer (PLCO input)`, `Maddams et al. 2012 cancer prevalence projections`, `history_of_cancer_for_plco.R`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Zenodo Upload Tests`** (3 nodes): `keyring (R package)`, `test_zenodo_upload.R`, `ZenodoAssetManager (R6 class)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `VSCode Debug Env Capture`** (2 nodes): `package .vscode/debug.R env capture`, `root .vscode/debug.R env capture`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Compression Helpers (deflate/inflate)`** (2 nodes): `deflate (helper)`, `inflate (helper)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `GAMLSS Prediction Helpers`** (2 nodes): `centile_predictAll (helper)`, `mean_predictAll (helper)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `ONS Incidence/Mortality Fetchers`** (2 nodes): `get_ons_incd (helper)`, `get_ons_mrtl (helper)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Smoking Cessation Sim`** (2 nodes): `simsmok_cessation (Rcpp export)`, `simsmok_complete_cessation (Rcpp export)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Smoking Policy Impact Sim`** (2 nodes): `simsmok_policy_impact_decr (Rcpp export)`, `simsmok_policy_impact_incr (Rcpp export)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Test Framework (tinytest)`** (2 nodes): `package: tinytest`, `tinytest.R runner`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `SynthPop Validation Plots`** (2 nodes): `plot_synthpop_val (helper)`, `package: cowplot`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Project Critical Rules`** (1 nodes): `Critical Rules (no delete, always test)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `R6 Hierarchy Doc`** (1 nodes): `R6 Class Hierarchy section`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Data Flow Diagram`** (1 nodes): `Data Flow (Inputs to Outputs)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Correlation Plot Helper`** (1 nodes): `plot_cor (helper)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Distribution Validation Helper`** (1 nodes): `distr_validation (helper)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `RNG Generation Helper`** (1 nodes): `generate_rns (RNG helper)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Output Directory Helper`** (1 nodes): `output_dir (helper)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Smoking Cigarette Sim`** (1 nodes): `simsmok_cig (Rcpp export)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Smoking Cigarette Scenario Sim`** (1 nodes): `simsmok_cig_sc (Rcpp export)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Smoking Cessation Cigarette Sim`** (1 nodes): `simsmok_complete_cessation_cig (Rcpp export)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Rcpp Auto-Generated Bindings`** (1 nodes): `RcppExports.cpp (auto-generated bindings)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `C++ Fortify Header (2)`** (1 nodes): `undef_fortify.h (header)`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Stratification (Age/Sex/DIMD)`** (1 nodes): `Age/sex/dimd/sha/ethnicity stratification`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Disease Burden Estimation`** (1 nodes): `Incidence/Prevalence/Fatality/Duration estimation pipeline`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Data Harmonisation`** (1 nodes): `harmonise() data preparation`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `Age-Sex-DIMD-SHA-Ethnicity Stratification`** (1 nodes): `Age-Sex-DIMD-SHA-Ethnicity Stratification`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.
- **Thin community `GAMLSS Methodology`** (1 nodes): `GAMLSS distribution fitting methodology`
  Too small to be a meaningful cluster - may be noise or needs more connections extracted.

## Suggested Questions
_Questions this graph is uniquely positioned to answer:_

- **Why does `BCT distribution family` connect `CPRD Validation & HFSS Energy Models` to `Cardiometabolic Risk Fitting (T2DM/HDL/TChol)`?**
  _High betweenness centrality (0.195) - this node is a cross-community bridge._
- **Why does `SynthPop (R6 class)` connect `SynthPop R6 Class Methods` to `CPRD Validation & HFSS Energy Models`, `Disease Epi Helpers (panel data)`?**
  _High betweenness centrality (0.137) - this node is a cross-community bridge._
- **Why does `R package: data.table` connect `Cancer Epidemiology Modeling` to `Cardiometabolic Risk Fitting (T2DM/HDL/TChol)`, `Simulation Testing Scripts`?**
  _High betweenness centrality (0.091) - this node is a cross-community bridge._
- **Are the 3 inferred relationships involving `ZenodoAssetManager (R6 class)` (e.g. with `package: zen4R` and `package: httr2`) actually correct?**
  _`ZenodoAssetManager (R6 class)` has 3 INFERRED edges - model-reasoned connections that need verification._
- **Are the 2 inferred relationships involving `SynthPop (R6 class)` (e.g. with `SynthPop$initialize` and `Vignette: How to test run`) actually correct?**
  _`SynthPop (R6 class)` has 2 INFERRED edges - model-reasoned connections that need verification._
- **What connects `Critical Rules (no delete, always test)`, `R6 Class Hierarchy section`, `Data Flow (Inputs to Outputs)` to the rest of the system?**
  _245 weakly-connected nodes found - possible documentation gaps or missing edges._
- **Should `Cardiometabolic Risk Fitting (T2DM/HDL/TChol)` be split into smaller, more focused modules?**
  _Cohesion score 0.06 - nodes in this community are weakly interconnected._