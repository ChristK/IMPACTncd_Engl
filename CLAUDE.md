# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Critical Rules

- **NEVER delete files** from this repository unless the user explicitly instructs you to do so.
- **Always test your implementations.** Verify that changes work as expected before presenting them. Keep iterating until they do. Never deliver untested solutions expecting the user to test them.

## Project Overview

IMPACTncd_Engl is a **microsimulation framework** for modeling chronic disease impacts in England. It simulates disease progression and risk factor effects using synthetic populations to evaluate "what-if" policy scenarios against baselines.

**Core workflow**: Baseline scenario → Policy scenario(s) → Export summaries → Analysis

## Architecture

### R6 Class Hierarchy
Main classes in `Rpackage/IMPACTncd_England_model_pkg/R/`:
- **Design**: Loads/validates YAML config from `inputs/sim_design.yaml`
- **SynthPop**: Manages synthetic population data (persisted as `.fst` files)
- **Simulation**: Orchestrates runs, scenarios, multicore execution (split across 4 files: `Simulation_class*.R`)
- **Disease/Exposure/ExposureEffect**: Handle disease and exposure-specific logic

### Performance-Critical C++ Code
Located in `Rpackage/IMPACTncd_England_model_pkg/src/`:
- `IMPACTncd_sim.cpp` - Core disease simulation engine
- `simsmok.cpp` - Smoking simulation engine

### Data Flow
- **Inputs**: YAML config, `.fst` exposure distributions, CSVY relative risks in `inputs/`
- **Processing**: SynthPop generation → Baseline run → Scenario modifications → Disease simulation
- **Outputs**: `.fst` files in `output_dir` (configured in YAML, typically external storage)

## Common Commands

### Initial Setup
```r
source("global.R")  # Initializes environment, installs deps, builds local package
```
Package rebuilds automatically when source changes are detected.

### Download model data (required on a fresh clone)
Large input/simulation data is not in git — it lives on Zenodo (concept DOI
`10.5281/zenodo.20812409`, CC-BY-SA-4.0). Download once; no Zenodo token needed
(published data is public):
```r
IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")
IMPACTncd$zenodo_connect()        # defaults to the published production record
IMPACTncd$zenodo_download_all()   # inputs + PARFs/RR (~13 GB)
```
See `vignette("zenodo_data_management")`. Data managers upload new versions with
`zenodo_upload_all(version = "...", publish = TRUE)` (requires `ZENODO_TOKEN`).

### Running Simulations
```r
source("global.R")
IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")
IMPACTncd$del_logs()$del_outputs()$run(1:2, multicore = TRUE, "sc0")$export_summaries(multicore = TRUE)
```

**Per-scenario comparators.** `comparator_scenario` accepts either a bare name
(historical: everything else vs it — byte-identical output) or a **named vector**
mapping intervention -> comparator, with at most one *unnamed* element as the
default: `c("sc0", mult_base_prs = "mult_base_qrisk")`. A scenario can be both an
intervention and someone else's comparator; chains are legal, cycles rejected.
Exactly **one comparator per scenario**, which is what keeps `scenario` a unique
row key and lets `comparator` be attached *after* quantiling — so it never enters
the key vector `x`, a group-by, one of the ~9 positional `setnames()` calls, or a
filename. The map gate `.comparator_is_map()` is **syntactic**, never data-derived,
so the schema can't vary between families or reruns. In map mode the contrast
tables gain a `comparator` column and `comparators.csv` is written; validation
(unknown names, cycles, orphans, `comparator`/`.cmp__` collisions) runs up front in
the parent. Gotcha: at all five join sites `cmp` must come from the FULL table, not
`d`'s complement — under a map the two sets overlap. Index the map with
`as.character(scenario)`: a factor uses level codes and mis-maps silently.
Tests: `test_comparator_map.R`.

**CEA cost perspectives** (`cea = TRUE`), one file each: `societal`
(`total_cost`, already net of economic output), `healthcare` (`healthcare_cost`)
and `healthcare_socialcare` (`healthcare_cost + socialcare_cost` — NICE's NHS+PSS
reference case, and usually the one to quote for a UK appraisal). Custom `*_costs`
always join societal, and join BOTH narrower perspectives when named in
`custom_costs_in_healthcare`. A perspective whose required built-in columns are
missing is **skipped**, never silently degraded — otherwise a run without
`socialcare_cost` would publish `healthcare_socialcare` tables identical to the
`healthcare` ones. Tests: `test_cea_perspectives.R`.

**Multi-level discounting (`export_tables()`).** `qaly_discount_rate` /
`cost_discount_rate` take **vectors** (default `c(0, 3.5)`) paired element-wise
into discount *levels*; a length-1 rate recycles against a longer one. Every
level is reported side by side, tagged in a `discount` column, so one export
yields the undiscounted **and** the NICE-3.5% figures instead of two runs into
separate directories. Differential rates need equal-length vectors —
`qaly = c(0, 1.5)` with `cost = c(0, 4)` is 2 levels, **not** 4.
Discounting is **no longer CEA-only**: it also reaches `qalys` / `net_qalys` /
`costs` / `net_costs`. Prevalence, incidence, mortality and exposure tables carry
no `discount` column (they are rates and counts, not flows), and the equity
tables are excluded too because their metric set mixes QALYs with plain event
counts — only the **main** and **CEA** tasks are passed `discount_levels`.
Four things are load-bearing:
**(a)** `expand_discount_levels()` scales the **per-year** values *before* any
`cumsum()`, so the `_cuml` columns accumulate present values rather than
discounting an already-summed total.
**(b)** Rates are **paired** with each other but **crossed** with `wtp`, giving
`length(wtp) × nrow(discount_levels)` NMB rows; a threshold is a ratio of present
values, so it is never itself discounted, and `ICER` appears once per level
rather than once per threshold.
**(c)** `discount` is excluded from the filename key `str6`, so the *set of
files written* is identical however many levels are requested — the levels live
in rows, not filenames.
**(d)** Every entry point defaults to `build_discount_levels(0, 0)` when called
without levels, so an internal caller cannot silently acquire discounting.
⚠️ **Breaking for downstream readers** of the qalys / costs / net_* tables: they
gained a `discount` column and now carry one block of rows per level, so a reader
that assumed one row per (scenario, year, stratum) will double-count unless it
filters on `discount`. The `0%` rows reproduce the old undiscounted figures.
**Year arguments take a two-digit shorthand** (`19` = 2019), because the
summaries store short years. `baseline_year_for_change_outputs` had always been
promoted; `discount_from_year` was **not**, so `discount_from_year = 19` was read
as the year 19 AD — the exponent became `year - 19` ≈ 2000, `1/1.035^2000`
underflowed, and every discounted figure was published as **zero** with no error
or warning. Both now go through the single helper `promote_short_year()`, so the
promotion cannot be added to one year argument and forgotten for the other.
Ported from IMPACTncd_Japan; the four helpers' executable code is byte-identical
to that model's, but the currency defaults are not carried over (`wtp` stays
`c(25000, 35000)` GBP). Tests: `test_discount_levels.R`.

**Integer summaries silently truncate rates — `.promote_integer_measures()`.** Every
rate in `Simulation_class_tables.R` is one summed column divided by another, written
back with `:=`. The form decides whether that is safe: `d[, v := v/n]` replaces the
whole column and **promotes** it to double (harmless), but a **join-update**
`d[cases, on = ..., v := v/i.v]` writes into a *subset* of the existing column and must
coerce to its type — so an **integer** column truncates, turning a case-fatality rate of
0.034 into `0`. data.table only warns, and a table of zeros looks like a result.
Every division site in `tbl_smmrs_core()` and `export_all_cause_mrtl_tables()` is the
join-update form. `tbl_smmrs_core()` had always guarded; `export_all_cause_mrtl_tables()`
had not, and was protected only by the stock summaries happening to store counts as
`double` — a property of whoever wrote the parquet, not of this code. Both now call the
shared `.promote_integer_measures(dt, keys)`. Tests:
`test_integer_measure_promotion.R` (incl. the plain-vs-join `:=` distinction, and an
integer fixture that must equal a double one end to end).

**`two_agegrps` + `two_agegrps_split_age`.** `two_agegrps = TRUE` reports `agegrp`
as two coarse groups and writes to `tables2agegrps/`. The split is now an argument:
`two_agegrps_split_age` (default `65L` = the historical `30-64` / `65-99`). The group
labels are derived from the **data**, not hard-coded 30..99, so a run reporting ages
20-89 gets `20-64` / `65-89` and an open-ended top band gives `65+`. Core helpers:
`two_agegrp_map()` / `collapse_two_agegrps()` (`Simulation_class_tables.R`, next to
`add_qimd_from_dimd()` — it is the same kind of pure relabel-before-the-group-by, and
therefore exact for counts and for summed-numerator/summed-denominator rates; `equity`
is the exception, since it fits a *model*).
Until 2026-08-05 the flag was plumbed into only the **main** and **equity** tasks, so
`export_all_cause_mrtl_tables()` silently ignored it and published 5-year bands inside
`tables2agegrps/` — files **md5-identical** to their `tables/` namesakes. A silent
no-op, not an error. It now reaches every family whose summaries carry `agegrp`: main,
all-cause-mortality-by-disease, CEA and equity. Out of scope by construction:
`dis_characteristics` has no `agegrp` column, and the `xps` family uses `agegrp20`
(`30-49`/`50-69`/`70-89`/`90+`/`All`) — a different variable that 65 would cut in half.
⚠️ **Breaking for downstream readers** of `tables2agegrps/all-cause mortality given
disease-*agegroup*.csv` written before 2026-08-05: those were 5-year-band tables.
The split must fall *between* bands, never inside one — `62` is refused because it
would cut `60-64` — and must leave both groups non-empty; it need **not** start a band
that is actually present, so a locality with an empty band at the split still works.
Both checks run in the **parent** (`validate_two_agegrps_split_age()` structurally,
`check_two_agegrps_split_age()` against a projected arrow scan of `agegrp`), so a typo
fails in seconds rather than from inside a forked worker.
Tests: `test_two_agegrps.R`; production verification:
`auxil/verify_two_agegrps_on_real_output.R` (pins the 18 main-family files
byte-identical to what the old code published).

**Baseline-only runs: turn off `cea` and `equity` explicitly.** Every contrast
(`cypp`/`cpp`/`dpp`/`net_qalys`/`net_costs`, all CEA tables, all equity indices)
needs an intervention arm, so with only `sc0` in the summaries none of those files
is written. The tasks are still *dispatched*, and — because the scenario list is a
property of the data, not of the arguments — each reads its summaries **before** it
can discover there is nothing to contrast: CEA loads `qalys` + `costs` (~3 GB), the
equity worker loads `incd`, `prvl`, `mrtl`, `qalys` in turn (~7 GB peak). Passing
`cea = FALSE, equity = FALSE` skips both reads and drops `length(tasks)` from 6 to
4, which also lowers the `min(length(tasks), clusternumber_export)` fork count.
Output is unchanged: those files never existed. In `export_main_tables()` the
per-metric working copy `tt <- copy(tt_base)` is taken **below** the same guard, so
a skipped comparison metric no longer duplicates the summary (it used to, ten times
per call). The copy must stay above the `two_agegrps` relabel and the ftlt
`absorb_dt()` — both mutate `tt` by reference, and dropping it makes the second
`ftlt` metric fail outright. Tests: `test_export_main_tables_working_copy.R`.
Memory profile: `auxil/diagnose_export_tables_memory.R`, `..._peak.R`,
`auxil/diagnose_skipped_workers.R`.

**Output column prefixes** come from `.tbl_col_prefixes()` (`Simulation_class_tables.R`),
one entry per metric, applied to every quantile column. The convention: a level metric
takes `<family>_rate_`/`_mean_`, and BOTH of its change variants take a prefix distinct
from the level and from each other, always containing `change` — that is what makes
`grep("change", names(dt))` a safe way to pick change columns downstream.
Every **relative** prefix is additionally unique across families. Two entries were
corrected on 2026-08-05, both publishing a name the values did not support:
**(a)** `ftlt_change_relative` mapped to `ftlt_rate_`, the *level* prefix — a
ratio-to-baseline (0.72 = 72% of the baseline year) under a name claiming to be a rate,
which also hid those six case-fatality tables from any name-pattern selection. Now
`ftlt_change_relative_`.
**(b)** `prvl_change_relative` and `incd_change_relative` both mapped to
`prct_change_relative_`, asserting a *percentage* where the value is a **ratio**
(`value / value[baseline_year]`: 1.0 = no change, 0.9 = a 10% fall). Taken literally that
understates every effect 100-fold while looking entirely plausible. Now
`prvl_change_relative_` / `incd_change_relative_`.
`abs_change_` is deliberately **not** renamed — an absolute change is exactly what that
name claims, so it is merely generic (shared by `prvl`/`incd`, disambiguated by the
filename), not wrong; the uniqueness test therefore covers relative only.
⚠️ **Breaking for downstream readers** of `*relative change by *.csv` written before
2026-08-05 — values are unchanged, only the column name.
Tests: `test_tbl_col_prefixes.R`.

**Life-expectancy tables (`life_expectancy = TRUE`, default) — new 2026-09-02.**
`export_summaries()` had always written `le`, `le60`, `hle_1st_cond` and
`hle_cmsmm1.5`, but **no table family read them**, so a run produced the summaries
and no CSV. `export_le_tables()` publishes all four in one file per
stratification (`life expectancy by <strata> (...).csv`), distinguished by a
`type` column, in both standardisations. Strata element names: `le_ons` /
`le_esp` (the vocabulary is now **ten**, not eight).
The values are **weight-averaged ages at an event**, not period life-table
expectancies: each summary is
`SELECT <strata>, SUM(w) AS popsize, SUM(age * w)/SUM(w) AS LE FROM lc_table
WHERE all_cause_mrtl > 0 [AND age > 60]`, so `popsize` is the weighted **count of
events**, not a population. Collapsing over strata is therefore **exact** — the
weights are the denominators the means were formed with — but only as a
`weighted.mean(value, popsize)`; averaging or summing them is wrong. Each figure
carries `n_events` (median weighted count over mc) so a thin stratum is visible.
`two_agegrps` and discounting do not apply (no age axis in either mode; an age at
death is not a flow). `qimd` is derived from `dimd` as everywhere else.
Tests: `test_le_tables.R` (25). Production verification:
`claude_process/ineq_diag/verify_le_tables_on_real_output.R` (recomputes a cell
independently from the parquet; on the Ineq run the 2043 most-vs-least-deprived
LE gap goes 11.14 y in `sc0` to 4.38 y under `PW_*_scn`).

**Exposure quantiles (`export_xps_quantiles: yes`, default `no`) — new 2026-09-02.**
`export_xps()` only ever stored `weighted.mean` per exposure, so a *shift in
distribution* was unrecoverable from output — two scenarios can share a mean BMI
while only one thins the upper tail. With the flag on it additionally writes
`xps20_qnt` / `xps5_qnt` (long format: one row per exposure per probability,
probabilities from `xps_quantile_probs`, default the 5th-95th centiles), and
`export_xps_tables()` publishes `exposure quantiles by <strata> (...).csv`.
Off by default because it is not free, and **`xps_quantile_strata` is the knob
that decides how much**. Cost is ~linear in (stratifications) x (columns) and
**not** in the number of groups, so the 8 stratifications the *means* need are
what made it expensive. Measured end to end through `export_xps()` at the stock
England shape (n = 100,000 x 31 years, 18 exposure columns of which 4 are
indicators), against 8.6 s for the means alone:

| config | per call | added |
|---|---|---|
| all 8 mean stratifications, every column | 26.5 s | +17.9 s |
| default `list("year", c("year","qimd"))`, informative columns only | 13.3 s | **+4.7 s** |

For the 17-scenario x 100-iteration Ineq design over 10 cores that is **~3.4 h
-> ~0.5 h**. Three things get it there, all exact:
**(a)** `xps_quantile_strata` (Design-level, default the two above) — every
entry must contain `year` and may name only `agegrp20`/`sex`/`qimd`, validated
in the **parent**; the *default* is intersected with what is available rather
than validated, because a design whose `strata_for_output` omits `dimd` has no
`qimd` and a default must never be able to fail.
**(b)** indicator and constant columns are skipped — a quantile of `ets` /
`statin_px` / `smok_never` / `smok_active` can only return 0 or 1, so it says
nothing the mean does not. Detected with `range()`, **not** `unique()`: hashing
3M distinct doubles once per exposure cost more than the quantiles it saved.
They keep their means.
**(c)** the `all_cause_mrtl >= 0 & year >= init_year & age >= ageL` subset is
taken **once** and shared — xps20, xps5 and both quantile calls used a
byte-identical filter, i.e. three or four full copies of `sp$pop`. This also
made the means ~8% faster; verified byte-identical to the inline-subset form.
Adding more `xps_quantile_probs` is nearly free (the sort dominates), so widen
the probabilities freely and the stratifications sparingly.
Benchmarks: `claude_process/ineq_diag/bench_xps_quantiles{,_v2,_v3}.R`
(v2 decomposes the cost: 26 of 33 s was the per-group sort, only 6.9 s
data.table machinery; a sort-once formulation was tried and is only 1.4x, so
the readable `groupingsets` form was kept — the two agree exactly).
End-to-end check incl. the unchanged means:
`claude_process/ineq_diag/verify_xps_quantile_speedup.R`.
`weighted_quantile()` (`aux_functions.R`, exported) is the **inverse-CDF (type-1)**
weighted quantile. It is a thin wrapper over the internal
`.weighted_quantile_pre()`, which takes the weight cleaning already done — the
one implementation, so there is no second copy to drift, and `export_xps()`
hoists `!is.na(w) & w > 0` out of its per-exposure loop instead of repeating it
~18 times over the same weight vector (~20% of the quantile cost). A 300-case
fuzz test pins the hoisted path identical to the public one. It is — no interpolation, deliberately, because several exposures are
discrete (`active_days` 0-7, integer pack-years) and interpolation would publish
values the variable cannot take. Verified against `stats::quantile(type = 1)` and
for replication invariance (integer weights == repeating the values).
⚠️ **Quantiles are NOT additive** — unlike every count in `summaries/` and unlike
the means. They cannot be averaged across mc to get "the" quantile and cannot be
combined across strata. In the CSVs the row-level `prob` is a centile of the
*exposure*; the `xps_qnt_*` column names are *Monte Carlo* uncertainty on it.
Also fixed alongside: `export_xps()`'s "All"-marginal relabel looped over **every**
column, so a legitimately `NaN` exposure mean (a stratum with no ex-smokers, say)
had `"All"` coerced into a double column — `NA` behind a "NAs introduced by
coercion" warning on every run. Now restricted to the strata columns, so a `NaN`
stays `NaN` (both are `is.na()`; on disk it is now a float NaN rather than a null).

**`export_xps_tables()` was materially broken — four fixes (2026-09-02).** Found
by running the exposure family end to end for the first time (the Ineq run's
`xps/` had been emptied by `del_outputs()`, so nobody had seen its output).
All four were **pre-existing** and affected the exposure *means*, not just the
new quantiles:
**(a)** `year` was promoted twice. Every reader under `summaries/` does
`year := year + 2000L`, and that line was copied here — but `export_xps()`
already promotes the year before writing, so the tables published **4013-4043**.
Now via `promote_xps_year()`, which promotes only if `max(year) < 200`, so the
convention lives in one place and a mixed tree cannot publish two.
**(b)** `detect_filterable_vars()` tested `is.character()`, but these datasets
come back from parquet **dictionary-encoded as factors**. So it found nothing,
`filter_expr` fell back to `quote(TRUE)`, and **every stratification's file got
the whole dataset** — the year-only marginal averaged together with the
`agegrp20 x sex x qimd` cells, i.e. a mean of means instead of the marginal.
**(c)** Same root cause: files were written for stratifications the data cannot
support (`... by year-agegroup ...` whose `agegrp20` is entirely `"All"`), and
**(d)** `standardised_vars` came out empty, naming files `( standardised)`.
Fixing (b) exposed a fifth: the ESP label omitted `age`. The xps5 dataset is
age-standardised *by construction* (weighted by `wt_esp`, no age axis), so
`agegrp20` is not among its filterable columns and `age` could never appear —
understating every ESP filename (`(sex-qimd standardised)` for what is
age-sex-qimd) and leaving the fully-stratified one labelled with nothing.
`age` is now **stated** in that branch, not derived: it can never be one of the
strata there, so deriving it with a `setdiff` against them implied otherwise,
and a stray `agegrp20` in `strata$xps_esp` removed it — labelling an
age-standardised table `(sex-qimd standardised)`. Only the remaining axes are
derived. Matches the main families' `... (age-sex-dimd standardised).csv`.
And a sixth, found by testing that stray config: naming a column the loaded xps
dataset lacks (`agegrp20` in `strata$xps_esp`) reached `keyby` and raised
``object 'agegrp20' not found`` **from inside a forked worker**, after the whole
dataset had been read. `xps_*` is exempt from `check_strata_columns()` (those
datasets are not built from `summaries/`, so the parent has nothing to check
against), which left this gap. `skip_missing_strata()` now warns and skips just
that stratification — the same answer the equity family gives.
⚠️ **Breaking for any reader of `exposures by *.csv`**: the years change by
-2000, the non-marginal files change value (they were averaging marginals in),
and most ESP filenames gain `age-`. Note `export_tables()` never wipes
`tables/`, so a re-export leaves the old wrongly-named files behind as orphans
— clear the directory or `del_outputs()` first.
Why the tests missed it: the fixture built its strata columns as **character**,
so `is.character()` succeeded there and failed in production; and its marginal
rows were the exact unweighted mean of their cells, so filtered and unfiltered
results coincided. The fixture now uses factors and a population-**weighted**
marginal, and asserts the two differ. Tests: `test_xps_quantiles.R` (387).

**Pack-years are reported among EVER-SMOKERS — `*_pop_curr_xps` (2026-09-02).**
`export_xps()` sets `smok_packyrs`/`smok_dur`/`smok_cig`/`smok_quit_yrs` to `NA`
for never-smokers (`smok_status_curr_xps == "1"`) and aggregates with
`weighted.mean(na.rm = TRUE)`, so those means are **conditional on ever having
smoked** — never-smokers are dropped from the denominator, not counted as zero.
Consequence, and the reason it looks like a paradox: a scenario that moves people
into never-smoking leaves a more heavily-smoking remainder, so **mean pack-years
RISES** while population-total pack-years falls. The conditional mean is also
`NaN` in a stratum with no ever-smokers, whereas a population mean is defined
everywhere. `smok_packyrs_pop_curr_xps`, `smok_dur_pop_curr_xps` and
`smok_cig_pop_curr_xps` now carry the population-weighted reading alongside
(never-smokers as 0), created before the `NA`s and picked up by the same
`_curr_xps$` grep, and dropped again before the lifecourse write. No companion
for `smok_quit_yrs`: a population mean over never- and current smokers alike,
neither of whom has years since quitting, would mean nothing.
Tests: `test_xps_quantiles.R` (52, incl. an end-to-end `export_xps()` run
asserting the conditional mean is blind to a scenario the population mean sees).

`export_tables()` turns the summaries into CSVs. Besides the main / mortality /
disease-characteristics / exposure / cost-effectiveness (`cea = TRUE`) tables, it also
writes **equity slope-index** tables (`equity = TRUE`, default): absolute (`AEI_total`,
`AEI_per100k`) and relative (`REI_rel`, `RII_ratio`) analogues of the Slope/Relative Index
of Inequality for cumulative CPP, CYPP, DPP and net QALYs across DIMD deciles (pro-poor sign
convention; MC iterations give the uncertainty), plus `total_benefit`, `fit_R2`, `n_mc` and
`p_pro_poor` columns. The estimator is the grouped-data weighted least squares of Kakwani,
Wagstaff & van Doorslaer (1997); `RII_ratio` is the Kunst–Mackenbach ratio of fitted extremes,
*not* Moreno-Betancur's RII (they criticise that formula; theirs is `exp(beta)` from a
log-linear model, impossible on a signed quantity). Three things bite:
**(a)** `AEI_total = SII * sum(PY)` is a *hypothetical* whole-population figure, **not** the
gap between the extreme deciles (~8x smaller) and not Cookson's published gap
(`* (J-1)/J^2`); it scales with total person-years, so it grows with the **horizon** as well
as with population — only `AEI_per100k`/`REI_rel` compare across localities, and even
`AEI_total` only compares at a fixed reporting year.
**(b)** Both relative indices are `NA` unless the fitted benefit is **positive** — with a
negative mean they invert, publishing a pro-poor result as pro-rich.
**(c)** `all_diseases_sum` excludes the design's umbrella conditions (incidence `type: 0`:
`dm`, `ctdra`, `cancer`) as well as `cms*`, via `umbrella_disease_names()` — summing an
umbrella beside its own components double-counts ~9.5% of events and distorts the slope.
The ridit reference population is selectable via
`equity_ridit_reference` (`"comparator"` default vs `"scenario"` — the latter violates
Cookson's equal-group-size precondition, so it is a sensitivity analysis; see the Renard 2019
decision support in the vignette). Core math lives in
`calc_equity_slope_indices()` / `export_equity_tables()` in `Simulation_class_tables.R`; see
`vignette("understanding_model_outputs")`.

**Person-years denominator (2026-08-13) — the numerator's exposure base must match.**
`B_k(T)` cumulates a flow over every year from `baseline_year`; it used to be divided by
`N_k(T)`, the **single reporting-year** population, and that failed the estimator's own
neutrality test. For a policy giving every living person the identical benefit each year,
`y_k = sum_t c(t) N_k(t) / N_k(T)` is systematically smaller for groups whose population
**grew** — and the more deprived deciles grow fastest in the ONS projection input (HFSS
England comparator, 2026→2043: +21.2% in the most deprived quintile vs +1.5% in the least).
So the estimator published a **pro-rich gradient out of demography alone**, growing with the
horizon, which distorted the *trend* as well as the level. `N` is now cumulative person-years
`cumsum(N)` over the same window, used for the ridit, the weights **and** the denominator —
all three must share the numerator's base or the invariance fails. Measured on the HFSS run
(`baseline_year = 2026`) the fix removes ~a third of the magnitude, sign intact: `AEI_total`
at 2043 −179,077 → −115,490 (35.5% artefact), and the share **grows with the horizon**
(15.5% at 2030 → 35.5% at 2043) — a trend distortion, not just a level shift.
Note the mirror-image trap the code had already avoided:
pinning `N` at `baseline_year` was rejected for reintroducing the same thing from the other
side. ⚠️ **Breaking**: every equity value changes, and `AEI_per100k` changes **units** (per
100k people → per **100k person-years** — which finally makes it comparable with a published
annual-rate SII, reversing the old vignette warning). Exact for an undiscounted stream;
second-order residual under discounting (a group whose person-years fall later genuinely
receives less present value per person-year — time preference, not artefact).
Tests: `test_equity_person_years_discount.R` (the neutral-stream fixture asserts exact zero
AND that the old formula gives −32 on the same fixture, so it can't pass on unfixed code).
Production verification: `claude_process/verify_equity_fix_on_nwl.R`.

**Equity tables now carry `discount` and `p_pro_poor`.** Discounting follows the same
per-metric rule as the main tables: `net_qalys` is a QALY flow → one block per
`qaly_discount_rate` level; `cpp`/`cypp`/`dpp` are counts the model never discounts → single
`"0%"` block. All four gain the column so the family shares one schema; levels live in rows,
not filenames; scaling is applied per-year **before** the cumsum. Gotcha: discounting is
**not** a uniform rescaling of a slope index — it re-weights the years making up the
cumulative benefit, so where the gradient changes sign over the horizon the discounted index
can in principle be *larger* in magnitude (it shrinks in all 18 years of the HFSS run, so
this is a caution, not a common case — but don't assume monotonicity).
`p_pro_poor` is the share of the same finite draws on the pro-poor side of each type's own
null: `> 0` for `AEI_*`/`REI_rel`, `> 1` for `RII_ratio`, `NA` for `total_benefit`/`fit_R2`.
It answers the *directional* question a percentile interval can't (HFSS 2043: median
−115,490 with a 95% interval of −269,921..+95,399, but `p_pro_poor = 0.07`).
⚠️ **Quote `p_pro_poor` from `AEI_total`/`AEI_per100k` ONLY.** On `REI_rel`/`RII_ratio` it is
biased **upwards (too pro-poor)** and is not a probability the policy is pro-poor. Not noise —
**selection**, so it does not shrink with more draws, and it is worst exactly when the
pro-rich signal is strongest. Mechanism: the relative indices are `NA` unless the fitted
benefit is positive, and the draws that fail that are overwhelmingly the steeply *pro-rich*
ones (a steep negative slope is what drags the fitted most-deprived end below zero) — so the
share is taken over a sample the pro-rich draws were filtered out of. Test fixture, 4 draws
(3 pro-poor, 1 pro-rich): `AEI_total` → `n_mc=4, p=0.75`; `RII_ratio` → `n_mc=3, p=1.00`.
Diagnosis: a relative index's `n_mc` **below** the absolute indices' `n_mc` is the tell, and
the shortfall is the number of dropped draws; equality is the only condition under which a
relative `p_pro_poor` is trustworthy. `AEI_total`/`AEI_per100k` never drop a draw.
`discount` is now in `.equity_reserved_vars`, so it cannot be a stratum.

**`strata_for_output` with both `dimd` and `qimd`.** The xps tables map `agegrp ->
agegrp20` and `dimd -> qimd` (exposures are reported in 20-year bands and quintiles).
That mapping is **many-to-one**, so listing both sides of a pair used to yield
`by = c(..., "qimd", "qimd")` and crash every parallel task of `run()` inside
`groupingsets()` with "Argument 'by' must have unique column names for grouping".
`xps_strata_from_output()` (`Simulation_class.R`) now de-duplicates it. Tests:
`test_xps_strata_from_output.R`.

**Deciles vs quintiles (any table family).** Every `strata` element accepts `qimd`
instead of `dimd`, and `qimd` need **not** be in `strata_for_output`: whenever a loaded
summary carries `dimd`, `add_qimd_from_dimd()` (`Simulation_class_tables.R`) relabels it
into quintiles and the caller's existing group-by does the collapsing. That is **exact** —
aggregations group-then-sum, so summing two deciles equals aggregating by quintile at
summary time, and rates follow since they are summed-numerator/summed-denominator *after*
the group-by. Putting `qimd` in `strata_for_output` changes no result, only file size
(it adds no grouping granularity). The converse is impossible: `dimd` strata need `dimd`
in the summaries, and are unavailable in the `xps_*` tables (export_xps writes those
quintile-based). Shared deprivation helpers: `.deprivation_vars`, `.deprivation_levels`,
`.dimd_to_qimd`, `deprivation_rank()` (errors on unrecognised labels rather than
silently mis-mapping). Tests: `test_qimd_derivation.R`.

**Single years of age (`export_summaries(single_year_of_age = TRUE)`).** Swaps
`agegrp` for `age` in the stratification (`strata <- strata_age` — one line,
`Simulation_class_summaries.R`), so **no** summary keeps `agegrp`: 20 of the 30
datasets gain `age`, the other 10 (`le`, `le60`, `hle_*`, `dis_characteristics`)
take `strata_noagegrp` and have neither in either mode. ~5x the rows (70 single
years vs 14 bands on the stock 30-99 design).
`export_tables()` merges `strata` over its defaults **by name only** — ten
independent keys (`ons`, `esp`, `mrtl_ons`, `mrtl_esp`, `disease_char`,
`xps_ons`, `xps_esp`, `le_ons`, `le_esp`, `equity`), whole-element replacement,
no coupling — and
two of them hard-code `agegrp`: `ons` (task 1 main, reused verbatim by task 5
CEA) and `mrtl_ons` (task 2). So overriding `ons` alone used to abort the export
with ``task 2 failed - "object 'agegrp' not found"`` — the message form
`keyby = eval(outstrata)` produces for a missing column, raised from inside a
forked worker after GBs were loaded, *after* the two agegrp-free defaults had
already been written.
Fixed 2026-08-18 by **`add_agegrp_from_age()`** (next to `add_qimd_from_dimd()`,
called at every one of the 9 summary-load sites): derives the bands from `age`
before any group-by, so the shipped `agegrp` defaults just work. Deriving, not
substituting `age` into the defaults — a substitution multiplies the group-by
~5x inside a worker already peaking in the tens of GB, and renames every file
whose strata name `agegrp`. Bands come from `CKutils::to_agegrp(dt, 5L, 99L, ...)`, the same
call `Simulation_class.R` used to build the column originally, with `min_age`
left data-derived so `ageL: 20` stays consistent; the column is **character**,
matching the on-disk type. Exact for counts and summed-num/summed-denom rates
(relabel-before-the-group-by, same argument as `qimd`); the converse is
impossible.
Also guards the top edge: `to_agegrp()` builds its age vector over
`max(actual, max_age)` but its LABELS from `max_age` alone, so on a design with
`ageH > 99` age 100 recycled onto the first band and was published as `30-34`
behind a "recycled with remainder" warning. The helper now extends `max_age` to
the data.
Three gotchas: **(a)** the denominator must be derived too — `pp_scaled` is
aggregated on the same strata before being joined to the numerator, so deriving
on one side only fails in *that* group-by with the same `object 'agegrp' not
found`. **(b)** the derivation must precede `collapse_two_agegrps()`, which
returns FALSE and does nothing when `agegrp` is absent. Every family had it right
except **equity**, which collapsed first — caught in review, fixed the same day.
That ordering is *silent* (a `tables2agegrps/` file md5-identical to its
`tables/` namesake — the 2026-08-05 signature) and worst in equity, because that
family fits a *model*: the two-group and 5-year-band fits are different
estimands, not a relabel. `check_two_agegrps_split_age()` also falls back to
deriving from `age`, so the parent-side probe is not silently skipped.
**(c)** `age` is a legal stratum but 5x the cells — `c("year","age")` is cheap,
`c("year","age","sex","dimd")` is ~347M rows post-melt.
Measured inflation of the summaries themselves: **exactly 4.99x the rows**,
2.90-5.05x the bytes (six main summaries, Gloucester single-year vs Wales banded;
`claude_process/measure_single_year_inflation.R`).
**None of this changes a banded run's output.** Verified by building the package
at the pre-change commit into a separate lib and diffing real output: full
`export_tables()` on Wales (banded, 204 files) and equity+CEA at both
`two_agegrps` settings on ineqCK (banded, 17 scenarios, 52 files) are **byte-
identical** old vs new. Scripts: `claude_process/verify_no_change_to_agegrp_runs.R`,
`..._equity_cea.R`, `claude_process/compare_agegrp_run_md5.R`.
Two parent-side guards were added alongside (2026-08-18), both for mistakes that
used to surface as ``task N failed - "object 'x' not found"`` from a fork:
`build_strata_config()` now **errors** on an element name outside the ten
(with a prefix-aware "did you mean" — `mrtl` -> `mrtl_ons`) instead of dropping
it silently, and `check_strata_columns()` errors on a stratum the summaries
cannot supply (the reverse case, `age` against a banded run, which is
underivable). The column check reads **schemas only**, applies the same
`qimd`-from-`dimd` / `agegrp`-from-`age` rules as the workers
(`.derivable_strata_columns()`), and flags a column only when it is missing from
**every** existing dataset that family reads — so it cannot fire where the
worker would have coped. Exempt by design: `xps_*` (not built from `summaries/`;
it filters pre-aggregated `"All"` rows) and `equity` (its per-table
warn-and-skip text is pinned by `test_export_equity_tables.R`).
Tests: `test_strata_validation.R` (43).

Known wart, not fixed: `cms` is passed both `strata` and `strata_age`, which the
flag makes the same vector, so `cms_score_*` and `cms_score_by_age_*` are
byte-identical (verified: 400 files, ~140 MB/run). The flag is recorded nowhere
on disk — the parquet schema is the only evidence of which mode ran.
Tests: `test_agegrp_derivation.R` (73 — mapping, no-ops, ordering vs
`two_agegrps` through the real main/acm/CEA/equity export paths, exactness vs a
natively banded fixture, ages above the top edge). Two meaningfulness proofs live
in the gitignored `claude_process/` and are therefore NOT in CI:
`verify_test_fails_without_fix.R` (neuters the helper; must reproduce
`object 'agegrp' not found`) and `verify_equity_ordering_regression.R` (rebuilds
the function with the pre-fix ordering; must publish 5-year bands).

**Equity strata + gradient axis.** Each entry of `strata$equity` is one table. Any
column the summaries carry (`strata_for_output`) can be an output stratum — the gradient
is fit *within* it, identically to filtering to that stratum and fitting without it.
Exceptions: `dimd`/`qimd` don't stratify, they *select the gradient axis* (deprivation is
consumed into the index); one token → that axis, both → one table per axis, neither →
`dimd` so the defaults keep their historical filenames. Enforced: `year` is required
(rows are cumulative-to-that-year against that year's reference population), `mc`/
`scenario` are reserved. The axis is echoed in a `gradient` column and (when the caller
wrote it) in the filename suffix, where `agegrp` is rewritten to `agegroup` as in every
other table family. `qimd` need not be in `strata_for_output` — it is
derived from `dimd` by the standard decile-pair collapse. The derived **counts** are exact
(counts are additive and the collapse preserves the ridit midpoint), but here — unlike every
other table family — a *model* is fit, so the **index** is invariant only under an exactly
linear gradient; under curvature the five-point and ten-point fits are different estimands
(a step-shaped gradient shifts the SII ~3.5%). The reverse is impossible. Unavailable
gradient/stratum columns and a missing comparator scenario raise `warning()`s, not
`logs`-gated messages, and skip only the affected table. `two_agegrps` collapses `agegrp`
here too. Tests: `test_equity_slope_index.R`, `test_export_equity_tables.R`.

### Testing
```r
tinytest::test_package("IMPACTncdEngland")
```

### Rebuild Documentation
```r
roxygen2::roxygenise("Rpackage/IMPACTncd_England_model_pkg")
```

### Docker Environment
```bash
./docker_setup/setup_user_docker_env.sh          # Interactive setup
./docker_setup/setup_dev_docker_env.sh           # Developer setup with mounted source
./docker_setup/docker_build_push.sh Dockerfile.IMPACTncdENGL --push  # Build/push images
```

**Three-layer images** (build in order): `prerequisite` (R env) → `data`
(prerequisite + ~13GB Zenodo data) → `impactncdengl` (data + model code). The
`data` layer is heavy (~33GB) and changes rarely, so it is built **locally** and
tagged by Zenodo version:
```bash
./docker_setup/build_push_data.sh   # auto-tags data.impactncdengl:<zenodo_ver> + :latest, pushes
```

**The image does NOT ship `scenarios/`.** Those drivers and their
`sim_design_*.yaml` are project-specific working files, and
`/IMPACTncd_England/scenarios` is where `setup_user_docker_env.sh -ScenariosDir`
bind-mounts the *user's* folder — a baked copy is shadowed exactly when it would
be used. The image keeps the path as an empty, world-writable mount point.
**Three** code-entry paths, **three** mechanisms — they are not interchangeable:
**(a)** the **copy** method is covered by `.dockerignore` (context filtering, so
the files never enter the `COPY` layer); **(b)** CI's **git** method clones
*inside* the build, where `.dockerignore` has no reach, so
`Dockerfile.IMPACTncdENGL` prunes `/tmp/repo/scenarios` in the **same `RUN`** as
the clone; **(c)** `build_model_via_commit.sh` (the containerd-GC-race fallback)
has neither a build context nor a Dockerfile — it pipes `git archive` into
`docker cp` — so the exclusion lives in its pathspec
(`git archive HEAD -- . ':(exclude)scenarios'`). Miss any one and that path
publishes an image unlike the other two.
Same-RUN matters in (b): an `rm -rf` in a later layer leaves the files in the
earlier one, still pushed and still extractable. Both (b) and (c) then assert the
folder is empty and **fail the build** if it is not, so a future overlay path
cannot silently re-add them. The CI workflow needs no filtering step of its own.
Verification (gitignored): `claude_process/test_docker_excludes_scenarios.sh` —
builds all three paths against `data.impactncdengl:local`, probes the copy-method
build *context*, and runs a negative control on the pre-change Dockerfile.

**That same asymmetry had already broken the copy method** (2026-07-02 →
2026-09-04, fixed): `.dockerignore` excludes `docker_setup/`, so the "Refresh the
container ENTRYPOINT" step could not find `docker_setup/entrypoint.sh` and every
`./docker_build_push.sh Dockerfile.IMPACTncdENGL` died there — while CI (clone)
and the commit fallback (`git archive`) sailed past, because neither filters
through `.dockerignore`. It stayed hidden because a failing layer-3 `docker
build` on this host is *documented* as the containerd GC race, whose fix is "use
`build_model_via_commit.sh`" — the same workaround absorbs both failures. Fixed
with a `!docker_setup/entrypoint.sh` negation. Rule of thumb: a `RUN` that reads
a repo file from inside the image must be checked against `.dockerignore`, not
just against git.

Downstream Dockerfiles parametrize their base via build-args
(`FROM ${DATA_IMAGE}` / `FROM ${PREREQ_IMAGE}`, default `:local`), so CI/other
builds point `FROM` at the pushed image (e.g.
`--build-arg DATA_IMAGE=chriskypri/data.impactncdengl:latest`) without
re-downloading from Zenodo. CI: `prerequisite` and `impactncdengl` build in
GitHub Actions; the model build pulls the data layers from Docker Hub.
If a layer-3 build trips the host containerd GC race
(`blob ... not found` at export), use `./docker_setup/build_model_via_commit.sh`.

## Code Conventions

### R Code Style
- **Always use data.table** for data manipulation (not dplyr/tidyverse)
- **Reference semantics**: Use `:=` for in-place modification, avoid copies
- **Memory efficiency**: Store integers as `integer` not `double`

### Column Naming
- `*_curr_xps` - Current exposure value
- `*_lagged` - Lagged/historical values
- `*_prvl` - Prevalence
- `*_incd` - Incidence
- `*_mrtl` - Mortality

### Key Exposure Variables (Primary Prevention)
In `synthpop$pop`:
- `bmi_curr_xps` - Body mass index (kg/m²)
- `sbp_curr_xps` - Systolic blood pressure (mmHg)
- `tchol_curr_xps` - Total cholesterol (mmol/L)
- `smok_status_curr_xps` - Smoking status (1=never, 2=occasional, 3=ex, 4=current)
- `smok_cig_curr_xps` - Cigarettes per day
- `alcohol_curr_xps` - Alcohol consumption (g/day)
- `fruit_curr_xps`, `veg_curr_xps` - Fruit/veg (g/day, 80g = 1 portion)
- `active_days_curr_xps` - Physical activity days per week

### Key Disease Variables (Secondary Prevention)
- `dm`, `t1dm`, `t2dm` - Diabetes
- `chd`, `stroke`, `af`, `hf` - Cardiovascular diseases
- `ckd`, `colorect_ca`, `prostate_ca`, `breast_ca`, `lung_ca`

### Scenario Definition Pattern
```r
IMPACTncd$update_primary_prevention_scn(
  function(sp) {
    sp$pop[year >= 23L, bmi_curr_xps := bmi_curr_xps * 0.9]  # 10% reduction
  }
)
IMPACTncd$run(1:n_runs, multicore = TRUE, "sc1")
```

**Deprivation-change scenarios**: if a scenario forces a simulant's `dimd`/`qimd`
(to model a deprivation-reduction policy), those forced values drive the simulation
(they are live covariates for disease incidence/duration and smoking), but ALL outputs
(lifecourse, summaries, tables, equity, xps) are stratified by each simulant's **original**
pre-scenario deprivation. The original is snapshotted before the scenario runs and restored
after `simcpp()` automatically (in `Simulation_class.R` `run_sim`); scenarios that don't
touch `dimd`/`qimd` are unaffected.

**`dimd` vs `qimd` are two different channels — `parf_deprivation` (2026-09-02).**
`qimd` keys **every** exposure distribution table (all 36 in
`inputs/exposure_distributions/`) and the whole smoking sub-model (`simsmok.cpp`
reads `qimd`, and `qimd_sc` for its scenario variant): it is the **modelled**
channel. `dimd` keys the cached PARF table (`simulation/parf/`), every
disease-burden input (`*_prvl`, `*_incd`, `*_ftlt`, `*_dur`) and
`IMPACTncd_sim.cpp`'s stochastic disease-duration predictor: it is the
**non-modelled** channel — the part of the deprivation gradient the modelled
risk factors do not explain.
`Disease$gen_parf()` recomputes nothing from a scenario's exposures; it reads the
cached table and **joins it with `on = .NATURAL`**, and `dimd` is one of those
keys. Its payload `p0 = mu * (1 - parf)` / `m0 = mu_mrtl * (1 - parf_mrtl)`
carries the *observed*, decile-specific rates, consumed positionally as
`private$parf$p0` in `set_incd_prb()` and `private$parf$m0` in `set_mrtl_prb()`
for the whole horizon. So **when** the join happens decides which decile's
background risk each simulant gets.
`update_primary_prevention_scn(fn, parf_deprivation = )` now states it:
`"original"` (default) joins at the pre-scenario `dimd` — only the regenerated
exposures act; `"scenario"` keeps the forced `dimd` and transfers the
non-modelled gradient too. Applied by `with_original_dimd()` (`aux_functions.R`)
around the `gen_parf()`+`set_init_prvl()` block in `run_sim`, with `on.exit()`
restoration so a throwing disease cannot leave the rest of the run on the wrong
deprivation. A **no-op** whenever `dimd == dimd_orig`, so every scenario that
does not model deprivation change is bit-unaffected; `scenarios/Ineq/` is the
only script in the repo that rewrites deprivation.
⚠️ **This was a silent regression.** Commit `adbf762` (2026-04-28, "move
set_init_prvl() after the scenario") moved that block from *before* the scenario
to *after* it, so the framework began doing the re-join unconditionally. The Ineq
paper's "Modeled" and "Non-Modeled" arms differed **only** in that the
Non-Modeled one called `gen_parf()` itself after rewriting `dimd`; that call
became redundant and **both arms became the Non-Modeled arm**. Measured: 13 of 24
arm-pairs bit-identical in the production summaries, the rest differing by a mean
signed 1.7e-5 of either sign; population-mean COPD `p0` falls 35.6% from the join
timing alone (most-deprived decile to 0.39x its own value), residual all-cause
`m0` 17.0%; extreme-decile `p0` ratios reach 3.0x (COPD, psychosis) and 12 of 21
diseases differ >25%. Net QALYs went 4.9M -> 22.9M on the Modeled arm while the
Non-Modeled arm reconciled with the old figures once discounting was matched.
Not covered, deliberately (it is pre-existing and unchanged): `set_incd_prb()`
computes its own `clbfctr = sum(mu)/sum(rp)` keyed on `.(dimd, sex)` at
`init_year` and looks `mu` up at the *current* `dimd`. `clbfctr` is a ratio of
observed to modelled (≈1), so for a **future**-change scenario
(`sc_year > init_year`, where `init_year`'s `dimd` is still original) the effect
is second-order. For a **past**-change scenario it is not: `mu` is then taken at
the forced decile, so incidence is partly re-calibrated to it whatever
`parf_deprivation` says.
Tests: `test_parf_deprivation.R` (36). Meaningfulness proof (gitignored, not in
CI): `claude_process/ineq_diag/verify_parf_test_fails_without_fix.R` neuters
`with_original_dimd()` and 8 assertions must fail.

**RNG restore in `run_sim` (2026-09-02).** Both scenario hooks saved
`rs <- .Random.seed` and restored with `set.seed(rs)`. `set.seed()` keeps only
`rs[1]` — the RNG **kind** code, a constant (10403) — so this silently reseeded
the base-R stream from one fixed value after every scenario and every iteration
instead of restoring it. Now `assign(".Random.seed", rs, envir = globalenv())`,
the idiom `Disease_class.R`'s `with_random()` already used. Result-neutral in
current code: every base-R consumer downstream sets its own seed
(`sample_int_expj` at `Disease_class.R:843`/`943`, `runif` inside
`with_random()`), and `dqrng` — which the C++ engine uses — was always restored
correctly via `dqrng_set_state()`. Audit:
`claude_process/ineq_diag/audit_base_rng_consumers.R`. One exception: if a PARF
cache is **missing**, `gen_parf_files()` now draws from a different stream, so
regenerate PARFs before a comparison intended to be bit-exact.

### File Paths
- **Absolute paths required** for `output_dir`, `synthpop_dir` in YAML
- **Relative paths** for inputs from project root (e.g., `./inputs/`)
- External storage typically at `/mnt/storage_fast/`

### Exposure table resolution (`exposure_definitions[*].file_name`)
`file_name` (the parquet dataset backing each exposure) resolves 3 ways, in
`Design$resolve_exposure_path()`:
- **bare name** (`bmi`) → stock `./inputs/exposure_distributions/<name>`
- **relative path with `/`** (`tables/bmi_local`) → resolved against the
  **sim_design YAML's own directory** (so Docker users can ship custom tables
  alongside their scenario scripts; the folder is bind-mounted)
- **absolute path** → used verbatim
Use `Design$xps_table_path("<file_name>")` for the reads that bypass `Exposure`
objects (PARF generation in `Disease_class.R`, `smok_relapse` in
`SynthPop_class.R`) so re-pointed tables stay consistent everywhere. See
`vignette("how_to_run_scenarios")` → "Bring your own exposure tables".

## Key YAML Parameters (inputs/sim_design.yaml)

- `n` - Population size per chunk
- `num_chunks` - Number of synthpops (total population = n × num_chunks)
- `clusternumber` - CPU cores for simulation (~10GB RAM per core)
- `clusternumber_export` - CPU cores for export (may need lower than clusternumber)
- `locality` - "England" or subnational area name
- `simulation_files_overwrite` - `yes` to overwrite existing outputs

**Important**: Copy `sim_design.yaml` before modifying to prevent merge conflicts.

## Key Dependencies

- **data.table** - Fast data manipulation (core throughout)
- **Rcpp/dqrng** - C++ integration and quality RNG
- **foreach/doParallel** - Multicore execution
- **fst** - Fast columnar storage for SynthPop
- **gamlss** - Distribution fitting for exposures
- **CKutils** - Custom utilities (installed from GitHub: `ChristK/CKutils`)

## Debugging

- Set `validation: yes` in YAML for additional checks (slower)
- For C++ debugging: Set `PKG_CXXFLAGS= -O0 -ggdb3` in Makevars and rebuild
- Set `dev_mode <- TRUE` in `global.R` for verbose startup logging
- VS Code: Use `.vscode/debug` task

### Run logs (`logs: yes`)

`<output_dir>/logs/` holds `times.txt` (phase + per-iteration timestamps),
`log<mc>.txt` (one per worker) and **`console.txt`** — everything the *parent*
process would have printed during `run()` / `export_summaries()` /
`export_tables()`, including `foreach`'s `.verbose` bookkeeping. Written by
`start_console_log()` / `stop_console_log()` (`Simulation_class.R`, next to
`time_mark()`); phases append behind a `===== <phase> at: <ts> =====` banner.
Tests: `test_console_log.R`.

**Why it's not just convenience.** `logs: yes` also sets `.verbose = TRUE` on
every `foreach`, and foreach's accumulator emits *all* per-task bookkeeping in
one burst when the loop returns (~21 KB for 200 tasks). Writing a burst to a
**terminal** can block forever: a stray `^S` stops the tty, the stop survives
detaching a `screen` session, and the process then sleeps in `write()` with no
CPU, no I/O and no error — indistinguishable from a hung simulation. A Wales run
lost **17h29m** to exactly this on 2026-07-28 (all 200 workers done at 15:59;
parent resumed 3 s after the screen was reattached at 09:28 next morning). Files
have no line discipline and cannot be stopped this way — and are ~20% faster to
write than a pty.

Caveats: it's an **R-level** redirect (captures `cat`/`print`/`message`/`Rcout`,
*not* writes a linked C library makes straight to fd 1/2), and it is released
via `on.exit()` so fatal errors still reach the terminal. For unattended runs
also redirect at the shell: `Rscript scenario.R > run.log 2>&1`.

Related: `run()`'s end-of-function sink cleanup unwinds only to the depth seen
on entry (it used to unwind to 0, destroying a caller's `capture.output()` /
knitr sink), and `sink.number(type = "message")` is a **connection number**
(2 = none, >2 = active), not a count — comparing it to 0 is always true.

**Diagnosing a slow run:** `auxil/profile_times_log.R <times.txt>` breaks one
run into phases, per-iteration durations, realised concurrency and any dead
gaps; `auxil/scan_times_gaps.R <times.txt>...` scans many runs for the tail gap
between the last worker finishing and `End of parallelisation` (should be ~0 s).

## Common Pitfalls

1. **Windows .Random.seed error** - `global.R` initializes RNG before Simulation creation
2. **BLAS threading** - Single-threaded mode forced to avoid fork issues with OpenMP
3. **Memory for exports** - `clusternumber_export` may need lower value than `clusternumber`
4. **Docker detection** - Code checks `file.exists("/.dockerenv")` for container-specific config

## Documentation

Vignettes in `Rpackage/IMPACTncd_England_model_pkg/vignettes/`:
- `how_to_test_run.Rmd` - Quick start guide
- `how_to_run_scenarios.Rmd` - Scenario setup
- `understanding_model_outputs.Rmd` - Output interpretation
- `inputs_manifest_system.Rmd` - Data asset management
- `zenodo_data_management.Rmd` - Zenodo integration

Access via: `vignette(package = "IMPACTncdEngland")`

## graphify

This project has a graphify knowledge graph at graphify-out/.

Rules:
- Before answering architecture or codebase questions, read graphify-out/GRAPH_REPORT.md for god nodes and community structure
- If graphify-out/wiki/index.md exists, navigate it instead of reading raw files
- For cross-module "how does X relate to Y" questions, prefer `graphify query "<question>"`, `graphify path "<A>" "<B>"`, or `graphify explain "<concept>"` over grep — these traverse the graph's EXTRACTED + INFERRED edges instead of scanning files
- After modifying code files in this session, run `graphify update .` to keep the graph current (AST-only, no API cost)
