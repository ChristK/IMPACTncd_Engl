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

**Output column prefixes** come from `.tbl_col_prefixes()` (`Simulation_class_tables.R`),
one entry per metric, applied to every quantile column. The convention: a level metric
takes `<family>_rate_`/`_mean_`, and BOTH of its change variants take a prefix distinct
from the level and from each other, always containing `change` — that is what makes
`grep("change", names(dt))` a safe way to pick change columns downstream.
`ftlt_change_relative` used to map to `ftlt_rate_`, the *level* prefix: it published a
ratio-to-baseline (0.72 = 72% of the baseline year) under a name claiming to be a rate,
and hid those six case-fatality tables from any name-pattern selection. Now
`ftlt_change_relative_`. ⚠️ **Breaking for downstream readers** of `case fatality relative
change by *.csv` written before 2026-08-05 — values are unchanged, only the column name.
Tests: `test_tbl_col_prefixes.R`.

`export_tables()` turns the summaries into CSVs. Besides the main / mortality /
disease-characteristics / exposure / cost-effectiveness (`cea = TRUE`) tables, it also
writes **equity slope-index** tables (`equity = TRUE`, default): absolute (`AEI_total`,
`AEI_per100k`) and relative (`REI_rel`, `RII_ratio`) analogues of the Slope/Relative Index
of Inequality for cumulative CPP, CYPP, DPP and net QALYs across DIMD deciles (pro-poor sign
convention; MC iterations give the uncertainty), plus `total_benefit` and `fit_R2` and an
`n_mc` column. The estimator is the grouped-data weighted least squares of Kakwani, Wagstaff
& van Doorslaer (1997); `RII_ratio` is the Kunst–Mackenbach ratio of fitted extremes, *not*
Moreno-Betancur's RII (they criticise that formula; theirs is `exp(beta)` from a log-linear
model, impossible on a signed quantity). Three things bite:
**(a)** `AEI_total = SII * n` is a *hypothetical* whole-population figure, **not** the gap
between the extreme deciles (~8x smaller) and not Cookson's published gap (`* (J-1)/J^2`);
it scales with `n`, so only `AEI_per100k`/`REI_rel` compare across localities.
**(b)** Both relative indices are `NA` unless the fitted benefit is **positive** — with a
negative mean they invert, publishing a pro-poor result as pro-rich.
**(c)** `all_diseases_sum` excludes the design's umbrella conditions (incidence `type: 0`:
`dm`, `ctdra`, `cancer`) as well as `cms*`, via `umbrella_disease_names()` — summing an
umbrella beside its own components double-counts ~9.5% of events and distorts the slope.
The ridit reference population is selectable via
`equity_ridit_reference` (`"comparator"` default vs `"scenario"` — the latter violates
Cookson's equal-group-size precondition, so it is a sensitivity analysis; see the Renard 2019
decision support in the vignette). The ridit is re-derived each **year** (N is the
comparator's reporting-year population) — deliberately, since pinning it at baseline would
manufacture a gradient out of differential population growth. Core math lives in
`calc_equity_slope_indices()` / `export_equity_tables()` in `Simulation_class_tables.R`; see
`vignette("understanding_model_outputs")`.

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
