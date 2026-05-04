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

### Running Simulations
```r
source("global.R")
IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")
IMPACTncd$del_logs()$del_outputs()$run(1:2, multicore = TRUE, "sc0")$export_summaries(multicore = TRUE)
```

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

### File Paths
- **Absolute paths required** for `output_dir`, `synthpop_dir` in YAML
- **Relative paths** for inputs from project root (e.g., `./inputs/`)
- External storage typically at `/mnt/storage_fast/`

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
