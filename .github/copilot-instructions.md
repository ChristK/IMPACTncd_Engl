# IMPACTncd_Engl AI Agent Instructions

## Project Overview

IMPACTncd_Engl is a **microsimulation framework** for modeling the impact of public health interventions on chronic disease outcomes in England. It uses synthetic populations to simulate disease progression, evaluating "what-if" policy scenarios against a baseline.

**Core Workflow**: Baseline scenario → Policy scenario(s) → Export summaries → Analysis

## Architecture

### R6 Class Hierarchy
The simulation is built on three main R6 classes in [Rpackage/IMPACTncd_Engl_model_pkg/R/](Rpackage/IMPACTncd_Engl_model_pkg/R/):

1. **`Design`** ([Design_class.R](Rpackage/IMPACTncd_Engl_model_pkg/R/Design_class.R)): Loads and validates YAML configuration from [inputs/sim_design.yaml](inputs/sim_design.yaml)
2. **`SynthPop`** ([SynthPop_class.R](Rpackage/IMPACTncd_Engl_model_pkg/R/SynthPop_class.R)): Manages synthetic population lifecourse data (stored as `.fst` files)
3. **`Simulation`** ([Simulation_class.R](Rpackage/IMPACTncd_Engl_model_pkg/R/Simulation_class.R)): Orchestrates simulation runs, scenarios, and multicore execution

**Additional Classes**: `Disease` and `ExposureEffect` classes handle disease and exposure-specific logic.

### Data Flow
- **Inputs**: 
  - YAML config ([inputs/sim_design.yaml](inputs/sim_design.yaml))
  - `.fst` exposure distributions in [inputs/exposure_distributions/](inputs/exposure_distributions/)
  - CSVY files with relative risks in [inputs/RR/](inputs/RR/)
  - Mortality data in [inputs/mortality/](inputs/mortality/)
  - Population estimates/projections in [inputs/pop_estimates_lsoa/](inputs/pop_estimates_lsoa/) and [inputs/pop_projections/](inputs/pop_projections/)
- **Processing**: Synthetic population generated/loaded → Baseline run → Scenario modifications → Disease progression simulation
- **Outputs**: Stored in `output_dir` (configured in YAML, typically `/mnt/storage_fast/output/`) with `.fst` files for efficient read/write
  - Subdirectories: `lifecourse/`, `logs/`, `plots/`, `summaries/`, `tables/`, `xps/`

### Performance-Critical Components
- **C++/Rcpp**: Performance-critical code in [src/](Rpackage/IMPACTncd_Engl_model_pkg/src/) 
  - [simsmok.cpp](Rpackage/IMPACTncd_Engl_model_pkg/src/simsmok.cpp) - Smoking simulation
  - [IMPACTncd_sim.cpp](Rpackage/IMPACTncd_Engl_model_pkg/src/IMPACTncd_sim.cpp) - Core disease simulation engine
  - Custom GAMLSS distributions (e.g., `NBI_distribution.cpp`, `BCPEo_distribution.cpp`)
- **Multicore**: Uses `foreach` + `data.table` threading (set via `clusternumber` in YAML)
- **File-based persistence**: `.fst` format for fast columnar read/write with compression
  - SynthPop files named as `synthpop_<checksum>_<mc>.fst` with companion `_meta.yaml` files

## Development Workflow

### Initial Setup
1. Run [global.R](global.R) to initialize environment (installs dependencies, builds local R package)
2. Package rebuild happens automatically via `installLocalPackageIfChanged()` when source changes detected
3. Dependencies listed in `docker_setup/r-packages.txt` are auto-installed if missing

### Running Simulations
```r
source("global.R")
IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")
IMPACTncd$del_logs()$del_outputs()$run(1:2, multicore = TRUE, "sc0")
```

See [simulation.R](simulation.R) for baseline template, [simulate_mrtl_paper.R](simulate_mrtl_paper.R) for complex scenario examples.

### Debugging
- Use `.vscode/debug` task (runs [.vscode/debug.R](.vscode/debug.R) to setup env vars for R debugging)
- Check [DEBUG-NOTES.txt](DEBUG-NOTES.txt) for known issues and C++ debugging setup
- `validation: yes` in YAML enables additional checks (slower execution)
- For C++ call stacks: Set `PKG_CXXFLAGS= -O0 -ggdb3` in Makevars and rebuild package

### Testing
- Package tests in [Rpackage/IMPACTncd_Engl_model_pkg/tests/](Rpackage/IMPACTncd_Engl_model_pkg/tests/) use `tinytest`
- Run test suite: `tinytest::test_package("IMPACTncdEngland")`

## Key Conventions

### R Code Style
- **data.table** for all data manipulation (not dplyr/tidyverse)
- **Reference semantics**: Prefer `:=` for in-place modification, avoid unnecessary copies
- **Column naming**: Snake_case with suffixes like `_curr_xps` (current exposure), `_lagged` (lagged values)
- **Memory efficiency**: Store integers as `integer` not `double` (see [TODO.Rmd](TODO.Rmd) line 12)

### Exposure Variables (Primary Prevention)
Key variables in `synthpop$pop` for primary prevention scenarios:
- `active_days_curr_xps` - Physical activity days per week
- `fruit_curr_xps`, `veg_curr_xps` - Fruit/veg consumption (g/day, 80g = 1 portion)
- `smok_status_curr_xps` - Smoking status (1=never, 2=occasional, 3=ex, 4=current)
- `smok_cig_curr_xps` - Cigarettes per day
- `ets_curr_xps` - Second-hand smoke exposure (0/1)
- `alcohol_curr_xps` - Alcohol consumption (g/day)
- `bmi_curr_xps` - Body mass index (kg/m²)
- `sbp_curr_xps` - Systolic blood pressure (mmHg)
- `tchol_curr_xps` - Total cholesterol (mmol/L)

### Disease Variables (Secondary Prevention)
Key disease status columns for secondary prevention:
- `dm`, `t1dm`, `t2dm` - Diabetes (any, type 1, type 2)
- `chd`, `stroke`, `af`, `hf` - Cardiovascular diseases
- `ckd` - Chronic kidney disease
- `colorect_ca`, `prostate_ca` - Cancers
- See [how_to_run_scenarios.Rmd](Rpackage/IMPACTncd_Engl_model_pkg/vignettes/how_to_run_scenarios.Rmd) for complete list

### Scenario Definition Pattern
Scenarios modify populations via `update_primary_prevention_scn()` or `update_secondary_prevention_scn()`:
```r
IMPACTncd$update_primary_prevention_scn(
  function(synthpop) {
    change_10pc <- 0.1
    sc_year <- 23L
    synthpop$pop[year >= sc_year, bmi_curr_xps := bmi_curr_xps * (1 - change_10pc)]
  }
)
```
**Pattern**: Define function that modifies `synthpop$pop` data.table in-place. See examples in [simulate_mrtl_paper.R](simulate_mrtl_paper.R).

### File Paths
- **Absolute paths required** for output/synthpop directories (e.g., `/mnt/storage_fast/`)
- **Relative paths** for inputs from project root (e.g., `./inputs/exposure_distributions/`)
- **External storage**: `/mnt/storage_fast/` and `/mnt/storage_fast4/` for large output files

### YAML Configuration
- Always copy [inputs/sim_design.yaml](inputs/sim_design.yaml) before modifying (prevents merge conflicts)
- Critical parameters: `n` (population size), `clusternumber` (parallelism), `output_dir`, `synthpop_dir`
- Disease definitions in `diseases:` array with `name`, `meta`, `rr_dependencies` fields

## Cross-Project Integration

- **CKutils**: Custom R package ([/home/ckyprid/GH_projects/CKutils](../../../GH_projects/CKutils)) with helper functions (installed from GitHub)
- **IMPACTncd_Japan**: Sister project ([/home/ckyprid/My_Models/IMPACTncd_Japan](../IMPACTncd_Japan)) with similar architecture but Japan-specific inputs
- Shared patterns: R6 class structure, `.fst` persistence, YAML configuration

## Documentation

- **Vignettes**: See [Rpackage/IMPACTncd_Engl_model_pkg/vignettes/](Rpackage/IMPACTncd_Engl_model_pkg/vignettes/) for:
  - [how_to_test_run.Rmd](Rpackage/IMPACTncd_Engl_model_pkg/vignettes/how_to_test_run.Rmd) - Quick start guide
  - [how_to_run_scenarios.Rmd](Rpackage/IMPACTncd_Engl_model_pkg/vignettes/how_to_run_scenarios.Rmd) - Scenario setup
- **roxygen2**: All R6 methods documented with roxygen2 comments; rebuild with `roxygen2::roxygenise()`

## Common Pitfalls

1. **Windows .Random.seed error**: Run `invisible(runif(1))` before creating Simulation object (see [global.R](global.R) line 93)
2. **Docker environment detection**: Code checks `file.exists("/.dockerenv")` to skip CRAN mirror selection in containers
3. **Seed management**: Each SynthPop has deterministic seeds based on `mc` (Monte Carlo iteration) parameter
4. **File overwrites**: Set `simulation_files_overwrite: yes` in YAML to overwrite existing outputs (default: no)

## Asset Management

- GitHub releases via `piggyback`: See [gh_upload.R](gh_upload.R) for uploading assets, [ghDeleteSurplusAssets.R](ghDeleteSurplusAssets.R) for cleanup
- Assets deployed automatically via `DeployGitHubAssets()` in Simulation initialization
