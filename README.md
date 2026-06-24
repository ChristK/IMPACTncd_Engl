# IMPACT<sub>NCD_Engl</sub> microsimulation

--------------------------------------------------------------------------------

IMPACT<sub>NCD_Engl</sub> is an implementation of the IMPACT<sub>NCD</sub> framework, developed by Chris
Kypridemos with contributions from Max Birkett, Karl Emmert-Fees, Anna Head, Brendan Collins, Martin O'Flaherty,
Peter Crowther (Melandra Ltd), Maria Guzman-Castillo, Amandine Robert, Piotr Bandosz, and Adithi R. Upadhya.

Several research grants have supported its development including grants from the Health Foundation,
NIHR, EU Horizon2020, Liverpool City Council, MRC, NIH, and the National Cerebral and Cardiovascular Center in Japan.

Copyright (C) 2018-2025 University of Liverpool, Chris Kypridemos.

IMPACT<sub>NCD_Engl</sub> is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version. This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along
with this program; if not, see <http://www.gnu.org/licenses/> or write to the
Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA.


## Overview

This microsimulation translates changes in the trends of established chronic disease risk factors to changes in the incidence and prevalence of common chronic
diseases, and changes in all-cause mortality, in England. A typical simulation consists of a *baseline* and one or more *what-if* policy scenarios. Policy scenarios usually comprise of changes in the exposure of a population segment to a risk factor, e.g. after some future health policy, 10% of the population show

- improved body mass index (BMI), blood pressure, or cholesterol level (each by 5%).
- increased fruit and vegetable intake (by 5%).
- reduced alcohol intake (by 5%).
- reduced smoking prevalence or passive smoking (each by 5%).
- increased physical activity (by having an additional active day).

The microsimulation can then estimate the effectiveness and equity of the modelled health policy.

## Installation

The model has two parts that are obtained separately:

1. **The code** — cloned (or downloaded) from the public [GitHub repository](https://github.com/ChristK/IMPACTncd_Engl). No GitHub account or access token is required.
2. **The data** — the large input data and pre-computed PARFs/RR tables are **not** in git. They are published on Zenodo ([concept DOI `10.5281/zenodo.20812409`](https://doi.org/10.5281/zenodo.20812409), CC-BY-SA-4.0) and downloaded with one command. The data is public, so **no Zenodo account or token is required**.

There are two ways to install and run the model: a **Docker** image (recommended — everything is bundled and reproducible) or a **native R** installation.

### Option A — Docker (recommended)

Pre-built Docker images bundle the full R environment, all model data (from Zenodo), and the model code, so there is nothing to compile or download by hand. After installing Docker:

```bash
cd docker_setup
./setup_user_docker_env.sh        # Linux/macOS   (Windows PowerShell: setup_user_docker_env.ps1)
```

See **[docker_setup/README.md](docker_setup/README.md)** for the full Docker guide: installing Docker, the `setup_user_docker_env`/`setup_dev_docker_env` helpers, the three-layer image architecture (prerequisite → data → model), and how to build the images locally.

### Option B — Native R installation

Run the model directly in R on your machine.

**Prerequisites** — install:

- [R](https://cran.r-project.org/) (≥ 4.1.0; the model is developed and tested against R 4.6.0) — guides for [Windows](installation_docs/installing_R_on_windows.md) and [Linux](installation_docs/installing_R_on_linux.md).
- [Rtools](installation_docs/installing_rtools_on_windows.md) — **Windows only**, required to compile the C++ simulation engine.
- [Git](installation_docs/installing_git.md) — recommended (eases future updates), and optionally [RStudio](installation_docs/installing_RStudio_on_windows.md).

Additional R packages (listed in [`docker_setup/r-packages.txt`](docker_setup/r-packages.txt), e.g. `data.table`, `foreach`) are installed automatically the first time you `source("global.R")`.

**1. Get the code.** Clone the public repository:

```bash
git clone https://github.com/ChristK/IMPACTncd_Engl.git
```

Or download a ZIP: on the [repository page](https://github.com/ChristK/IMPACTncd_Engl), click **Code ▸ Download ZIP**, then unzip it to your chosen location.

**2. Download the model data from Zenodo** (required on a fresh clone — a freshly cloned repository has no data). The data is public, so no account or token is needed:

```r
source("global.R")
IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")
IMPACTncd$zenodo_connect()        # defaults to the published record, anonymous
IMPACTncd$zenodo_download_all()   # inputs + pre-computed PARFs/RR tables (~13 GB)
```

See the **`zenodo_data_management`** vignette for selective downloads, updating to a new data version, and (for data managers) publishing new data.

> **Tip:** It is highly recommended to create your own copies of `inputs/sim_design.yaml` and `simulation.R`, edit those, and leave the originals untouched. This minimises merge conflicts when pulling updates from GitHub.

## Documentation and tutorials

The package ships tutorial vignettes. List them all with:

```r
vignette(package = "IMPACTncdEngland")
```

and open one with `vignette("<name>", package = "IMPACTncdEngland")`. The source `.Rmd` files live in [`Rpackage/IMPACTncd_England_model_pkg/vignettes/`](Rpackage/IMPACTncd_England_model_pkg/vignettes).

| Vignette | Topic |
|---|---|
| `how_to_test_run` | **Start here** — run a test simulation (also covers how to uninstall the package). |
| `how_to_run_scenarios` | Define and run a baseline plus what-if policy scenarios. |
| `understanding_model_outputs` | Interpret the model's output files. |
| `custom-scenario-columns` | Export custom columns created within scenarios. |
| `inputs_manifest_system` | The inputs manifest / data-asset tracking system. |
| `zenodo_data_management` | Download, manage, and publish the model's input data on Zenodo. |

Other help files:

- **[docker_setup/README.md](docker_setup/README.md)** — Docker images, the three-layer build, and user/developer environment setup.
- **[installation_docs/](installation_docs)** — step-by-step guides for installing R, Rtools, Git, and RStudio.
- **CLAUDE.md** — architecture notes and common commands for the repository.

## Further notes and references

You can find a detailed technical description and validation of the model [here](https://www.health.org.uk/sites/default/files/2023-07/REAL_Insights_Technical%20appendix.pdf).

This version of the model has been used to produce two REAL Centre reports, named [Health in 2040: projected patterns of illness in England](https://www.health.org.uk/publications/health-in-2040) and [Health inequalities in 2040: current and projected patterns of illness by level of deprivation in England](https://www.health.org.uk/publications/reports/health-inequalities-in-2040).

Previous versions of the IMPACT<sub>NCD</sub> microsimulation are available in [Chris Kypridemos Github repository](https://github.com/ChristK) and have been published in various academic journals, such us [BMJ](https://doi.org/10.1136/bmj.i2793), [PLOS Med](https://doi.org/10.1371/journal.pmed.1002573), and [Circulation](https://doi.org/10.1161/CIRCULATIONAHA.118.036751).
