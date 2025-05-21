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

## Requirements

The IMPACT<sub>NCD_Engl</sub> distribution is usually installed directly from [GitHub](https://github.com/ChristK/IMPACTncd_Engl/) via the [Git](https://git-scm.com/) version control system, which should be installed on the target system. While not absolutely required, Git eases installation of future updates and previous releases. While earlier versions will suffice, Git 2.40.0 was the latest at release time.

The IMPACT<sub>NCD_Engl</sub> model is written primarily in, and so requires installation of, the [R programming language](https://cran.r-project.org/), for which version 4.2.3 was the latest at release time. Additional R packages (listed in `dependencies.yaml`) will be installed automatically if missing at execution time, e.g. `data.table`, `piggyback`, `foreach`.

## Installation

IMPACT<sub>NCD_Engl</sub> is installed directly from GitHub, after which a model-specific configuration is set. With multiple configurations, different policy scenarios may be tested on the same server. 

As the installation requires downloading numerous files from GitHub, you will need a Personal Access Token. Instructions [here](https://github.com/ChristK/IMPACTncd_Engl/blob/main/installation_docs/installing_git.md#setting-a-personal-access-token) 

### Download without GitHub

##### To install and run IMPACT<sub>NCD_Engl</sub> we need to install R on [Windows](installation_docs/installing_R_on_windows.md) / [Linux](installation_docs/installing_R_on_linux.md); and [rtools](installation_docs/installing_rtools_on_windows.md) (only for Windows) before these following steps

1. To use IMPACTncd England navigate here: https://github.com/ChristK/IMPACTncd_Engl 

2. Select Code on the right side corner, this will show an option of **Download ZIP**

![](installation_docs/img/Github_4.jpeg)

3. Once downloaded unzip the folder and save it in desired location


### Download with GitHub

##### To install and run IMPACT<sub>NCD_Engl</sub> we need to install R on [Windows](installation_docs/installing_R_on_windows.md) / [Linux](installation_docs/installing_R_on_linux.md); [Git](installation_docs/installing_git.md); and [rtools](installation_docs/installing_rtools_on_windows.md) (only for Windows) before these following steps

1. To use IMPACT<sub>NCD_Engl</sub> navigate here: https://github.com/ChristK/IMPACTncd_Engl 

2. Select Code on the right side corner, this will present a link which can be used to download the package in local system. Copy this link

![](installation_docs/img/Github_4.jpeg)


3. In terminal, type git clone https://github.com/ChristK/IMPACTncd_Engl.git (which is the copied link) 


### Modelling with IMPACT<sub>NCD_Engl</sub>


##### Important Note : It is highly recommended to create your own copies of `sim_design.yaml` present in `inputs` folder and `simulation.R` in the root folder to run this model and make no changes to the original `sim_design.yaml` and `simulation.R`, as this would help reducing conflicts while pulling changes from [GitHub](https://github.com/ChristK/IMPACTncd_Engl). 

1. All vignettes can be viewed using the code below 

```{r}
vignette(package = "IMPACTncdEngl")
```

2. The next steps can be found in the vignette which can be accessed using the code below

```{r}
vignette("how_to_test_run", package = "IMPACTncdEngl")
```

3. Use the following code to open the vignette to run different policy scenarios

```{r}
vignette("how_to_run_scenarios", package = "IMPACTncdEngl")
```

4. Use the following code to open the vignette to understand model outputs 

```{r}
vignette("understanding_model_outputs", package = "IMPACTncdEngl")
```

5. To completely delete the package IMPACTncdEngl is explained in a section called **How to remove IMPACTncdEngl installed package** in the vignette mentioned below

```{r}
vignette("how_to_test_run", package = "IMPACTncdEngl")
```

## Further notes and references

You can find a detailed technical description and validation of the model [here](https://www.health.org.uk/sites/default/files/2023-07/REAL_Insights_Technical%20appendix.pdf).

This version of the model has been used to produce two REAL Centre reports, named [Health in 2040: projected patterns of illness in England](https://www.health.org.uk/publications/health-in-2040) and [Health inequalities in 2040: current and projected patterns of illness by level of deprivation in England](https://www.health.org.uk/publications/reports/health-inequalities-in-2040).

Previous versions of the IMPACT<sub>NCD</sub> microsimulation are available in [Chris Kypridemos Github repository](https://github.com/ChristK) and have been published in various academic journals, such us [BMJ](https://doi.org/10.1136/bmj.i2793), [PLOS Med](https://doi.org/10.1371/journal.pmed.1002573), and [Circulation](https://doi.org/10.1161/CIRCULATIONAHA.118.036751).
