# IMPACTncd_Engl microsimulation

--------------------------------------------------------------------------------

IMPACTncd_Engl is an implementation of the IMPACTncd framework, developed by Chris
Kypridemos with contributions from Max Birkett, Karl Emmert-Fees, Anna Head, Brendan Collins, Martin O'Flaherty,
Peter Crowther (Melandra Ltd), Maria Guzman-Castillo, Amandine Robert, and Piotr Bandosz. 

Several research grants have supported its development including grants from the Health Foundation,
NIHR, EU Horizon2020, Liverpool City Council, MRC, NIH, and the National Cerebral and Cardiovascular Center in Japan.  

Copyright (C) 2018-2023 University of Liverpool, Chris Kypridemos.

IMPACTncd_Engl is free software; you can redistribute it and/or modify it under the
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

A simulation consists of *baseline* and *what-if* policy scenarios (or hypotheses). Policy scenarios typically comprise populations with *improved* risk factors, e.g. after some future health policy, 20% of the population show 

- improved body mass index (BMI), blood pressure, or cholesterol level (each by 20%).
- increased fruit and vegetable intake (by 20%).
- reduced alcohol intake (by 20%).
- reduced smoking prevalence or passive smoking (each by 20%).
- increased physical activity (by a single day).

Output data is then analysed for health improvements, e.g. in disease prevalence, incidence, or mortality. 

## Requirements

The IMPACTncd_Engl distribution is usually installed directly from [GitHub](https://github.com/ChristK/IMPACTncd_Engl/) via the [Git](https://git-scm.com/) version control system, which should be installed on the target system. While not absolutely required, Git eases installation of future updates and previous releases. While earlier versions will suffice, Git 2.40.0 was the latest at release time.

The IMPACTncd_Engl model is written primarily in, and so requires installation of, the [R programming language](https://cran.r-project.org/), for which version 4.2.3 was the latest at release time. Additional R packages (listed in `dependencies.yaml`) will be installed automatically if missing at execution time, e.g. `data.table`, `piggyback`, `foreach`.

## Installation

IMPACTncd_Engl is installed directly from GitHub, after which a model-specific configuration is set. With multiple configurations, different policy scenarios may be tested on the same server. The following notes describe installation and execution of a test model:

1. create a folder in which to install your IMPACTncd_Engl distribution and change the current directory to this folder, e.g. from a Linux terminal:

	```
	$ mkdir IMPACTncdEngl
	$ cd IMPACTncdEngl
	```
  
2. install the IMPACTncd_Engl repository ('repo') from GitHub. As IMPACTncd_Engl is a publicly visible repo, no GitHub account is needed to access the source code or other repo assets:

	```$ git clone https://github.com/ChristK/IMPACTncd_Engl.git```

	The above installs the most recent IMPACTncd_Engl development version. If a particular earlier version is required, perhaps a 
	specific release ('tag') or code change ('commit'), select this now with:
	
	```$ git checkout <release-or-commit-id>```	
	
	where `<release-or-commit-id>` is a particular release or commit ID; `git tag` or `git log` will show a complete list of tag or commit options. After examining this earlier version, return to the most recent code with:
	
	```$ git checkout main```	
	
3. now to install certain large data files, stored as GitHub repo *assets*. A config file `auxil/ghAssetConfig.yaml` holds details which ease asset installation. In here, set [`id`] to a suitable name for this particular IMPACTncd_Engl simulation, and [`uploadSrcDirectory`] and [`deployToRootDirectory`] to a desired target IMPACTncd_Engl folder path/name. Download the first five assets with (from a Linux terminal):

	```
	$ Rscript gh_deploy.R auxil/ghAssetConfig.yaml <simulation-id>
	```
	
	or, on Windows:
	
	```
	<R-installation-directory>\bin\Rscript gh_deploy.R auxil/ghAssetConfig.yaml <simulation-id>
	```  
	
	where `<R-installation-directory>` is likely to be something like `C:\Program Files\R\R-4.2.0`, `<simulation-id>` is your name (from `auxil/ghAssetConfig.yaml`) for this IMPACTncd_Engl simulation. Alternatively, if wanting to avoid the Linux terminal, assets may be downloaded from within RStudio as follows:
	
	```
	setwd("<path-to-your-IMPACTncd_Engl-distribution>")
	source("gh_deploy.R")
	DeployGitHubAssets(sAssetConfigFilePath="auxil/ghAssetConfig.yaml",sGitHubAssetRouteId="<simulation-id>")
	```
	
	where `<path-to-your-IMPACTncd_Engl-distribution>` is the directory in which the IMPACTncd_Engl simulation has been installed, something like `C:/.../IMPACTncdEngl/IMPACTncd_Engl` on Windows.

4. after this, the first five assets in `auxil/filindx.csv` should have downloaded:

	```
	inputs/pop_estimates_lsoa/LSOA_1st_April_population_estimates.fst
	inputs/pop_estimates_lsoa/LSOA_mid_year_population_estimates.fst
	simulation/parf/PARF_af_9dceba2c7567c3e48a443acdffbe19d8.fst
	simulation/parf/PARF_andep_33c23b1449ddd6fb61bad0aa8cfe135c.fst
	simulation/parf/PARF_asthma_5ba32a1a5d21f638e4ab8afda42befc3.fst
	```

5. if the last step was successful, the full set of 631 assets listed in `auxil/filindx.csv` (as of 2023-03-27) should now be downloaded. First edit `auxil/ghAssetConfig.yaml` and comment-out the [`testWithFirstNAssets`] line. Now save `ghAssetConfig.yaml` and re-run the above step 3.

6. finally, a model and configuration is set. The model defines the *baseline* and *what-if* scenarios, while the configuration specifies other computational and model settings, e.g. output directories, parallelisation, population and disease parameters (see the **Model design and execution** section below). Using the default model, set the below configuration changes:

	a. in `auxil/simulation.R`, set the current `run(1:200,` command (lines 16 and 132), which describe a number of repetitions, to something appropriate - perhaps `run(1:2,` for testing.
	
	b. in `auxil/sim_design_parf.yaml` and `inputs/sim_design.yaml`, set:

		* paths to pre-existing output directories for [`output_dir`] and [`synthpop_dir`], as otherwise any prior data here may be deleted. 
		* the [`validation`] value so inactive, i.e. `validation: no`.
		* the [`clusternumber`] parameter to an appropriate number of cores on which to run the model in parallel, e.g. perhaps only 2-4 cores for either a personal laptop or if running on a shared server. The default 15 cores will impact other users using the same shared server.
		
7. the model is now ready for execution; see the execution notes in the following **Model design and execution** section.

## Model design and execution

Presently, `auxil/simulation.R` tool both prepares a default model,[^policyScenario] and then executes the *baseline* and *what-if* scenarios after loading a configuration set in `inputs/sim_design.yaml`. The what-if dataset is created by code in `auxil/simulation.R:scenario_fn()`, which enables evaluation of detailed models. Multiple models may be examined from one IMPACTncd_Engl directory by creating copies of `auxil/simulation.R` and `inputs/sim_design.yaml` - each with an appropriate `scenario_fn()` and relevant parameters. In a subsequent release, a tool will create a range of what-if datasets without such source code changes.

Due to RStudio parallelisation issues the simulation is best executed from within a Linux terminal. Within an R session, models may be executed with `source("auxil/simulation.R")`. Alternatively, directly within a Linux terminal:

```
$ Rscript auxil/simulation.R > myLogFile230321a.txt 2>&1 &
```

where the optional piping (`> myLogFile230321a.txt 2>&1 &`) commands direct the *stdout* and *stderr* outputs to a `myLogFile230321a.txt` log file, useful for examining behaviour later, while the final `&` runs the simulation as a background task, so the terminal remains usable. As models may take many hours to finish, depending on their configuration and the host machine's capacity, it may be preferred to run these on a remote server within a window manager such as $screen$, so that users can disconnect from the terminal and reconnect at a convenient time later.
	
## Further notes and references

[^policyScenario]: the default what-if scenario is described further elsewhere.
