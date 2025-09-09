# IMPACTncd_Engl prerequisite dockerfile

This repository contains the Dockerfile used to build the prerequisite container for the IMPACTncd model. The container is based on Ubuntu and includes R version 4.4.3 with package versions frozen as of 31/03/2025, using the [RStudio Package Manager](https://packagemanager.posit.co/client/#/). You can update the required R packages by editing the file `r-packages.txt` and running the build script. Similarly, you can update the required operating system libraries by editing the file `apt-packages.txt` and running the build script. You can check the system prerequisites for Ubuntu 24.04 (Noble)  for any additional package you may want to add `r-packages.txt` [here](https://packagemanager.posit.co/client/#/repos/cran/packages/overview?search=) and then if necessary add the required system library to `apt-packages.txt`. To find the version of the system library, for full reproducibility, you can use the command `apt-cache policy <package_name>` from a terminal within the container.

This Docker container is required for the branch 'master' of the IMPACTncd_Engl model.

## 🐳 Docker Setup for IMPACTncd Engl

This directory contains the Docker configuration and scripts needed to build and run a containerized environment for the **IMPACTncd Engl** project.

The main scripts, [`create_env.sh`](./create_env.sh) and [`create_env.ps1`](./create_env.ps1), build a Docker image with the necessary system and R packages, and run a container with your project directory and additional folders (defined in `sim_design.yaml`) mounted inside the container for immediate use.

---

### 🚀 Quick Start (Linux/macOS)

```bash
cd docker_setup
./create_env.sh [optional_path_to_sim_design.yaml]
```

If no YAML file is specified, the default path `../inputs/sim_design.yaml` will be used.

> On macOS, install `coreutils` first (if needed):
> ```bash
> brew install coreutils
> ```

---

### 🪟 Quick Start (Windows - PowerShell)

Use the PowerShell script [`create_env.ps1`](./create_env.ps1):

```powershell
cd docker_setup
.\create_env.ps1 [optional_path_to_sim_design.yaml]
```

If no path is specified, the script uses the default: `..\inputs\sim_design.yaml`.

> If you encounter an execution policy error, run:
> ```powershell
> Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
> ```

---

## 📦 What the `create_env.*` Scripts Do

1. Accept an optional path to a `sim_design.yaml` file.
2. Extract the following from the YAML file:
   - `output_dir`: a local directory to mount
   - `synthpop_dir`: another directory to mount
3. Compute a hash from:
   - `Dockerfile.prerequisite`
   - `apt-packages.txt`
   - `r-packages.txt`
4. Compare the hash to a saved version in `.docker_build_hash`.
5. Rebuild the Docker image **only if** any of those files changed.
6. Mount the following directories into the Docker container:
   - The **project root** to `/IMPACTncd_Engl`
   - The `output_dir` to `/IMPACTncd_Engl/outputs`
   - The `synthpop_dir` to `/IMPACTncd_Engl/synthpop`
7. Launch the container with an interactive shell.

---

## 🔍 Directory Mounting Summary

| Host Path (automatically resolved) | Mounted to inside container              |
|-----------------------------------|------------------------------------------|
| Project root (`../`)              | `/IMPACTncd_Engl`                       |
| `output_dir` from YAML            | `/IMPACTncd_Engl/outputs`              |
| `synthpop_dir` from YAML          | `/IMPACTncd_Engl/synthpop`             |

This allows you to run R scripts inside the container like:

```bash
Rscript /IMPACTncd_Engl/scripts/my_analysis.R
```

---

## 🐳 Docker Image

- **Image name:** `impactncd-engl-r-prerequisite:latest`
- **Base image:** [`rocker/r-ver`](https://hub.docker.com/r/rocker/r-ver)
- **System packages:** listed in [`apt-packages.txt`](./apt-packages.txt)
- **R packages:** listed in [`r-packages.txt`](./r-packages.txt)

---

## 🧼 Clean-Up

To remove the Docker image:

```bash
docker rmi impactncd-engl-r-prerequisite:latest
```

To prune unused containers/images:

```bash
docker system prune
```

---

## 🛠 Build and Push Using Provided Scripts

You can use the provided scripts:

- **Prerequisite Container:**
  - On Linux or macOS: `./build_and_push_prerequisite.sh`
    - Use the `--push` argument to push the Docker image to Docker Hub after building.
  - On Windows: `build_and_push_prerequisite.ps1`
    - Use the `-Push` argument to push the Docker image to Docker Hub after building.

- **IMPACTncd Container:**
  - On Linux or macOS: `./build_and_push_IMPACTncd.sh`
    - Use the `--push` argument to push the Docker image to Docker Hub after building.
  - On Windows: `build_and_push_IMPACTncd.ps1`
    - Use the `-Push` argument to push the Docker image to Docker Hub after building.

---

## ❓ Troubleshooting

- **Docker not found?** Make sure Docker Desktop is running.
- **On macOS:** You may need `coreutils` (for `gsha256sum`).
- **On Windows:** Use PowerShell, not Command Prompt.

---

## 📬 Need Help?

If you encounter issues, please contact the project maintainers or raise an issue in this repository.
