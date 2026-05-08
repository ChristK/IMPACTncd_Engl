# IMPACTncd_England Docker Setup

Cross-platform Docker configuration and setup scripts for the **IMPACTncd England** microsimulation project. Two entry points cover the common workflows:

- **`setup_user_docker_env.{sh,ps1}`** — pulls a pre-built image from Docker Hub and runs the model. Most users want this.
- **`setup_dev_docker_env.{sh,ps1}`** — builds the prerequisite image locally from `Dockerfile.prerequisite.IMPACTncdENGL` and mounts the project source into the container. Use this when you're modifying the model code or its dependency lists.

---

## 📋 Table of Contents

- [Prerequisites](#-prerequisites)
- [Quick Start — running simulations](#-quick-start--running-simulations)
- [Mount points](#-mount-points)
- [Docker images](#-docker-images)
- [Troubleshooting](#-troubleshooting)
- [Cleanup](#-cleanup)
- [Developer documentation](#-developer-documentation)

---

## 💾 Prerequisites

### Install Docker

- **Windows:** [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/)
- **macOS:** [Docker Desktop for Mac](https://docs.docker.com/desktop/install/mac-install/)
- **Linux:** [Ubuntu](https://docs.docker.com/engine/install/ubuntu/) · [Debian](https://docs.docker.com/engine/install/debian/) · [Fedora](https://docs.docker.com/engine/install/fedora/) · [CentOS / Rocky Linux](https://docs.docker.com/engine/install/centos/)

### Platform-specific requirements

**Windows**
- PowerShell 5.1+ (preinstalled on Windows 10/11).
- WSL2 backend strongly recommended for Docker Desktop.
- If your YAML config uses Windows-style absolute paths (e.g. `output_dir: C:/data/...`), the scripts auto-translate them to the correct WSL mount path via `wsl.exe wslpath -u`. This works correctly even when `/etc/wsl.conf` sets a non-default `[automount] root` (e.g. `/mnt/host/`).

**macOS**
```bash
brew install coreutils   # provides gsha256sum
```

**Linux**
```bash
sudo usermod -aG docker $USER   # then log out and back in
```

---

## 🚀 Quick Start — running simulations

Both scripts accept the same parameters. Note that `--UseVolumes` is a double-dash flag in bash but `-UseVolumes` in PowerShell — that's a PowerShell convention, not a typo.

### Linux / macOS (Bash)

The bash script uses **flag-style** arguments (`-Tag value`), not positional.

```bash
cd docker_setup

# Default: pulls chriskypri/impactncdengl:main from Docker Hub
./setup_user_docker_env.sh

# Specific Docker Hub tag
./setup_user_docker_env.sh -Tag v1.0.0

# Locally built image (build first — see "Building images locally" below)
./setup_user_docker_env.sh -Tag local

# Custom YAML
./setup_user_docker_env.sh -SimDesignYaml /path/to/my_sim_design.yaml

# Mount a scenarios directory at /IMPACTncd_England/scenarios in the container
./setup_user_docker_env.sh -ScenariosDir /path/to/my_scenarios

# Use Docker volumes instead of bind mounts (recommended on Windows / macOS)
./setup_user_docker_env.sh --UseVolumes

# Combine options
./setup_user_docker_env.sh -Tag v1.0.0 -ScenariosDir /path/to/scenarios --UseVolumes
```

### Windows (PowerShell)

```powershell
cd docker_setup

# Default: pulls chriskypri/impactncdengl:main from Docker Hub
.\setup_user_docker_env.ps1

# Specific Docker Hub tag
.\setup_user_docker_env.ps1 -Tag "v1.0.0"

# Locally built image
.\setup_user_docker_env.ps1 -Tag "local"

# Custom YAML
.\setup_user_docker_env.ps1 -SimDesignYaml "C:\path\to\my_sim_design.yaml"

# Mount a scenarios directory at /IMPACTncd_England/scenarios in the container
.\setup_user_docker_env.ps1 -ScenariosDir "C:\path\to\my_scenarios"

# Use Docker volumes
.\setup_user_docker_env.ps1 -UseVolumes

# Combine options
.\setup_user_docker_env.ps1 -Tag "v1.0.0" -ScenariosDir "..\scenarios" -UseVolumes
```

### Parameters

| PowerShell | Bash | Description |
|---|---|---|
| `-Tag <name>` | `-Tag <name>` | Image tag. Default: `main`. |
| `-ScenariosDir <path>` | `-ScenariosDir <path>` | Optional. Mounted at `/IMPACTncd_England/scenarios`. |
| `-SimDesignYaml <path>` | `-SimDesignYaml <path>` | Path to YAML config. Default: `../inputs/sim_design.yaml`. |
| `-UseVolumes` | `--UseVolumes` | Use Docker volumes instead of bind mounts. |

### What happens when you run the script

1. **Pulls the Docker image** from Docker Hub (if not already cached).
2. **Reads your `sim_design.yaml`** to find `output_dir` and `synthpop_dir`.
3. **Creates those host directories** if they don't exist.
4. **Starts an interactive container** with those directories mounted.
5. **Runs as your user** (UID/GID auto-detected) — not root — to avoid permission issues.

Once inside the container:
```r
source("global.R")
IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")
IMPACTncd$run(1:10, multicore = TRUE, "sc0")$export_summaries(multicore = TRUE)
```

---

## 📁 Mount points

| Host | Container | Description |
|---|---|---|
| (built into image) | `/IMPACTncd_England` | Project source. The image already contains it; no host project mount is needed. |
| `output_dir` from YAML | `/outputs` | Simulation outputs (lifecourse, summaries, etc.). |
| `synthpop_dir` from YAML | `/synthpop` | Synthetic population cache. |
| `-ScenariosDir` (optional) | `/IMPACTncd_England/scenarios` | Custom scenario scripts. |

### Bind mount mode (default)
Host directories are mounted directly. Real-time visibility, lower overhead. Recommended on Linux.

### Volume mode (`--UseVolumes` / `-UseVolumes`)
The script creates Docker-managed volumes for the output and synthpop directories, runs the container against them, then rsyncs results back to the host directories on exit and removes the volumes. Recommended on Windows and macOS where bind mounts have higher I/O overhead.

---

## 🐳 Docker images

### Image selection logic

| `-Tag` value | Image used | Source |
|---|---|---|
| `main` (default) | `chriskypri/impactncdengl:main` | Docker Hub |
| `local` | `impactncdengl:local` | Local Docker registry — build first (see [Developer documentation](#-developer-documentation)) |
| anything else, e.g. `v1.0.0` | `chriskypri/impactncdengl:<tag>` | Docker Hub |

### Available remote tags

https://hub.docker.com/r/chriskypri/impactncdengl/tags

---

## ❓ Troubleshooting

**"Cannot connect to the Docker daemon"**
- Windows / macOS: ensure Docker Desktop is running (check the system tray).
- Linux: `sudo systemctl start docker`.
- Verify with `docker info`.

**"Failed to pull Docker image"**
- Confirm the tag exists at https://hub.docker.com/r/chriskypri/impactncdengl/tags.
- For private repositories: `docker login` first.
- To verify a *local* image exists, use `docker image inspect impactncdengl:local` — `docker pull` only works for remote images.

**Linux — "permission denied" on the Docker socket**
```bash
sudo usermod -aG docker $USER   # then log out and back in
groups                          # should show 'docker'
```

**Windows — "Execution policy" error**
```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

**Windows + WSL — bind mount paths**
The scripts auto-translate Windows-style absolute paths (`C:/...`) via `wsl.exe wslpath -u`. If you see a warning that wsl.exe wasn't found and the script is falling back to a legacy `/c/...` form, your daemon is likely Hyper-V-backed or running on a host without WSL — either install WSL or switch to POSIX paths in your YAML.

**Container runs out of memory**
Edit your `sim_design.yaml` and reduce:
- `clusternumber` — fewer parallel cores, less RAM (~10 GB per core).
- `n` and `num_chunks` — smaller synthetic population.

**Specific apt versions failing**
The intelligent installer in `Dockerfile.prerequisite.IMPACTncdENGL` automatically substitutes the closest available version and prints a `PACKAGE VERSION UPDATES DETECTED` block at the end of the build. Run `update-apt-packages.{sh,ps1} --interactive` to fold those substitutions back into `apt-packages.txt`.

---

## 🧹 Cleanup

```bash
# User-facing images
docker rmi chriskypri/impactncdengl:main impactncdengl:local

# Dev-facing image (built by setup_dev_docker_env.{sh,ps1})
docker rmi prerequisite.impactncdengl:local

# Remove all unused images and stopped containers
docker system prune

# Remove everything including volumes (⚠️ destructive)
docker system prune -a --volumes
```

---

## 📬 Need help?

- Check the [main project README](../README.md).
- Review the [vignettes](../Rpackage/IMPACTncd_England_model_pkg/vignettes/) for simulation guidance.
- Open an issue on GitHub or contact the project maintainers.

---

# 🔧 Developer documentation

The following sections are for developers who need to build Docker images or maintain the package system.

---

## 🛠 Building images locally

### Developer workflow — `setup_dev_docker_env.{sh,ps1}`

This script builds `prerequisite.impactncdengl:local` from `Dockerfile.prerequisite.IMPACTncdENGL` and runs the container with your project source mounted at `/IMPACTncd_England`. It auto-rebuilds the image when any of `Dockerfile.prerequisite.IMPACTncdENGL`, `apt-packages.txt`, `r-packages.txt`, or `entrypoint.sh` changes (tracked via a hash file in this directory).

```bash
# Linux/macOS
./setup_dev_docker_env.sh                              # bind-mount mode
./setup_dev_docker_env.sh ../inputs/sim_design_test.yaml
./setup_dev_docker_env.sh --UseVolumes

# Windows
.\setup_dev_docker_env.ps1
.\setup_dev_docker_env.ps1 -SimDesignYaml "..\inputs\sim_design_test.yaml"
.\setup_dev_docker_env.ps1 -UseVolumes
```

The dev script does not accept `-Tag` — its image name is fixed.

### Building images with `docker_build_push.{sh,ps1}`

The project uses a layered architecture:

1. **Prerequisite image** (`prerequisite.impactncdengl:local`) — base R environment with all dependencies.
2. **Main image** (`impactncdengl:local`) — prerequisite + project code + input data (~20 GB).

```bash
cd docker_setup

# Build prerequisite image
./docker_build_push.sh Dockerfile.prerequisite.IMPACTncdENGL

# Build main image
./docker_build_push.sh Dockerfile.IMPACTncdENGL

# Build and push to Docker Hub (requires DOCKERHUB_USERNAME / DOCKERHUB_TOKEN)
./docker_build_push.sh Dockerfile.IMPACTncdENGL --push
```

To produce the `impactncdengl:local` image expected by `setup_user_docker_env.sh -Tag local`, override the default `--image-name`:

```bash
./docker_build_push.sh --image-name impactncdengl --image-tag local Dockerfile.IMPACTncdENGL
```

### Build script options

```
docker_build_push.sh <Dockerfile> [--image-name <name>] [--image-tag <tag>] [--push]
docker_build_push.ps1 <Dockerfile> [-ImageName <name>] [-ImageTag <tag>] [-Push]
```

| Option | Description |
|---|---|
| `<Dockerfile>` | Required. Dockerfile to build. |
| `--image-name` / `-ImageName` | Image name (default: derived from Dockerfile filename, lowercased). |
| `--image-tag` / `-ImageTag` | Image tag (default: `local`). |
| `--push` / `-Push` | Push to Docker Hub after a successful build. |

Push credentials come from the environment or a `.env` file alongside the script:
```
export DOCKERHUB_USERNAME=yourusername
export DOCKERHUB_TOKEN=youraccesstoken
```

### Disk space

The main image is large (~20 GB) due to bundled input data. Make sure Docker has room:
- **Linux:** check the `/var` partition or set `data-root` in `/etc/docker/daemon.json`.
- **Windows / macOS:** check Docker Desktop disk allocation in settings.

### Verbose builds (see package version substitutions)

```bash
# Bash
docker build -f Dockerfile.prerequisite.IMPACTncdENGL --progress=plain --no-cache -t my-image . 2>&1 \
  | grep -E "(Processing package|Successfully installed|Version.*not available|Installing available|PACKAGE VERSION)"
```

```powershell
# PowerShell
docker build -f Dockerfile.prerequisite.IMPACTncdENGL --progress=plain --no-cache -t my-image . 2>&1 `
  | Select-String -Pattern "(Processing package|Successfully installed|Version.*not available|Installing available|PACKAGE VERSION)"
```

---

## 📦 Package management

### System packages — `apt-packages.txt`

Pinned Ubuntu / Debian package versions:
```
automake=1:1.16.5-1.3ubuntu1
cmake=3.28.3-1build7
git=1:2.43.0-1ubuntu7.3
```

The intelligent installer (`install_packages.sh`) falls back to the nearest available version when a pinned version has been rotated out of the apt repository, and reports what it substituted.

### R packages — `r-packages.txt`

Installed from the Posit Package Manager snapshot whose date is the first line of `r-packages.txt`. The snapshot URL pins the versions; this file just lists the package names.

```
# CRAN snapshot date: 2026-02-18
data.table
ggplot2
fst
```

### Updating pinned package versions

When a build reports unavailable apt versions, refresh `apt-packages.txt`:

```bash
# Linux/macOS
./update-apt-packages.sh -i              # interactive
./update-apt-packages.sh -f build.log    # from a build log file

# Long-form aliases also work: --interactive, --file
```

```powershell
# Windows
.\update-apt-packages.ps1 -Interactive
.\update-apt-packages.ps1 -BuildLogFile "build.log"
```

To bump the R package snapshot, edit the first line of `r-packages.txt` and rebuild the prerequisite image.

---

## 🔄 Development workflow

1. **Make code changes** in the R package or project files.
2. **Rebuild the main image:**
   ```bash
   ./docker_build_push.sh Dockerfile.IMPACTncdENGL
   ```
3. **Test locally:**
   ```bash
   ./setup_user_docker_env.sh -Tag local
   ```
4. **Push to Docker Hub when ready:**
   ```bash
   ./docker_build_push.sh Dockerfile.IMPACTncdENGL --push
   ```

### Adding R or system dependencies

1. Edit `r-packages.txt` or `apt-packages.txt`.
2. Rebuild the prerequisite image:
   ```bash
   ./docker_build_push.sh Dockerfile.prerequisite.IMPACTncdENGL
   ```
3. Rebuild the main image:
   ```bash
   ./docker_build_push.sh Dockerfile.IMPACTncdENGL
   ```

---

## 🏷 Project details

| Property | Value | Source of truth |
|---|---|---|
| Base image | `rocker/r-ver:4.5.2` | `Dockerfile.prerequisite.IMPACTncdENGL` line 1 |
| R version | 4.5.2 | implied by base image |
| CRAN snapshot date | see first line of `r-packages.txt` (currently **2026-02-18**) | `r-packages.txt` |
| Package manager | [Posit Package Manager](https://packagemanager.posit.co/) | — |
| Docker Hub | [chriskypri/impactncdengl](https://hub.docker.com/r/chriskypri/impactncdengl) | — |
| Supported platforms | Windows 10/11 (incl. PowerShell-in-WSL), macOS, Linux (Ubuntu / Debian / Fedora / CentOS / Rocky) | — |

> The "source of truth" column is intentional: pinning version numbers in this README directly invites drift. If you change R or the snapshot date, update those files — this table reads from them.
