# IMPACTncd_England Docker Setup - Complete Guide

This directory contains Docker configuration and cross-platform setup scripts for the **IMPACTncd England** microsimulation project. The system provides reproducible containerized environments using pre-built Docker images from Docker Hub.

---

## 📋 Table of Contents

- [Prerequisites](#-prerequisites)
- [Quick Start](#-quick-start)
- [Directory Structure](#-directory-structure)
- [Docker Images](#-docker-images)
- [Troubleshooting](#-troubleshooting)
- [Developer Documentation](#-developer-documentation)

---

## 💾 Prerequisites

### Install Docker

- **Windows:** [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/)
- **macOS:** [Docker Desktop for Mac](https://docs.docker.com/desktop/install/mac-install/)
- **Linux:** Follow instructions for your distribution:
  - [Ubuntu](https://docs.docker.com/engine/install/ubuntu/)
  - [Debian](https://docs.docker.com/engine/install/debian/)
  - [Fedora](https://docs.docker.com/engine/install/fedora/)
  - [CentOS/Rocky Linux](https://docs.docker.com/engine/install/centos/)

### Platform-Specific Requirements

**Windows:**
- PowerShell 5.1+ (pre-installed on Windows 10/11)
- WSL2 backend recommended for Docker Desktop

**macOS:**
```bash
# Install coreutils for gsha256sum (optional, for certain operations)
brew install coreutils
```

**Linux:**
```bash
# Add user to docker group (avoids needing sudo for docker commands)
sudo usermod -aG docker $USER
# Log out and back in for the change to take effect
```

---

## 🚀 Quick Start

### For End Users (Running Simulations)

The easiest way to get started is using the pre-built Docker images:

**Linux/macOS:**
```bash
cd docker_setup

# Run with default settings (pulls chriskypri/impactncdengl:main from Docker Hub)
./setup_user_docker_env.sh

# Use a specific version tag
./setup_user_docker_env.sh -Tag v1.0.0

# Use a locally built image
./setup_user_docker_env.sh -Tag local

# Use Docker volumes for better I/O performance (recommended for macOS/Windows)
./setup_user_docker_env.sh --UseVolumes

# Specify custom simulation design YAML
./setup_user_docker_env.sh -SimDesignYaml /path/to/my_sim_design.yaml

# Mount a custom scenarios directory
./setup_user_docker_env.sh -ScenariosDir /path/to/my_scenarios
```

**Windows (PowerShell):**
```powershell
cd docker_setup

# Run with default settings
.\setup_user_docker_env.ps1

# Use a specific version tag
.\setup_user_docker_env.ps1 -Tag v1.0.0

# Use Docker volumes for better performance
.\setup_user_docker_env.ps1 -UseVolumes

# Specify custom simulation design YAML
.\setup_user_docker_env.ps1 -SimDesignYaml "C:\path\to\my_sim_design.yaml"
```

### What Happens When You Run the Script

1. **Pulls the Docker image** from Docker Hub (if not already cached)
2. **Reads your `sim_design.yaml`** to find output and synthpop directories
3. **Creates those directories** on your host if they don't exist
4. **Starts an interactive container** with your directories mounted
5. **Runs as your user** (not root) to avoid permission issues

Once inside the container, you can run simulations:
```r
source("global.R")
IMPACTncd <- Simulation$new("./inputs/sim_design.yaml")
IMPACTncd$run(1:10, multicore = TRUE, "sc0")$export_summaries(multicore = TRUE)
```

---

## 📁 Directory Structure

### Files in This Directory

| File | Description |
|------|-------------|
| `setup_user_docker_env.sh` | **Main script** for Linux/macOS users to run simulations |
| `setup_user_docker_env.ps1` | **Main script** for Windows users to run simulations |
| `docker_build_push.sh` | Build and push Docker images (developers) |
| `docker_build_push.ps1` | Build and push Docker images (Windows developers) |
| `Dockerfile.prerequisite.IMPACTncdENGL` | Base image with R and dependencies |
| `Dockerfile.IMPACTncdENGL` | Main image with project code and data |
| `apt-packages.txt` | System packages with pinned versions |
| `r-packages.txt` | R packages to install |
| `install_packages.sh` | Intelligent package installer |
| `entrypoint.sh` | Container entrypoint script |
| `update-apt-packages.sh` | Update package versions from build logs |

### Container Directory Mounts

When you run the setup script, these directories are mounted:

| Host Path | Container Path | Description |
|-----------|----------------|-------------|
| `output_dir` from YAML | `/output` | Simulation outputs (lifecourse, summaries, etc.) |
| `synthpop_dir` from YAML | `/synthpop` | Synthetic population cache |
| `scenarios/` (if specified) | `/IMPACTncd_England/scenarios` | Custom scenario scripts |

> **Note:** The project code (`/IMPACTncd_England`) and input data are already baked into the Docker image. Only output directories need to be mounted.

---

## 🐳 Docker Images

### Image Selection Logic

| Tag | Image Used | Source |
|-----|------------|--------|
| `main` (default) | `chriskypri/impactncdengl:main` | Docker Hub |
| `local` | `impactncdengl:local` | Local Docker registry |
| `v1.0.0`, etc. | `chriskypri/impactncdengl:<tag>` | Docker Hub |

### Available Tags

Check Docker Hub for available tags: https://hub.docker.com/r/chriskypri/impactncdengl/tags

### Bind Mount vs Volume Mode

| Mode | Flag | Best For | Description |
|------|------|----------|-------------|
| **Bind Mount** (default) | none | Linux | Direct filesystem access, real-time file visibility |
| **Volume Mode** | `--UseVolumes` | macOS, Windows | Better I/O performance, syncs data after container exit |

---

## ❓ Troubleshooting

### "Cannot connect to the Docker daemon"

**Windows/macOS:**
- Ensure Docker Desktop is running (check system tray icon)

**Linux:**
```bash
# Check Docker status
sudo systemctl status docker

# Start Docker if not running
sudo systemctl start docker

# Verify Docker works
docker info
```

### "Permission denied" (Linux)

```bash
# Add your user to the docker group
sudo usermod -aG docker $USER

# Log out and back in, then verify
groups  # should show 'docker'
```

### "Failed to pull Docker image"

- **Check internet connection** for remote images
- **Verify the tag exists:** `docker pull chriskypri/impactncdengl:main`
- **For private repos:** Run `docker login` first

### Windows: "Execution policy" error

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

### Container runs out of memory

Edit your `sim_design.yaml` and reduce:
- `clusternumber` (fewer parallel cores = less RAM needed)
- `n` and `num_chunks` (smaller synthetic population)

---

## 🧹 Cleanup

```bash
# Remove a specific image
docker rmi chriskypri/impactncdengl:main

# Remove all unused images and containers
docker system prune

# Remove everything including volumes (⚠️ destructive)
docker system prune -a --volumes
```

---

## 📬 Need Help?

- Check the [main project README](../README.md)
- Review the [vignettes](../Rpackage/IMPACTncd_England_model_pkg/vignettes/) for simulation guidance
- Open an issue on GitHub
- Contact project maintainers

---

# 🔧 Developer Documentation

The following sections are for developers who need to build Docker images or maintain the package system.

---

## 🏗️ Building Docker Images

There are two Docker images in a layered architecture:

1. **Prerequisite image** (`prerequisite.impactncdengl:local`): Base R environment with all dependencies
2. **Main image** (`impactncdengl:local`): Prerequisite + project code + input data (~15GB)

### Build Commands

```bash
cd docker_setup

# Build prerequisite image first (base R environment with packages)
./docker_build_push.sh Dockerfile.prerequisite.IMPACTncdENGL

# Build main image (includes project code and 14GB of input data)
./docker_build_push.sh Dockerfile.IMPACTncdENGL

# Build and push to Docker Hub
./docker_build_push.sh Dockerfile.IMPACTncdENGL --push
```

### Build Script Options

```bash
./docker_build_push.sh <Dockerfile> [--image-name <name>] [--image-tag <tag>] [--push]
```

| Option | Description |
|--------|-------------|
| `<Dockerfile>` | Required. The Dockerfile to build |
| `--image-name` | Custom image name (default: derived from Dockerfile) |
| `--image-tag` | Custom tag (default: `local`) |
| `--push` | Push to Docker Hub after building |

### Docker Storage Requirements

The main image is large (~20GB) due to input data. Ensure you have sufficient disk space:
- **Linux:** Check `/var` partition or configure Docker's `data-root` in `/etc/docker/daemon.json`
- **Windows/macOS:** Check Docker Desktop disk allocation in settings

---

## 📦 Package Management

### Container Specifications

- **Base Image:** [rocker/r-ver:4.5.1](https://hub.docker.com/r/rocker/r-ver)
- **R Version:** 4.5.1
- **CRAN Snapshot:** July 20, 2025
- **Package Manager:** [Posit Package Manager](https://packagemanager.posit.co/)

### System Packages (`apt-packages.txt`)

Ubuntu/Debian packages with pinned versions:
```
automake=1:1.16.5-1.3ubuntu1
cmake=3.28.3-1build7
git=1:2.43.0-1ubuntu7.3
```

### R Packages (`r-packages.txt`)

R packages installed from the frozen CRAN snapshot:
```
# CRAN snapshot date: 2025-07-20
data.table
ggplot2
fst
```

### Updating Package Versions

When Docker builds detect version mismatches, the intelligent installer (`install_packages.sh`) will:
1. Try to install the pinned version
2. Fall back to the latest available version if pinned version is unavailable
3. Report which packages were updated

To update `apt-packages.txt` with new versions:

```bash
# Interactive mode
./update-apt-packages.sh -i

# From a build log file
./update-apt-packages.sh -f build.log
```

---

## 🔄 Development Workflow

1. **Make code changes** in the R package or project files
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

### Updating R or System Dependencies

If you need to add new packages:

1. Edit `r-packages.txt` or `apt-packages.txt`
2. Rebuild the prerequisite image:
   ```bash
   ./docker_build_push.sh Dockerfile.prerequisite.IMPACTncdENGL
   ```
3. Rebuild the main image:
   ```bash
   ./docker_build_push.sh Dockerfile.IMPACTncdENGL
   ```

---

## 🏷️ Project Details

| Property | Value |
|----------|-------|
| R Version | 4.5.1 |
| CRAN Snapshot | January 20, 2026 |
| Base Image | rocker/r-ver:4.5.2 |
| Docker Hub | [chriskypri/impactncdengl](https://hub.docker.com/r/chriskypri/impactncdengl) |
| Platforms | Windows 10/11, macOS, Linux |

---

*Last updated: January 2026*
