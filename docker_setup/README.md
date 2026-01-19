# IMPACTncd_England Docker Setup - Complete Guide

This repository contains the Docker configuration and cross-platform setup scripts for the **IMPACTncd England** project. The system provides robust, reproducible containerized environments using pre-built Docker images.

---

## 💾 Prerequisites

First, install Docker for your operating system:

- **Windows:** [Docker Desktop for Windows](https://docs.docker.com/desktop/install/windows-install/)
- **macOS:** [Docker Desktop for Mac](https://docs.docker.com/desktop/install/mac-install/)
- **Linux:** Follow instructions for your distribution:
  - [Ubuntu](https://docs.docker.com/engine/install/ubuntu/)
  - [Debian](https://docs.docker.com/engine/install/debian/)
  - [Fedora](https://docs.docker.com/engine/install/fedora/)
  - [CentOS](https://docs.docker.com/engine/install/centos/)

**Additional Platform Requirements:**

**Windows:**
- PowerShell 5.1+ (usually pre-installed)

**macOS:**
```bash
# Install coreutils which provides gsha256sum (only needed for certain operations)
brew install coreutils
```

**Linux:**
```bash
# Ubuntu/Debian - ensure Docker is installed
sudo apt-get update && sudo apt-get install docker.io

# Add user to docker group (optional, requires logout/login)
sudo usermod -aG docker $USER
```

---

## 🚀 Quick Start

### Container Selection

The scripts now use pre-built Docker images instead of building locally:

- **Default (main):** `chriskypri/impactncdengl:main` - pulls from Docker Hub
- **Local images:** `impactncdengl:local` - uses local Docker registry (when tag="local")
- **Remote images:** `chriskypri/impactncdengl:<tag>` - pulls from Docker Hub

### Basic Usage

Choose your platform and run the setup:

### Windows (PowerShell)
```powershell
# Default main image
.\setup_user_docker_env.ps1

# Use local image
.\setup_user_docker_env.ps1 -Tag "local"

# Use specific remote tag
.\setup_user_docker_env.ps1 -Tag "v1.2.3"

# Use volumes for better performance
.\setup_user_docker_env.ps1 -Tag "latest" -UseVolumes

# Custom YAML with remote image
.\setup_user_docker_env.ps1 -Tag "v1.2.3" -SimDesignYaml "C:\path\to\custom_sim_design.yaml"
```

### macOS/Linux (Bash)
```bash
# Default main image
./setup_user_docker_env.sh

# Use local image
./setup_user_docker_env.sh local

# Use specific remote tag
./setup_user_docker_env.sh v1.2.3

# Use volumes for better performance
./setup_user_docker_env.sh latest --UseVolumes

# Custom YAML with remote image
./setup_user_docker_env.sh v1.2.3 /path/to/custom_sim_design.yaml
```

### Script Parameters

**Bash Script (`setup_user_docker_env.sh`):**
- `[tag]`: Optional Docker image tag (default: "main") - first positional argument
- `[path_to_yaml]`: Optional path to simulation design YAML file
- `--UseVolumes`: Use Docker volumes for enhanced I/O performance

**PowerShell Script (`setup_user_docker_env.ps1`):**
- `-Tag <tag>`: Docker image tag (default: "main") - first parameter
- `-SimDesignYaml <path>`: Optional path to simulation design YAML file
- `-UseVolumes`: Switch for Docker volumes

---

## 📁 Directory Structure

When you run the setup, the following directories are mounted:

| Host Path | Container Mount | Description |
|-----------|-----------------|-------------|
| **Pre-built in image** | `/IMPACTncd_England` | Main project directory (already in Docker image) |
| `output_dir` from `sim_design.yaml` | `/output` | Simulation outputs |
| `synthpop_dir` from `sim_design.yaml` | `/synthpop` | Synthetic population data |

**Note:** The Docker images already contain the `/IMPACTncd_England` project folder, so no project directory mounting or copying is required.

## 🐳 Docker Images

### Image Selection Logic
- **Tag = "local"** → Uses `impactncdengl:local`
- **Tag = "main" (default)** → Uses `chriskypri/impactncdengl:main`
- **Tag = anything else** → Uses `chriskypri/impactncdengl:<tag>`

### Available Tags
Check Docker Hub for available tags: https://hub.docker.com/r/chriskypri/impactncdengl/tags

### Volume vs Bind Mount Modes

**Volume Mode (`--UseVolumes` / `-UseVolumes`):**
- Recommended for Windows and macOS
- Better I/O performance
- Creates temporary Docker volumes
- Syncs data back to host after container exit

**Bind Mount Mode (default):**
- Recommended for Linux
- Direct filesystem access
- Real-time file visibility
- Lower overhead on Linux

---

## ❓ Troubleshooting

### Common Issues

**"Cannot connect to the Docker daemon"**
- **Windows/macOS:** Make sure Docker Desktop is running (check system tray)
- **Linux:** Start Docker service: `sudo systemctl start docker`
- **Test:** Run `docker info` to verify Docker is accessible

**"Failed to pull Docker image"**
- **Check image exists:** Verify the tag exists on Docker Hub or locally
- **Network connectivity:** Ensure internet connection for remote images
- **Authentication:** For private repositories, run `docker login`
- **Manual test:** Try `docker pull impactncdengl:local` or `docker pull chriskypri/impactncdengl:latest`

**Windows: "Execution policy" error**  
- Run in PowerShell: `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`

**Linux: Permission denied**
- Add user to docker group: `sudo usermod -aG docker $USER` (then logout/login)
- Or use `sudo docker` commands

---

## 🧼 Cleanup

When you're done, you can clean up Docker resources:

```bash
# Remove the Docker image
docker rmi impactncd-england-r-prerequisite:latest

# Clean up unused containers/images
docker system prune
```

---

## 📬 Need Help?

For assistance or issues:
- Check the troubleshooting section above
- Open an issue in this repository
- Contact project maintainers

---

# 🔧 Technical Details & Developer Documentation

The following sections contain technical details, advanced usage, and developer information.

## 🌟 Intelligent Package Management System

---

## 🐳 Docker Configuration

### Container Specifications
- **Image Name:** `impactncd-england-r-prerequisite:latest`
- **Base Image:** [rocker/r-ver:4.5.1](https://hub.docker.com/r/rocker/r-ver)
- **R Version:** 4.5.1 with packages frozen as of July 20, 2025
- **Package Manager:** [RStudio Package Manager](https://packagemanager.posit.co/client/#/)

### File Structure
- **`Dockerfile.prerequisite.IMPACTncdENGL`**: Main Docker configuration with intelligent installer
- **`apt-packages.txt`**: System packages with pinned versions
- **`r-packages.txt`**: R packages list
- **`install_packages.sh`**: Intelligent package installer script
- **`entrypoint.sh`**: Container entrypoint for dynamic user creation

---

## � Usage Guide

### Basic Docker Build
```bash
# All platforms
docker build -f Dockerfile.prerequisite.IMPACTncdENGL -t impactncd-england .
```

### Verbose Build (See Package Updates)
```bash
# Windows PowerShell
docker build -f Dockerfile.prerequisite.IMPACTncdENGL --progress=plain --no-cache -t my-image . 2>&1 | Select-String -Pattern "(Processing package|Successfully installed|Version.*not available|Installing available|PACKAGE VERSION)"

# macOS/Linux Bash  
docker build -f Dockerfile.prerequisite.IMPACTncdENGL --progress=plain --no-cache -t my-image . 2>&1 | grep -E "(Processing package|Successfully installed|Version.*not available|Installing available|PACKAGE VERSION)"
```

### Update Package Versions

When you see package version updates in build output:

**Windows:**
```powershell
# Interactive update
.\update-apt-packages.ps1 -Interactive

# From build log file
.\update-apt-packages.ps1 -BuildLogFile "build.log"
```

**macOS/Linux:**
```bash
# Make executable (first time only)
chmod +x update-apt-packages.sh

# Interactive update
./update-apt-packages.sh --interactive

# From build log file  
./update-apt-packages.sh --file build.log

# Pipe directly from build
docker build ... 2>&1 | ./update-apt-packages.sh
```

## 🏗 Environment Setup Modes

### Bind Mount Mode (Default)
- Direct bind mounts between host and container
- Real-time interaction between host and container
- More flexible but potentially slower

### Volume Mode (Recommended for Large Simulations)
Use `--UseVolumes` (Bash) or `-UseVolumes` (PowerShell):
- Project directory copied to Docker-managed volume
- Better performance for large datasets
- Includes post-simulation sync back to host

## 🔄 Development Workflow

The same workflow applies across all platforms:

1. **Edit Code**: Modify R code or dependencies
2. **Build Image**: Run `docker build` command  
3. **Check Updates**: Look for package version update messages
4. **Update Versions**: Use appropriate update script for your platform
5. **Rebuild**: Build again with updated versions for reproducibility

### Multi-Platform Team Example
```bash
# Developer on macOS
./update-apt-packages.sh --interactive
git add apt-packages.txt
git commit -m "Update git package version"
git push

# Developer on Windows
git pull
.\update-apt-packages.ps1 -Interactive  # Shows no updates needed

# Developer on Linux  
git pull
./update-apt-packages.sh --interactive  # Shows no updates needed
```

## 🛠 Build and Push Scripts

Build and optionally push Docker images:

**Prerequisite Container:**
- Linux/macOS: `./build_push_prerequisite.sh [--push]`
- Windows: `build_push_prerequisite.ps1 [-Push]`

**IMPACTncd Container:**
- Linux/macOS: `./build_push_IMPACTncdENGL.sh [--push]`
- Windows: `build_push_IMPACTncdENGL.ps1 [-Push]`

## 🔍 Advanced Troubleshooting

### General Docker Issues
- **Check Docker Status:** Run `docker info`
- **Expected:** Detailed Docker installation info without errors
- **Common Error:** "Cannot connect to the Docker daemon"

**Solutions:**
- **Windows/macOS:** Ensure Docker Desktop is running
- **Linux:** Check service with `sudo systemctl status docker`
  - Start if needed: `sudo systemctl start docker`
  - Add user to docker group: `sudo usermod -aG docker $USER` (then logout/login)

### Platform-Specific Issues

**Windows:**
- Use PowerShell (not Command Prompt) for best compatibility
- May need to adjust execution policy: `Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass`
- Check WSL2 configuration for Docker Desktop

**macOS:**
- Install required tools: `brew install coreutils`
- Make scripts executable: `chmod +x update-apt-packages.sh`
- Check Docker Desktop permissions in System Preferences

**Linux:**
- Ensure Docker service is running: `sudo systemctl start docker`
- Make scripts executable: `chmod +x update-apt-packages.sh`
- May need `sudo docker` if user not in docker group

### Package Update Issues
- If specific versions fail, the intelligent installer automatically uses available versions
- Check build output for "PACKAGE VERSION UPDATES DETECTED" section
- Update `apt-packages.txt` with new versions for reproducibility

## 📦 Package Management Details

### System Packages (`apt-packages.txt`)
Contains Ubuntu packages with pinned versions:
```
automake=1:1.16.5-1.3ubuntu1
cmake=3.28.3-1build7
git=1:2.43.0-1ubuntu7.3
# ... more packages
```

### R Packages (`r-packages.txt`)  
R packages from CRAN snapshot (July 20, 2025):
```
data.table
ggplot2
dplyr
# ... more packages
```

### Version Update Workflow
1. **Automatic Detection**: Build detects unavailable versions
2. **Clear Reporting**: Shows which packages were updated
3. **Manual Update**: Update package files with new versions
4. **Automation**: Use provided scripts to streamline updates

## 🎯 Benefits

1. **Reproducible Builds**: Version pinning where possible, intelligent fallback when needed
2. **Cross-Platform**: Identical behavior on Windows, macOS, and Linux
3. **Resilient**: Builds don't fail due to outdated package versions
4. **Clear Feedback**: Visual indicators show exactly what happened
5. **Easy Maintenance**: Automated tools for version management
6. **Team Collaboration**: Consistent environment across different platforms

## 📝 Reproducibility Notes

While pinning versions with `<package>=<version>` is the standard apt method, package repositories don't guarantee indefinite availability of historical versions. Security updates and repository cleanup can affect older versions.

**Advanced Options for Maximum Reproducibility:**

1. **Repository Mirroring**: Create local mirrors of Ubuntu repositories at specific points in time (complex but complete control)

2. **Multi-Stage Builds with .deb Files**: Store specific .deb files alongside Dockerfile and install using `dpkg -i` (harder dependency management)

The intelligent package system balances reproducibility with practical build reliability by maintaining version pinning where possible while gracefully handling unavailable versions.

## 🏷 Project Details

- **Supported Branch:** master of IMPACTncd England model
- **R Version:** 4.5.1
- **Package Snapshot:** July 20, 2025
- **Base Image:** rocker/r-ver:4.5.1
- **Supported Platforms:** Windows 10/11, macOS, Linux (Ubuntu/Debian/CentOS/Fedora)

---

*This documentation combines intelligent package management, cross-platform compatibility, and comprehensive setup instructions for the IMPACTncd_England Docker environment.*
