# -----------------------------------------------------------------------------
# create_env.ps1
#
# PowerShell script for building and running a Docker container for the
# IMPACTncd Engl project on Windows.
#
# Features:
# - Accepts optional path to sim_design.yaml as argument
# - Extracts `output_dir` and `synthpop_dir` from YAML
# - Rebuilds Docker image only if build inputs have changed
# - Mounts project root, output_dir, and synthpop_dir into container
#
# Usage:
#   .\create_env.ps1 [path\to\sim_design.yaml]
#
# PowerShell script for building and running a Docker container for the
# IMPACTncd Engl project on Windows.
# Rebuilds the Docker image only if Dockerfile.prerequisite, apt-packages.txt,
# or r-packages.txt have changed.
# This script is designed to be run from the docker_setup directory, and binds 
# its parent directory to /IMPACTncd_Engl in the container.
# If you get an execution policy error, run: 
# Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
# -----------------------------------------------------------------------------
param (
    [string]$SimDesignYaml = "..\inputs\sim_design.yaml"
)

$ImageName = "impactncd-engl-r-prerequisite:latest"
$Dockerfile = "Dockerfile.prerequisite"
$HashFile = ".docker_build_hash"

# Resolve script directory
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $ScriptDir

# Validate YAML path
if (-not (Test-Path $SimDesignYaml)) {
    Write-Host "Error: YAML file not found at '$SimDesignYaml'"
    Exit 1
}

Write-Host "Using configuration file: $SimDesignYaml"

# Compute build hash from key input files with normalized newlines
$FilesToHash = @(
    (Get-Content "$Dockerfile" -Raw).Replace("`r`n", "`n")
    (Get-Content "apt-packages.txt" -Raw).Replace("`r`n", "`n")
    (Get-Content "r-packages.txt" -Raw).Replace("`r`n", "`n")
)
$HashAlgorithm = [System.Security.Cryptography.SHA256]::Create()
$Bytes = [System.Text.Encoding]::UTF8.GetBytes($FilesToHash -join "`n")
$HashBytes = $HashAlgorithm.ComputeHash($Bytes)
$BuildHash = [BitConverter]::ToString($HashBytes) -replace "-", ""

if (Test-Path $HashFile) {
    $LastHash = Get-Content $HashFile -Raw
    Write-Host "Previous build hash: $LastHash"
    Write-Host "Current build hash:  $BuildHash"
    if ($LastHash -ne $BuildHash) {
        Write-Host "Build hash mismatch — rebuild required."
    } else {
        Write-Host "Build hash match — no rebuild needed."
    }
}

# Determine if rebuild is needed
$NeedsBuild = $false
if (-not (docker image inspect $ImageName > $null 2>&1)) {
    Write-Host "Docker image does not exist. Need to build."
    $NeedsBuild = $true
} elseif (-not (Test-Path $HashFile)) {
    Write-Host "No previous build hash found. Need to build."
    $NeedsBuild = $true
} else {
    $LastHash = Get-Content $HashFile -Raw
    if ($LastHash -ne $BuildHash) {
        Write-Host "Detected changes in build inputs. Rebuilding Docker image..."
        $NeedsBuild = $true
    } else {
        Write-Host "No changes detected. Skipping Docker build."
    }
}

# Build Docker image if needed
if ($NeedsBuild) {
    Write-Host "Building Docker image using --no-cache..."
    docker build --no-cache -f $Dockerfile -t $ImageName .
    $BuildHash | Set-Content $HashFile
} else {
    Write-Host "Docker image is up to date. Skipping build."
}

# Extract output_dir and synthpop_dir from YAML
function Get-YamlPathValue {
    param (
        [string]$YamlPath,
        [string]$Key
    )

    $line = Select-String "^$Key\s*:" $YamlPath | Select-Object -First 1
    if ($line) {
        $value = ($line.Line -split ":\s*", 2)[1].Split("#")[0].Trim()

        # Resolve absolute path from YAML value
        $resolved = Resolve-Path $value -ErrorAction Stop

        # Convert to Docker-compatible mount path
        $dockerPath = $resolved.Path -replace '^([A-Za-z]):', {
            "/run/desktop/mnt/host/$($args[0].Groups[1].Value.ToLower())"
        } -replace '\\', '/'

        return $dockerPath
    }

    return $null
}

$outputDir = Get-YamlPathValue -YamlPath $SimDesignYaml -Key "output_dir"
$synthpopDir = Get-YamlPathValue -YamlPath $SimDesignYaml -Key "synthpop_dir"

if (-not $outputDir -or -not (Test-Path $outputDir)) {
    Write-Host "Error: output_dir path not found or invalid: $outputDir"
    Exit 1
}

if (-not $synthpopDir -or -not (Test-Path $synthpopDir)) {
    Write-Host "Error: synthpop_dir path not found or invalid: $synthpopDir"
    Exit 1
}

Write-Host "Mounting output_dir:    $outputDir"
Write-Host "Mounting synthpop_dir: $synthpopDir"

# Resolve project root directory (one level up)
$ProjectRoot = Resolve-Path "$ScriptDir\.."
$ProjectRoot = $ProjectRoot.Path -replace '\\', '/'

# Format other mount paths for Docker
$outputDir = $outputDir -replace '\\', '/'
$synthpopDir = $synthpopDir -replace '\\', '/'

# Run Docker container
docker run -it `
  --mount type=bind,source="$ProjectRoot",target=/IMPACTncd_Engl `
  --mount type=bind,source="$outputDir",target=/IMPACTncd_Engl/outputs `
  --mount type=bind,source="$synthpopDir",target=/IMPACTncd_Engl/synthpop `
  $ImageName `
  bash

Pop-Location