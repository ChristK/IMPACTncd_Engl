# -----------------------------------------------------------------------------
# setup_user_docker_env.ps1
#
# PowerShell script for pulling and running a Docker container for the
# IMPACTncd England project. This script supports two operation modes:
#
# USAGE:
#   .\setup_user_docker_env.ps1 [[-Tag] <string>] [[-ScenariosDir] <string>] [[-SimDesignYaml] <string>] [-UseVolumes]
#
# EXAMPLES:
#   # Basic usage with default settings (uses main tag, bind mounts)
#   .\setup_user_docker_env.ps1
#   
#   # Use a specific Docker image tag
#   .\setup_user_docker_env.ps1 -Tag "v1.2.3"
#   .\setup_user_docker_env.ps1 -Tag "local"
#   
#   # Mount custom scenarios directory
#   .\setup_user_docker_env.ps1 -ScenariosDir "C:\path\to\my\scenarios"
#   .\setup_user_docker_env.ps1 -ScenariosDir "..\custom_scenarios"
#   
#   # Use a custom YAML configuration file
#   .\setup_user_docker_env.ps1 -SimDesignYaml "..\inputs\sim_design_test.yaml"
#   
#   # Use Docker volumes for better performance (recommended for macOS/Windows)
#   .\setup_user_docker_env.ps1 -UseVolumes
#   
#   # Combine multiple options
#   .\setup_user_docker_env.ps1 -Tag "v1.2.3" -ScenariosDir "..\my_scenarios" -UseVolumes
#   .\setup_user_docker_env.ps1 -Tag "local" -SimDesignYaml "..\inputs\sim_design_clbr.yaml" -UseVolumes
#   
#   # Using positional parameters
#   .\setup_user_docker_env.ps1 "local" "..\scenarios" "..\inputs\sim_design.yaml" -UseVolumes
#
# PARAMETERS:
#   -Tag <string>
#       Docker image tag to use. Options:
#       • "main" (default): pulls chriskypri/impactncdengl:main
#       • "local": uses locally built impactncdengl:local  
#       • Any other value: pulls chriskypri/impactncdengl:<tag>
#       
#   -ScenariosDir <string>
#       Path to scenarios directory to mount in container.
#       Will be available as /IMPACTncd_England/scenarios inside container.
#       Default: None (no scenarios mounted)
#       
#   -SimDesignYaml <string>
#       Path to the simulation design YAML file. Can be relative or absolute.
#       Default: "..\inputs\sim_design.yaml"
#       
#   -UseVolumes [<SwitchParameter>]
#       Use Docker-managed volumes instead of direct bind mounts.
#       Recommended for macOS and Windows for better I/O performance.
#       Default: False (uses bind mounts)
#
# Container Selection:
#   - If Tag is "main" (default): pulls and uses "chriskypri/impactncdengl:main"
#   - If Tag is "local": uses "impactncdengl:local" (built locally)
#   - If Tag is any other value: pulls and uses "chriskypri/impactncdengl:<tag>"
#
# Scenarios Directory:
#   - If ScenariosDir is provided, that directory will be mounted as
#     /IMPACTncd_England/scenarios inside the container, making the scenarios
#     files available at runtime.
#
# Operation Modes:
# 1. Using Docker-managed volumes (with -UseVolumes):
#      - Creates Docker volumes for output_dir and synthpop_dir (defined in YAML).
#      - Pre-populates volumes from local folders.
#      - Synchronizes volumes back to local folders after container exits.
#      - Removes volumes after synchronization.
#
# 2. Using direct bind mounts (default):
#      - Mounts local directories directly into the container.
#      - Changes are immediately visible on the host filesystem.
#
# Security:
#   - Containers run as non-root users to prevent permission issues.
#   - On Windows, Docker Desktop runs containers in a Linux VM, so UID/GID 1000:1000
#     is used for compatibility.
#
# Key Features:
# - Validates YAML and scenarios directory paths.
# - Extracts output_dir and synthpop_dir from YAML and ensures they exist.
# - Supports path normalization for Docker compatibility.
# - Provides detailed error messages for Docker connectivity and permission issues.
# - Automatically builds and uses an rsync-alpine image for volume synchronization.
#
# Notes:
# - If you encounter an execution policy error, run:
#     Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
# - For macOS and Windows, using Docker volumes (-UseVolumes) is recommended for better performance.
# - For Linux, ensure your user has Docker permissions (e.g., part of the "docker" group).
# -----------------------------------------------------------------------------

param (
    [string]$Tag = "main",
    [string]$ScenariosDir = "",
    [string]$SimDesignYaml = "..\inputs\sim_design.yaml",
    [switch]$UseVolumes
)

# Resolve script directory and switch to it
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $ScriptDir

# Resolve project root directory (one level above the current script directory)
$ProjectRoot = (Resolve-Path "$ScriptDir/..").Path -replace '\\', '/'

# If SimDesignYaml is a relative path, resolve it relative to the project root
if (-not [System.IO.Path]::IsPathRooted($SimDesignYaml)) {
    # Normalize path separators to forward slashes for cross-platform compatibility
    $SimDesignYamlNormalized = $SimDesignYaml -replace '\\', '/'
    $TempPath = "$ProjectRoot/$SimDesignYamlNormalized" -replace '/+', '/'
    # Resolve the path to handle .. components properly
    $SimDesignYaml = (Resolve-Path $TempPath -ErrorAction SilentlyContinue).Path
    if (-not $SimDesignYaml) {
        # If Resolve-Path fails, try manual construction (for the actual inputs directory)
        if ($SimDesignYamlNormalized -eq "../inputs/sim_design.yaml") {
            $SimDesignYaml = "$ProjectRoot/inputs/sim_design.yaml"
        } else {
            $SimDesignYaml = $TempPath
        }
    }
}

# Validate that the YAML file exists
if (-not (Test-Path $SimDesignYaml)) {
    Write-Host "Error: YAML file not found at '$SimDesignYaml'"
    Write-Host "Original path provided: '..\inputs\sim_design.yaml'"
    Write-Host "Project root: '$ProjectRoot'"
    Pop-Location
    Exit 1
}

Write-Host "Using configuration file: $SimDesignYaml"

# Validate scenarios directory if provided
if ($ScenariosDir -and -not (Test-Path $ScenariosDir -PathType Container)) {
    Write-Host "Error: Scenarios directory not found at '$ScenariosDir'"
    Pop-Location
    Exit 1
}

if ($ScenariosDir) {
    # Resolve to absolute path with forward slashes
    $ScenariosDir = (Resolve-Path $ScenariosDir).Path -replace '\\', '/'
    Write-Host "Using scenarios from: $ScenariosDir"
}

# Variable definitions
# Determine the Docker image name based on the tag
if ($Tag -eq "local") {
    $ImageName = "impactncdengl:local"
} else {
    $ImageName = "chriskypri/impactncdengl:$Tag"
}

Write-Host "Using Docker image: $ImageName"

# Get current user (for user-specific volume names)
$CurrentUser = $env:USERNAME

# Sanitize username for Docker volume names (replace spaces and special characters with underscores)
$SafeCurrentUser = $CurrentUser -replace '[^a-zA-Z0-9]', '_' -replace '__+', '_' -replace '^_|_$', ''
if ([string]::IsNullOrEmpty($SafeCurrentUser)) {
    $SafeCurrentUser = "dockeruser"
}

# Get user identity information for non-root Docker execution
# Note: On Windows, Docker Desktop runs containers in a Linux VM, so we use
# default UID/GID (1000:1000) which works well for most cases
$UserId = "1000"
$GroupId = "1000"
$UserName = $CurrentUser
# Use a safe group name - if it conflicts, the entrypoint will create a fallback
$GroupName = "dockergroup"

# Define user-specific Docker volume names using sanitized username (only for output and synthpop)
$VolumeOutput    = "impactncd_england_output_$SafeCurrentUser"
$VolumeSynthpop  = "impactncd_england_synthpop_$SafeCurrentUser"

# --- Docker Permission Check ---
# Check if the user can connect to the Docker daemon
Write-Host "Checking Docker daemon connectivity..."
docker info > $null 2>&1
if ($LASTEXITCODE -ne 0) {
    Write-Host "---------------------------------------------------------------------" -ForegroundColor Red
    Write-Host "Error: Cannot connect to the Docker daemon." -ForegroundColor Red
    Write-Host "This usually means Docker Desktop (Windows/macOS) or the Docker service (Linux) is not running or your user account lacks permission."
    Write-Host "Please ensure Docker is running and accessible before proceeding."
    Write-Host "" # Blank line for spacing
    Write-Host "How to check/fix:" -ForegroundColor Yellow
    Write-Host "  1. Run 'docker info' in your terminal. If it fails with a similar error, Docker is not accessible."
    Write-Host "  2. Windows/macOS: Make sure Docker Desktop is running (check the system tray or application list)."
    Write-Host "  3. Linux: Check service status with 'sudo systemctl status docker'. If inactive, start it with 'sudo systemctl start docker'."
    Write-Host "     You might also need to add your user to the 'docker' group ('sudo usermod -aG docker $env:USERNAME') and then log out and back in."
    Write-Host "  4. If running in WSL (Windows Subsystem for Linux), ensure Docker Desktop's WSL integration is enabled for your distribution."
    Write-Host "---------------------------------------------------------------------" -ForegroundColor Red
    Pop-Location # Restore original location before exiting
    Exit 1
} else {
    Write-Host "Docker daemon connection successful."
}
# --- End Docker Permission Check ---

# -----------------------------
# Pull the Docker image
# -----------------------------
Write-Host "Pulling Docker image: $ImageName"
docker pull $ImageName
if ($LASTEXITCODE -ne 0) {
    Write-Host "Error: Failed to pull Docker image: $ImageName" -ForegroundColor Red
    Write-Host "Please check:" -ForegroundColor Yellow
    Write-Host "  1. The image exists and is accessible"
    Write-Host "  2. You have the correct permissions"
    Write-Host "  3. Your internet connection is working"
    if ($Tag -ne "local") {
        Write-Host "  4. The tag '$Tag' exists in the chriskypri/impactncdengl repository"
    }
    Pop-Location
    Exit 1
}

# -----------------------------
# Extract paths from YAML
# -----------------------------
# Helper function to extract and construct potential paths from the YAML file
function Get-YamlPathValue {
    param (
        [string]$YamlPath,
        [string]$Key,
        [string]$BaseDir # Pass ProjectRoot here (already uses forward slashes)
    )
    $line = Select-String -Path $YamlPath -Pattern "^$Key\s*:" | Select-Object -First 1
    if ($line) {
        $value = ($line.Line -split ":\s*", 2)[1].Split("#")[0].Trim()
        $constructedPath = $null

        # Check if the path from YAML is absolute (Windows or Unix-like)
        if ([System.IO.Path]::IsPathRooted($value) -or $value.StartsWith('/')) {
            $constructedPath = $value
            Write-Host "Path '$value' for key '$Key' is absolute."
        } else {
            # Construct path relative to the specified BaseDir (ProjectRoot)
            # Ensure BaseDir and value use consistent slashes for joining
            $valueNormalized = $value -replace '\\', '/'
            $constructedPath = "$BaseDir/$valueNormalized" # Simple string concatenation with forward slashes
            # Clean up potential double slashes, except after protocol like C://
            $constructedPath = $constructedPath -replace '(?<!:)/{2,}', '/'
            Write-Host "Path '$value' for key '$Key' is relative. Constructed as '$constructedPath'."
        }

        # Normalize to forward slashes for consistency before returning
        $normalizedPath = $constructedPath -replace '\\', '/'
        return $normalizedPath
    }
    Write-Host "Warning: No matching line found for key: $Key in '$YamlPath'"
    return $null
}

# Call the function passing $ProjectRoot
$outputDir    = Get-YamlPathValue -YamlPath $SimDesignYaml -Key "output_dir" -BaseDir $ProjectRoot
$synthpopDir  = Get-YamlPathValue -YamlPath $SimDesignYaml -Key "synthpop_dir" -BaseDir $ProjectRoot

# --- Path Validation and Creation ---
# Helper function to check and create directory
function Test-AndCreateDirectory {
    param(
        [string]$Path,
        [string]$PathKey # For logging purposes (e.g., "output_dir")
    )
    if (-not $Path) {
        Write-Host "Error: Could not determine $PathKey path from YAML."
        return $false
    }

    # Use native path format for Test-Path and New-Item
    $NativePath = $Path -replace '/', '\\'

    if (-not (Test-Path $NativePath)) {
        Write-Host "Warning: $PathKey path not found: $NativePath. Creating directory..."
        try {
            New-Item -ItemType Directory -Path $NativePath -Force -ErrorAction Stop | Out-Null
            Write-Host "Successfully created $PathKey directory: $NativePath"
            return $true
        } catch {
            Write-Host "Error: Failed to create $PathKey directory: $NativePath - $($_.Exception.Message)"
            # Attempt to resolve the path to see if it exists now, maybe a race condition or delay
            if(Test-Path $NativePath) {
                 Write-Host "Info: Directory $NativePath seems to exist now despite previous error."
                 return $true
            }
            return $false
        }
    } elseif (-not (Get-Item $NativePath).PSIsContainer) {
        Write-Host "Error: The path specified for $PathKey exists but is a file, not a directory: $NativePath"
        return $false
    } else {
         # Directory exists
         return $true
    }
}

# Validate or create output directory
if (-not (Test-AndCreateDirectory -Path $outputDir -PathKey "output_dir")) {
    Pop-Location
    Exit 1
}

# Validate or create synthpop directory
if (-not (Test-AndCreateDirectory -Path $synthpopDir -PathKey "synthpop_dir")) {
    Pop-Location
    Exit 1
}
# --- End Path Validation and Creation ---


Write-Host "Mounting output_dir:    $outputDir"       # Keep using forward slashes for Docker mounts
Write-Host "Mounting synthpop_dir:  $synthpopDir"      # Keep using forward slashes for Docker mounts

# Helper function to convert Windows path to Docker Desktop/WSL format
function Convert-PathToDockerFormat {
    param([string]$Path)
    # Input example: P:/My_Models/IMPACTncd_England
    # Match drive letter (e.g., P) and the rest of the path
    if ($Path -match '^([A-Za-z]):/(.*)') {
        $driveLetter = $matches[1].ToLower()
        $restOfPath = $matches[2]
        # Construct the Docker path: /<drive_letter>/<rest_of_path>
        $dockerPath = "/$driveLetter/$restOfPath"
        # Remove trailing slash if present
        $dockerPath = $dockerPath -replace '/$', ''
        return $dockerPath
    } else {
        Write-Warning "Path '$Path' did not match expected Windows format (e.g., C:/path/to/dir)"
        return $Path # Return original path if format is unexpected
    }
}

# -----------------------------
# Run Docker container
# -----------------------------
if ($UseVolumes) {
    Write-Host "`nUsing Docker volumes for outputs and synthpop..."

    # Build rsync-alpine image if it doesn't already exist.
    $rsyncImage = "rsync-alpine"
    docker image inspect $rsyncImage > $null 2>&1
    if ($LASTEXITCODE -ne 0) {
        Write-Host "Building rsync-alpine image..."
        
        # Check if Dockerfile.rsync exists
        $DockerfileRsync = Join-Path $ScriptDir "Dockerfile.rsync"
        if (Test-Path $DockerfileRsync) {
            Write-Host "Using Dockerfile.rsync..."
            docker build -f "$DockerfileRsync" -t $rsyncImage $ScriptDir
        } else {
            Write-Host "Dockerfile.rsync not found, creating rsync image inline..."
            $InlineDockerfile = @"
FROM alpine:latest
RUN apk add --no-cache rsync
"@
            $InlineDockerfile | docker build -t $rsyncImage -
        }
    } else {
        Write-Host "Using existing rsync-alpine image."
    }

    # Ensure local output directories exist
    if (-not (Test-Path $outputDir)) { New-Item -ItemType Directory -Path $outputDir | Out-Null }
    if (-not (Test-Path $synthpopDir)) { New-Item -ItemType Directory -Path $synthpopDir | Out-Null }

    # Remove any existing volumes (ignore errors if not removable)
    Write-Host "Removing any existing volumes (if possible)..."
    docker volume rm $VolumeOutput -f 2>$null
    docker volume rm $VolumeSynthpop -f 2>$null

    # Create fresh Docker-managed volumes
    docker volume create $VolumeOutput | Out-Null
    docker volume create $VolumeSynthpop | Out-Null

    # Fix volume ownership and pre-populate volumes:
    # Docker volumes are created with root ownership by default. We need to fix
    # the ownership before we can populate them as the calling user.
    Write-Host "Setting correct ownership for Docker volumes..."
    docker run --rm -v "${VolumeOutput}:/volume" alpine sh -c "chown ${UserId}:${GroupId} /volume"
    docker run --rm -v "${VolumeSynthpop}:/volume" alpine sh -c "chown ${UserId}:${GroupId} /volume"

    # Pre-populate volumes:
    # The output and synthpop volumes are populated from the respective local folders.
    Write-Host "Populating output volume from local folder..."
    # Use permission-tolerant copy with fallback logic
    docker run --rm --user "${UserId}:${GroupId}" -v "${outputDir}:/source" -v "${VolumeOutput}:/volume" alpine sh -c "cp -r /source/. /volume/ 2>/dev/null || cp -a /source/. /volume/ 2>/dev/null || true"
    Write-Host "Populating synthpop volume from local folder..."
    # Use permission-tolerant copy with fallback logic
    docker run --rm --user "${UserId}:${GroupId}" -v "${synthpopDir}:/source" -v "${VolumeSynthpop}:/volume" alpine sh -c "cp -r /source/. /volume/ 2>/dev/null || cp -a /source/. /volume/ 2>/dev/null || true"

    # Run the main container with volumes mounted.
    Write-Host "Running the main container using Docker volumes..."
    # Construct arguments as an array for reliable passing
    $dockerArgs = @(
        "run", "-it", "--rm",
        # User identity environment variables
        "-e", "USER_ID=$UserId",
        "-e", "GROUP_ID=$GroupId", 
        "-e", "USER_NAME=$UserName",
        "-e", "GROUP_NAME=$GroupName",
        # Use -v syntax within the array elements (no project volume needed)
        "-v", "${VolumeOutput}:/output",
        "-v", "${VolumeSynthpop}:/synthpop"
    )
    
    # Add scenarios mount if provided
    if ($ScenariosDir) {
        $DockerScenariosDir = Convert-PathToDockerFormat -Path $ScenariosDir
        $dockerArgs += "--mount"
        $dockerArgs += "type=bind,source=$DockerScenariosDir,target=/IMPACTncd_England/scenarios"
    }
    
    # Add final arguments
    $dockerArgs += "--workdir"
    $dockerArgs += "/IMPACTncd_England"
    $dockerArgs += $ImageName
    $dockerArgs += "bash"
    
    # Execute docker with the arguments array
    & docker $dockerArgs

    # After the container exits:
    # Synchronize the output and synthpop volumes back to the local directories using rsync.
    Write-Host "Container exited. Syncing volumes back to local directories using rsync (checksum mode)..."
    # Use ${} to delimit variable name before the colon and add permission flags
    # Added --no-perms and --chmod=ugo=rwX to prevent permission issues on Windows
    docker run --rm --user "${UserId}:${GroupId}" -v "${VolumeOutput}:/volume" -v "${outputDir}:/backup" $rsyncImage rsync -avc --no-owner --no-group --no-times --no-perms --chmod=ugo=rwX /volume/ /backup/
    docker run --rm --user "${UserId}:${GroupId}" -v "${VolumeSynthpop}:/volume" -v "${synthpopDir}:/backup" $rsyncImage rsync -avc --no-owner --no-group --no-times --no-perms --chmod=ugo=rwX /volume/ /backup/

    # Clean up all the Docker volumes used for the simulation.
    Write-Host "Cleaning up Docker volumes..."
    docker volume rm $VolumeOutput | Out-Null
    docker volume rm $VolumeSynthpop | Out-Null

} else {
    Write-Host "`nUsing direct bind mounts for outputs and synthpop..."

    # Convert paths for Docker bind mount
    $DockerOutputDir = Convert-PathToDockerFormat -Path $outputDir
    $DockerSynthpopDir = Convert-PathToDockerFormat -Path $synthpopDir

    Write-Host "Docker Output Dir:   $DockerOutputDir"
    Write-Host "Docker Synthpop Dir: $DockerSynthpopDir"

    # Pass mount arguments correctly to docker run (no project mount needed)
    if ($ScenariosDir) {
        $DockerScenariosDir = Convert-PathToDockerFormat -Path $ScenariosDir
        docker run -it --rm `
            -e "USER_ID=$UserId" `
            -e "GROUP_ID=$GroupId" `
            -e "USER_NAME=$UserName" `
            -e "GROUP_NAME=$GroupName" `
            --mount "type=bind,source=$DockerOutputDir,target=/output" `
            --mount "type=bind,source=$DockerSynthpopDir,target=/synthpop" `
            --mount "type=bind,source=$DockerScenariosDir,target=/IMPACTncd_England/scenarios" `
            --workdir /IMPACTncd_England `
            $ImageName `
            bash
    } else {
        docker run -it --rm `
            -e "USER_ID=$UserId" `
            -e "GROUP_ID=$GroupId" `
            -e "USER_NAME=$UserName" `
            -e "GROUP_NAME=$GroupName" `
            --mount "type=bind,source=$DockerOutputDir,target=/output" `
            --mount "type=bind,source=$DockerSynthpopDir,target=/synthpop" `
            --workdir /IMPACTncd_England `
            $ImageName `
            bash
    }
}

# Restore the original directory
Pop-Location