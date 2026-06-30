# -----------------------------------------------------------------------------
# setup_dev_docker_env.ps1
#
# PowerShell script for building and running a Docker container for the
# IMPACTncd England project. This script supports two operation modes:
#
# USAGE:
#   .\setup_dev_docker_env.ps1 [[-SimDesignYaml] <string>] [-UseVolumes]
#
# EXAMPLES:
#   # Basic usage with default settings (uses bind mounts)
#   .\setup_dev_docker_env.ps1
#   
#   # Use a custom YAML configuration file
#   .\setup_dev_docker_env.ps1 -SimDesignYaml "C:\path\to\custom_config.yaml"
#   .\setup_dev_docker_env.ps1 "..\inputs\sim_design_test.yaml"
#   
#   # Use Docker volumes for better performance (recommended for macOS/Windows)
#   .\setup_dev_docker_env.ps1 -UseVolumes
#   
#   # Combine custom YAML with volumes
#   .\setup_dev_docker_env.ps1 -SimDesignYaml "..\inputs\sim_design_clbr.yaml" -UseVolumes
#   
#   # Using positional parameter (YAML path as first argument)
#   .\setup_dev_docker_env.ps1 "..\inputs\sim_design_test.yaml" -UseVolumes
#
# PARAMETERS:
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
#   - Builds a local development image: prerequisite.impactncdengl:local (from Dockerfile.prerequisite.IMPACTncdENGL).
#   - Automatically detects changes in build inputs (Dockerfile, apt-packages.txt, etc.) and rebuilds the image if needed.
#
# Operation Modes:
# 1. Using Docker-managed volumes (with -UseVolumes):
#      - Copies the project directory into a Docker volume for faster I/O.
#      - Creates separate volumes for output_dir and synthpop_dir (defined in YAML).
#      - Synchronizes volumes back to local folders after container exits.
#      - Removes volumes after synchronization.
#
# 2. Using direct bind mounts (default):
#      - Mounts local directories directly into the container.
#      - Changes are immediately visible on the host filesystem.
#
# Security:
#   - Containers run as the calling user (non-root) to prevent permission issues.
#   - Automatically detects the current user's UID and GID and passes them to Docker.
#
# Notes:
# - Compatible with Windows PowerShell and Linux/macOS (requires Docker).
# - For macOS and Windows, using Docker volumes (-UseVolumes) is recommended for better performance.
# - For Linux, ensure your user has Docker permissions (e.g., part of the "docker" group).
# - If you encounter permission issues, ensure the output_dir and synthpop_dir exist and are writable.
# - If you get an execution policy error, run:
#     Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
# -----------------------------------------------------------------------------

param (
    [string]$SimDesignYaml = "..\inputs\sim_design.yaml",
    [switch]$UseVolumes
)

# Resolve script directory and switch to it
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $ScriptDir

# Resolve project root directory (one level above the current script directory)
$ProjectRoot = (Resolve-Path "$ScriptDir/..").Path -replace '\\', '/'

# Resolve the YAML file path relative to the project root if it's a relative path
if (-not [System.IO.Path]::IsPathRooted($SimDesignYaml)) {
    # Normalize path separators to forward slashes for cross-platform compatibility
    $SimDesignYamlNormalized = $SimDesignYaml -replace '\\', '/'
    $TempPath = "$ProjectRoot/$SimDesignYamlNormalized" -replace '/+', '/'
    # Resolve the path to handle .. components properly
    $SimDesignYamlResolved = (Resolve-Path $TempPath -ErrorAction SilentlyContinue).Path
    if (-not $SimDesignYamlResolved) {
        # If Resolve-Path fails, try manual construction (for the actual inputs directory)
        if ($SimDesignYamlNormalized -eq "../inputs/sim_design.yaml") {
            $SimDesignYamlResolved = "$ProjectRoot/inputs/sim_design.yaml"
        } else {
            $SimDesignYamlResolved = $TempPath
        }
    }
} else {
    $SimDesignYamlResolved = $SimDesignYaml
}

# Validate that the YAML file exists
if (-not (Test-Path $SimDesignYamlResolved)) {
    Write-Host "Error: YAML file not found at '$SimDesignYamlResolved'"
    Write-Host "Original path provided: '$SimDesignYaml'"
    Write-Host "Project root: '$ProjectRoot'"
    Pop-Location
    Exit 1
}

Write-Host "Using configuration file: $SimDesignYamlResolved"

# Variable definitions
$ImageName   = "prerequisite.impactncdengl:local"
$Dockerfile  = "Dockerfile.prerequisite.IMPACTncdENGL"
$HashFile    = ".docker_build_hash"

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

# Define user-specific Docker volume names using sanitized username
$VolumeProject   = "impactncd_england_project_$SafeCurrentUser"
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
# Build hash and rebuild logic
# -----------------------------
# Compute build hash from Dockerfile, package lists, and entrypoint script
# Use the same method as the bash script for consistency
$TempFile = [System.IO.Path]::GetTempFileName()
try {
    # Concatenate files in the same way as the bash script
    $FilesToConcat = @($Dockerfile, "apt-packages.txt", "r-packages.txt", "entrypoint.sh")
    $CombinedContent = ""
    foreach ($File in $FilesToConcat) {
        $CombinedContent += Get-Content -Raw -Path $File -Encoding UTF8
    }
    
    # Write to temp file and compute hash using external command for consistency
    $CombinedContent | Set-Content -Path $TempFile -Encoding UTF8 -NoNewline
    $BuildHash = (Get-FileHash -Path $TempFile -Algorithm SHA256).Hash.ToLower()
} finally {
    if (Test-Path $TempFile) { Remove-Item $TempFile }
}

$NeedsBuild = $false
Write-Host "Checking for Docker image: '$ImageName'"
docker image inspect $ImageName > $null 2>&1
$InspectExitCode = $LASTEXITCODE
Write-Host "docker image inspect exit code: $InspectExitCode"

if ($InspectExitCode -ne 0) {
    Write-Host "Docker image does not exist or inspect failed. Need to build."
    $NeedsBuild = $true
} elseif (-not (Test-Path $HashFile)) {
    Write-Host "No previous build hash found. Need to build."
    $NeedsBuild = $true
} else {
    $LastHash = (Get-Content -Raw -Encoding UTF8 $HashFile).Trim()
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
    $BuildHash | Set-Content -NoNewline -Encoding UTF8 $HashFile
} else {
    Write-Host "Docker image is up to date. Skipping build."
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

        # Check if the path from YAML is absolute (Windows or Unix-like).
        # On Linux PowerShell, [Path]::IsPathRooted("C:/foo") returns false, so we
        # also test the Windows-drive pattern explicitly to handle YAMLs with
        # Windows absolute paths read by PowerShell running inside WSL.
        if ([System.IO.Path]::IsPathRooted($value) -or $value.StartsWith('/') -or $value -match '^[A-Za-z]:[/\\]') {
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
$outputDir    = Get-YamlPathValue -YamlPath $SimDesignYamlResolved -Key "output_dir" -BaseDir $ProjectRoot
$synthpopDir  = Get-YamlPathValue -YamlPath $SimDesignYamlResolved -Key "synthpop_dir" -BaseDir $ProjectRoot

# --- Path Validation and Creation ---
# Helper function to check and create directory
function Ensure-DirectoryExists {
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
if (-not (Ensure-DirectoryExists -Path $outputDir -PathKey "output_dir")) {
    Pop-Location
    Exit 1
}

# Validate or create synthpop directory
if (-not (Ensure-DirectoryExists -Path $synthpopDir -PathKey "synthpop_dir")) {
    Pop-Location
    Exit 1
}
# --- End Path Validation and Creation ---


Write-Host "Mounting output_dir:    $outputDir"       # Keep using forward slashes for Docker mounts
Write-Host "Mounting synthpop_dir:  $synthpopDir"      # Keep using forward slashes for Docker mounts

# -----------------------------
# Run Docker container
# -----------------------------
if ($UseVolumes) {
    Write-Host "`nUsing Docker volumes for project, outputs, and synthpop..."

    # Build rsync-alpine image if it doesn't already exist.
    $rsyncImage = "rsync-alpine"
    docker image inspect $rsyncImage > $null 2>&1
    if ($LASTEXITCODE -ne 0) {
        Write-Host "Building rsync-alpine image..."
        docker build -f "Dockerfile.rsync" -t $rsyncImage .
    } else {
        Write-Host "Using existing rsync-alpine image."
    }

    # Ensure local output directories exist
    if (-not (Test-Path $outputDir)) { New-Item -ItemType Directory -Path $outputDir | Out-Null }
    if (-not (Test-Path $synthpopDir)) { New-Item -ItemType Directory -Path $synthpopDir | Out-Null }

    # Remove any existing volumes (ignore errors if not removable)
    Write-Host "Removing any existing volumes (if possible)..."
    docker volume rm $VolumeProject -f 2>$null
    docker volume rm $VolumeOutput -f 2>$null
    docker volume rm $VolumeSynthpop -f 2>$null

    # Create fresh Docker-managed volumes
    docker volume create $VolumeProject | Out-Null
    docker volume create $VolumeOutput | Out-Null
    docker volume create $VolumeSynthpop | Out-Null

    # Fix volume ownership and pre-populate volumes:
    # Docker volumes are created with root ownership by default. We need to fix
    # the ownership before we can populate them as the calling user.
    Write-Host "Setting correct ownership for Docker volumes..."
    docker run --rm -v "${VolumeProject}:/destination" alpine sh -c "chown ${UserId}:${GroupId} /destination"
    docker run --rm -v "${VolumeOutput}:/volume" alpine sh -c "chown ${UserId}:${GroupId} /volume"
    docker run --rm -v "${VolumeSynthpop}:/volume" alpine sh -c "chown ${UserId}:${GroupId} /volume"

    # Pre-populate volumes:
    # 1. The project volume is populated with the entire project directory using tar, excluding dot files/folders.
    # 2. The output and synthpop volumes are populated from the respective local folders.
    Write-Host "Populating project volume from host project directory (excluding dot files/folders)..."
    # Use tar to copy project files, excluding folders starting with .
    docker run --rm --user "${UserId}:${GroupId}" -v "${ProjectRoot}:/source" -v "${VolumeProject}:/destination" alpine sh -c "cd /source && tar cf - --exclude='./.*' . | (cd /destination && tar xf -)"

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
        # Use -v syntax within the array elements
        "-v", "${VolumeProject}:/IMPACTncd_England",
        "-v", "${VolumeOutput}:/outputs",
        "-v", "${VolumeSynthpop}:/synthpop",
        "--workdir", "/IMPACTncd_England",
        $ImageName,
        "bash"
    )
    # Execute docker with the arguments array
    & docker $dockerArgs

    # After the container exits:
    # Synchronize the output, synthpop, and simulation volumes back to the local directories using rsync.
    Write-Host "Container exited. Syncing volumes back to local directories using rsync (checksum mode)..."
    # Use ${} to delimit variable name before the colon and add permission flags
    # Added --no-perms and --chmod=ugo=rwX to prevent permission issues on Windows
    docker run --rm --user "${UserId}:${GroupId}" -v "${VolumeOutput}:/volume" -v "${outputDir}:/backup" $rsyncImage rsync -avc --no-owner --no-group --no-times --no-perms --chmod=ugo=rwX /volume/ /backup/
    docker run --rm --user "${UserId}:${GroupId}" -v "${VolumeSynthpop}:/volume" -v "${synthpopDir}:/backup" $rsyncImage rsync -avc --no-owner --no-group --no-times --no-perms --chmod=ugo=rwX /volume/ /backup/
    # Sync simulation folder back to the project directory
    $SimulationDir = "$ProjectRoot/simulation" -replace '\\', '/'
    Write-Host "Syncing simulation folder back to: $SimulationDir"
    docker run --rm --user "${UserId}:${GroupId}" -v "${VolumeProject}:/project" -v "${SimulationDir}:/backup" $rsyncImage rsync -avc --no-owner --no-group --no-times --no-perms --chmod=ugo=rwX /project/simulation/ /backup/

    # Clean up all the Docker volumes used for the simulation.
    Write-Host "Cleaning up Docker volumes..."
    docker volume rm $VolumeProject | Out-Null
    docker volume rm $VolumeOutput | Out-Null
    docker volume rm $VolumeSynthpop | Out-Null

} else {
    # Helper function to translate a path into the form the local docker daemon
    # expects for bind-mount sources. The docker daemon on Windows is almost
    # always Linux-side these days (Docker Desktop's WSL2 backend, or docker
    # installed directly inside a WSL distro), so bind mounts must use a POSIX
    # path that's valid in that distro's filesystem — and the mount root for
    # Windows drives is configurable via /etc/wsl.conf [automount] root
    # (commonly /mnt, but /mnt/host/ on some setups).
    #
    # Strategy: ask `wslpath` (canonical, honours per-distro config) whether
    # we're running inside WSL or on native Windows with WSL installed. Fall
    # back to /etc/wsl.conf parsing or the legacy /c/foo Docker-Desktop short
    # form only if wslpath is unreachable.
    function Convert-PathToDockerFormat {
        param([string]$Path)

        # Normalize backslashes for parsing (preserve drive-letter casing)
        $normalized = $Path -replace '\\', '/'

        # Already a POSIX absolute path (e.g. /home/user, /mnt/c/foo) — pass through.
        if ($normalized.StartsWith('/')) {
            return ($normalized -replace '/$', '')
        }

        # Detect Windows-drive-rooted path (e.g. C:/foo, P:\bar)
        if ($normalized -match '^([A-Za-z]):/(.*)$') {
            $driveLetter = $matches[1].ToLower()
            $restOfPath  = $matches[2] -replace '/$', ''
            $winStyle    = "${driveLetter}:\" + ($restOfPath -replace '/', '\')

            $inWsl = $false
            if ($IsLinux) {
                if ($env:WSL_DISTRO_NAME -or (Test-Path '/proc/sys/fs/binfmt_misc/WSLInterop')) {
                    $inWsl = $true
                }
            }

            # Attempt wslpath resolution. Inside WSL this is the binary itself;
            # on native Windows we delegate to `wsl.exe wslpath -u` which runs
            # in the user's default WSL distro. That's the same distro Docker
            # Desktop's WSL2 backend integrates with, so the mount root it
            # reports matches what the daemon will see.
            $wslResolved = $null
            try {
                if ($inWsl) {
                    $out = (& wslpath -u $winStyle 2>$null)
                    if ($LASTEXITCODE -eq 0 -and $out) { $wslResolved = ($out -join "`n").Trim() }
                } elseif (Get-Command 'wsl.exe' -ErrorAction SilentlyContinue) {
                    $out = (& wsl.exe wslpath -u $winStyle 2>$null)
                    if ($LASTEXITCODE -eq 0 -and $out) {
                        # wsl.exe can emit UTF-16 BOMs / CRs on older Windows
                        $wslResolved = (($out -join "`n") -replace "`r", '').Trim()
                    }
                }
            } catch { }

            if ($wslResolved) { return ($wslResolved -replace '/$', '') }

            # Fallback for WSL when wslpath is unavailable: parse [automount] root.
            if ($inWsl) {
                $mountRoot = '/mnt'
                if (Test-Path '/etc/wsl.conf') {
                    $conf = Get-Content '/etc/wsl.conf' -Raw -ErrorAction SilentlyContinue
                    if ($conf -match '(?ms)^\s*\[automount\][^\[]*?^\s*root\s*=\s*(\S+)') {
                        $mountRoot = $matches[1].TrimEnd('/')
                    }
                }
                return ("$mountRoot/$driveLetter/$restOfPath" -replace '/$', '')
            }

            # Last resort on native Windows with no WSL: legacy Docker Desktop
            # short form. Works on Hyper-V-backend / pre-WSL2 setups.
            Write-Warning "wsl.exe not found; falling back to Docker Desktop legacy '/${driveLetter}/...' bind-mount form. If your daemon is WSL-backed this may not resolve correctly."
            return ("/$driveLetter/$restOfPath" -replace '/$', '')
        }

        Write-Warning "Path '$Path' did not match an expected Windows or POSIX absolute path."
        return $Path
    }

    Write-Host "`nUsing direct bind mounts for project, outputs, and synthpop..."

    # Convert paths for Docker bind mount
    $DockerProjectRoot = Convert-PathToDockerFormat -Path $ProjectRoot
    $DockerOutputDir = Convert-PathToDockerFormat -Path $outputDir
    $DockerSynthpopDir = Convert-PathToDockerFormat -Path $synthpopDir

    Write-Host "Docker Project Root: $DockerProjectRoot"
    Write-Host "Docker Output Dir:   $DockerOutputDir"
    Write-Host "Docker Synthpop Dir: $DockerSynthpopDir"

    # Pass mount arguments correctly to docker run
    docker run -it --rm `
        -e "USER_ID=$UserId" `
        -e "GROUP_ID=$GroupId" `
        -e "USER_NAME=$UserName" `
        -e "GROUP_NAME=$GroupName" `
        --mount "type=bind,source=$DockerProjectRoot,target=/IMPACTncd_England" `
        --mount "type=bind,source=$DockerOutputDir,target=/outputs" `
        --mount "type=bind,source=$DockerSynthpopDir,target=/synthpop" `
        --workdir /IMPACTncd_England `
        $ImageName `
        bash
}

# Restore the original directory
Pop-Location