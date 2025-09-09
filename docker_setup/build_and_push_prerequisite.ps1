<#
    This script automates building and optionally pushing a Docker image to Docker Hub.

    Usage:
    1. (Optional) Create a `.env` file in the same directory with the following variables:
           DOCKERHUB_USERNAME=yourusername
           DOCKERHUB_TOKEN=youraccesstoken

    2. Run the script as follows:
           .\build_and_push_prerequisite.ps1 [-Push]

       Use the -Push switch to push the Docker image to Docker Hub. Without the switch, the script will only build the image.
#>
param (
    [string]$DockerHubUsername = $env:DOCKERHUB_USERNAME,
    [string]$DockerHubToken = $env:DOCKERHUB_TOKEN,
    [switch]$Push = $false
)

# Load environment variables from .env if not already set
$envFile = ".env"
if (-not $DockerHubUsername -or -not $DockerHubToken) {
    if (Test-Path $envFile) {
        $lines = Get-Content $envFile | Where-Object { $_ -notmatch '^(\s*#|\s*$)' }
        if ($lines.Count -eq 0) {
            Write-Host ".env file is empty or contains only comments."
        } else {
            $lines | ForEach-Object {
                $parts = $_ -split '=', 2
                if ($parts.Length -eq 2) {
                    $key = $parts[0].Trim()
                    $value = $parts[1].Trim()
                    if (-not (Get-Item -Path Env:$key -ErrorAction SilentlyContinue)) {
                        [Environment]::SetEnvironmentVariable($key, $value, "Process")
                        Write-Host "Set $key from .env"
                    }
                }
            }
            Write-Host "Environment variables loaded from .env"
        }
    } else {
        Write-Host ".env file not found."
    }
}

function Log {
    param ($Message)
    Write-Output "$(Get-Date -Format 'yyyy-MM-dd HH:mm:ss') - $Message"
}

# Prompt for username if not set
if (-not $DockerHubUsername) {
    $DockerHubUsername = Read-Host "Enter your Docker Hub username"
}

$ImageName = "$DockerHubUsername/impactncd-engl-r-prerequisite:latest"

# Login
if (-not $DockerHubToken) {
    Log "DOCKERHUB_TOKEN not set. Prompting for manual login..."
    docker login
} else {
    Log "Logging into Docker Hub..."
    $DockerHubToken | docker login -u $DockerHubUsername --password-stdin
    if ($LASTEXITCODE -ne 0) {
        Log "Docker login failed."
        exit 1
    }
}

# Build
Log "Building Docker image..."
docker build --no-cache -f Dockerfile.prerequisite -t $ImageName .
if ($LASTEXITCODE -eq 0) {
    Log "Docker image built successfully."
} else {
    Log "Docker image build failed."
    exit 1
}

if ($Push) {
    Log "Pushing Docker image to Docker Hub..."
    docker push $ImageName
    if ($LASTEXITCODE -eq 0) {
        Log "Docker image pushed successfully."
    } else {
        Log "Docker push failed."
        exit 1
    }
} else {
    Log "Skipping Docker push as -Push switch was not provided."
}
