<file name=build_and_push_IMPACTncd.ps1 path=/home/ckyprid/My_Models/IMPACTncd/docker_setup><# 
.SYNOPSIS 
    This script automates building and pushing a Docker image to Docker Hub on Windows. 

.DESCRIPTION 
    The script performs the following steps: 
      - Loads environment variables from a .env file if present. 
      - Prompts for Docker Hub username if not set. 
      - Logs into Docker Hub using the DOCKERHUB_TOKEN, or asks for manual login if not provided. 
      - Builds the Docker image from Dockerfile.IMPACTncd using the parent directory as context. 
      - Optionally pushes the image to Docker Hub if the --push flag is provided. 

.PARAMETER --push 
    Optional flag to push the Docker image after building it. 

.EXAMPLE 
    .\build_and_push_IMPACTncd.ps1 --push 
    This will build the image and push it to Docker Hub. 
#> 

# Enable strict mode for better scripting practices 
Set-StrictMode -Version Latest 

# Determine if the --push flag is provided 
$PushImage = $false 
if ($args.Length -gt 0 -and $args[0] -eq '--push') { 
    $PushImage = $true 
} 

# Function for timestamped logging 
function Log($message) { 
    $timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss" 
    Write-Host "$timestamp - $message" 
} 

# Function to load environment variables from a .env file 
function Load-DotEnv { 
    $envPath = Join-Path -Path (Get-Location) -ChildPath ".env" 
    if (Test-Path $envPath) { 
        Log "Sourcing environment variables from .env" 
        Get-Content $envPath | ForEach-Object { 
            # Skip comments and empty lines 
            if ($_ -match '^\s*#' -or $_.Trim() -eq '') { return } 
            if ($_ -match '^(\w+)=(.*)$') { 
                $key = $matches[1].Trim() 
                $value = $matches[2].Trim() 
                # Remove any surrounding quotes from the value 
                $value = $value.Trim('"').Trim("'") 
                [System.Environment]::SetEnvironmentVariable($key, $value) 
            } 
        } 
    } 
} 

# Load environment variables if .env exists 
Load-DotEnv 

# Prompt for Docker Hub username if not already set 
if (-not $env:DOCKERHUB_USERNAME) { 
    $env:DOCKERHUB_USERNAME = Read-Host "Enter your Docker Hub username" 
} 

# Define the image name 
$IMAGE_NAME = "$($env:DOCKERHUB_USERNAME)/impactncd_engl:latest" 

# Docker login 
Log "Logging into Docker Hub..." 
if (-not $env:DOCKERHUB_USERNAME -or -not $env:DOCKERHUB_TOKEN) { 
    Log "Environment variables DOCKERHUB_USERNAME and/or DOCKERHUB_TOKEN not set." 
    Log "Please log in manually." 
    docker login 
} else { 
    try { 
        # Pipe the token into the docker login command 
        $env:DOCKERHUB_TOKEN | docker login -u $env:DOCKERHUB_USERNAME --password-stdin 
        if ($LASTEXITCODE -eq 0) { 
            Log "Docker login successful." 
        } else { 
            Log "Docker login failed." 
            exit 1 
        } 
    } catch { 
        Log "Docker login encountered an error: $_" 
        exit 1 
    } 
} 

# Build the Docker image using the parent directory as the build context 
Log "Building Docker image from Dockerfile.IMPACTncd..." 
try { 
    docker build --no-cache -f Dockerfile.IMPACTncd -t $IMAGE_NAME .. 
    if ($LASTEXITCODE -eq 0) { 
        Log "Docker image built successfully." 
    } else { 
        Log "Docker image build failed." 
        exit 1 
    } 
} catch { 
    Log "Docker build encountered an error: $_" 
    exit 1 
} 

# Push the Docker image if the --push flag was provided 
if ($PushImage) { 
    Log "Pushing Docker image to Docker Hub..." 
    try { 
        docker push $IMAGE_NAME 
        if ($LASTEXITCODE -eq 0) { 
            Log "Docker image pushed successfully." 
        } else { 
            Log "Docker push failed." 
            exit 1 
        } 
    } catch { 
        Log "Docker push encountered an error: $_" 
        exit 1 
    } 
} else { 
    Log "Skipping Docker push as --push argument not provided." 
} 
</file>