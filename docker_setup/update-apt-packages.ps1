#!/usr/bin/env pwsh

<#
.SYNOPSIS
    Updates apt-packages.txt with the actual installed versions from Docker build output.

.DESCRIPTION
    This script parses Docker build output to find package version updates and 
    automatically updates the apt-packages.txt file with the correct versions.

.PARAMETER BuildLogFile
    Path to a file containing Docker build output (optional)

.PARAMETER Interactive
    Run in interactive mode to manually select which versions to update

.EXAMPLE
    .\update-apt-packages.ps1
    Parse Docker build output from pipeline and update apt-packages.txt

.EXAMPLE
    .\update-apt-packages.ps1 -BuildLogFile "build.log"
    Parse specific build log file and update apt-packages.txt

.EXAMPLE
    .\update-apt-packages.ps1 -Interactive
    Interactively choose which versions to update
#>

param(
    [string]$BuildLogFile,
    [switch]$Interactive
)

$AptPackagesFile = Join-Path $PSScriptRoot "apt-packages.txt"

function Parse-DockerBuildOutput {
    param(
        [string[]]$BuildOutput
    )
    
    $updates = @{}
    $inUpdateSection = $false
    
    foreach ($line in $BuildOutput) {
        if ($line -match "PACKAGE VERSION UPDATES DETECTED:") {
            $inUpdateSection = $true
            continue
        }
        
        if ($inUpdateSection) {
            if ($line -match "=============") {
                if ($line -match "RECOMMENDED ACTION:") {
                    continue
                } else {
                    $inUpdateSection = $false
                    continue
                }
            }
            
            if ($line -match "^([a-zA-Z0-9\-\.]+)=(.+)$") {
                $packageName = $matches[1]
                $version = $matches[2]
                $updates[$packageName] = $version
                Write-Host "Found update: $packageName = $version" -ForegroundColor Green
            }
        }
    }
    
    return $updates
}

function Update-AptPackagesFile {
    param(
        [hashtable]$Updates
    )
    
    if (-not (Test-Path $AptPackagesFile)) {
        Write-Error "apt-packages.txt not found at: $AptPackagesFile"
        return
    }
    
    $content = Get-Content $AptPackagesFile
    $newContent = @()
    $updated = $false
    
    foreach ($line in $content) {
        if ($line -match "^([a-zA-Z0-9\-\.]+)=") {
            $packageName = $matches[1]
            if ($Updates.ContainsKey($packageName)) {
                $newVersion = $Updates[$packageName]
                $newLine = "$packageName=$newVersion"
                Write-Host "Updating $packageName to version $newVersion" -ForegroundColor Yellow
                $newContent += $newLine
                $updated = $true
            } else {
                $newContent += $line
            }
        } else {
            $newContent += $line
        }
    }
    
    if ($updated) {
        $backupFile = "$AptPackagesFile.backup.$(Get-Date -Format 'yyyyMMdd_HHmmss')"
        Copy-Item $AptPackagesFile $backupFile
        Write-Host "Backup created: $backupFile" -ForegroundColor Cyan
        
        $newContent | Out-File -FilePath $AptPackagesFile -Encoding UTF8
        Write-Host "Updated $AptPackagesFile successfully!" -ForegroundColor Green
    } else {
        Write-Host "No updates were needed." -ForegroundColor Blue
    }
}

# Main execution
try {
    Write-Host "APT Packages Version Updater" -ForegroundColor Magenta
    Write-Host "=============================" -ForegroundColor Magenta
    
    if ($BuildLogFile) {
        if (-not (Test-Path $BuildLogFile)) {
            Write-Error "Build log file not found: $BuildLogFile"
            exit 1
        }
        $buildOutput = Get-Content $BuildLogFile
        Write-Host "Reading from build log file: $BuildLogFile" -ForegroundColor Cyan
    } else {
        Write-Host "Reading from pipeline input..." -ForegroundColor Cyan
        $buildOutput = @()
        while ($input = Read-Host -Prompt "Enter Docker build output (or 'END' to finish)") {
            if ($input -eq "END") { break }
            $buildOutput += $input
        }
    }
    
    $updates = Parse-DockerBuildOutput -BuildOutput $buildOutput
    
    if ($updates.Count -eq 0) {
        Write-Host "No package version updates found in the build output." -ForegroundColor Yellow
        exit 0
    }
    
    Write-Host "`nFound $($updates.Count) package updates:" -ForegroundColor Cyan
    foreach ($pkg in $updates.Keys) {
        Write-Host "  $pkg = $($updates[$pkg])" -ForegroundColor White
    }
    
    if ($Interactive) {
        Write-Host "`nReview the updates above." -ForegroundColor Cyan
        $confirm = Read-Host "Do you want to apply these updates to apt-packages.txt? (y/N)"
        if ($confirm -ne "y" -and $confirm -ne "Y") {
            Write-Host "Updates cancelled by user." -ForegroundColor Yellow
            exit 0
        }
    }
    
    Update-AptPackagesFile -Updates $updates
    
} catch {
    Write-Error "An error occurred: $($_.Exception.Message)"
    exit 1
}
