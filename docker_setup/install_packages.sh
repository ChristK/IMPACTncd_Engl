#!/bin/bash
set -e

updated_packages=""
failed_packages=""

echo "Starting intelligent package installation..."

# Read packages from apt-packages.txt
while IFS= read -r line || [[ -n "$line" ]]; do
  # Skip empty lines and comments
  [[ -z "$line" || "$line" =~ ^[[:space:]]*# ]] && continue
  
  package_spec="$line"
  package_name=$(echo "$package_spec" | cut -d'=' -f1)
  requested_version=$(echo "$package_spec" | cut -d'=' -f2-)
  
  echo "Processing package: $package_name (requested version: $requested_version)"
  
  # Try to install the specific version first
  if apt-get install -y --no-install-recommends "$package_spec" 2>/dev/null; then
    echo "✓ Successfully installed $package_spec"
  else
    echo "⚠ Version $requested_version for $package_name not available"
    
    # Get the available version
    available_version=$(apt-cache policy "$package_name" | grep "Candidate:" | awk '{print $2}')
    
    if [[ -n "$available_version" && "$available_version" != "(none)" ]]; then
      echo "→ Installing available version: $package_name=$available_version"
      
      # Install the available version
      if apt-get install -y --no-install-recommends "$package_name=$available_version"; then
        echo "✓ Successfully installed $package_name=$available_version"
        updated_packages="$updated_packages$package_name=$available_version\n"
      else
        echo "✗ Failed to install $package_name=$available_version"
        failed_packages="$failed_packages$package_name\n"
      fi
    else
      echo "→ Installing latest version of $package_name"
      if apt-get install -y --no-install-recommends "$package_name"; then
        installed_version=$(dpkg-query -W -f='${Version}' "$package_name" 2>/dev/null || echo "unknown")
        echo "✓ Successfully installed $package_name (version: $installed_version)"
        updated_packages="$updated_packages$package_name=$installed_version\n"
      else
        echo "✗ Failed to install $package_name"
        failed_packages="$failed_packages$package_name\n"
      fi
    fi
  fi
done < /tmp/apt-packages.txt

# Install gosu
echo "Installing gosu..."
apt-get install -y --no-install-recommends gosu

# Report results
if [[ -n "$updated_packages" ]]; then
  echo ""
  echo "=================================================="
  echo "PACKAGE VERSION UPDATES DETECTED:"
  echo "=================================================="
  echo "The following packages were installed with different versions than requested:"
  echo -e "$updated_packages"
  echo ""
  echo "RECOMMENDED ACTION: Update your apt-packages.txt file with these versions:"
  echo -e "$updated_packages"
  echo "=================================================="
fi

if [[ -n "$failed_packages" ]]; then
  echo ""
  echo "=================================================="
  echo "PACKAGE INSTALLATION FAILURES:"
  echo "=================================================="
  echo "The following packages failed to install:"
  echo -e "$failed_packages"
  echo "=================================================="
  exit 1
fi

echo "Package installation completed successfully!"
