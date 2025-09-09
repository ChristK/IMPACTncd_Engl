#!/bin/bash

# -----------------------------------------------------------------------------
# create_env.sh
#
# Optional argument: path to YAML file
# Default: ../inputs/sim_design.yaml
# Usage: ./create_env.sh [path_to_yaml]
#
# Builds and runs a Docker container for the IMPACTncd Engl project.
# Rebuilds the image only if relevant input files have changed.
# This script is designed to be run from the docker_setup directory, and binds 
# its parent directory to /IMPACTncd_Engl in the container.
# Compatible with Linux and macOS (requires coreutils on macOS).
# -----------------------------------------------------------------------------

IMAGE_NAME="impactncd-engl-r-prerequisite:latest"
DOCKERFILE="Dockerfile.prerequisite"
HASH_FILE=".docker_build_hash"

YAML_FILE="${1:-../inputs/sim_design.yaml}"

if [ ! -f "$YAML_FILE" ]; then
  echo "Error: YAML file not found at $YAML_FILE"
  exit 1
fi

echo "Using configuration from: $YAML_FILE"

# Path to the YAML file (adjust if needed)
SIM_DESIGN_FILE="$YAML_FILE"

# Extract output_dir and synthpop_dir using grep + sed
OUTPUT_DIR=$(grep '^output_dir:' "$SIM_DESIGN_FILE" | sed -E 's/output_dir:[[:space:]]*([^#]*).*/\1/' | xargs)
SYNTHPOP_DIR=$(grep '^synthpop_dir:' "$SIM_DESIGN_FILE" | sed -E 's/synthpop_dir:[[:space:]]*([^#]*).*/\1/' | xargs)
# Convert Windows-style or escaped paths if needed (only if using this in mixed environments)
# Make sure the paths are absolute
if [[ "$OUTPUT_DIR" != /* ]]; then
  OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"
fi

if [[ "$SYNTHPOP_DIR" != /* ]]; then
  SYNTHPOP_DIR="$(realpath "$SYNTHPOP_DIR")"
fi

echo "Mounting output_dir: $OUTPUT_DIR"
echo "Mounting synthpop_dir: $SYNTHPOP_DIR"

# Detect OS and choose hash command
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  HASH_CMD="sha256sum"
elif [[ "$OSTYPE" == "darwin"* ]]; then
  # macOS: Use gsha256sum from coreutils (brew install coreutils)
  if command -v gsha256sum > /dev/null; then
    HASH_CMD="gsha256sum"
  else
    echo "Error: gsha256sum not found. Please install coreutils with:"
    echo "  brew install coreutils"
    exit 1
  fi
else
  echo "Unsupported OS: $OSTYPE"
  exit 1
fi

# Compute content hash of Dockerfile and package lists
BUILD_HASH=$(cat "$DOCKERFILE" apt-packages.txt r-packages.txt | $HASH_CMD | cut -d ' ' -f1)

# Determine whether rebuild is needed
NEEDS_BUILD=false
if ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
  echo "Docker image does not exist. Need to build."
  NEEDS_BUILD=true
elif [ ! -f "$HASH_FILE" ]; then
  echo "No previous build hash found. Need to build."
  NEEDS_BUILD=true
else
  LAST_HASH=$(cat "$HASH_FILE")
  if [ "$BUILD_HASH" != "$LAST_HASH" ]; then
    echo "Detected changes in build inputs. Rebuilding Docker image..."
    NEEDS_BUILD=true
  else
    echo "No changes detected. Skipping Docker build."
  fi
fi

# Build Docker image if needed
if [ "$NEEDS_BUILD" = true ]; then
  docker build --no-cache -f "$DOCKERFILE" -t "$IMAGE_NAME" .
  echo "$BUILD_HASH" > "$HASH_FILE"
fi

# Run Docker container interactively with bind mount
docker run -it \
  --mount type=bind,source="$(pwd)/..",target=/IMPACTncd_Engl \
  --mount type=bind,source="$OUTPUT_DIR",target=/IMPACTncd_Engl/outputs \
  --mount type=bind,source="$SYNTHPOP_DIR",target=/IMPACTncd_Engl/synthpop \
  "$IMAGE_NAME" \
  bash