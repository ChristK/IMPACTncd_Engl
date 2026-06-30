#!/bin/bash
# -----------------------------------------------------------------------------
# setup_dev_docker_env.sh
#
# Usage:
#   ./setup_dev_docker_env.sh [path_to_yaml] [--UseVolumes]
#
# Description:
#   This script builds and runs a Docker container for the IMPACTncd England project.
#   It rebuilds the Docker image only if build inputs have changed.
#
# Container Selection:
#   - Builds a local development image: prerequisite.impactncdengl:local (from Dockerfile.prerequisite.IMPACTncdENGL).
#   - Automatically detects changes in build inputs (Dockerfile, apt-packages.txt, etc.) and rebuilds the image if needed.
#
# Operation Modes:
# 1. Using Docker-managed volumes (recommended for macOS and Windows):
#      - Copies the project directory into a Docker volume for faster I/O.
#      - Creates separate volumes for output_dir and synthpop_dir (defined in YAML).
#      - Synchronizes volumes back to local folders after container exits.
#      - Removes volumes after synchronization.
#
# 2. Using direct bind mounts (less efficient, but useful for interactive access):
#      - Mounts local directories directly into the container.
#
# Security:
#   - Containers run as the calling user (non-root) to prevent permission issues.
#   - Automatically detects the current user's UID and GID and passes them to Docker.
#
# Notes:
# - Compatible with Linux and macOS (requires coreutils on macOS).
# - For macOS and Windows, using Docker volumes is recommended for better performance.
# - For Linux, ensure your user has Docker permissions (e.g., part of the "docker" group).
# - If you encounter permission issues, ensure the output_dir and synthpop_dir exist and are writable.
# -----------------------------------------------------------------------------

# Get the directory where the script is located
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the project root directory (one level above the script directory)
PROJECT_ROOT=$(realpath "$SCRIPT_DIR/..")

# Variable definitions
IMAGE_NAME="prerequisite.impactncdengl:local"
DOCKERFILE="Dockerfile.prerequisite.IMPACTncdENGL"
HASH_FILE="$SCRIPT_DIR/.docker_build_hash" # Store hash file in script directory
YAML_FILE="$PROJECT_ROOT/inputs/sim_design.yaml" # Default YAML path relative to project root
CURRENT_USER=$(whoami)
# Get current user's UID and GID for running containers as non-root
USER_ID=$(id -u)
GROUP_ID=$(id -g)
USER_NAME=$(whoami)
GROUP_NAME=$(id -gn)
# User-specific Docker volume names to avoid conflicts
VOLUME_PROJECT="impactncd_england_project_${CURRENT_USER}"
VOLUME_OUTPUT_NAME="impactncd_england_output_${CURRENT_USER}"
VOLUME_SYNTHPOP_NAME="impactncd_england_synthpop_${CURRENT_USER}"

# --- Docker Permission Check ---
# Check if the user can connect to the Docker daemon
if ! docker info > /dev/null 2>&1; then
  echo "---------------------------------------------------------------------"
  echo "Error: Cannot connect to the Docker daemon."
  echo "Please ensure Docker is running and you have the necessary permissions."
  echo "You might need to:"
  echo "  1. Start the Docker daemon."
  echo "  2. Add your user to the 'docker' group:"
  echo "     sudo usermod -aG docker $USER"
  echo "     (You'll need to log out and back in for this change to take effect)"
  echo "  3. Or run this script using 'sudo':"
  echo "     sudo ./setup_dev_docker_env.sh [options]"
  echo "---------------------------------------------------------------------"
  exit 1
fi
# --- End Docker Permission Check ---
# Process command-line arguments for YAML file and volume usage flag
USE_VOLUMES=false # Default to not using volumes
for arg in "$@"; do
  if [[ "$arg" == --UseVolumes ]]; then
    USE_VOLUMES=true
  elif [[ "$arg" == *.yaml || "$arg" == *.yml ]]; then
    # If YAML path is absolute, use it as-is
    if [[ "$arg" == /* || "$arg" == ~* ]]; then
        YAML_FILE="$arg"
    else
        # If YAML path is relative, resolve it relative to the project root
        YAML_FILE="$(realpath "$PROJECT_ROOT/$arg")"
    fi
  fi
done

# Default to not using volumes if flag is not provided
USE_VOLUMES=false

# Check again for --UseVolumes flag
for arg in "$@"; do
  if [[ "$arg" == "--UseVolumes" ]]; then
    USE_VOLUMES=true
  fi
done

if [ ! -f "$YAML_FILE" ]; then
  echo "Error: YAML file not found at $YAML_FILE"
  echo "Project root: $PROJECT_ROOT"
  exit 1
fi

echo "Using configuration from: $YAML_FILE"

# Set simulation design file and extract output directories from YAML
SIM_DESIGN_FILE="$YAML_FILE"
OUTPUT_DIR_RAW=$(grep '^output_dir:' "$SIM_DESIGN_FILE" | sed -E 's/output_dir:[[:space:]]*([^#]*).*/\1/' | xargs)
SYNTHPOP_DIR_RAW=$(grep '^synthpop_dir:' "$SIM_DESIGN_FILE" | sed -E 's/synthpop_dir:[[:space:]]*([^#]*).*/\1/' | xargs)

# Resolve paths relative to the PROJECT_ROOT if they are not absolute
if [[ "$OUTPUT_DIR_RAW" != /* && "$OUTPUT_DIR_RAW" != ~* ]]; then
  OUTPUT_DIR_TEMP="$PROJECT_ROOT/$OUTPUT_DIR_RAW"
else
  OUTPUT_DIR_TEMP="$OUTPUT_DIR_RAW"
fi
if [[ "$SYNTHPOP_DIR_RAW" != /* && "$SYNTHPOP_DIR_RAW" != ~* ]]; then
  SYNTHPOP_DIR_TEMP="$PROJECT_ROOT/$SYNTHPOP_DIR_RAW"
else
  SYNTHPOP_DIR_TEMP="$SYNTHPOP_DIR_RAW"
fi

# Create directories if they don't exist, then resolve with realpath
mkdir -p "$OUTPUT_DIR_TEMP"
mkdir -p "$SYNTHPOP_DIR_TEMP"
OUTPUT_DIR="$(realpath "$OUTPUT_DIR_TEMP")"
SYNTHPOP_DIR="$(realpath "$SYNTHPOP_DIR_TEMP")"

echo "Mounting output_dir: $OUTPUT_DIR"
echo "Mounting synthpop_dir: $SYNTHPOP_DIR"

# Detect OS and choose appropriate hash command
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
  HASH_CMD="sha256sum"
elif [[ "$OSTYPE" == "darwin"* ]]; then
  if command -v gsha256sum > /dev/null; then
    HASH_CMD="gsha256sum"
  else
    echo "Error: gsha256sum not found. Please install coreutils."
    echo "If you use Homebrew, run: brew install coreutils"
    exit 1
  fi
else
  echo "Unsupported OS: $OSTYPE"
  exit 1
fi

# Compute hash of build inputs (Dockerfile, apt-packages.txt, r-packages.txt, entrypoint.sh)
BUILD_HASH=$(cat "$SCRIPT_DIR/$DOCKERFILE" "$SCRIPT_DIR/apt-packages.txt" "$SCRIPT_DIR/r-packages.txt" "$SCRIPT_DIR/entrypoint.sh" | $HASH_CMD | cut -d ' ' -f1)

# Determine whether rebuild of the Docker image is needed
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

# -----------------------------------------------------------------------------
# Optionally create and use Docker volumes for simulation
#
# When using volumes:
#
#   - The entire project directory (parent of docker_setup) is copied to a Docker
#     volume (VOLUME_PROJECT). This volume acts as the main drive during simulation.
#
#   - Separate volumes (VOLUME_OUTPUT_NAME and VOLUME_SYNTHPOP_NAME) for the outputs 
#     and synthpop directories (as specified in the YAML file) are created.
#
#   - Prior to simulation, the local outputs and synthpop folders are copied into these volumes.
#
#   - The container runs with these Docker volumes mounted. This improves I/O performance.
#
#   - After the container exits, the content of the output and synthpop volumes is 
#     synchronized back to the corresponding local folders using rsync.
#
#   - Finally, all these Docker volumes are removed to clean up.
#
# When not using volumes, the script uses direct bind mounts.
# -----------------------------------------------------------------------------
if [ "$USE_VOLUMES" = true ]; then
  echo "Using Docker volumes for project, outputs, and synthpop..."

  # Build rsync-alpine image (only if it doesn't already exist)
  if ! docker image inspect rsync-alpine > /dev/null 2>&1; then
    echo "Building rsync-alpine image..."
    docker build -f Dockerfile.rsync -t rsync-alpine .
  else
    echo "Using existing rsync-alpine image."
  fi

  # Ensure local output directories exist
  mkdir -p "$OUTPUT_DIR"
  mkdir -p "$SYNTHPOP_DIR"

  # Remove any existing volumes (ignore errors if not removable)
  echo "Removing any existing volumes (if possible)..."
  docker volume rm "$VOLUME_PROJECT" 2>/dev/null
  docker volume rm "$VOLUME_OUTPUT_NAME" 2>/dev/null
  docker volume rm "$VOLUME_SYNTHPOP_NAME" 2>/dev/null

  # Create fresh Docker-managed volumes
  docker volume create "$VOLUME_PROJECT"
  docker volume create "$VOLUME_OUTPUT_NAME"
  docker volume create "$VOLUME_SYNTHPOP_NAME"

  # --------------------------------------------------------------------------
  # Fix volume ownership and pre-populate volumes:
  #
  # Docker volumes are created with root ownership by default. We need to fix
  # the ownership before we can populate them as the calling user.
  #
  # 1. The project volume is populated with the entire project directory (from
  #    one level above the docker_setup folder), excluding dot files/folders.
  #    This volume serves as the main drive.
  #
  # 2. The output and synthpop volumes are populated from the respective local folders.
  # --------------------------------------------------------------------------
  
  # Fix ownership of volume directories first (run as root, then change ownership)
  echo "Setting correct ownership for Docker volumes..."
  docker run --rm \
    -v "${VOLUME_PROJECT}":/destination \
    alpine sh -c "chown ${USER_ID}:${GROUP_ID} /destination"
  docker run --rm \
    -v "$VOLUME_OUTPUT_NAME":/volume \
    alpine sh -c "chown ${USER_ID}:${GROUP_ID} /volume"
  docker run --rm \
    -v "$VOLUME_SYNTHPOP_NAME":/volume \
    alpine sh -c "chown ${USER_ID}:${GROUP_ID} /volume"

  echo "Populating project volume from host project directory (excluding dot files/folders)..."
  # Use tar to copy, excluding dot files/folders at the root of the source
  docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" \
    -v "${PROJECT_ROOT}":/source \
    -v "${VOLUME_PROJECT}":/destination \
    alpine sh -c "tar -C /source --exclude='./.*' -cf - . | tar -C /destination -xf -"

  echo "Populating output and synthpop volumes from local folders..."
  docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" \
    -v "$OUTPUT_DIR":/source \
    -v "$VOLUME_OUTPUT_NAME":/volume \
    alpine sh -c "cp -r /source/. /volume/ 2>/dev/null || cp -a /source/. /volume/ 2>/dev/null || true"
  docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" \
    -v "$SYNTHPOP_DIR":/source \
    -v "$VOLUME_SYNTHPOP_NAME":/volume \
    alpine sh -c "cp -r /source/. /volume/ 2>/dev/null || cp -a /source/. /volume/ 2>/dev/null || true"

  # Run the main container with the project volume mounted in place of the project bind mount.
  echo "Running the main container using Docker volumes..."
  docker run -it --rm \
    -e USER_ID="${USER_ID}" \
    -e GROUP_ID="${GROUP_ID}" \
    -e USER_NAME="${USER_NAME}" \
    -e GROUP_NAME="${GROUP_NAME}" \
    --mount type=volume,source="$VOLUME_PROJECT",target=/IMPACTncd_England \
    --mount type=volume,source="$VOLUME_OUTPUT_NAME",target=/outputs \
    --mount type=volume,source="$VOLUME_SYNTHPOP_NAME",target=/synthpop \
    --workdir /IMPACTncd_England \
    "$IMAGE_NAME" \
    bash

  # After the container exits:
  # - Synchronize the volumes back to the local directories using rsync (checksum mode).
  echo "Container exited. Syncing volumes back to local directories using rsync (checksum mode)..."
  docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" \
    -v "$VOLUME_OUTPUT_NAME":/volume \
    -v "$OUTPUT_DIR":/backup \
    rsync-alpine rsync -avc --no-owner --no-group --no-times /volume/ /backup/
  docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" \
    -v "$VOLUME_SYNTHPOP_NAME":/volume \
    -v "$SYNTHPOP_DIR":/backup \
    rsync-alpine rsync -avc --no-owner --no-group --no-times /volume/ /backup/
  # Sync simulation folder back to the project directory
  SIMULATION_DIR="$PROJECT_ROOT/simulation"
  echo "Syncing simulation folder back to: $SIMULATION_DIR"
  docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" \
    -v "$VOLUME_PROJECT":/project \
    -v "$SIMULATION_DIR":/backup \
    rsync-alpine rsync -avc --no-owner --no-group --no-times /project/simulation/ /backup/

  # Clean up all the Docker volumes used for the simulation.
  echo "Cleaning up Docker volumes..."
  docker volume rm "$VOLUME_PROJECT"
  docker volume rm "$VOLUME_OUTPUT_NAME"
  docker volume rm "$VOLUME_SYNTHPOP_NAME"
else
  echo "Using direct bind mounts for project, outputs, and synthpop..."

  docker run -it --rm \
    -e USER_ID="${USER_ID}" \
    -e GROUP_ID="${GROUP_ID}" \
    -e USER_NAME="${USER_NAME}" \
    -e GROUP_NAME="${GROUP_NAME}" \
    --mount type=bind,source="$PROJECT_ROOT",target=/IMPACTncd_England \
    --mount type=bind,source="$OUTPUT_DIR",target=/outputs \
    --mount type=bind,source="$SYNTHPOP_DIR",target=/synthpop \
    --workdir /IMPACTncd_England \
    "$IMAGE_NAME" \
    bash
fi