#!/bin/bash
# -----------------------------------------------------------------------------
# setup_user_docker_env.sh
#
# Usage:
#   ./setup_user_docker_env.sh [-Tag <tag>] [-ScenariosDir <path_to_scenarios>] [-SimDesignYaml <path_to_yaml>] [--UseVolumes]
#
# Description:
#   This script pulls and runs a Docker container for the IMPACTncd England project.
#
# Container Selection:
#   - If <tag> is "main" (default): pulls and uses "chriskypri/impactncdengl:main".
#   - If <tag> is "local": uses "impactncdengl:local" (built locally).
#   - If <tag> is any other value: pulls and uses "chriskypri/impactncdengl:<tag>".
#
# Scenarios Directory:
#   - If [path_to_scenarios] is provided, it will be mounted as /IMPACTncd_England/scenarios inside the container.
#
# Operation Modes:
# 1. Using Docker-managed volumes (recommended for macOS and Windows):
#      - Creates Docker volumes for output_dir and synthpop_dir (defined in YAML).
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
#
# Examples:
#
# 1. Run with the default tag ("main") and default YAML file:
#    ./setup_user_docker_env.sh
#
# 2. Run with a specific tag (e.g., "v1.2.3") and default YAML file:
#    ./setup_user_docker_env.sh -Tag v1.2.3
#
# 3. Run with a custom scenarios directory and default YAML file:
#    ./setup_user_docker_env.sh -Tag main -ScenariosDir /path/to/scenarios
#
# 4. Run with a custom YAML file:
#    ./setup_user_docker_env.sh -Tag main -ScenariosDir /path/to/scenarios -SimDesignYaml /path/to/custom_sim_design.yaml
#
# 5. Use Docker volumes for better performance:
#    ./setup_user_docker_env.sh -Tag main -ScenariosDir /path/to/scenarios -SimDesignYaml /path/to/custom_sim_design.yaml --UseVolumes
# -----------------------------------------------------------------------------

# Get the directory where the script is located
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
# Get the project root directory (one level above the script directory)
PROJECT_ROOT=$(realpath "$SCRIPT_DIR/..")

# Variable definitions
DOCKER_TAG="main"  # Default tag
YAML_FILE="$PROJECT_ROOT/inputs/sim_design.yaml" # Default YAML path relative to project root
SCENARIOS_DIR=""  # Scenarios directory to copy into container
CURRENT_USER=$(whoami)
# Get current user's UID and GID for running containers as non-root
USER_ID=$(id -u)
GROUP_ID=$(id -g)
USER_NAME=$(whoami)
GROUP_NAME=$(id -gn)
# User-specific Docker volume names to avoid conflicts (only for output and synthpop)
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

# Update argument parsing to match PowerShell version
USE_VOLUMES=false # Default to not using volumes

# Process command-line arguments for scenarios folder, YAML file, volume usage flag, and Docker tag
while [[ $# -gt 0 ]]; do
  case $1 in
    -Tag)
      DOCKER_TAG="$2"
      shift 2
      ;;
    -ScenariosDir)
      if [ ! -d "$2" ]; then
        echo "Error: Scenarios directory not found at '$2'"
        exit 1
      fi
      SCENARIOS_DIR="$(realpath "$2")"
      shift 2
      ;;
    -SimDesignYaml)
      if [ ! -f "$2" ]; then
        echo "Error: YAML file not found at '$2'"
        exit 1
      fi
      YAML_FILE="$(realpath "$2")"
      shift 2
      ;;
    --UseVolumes)
      USE_VOLUMES=true
      shift
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

# Ensure default values if arguments are not provided
DOCKER_TAG=${DOCKER_TAG:-"main"}
YAML_FILE=${YAML_FILE:-"$PROJECT_ROOT/inputs/sim_design.yaml"}

# Determine the Docker image name based on the tag
if [[ "$DOCKER_TAG" == "local" ]]; then
  IMAGE_NAME="impactncdengl:local"
else
  IMAGE_NAME="chriskypri/impactncdengl:${DOCKER_TAG}"
fi

echo "Using Docker image: $IMAGE_NAME"

if [ ! -f "$YAML_FILE" ]; then
  echo "Error: YAML file not found at $YAML_FILE"
  exit 1
fi

echo "Using configuration from: $YAML_FILE"

if [[ -n "$SCENARIOS_DIR" ]]; then
  if [ ! -d "$SCENARIOS_DIR" ]; then
    echo "Error: Scenarios directory not found at $SCENARIOS_DIR"
    exit 1
  fi
  echo "Using scenarios from: $SCENARIOS_DIR"
fi

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

# Pull the Docker image (skip for local builds)
if [[ "$DOCKER_TAG" == "local" ]]; then
  echo "Using local Docker image: $IMAGE_NAME"
  if ! docker image inspect "$IMAGE_NAME" > /dev/null 2>&1; then
    echo "Error: Local image '$IMAGE_NAME' not found."
    echo "Build it first with: docker build -t $IMAGE_NAME -f docker_setup/Dockerfile.IMPACTncdENGL docker_setup"
    exit 1
  fi
else
  echo "Pulling Docker image: $IMAGE_NAME"
  if ! docker pull "$IMAGE_NAME"; then
    echo "Error: Failed to pull Docker image: $IMAGE_NAME"
    echo "Please check:"
    echo "  1. The image exists and is accessible"
    echo "  2. You have the correct permissions"
    echo "  3. Your internet connection is working"
    echo "  4. The tag '$DOCKER_TAG' exists in the chriskypri/impactncdengl repository"
    exit 1
  fi
fi

# -----------------------------------------------------------------------------
# Optionally create and use Docker volumes for simulation
#
# When using volumes:
#   - Separate volumes (VOLUME_OUTPUT_NAME and VOLUME_SYNTHPOP_NAME) for the outputs 
#     and synthpop directories (as specified in the YAML file) are created.
#   - Prior to simulation, the local outputs and synthpop folders are copied into these volumes.
#   - The container runs with these Docker volumes mounted. This improves I/O performance.
#   - After the container exits, the content of the output and synthpop volumes is 
#     synchronized back to the corresponding local folders using rsync.
#   - Finally, all these Docker volumes are removed to clean up.
#
# When not using volumes, the script uses direct bind mounts for output and synthpop directories.
#
# Note: The Docker image already contains the /IMPACTncd_England project, so no project
# volume or bind mount is needed.
# -----------------------------------------------------------------------------
if [ "$USE_VOLUMES" = true ]; then
  echo "Using Docker volumes for outputs and synthpop..."

  # Build rsync-alpine image (only if it doesn't already exist)
  if ! docker image inspect rsync-alpine > /dev/null 2>&1; then
    echo "Building rsync-alpine image..."
    
    # Check if Dockerfile.rsync exists
    if [ -f "Dockerfile.rsync" ]; then
      echo "Using Dockerfile.rsync..."
      docker build -f Dockerfile.rsync -t rsync-alpine .
    else
      echo "Dockerfile.rsync not found, creating rsync image inline..."
      cat << 'EOF' | docker build -t rsync-alpine -
FROM alpine:latest
RUN apk add --no-cache rsync
EOF
    fi
  else
    echo "Using existing rsync-alpine image."
  fi

  # Ensure local output directories exist
  mkdir -p "$OUTPUT_DIR"
  mkdir -p "$SYNTHPOP_DIR"

  # Remove any existing volumes (ignore errors if not removable)
  echo "Removing any existing volumes (if possible)..."
  docker volume rm "$VOLUME_OUTPUT_NAME" 2>/dev/null
  docker volume rm "$VOLUME_SYNTHPOP_NAME" 2>/dev/null

  # Create fresh Docker-managed volumes
  docker volume create "$VOLUME_OUTPUT_NAME"
  docker volume create "$VOLUME_SYNTHPOP_NAME"

  # --------------------------------------------------------------------------
  # Fix volume ownership and pre-populate volumes:
  #
  # Docker volumes are created with root ownership by default. We need to fix
  # the ownership before we can populate them as the calling user.
  #
  # The output and synthpop volumes are populated from the respective local folders.
  # --------------------------------------------------------------------------
  
  # Fix ownership of volume directories first (run as root, then change ownership)
  echo "Setting correct ownership for Docker volumes..."
  docker run --rm \
    -v "$VOLUME_OUTPUT_NAME":/volume \
    alpine sh -c "chown ${USER_ID}:${GROUP_ID} /volume"
  docker run --rm \
    -v "$VOLUME_SYNTHPOP_NAME":/volume \
    alpine sh -c "chown ${USER_ID}:${GROUP_ID} /volume"

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

  # Run the main container using the pre-built image
  echo "Running the main container using Docker volumes..."
  
  # Prepare docker run command with scenarios mount if provided
  if [[ -n "$SCENARIOS_DIR" ]]; then
    docker run -it --rm \
      -e USER_ID="${USER_ID}" \
      -e GROUP_ID="${GROUP_ID}" \
      -e USER_NAME="${USER_NAME}" \
      -e GROUP_NAME="${GROUP_NAME}" \
      --mount type=volume,source="$VOLUME_OUTPUT_NAME",target=/outputs \
      --mount type=volume,source="$VOLUME_SYNTHPOP_NAME",target=/synthpop \
      --mount type=bind,source="$SCENARIOS_DIR",target=/IMPACTncd_England/scenarios \
      --workdir /IMPACTncd_England \
      "$IMAGE_NAME" \
      bash
  else
    docker run -it --rm \
      -e USER_ID="${USER_ID}" \
      -e GROUP_ID="${GROUP_ID}" \
      -e USER_NAME="${USER_NAME}" \
      -e GROUP_NAME="${GROUP_NAME}" \
      --mount type=volume,source="$VOLUME_OUTPUT_NAME",target=/outputs \
      --mount type=volume,source="$VOLUME_SYNTHPOP_NAME",target=/synthpop \
      --workdir /IMPACTncd_England \
      "$IMAGE_NAME" \
      bash
  fi

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

  # Clean up all the Docker volumes used for the simulation.
  echo "Cleaning up Docker volumes..."
  docker volume rm "$VOLUME_OUTPUT_NAME"
  docker volume rm "$VOLUME_SYNTHPOP_NAME"
else
  echo "Using direct bind mounts for outputs and synthpop..."

  # Prepare docker run command with scenarios mount if provided
  if [[ -n "$SCENARIOS_DIR" ]]; then
    docker run -it --rm \
      -e USER_ID="${USER_ID}" \
      -e GROUP_ID="${GROUP_ID}" \
      -e USER_NAME="${USER_NAME}" \
      -e GROUP_NAME="${GROUP_NAME}" \
      --mount type=bind,source="$OUTPUT_DIR",target=/outputs \
      --mount type=bind,source="$SYNTHPOP_DIR",target=/synthpop \
      --mount type=bind,source="$SCENARIOS_DIR",target=/IMPACTncd_England/scenarios \
      --workdir /IMPACTncd_England \
      "$IMAGE_NAME" \
      bash
  else
    docker run -it --rm \
      -e USER_ID="${USER_ID}" \
      -e GROUP_ID="${GROUP_ID}" \
      -e USER_NAME="${USER_NAME}" \
      -e GROUP_NAME="${GROUP_NAME}" \
      --mount type=bind,source="$OUTPUT_DIR",target=/outputs \
      --mount type=bind,source="$SYNTHPOP_DIR",target=/synthpop \
      --workdir /IMPACTncd_England \
      "$IMAGE_NAME" \
      bash
  fi
fi