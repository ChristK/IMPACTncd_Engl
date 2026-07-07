#!/bin/bash
# -----------------------------------------------------------------------------
# setup_user_docker_env.sh
#
# Usage:
#   ./setup_user_docker_env.sh [-Tag <tag>] [-ScenariosDir <path_to_scenarios>] [-SimDesignYaml <path_to_yaml>] [-Run] [--UseVolumes]
#
# Description:
#   This script pulls and runs a Docker container for the IMPACTncd England project.
#
# Standalone ("drop-in folder") mode:
#   Download just this script into a folder that also holds your own
#   sim_design*.yaml (plus any scenarios / a driver R script), then run it with
#   NO arguments. It will:
#     - auto-detect the single sim_design*.yaml sitting next to it (if exactly
#       one YAML is present; otherwise pass -SimDesignYaml to choose), and
#     - merge that folder's files into the model folder /IMPACTncd_England
#       (via symlinks) so your config + scripts sit alongside global.R — you
#       reference them with no subfolder prefix and pass no paths.
#   Relative output_dir/synthpop_dir in the YAML are resolved relative to the
#   YAML's own folder.
#   Add -Run to auto-run the single simulate*.R driver in the folder via Rscript
#   and exit (otherwise you land in an interactive bash shell, as before).
#
#   If no YAML sits next to the script, it falls back to the repo default
#   inputs/sim_design.yaml and the classic behaviour (nothing extra mounted).
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
YAML_FILE=""       # Sim-design YAML; stays empty until set by -SimDesignYaml or auto-detected
SCENARIOS_DIR=""  # Scenarios directory to copy into container
RUN_MODE=false     # -Run: auto-run the single driver R script in the working folder, then exit
CURRENT_USER=$(whoami)
# Get current user's UID and GID for running containers as non-root
USER_ID=$(id -u)
GROUP_ID=$(id -g)
USER_NAME=$(whoami)
GROUP_NAME=$(id -gn)
# Supplementary host GIDs (space-separated). The entrypoint grants these inside
# the container so group-based access to shared bind mounts (e.g. the team's
# shared synthpop cache) works — `gosu uid:gid` alone strips supplementary
# groups. Older images ignore the variable; older scripts omit it harmlessly.
ADDITIONAL_GIDS=$(id -G)
# The same GIDs as --group-add flags for the HELPER containers (volume
# populate / rsync sync-back) that run with a bare `--user uid:gid` and no
# entrypoint — Docker grants supplementary groups to those directly.
# GID 0 (root group) is never granted; non-numeric tokens are dropped.
GROUP_ADD_FLAGS=""
for gid in $ADDITIONAL_GIDS; do
  case "$gid" in ''|*[!0-9]*) continue ;; esac
  [ "$gid" = "$GROUP_ID" ] && continue
  [ "$gid" -eq 0 ] && continue
  GROUP_ADD_FLAGS="$GROUP_ADD_FLAGS --group-add $gid"
done
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
    -Run)
      RUN_MODE=true
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

# -----------------------------------------------------------------------------
# Resolve the simulation-design YAML and, in "standalone" mode, the working
# folder to mount.
#
# Standalone mode = this launcher was downloaded into a folder that also holds
# the user's own sim_design*.yaml (plus scenarios / a driver script). When no
# -SimDesignYaml is given and exactly one YAML sits next to the script, it is
# used automatically and its files are merged (symlinked) into the model folder
# /IMPACTncd_England — so the user's config and scripts sit alongside the model
# and are referenced with no subfolder prefix, without passing any paths.
#
# If no YAML sits next to the script we fall back to the repo default
# (inputs/sim_design.yaml) and the classic behaviour, with nothing extra
# mounted. Passing -SimDesignYaml explicitly also mounts its containing folder,
# unless that folder is inside the repo tree (the developer workflow).
# -----------------------------------------------------------------------------
WORKDIR_MOUNT=""   # host folder merged into /IMPACTncd_England ("" = none)

if [ -n "$YAML_FILE" ]; then
  # Explicit -SimDesignYaml (already realpath'd and existence-checked above).
  YAML_DIR="$(cd "$(dirname "$YAML_FILE")" && pwd)"
  case "$YAML_DIR/" in
    "$PROJECT_ROOT"/*) : ;;            # inside the repo tree -> classic, no work mount
    *) WORKDIR_MOUNT="$YAML_DIR" ;;    # a user folder outside the repo -> mount it
  esac
else
  # No -SimDesignYaml: look for a single YAML next to this script. Prefer the
  # sim_design* naming convention; fall back to any single generic *.yaml/*.yml.
  shopt -s nullglob
  yaml_candidates=( "$SCRIPT_DIR"/sim_design*.yaml "$SCRIPT_DIR"/sim_design*.yml )
  if [ ${#yaml_candidates[@]} -eq 0 ]; then
    yaml_candidates=( "$SCRIPT_DIR"/*.yaml "$SCRIPT_DIR"/*.yml )
  fi
  shopt -u nullglob

  if [ ${#yaml_candidates[@]} -eq 1 ]; then
    YAML_FILE="$(realpath "${yaml_candidates[0]}")"
    WORKDIR_MOUNT="$SCRIPT_DIR"
    echo "Auto-detected sim-design YAML next to the script: $(basename "$YAML_FILE")"
  elif [ ${#yaml_candidates[@]} -gt 1 ]; then
    echo "Error: multiple YAML files sit next to the script; cannot choose one automatically:"
    for y in "${yaml_candidates[@]}"; do echo "  - $(basename "$y")"; done
    echo "Re-run with -SimDesignYaml <file> to pick one."
    exit 1
  else
    # Nothing next to the script -> classic repo default.
    YAML_FILE="$PROJECT_ROOT/inputs/sim_design.yaml"
  fi
fi

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

# Set simulation design file and extract output directories from YAML.
# `tr -d '\r'` strips any trailing carriage return so a YAML saved with
# Windows/CRLF line endings can't smuggle a \r into the mount path — that
# makes `docker run --mount` fail with "value should not have whitespace".
SIM_DESIGN_FILE="$YAML_FILE"
OUTPUT_DIR_RAW=$(grep '^output_dir:' "$SIM_DESIGN_FILE" | sed -E 's/output_dir:[[:space:]]*([^#]*).*/\1/' | tr -d '\r' | xargs)
SYNTHPOP_DIR_RAW=$(grep '^synthpop_dir:' "$SIM_DESIGN_FILE" | sed -E 's/synthpop_dir:[[:space:]]*([^#]*).*/\1/' | tr -d '\r' | xargs)

# Base for resolving RELATIVE output_dir/synthpop_dir from the YAML. In
# standalone mode a relative path (e.g. ./outputs) is relative to the YAML's
# own folder — i.e. inside the mounted working folder; in classic mode it is
# relative to the repo root, preserving the previous behaviour.
if [ -n "$WORKDIR_MOUNT" ]; then
  REL_BASE="$(cd "$(dirname "$YAML_FILE")" && pwd)"
else
  REL_BASE="$PROJECT_ROOT"
fi

# Resolve paths relative to REL_BASE if they are not absolute
if [[ "$OUTPUT_DIR_RAW" != /* && "$OUTPUT_DIR_RAW" != ~* ]]; then
  OUTPUT_DIR_TEMP="$REL_BASE/$OUTPUT_DIR_RAW"
else
  OUTPUT_DIR_TEMP="$OUTPUT_DIR_RAW"
fi
if [[ "$SYNTHPOP_DIR_RAW" != /* && "$SYNTHPOP_DIR_RAW" != ~* ]]; then
  SYNTHPOP_DIR_TEMP="$REL_BASE/$SYNTHPOP_DIR_RAW"
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

# -----------------------------------------------------------------------------
# Assemble the pieces shared by both mount strategies (volumes and bind mounts)
# so the actual `docker run` is written once each, instead of being duplicated
# for every scenarios / -Run combination.
# -----------------------------------------------------------------------------
RUN_ENV_ARGS=(
  -e USER_ID="${USER_ID}"
  -e GROUP_ID="${GROUP_ID}"
  -e USER_NAME="${USER_NAME}"
  -e GROUP_NAME="${GROUP_NAME}"
  -e ADDITIONAL_GIDS="${ADDITIONAL_GIDS}"
)

# Optional bind mounts: a user scenarios dir and/or the standalone working
# folder. In standalone mode the working folder is bind-mounted at a private
# STAGING path, and MERGE_SCRIPT (below) symlinks its files into the model root
# /IMPACTncd_England at container start — so they sit alongside global.R and the
# user references them with NO subfolder prefix. SIM_DESIGN_YAML points at that
# post-symlink location so a driver can do Simulation$new(Sys.getenv("SIM_DESIGN_YAML")).
WORK_STAGE="/mnt/impactncd_project"   # internal staging mount (users never see it)
EXTRA_MOUNTS=()
if [[ -n "$SCENARIOS_DIR" ]]; then
  EXTRA_MOUNTS+=( --mount "type=bind,source=$SCENARIOS_DIR,target=/IMPACTncd_England/scenarios" )
fi
if [[ -n "$WORKDIR_MOUNT" ]]; then
  EXTRA_MOUNTS+=( --mount "type=bind,source=$WORKDIR_MOUNT,target=$WORK_STAGE" )
  RUN_ENV_ARGS+=( -e SIM_DESIGN_YAML="/IMPACTncd_England/$(basename "$YAML_FILE")" )
fi

# Wrapper run at container start in standalone mode: symlink each file from the
# staging mount into the model root, skipping the launcher itself and any name
# that ALREADY exists there — so a user file never shadows a model file (global.R,
# Rpackage, …). The skip check runs against the real in-container tree, so there
# is no hardcoded reserved list to keep in sync. Then cd to the model root so the
# relative paths in global.R and in the user's driver resolve.
MERGE_SCRIPT='for f in '"$WORK_STAGE"'/*; do
  [ -e "$f" ] || continue
  n=$(basename "$f")
  case "$n" in setup_*_docker_env.sh|setup_*_docker_env.ps1) continue ;; esac
  if [ -e "/IMPACTncd_England/$n" ]; then
    echo "NOTE: \"$n\" already exists in the model folder; keeping the model version (rename your copy if you need it)." >&2
    continue
  fi
  ln -sfn "$f" "/IMPACTncd_England/$n"
done
cd /IMPACTncd_England'

# -Run: locate the single driver R script in the working folder to execute
# non-interactively (requires standalone mode, i.e. a mounted working folder).
DRIVER_NAME=""
if [ "$RUN_MODE" = true ]; then
  if [ -z "$WORKDIR_MOUNT" ]; then
    echo "Error: -Run needs a working folder (standalone mode). Put your sim_design*.yaml"
    echo "       and a simulate*.R driver next to this script, then re-run with -Run."
    exit 1
  fi
  shopt -s nullglob
  drivers=( "$WORKDIR_MOUNT"/simulate*.R "$WORKDIR_MOUNT"/simulate*.r )
  if [ ${#drivers[@]} -eq 0 ]; then
    drivers=( "$WORKDIR_MOUNT"/*.R "$WORKDIR_MOUNT"/*.r )
  fi
  shopt -u nullglob
  if [ ${#drivers[@]} -eq 1 ]; then
    DRIVER_NAME="$(basename "${drivers[0]}")"
  elif [ ${#drivers[@]} -gt 1 ]; then
    echo "Error: -Run found multiple R scripts in the working folder; cannot choose one:"
    for d in "${drivers[@]}"; do echo "  - $(basename "$d")"; done
    echo "Keep a single simulate*.R driver, or run interactively (omit -Run)."
    exit 1
  else
    echo "Error: -Run found no R script in the working folder to execute."
    exit 1
  fi
fi

# Container command + TTY flags. In standalone mode the real command (bash, or
# Rscript for -Run) is wrapped with MERGE_SCRIPT, which symlinks the user's files
# into the model root first. The whole in-container script is base64-encoded and
# passed as a single blob, then decoded to a file and exec'd inside. Encoding is
# essential for cross-shell safety: PowerShell 5.1 corrupts `docker run`
# arguments that contain embedded quotes/newlines, which mangled the raw script.
# A base64 blob has none of those characters, so both bash and PowerShell pass it
# intact. Running the decoded FILE (not a pipe) lets `exec bash` keep the TTY.
LAUNCH_DECODE='echo $0 | base64 -d > /tmp/impactncd_launch.sh; exec bash /tmp/impactncd_launch.sh'
if [ -n "$WORKDIR_MOUNT" ]; then
  if [ "$RUN_MODE" = true ]; then
    TTY_FLAGS=()                                 # non-interactive (batch)
    FULL_SCRIPT="$MERGE_SCRIPT"$'\n'"exec Rscript '$DRIVER_NAME'"
    echo "Auto-run mode: will run $DRIVER_NAME via Rscript (merged into the model folder), then exit."
  else
    TTY_FLAGS=( -it )
    FULL_SCRIPT="$MERGE_SCRIPT"$'\n'"exec bash"
  fi
  LAUNCH_B64=$(printf '%s' "$FULL_SCRIPT" | base64 | tr -d '\n')
  CONTAINER_CMD=( bash -c "$LAUNCH_DECODE" "$LAUNCH_B64" )
else
  TTY_FLAGS=( -it )                              # classic mode: nothing to merge
  CONTAINER_CMD=( bash )
fi

# Friendly hint for interactive standalone runs.
if [ "$RUN_MODE" != true ] && [ -n "$WORKDIR_MOUNT" ]; then
  echo "-----------------------------------------------------------------------"
  echo "Your files are available directly in the model folder (no subfolder to"
  echo "remember). Inside the container, run your simulation, e.g.:"
  echo "    R"
  echo "    > source(\"global.R\")"
  echo "    > IMPACTncd <- Simulation\$new(\"$(basename "$YAML_FILE")\")"
  echo "    > IMPACTncd\$run(1:10, multicore = TRUE, \"sc0\")\$export_summaries(multicore = TRUE)"
  echo "  ...or source your driver:  source(\"<your-driver>.R\")"
  echo "  (outputs are written to /outputs -> your host output_dir)"
  echo "-----------------------------------------------------------------------"
fi

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
    --user "${USER_ID}:${GROUP_ID}" $GROUP_ADD_FLAGS \
    -v "$OUTPUT_DIR":/source \
    -v "$VOLUME_OUTPUT_NAME":/volume \
    alpine sh -c "cp -r /source/. /volume/ 2>/dev/null || cp -a /source/. /volume/ 2>/dev/null || true"
  docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" $GROUP_ADD_FLAGS \
    -v "$SYNTHPOP_DIR":/source \
    -v "$VOLUME_SYNTHPOP_NAME":/volume \
    alpine sh -c "cp -r /source/. /volume/ 2>/dev/null || cp -a /source/. /volume/ 2>/dev/null || true"

  # Run the main container using the pre-built image
  echo "Running the main container using Docker volumes..."
  docker run "${TTY_FLAGS[@]}" --rm \
    "${RUN_ENV_ARGS[@]}" \
    --mount type=volume,source="$VOLUME_OUTPUT_NAME",target=/outputs \
    --mount type=volume,source="$VOLUME_SYNTHPOP_NAME",target=/synthpop \
    "${EXTRA_MOUNTS[@]}" \
    --workdir /IMPACTncd_England \
    "$IMAGE_NAME" \
    "${CONTAINER_CMD[@]}"

  # After the container exits:
  # - Synchronize the volumes back to the local directories using rsync (checksum mode).
  echo "Container exited. Syncing volumes back to local directories using rsync (checksum mode)..."
  docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" $GROUP_ADD_FLAGS \
    -v "$VOLUME_OUTPUT_NAME":/volume \
    -v "$OUTPUT_DIR":/backup \
    rsync-alpine rsync -avc --no-owner --no-group --no-times /volume/ /backup/
  docker run --rm \
    --user "${USER_ID}:${GROUP_ID}" $GROUP_ADD_FLAGS \
    -v "$VOLUME_SYNTHPOP_NAME":/volume \
    -v "$SYNTHPOP_DIR":/backup \
    rsync-alpine rsync -avc --no-owner --no-group --no-times /volume/ /backup/

  # Clean up all the Docker volumes used for the simulation.
  echo "Cleaning up Docker volumes..."
  docker volume rm "$VOLUME_OUTPUT_NAME"
  docker volume rm "$VOLUME_SYNTHPOP_NAME"
else
  echo "Using direct bind mounts for outputs and synthpop..."
  docker run "${TTY_FLAGS[@]}" --rm \
    "${RUN_ENV_ARGS[@]}" \
    --mount type=bind,source="$OUTPUT_DIR",target=/outputs \
    --mount type=bind,source="$SYNTHPOP_DIR",target=/synthpop \
    "${EXTRA_MOUNTS[@]}" \
    --workdir /IMPACTncd_England \
    "$IMAGE_NAME" \
    "${CONTAINER_CMD[@]}"
fi