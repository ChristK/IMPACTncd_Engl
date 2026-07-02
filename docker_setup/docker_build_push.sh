#!/bin/bash

# This script automates building and pushing a Docker image to Docker Hub.
#
# Features:
# - Sources config from `.env` (DOCKERHUB_USERNAME, etc.)
# - Constructs the image name as: <image-name>:<image-tag>.
# - If pushed to Docker Hub, tagged as <DOCKERHUB_USERNAME>/<image-name>:<image-tag>.
# - Image name and tag can be specified via arguments `--image-name` and `--image-tag`
# - Uses token-based login if available; otherwise prompts for manual login.
# - For Dockerfile.IMPACTncdENGL: builds a clean context from git-tracked files
#   only. The model data is NOT bundled — it is downloaded from Zenodo during
#   the build (see Dockerfile.IMPACTncdENGL / download_zenodo_data.R).
#
# Config:
#   Non-secret config is read from `.env`:
#     DOCKERHUB_USERNAME
#
# Usage:
#   ./docker_build_push.sh <dockerfile> [--push]
#   ./docker_build_push.sh <dockerfile> [--image-name <name>] [--image-tag <tag>] [--push]
#
# Examples:
#   ./docker_build_push.sh Dockerfile.IMPACTncdENGL                 # build only
#   ./docker_build_push.sh Dockerfile.IMPACTncdENGL --push          # build + push
#   ./docker_build_push.sh Dockerfile.prerequisite.IMPACTncdENGL    # prerequisite image

set -euo pipefail

# Parse arguments for image name, tag, and dockerfile
IMAGE_TAG="local"
PUSH_IMAGE=false
DOCKERFILE=""

while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --image-name)
      IMAGE_NAME="$2"
      shift 2
      ;;
    --image-tag)
      IMAGE_TAG="$2"
      shift 2
      ;;
    --push)
      PUSH_IMAGE=true
      shift
      ;;
    *)
      DOCKERFILE="$1"
      shift
      ;;
  esac
done

# Ensure DOCKERFILE is set before IMAGE_NAME is derived
if [[ -z "$DOCKERFILE" ]]; then
  echo "Error: Dockerfile argument is required."
  exit 1
fi

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

# Source .env file for non-secret config (DOIs, usernames, flags).
# The .env file is gitignored. Permissions should be restricted: chmod 600 .env
# CALLER-ENV PRECEDENCE for the base-image overrides: a caller that explicitly
# sets PREREQ_IMAGE/DATA_IMAGE (build_push_data.sh exports the branch-specific
# prerequisite on non-main branches) must not have it silently clobbered by a
# stray line in .env — an explicit environment beats a config-file default.
_CALLER_PREREQ_IMAGE="${PREREQ_IMAGE:-}"
_CALLER_DATA_IMAGE="${DATA_IMAGE:-}"
if [[ -f ".env" ]]; then
  echo "Sourcing config from .env"
  set -a
  source .env
  set +a
fi
if [[ -n "$_CALLER_PREREQ_IMAGE" && "$_CALLER_PREREQ_IMAGE" != "${PREREQ_IMAGE:-}" ]]; then
  echo "NOTE: keeping caller-provided PREREQ_IMAGE=${_CALLER_PREREQ_IMAGE} (overrides .env)"
  PREREQ_IMAGE="$_CALLER_PREREQ_IMAGE"
fi
if [[ -n "$_CALLER_DATA_IMAGE" && "$_CALLER_DATA_IMAGE" != "${DATA_IMAGE:-}" ]]; then
  echo "NOTE: keeping caller-provided DATA_IMAGE=${_CALLER_DATA_IMAGE} (overrides .env)"
  DATA_IMAGE="$_CALLER_DATA_IMAGE"
fi


# Remove hardcoded default for IMAGE_NAME
# Respect user-defined name if provided, otherwise derive from Dockerfile
if [[ -z "${IMAGE_NAME:-}" ]]; then
  IMAGE_NAME="$(basename "$DOCKERFILE" | sed 's/^[Dd]ockerfile\.//' | tr '[:upper:]' '[:lower:]')"
fi

# Convert the image NAME to lowercase — Docker requires lowercase repository
# names. The TAG is deliberately case-PRESERVED: Docker tags are case-sensitive
# and may contain uppercase, and callers align tags with git branch names
# (e.g. data.impactncdengl:Report_1, prerequisite.impactncdengl:update_from_JPN
# — the CI workflows push case-preserved branch tags). Lowercasing the tag here
# would push data:report_1 while CI resolves data:Report_1 (silent fallback to
# :latest) and then crash the caller's retag step after the full build.
if [[ "$IMAGE_NAME" != "${IMAGE_NAME,,}" ]]; then
    echo "Image name '$IMAGE_NAME' converted to lowercase."
    IMAGE_NAME="${IMAGE_NAME,,}"
fi

# Construct the image name for building
BUILD_IMAGE_NAME="${IMAGE_NAME}:${IMAGE_TAG}"

# Function for timestamped log messages
log() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Skip credential checks if push argument is false
if [ "$PUSH_IMAGE" = true ]; then
  if [[ -z "${DOCKERHUB_USERNAME:-}" ]]; then
    read -p "Enter your Docker Hub username: " DOCKERHUB_USERNAME
  fi

  # Construct the full image name for pushing
  FULL_IMAGE_NAME="${DOCKERHUB_USERNAME}/${IMAGE_NAME}:${IMAGE_TAG}"

  log "Logging into Docker Hub..."
  if [[ -z "${DOCKERHUB_USERNAME:-}" || -z "${DOCKERHUB_TOKEN:-}" ]]; then
    log "Environment variables DOCKERHUB_USERNAME and/or DOCKERHUB_TOKEN not set."
    log "Please log in manually:"
    docker login
  else
    if echo "$DOCKERHUB_TOKEN" | docker login -u "$DOCKERHUB_USERNAME" --password-stdin; then
      log "Docker login successful."
    else
      log "Docker login failed."
      exit 1
    fi
  fi
fi

# Build the Docker image
log "Building Docker image..."
# For the data and model images, create a clean build context from git-tracked
# files only (via `git archive`). The model DATA is NOT bundled — it is
# downloaded from Zenodo in the data image. Using only git-tracked files keeps
# sensitive files (.env, tokens) and other untracked content out of the image
# and keeps the build context small. The prerequisite image uses the current
# directory (docker_setup) as its build context.
BUILD_CONTEXT_DIR=""
DOCKERFILE_BASE="$(basename "$DOCKERFILE")"
if [[ "$DOCKERFILE_BASE" == "Dockerfile.IMPACTncdENGL" || "$DOCKERFILE_BASE" == "Dockerfile.data.IMPACTncdENGL" ]]; then
  BUILD_CONTEXT_DIR=$(mktemp -d)
  trap 'rm -rf "$BUILD_CONTEXT_DIR"' EXIT

  log "Creating build context from git-tracked files..."
  git -C ".." archive HEAD | tar -x -C "$BUILD_CONTEXT_DIR"

  BUILD_CONTEXT="$BUILD_CONTEXT_DIR"
  log "Build context ready at $BUILD_CONTEXT_DIR"
else
  BUILD_CONTEXT="."
  log "Using current directory (.) as build context for $(basename "$DOCKERFILE")"
fi

# Pass the Zenodo concept DOI (and download flag) from .env to the build, so the
# image downloads the matching data record. The Dockerfile has sensible
# defaults (the published record, DOWNLOAD_DATA=true) if these are unset.
BUILD_ARGS=()
[[ -n "${ZENODO_CONCEPT_DOI:-}" ]] && BUILD_ARGS+=(--build-arg "ZENODO_CONCEPT_DOI=${ZENODO_CONCEPT_DOI}")
[[ -n "${DOWNLOAD_DATA:-}" ]] && BUILD_ARGS+=(--build-arg "DOWNLOAD_DATA=${DOWNLOAD_DATA}")
# Base-image overrides (the Dockerfiles declare ARG PREREQ_IMAGE / DATA_IMAGE
# with :local defaults). Lets a caller point the build at e.g. a branch-tagged
# base from Docker Hub, as build_push_data.sh does for non-main branches:
#   PREREQ_IMAGE=chriskypri/prerequisite.impactncdengl:<branch>
[[ -n "${PREREQ_IMAGE:-}" ]] && BUILD_ARGS+=(--build-arg "PREREQ_IMAGE=${PREREQ_IMAGE}")
[[ -n "${DATA_IMAGE:-}" ]] && BUILD_ARGS+=(--build-arg "DATA_IMAGE=${DATA_IMAGE}")

# Optional build network. Some hosts use an internal DNS that the Docker bridge
# network cannot reach, which makes apt (prerequisite image) and the in-build
# Zenodo download (main image) fail to resolve names. Set
# DOCKER_BUILD_NETWORK=host (e.g. in .env) to build with host networking.
# Unset = Docker's default bridge network.
NETWORK_ARG=()
[[ -n "${DOCKER_BUILD_NETWORK:-}" ]] && NETWORK_ARG=(--network "${DOCKER_BUILD_NETWORK}")

if docker build --no-cache ${NETWORK_ARG[@]+"${NETWORK_ARG[@]}"} ${BUILD_ARGS[@]+"${BUILD_ARGS[@]}"} -f "$DOCKERFILE" -t "$BUILD_IMAGE_NAME" "$BUILD_CONTEXT"; then
  log "Docker image built successfully."
else
  log "Docker image build failed."
  exit 1
fi

# Push the Docker image to Docker Hub
if [ "$PUSH_IMAGE" = true ]; then
  log "Tagging Docker image for Docker Hub..."
  docker tag "$BUILD_IMAGE_NAME" "$FULL_IMAGE_NAME"

  log "Pushing Docker image to Docker Hub..."
  if docker push "$FULL_IMAGE_NAME"; then
    log "Docker image pushed successfully."
  else
    log "Docker push failed."
    exit 1
  fi
else
  log "Skipping Docker push as --push argument not provided."
fi

# Display the final Docker image name
echo "Final Docker image name: ${IMAGE_NAME}:${IMAGE_TAG}"