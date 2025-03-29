#!/bin/bash

# This script automates building and pushing a Docker image to Docker Hub.
#
# Features:
# - Sources environment variables from `.env` if present.
# - Prompts for Docker Hub username if not set.
# - Builds from Dockerfile.IMPACTncd, with build context one folder up.
# - Tags image as <DOCKERHUB_USERNAME>/impactncd_engl:latest
# - Pushes the image to Docker Hub with logging and error handling.
#
# Usage:
# 1. (Optional) Create a `.env` file in the same directory with the following variables:
#       DOCKERHUB_USERNAME=yourusername
#       DOCKERHUB_TOKEN=youraccesstoken
#
# 2. Run the script as follows:
#       ./build_and_push_IMPACTncd.sh [--push]
#
#    Use the '--push' flag to push the Docker image to Docker Hub. Without the flag, the script will only build the image.

set -euo pipefail

# Check if push to Docker should be executed
PUSH_IMAGE=false
if [[ "$#" -gt 0 && "$1" == "--push" ]]; then
  PUSH_IMAGE=true
  shift
fi

# Source .env file if it exists
if [[ -f ".env" ]]; then
  echo "Sourcing environment variables from .env"
  source .env
fi

# Prompt for Docker Hub username if not set
if [[ -z "${DOCKERHUB_USERNAME:-}" ]]; then
  read -p "Enter your Docker Hub username: " DOCKERHUB_USERNAME
fi

IMAGE_NAME="${DOCKERHUB_USERNAME}/impactncd_engl:latest"

# Timestamped logging function
log() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"
}

# Docker login
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

# Build the Docker image using the parent folder as build context
log "Building Docker image from Dockerfile.IMPACTncd..."
if docker build --no-cache -f Dockerfile.IMPACTncd -t "$IMAGE_NAME" ..; then
  log "Docker image built successfully."
else
  log "Docker image build failed."
  exit 1
fi

# Push the Docker image
if [ "$PUSH_IMAGE" = true ]; then
  log "Pushing Docker image to Docker Hub..."
  if docker push "$IMAGE_NAME"; then
    log "Docker image pushed successfully."
  else
    log "Docker push failed."
    exit 1
  fi
else
  log "Skipping Docker push as --push argument not provided."
fi