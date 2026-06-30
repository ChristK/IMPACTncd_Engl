#!/bin/bash
# -----------------------------------------------------------------------------
# build_push_data.sh
#
# Build the DATA image (layer 2 of 3) LOCALLY, tag it with the Zenodo data
# version (auto-derived) plus a moving :latest, and push both to Docker Hub.
#
# WHY LOCAL (not GitHub CI):
#   The data image is ~33GB unpacked (~13GB Zenodo download). Building it on
#   GitHub-hosted runners is slow and disk-tight. The data also changes rarely.
#   So it is built here, on a machine with disk + bandwidth, and the downstream
#   images (model, etc.) build FROM the pushed image — they pull the data layers
#   from Docker Hub and never re-download from Zenodo.
#
# TAGS PRODUCED (DOCKERHUB_USERNAME from docker_setup/.env, e.g. chriskypri):
#   <user>/data.impactncdengl:<zenodo_version>   immutable, the data's identity
#   <user>/data.impactncdengl:latest             moving alias (downstream default)
#   data.impactncdengl:local                     local alias for default local
#                                                 model builds (ARG DATA_IMAGE)
#
# PREREQUISITE: build the prerequisite image first (its :local tag is the base):
#   ./docker_build_push.sh Dockerfile.prerequisite.IMPACTncdENGL
#
# USAGE (run from anywhere):
#   docker_setup/build_push_data.sh
#   ZENODO_CONCEPT_DOI=10.5281/zenodo.XXXX docker_setup/build_push_data.sh
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # docker_setup
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
IMAGE="data.impactncdengl"

# Non-secret config (DOCKERHUB_USERNAME, ZENODO_CONCEPT_DOI) from .env.
if [[ -f "$SCRIPT_DIR/.env" ]]; then
  set -a; source "$SCRIPT_DIR/.env"; set +a
fi
CONCEPT_DOI="${ZENODO_CONCEPT_DOI:-10.5281/zenodo.20812409}"
DH_USER="${DOCKERHUB_USERNAME:?Set DOCKERHUB_USERNAME (in docker_setup/.env)}"

log() { echo "$(date '+%Y-%m-%d %H:%M:%S') - $1"; }

# 1) Resolve the Zenodo data version (drives the immutable tag).
log "Resolving Zenodo data version for concept DOI ${CONCEPT_DOI} ..."
VERSION="$(cd "$REPO_ROOT" && ZENODO_CONCEPT_DOI="$CONCEPT_DOI" \
  Rscript docker_setup/print_zenodo_version.R 2>/dev/null \
  | grep -oE 'ZENODO_VERSION=[^[:space:]]+' | head -n1 | cut -d= -f2)"
if [[ -z "${VERSION}" ]]; then
  echo "ERROR: could not resolve the Zenodo data version for ${CONCEPT_DOI}." >&2
  echo "       Check network access and the concept DOI." >&2
  exit 1
fi
log "Zenodo data version: ${VERSION}"

# 2) Build + push the version-tagged data image. docker_build_push.sh handles
#    the clean git-archive build context, Docker Hub login, and the build-args
#    (ZENODO_CONCEPT_DOI / DOWNLOAD_DATA) from .env. It must run from docker_setup.
log "Building + pushing ${DH_USER}/${IMAGE}:${VERSION} ..."
( cd "$SCRIPT_DIR" && ./docker_build_push.sh Dockerfile.data.IMPACTncdENGL --image-tag "${VERSION}" --push )

# 3) Add + push the moving :latest alias; keep a :local alias for local model builds.
log "Tagging + pushing ${DH_USER}/${IMAGE}:latest ..."
docker tag "${DH_USER}/${IMAGE}:${VERSION}" "${DH_USER}/${IMAGE}:latest"
docker push "${DH_USER}/${IMAGE}:latest"
docker tag "${DH_USER}/${IMAGE}:${VERSION}" "${IMAGE}:local"

log "Done."
echo "  pushed  ${DH_USER}/${IMAGE}:${VERSION}   (immutable Zenodo version)"
echo "  pushed  ${DH_USER}/${IMAGE}:latest        (moving alias; downstream default)"
echo "  local   ${IMAGE}:local                    (for default local model builds)"
