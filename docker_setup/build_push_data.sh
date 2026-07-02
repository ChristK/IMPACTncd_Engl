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

# 0) Warn if the local clone is not at the tip of origin/main. The prerequisite
#    builds from the live docker_setup dir and the data image from
#    `git archive HEAD`, so a stale or wrong-branch checkout would bake OLD
#    content into the pushed image — while still turning the CI freshness
#    guard green (it only compares timestamps). Non-fatal: offline/detached
#    use is allowed, but the operator must know what they are shipping.
REMOTE_MAIN="$(git -C "$REPO_ROOT" ls-remote origin refs/heads/main 2>/dev/null | cut -f1 || true)"
LOCAL_HEAD="$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || true)"
if [[ -n "$REMOTE_MAIN" && -n "$LOCAL_HEAD" && "$LOCAL_HEAD" != "$REMOTE_MAIN" ]]; then
  log "WARNING: HEAD (${LOCAL_HEAD:0:8}) is not the tip of origin/main (${REMOTE_MAIN:0:8})."
  log "         The data image will bake THIS checkout's content. If that is not"
  log "         intended, Ctrl-C now and run: git checkout main && git pull --ff-only"
fi

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

# 4) Close the rebuild chain: trigger the model image CI build on main so the
#    fresh data:latest (and the prerequisite baked beneath it) reaches
#    ${CI_DATA_NS}/impactncdengl:main without waiting for the next code push.
#    Uses the GitHub CLI if available + authenticated; otherwise prints the
#    manual command. Failure here is non-fatal — the data push above succeeded.
#    NOTE: CI consumes ${CI_DATA_NS}/data.impactncdengl:latest (the repo's
#    DOCKERHUB_USERNAME secret). If this script pushed to a DIFFERENT
#    namespace, a rebuild now would NOT contain the data just pushed — skip
#    the trigger and say so, instead of giving false assurance.
CI_DATA_NS="chriskypri"   # must match secrets.DOCKERHUB_USERNAME of ChristK/IMPACTncd_Engl
if [[ "${DH_USER}" != "${CI_DATA_NS}" ]]; then
  log "WARNING: data image was pushed to ${DH_USER}/, but CI on ChristK/IMPACTncd_Engl"
  log "         builds FROM ${CI_DATA_NS}/data.impactncdengl:latest. Skipping the automatic"
  log "         model-rebuild trigger — a rebuild now would NOT contain the data you just"
  log "         pushed. Either push as ${CI_DATA_NS}, or update the DOCKERHUB_USERNAME"
  log "         secret / DATA_IMAGE build-arg in the repo."
elif command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
  log "Triggering the model image rebuild on main (GitHub Actions) ..."
  if gh workflow run build_push_impactncdengl.yml --ref main --repo ChristK/IMPACTncd_Engl; then
    log "Model rebuild triggered. Watch it with:"
    echo "      gh run list --workflow build_push_impactncdengl.yml --repo ChristK/IMPACTncd_Engl"
  else
    log "WARNING: failed to trigger the model workflow — trigger it manually:"
    echo "      gh workflow run build_push_impactncdengl.yml --ref main --repo ChristK/IMPACTncd_Engl"
    echo "      (needs the workflow_dispatch trigger on main; present since 2026-07-02)"
  fi
else
  log "NOTE: GitHub CLI unavailable or unauthenticated — trigger the model rebuild manually:"
  echo "      gh workflow run build_push_impactncdengl.yml --ref main --repo ChristK/IMPACTncd_Engl"
  echo "      (or use the 'Run workflow' button on the repo's Actions tab)"
fi

log "Done."
echo "  pushed  ${DH_USER}/${IMAGE}:${VERSION}   (immutable Zenodo version)"
echo "  pushed  ${DH_USER}/${IMAGE}:latest        (moving alias; downstream default)"
echo "  local   ${IMAGE}:local                    (for default local model builds)"
