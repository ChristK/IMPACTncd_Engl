#!/bin/bash
# -----------------------------------------------------------------------------
# build_push_data.sh
#
# Build the DATA image (layer 2 of 3) LOCALLY and push it to Docker Hub.
#
# WHY LOCAL (not GitHub CI):
#   The data image is ~33GB unpacked (~13GB Zenodo download). Building it on
#   GitHub-hosted runners is slow and disk-tight. The data also changes rarely.
#   So it is built here, on a machine with disk + bandwidth, and the downstream
#   images (model, etc.) build FROM the pushed image — they pull the data layers
#   from Docker Hub and never re-download from Zenodo.
#
# BRANCH BEHAVIOR (the script keys off the current git branch):
#   on main            -> maintains the CANONICAL images:
#       <user>/data.impactncdengl:<zenodo_version>   immutable, the data's identity
#       <user>/data.impactncdengl:latest             moving alias (downstream default)
#   on any other branch -> produces a BRANCH-scoped image instead:
#       <user>/data.impactncdengl:<branch>           for branches whose prereq/data
#                                                    inputs diverge from main
#     It builds FROM the branch's CI-built prerequisite
#     (<user>/prerequisite.impactncdengl:<branch>) and NEVER touches :latest or
#     the version tag — those belong to main. The model CI resolves
#     data:<branch> first and falls back to :latest, so branches without a
#     branch data image keep using main's environment automatically.
#   Both paths leave a local data.impactncdengl:local alias for local model builds.
#
# PREREQUISITE base (main builds): prerequisite.impactncdengl:local — pull the
# CI-built one (fast) rather than rebuilding:
#   docker pull <user>/prerequisite.impactncdengl:main
#   docker tag  <user>/prerequisite.impactncdengl:main prerequisite.impactncdengl:local
#
# USAGE (run from anywhere; the checked-out branch decides the behavior):
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

# Branch awareness: which branch is checked out decides what this run produces
# (see BRANCH BEHAVIOR above). The tag is sanitized the same way the CI
# workflows sanitize github.ref_name, so the tags line up.
BRANCH="$(git -C "$REPO_ROOT" rev-parse --abbrev-ref HEAD 2>/dev/null || true)"
if [[ -z "$BRANCH" || "$BRANCH" == "HEAD" ]]; then
  # rev-parse prints the literal "HEAD" on a detached checkout — pushing a
  # data image tagged :HEAD (or misclassifying the run) helps nobody.
  echo "ERROR: not on a branch (detached HEAD or not a git checkout)." >&2
  echo "       Check out the branch this data image is for, then re-run:" >&2
  echo "       git checkout main   (or your feature branch)" >&2
  exit 1
fi
SAFE_TAG="$(printf '%s' "$BRANCH" | sed -e 's|[^A-Za-z0-9._-]|-|g' -e 's|^[.-]|_|')"
SAFE_TAG="${SAFE_TAG:0:128}"
ON_MAIN=true
[[ "$BRANCH" != "main" ]] && ON_MAIN=false
if ! $ON_MAIN; then
  log "Branch build: '$BRANCH' -> will push ${DH_USER}/${IMAGE}:${SAFE_TAG} (NOT :latest)."
fi

# 0) Warn if the local clone is not at the tip of its remote branch. The
#    prerequisite builds from the live docker_setup dir and the data image from
#    `git archive HEAD`, so a stale checkout would bake OLD content into the
#    pushed image — while still turning the CI freshness guard green (it only
#    compares timestamps). Non-fatal: offline/detached use is allowed, but the
#    operator must know what they are shipping.
REMOTE_TIP="$(git -C "$REPO_ROOT" ls-remote origin "refs/heads/$BRANCH" 2>/dev/null | cut -f1 || true)"
LOCAL_HEAD="$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || true)"
if [[ -n "$REMOTE_TIP" && -n "$LOCAL_HEAD" && "$LOCAL_HEAD" != "$REMOTE_TIP" ]]; then
  log "WARNING: HEAD (${LOCAL_HEAD:0:8}) is not the tip of origin/$BRANCH (${REMOTE_TIP:0:8})."
  log "         The data image will bake THIS checkout's content. If that is not"
  log "         intended, Ctrl-C now and run: git pull --ff-only"
fi

# 0b) Resolve/validate the prerequisite base.
# Make any ambient base-image override VISIBLE: .env is sourced with `set -a`,
# so a stray PREREQ_IMAGE line there (or in the calling shell) silently reaches
# docker_build_push.sh as a --build-arg and changes what the data image builds
# FROM. Legitimate for experiments — but never silent.
if [[ -n "${PREREQ_IMAGE:-}" ]]; then
  log "NOTE: PREREQ_IMAGE override is set in the environment: ${PREREQ_IMAGE}"
  log "      The data image will build FROM it (not the default base). Unset it if unintended."
fi
if $ON_MAIN; then
  # main: the base is prerequisite.impactncdengl:local. Warn if it is missing
  # or older than the CI-built one — a stale :local would silently bake an OLD
  # prerequisite into the data image while the freshness guard (timestamps
  # only) still turns green. Pull+retag is the fix: faster than a rebuild and
  # the exact image CI validated. Non-fatal: offline use and intentional
  # custom local prereqs are allowed.
  PREREQ_LOCAL_CREATED="$(docker image inspect -f '{{.Created}}' prerequisite.impactncdengl:local 2>/dev/null || true)"
  PREREQ_HUB_JSON="$(curl -fsSL --max-time 15 \
    "https://hub.docker.com/v2/repositories/${DH_USER}/prerequisite.impactncdengl/tags/main" 2>/dev/null || true)"
  if [[ -z "$PREREQ_LOCAL_CREATED" ]]; then
    log "WARNING: prerequisite.impactncdengl:local not found. Pull the CI-built one first:"
    echo "      docker pull ${DH_USER}/prerequisite.impactncdengl:main"
    echo "      docker tag ${DH_USER}/prerequisite.impactncdengl:main prerequisite.impactncdengl:local"
  elif [[ -n "$PREREQ_HUB_JSON" ]]; then
    # Exact check first: a PULLED image carries a RepoDigest — if it matches any
    # digest the Hub lists for :main, the local base IS the CI-built image.
    # (.Created-vs-last_updated alone would false-alarm: build time always
    # precedes push time, even for the identical image.)
    HUB_DIGESTS="$(grep -oE '"digest":"sha256:[a-f0-9]+"' <<<"$PREREQ_HUB_JSON" | cut -d'"' -f4 | sort -u || true)"
    LOCAL_DIGESTS="$(docker image inspect -f '{{join .RepoDigests " "}}' prerequisite.impactncdengl:local 2>/dev/null || true)"
    DIGEST_MATCH=false
    for d in $HUB_DIGESTS; do
      [[ "$LOCAL_DIGESTS" == *"$d"* ]] && DIGEST_MATCH=true && break
    done
    if [[ "$DIGEST_MATCH" == false ]]; then
      # Fallback heuristic for locally-BUILT prereqs (no matching RepoDigest by
      # construction): only warn when the local image predates the Hub push.
      PREREQ_HUB_PUSHED="$(grep -oE '"last_updated":"[^"]+"' <<<"$PREREQ_HUB_JSON" | head -n1 | cut -d'"' -f4 || true)"
      T_LOCAL=$(date -d "$PREREQ_LOCAL_CREATED" +%s 2>/dev/null || echo 0)
      T_HUB=$(date -d "${PREREQ_HUB_PUSHED:-}" +%s 2>/dev/null || echo 0)
      if [[ "$T_LOCAL" -gt 0 && "$T_HUB" -gt 0 && "$T_LOCAL" -lt "$T_HUB" ]]; then
        log "WARNING: local prerequisite.impactncdengl:local (created $PREREQ_LOCAL_CREATED) looks"
        log "         OLDER than the CI-built ${DH_USER}/prerequisite.impactncdengl:main (pushed $PREREQ_HUB_PUSHED)."
        log "         The data image would bake the OLD prerequisite. If that is not intended,"
        log "         Ctrl-C now and refresh it with:"
        echo "      docker pull ${DH_USER}/prerequisite.impactncdengl:main"
        echo "      docker tag ${DH_USER}/prerequisite.impactncdengl:main prerequisite.impactncdengl:local"
      fi
    fi
  fi
else
  # Branch: build FROM the branch's CI-built prerequisite on Docker Hub (the
  # prerequisite workflow builds it on every push to the branch, including a
  # new-branch bootstrap). PREREQ_IMAGE is exported for docker_build_push.sh,
  # which forwards it as a --build-arg (Dockerfile: ARG PREREQ_IMAGE / FROM).
  BRANCH_PREREQ="${DH_USER}/prerequisite.impactncdengl:${SAFE_TAG}"
  if docker manifest inspect "$BRANCH_PREREQ" >/dev/null 2>&1; then
    export PREREQ_IMAGE="$BRANCH_PREREQ"
    log "Data image will build FROM the branch prerequisite: ${BRANCH_PREREQ}"
  else
    log "WARNING: no CI-built prerequisite for branch '$BRANCH' on Docker Hub"
    log "         (${BRANCH_PREREQ} not found). Falling back to the LOCAL"
    log "         prerequisite.impactncdengl:local — make sure it matches this"
    log "         branch's prerequisite inputs (push the branch so CI builds it)."
  fi
fi

# 1) Resolve the Zenodo data version (drives the immutable tag on main; logged
#    for provenance on branch builds).
#    Every stage is failure-guarded: under `set -e`, an unguarded pipeline
#    (Rscript stop()s when Zenodo is unreachable; grep exits 1 on no match)
#    would kill the script BEFORE the error handler below, silently. stderr is
#    captured, not discarded, so a failure shows WHY.
log "Resolving Zenodo data version for concept DOI ${CONCEPT_DOI} ..."
VERSION_RAW="$(cd "$REPO_ROOT" && ZENODO_CONCEPT_DOI="$CONCEPT_DOI" \
  Rscript docker_setup/print_zenodo_version.R 2>&1 || true)"
VERSION="$(grep -oE 'ZENODO_VERSION=[^[:space:]]+' <<<"$VERSION_RAW" | head -n1 | cut -d= -f2 || true)"
if [[ -z "${VERSION}" ]]; then
  if $ON_MAIN; then
    # main: the version IS the immutable tag — cannot proceed without it.
    echo "ERROR: could not resolve the Zenodo data version for ${CONCEPT_DOI}." >&2
    echo "       Check network access and the concept DOI. Rscript output (tail):" >&2
    tail -n 5 <<<"$VERSION_RAW" | sed 's/^/       /' >&2
    exit 1
  else
    # branch: the tag is the branch name; the version is provenance-only.
    # A host-side R problem must not kill a build whose in-container download
    # resolves the record independently anyway.
    log "WARNING: could not resolve the Zenodo data version (host R issue or Zenodo unreachable)."
    log "         Continuing — the branch tag is '${SAFE_TAG}'; version is provenance-only here."
    VERSION="unknown"
  fi
fi
log "Zenodo data version: ${VERSION}"

# 2) Build + push the data image. docker_build_push.sh handles the clean
#    git-archive build context (of the CURRENT branch's HEAD), Docker Hub
#    login, and the build-args (ZENODO_CONCEPT_DOI / DOWNLOAD_DATA /
#    PREREQ_IMAGE) from the environment. It must run from docker_setup.
#    main   -> immutable Zenodo-version tag;  branch -> the branch tag.
if $ON_MAIN; then
  TAG_TO_BUILD="${VERSION}"
else
  TAG_TO_BUILD="${SAFE_TAG}"
fi
log "Building + pushing ${DH_USER}/${IMAGE}:${TAG_TO_BUILD} ..."
( cd "$SCRIPT_DIR" && ./docker_build_push.sh Dockerfile.data.IMPACTncdENGL --image-tag "${TAG_TO_BUILD}" --push )

# 3) main only: add + push the moving :latest alias. Both paths keep a :local
#    alias for local model builds. A branch build must NEVER move :latest —
#    that would swap the whole team's environment for branch content.
if $ON_MAIN; then
  log "Tagging + pushing ${DH_USER}/${IMAGE}:latest ..."
  docker tag "${DH_USER}/${IMAGE}:${TAG_TO_BUILD}" "${DH_USER}/${IMAGE}:latest"
  docker push "${DH_USER}/${IMAGE}:latest"
else
  log "Branch build: :latest deliberately NOT updated (belongs to main)."
fi
docker tag "${DH_USER}/${IMAGE}:${TAG_TO_BUILD}" "${IMAGE}:local"

# 4) Close the rebuild chain: trigger the model image CI build on this branch
#    so the fresh data image reaches ${CI_DATA_NS}/impactncdengl:<branch>
#    without waiting for the next code push. Uses the GitHub CLI if available +
#    authenticated; otherwise prints the manual command. Failure here is
#    non-fatal — the data push above succeeded.
#    NOTE: CI consumes images from the ${CI_DATA_NS} namespace (the repo's
#    DOCKERHUB_USERNAME secret). If this script pushed to a DIFFERENT
#    namespace, a rebuild now would NOT contain the data just pushed — skip
#    the trigger and say so, instead of giving false assurance.
CI_DATA_NS="chriskypri"   # must match secrets.DOCKERHUB_USERNAME of ChristK/IMPACTncd_Engl
if [[ "${DH_USER}" != "${CI_DATA_NS}" ]]; then
  log "WARNING: data image was pushed to ${DH_USER}/, but CI on ChristK/IMPACTncd_Engl"
  log "         builds FROM the ${CI_DATA_NS}/ namespace. Skipping the automatic"
  log "         model-rebuild trigger — a rebuild now would NOT contain the data you just"
  log "         pushed. Either push as ${CI_DATA_NS}, or update the DOCKERHUB_USERNAME"
  log "         secret / DATA_IMAGE build-arg in the repo."
elif command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
  log "Triggering the model image rebuild on '$BRANCH' (GitHub Actions) ..."
  if gh workflow run build_push_impactncdengl.yml --ref "$BRANCH" --repo ChristK/IMPACTncd_Engl; then
    log "Model rebuild triggered. Watch it with:"
    echo "      gh run list --workflow build_push_impactncdengl.yml --repo ChristK/IMPACTncd_Engl"
    if ! $ON_MAIN; then
      log "NOTE: the dispatched run executes the BRANCH's copy of the workflow. If"
      log "      '$BRANCH' predates the branch-aware data resolution, that run still"
      log "      builds FROM data:latest — merge main into the branch to pick it up."
    fi
  else
    log "WARNING: failed to trigger the model workflow — trigger it manually:"
    echo "      gh workflow run build_push_impactncdengl.yml --ref '$BRANCH' --repo ChristK/IMPACTncd_Engl"
    echo "      (the branch's workflow file needs the workflow_dispatch trigger AND the"
    echo "       branch-aware data resolution — merge main into the branch if missing)"
  fi
else
  log "NOTE: GitHub CLI unavailable or unauthenticated — trigger the model rebuild manually:"
  echo "      gh workflow run build_push_impactncdengl.yml --ref '$BRANCH' --repo ChristK/IMPACTncd_Engl"
  echo "      (or use the 'Run workflow' button on the repo's Actions tab)"
fi

log "Done."
if $ON_MAIN; then
  echo "  pushed  ${DH_USER}/${IMAGE}:${VERSION}   (immutable Zenodo version)"
  echo "  pushed  ${DH_USER}/${IMAGE}:latest        (moving alias; downstream default)"
else
  echo "  pushed  ${DH_USER}/${IMAGE}:${SAFE_TAG}   (branch-scoped; consumed by the"
  echo "          branch's model CI via its data-resolution step; :latest untouched)"
  echo "  data version baked in: ${VERSION} (concept DOI ${CONCEPT_DOI})"
fi
echo "  local   ${IMAGE}:local                    (for default local model builds)"
