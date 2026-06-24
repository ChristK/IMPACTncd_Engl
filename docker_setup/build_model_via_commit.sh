#!/bin/bash
# -----------------------------------------------------------------------------
# build_model_via_commit.sh
#
# Build the model image (layer 3) on top of the data image WITHOUT buildkit,
# using `docker commit`. This mirrors Dockerfile.IMPACTncdENGL exactly (merge
# code over the data layer, build + install the package, snapshot, chmod) but
# avoids buildkit's image export.
#
# WHY THIS EXISTS:
#   On hosts that use the containerd image store with an aggressive garbage
#   collector, containerd can evict a base layer's *content* blob WHILE a build
#   is running. The image still runs (running uses the unpacked snapshot store),
#   but `docker build` of a child image fails at export with:
#       failed to extract layer sha256:...: blob sha256:... not found
#   `docker commit` never re-reads base content blobs (it references lower
#   layers by digest and only writes one new top layer), so it is immune to the
#   GC race and builds the model image reliably.
#
#   The proper long-term fix is on the host: tune the containerd GC (e.g.
#   io.containerd.gc.v1.scheduler thresholds) or ensure buildkit holds a lease
#   on the base image for the build's duration. Until then, use this script.
#
# PREREQUISITE: build the data image first (this works via buildkit):
#   ./docker_build_push.sh Dockerfile.data.IMPACTncdENGL
#
# USAGE:
#   ./build_model_via_commit.sh [DATA_IMAGE] [MODEL_IMAGE]
#   (defaults: data.impactncdengl:local -> impactncdengl:local)
# -----------------------------------------------------------------------------
set -euo pipefail

DATA_IMAGE="${1:-data.impactncdengl:local}"
MODEL_IMAGE="${2:-impactncdengl:local}"
CONTAINER="model_build_$$"
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

log(){ echo "$(date '+%H:%M:%S') $*"; }

if ! docker image inspect "$DATA_IMAGE" >/dev/null 2>&1; then
  echo "ERROR: data image '$DATA_IMAGE' not found."
  echo "Build it first:  ./docker_build_push.sh Dockerfile.data.IMPACTncdENGL"
  exit 1
fi

cleanup(){ docker rm -f "$CONTAINER" >/dev/null 2>&1 || true; }
trap cleanup EXIT

log "Starting build container from $DATA_IMAGE ..."
docker run -d --name "$CONTAINER" --entrypoint sleep "$DATA_IMAGE" infinity >/dev/null

log "Merging code (git-tracked HEAD) over the data layer ..."
git -C "$PROJECT_ROOT" archive HEAD | docker cp - "$CONTAINER":/IMPACTncd_England

log "Building the package inside the container (roxy + build + install) ..."
docker exec -i "$CONTAINER" bash -s <<'BUILD'
set -e
export PATH=/usr/local/lib/R/site-library/littler/bin:/usr/local/lib/R/site-library/littler/examples:$PATH
mkdir -p /outputs /synthpop
cd /IMPACTncd_England/Rpackage/IMPACTncd_England_model_pkg
roxy.r
build.r
PACKAGE=$(grep "^Package:" DESCRIPTION | awk '{print $2}')
VERSION=$(grep "^Version:" DESCRIPTION | awk '{print $2}')
install2.r "${PACKAGE}_${VERSION}.tar.gz"
cd /IMPACTncd_England/
Rscript -e 'snapshot <- utils::fileSnapshot("Rpackage/IMPACTncd_England_model_pkg/", timestamp = NULL, md5sum = TRUE, recursive = TRUE); saveRDS(snapshot, "Rpackage/.IMPACTncd_England_model_pkg_snapshot.rds")'
# Make code + runtime dirs writable (data is already a+rw from the data layer).
chmod -R a+rw /IMPACTncd_England/Rpackage /outputs /synthpop
chmod a+rw /IMPACTncd_England /IMPACTncd_England/* 2>/dev/null || true
Rscript -e 'stopifnot(requireNamespace("IMPACTncdEngland", quietly = TRUE)); cat("package installs OK\n")'
BUILD

log "Committing container -> $MODEL_IMAGE (restoring the prerequisite entrypoint/PATH) ..."
docker commit \
  --change 'ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]' \
  --change 'CMD ["/bin/bash"]' \
  --change 'WORKDIR /IMPACTncd_England' \
  --change 'ENV PATH=/usr/local/lib/R/site-library/littler/bin:/usr/local/lib/R/site-library/littler/examples:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin' \
  "$CONTAINER" "$MODEL_IMAGE" >/dev/null

log "Done. Built model image:"
docker images "$MODEL_IMAGE" --format "  {{.Repository}}:{{.Tag}}  {{.ID}}  {{.Size}}"
