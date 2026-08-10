#!/usr/bin/env bash
# =============================================================================
# Regression test for .github/workflows/data_image_freshness_guard.yml
# =============================================================================
#
# Runs the workflow's OWN shell script (extracted from the YAML block scalar),
# so the logic under test cannot drift from a copy kept here. Only the network
# is stubbed: the two `hub_tag_epoch` calls are replaced by fixture epochs, and
# the substitution is verified to have taken effect — a silently un-stubbed
# script would reach out to Docker Hub and quietly test nothing.
#
# Covers the property added on 2026-08-10: a line-ending-only commit must not
# trip CHECK 1 for files whose EOLs the build cannot see (Dockerfiles, R
# sources), while staying STRICT for shell scripts and shell-parsed package
# lists, where a stray CR is a real, breaking change.
#
# Usage:  auxil/test_data_image_freshness_guard.sh          # offline test suite
#         auxil/test_data_image_freshness_guard.sh --live   # run the real guard
#
# --live runs the unstubbed guard against Docker Hub and this checkout, i.e.
# exactly what CI will say about HEAD. Run it before pushing to main to find out
# whether the data image needs a rebuild without waiting for the red X.
# Override the namespace with HUB_USER=... (default chriskypri).
#
# Exit:   0 = all cases pass, 1 = at least one failed.
# =============================================================================

set -uo pipefail

LIVE=0
[ "${1:-}" = "--live" ] && LIVE=1

REPO_ROOT="$(git -C "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" rev-parse --show-toplevel)"
WF="$REPO_ROOT/.github/workflows/data_image_freshness_guard.yml"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

PASS=0
FAIL=0

# --- extract the `run: |` block scalar -------------------------------------
# Generic enough to survive the step moving within the file: keyed on the
# `run: |` line's indent, ends at the first line indented no deeper.
extract_run_block() {
  awk '
    !seen && /^[[:space:]]*run:[[:space:]]*\|[[:space:]]*$/ {
      seen = 1; match($0, /^[[:space:]]*/); key = RLENGTH; next
    }
    seen {
      if ($0 ~ /^[[:space:]]*$/) { print ""; next }
      match($0, /^[[:space:]]*/)
      if (RLENGTH <= key) exit
      if (!body) body = RLENGTH
      print substr($0, body + 1)
    }
  ' "$1"
}

SCRIPT="$WORK/guard.sh"
extract_run_block "$WF" > "$SCRIPT.raw"
if [ ! -s "$SCRIPT.raw" ]; then
  echo "FATAL: could not extract the run: block from $WF" >&2
  exit 1
fi

# Cross-check the awk extractor against a real YAML parser (R's yaml is already
# a hard dependency of this project). This validates the workflow's syntax AND
# proves the block scalar we are testing is byte-identical to the one CI runs —
# an extractor that quietly dropped the tail would still "pass" every case.
if Rscript -e 'quit(status = !requireNamespace("yaml", quietly = TRUE))' 2>/dev/null; then
  cat > "$WORK/extract.R" <<'EOF'
a <- commandArgs(TRUE)
wf <- yaml::read_yaml(a[1])
runs <- Filter(Negate(is.null), lapply(wf$jobs$check_freshness$steps, `[[`, "run"))
if (length(runs) != 1L) stop("expected exactly one run: step, found ", length(runs))
cat(runs[[1]], file = a[2])
EOF
  if ! Rscript "$WORK/extract.R" "$WF" "$WORK/yaml.raw" 2>"$WORK/yaml.err"; then
    echo "FATAL: $WF is not valid YAML (or the job/step layout changed):" >&2
    cat "$WORK/yaml.err" >&2
    exit 1
  fi
  # `awk 1` normalises the trailing newline on both sides: R's yaml drops the
  # clip-chomped final newline of a `|` block ("one\ntwo", not "one\ntwo\n"),
  # which is a parser quirk and not a difference in the file.
  awk 1 "$WORK/yaml.raw" > "$WORK/yaml.norm"
  awk 1 "$SCRIPT.raw"    > "$WORK/awk.norm"
  if ! diff -u "$WORK/yaml.norm" "$WORK/awk.norm" > "$WORK/yaml.diff"; then
    echo "FATAL: awk extraction differs from the YAML parser's block scalar:" >&2
    head -40 "$WORK/yaml.diff" >&2
    exit 1
  fi
  echo "ok   workflow YAML parses; extracted script matches the parser byte for byte"
  PASS=$((PASS + 1))
else
  echo "SKIP YAML cross-check (R package 'yaml' unavailable)"
fi

if [ "$LIVE" -eq 1 ]; then
  echo "== live run against Docker Hub (namespace: ${HUB_USER:-chriskypri}) =="
  ( cd "$REPO_ROOT" && HUB_USER="${HUB_USER:-chriskypri}" bash "$SCRIPT.raw" )
  rc=$?
  echo
  if [ "$rc" -eq 0 ]; then
    echo "guard is GREEN for $(git -C "$REPO_ROOT" rev-parse --short HEAD)"
  else
    echo "guard is RED for $(git -C "$REPO_ROOT" rev-parse --short HEAD) — see above"
  fi
  exit "$rc"
fi

sed -e 's|^T_IMAGE=\$(hub_tag_epoch data\.impactncdengl latest)$|T_IMAGE=$FAKE_T_IMAGE|' \
    -e 's|^T_MODEL=\$(hub_tag_epoch impactncdengl main)$|T_MODEL=$FAKE_T_MODEL|' \
    "$SCRIPT.raw" > "$SCRIPT"

for stub in 'T_IMAGE=$FAKE_T_IMAGE' 'T_MODEL=$FAKE_T_MODEL'; do
  if ! grep -qF "$stub" "$SCRIPT"; then
    echo "FATAL: network stub '$stub' did not apply — the hub_tag_epoch call site changed." >&2
    echo "       Update the sed patterns in $0 before trusting this test." >&2
    exit 1
  fi
done
if grep -q 'hub_tag_epoch data\.impactncdengl\|hub_tag_epoch impactncdengl' "$SCRIPT"; then
  echo "FATAL: a live hub_tag_epoch call survived stubbing." >&2
  exit 1
fi

# --- harness ---------------------------------------------------------------
# run_guard <repo-dir> <T_IMAGE epoch> <T_MODEL epoch>  -> sets OUT and RC.
# Deliberately NOT `OUT=$(run_guard ...)`: an exit status assigned inside a
# command substitution dies with the subshell, so RC would silently stay at the
# previous case's value and every rc assertion would be vacuous.
RC=0
OUT=""
run_guard() {
  ( cd "$1" && HUB_USER=stub FAKE_T_IMAGE="$2" FAKE_T_MODEL="$3" bash "$SCRIPT" ) \
    > "$WORK/out.txt" 2>&1
  RC=$?
  OUT=$(cat "$WORK/out.txt")
}

# check <name> <expected rc> <expected substring|-> <forbidden substring|->
check() {
  local name="$1" want_rc="$2" want="$3" forbid="$4" out="$5" ok=1
  [ "$RC" -eq "$want_rc" ] || { ok=0; echo "    rc=$RC, wanted $want_rc"; }
  if [ "$want" != "-" ] && ! printf '%s' "$out" | grep -qF "$want"; then
    ok=0; echo "    missing expected text: $want"
  fi
  if [ "$forbid" != "-" ] && printf '%s' "$out" | grep -qF "$forbid"; then
    ok=0; echo "    found forbidden text: $forbid"
  fi
  if [ "$ok" -eq 1 ]; then PASS=$((PASS + 1)); echo "ok   $name"
  else FAIL=$((FAIL + 1)); echo "FAIL $name"; printf '%s\n' "$out" | sed 's/^/     | /'; fi
}

epoch() { date -u -d "$1" +%s; }

# ===========================================================================
# Part 1 — against the REAL repository history
# ===========================================================================
echo "== real repo history =="

# The actual Docker Hub push time of data.impactncdengl:latest that red-X'd
# main, and the EOL-only commit that did it.
T_REAL_IMAGE=$(epoch '2026-07-23T22:21:04Z')
EOL_COMMIT=96524fa   # "chore: enforce LF repo-wide via `* text=auto eol=lf`"
CKUTILS_COMMIT=121129a  # last real content change to the prerequisite Dockerfile

run_guard "$REPO_ROOT" "$T_REAL_IMAGE" "$(epoch '2026-08-08T23:22:05Z')"
check "EOL-only renormalisation no longer trips CHECK 1" 0 "CHECK 1 OK" "is STALE" "$OUT"
check "  ...and says why it was ignored" 0 "line-ending-only change" - "$OUT"

# Pin the image to one second before the CKutils bump landed: a genuine
# content change to the prerequisite Dockerfile must still fail.
T_BEFORE_CKUTILS=$(( $(git -C "$REPO_ROOT" log -1 --format=%ct "$CKUTILS_COMMIT") - 1 ))
run_guard "$REPO_ROOT" "$T_BEFORE_CKUTILS" "$T_BEFORE_CKUTILS"
check "genuine content change still trips CHECK 1" 1 \
      "Dockerfile.prerequisite.IMPACTncdENGL" - "$OUT"
check "  ...and reports the content commit, not the EOL one" 1 \
      "$(git -C "$REPO_ROOT" log -1 --format=%h "$CKUTILS_COMMIT")" \
      "$(git -C "$REPO_ROOT" log -1 --format=%h "$EOL_COMMIT") 2026-08-08" "$OUT"

# ===========================================================================
# Part 2 — synthetic histories (the opt-in-vs-blanket property)
# ===========================================================================
echo "== synthetic histories =="

HARD_FILES="
docker_setup/Dockerfile.data.IMPACTncdENGL
docker_setup/download_zenodo_data.R
docker_setup/Dockerfile.prerequisite.IMPACTncdENGL
docker_setup/apt-packages.txt
docker_setup/r-packages.txt
docker_setup/install_packages.sh
docker_setup/entrypoint.sh
Rpackage/IMPACTncd_England_model_pkg/R/ZenodoManager_class.R
"

# new_repo <dir> — all guarded files present, committed CRLF at a fixed old date
new_repo() {
  local d="$1" f
  mkdir -p "$d" && git -C "$d" init -q -b main
  git -C "$d" config user.email t@t.invalid && git -C "$d" config user.name t
  # No .gitattributes: keep git out of the way so blobs are exactly as written.
  for f in $HARD_FILES; do
    mkdir -p "$d/$(dirname "$f")"
    printf 'line one\r\nline two\r\n' > "$d/$f"
  done
  git -C "$d" add -A
  GIT_AUTHOR_DATE='2026-01-01T00:00:00Z' GIT_COMMITTER_DATE='2026-01-01T00:00:00Z' \
    git -C "$d" commit -qm 'initial'
}

# commit_in <dir> <date> <msg> — commit whatever is staged/modified
commit_in() {
  local d="$1" when="$2" msg="$3"
  git -C "$d" add -A
  GIT_AUTHOR_DATE="$when" GIT_COMMITTER_DATE="$when" git -C "$d" commit -qm "$msg"
}

# The whole synthetic timeline must sit in the PAST, or every commit trips the
# future-clock branch and is excluded from the gate — which would turn each
# "must fail" case green for entirely the wrong reason.
#   initial commit 2026-01-01 < image 2026-03-01 < branch 2026-04-01
#   < later commit 2026-05-01 < model image 2026-06-01 < today
T_SYNTH_IMAGE=$(epoch '2026-03-01T00:00:00Z')
T_LATE=$(epoch '2026-06-01T00:00:00Z')          # keeps CHECK 2 green
if [ "$T_LATE" -ge "$(date -u +%s)" ]; then
  echo "FATAL: the synthetic timeline has aged into the future — move the fixture dates back." >&2
  exit 1
fi

# (C) EOL-only change to an EOL-INSENSITIVE file -> must pass
R="$WORK/c"; new_repo "$R"
printf 'line one\nline two\n' > "$R/docker_setup/Dockerfile.prerequisite.IMPACTncdENGL"
commit_in "$R" '2026-05-01T00:00:00Z' 'renormalise EOLs'
run_guard "$R" "$T_SYNTH_IMAGE" "$T_LATE"
check "EOL-only on Dockerfile (insensitive) -> pass" 0 "CHECK 1 OK" "is STALE" "$OUT"

# (D) the same EOL-only change to a STRICT file -> must still fail
R="$WORK/d"; new_repo "$R"
printf 'line one\nline two\n' > "$R/docker_setup/install_packages.sh"
commit_in "$R" '2026-05-01T00:00:00Z' 'renormalise EOLs'
run_guard "$R" "$T_SYNTH_IMAGE" "$T_LATE"
check "EOL-only on install_packages.sh (strict) -> fail" 1 "install_packages.sh" - "$OUT"

# (E) real content change to an EOL-insensitive file -> must fail
R="$WORK/e"; new_repo "$R"
printf 'line one\r\nline two\r\nRUN something new\r\n' \
  > "$R/docker_setup/Dockerfile.data.IMPACTncdENGL"
commit_in "$R" '2026-05-01T00:00:00Z' 'real change'
run_guard "$R" "$T_SYNTH_IMAGE" "$T_LATE"
check "content change on Dockerfile (insensitive) -> fail" 1 \
      "Dockerfile.data.IMPACTncdENGL" - "$OUT"

# (F) EOL-only change arriving via a --no-ff merge: the walk is --first-parent,
#     so the diff must be taken against the FIRST parent (old main), not the
#     feature branch, or the merge looks like a no-op and is skipped wrongly.
R="$WORK/f"; new_repo "$R"
git -C "$R" checkout -q -b feature
printf 'line one\nline two\n' > "$R/docker_setup/Dockerfile.prerequisite.IMPACTncdENGL"
commit_in "$R" '2026-04-01T00:00:00Z' 'renormalise EOLs on a branch'
git -C "$R" checkout -q main
GIT_AUTHOR_DATE='2026-05-01T00:00:00Z' GIT_COMMITTER_DATE='2026-05-01T00:00:00Z' \
  git -C "$R" merge -q --no-ff -m 'Merge feature' feature
run_guard "$R" "$T_SYNTH_IMAGE" "$T_LATE"
check "EOL-only merged via --no-ff -> pass" 0 "CHECK 1 OK" "is STALE" "$OUT"

# (G) same shape, but the merge carries a REAL change -> must fail
R="$WORK/g"; new_repo "$R"
git -C "$R" checkout -q -b feature
printf 'line one\r\nRUN new\r\n' > "$R/docker_setup/Dockerfile.prerequisite.IMPACTncdENGL"
commit_in "$R" '2026-04-01T00:00:00Z' 'real change on a branch'
git -C "$R" checkout -q main
GIT_AUTHOR_DATE='2026-05-01T00:00:00Z' GIT_COMMITTER_DATE='2026-05-01T00:00:00Z' \
  git -C "$R" merge -q --no-ff -m 'Merge feature' feature
run_guard "$R" "$T_SYNTH_IMAGE" "$T_LATE"
check "real change merged via --no-ff -> fail" 1 "Dockerfile.prerequisite" - "$OUT"

# (H) an EOL-only commit must not mask an OLDER real change that is also newer
#     than the image: the walk has to keep going past the skipped commit.
R="$WORK/h"; new_repo "$R"
printf 'line one\r\nRUN new\r\n' > "$R/docker_setup/Dockerfile.prerequisite.IMPACTncdENGL"
commit_in "$R" '2026-04-01T00:00:00Z' 'real change'
printf 'line one\nRUN new\n' > "$R/docker_setup/Dockerfile.prerequisite.IMPACTncdENGL"
commit_in "$R" '2026-05-01T00:00:00Z' 'renormalise EOLs'
run_guard "$R" "$T_SYNTH_IMAGE" "$T_LATE"
check "EOL-only does not mask an older real change" 1 "Dockerfile.prerequisite" - "$OUT"

# (I) all files older than the image -> pass, and CHECK 2 still evaluated
R="$WORK/i"; new_repo "$R"
run_guard "$R" "$T_SYNTH_IMAGE" "$T_LATE"
check "everything older than the image -> pass" 0 "CHECK 2 OK" "is STALE" "$OUT"

# (J) CHECK 2 unaffected: model image older than the data image -> fail
R="$WORK/j"; new_repo "$R"
run_guard "$R" "$T_LATE" "$T_SYNTH_IMAGE"
check "CHECK 2 still fails on a stale model image" 1 "is OLDER than" - "$OUT"

# (K) future committer time is still excluded, not fatal
R="$WORK/k"; new_repo "$R"
printf 'line one\r\nRUN new\r\n' > "$R/docker_setup/Dockerfile.data.IMPACTncdENGL"
commit_in "$R" '2099-01-01T00:00:00Z' 'clock skew'
run_guard "$R" "$T_SYNTH_IMAGE" "$T_LATE"
check "future committer time excluded from the gate" 0 "FUTURE committer time" "is STALE" "$OUT"

echo
echo "passed: $PASS   failed: $FAIL"
[ "$FAIL" -eq 0 ]
