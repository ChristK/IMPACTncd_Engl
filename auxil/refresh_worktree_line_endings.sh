#!/usr/bin/env bash
# Rewrite every tracked file in a worktree from the index, so what is ON DISK
# matches the (normalised) blob.
#
# Why this is a separate, non-obvious step:
#
#   `git add --renormalize .` updates the INDEX and the resulting commit, but it
#   does NOT rewrite the files in your working tree. After committing, the blob
#   is LF while the file on disk is still CRLF -- and `git status` says clean,
#   because the clean filter strips the CRs when it compares. So the repository
#   looks fixed while the CRLF is still sitting in the file that docker, bash and
#   R actually read. That is precisely the bug the normalisation was meant to
#   cure: a stray \r in a YAML value breaks `docker run --mount`, and a CRLF
#   shebang fails as `bad interpreter: ...^M`.
#
#   `git checkout -- .` will NOT fix it either, because git considers those files
#   unchanged and skips them. `git checkout-index -f -a` forces the rewrite.
#
# Refuses to run on a dirty worktree, because forcing a rewrite from the index
# would discard uncommitted work.
#
# Usage: auxil/refresh_worktree_line_endings.sh [worktree ...]

set -uo pipefail
WTS=("$@")
if [ ${#WTS[@]} -eq 0 ]; then
  mapfile -t WTS < <(git worktree list --porcelain | awk '/^worktree /{print $2}')
fi

crcount() { tr -dc '\r' < "$1" 2>/dev/null | wc -c; }

for WT in "${WTS[@]}"; do
  [ -d "$WT" ] || { echo "!! missing: $WT"; continue; }
  echo "=== $WT ($(git -C "$WT" branch --show-current 2>/dev/null)) ==="

  dirty=$(git -C "$WT" status --porcelain | wc -l)
  if [ "$dirty" -ne 0 ]; then
    echo "  REFUSING: $dirty uncommitted change(s); commit or discard them first."
    git -C "$WT" status --porcelain | head -5 | sed 's/^/    /'
    continue
  fi

  before=0
  while IFS= read -r -d '' f; do
    [ -f "$WT/$f" ] || continue
    case "$f" in *.jpg|*.jpeg|*.png|*.pdf|*.docx|*.xlsx|*.zip|*.gz) continue;; esac
    n=$(crcount "$WT/$f")
    [ "$n" -gt 0 ] && before=$(( before + 1 ))
  done < <(git -C "$WT" ls-files -z)

  if [ "$before" -eq 0 ]; then
    echo "  already LF on disk; nothing to do"
    continue
  fi

  ( cd "$WT" && git checkout-index -f -a )

  after=0
  while IFS= read -r -d '' f; do
    [ -f "$WT/$f" ] || continue
    case "$f" in *.jpg|*.jpeg|*.png|*.pdf|*.docx|*.xlsx|*.zip|*.gz) continue;; esac
    n=$(crcount "$WT/$f")
    [ "$n" -gt 0 ] && { after=$(( after + 1 )); echo "    STILL CRLF: $f"; }
  done < <(git -C "$WT" ls-files -z)

  echo "  files with CR on disk: $before -> $after"
  echo "  worktree still clean: $([ "$(git -C "$WT" status --porcelain | wc -l)" -eq 0 ] && echo yes || echo NO)"
done
