#!/usr/bin/env bash
# Audit line endings across the repository.
#
# Reports, per branch, which TRACKED blobs contain CR bytes, grouped by
# extension, and whether .gitattributes currently claims that extension is
# text/eol=lf. A file that holds CRLF in the blob while an attribute says
# eol=lf is the pathological case: `eol=lf` governs only checkout, so the CRs
# survive it, while the clean filter strips them on staging -- so git reports
# the file as modified forever, in every worktree, even a pristine checkout.
#
# Reads blobs from the branch (git cat-file), NOT the working tree, so the
# result is a property of history rather than of whoever last checked out.
#
# Binary files are excluded by asking git itself (`git grep -I`-style check via
# the `binary` attribute and a NUL-byte probe), so a .parquet or .fst that
# happens to contain 0x0D is never reported as a line-ending problem.
#
# Usage: auxil/audit_line_endings.sh [branch ...]     (default: main v2026 polygenic)

set -uo pipefail
BRANCHES=("$@")
[ ${#BRANCHES[@]} -eq 0 ] && BRANCHES=(main v2026 polygenic)

cd "$(git rev-parse --show-toplevel)" || exit 1

for BR in "${BRANCHES[@]}"; do
  git rev-parse --verify -q "$BR" >/dev/null || { echo "!! no such branch: $BR"; continue; }
  echo "================ $BR ($(git rev-parse --short "$BR")) ================"

  total=0
  declare -A crlf_by_ext=()
  declare -A all_by_ext=()

  while IFS=$'\t' read -r mode type _ path; do
    [ "$type" = blob ] || continue
    case "$path" in */) continue;; esac

    ext="${path##*/}"
    if [[ "$ext" == *.* ]]; then ext=".${ext##*.}"; else ext="(noext)"; fi

    # Skip binaries: probe the first 8000 bytes for a NUL, as git does.
    if git cat-file blob "$BR:$path" 2>/dev/null | head -c 8000 | LC_ALL=C grep -qU $'\x00'; then
      continue
    fi

    all_by_ext["$ext"]=$(( ${all_by_ext["$ext"]:-0} + 1 ))
    if git cat-file blob "$BR:$path" 2>/dev/null | LC_ALL=C grep -qU $'\r'; then
      crlf_by_ext["$ext"]=$(( ${crlf_by_ext["$ext"]:-0} + 1 ))
      total=$(( total + 1 ))
    fi
  done < <(git ls-tree -r --full-tree "$BR")

  if [ "$total" -eq 0 ]; then
    echo "  clean: no tracked text blob contains CR"
  else
    printf '  %-12s %8s %8s   %s\n' EXT WITH_CR TOTAL ".gitattributes says"
    for ext in $(printf '%s\n' "${!crlf_by_ext[@]}" | sort); do
      # What does git think this extension is, per .gitattributes on that branch?
      probe="probe${ext}"
      [ "$ext" = "(noext)" ] && probe="probe"
      attr=$(git check-attr text eol -- "$probe" 2>/dev/null \
             | sed 's/^[^:]*: //' | paste -sd' ' - | tr -s ' ')
      printf '  %-12s %8s %8s   %s\n' "$ext" "${crlf_by_ext[$ext]}" "${all_by_ext[$ext]:-0}" "${attr:-unspecified}"
    done
    echo "  ---- $total tracked text blob(s) contain CR"
  fi
  unset crlf_by_ext all_by_ext
  echo
done

echo "=== current .gitattributes (working tree) ==="
sed 's/^/  /' .gitattributes 2>/dev/null || echo "  (none)"
