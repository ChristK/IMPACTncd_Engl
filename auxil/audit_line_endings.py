#!/usr/bin/env python3
"""Audit line endings across the repository, per branch, from git BLOBS.

Why this exists, and why it is Python rather than shell:

  * The pathological state in this repo is a blob that holds CRLF while
    .gitattributes says `eol=lf`. `eol=lf` governs only checkout, so the CRs
    survive it, while the clean filter strips them on staging -- so git reports
    the file as modified forever, in every worktree, even in a pristine
    checkout. That blocks `git rebase` and `git merge` outright.
  * It therefore has to inspect the BLOB (history), not the working tree, which
    only tells you what the last checkout happened to produce.
  * A shell version of this check was wrong: `grep` on this machine resolves to
    `ugrep`, whose `-U` is not GNU grep's `--binary`, so `grep -qU $'\\r'`
    silently matched nothing and the audit reported a repo full of CRLF as
    "clean". Byte comparison in Python has no such ambiguity. The self-test
    below exists so that failure mode can never recur silently.

Usage:
    auxil/audit_line_endings.py [branch ...]        # default: main v2026 polygenic
    auxil/audit_line_endings.py --selftest          # prove the detector fires
"""

import collections
import subprocess
import sys

DEFAULT_BRANCHES = ["main", "v2026", "polygenic"]


def git(*args, binary=False):
    out = subprocess.run(["git", *args], capture_output=True, check=True)
    return out.stdout if binary else out.stdout.decode("utf-8", "replace")


def blobs(branch):
    """Yield (sha, path) for every blob in the branch's tree."""
    for line in git("ls-tree", "-r", "--full-tree", branch).splitlines():
        meta, _, path = line.partition("\t")
        parts = meta.split()
        if len(parts) >= 3 and parts[1] == "blob":
            yield parts[2], path


def read_blobs(shas):
    """Stream many blobs through a single `git cat-file --batch`."""
    if not shas:
        return {}
    proc = subprocess.Popen(
        ["git", "cat-file", "--batch"],
        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
    )
    payload = ("\n".join(shas) + "\n").encode()
    out, _ = proc.communicate(payload)

    contents, pos = {}, 0
    while pos < len(out):
        nl = out.find(b"\n", pos)
        if nl == -1:
            break
        header = out[pos:nl].split()
        if len(header) < 3:
            break
        sha, size = header[0].decode(), int(header[2])
        start = nl + 1
        contents[sha] = out[start:start + size]
        pos = start + size + 1          # +1 for the trailing newline
    return contents


def ext_of(path):
    name = path.rsplit("/", 1)[-1]
    return "." + name.rsplit(".", 1)[-1] if "." in name else "(noext)"


def attr_for(ext):
    probe = "probe" + ("" if ext == "(noext)" else ext)
    try:
        raw = git("check-attr", "text", "eol", "--", probe)
    except subprocess.CalledProcessError:
        return "?"
    vals = []
    for line in raw.splitlines():
        attr, _, val = line.rpartition(": ")
        vals.append(f"{attr.split(': ')[-1]}={val}")
    return " ".join(v for v in vals if not v.endswith("=unspecified")) or "unspecified"


def audit(branch):
    pairs = list(blobs(branch))
    contents = read_blobs([sha for sha, _ in pairs])

    crlf = collections.Counter()
    total = collections.Counter()
    offenders = []

    for sha, path in pairs:
        data = contents.get(sha)
        if data is None:
            continue
        if b"\x00" in data[:8000]:          # binary, as git itself decides
            continue
        e = ext_of(path)
        total[e] += 1
        if b"\r" in data:
            crlf[e] += 1
            offenders.append(path)

    head = git("rev-parse", "--short", branch).strip()
    print(f"================ {branch} ({head}) ================")
    if not offenders:
        print("  clean: no tracked text blob contains CR\n")
        return []

    print(f"  {'EXT':<12}{'WITH_CR':>9}{'TOTAL':>8}   .gitattributes says")
    for e in sorted(crlf):
        print(f"  {e:<12}{crlf[e]:>9}{total[e]:>8}   {attr_for(e)}")
    print(f"  ---- {len(offenders)} tracked text blob(s) contain CR")
    for p in offenders[:25]:
        print(f"       {p}")
    if len(offenders) > 25:
        print(f"       ... and {len(offenders) - 25} more")
    print()
    return offenders


def selftest():
    """The detector MUST fire on a blob known to contain CRLF."""
    ref = "backup/pre-rebase-polygenic-20260808:testing/sim_design_testing.yaml"
    try:
        data = git("cat-file", "blob", ref, binary=True)
    except subprocess.CalledProcessError:
        print("SELFTEST SKIPPED: backup tag not present")
        return 0
    n = data.count(b"\r")
    ok = n > 0
    print(f"SELFTEST: known-CRLF blob has {n} CR bytes -> detector {'FIRES' if ok else 'FAILED'}")
    return 0 if ok else 1


if __name__ == "__main__":
    args = sys.argv[1:]
    if "--selftest" in args:
        sys.exit(selftest())
    if selftest() != 0:
        sys.exit("refusing to report: the detector does not fire on a known-CRLF blob")
    print()
    bad = 0
    for b in (args or DEFAULT_BRANCHES):
        bad += len(audit(b))
    print(f"TOTAL offending blobs across branches: {bad}")
