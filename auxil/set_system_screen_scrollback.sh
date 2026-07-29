#!/usr/bin/env bash
# set_system_screen_scrollback.sh -------------------------------------------
# Raise the system-wide GNU screen scrollback cap.
#
# /etc/screenrc ships `defscrollback 1000`. That is too small to diagnose a
# stalled run: on 2026-07-28 a Wales run lost 17h29m to a wedged tty and the
# console history explaining it had already been overwritten by the time we
# looked. 20000 lines costs a few MB of RAM per window.
#
# A user with their own ~/.screenrc is already covered (it is read after
# /etc/screenrc, so it wins). This script is for the other model users on a
# shared box who have none.
#
# Run as root:   sudo ./auxil/set_system_screen_scrollback.sh
# Takes effect for NEW screen sessions only.
set -euo pipefail

RC=/etc/screenrc
WANT=20000

[[ $EUID -eq 0 ]] || { echo "must run as root: sudo $0" >&2; exit 1; }
[[ -f $RC ]] || { echo "$RC not found" >&2; exit 1; }

BACKUP="${RC}.bak.$(date +%Y%m%d%H%M%S)"
cp -p "$RC" "$BACKUP"
echo "backed up $RC -> $BACKUP"

echo "before: $(grep -n '^defscrollback' "$RC" || echo '(no defscrollback line)')"

if grep -q '^defscrollback' "$RC"; then
  sed -i "s/^defscrollback .*/defscrollback $WANT/" "$RC"
else
  printf '\n# raised for post-mortem diagnosis of long-running jobs\ndefscrollback %s\n' "$WANT" >> "$RC"
fi

echo "after : $(grep -n '^defscrollback' "$RC")"
echo "done. New screen sessions will use the new value; existing ones are unchanged."
