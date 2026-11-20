#!/usr/bin/env bash
# Move stale outputs to $W/pre_repair/<repo path> (mv), or copy (cp -p) the in-place-updated ones.
#   usage: r2_moveaside.sh mv|cp <repo-relative path> ...
set -euo pipefail
R=/home/freiburger/Documents/prFBA; W=/home/freiburger/Documents/prFBA_rescore_2026-11
op=$1; shift
for p in "$@"; do
  src=$R/$p; dst=$W/pre_repair/$p
  [ -e "$src" ] || { echo "[aside] MISSING $p"; continue; }
  [ -e "$dst" ] && { echo "[aside] EXISTS (kept) $dst"; [ "$op" = mv ] && { echo "[aside] refusing to mv over an existing copy" >&2; exit 3; }; continue; }
  mkdir -p "$(dirname "$dst")"
  if [ "$op" = mv ]; then mv "$src" "$dst"; else cp -p "$src" "$dst"; fi
  echo "[aside] $op $p ($(stat -c %s "$dst") B, sha256 $(sha256sum "$dst" | cut -c1-16))"
done
