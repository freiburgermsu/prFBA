#!/usr/bin/env bash
# POST stages for all six depth-100 trees, strictly one at a time: they all write the SHARED
# PGFam cache (load -> mutate -> whole-file replace; the last writer wins), so a second
# concurrent POST would lose fetches. An flock enforces this across invocations.
# Order: original union, original legacy, k-means union, k-means legacy, random union,
# random legacy. Each tree waits for its panel's alignment (polls align.done every 5 min),
# so this can be started while later panels are still aligning. A tree whose POST finished
# (data/post.done) is skipped; a partial POST resumes stage by stage.
#   POST_TREES="..." restricts the list (same names, same order rules).
set -euo pipefail
R=/home/freiburger/Documents/prFBA
LOCK="$R/region_validation_d100/.post_cache_writer.lock"
exec 9>"$LOCK"
flock -n 9 || { echo "[post-all] another POST driver holds $LOCK -- refusing to run two cache writers" >&2; exit 3; }
FAILED=()
TREES=(${POST_TREES:-region_validation_d100 region_validation_d100_legacy region_validation_kmeans_d100 region_validation_kmeans_d100_legacy region_validation_random_d100 region_validation_random_d100_legacy})
for T in "${TREES[@]}"; do
  U="${T%_legacy}"
  if [ -s "$R/$T/data/post.done" ]; then echo "[post-all] skip $T (post.done present)"; continue; fi
  while [ ! -s "$R/$U/data/hits/raw/align.done" ]; do
    echo "[post-all] $T waits for $U alignment ($(date -Is))"; sleep 300
  done
  if [ "$T" != "$U" ]; then bash "$R/region_validation_d100/seed_legacy_hits_d100.sh" "$U"; fi
  echo "[post-all] $T POST start $(date -Is)"
  STAGE=post NO_TMUX=1 POST_LOCK_HELD=1 bash "$R/$T/scripts/run_validation.sh" full
  if "/home/freiburger/Documents/py_venv/bin/python" "$R/region_validation_d100/check_d100.py" post "$T"; then
    date -Is > "$R/$T/data/post.done"
    echo "[post-all] $T POST done $(date -Is)"
  else
    FAILED+=("$T")
    # let the next launch re-close the representative gap and re-score gene capture
    # (both stages are otherwise skipped because their outputs exist)
    rm -f "$R/$T/data/close_exclude_rep_gap.done" "$R/$T/data/gene_capture_by_region.json"
    echo "[post-all] $T POST ran but its checks FAILED (no post.done written; see above) $(date -Is)"
  fi
done
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "[post-all] finished with failed checks: ${FAILED[*]} $(date -Is)"; exit 1
fi
echo "[post-all] finished $(date -Is)"
