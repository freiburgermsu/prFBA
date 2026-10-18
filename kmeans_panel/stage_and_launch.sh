#!/bin/bash
# Stage the k-means panel into the isolated pipeline tree and launch the 16S arm.
# The driver's first stage (select_10k.py) is skipped because its outputs already
# exist -- that is how the alternative panel is injected without touching the
# original selector.
set -euo pipefail
KM=~/Documents/prFBA/kmeans_panel
RVK=~/Documents/prFBA/region_validation_kmeans
test -s "$KM/selection_kmeans.json" || { echo "panel not built yet"; exit 1; }
mkdir -p "$RVK/data"
cp "$KM/selection_kmeans.json" "$RVK/data/selection_10k.json"
cp "$KM/src16s_kmeans.fasta"   "$RVK/data/src16s_10k.fasta"
echo "staged: $(grep -c '^>' "$RVK/data/src16s_10k.fasta") genomes into $RVK/data"
# safety: the original tree must be untouched
grep -q 'region_validation_kmeans' "$RVK/scripts/_common.py" || { echo "BASE not redirected!"; exit 1; }
grep -q '"/home/freiburger/Documents/prFBA/region_validation"' ~/Documents/prFBA/region_validation/scripts/_common.py \
  || { echo "ORIGINAL _common.py was modified - aborting"; exit 1; }
echo "safety checks passed; launching 16S arm (full mode)"
cd "$RVK/scripts"
exec bash run_validation.sh full
