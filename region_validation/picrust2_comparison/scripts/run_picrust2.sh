#!/bin/bash
# Run PICRUSt2 2.6.3 (default bacteria + archaea references, EC and KO traits,
# max-parsimony HSP) over every batch prepared by prepare_inputs.py --batches.
# Resumable: a batch with a DONE marker is skipped; a failed batch is re-run.
#
#   bash run_picrust2.sh <inputs_dir> <outputs_dir> [parallel_batches] [procs_per_batch]
set -euo pipefail
IN=$(realpath "$1"); OUT=$(realpath -m "$2"); PAR=${3:-2}; PROCS=${4:-30}
source "$HOME/opt/picrust2/activate.sh"
mkdir -p "$OUT"
{
  echo "picrust2 $(python -c 'import picrust2; print(picrust2.__version__)')"
  epa-ng --version; echo "gappa $(gappa --version 2>&1 | head -1)"
  hmmalign -h | sed -n 2p; Rscript -e 'cat("castor", as.character(packageVersion("castor")), "\n")'
} > "$OUT/versions.txt" 2>&1

run_batch() {
  local b=$1 in_dir="$IN/$1" out_dir="$OUT/$1"
  [ -f "$out_dir/DONE" ] && { echo "[run_picrust2] $b already done"; return 0; }
  rm -rf "$out_dir"
  echo "[run_picrust2] $b start $(date +%T)"
  if picrust2_pipeline.py -s "$in_dir/study_seqs.fna" -i "$in_dir/study_seqs_abundance.tsv" \
       -o "$out_dir" -p "$PROCS" --in_traits EC,KO --no_pathways --remove_intermediate \
       --verbose > "$OUT/$b.log" 2>&1; then
    touch "$out_dir/DONE"; echo "[run_picrust2] $b done $(date +%T)"
  else
    echo "[run_picrust2] $b FAILED (see $OUT/$b.log)"; return 1
  fi
}
export -f run_batch; export IN OUT PROCS

ls -d "$IN"/batch_* | xargs -n1 basename | xargs -P "$PAR" -I{} bash -c 'run_batch {}'
