#!/usr/bin/env bash
# d100_pipeline.sh -- the depth-100 real-amplicon deployment pipeline, shared by
# run_codif_d100.sh (deployment A) and run_emilykin_d100.sh (deployment B). Not run directly.
#
# Required env: NAME FASTA TAXCSV OUT SYN     Optional: STAGE=all|align|build (default all),
#   WORKERS=12 (align pool), PROCS=8 (merge pool), IOW=8 (fetch threads), MIN_AVAIL_GB=32,
#   GENOME_CACHE_DIR (default prFBA/bvbrc_cache/kbase_genome_cache; a smoke test points it at scratch)
#
# Stages (each resumable; reruns only what is missing):
#   1 align     align_hits.py --force-backend cuda --k2 100 --prefilter-k 500 --gpu-cand 200
#               --workers 12 --save-shortlist -> $OUT/asv_top20_alignment_hits.json (key `top20`
#               holds up to 100 hits; the name is the pipeline's fixed schema). A crash after
#               stage 1/2 resumes from the saved shortlist (--shortlist-in $OUT).
#   2 tieset    tieset_gene_cache.py: NEW small PGFam cache $OUT/gene_families_tieset.json,
#               streamed read-only from the shared bvbrc_cache/genome_gene_families.json
#   3 prefetch  prefetch_gene_families.py --cache <that small cache> (BV-BRC read-only GETs)
#   4 build     build_d100.py == build_synthetic_genomes_parallel.py + BV-BRC User-Agent fix + atomic genome-cache writes +
#               tie-set guard + gc.freeze -> $SYN (synthetic genomes) + $OUT/asv_reference_selection_d100.json
#               (shipped exact-tie union, P3 gene tie-break on: the parallel builder's defaults)
#   5 stats     deployment_stats.py -> $OUT/deployment_stats_d100.json
# Nothing outside $OUT and $SYN is written except NEW per-genome files in
# bvbrc_cache/kbase_genome_cache (the build's additive genome cache).
set -uo pipefail
PY=/home/freiburger/Documents/py_venv/bin/python
PRFBA=/home/freiburger/Documents/prFBA
D=$PRFBA/deployment_d100
cd "$PRFBA" || exit 2
: "${NAME:?}" "${FASTA:?}" "${TAXCSV:?}" "${OUT:?}" "${SYN:?}"
STAGE=${STAGE:-all}; WORKERS=${WORKERS:-12}; PROCS=${PROCS:-8}; IOW=${IOW:-8}; MIN_AVAIL_GB=${MIN_AVAIL_GB:-32}
GENOME_CACHE_DIR=${GENOME_CACHE_DIR:-$PRFBA/bvbrc_cache/kbase_genome_cache}
HITS=$OUT/asv_top20_alignment_hits.json
GC=$OUT/gene_families_tieset.json
SEL=$OUT/asv_reference_selection_d100.json

log() { echo "[$(date -Is)] [$NAME] $*"; }
for p in "$OUT" "$SYN"; do
  case "$p" in *_d100) ;; *) log "[abort] output path $p does not end in _d100 (refusing to touch 20-hit outputs)"; exit 2;; esac
done
mempeak() {  # the enclosing cgroup's peak memory so far (the systemd unit when launched by launch_dep100.sh)
  local cg; cg=$(cut -d: -f3 /proc/self/cgroup 2>/dev/null)
  if [ -r "/sys/fs/cgroup$cg/memory.peak" ]; then
    log "cgroup $cg memory.peak=$(( $(cat /sys/fs/cgroup$cg/memory.peak) / 1048576 )) MiB current=$(( $(cat /sys/fs/cgroup$cg/memory.current) / 1048576 )) MiB"
  fi
}
wait_mem() {  # cooperative: do not start a heavy stage while the machine is short of memory
  local n=0 avail
  while :; do
    avail=$(awk '/^MemAvailable:/{print int($2/1048576)}' /proc/meminfo)
    [ "$avail" -ge "$MIN_AVAIL_GB" ] && { log "MemAvailable=${avail} GB (>= ${MIN_AVAIL_GB}) -> start $1"; return 0; }
    [ $((n % 10)) -eq 0 ] && log "waiting: MemAvailable=${avail} GB < ${MIN_AVAIL_GB} GB before $1"
    n=$((n + 1)); sleep 60
  done
}
step() { local t0=$SECONDS; log "+ $*"; "$@"; local rc=$?; log "rc=$rc ($((SECONDS - t0)) s)"; mempeak; return $rc; }

align() {
  if [ -s "$HITS" ] && $PY -c "import json,sys; sys.exit(0 if 'run_stats' in json.load(open('$OUT/method_selection.json')) else 1)" 2>/dev/null; then
    log "align: complete hits already in $OUT -> skip"; return 0
  fi
  mkdir -p "$OUT"
  local resume=()
  if [ -s "$OUT/shortlist_meta.json" ] && [ -s "$OUT/edlib_top500_shortlist.npz" ]; then
    resume=(--shortlist-in "$OUT"); log "align: resuming from the saved shortlist in $OUT"
  fi
  wait_mem align
  step $PY align_hits.py --force-backend cuda --k2 100 --prefilter-k 500 --gpu-cand 200 \
       --workers "$WORKERS" --save-shortlist "${resume[@]}" \
       --fasta "$FASTA" --taxonomy-csv "$TAXCSV" --outdir "$OUT" || return $?
  # every record must hold up to 100 hits, ranks 1..n
  step $PY -c "
import json, collections, sys
d = json.load(open('$HITS'))
L = collections.Counter(len(r['top20']) for r in d.values())
ok = all([h['rank'] for h in r['top20']] == list(range(1, len(r['top20']) + 1)) for r in d.values())
print('[check] records', len(d), '| hits/record', dict(sorted(L.items())), '| ranks contiguous', ok)
sys.exit(0 if ok and max(L) <= 100 and max(L) == 100 else 3)"
}

build() {
  [ -s "$HITS" ] || { log "[abort] build: no hits at $HITS (run the align stage first)"; return 2; }
  if [ ! -s "$GC" ] || [ "$HITS" -nt "$GC" ]; then
    step $PY $D/tieset_gene_cache.py --hits "$HITS" --out "$GC" --force || return $?
  else
    log "tieset: $GC is newer than the hits -> reuse"
  fi
  step $PY prefetch_gene_families.py --hits "$HITS" --cache "$GC" --workers "$IOW" || return $?
  wait_mem build
  step $PY $D/build_d100.py --tieset "$GC.tieset.json" --hits "$HITS" --gene-cache "$GC" \
       --genome-cache-dir "$GENOME_CACHE_DIR" \
       --out-dir "$SYN" --selection-out "$SEL" --procs "$PROCS" --io-workers "$IOW" || return $?
  step $PY $D/deployment_stats.py --selection "$SEL" --hits "$HITS" --provenance-dir "$SYN/provenance" \
       --out "$OUT/deployment_stats_d100.json" --force
}

log "START stage=$STAGE FASTA=$FASTA OUT=$OUT SYN=$SYN genome_cache=$GENOME_CACHE_DIR workers=$WORKERS procs=$PROCS io=$IOW host=$(hostname) prFBA=$(git -C $PRFBA rev-parse --short HEAD)"
log "merge module: codiffusion_bioreactor/bvbrc_to_kbase_genome.py @ $(git -C /home/freiburger/Documents/codiffusion_bioreactor log -1 --format=%h -- bvbrc_to_kbase_genome.py) (dirty: $(git -C /home/freiburger/Documents/codiffusion_bioreactor status --porcelain -- bvbrc_to_kbase_genome.py | wc -l))"
rc=0
case "$STAGE" in
  align) align; rc=$? ;;
  build) build; rc=$? ;;
  all)   align && build; rc=$? ;;
  *) log "[abort] unknown STAGE=$STAGE"; rc=2 ;;
esac
log "END stage=$STAGE rc=$rc"
exit $rc
