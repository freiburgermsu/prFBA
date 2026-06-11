#!/usr/bin/env bash
# =============================================================================
# run_validation.sh -- orchestrator for the 16S hypervariable-region validation
# =============================================================================
# DESIGN: /home/freiburger/Documents/prFBA/region_validation/DESIGN.md
#   USER DECISIONS honored (DESIGN banner):
#     * NO self-exclusion -- include-self ONLY (one selection scenario).
#     * archaea KEPT & stratified by domain (only V4 / V4-V5 amplify them).
#     * ALL intermediates live under region_validation/data/ -- every path is
#       sourced from _common.PATHS, NEVER hardcoded in this orchestrator.
#   Every stage is RESUMABLE: it short-circuits when its output already exists.
#   Pass --force to ignore checkpoints and rebuild every stage from scratch.
#
# Usage:
#   run_validation.sh pilot          # MODE=pilot -> --n 500   (whole-pipeline smoke test)
#   run_validation.sh full           # MODE=full  -> --n 10000 (production run)
#   run_validation.sh pilot --force  # ignore checkpoints, rebuild everything
#
# Stage order (DESIGN §9):
#   FRONT-LOAD (sequential, resumable -- everything before alignment + truth):
#     select_10k.py --n N
#       -> insilico_pcr.py --fasta <src16s> --combined-fasta <combined> --emax 3
#       -> dedup_amplicons.py
#       -> build_truth.py
#   PARALLEL (two detached tmux jobs launched together, joined before POST):
#     align_hits.py (GPU + 60 CPU workers, --force-backend cuda)   [tmux]
#       ||  prefetch_pgfams.py (BV-BRC API: truth pass now, reps pass in POST)  [tmux]
#   POST (pure, fast, sequential -- no GPU):
#     fan hits md5->{genome_id}__{region} via amplicon_md5_fanout
#       -> select_references.py --gene-provider auto
#       -> prefetch_pgfams.py   (reps pass, now that selection exists)
#       -> score_taxacc.py
#       -> score_genecap.py
#       -> build_region_meta.py
#       -> render_accuracy.py
#       -> render_genecap.py
#
# Long jobs run detached in tmux (per project convention) so they survive the
# controlling session closing; this script waits on their completion sentinels.
#
# NOTE on the fan-out (DESIGN §5.0 + §9):  dedup keys the unique amplicon FASTA on
# the *insert md5*, so align_hits.py emits an md5-keyed hits JSON.  Every POST
# scorer (select_references / score_taxacc / score_genecap / build_truth) keys on
# the provenance id "{genome_id}__{region}".  The POST stage therefore expands the
# md5-keyed raw hits back out to every source (genome,region) via
# amplicon_md5_fanout.json BEFORE selection -- this is the join that restores
# per-source provenance after the dedup speedup.
# =============================================================================

set -euo pipefail

# --- interpreter (sole uv-managed venv; per global rules) -------------------
PY="/home/freiburger/Documents/py_venv/bin/python"

# --- repos ------------------------------------------------------------------
PRFBA="/home/freiburger/Documents/prFBA"                       # reused modules
RV="/home/freiburger/Documents/prFBA/region_validation"        # this study
SCRIPTS="${RV}/scripts"                                         # stage scripts

# --- reused prFBA modules (run as-is, DESIGN §9) ----------------------------
INSILICO_PCR="${PRFBA}/insilico_pcr.py"
ALIGN_HITS="${PRFBA}/align_hits.py"
SELECT_REFERENCES="${PRFBA}/select_references.py"

# =============================================================================
# Parse args: MODE in {pilot,full}; optional --force
# =============================================================================
MODE="${1:-}"
FORCE=0
shift || true
for a in "$@"; do
  case "$a" in
    --force) FORCE=1 ;;
    *) echo "run_validation.sh: unknown arg '$a'" >&2; exit 2 ;;
  esac
done

case "$MODE" in
  pilot) N=500 ;;
  full)  N=10000 ;;
  *)
    echo "Usage: run_validation.sh {pilot|full} [--force]" >&2
    echo "  pilot -> --n 500   full -> --n 10000" >&2
    exit 2
    ;;
esac

FORCE_FLAG=""
[ "${FORCE}" -eq 1 ] && FORCE_FLAG="--force"

# =============================================================================
# Resolve ALL absolute paths from _common.PATHS (single source of truth).
# Printed as `KEY="value"` shell assignments and eval'd -- never hardcode a data
# path in this orchestrator (DESIGN USER DECISIONS: all intermediates under data/).
# =============================================================================
eval "$(
  PYTHONPATH="${SCRIPTS}" "${PY}" - <<'PYEOF'
from _common import DATA, FIGS, PATHS  # noqa: E402
emit = {
    "DATA": DATA,
    "FIGS": FIGS,
    "SELECTION_JSON": PATHS.selection_json,     # S1 10k manifest (truth genomes)
    "SRC16S_FASTA": PATHS.src16s_fasta,         # S1 per-operon full 16S FASTA
    "COMBINED_FASTA": PATHS.combined_fasta,     # S3 combined amplicons
    "COMBINED_STATS": PATHS.combined_stats,     # S3 sibling stats.json
    "UNIQUE_FASTA": PATHS.unique_fasta,         # S4 dedup'd unique inserts (md5-keyed)
    "AMPLICON_FANOUT": PATHS.amplicon_fanout,   # S4 {insert_md5: [genome_id__region]}
    "HITS_DIR": PATHS.hits_dir,                 # S5 hits dir (fanned hits land here)
    "HITS_JSON": PATHS.hits_json,               # the FANNED {genome_id__region}-keyed hits
    "SELECTION_OUT": PATHS.selection_out,       # S6 selection.json
    "TRUTH_JSON": PATHS.truth_json,             # S7 benchmark_truth.json
    "PGFAM_CACHE": PATHS.pgfam_cache,           # S8 prFBA-wide PGFam cache
    "TAXACC_SUMMARY_JSON": PATHS.taxacc_summary_json,    # S9
    "GENECAP_SUMMARY_JSON": PATHS.genecap_summary_json,  # S10
    "REGION_META_CSV": PATHS.region_meta_csv,            # S11
}
for k, v in emit.items():
    print(f'{k}="{v}"')
PYEOF
)"

# Raw md5-keyed hits land in a subdir; we fan them out into HITS_JSON for POST.
RAW_HITS_DIR="${HITS_DIR}/raw"
RAW_HITS_JSON="${RAW_HITS_DIR}/asv_top20_alignment_hits.json"

# tmux sentinels + logs (under DATA/hits so they are resumability checkpoints too)
ALIGN_DONE="${RAW_HITS_DIR}/align.done"
ALIGN_LOG="${RAW_HITS_DIR}/align.log"
ALIGN_TMUX="rv_align_${MODE}"
PGFAM_DONE="${DATA}/prefetch_pgfams.truth.done"
PGFAM_LOG="${DATA}/prefetch_pgfams.log"
PGFAM_TMUX="rv_pgfam_${MODE}"

mkdir -p "${DATA}" "${FIGS}" "${HITS_DIR}" "${RAW_HITS_DIR}"

# =============================================================================
# Helpers
# =============================================================================
banner() {
  echo ""
  echo "============================================================================="
  echo ">>> $*"
  echo "============================================================================="
}

# done_file PATH -> true (0) if checkpoint usable (non-empty) AND not forcing.
done_file() {
  [ "${FORCE}" -eq 0 ] && [ -s "$1" ]
}

# run_stage "label" OUTPUT_PATH -- <command...>
# Skips the command if OUTPUT_PATH already exists (resumable), unless --force.
run_stage() {
  local label="$1"; shift
  local out="$1"; shift
  if done_file "${out}"; then
    echo "[skip] ${label}: checkpoint present -> ${out}"
    return 0
  fi
  echo "[run ] ${label}"
  "$@"
  echo "[ok  ] ${label} -> ${out}"
}

echo ""
echo "#############################################################################"
echo "#  16S region validation  |  MODE=${MODE}  N=${N}  FORCE=${FORCE}"
echo "#  interpreter : ${PY}"
echo "#  data dir    : ${DATA}"
echo "#  figures dir : ${FIGS}"
echo "#  resumability: every stage skips when its output exists; --force overrides"
echo "#  scope       : include-self ONLY; archaea kept & stratified (DESIGN USER DECISIONS)"
echo "#############################################################################"

cd "${SCRIPTS}"

# =============================================================================
# FRONT-LOAD (sequential, resumable) -- everything before alignment + truth
# =============================================================================
banner "FRONT-LOAD  (select_10k | insilico_pcr | dedup | build_truth)"
echo "Resumability: each stage below is skipped if its output file already exists."

# --- select_10k.py: deterministic N-genome selection -> manifest + src16s FASTA
run_stage "select_10k.py (--n ${N})" "${SRC16S_FASTA}" \
  "${PY}" "${SCRIPTS}/select_10k.py" --n "${N}" ${FORCE_FLAG}

# --- insilico_pcr.py: panel mode (both-strand, edits) -> combined amplicons + stats
#     emits <combined>.stats.json as a sibling automatically; --emax 3.
run_stage "insilico_pcr.py (panel mode, --emax 3)" "${COMBINED_FASTA}" \
  "${PY}" "${INSILICO_PCR}" \
    --fasta "${SRC16S_FASTA}" \
    --combined-fasta "${COMBINED_FASTA}" \
    --emax 3

# --- dedup_amplicons.py: collapse identical inserts -> unique fasta + fan-out
run_stage "dedup_amplicons.py" "${UNIQUE_FASTA}" \
  "${PY}" "${SCRIPTS}/dedup_amplicons.py" ${FORCE_FLAG}

# --- build_truth.py: benchmark truth (lineage_for over provenance) ----------
run_stage "build_truth.py" "${TRUTH_JSON}" \
  "${PY}" "${SCRIPTS}/build_truth.py" ${FORCE_FLAG}

# =============================================================================
# PARALLEL STAGE -- two detached tmux jobs launched together
#   align_hits.py     (GPU + 60 CPU workers, --force-backend cuda; the one heavy job)
#   prefetch_pgfams.py (BV-BRC API; truth pass now, reps pass deferred to POST)
# Each writes a `.done` sentinel on success; we then wait on both.
# =============================================================================
banner "PARALLEL  (align_hits.py [tmux]  ||  prefetch_pgfams.py [tmux])"
echo "Resumability: a job is NOT relaunched if its .done sentinel already exists."
echo "Long jobs run detached in tmux (survive session close, per project convention)."

# --- align_hits.py on the UNIQUE (md5-keyed) amplicons (detached tmux) -------
#     Writes <RAW_HITS_DIR>/asv_top20_alignment_hits.json keyed by insert md5;
#     POST fans it out to the provenance-keyed HITS_JSON.
if done_file "${ALIGN_DONE}" && done_file "${RAW_HITS_JSON}"; then
  echo "[skip] align_hits.py: sentinel + raw hits present -> ${RAW_HITS_JSON}"
elif tmux has-session -t "${ALIGN_TMUX}" 2>/dev/null; then
  echo "[wait] align_hits.py: tmux session '${ALIGN_TMUX}' already running -> joined below"
else
  echo "[run ] align_hits.py in tmux '${ALIGN_TMUX}' (GPU/60-worker)"
  [ "${FORCE}" -eq 1 ] && rm -f "${ALIGN_DONE}"
  tmux new-session -d -s "${ALIGN_TMUX}" \
    "set -o pipefail; \
     '${PY}' '${ALIGN_HITS}' \
       --fasta   '${UNIQUE_FASTA}' \
       --outdir  '${RAW_HITS_DIR}' \
       --prefilter-k 500 --k2 20 --gpu-cand 100 --workers 60 \
       --min-ref-len 200 --force-backend cuda \
       2>&1 | tee '${ALIGN_LOG}' \
     && echo done > '${ALIGN_DONE}'"
fi

# --- prefetch_pgfams.py truth pass (10k/N source genomes) (detached tmux) ----
#     prefetch_pgfams.py runs BOTH passes in one call (truth then reps); the reps
#     pass is skipped gracefully while SELECTION_OUT is absent, so this launch
#     warms the source/truth gene sets alongside alignment.  It is re-run in POST
#     (after selection) to fetch the distinct selected reps -- resumable against
#     the shared cache, so it only fetches the new reps then.
if done_file "${PGFAM_DONE}"; then
  echo "[skip] prefetch_pgfams.py (truth pass): sentinel present -> ${PGFAM_CACHE}"
elif tmux has-session -t "${PGFAM_TMUX}" 2>/dev/null; then
  echo "[wait] prefetch_pgfams.py: tmux session '${PGFAM_TMUX}' already running -> joined below"
else
  echo "[run ] prefetch_pgfams.py (truth pass) in tmux '${PGFAM_TMUX}'"
  [ "${FORCE}" -eq 1 ] && rm -f "${PGFAM_DONE}"
  tmux new-session -d -s "${PGFAM_TMUX}" \
    "set -o pipefail; \
     '${PY}' '${SCRIPTS}/prefetch_pgfams.py' ${FORCE_FLAG} \
       2>&1 | tee '${PGFAM_LOG}' \
     && echo done > '${PGFAM_DONE}'"
fi

# --- join: block until both detached jobs finish ----------------------------
echo ""
echo "[wait] joining parallel jobs (align + pgfam-truth) -- polling every 30s..."
while tmux has-session -t "${ALIGN_TMUX}" 2>/dev/null \
   || tmux has-session -t "${PGFAM_TMUX}" 2>/dev/null; do
  sleep 30
done

# Fail loudly if the heavy job did not produce its raw hits output.
if ! done_file "${RAW_HITS_JSON}"; then
  echo "[FAIL] align_hits.py produced no raw hits JSON (${RAW_HITS_JSON}); see ${ALIGN_LOG}" >&2
  exit 1
fi
echo "[ok  ] parallel stage complete: raw hits=${RAW_HITS_JSON}"

# =============================================================================
# POST (sequential, pure, fast -- no GPU)
# =============================================================================
banner "POST  (fan-out | select | reps prefetch | score | meta | render)"
echo "Resumability: each stage below skips when its output file already exists."

# --- fan hits md5 -> {genome_id}__{region} via amplicon_md5_fanout ----------
#     The dedup'd align output is keyed by insert md5; every downstream scorer
#     keys on the provenance id "{genome_id}__{region}".  Expand each md5 record
#     to all of its source (genome,region) members so selection + scoring see one
#     record per amplifying (genome,region).  Writes HITS_JSON (what POST reads).
run_stage "fan_out hits (md5 -> genome_id__region)" "${HITS_JSON}" \
  "${PY}" - "${RAW_HITS_JSON}" "${AMPLICON_FANOUT}" "${HITS_JSON}" <<'PYEOF'
import json, os, sys
raw_path, fanout_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
raw = json.load(open(raw_path))            # {insert_md5: hit_record}
fanout = json.load(open(fanout_path))      # {insert_md5: [genome_id__region, ...]}
fanned, n_md5, n_missing = {}, 0, 0
for md5, members in fanout.items():
    rec = raw.get(md5)
    if rec is None:
        n_missing += 1
        continue
    n_md5 += 1
    for member in members:                 # member == "{genome_id}__{region}"
        fanned[member] = rec               # share the same top20 record by provenance
tmp = out_path + ".tmp"
with open(tmp, "w") as fh:
    json.dump(fanned, fh)
os.replace(tmp, out_path)
print(f"[fanout] {n_md5} md5 hit-records -> {len(fanned)} (genome,region) amplicons "
      f"({n_missing} md5 in fan-out had no aligned record)")
PYEOF

# --- PREFETCH selection-candidate gene-sets FIRST (pilot lesson) -------------
#     select_references --gene-provider auto otherwise lazy-fetches every uncached
#     stage-7 candidate genome AND rewrites the whole PGFam cache per miss -> an
#     O(N^2) death-spiral on diverse references (37 min and ballooning in the pilot).
#     Prefetch them (threaded, incremental save-every-200) so selection has 0 misses.
echo "[run ] prefetch selection candidates (stage-7 gene-families) -> ${PGFAM_CACHE}"
"${PY}" "${PRFBA}/prefetch_gene_families.py" --hits "${HITS_JSON}" --cache "${PGFAM_CACHE}" --workers 16
echo "[ok  ] selection candidates prefetched"

# --- select_references.py: include-self selection (single scenario) ----------
#     --gene-provider auto, reading the shared PGFam cache.  Output keyed by the
#     same "{genome_id}__{region}" provenance ids as HITS_JSON.
run_stage "select_references.py (--gene-provider auto)" "${SELECTION_OUT}" \
  "${PY}" "${SELECT_REFERENCES}" \
    --hits "${HITS_JSON}" \
    --out  "${SELECTION_OUT}" \
    --gene-provider auto \
    --gene-cache "${PGFAM_CACHE}"

# --- prefetch_pgfams.py reps pass: fetch the distinct selected reps ----------
#     Now that SELECTION_OUT exists, the same script fetches the reps gene sets
#     (truth genomes are already cached -> only the new reps are fetched).
#     Resumable: idempotent against the shared cache.
echo "[run ] prefetch_pgfams.py (reps pass)"
"${PY}" "${SCRIPTS}/prefetch_pgfams.py" ${FORCE_FLAG}
echo "[ok  ] reps PGFam cache warmed -> ${PGFAM_CACHE}"

# --- prefetch the FULL rep universe (every top-20 hit genome) so score_genecap's
#     exclude-self re-selection + gene-capture unions never lazy-fetch + rewrite the
#     whole cache per miss (the pilot O(N^2) trap, one level deeper than selection).
echo "[run ] prefetch_universe.py (all top-20 hit genomes -> ${PGFAM_CACHE})"
"${PY}" "${SCRIPTS}/prefetch_universe.py"
echo "[ok  ] rep universe prefetched"

# --- score_taxacc.py: taxonomic-accuracy scoring -> summary + long tables ----
run_stage "score_taxacc.py" "${TAXACC_SUMMARY_JSON}" \
  "${PY}" "${SCRIPTS}/score_taxacc.py" ${FORCE_FLAG}

# --- score_genecap.py: gene-capture scoring -> summary + long tables ---------
run_stage "score_genecap.py" "${GENECAP_SUMMARY_JSON}" \
  "${PY}" "${SCRIPTS}/score_genecap.py" ${FORCE_FLAG}

# --- build_region_meta.py: region metadata (+ domain, extract_rate, stamps) --
run_stage "build_region_meta.py" "${REGION_META_CSV}" \
  "${PY}" "${SCRIPTS}/build_region_meta.py" ${FORCE_FLAG}

# --- render: figures (cheap, self-overwriting; always run so they refresh) ---
banner "RENDER  (accuracy + gene-capture figures)"
echo "[run ] render_accuracy.py"
"${PY}" "${SCRIPTS}/render_accuracy.py" ${FORCE_FLAG}
echo "[run ] render_genecap.py"
"${PY}" "${SCRIPTS}/render_genecap.py" ${FORCE_FLAG}

# =============================================================================
banner "DONE  (MODE=${MODE}, N=${N})"
echo "Tables : ${DATA}"
echo "Figures: ${FIGS}"
echo "Hits   : ${HITS_JSON}  (fanned; raw md5-keyed: ${RAW_HITS_JSON})"
echo "Re-run is fully resumable; add --force to rebuild from scratch."
