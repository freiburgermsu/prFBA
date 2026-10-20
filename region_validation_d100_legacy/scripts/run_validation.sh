#!/usr/bin/env bash
# =============================================================================
# run_validation.sh -- orchestrator for the 16S hypervariable-region validation
#                      DEPTH-100 RERUN TREE: region_validation_d100_legacy
# =============================================================================
# DESIGN: /home/freiburger/Documents/prFBA/region_validation/DESIGN.md
#   USER DECISIONS honored (DESIGN banner):
#     * NO self-exclusion -- include-self ONLY (one selection scenario).
#     * archaea KEPT & stratified by domain (only V4 / V4-V5 amplify them).
#     * ALL intermediates live under <tree>/data/ -- every path is
#       sourced from _common.PATHS, NEVER hardcoded in this orchestrator.
#   Every stage is RESUMABLE: it short-circuits when its output already exists.
#
# DEPTH-100 TREE (copied from region_validation/scripts, see
# region_validation_d100/RUNBOOK_d100.md). Differences from the depth-20 script:
#   * RV / _common.BASE point at THIS tree; tmux sessions are suffixed _d100_legacy.
#   * alignment: --prefilter-k 500 --k2 100 --gpu-cand 500 --workers ${ALIGN_WORKERS}
#     --min-ref-len 200 --force-backend cuda --save-shortlist   (record depth 100; the
#     whole 500-deep shortlist is Biopython-rescored, and the shortlist + GPU scores are
#     persisted so a re-run or another depth skips the edlib stage via --shortlist-in).
#   * SELECT_MODE (exact_tie | legacy) is FIXED per tree (this tree: legacy)
#     and threaded to prefetch_gene_families / select_references / prefetch_universe /
#     close_exclude_rep_gap / score_taxacc / score_genecap.
#   * kmeans_panel/close_exclude_rep_gap.py runs after the universe prefetch.
#   * FRONT-LOAD is check-only: the seven panel files are COPIES of the depth-20 panel;
#     a missing one aborts (it is never rebuilt here -- select_10k.py would silently
#     build a different panel).
#   * --force is refused (it would rebuild/re-fetch shared state); delete the specific
#     output file to re-run one stage.
#
# Environment knobs:
#   STAGE=all|align|post   all (default) = alignment then POST; align = alignment ONLY
#                          (never touches the shared PGFam cache); post = POST only
#                          (requires align.done).
#   NO_TMUX=1              run every long job in the FOREGROUND (for systemd-run units)
#                          instead of detached tmux sessions; also skips the parallel
#                          truth-pass prefetch (the POST prefetch_pgfams call runs the truth
#                          pass AND the reps pass).
#   ALIGN_WORKERS=N        CPU workers for align_hits.py (default 40).
#   SELECT_MODE=...        optional; must equal this tree's rule (legacy).
#   ALLOW_GPU_SHARE=1      skip the "another align_hits/run_prfba_fold is running" guard.
#
# Usage:
#   STAGE=align NO_TMUX=1 bash run_validation.sh full     # alignment only, foreground
#   STAGE=post  NO_TMUX=1 bash run_validation.sh full     # POST only, foreground
#
# Stage order (DESIGN §9):
#   FRONT-LOAD (check only in this tree)
#   ALIGN   align_hits.py (GPU + CPU workers, --force-backend cuda)
#           [|| prefetch_pgfams.py truth pass, tmux mode only]
#   POST (sequential -- no GPU):
#     fan hits md5->{genome_id}__{region} via amplicon_md5_fanout
#       -> prefetch_gene_families.py --select-mode   (stage-7 candidates)
#       -> select_references.py --gene-provider auto --select-mode
#       -> prefetch_pgfams.py   (truth + reps pass)
#       -> prefetch_universe.py --select-mode
#       -> close_exclude_rep_gap.py --tree <tree> [--legacy]
#       -> score_taxacc.py [--legacy]
#       -> score_genecap.py [--legacy]
#       -> mcnemar_regions.py -> build_region_meta.py -> render_accuracy.py -> render_genecap.py
#
# NOTE on the fan-out (DESIGN §5.0 + §9):  dedup keys the unique amplicon FASTA on
# the *insert md5*, so align_hits.py emits an md5-keyed hits JSON.  Every POST
# scorer keys on the provenance id "{genome_id}__{region}".  The POST stage expands
# the md5-keyed raw hits back out to every source (genome,region) via
# amplicon_md5_fanout.json BEFORE selection.  The file keeps its depth-20 name
# (asv_top20_alignment_hits.json) and list key ("top20") by design: >50 consumers read
# them; the depth is recorded in hits/raw/method_selection.json params.k2.
# =============================================================================

set -euo pipefail

# --- interpreter (sole uv-managed venv; per global rules) -------------------
PY="/home/freiburger/Documents/py_venv/bin/python"

# --- repos ------------------------------------------------------------------
PRFBA="/home/freiburger/Documents/prFBA"                       # reused modules
RV="/home/freiburger/Documents/prFBA/region_validation_d100_legacy"   # depth-100 tree (original panel, legacy)
SCRIPTS="${RV}/scripts"                                         # stage scripts

# --- reused prFBA modules (run as-is, DESIGN §9) ----------------------------
ALIGN_HITS="${PRFBA}/align_hits.py"
SELECT_REFERENCES="${PRFBA}/select_references.py"
CLOSE_GAP="${PRFBA}/kmeans_panel/close_exclude_rep_gap.py"

# --- depth-100 tree identity (fixed) ----------------------------------------
TREE_SELECT_MODE="legacy"     # exact_tie (union trees) | legacy (*_legacy trees)
ALIGN_ALLOWED=0             # 1 = union tree aligns; 0 = legacy tree copies hits
UNION_TREE="/home/freiburger/Documents/prFBA/region_validation_d100"                 # where a legacy tree's hits come from
K2=100; GPU_CAND=500; PREFILTER_K=500

# =============================================================================
# Parse args: MODE in {pilot,full}; --force refused in this tree
# =============================================================================
MODE="${1:-}"
FORCE=0
shift || true
for a in "$@"; do
  case "$a" in
    --force)
      echo "run_validation.sh: --force is refused in the depth-100 trees (it would rebuild" >&2
      echo "  the copied panel and re-fetch shared PGFam state). Delete the one output you" >&2
      echo "  want rebuilt and re-run instead." >&2
      exit 2 ;;
    *) echo "run_validation.sh: unknown arg '$a'" >&2; exit 2 ;;
  esac
done

case "$MODE" in
  pilot) N=500 ;;
  full)  N=10000 ;;
  *)
    echo "Usage: [STAGE=all|align|post] [NO_TMUX=1] run_validation.sh {pilot|full}" >&2
    exit 2
    ;;
esac
FORCE_FLAG=""

STAGE="${STAGE:-all}"
NO_TMUX="${NO_TMUX:-0}"
ALIGN_WORKERS="${ALIGN_WORKERS:-40}"
SELECT_MODE="${SELECT_MODE:-${TREE_SELECT_MODE}}"
case "${STAGE}" in all|align|post) ;; *) echo "STAGE must be all|align|post (got '${STAGE}')" >&2; exit 2 ;; esac
case "${SELECT_MODE}" in exact_tie|legacy) ;; *) echo "SELECT_MODE must be exact_tie|legacy (got '${SELECT_MODE}')" >&2; exit 2 ;; esac
if [ "${SELECT_MODE}" != "${TREE_SELECT_MODE}" ]; then
  echo "[FAIL] ${RV} is the ${TREE_SELECT_MODE} tree; SELECT_MODE=${SELECT_MODE} would mix two" >&2
  echo "       selection rules in one data/ dir. Use the matching *_d100 / *_d100_legacy tree." >&2
  exit 2
fi
if [ "${SELECT_MODE}" = "legacy" ]; then
  SR_MODE="legacy"; PF_MODE="legacy"; LEGACY_FLAG="--legacy"
else
  SR_MODE="exact-tie"; PF_MODE="exact-tie"; LEGACY_FLAG=""
fi

# --- fd limit: the random panel's depth-20 reps pass died on Errno 24 ---------
ulimit -n 1048576 2>/dev/null || ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
if [ "$(ulimit -n)" != "unlimited" ] && [ "$(ulimit -n)" -lt 65536 ]; then
  echo "[WARN] open-file limit is $(ulimit -n) (< 65536); launch under systemd-run -p LimitNOFILE=1048576" >&2
fi

# =============================================================================
# Resolve ALL absolute paths from _common.PATHS (single source of truth).
# =============================================================================
eval "$(
  PYTHONPATH="${SCRIPTS}" "${PY}" - <<'PYEOF'
from _common import BASE, DATA, FIGS, PATHS  # noqa: E402
emit = {
    "COMMON_BASE": BASE,
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

# Hard isolation check: _common.BASE must be this tree, and DATA must sit under it.
if [ "${COMMON_BASE}" != "${RV}" ] || [ "${DATA}" != "${RV}/data" ]; then
  echo "[FAIL] _common.BASE=${COMMON_BASE} DATA=${DATA} -- not this tree (${RV}); refusing to run" >&2
  exit 2
fi

# Raw md5-keyed hits land in a subdir; we fan them out into HITS_JSON for POST.
RAW_HITS_DIR="${HITS_DIR}/raw"
RAW_HITS_JSON="${RAW_HITS_DIR}/asv_top20_alignment_hits.json"
SHORTLIST_NPZ="${RAW_HITS_DIR}/edlib_top${PREFILTER_K}_shortlist.npz"
SHORTLIST_META="${RAW_HITS_DIR}/shortlist_meta.json"

# sentinels + logs (under DATA/hits so they are resumability checkpoints too)
ALIGN_DONE="${RAW_HITS_DIR}/align.done"
ALIGN_LOG="${RAW_HITS_DIR}/align.log"
ALIGN_TMUX="rv_align_${MODE}_d100_legacy"
PGFAM_DONE="${DATA}/prefetch_pgfams.truth.done"
PGFAM_LOG="${DATA}/prefetch_pgfams.log"
PGFAM_TMUX="rv_pgfam_${MODE}_d100_legacy"
GAP_DONE="${DATA}/close_exclude_rep_gap.done"
GAP_LOG="${DATA}/close_exclude_rep_gap.log"

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
# Skips the command if OUTPUT_PATH already exists (resumable).
take_cache_lock() {
  # one writer of the shared PGFam cache at a time (load -> mutate -> whole-file replace)
  [ "${POST_LOCK_HELD:-0}" = "1" ] && return 0
  [ -n "${CACHE_LOCK_TAKEN:-}" ] && return 0
  exec 9>/home/freiburger/Documents/prFBA/region_validation_d100/.post_cache_writer.lock
  flock -n 9 || { echo "[FAIL] another cache writer holds the lock" >&2; exit 3; }
  CACHE_LOCK_TAKEN=1
}

run_stage() {
  local label="$1"; shift
  local out="$1"; shift
  if done_file "${out}"; then
    echo "[skip] ${label}: checkpoint present -> ${out}"
    return 0
  fi
  echo "[run ] ${label}  ($(date -Is))"
  "$@"
  echo "[ok  ] ${label} -> ${out}  ($(date -Is))"
}

# check_alignment: method_selection.json must say k2=100 / gpu_cand=500 / prefilter 500
check_alignment() {
  "${PY}" - "${RAW_HITS_DIR}/method_selection.json" "${K2}" "${GPU_CAND}" "${PREFILTER_K}" <<'PYEOF'
import json, sys
p, k2, gc, pk = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
m = json.load(open(p)); prm = m["params"]
bad = [(k, prm.get(k), v) for k, v in (("k2", k2), ("gpu_cand", gc), ("prefilter_k", pk)) if prm.get(k) != v]
if bad or "run_stats" not in m:
    sys.exit(f"[FAIL] {p}: params mismatch {bad} or no run_stats (alignment incomplete)")
print(f"[ok  ] alignment params k2={prm['k2']} gpu_cand={prm['gpu_cand']} prefilter_k={prm['prefilter_k']} "
      f"wall={m.get('wall_seconds_total')}s n_asvs={m.get('n_asvs')}")
PYEOF
}

echo ""
echo "#############################################################################"
echo "#  16S region validation  |  DEPTH-100 tree ${RV}"
echo "#  MODE=${MODE}  N=${N}  STAGE=${STAGE}  NO_TMUX=${NO_TMUX}  SELECT_MODE=${SELECT_MODE}"
echo "#  alignment   : --prefilter-k ${PREFILTER_K} --k2 ${K2} --gpu-cand ${GPU_CAND} --workers ${ALIGN_WORKERS}"
echo "#  interpreter : ${PY}"
echo "#  data dir    : ${DATA}"
echo "#  figures dir : ${FIGS}"
echo "#  fd limit    : $(ulimit -n)"
echo "#  started     : $(date -Is)"
echo "#############################################################################"

cd "${SCRIPTS}"

# =============================================================================
# FRONT-LOAD (check only) -- the seven panel files are copies of the depth-20 panel
# =============================================================================
banner "FRONT-LOAD  (check only: the copied panel must be complete)"
missing=0
for f in "${SELECTION_JSON}" "${SRC16S_FASTA}" "${COMBINED_FASTA}" "${COMBINED_STATS}" \
         "${UNIQUE_FASTA}" "${AMPLICON_FANOUT}" "${TRUTH_JSON}"; do
  if [ ! -s "$f" ]; then echo "[FAIL] front-load file missing: $f" >&2; missing=1; fi
  if [ -L "$f" ]; then echo "[FAIL] front-load file is a symlink (must be a copy): $f" >&2; missing=1; fi
done
if [ "${missing}" -ne 0 ]; then
  echo "[FAIL] copy the seven front-load files from region_validation/data (cp, never ln -s)." >&2
  exit 2
fi
echo "[ok  ] front-load: 7 panel files present (copies)"

# =============================================================================
# ALIGN
# =============================================================================
if [ "${STAGE}" != "post" ]; then
  banner "ALIGN  (align_hits.py  k2=${K2}, $( [ "${NO_TMUX}" = 1 ] && echo foreground || echo tmux ))"
  if done_file "${ALIGN_DONE}" && done_file "${RAW_HITS_JSON}"; then
    echo "[skip] align_hits.py: sentinel + raw hits present -> ${RAW_HITS_JSON}"
  elif [ "${ALIGN_ALLOWED}" -ne 1 ]; then
    echo "[FAIL] ${RV} is a legacy-rule tree: it never aligns. Copy the hits from" >&2
    echo "       ${UNION_TREE}/data/hits/raw (asv_top20_alignment_hits.json," >&2
    echo "       asv_alignment_summary.csv, method_selection.json, align.log, align.done) with cp." >&2
    exit 2
  elif [ "${NO_TMUX}" != "1" ] && tmux has-session -t "${ALIGN_TMUX}" 2>/dev/null; then
    echo "[wait] align_hits.py: tmux session '${ALIGN_TMUX}' already running -> joined below"
  else
    if [ "${ALLOW_GPU_SHARE:-0}" != "1" ]; then
      others="$(pgrep -af 'align_hits\.py|run_prfba_fold\.py' | awk '$2 ~ /(^|\/)python[0-9.]*$/' || true)"
      if [ -n "${others}" ]; then
        echo "[FAIL] another GPU/CPU-heavy alignment is running (set ALLOW_GPU_SHARE=1 to override):" >&2
        echo "${others}" >&2
        exit 3
      fi
    fi
    ALIGN_FLAGS=(--prefilter-k "${PREFILTER_K}" --k2 "${K2}" --gpu-cand "${GPU_CAND}"
                 --workers "${ALIGN_WORKERS}" --min-ref-len 200 --force-backend cuda --save-shortlist)
    RESUME_FLAGS=()
    if [ -s "${SHORTLIST_NPZ}" ] && [ -s "${SHORTLIST_META}" ]; then
      RESUME_FLAGS=(--shortlist-in "${RAW_HITS_DIR}")
      echo "[info] saved edlib shortlist found -> resuming with --shortlist-in (skips the edlib stage)"
    fi
    rm -f "${ALIGN_DONE}"
    if [ "${NO_TMUX}" = "1" ]; then
      echo "[run ] align_hits.py in the foreground ($(date -Is)) -> ${ALIGN_LOG}"
      set -o pipefail
      "${PY}" "${ALIGN_HITS}" --fasta "${UNIQUE_FASTA}" --outdir "${RAW_HITS_DIR}" \
        "${ALIGN_FLAGS[@]}" "${RESUME_FLAGS[@]}" 2>&1 | tee -a "${ALIGN_LOG}"
      check_alignment
      echo done > "${ALIGN_DONE}"
      echo "[ok  ] align_hits.py ($(date -Is))"
    else
      echo "[run ] align_hits.py in tmux '${ALIGN_TMUX}'"
      tmux new-session -d -s "${ALIGN_TMUX}" \
        "set -o pipefail; \
         '${PY}' '${ALIGN_HITS}' \
           --fasta   '${UNIQUE_FASTA}' \
           --outdir  '${RAW_HITS_DIR}' \
           ${ALIGN_FLAGS[*]} ${RESUME_FLAGS[*]} \
           2>&1 | tee -a '${ALIGN_LOG}' \
         && echo done > '${ALIGN_DONE}'"
    fi
  fi

  # --- prefetch_pgfams.py truth pass alongside alignment: tmux, STAGE=all only ----
  #     (STAGE=align must never write the shared cache; in NO_TMUX mode the POST
  #     prefetch_pgfams call runs the truth pass and the reps pass in one go.)
  if [ "${STAGE}" = "all" ] && [ "${NO_TMUX}" != "1" ]; then
    take_cache_lock
    if done_file "${PGFAM_DONE}"; then
      echo "[skip] prefetch_pgfams.py (truth pass): sentinel present -> ${PGFAM_CACHE}"
    elif tmux has-session -t "${PGFAM_TMUX}" 2>/dev/null; then
      echo "[wait] prefetch_pgfams.py: tmux session '${PGFAM_TMUX}' already running -> joined below"
    else
      echo "[run ] prefetch_pgfams.py (truth pass) in tmux '${PGFAM_TMUX}'"
      tmux new-session -d -s "${PGFAM_TMUX}" \
        "set -o pipefail; \
         '${PY}' '${SCRIPTS}/prefetch_pgfams.py' \
           2>&1 | tee '${PGFAM_LOG}' \
         && echo done > '${PGFAM_DONE}'"
    fi
  fi

  # --- join (tmux mode only) ---------------------------------------------------
  if [ "${NO_TMUX}" != "1" ]; then
    echo "[wait] joining tmux jobs -- polling every 30s..."
    while tmux has-session -t "${ALIGN_TMUX}" 2>/dev/null \
       || tmux has-session -t "${PGFAM_TMUX}" 2>/dev/null; do
      sleep 30
    done
  fi
fi

# Fail loudly if the heavy job did not produce its raw hits output.
if ! done_file "${RAW_HITS_JSON}" || ! done_file "${ALIGN_DONE}"; then
  echo "[FAIL] no complete alignment (${RAW_HITS_JSON} + ${ALIGN_DONE}); see ${ALIGN_LOG}" >&2
  exit 1
fi
check_alignment
echo "[ok  ] alignment present: raw hits=${RAW_HITS_JSON}"

if [ "${STAGE}" = "align" ]; then
  banner "DONE  (STAGE=align; POST not run)"
  exit 0
fi

# =============================================================================
# POST (sequential -- no GPU). Writes the SHARED PGFam cache: never run two trees'
# POST stages at the same time (the cache is load -> mutate -> whole-file replace).
# =============================================================================
take_cache_lock
banner "POST  (fan-out | select [${SR_MODE}] | prefetch | gap close | score | meta | render)"
echo "Resumability: each stage below skips when its output file already exists."

# --- fan hits md5 -> {genome_id}__{region} via amplicon_md5_fanout ----------
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
#     Under the SAME reducer select_references will run (--select-mode), so the
#     selection below has 0 lazy misses (each miss rewrites the whole cache).
echo "[run ] prefetch selection candidates (stage-7, ${PF_MODE}) -> ${PGFAM_CACHE}  ($(date -Is))"
"${PY}" "${PRFBA}/prefetch_gene_families.py" --hits "${HITS_JSON}" --cache "${PGFAM_CACHE}" \
  --workers 16 --select-mode "${PF_MODE}"
echo "[ok  ] selection candidates prefetched"

# --- select_references.py: include-self selection (single scenario) ----------
run_stage "select_references.py (--gene-provider auto --select-mode ${SR_MODE})" "${SELECTION_OUT}" \
  bash -c '"$0" "$1" --hits "$2" --out "$3.tmp" --gene-provider auto --gene-cache "$4" \
             --select-mode "$5" && mv "$3.tmp" "$3"' \
    "${PY}" "${SELECT_REFERENCES}" "${HITS_JSON}" "${SELECTION_OUT}" "${PGFAM_CACHE}" "${SR_MODE}"

# --- prefetch_pgfams.py: truth pass + reps pass (distinct selected reps) -------
echo "[run ] prefetch_pgfams.py (truth + reps pass)  ($(date -Is))"
"${PY}" "${SCRIPTS}/prefetch_pgfams.py"
echo "[ok  ] truth + reps PGFam cache warmed -> ${PGFAM_CACHE}"

# --- prefetch the bounded rep universe (include + exclude-self candidates) -----
echo "[run ] prefetch_universe.py (--select-mode ${PF_MODE}) -> ${PGFAM_CACHE}  ($(date -Is))"
"${PY}" "${SCRIPTS}/prefetch_universe.py" --select-mode "${PF_MODE}"
echo "[ok  ] rep universe prefetched"

# --- close the exclude-mode rep gap (fixed point; k-means depth-20 lesson) -----
#     score_genecap's exclude-self RE-selection can name reps no prefetch enumerated
#     (they would score pgfam_missing). Iterated to a fixed point; never caches [].
run_stage "close_exclude_rep_gap.py (--tree $(basename "${RV}") ${LEGACY_FLAG})" "${GAP_DONE}" \
  bash -c "set -o pipefail; '${PY}' '${CLOSE_GAP}' --tree '$(basename "${RV}")' ${LEGACY_FLAG} 2>&1 | tee -a '${GAP_LOG}' && date -Is > '${GAP_DONE}'"

# --- score_taxacc.py: taxonomic-accuracy scoring -> summary + long tables ----
run_stage "score_taxacc.py ${LEGACY_FLAG}" "${TAXACC_SUMMARY_JSON}" \
  "${PY}" "${SCRIPTS}/score_taxacc.py" ${LEGACY_FLAG}

# --- score_genecap.py: gene-capture scoring -> summary + long tables ---------
run_stage "score_genecap.py ${LEGACY_FLAG}" "${GENECAP_SUMMARY_JSON}" \
  "${PY}" "${SCRIPTS}/score_genecap.py" ${LEGACY_FLAG}

# --- mcnemar_regions.py: paired McNemar across region pairs (DESIGN Section 6) -
run_stage "mcnemar_regions.py" "${DATA}/mcnemar_region_pairs.csv" \
  "${PY}" "${SCRIPTS}/mcnemar_regions.py"

# --- build_region_meta.py: region metadata (+ domain, extract_rate, stamps) --
run_stage "build_region_meta.py" "${REGION_META_CSV}" \
  "${PY}" "${SCRIPTS}/build_region_meta.py"

# --- render: figures (cheap, self-overwriting; always run so they refresh) ---
banner "RENDER  (accuracy + gene-capture figures)"
echo "[run ] render_accuracy.py"
"${PY}" "${SCRIPTS}/render_accuracy.py"
echo "[run ] render_genecap.py"
"${PY}" "${SCRIPTS}/render_genecap.py"

# =============================================================================
banner "DONE  (MODE=${MODE}, N=${N}, SELECT_MODE=${SELECT_MODE}, $(date -Is))"
echo "Tables : ${DATA}"
echo "Figures: ${FIGS}"
echo "Hits   : ${HITS_JSON}  (fanned; raw md5-keyed: ${RAW_HITS_JSON})"
echo "Re-run is fully resumable."
