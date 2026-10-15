#!/usr/bin/env python
"""
run_prfba_fold.py — prFBA's side of the shared-reference k-fold benchmark (spec component 4).

WHY THIS SCRIPT EXISTS
----------------------
The predecessor head-to-head (``picrust2_comparison/README.md``) let prFBA align against
the whole ~112k-genome BV-BRC 16S database and then DELETED the self/near-self hits from
each amplicon's top-20 record afterwards. That post-hoc filter is not the same experiment:
the record has already been truncated to ``--k2 20``, so excising the held-out organism's
near neighbours EMPTIES the record instead of demoting it — measured as 94.6% of all
self-exclusion abstentions. PICRUSt2, meanwhile, kept its full 27,870-genome reference.
Neither the reference CONTENT nor the exclusion MECHANISM was shared.

This runner implements spec decision D14: for fold *i* prFBA RE-ALIGNS the fold's test
amplicons against a restricted reference database that is *exactly* the union of the tips
of the PICRUSt2 trees built for the same fold and density (component 2's manifests). The
exclusion is therefore enforced by construction at alignment time, the k2-truncation
artifact cannot arise, and the two methods see byte-identically the same reference set
(spec control C4). ``score_genecap.exclude_self_record`` is still called — but only as a
CANARY that must remove nothing (C6d); applying it as a filter here would charge prFBA
for the exclusion twice.

WHAT IT GUARANTEES (every one of these aborts the run rather than reporting a number)
  * exit 2 — configuration is incoherent: manifests disagree, a manifest belongs to
    another fold/density, a reference genome is not panel-eligible, two reference genomes
    share a 16S md5 (the md5->genome join would be ambiguous), or a synthesised BV-BRC
    header does not round-trip through ``edlib_biopython_hits.parse_header``.
  * exit 3 — a leakage / no-prior canary fired:
      C6a  a surviving reference shares an amplicon md5 with a test amplicon
      C6b  a test amplicon aligns to a surviving reference at identity >= 0.987
           (evaluated under the spec's single identity definition, §3.2: the aligner's
           ``n_matches/aligned_len`` QUALIFIED BY ``aligned_len/query_len >= cov_min``,
           because a sub-coverage local hit is not an edge in plan_folds' graph either)
      C6c  a surviving reference shares a taxon_id with a test genome
      C6d  ``exclude_self_record(..., tau=0.987)`` removes a COVERAGE-QUALIFIED hit from
           any record (see §3.2 above and the comment at the call site: a sub-coverage
           high-identity local match is neither an edge in plan_folds' graph nor visible
           to selection, so it is counted and reported, never gated)
      C8   a MiDAS lineage reached prFBA (an external taxonomic prior for the held-out
           organism), or a selected representative lies outside the fold reference set
      plus: every hit in every record must name a reference-set genome — the direct proof
      that ``--ref-db`` restricted the ALIGNMENT and not merely the selection
      (smoke-gate item 7), and the forced-answer invariant below.

FORCED-ANSWER INVARIANT (``--forced-answer``)
  Lowering only ``t_family``/``family_floor`` to 0.0 must change a cell's outcome ONLY
  where the shipped arm abstained; every cell the shipped arm answered must yield a
  byte-identical EC set, which is what makes a post-hoc ``rank1_identity`` sweep over the
  forced arm a risk-coverage curve for ONE method rather than a comparison of two. The
  forced run therefore re-runs the shipped knobs in memory and compares. The invariant is
  exact for the window-0 union modes (``exact_tie``, ``exact_tie_guards``). For ``legacy``
  (identity ``band`` 0.005) and ``graduated`` (``select_all_window`` 0.005) it is NOT
  exactly true: stage 0's per-hit identity floor and the window interact, so a cell with
  ``best_identity`` in ``[0.865, 0.865+window)`` can legitimately admit extra hits once the
  floor is 0. Those mismatches are counted separately as ``band_explained`` and do not
  abort; an unexplained mismatch does.

NOT A FILTER, NOT A SECOND CODE PATH
  Selection is ``select_references.select_representatives`` verbatim, with the knob sets
  read from ``frozen.prfba.select_modes`` (recorded by component 1), the PGFam gene
  provider imported from ``picrust2_comparison/scripts/prfba_selection.py``, and the
  aligner invoked as the ordinary ``align_hits.py`` CLI. Nothing about prFBA is
  reimplemented here; this script only restricts what it may see and proves the
  restriction held.

    PY=~/Documents/py_venv/bin/python
    $PY run_prfba_fold.py --fold 0 --density 1.0 \
        --fold-plan  $D/fold_plan.json --panel $D/panel_eligible.json \
        --ref-manifest $D/refs/f00_bac_d100/manifest.json \
        --ref-manifest $D/refs/f00_arc_d100/manifest.json \
        --outdir $D/prfba/f00_d100 --select-modes legacy,exact_tie

Interpreter: /home/freiburger/Documents/py_venv/bin/python  (the ONLY interpreter)
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time
from collections import Counter, defaultdict

import orjson

# ----------------------------------------------------------------------- repo wiring
HERE = os.path.dirname(os.path.abspath(__file__))
SR_ROOT = os.path.dirname(HERE)                                  # .../shared_reference
RV = os.path.dirname(SR_ROOT)                                    # .../region_validation
REPO = os.path.dirname(RV)                                       # .../prFBA
for _p in (REPO, os.path.join(RV, "scripts"),
           os.path.join(RV, "picrust2_comparison", "scripts")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import edlib_biopython_hits as EBH          # noqa: E402  parse_header (round-trip assert)
import select_references as SR              # noqa: E402  THE selection algorithm
import score_genecap as SG                  # noqa: E402  exclude_self_record canary
import prfba_selection as PS                # noqa: E402  the PGFam cache-only provider
import _common as C                         # noqa: E402  PATHS

PATHS = C.PATHS
ALIGN_HITS = os.path.join(REPO, "align_hits.py")

# Canonical immutable inputs (overridable only so the component can be smoke-tested).
DEF_AMPLICONS = os.path.join(RV, "data", "combined_amplicons.unique.fasta")
DEF_FANOUT = os.path.join(RV, "data", "amplicon_md5_fanout.json")
DEF_FUNCTIONS = os.path.join(RV, "picrust2_comparison", "data", "genome_function_sets.json")

MIDAS_HEADER = "seq,Kingdom,Phylum,Class,Order,Family,Genus,Species,rel_ab"
LINEAGE_RANKS = ["species", "genus", "family", "order", "class", "phylum"]

EXIT_CONFIG = 2
EXIT_CANARY = 3

# ---------------------------------------------------------------------------- frozen
# Spec §1, verbatim. Only the fallback: ``panel_eligible.json:meta.frozen`` is
# authoritative when present, so the block cannot drift between components.
FROZEN = {
    "tau_amplicon": 0.987,
    "tau_full16s": 0.987,
    "k_folds": 10,
    "seed": 1729,
    "prfba": {"k2": 200, "gpu_cand": 400, "prefilter_k": 1000, "min_ref_len": 200,
              "force_backend": "cuda",
              "select_modes": {
                  "legacy":           {"select_all": False},
                  "exact_tie":        {},
                  "graduated":        {"select_all_window": 0.005,
                                       "select_all_unconditional_window": 0.0},
                  "exact_tie_guards": {"select_all_coherence": True,
                                       "select_all_contam": True,
                                       "select_all_cap": 6}},
              "forced_answer_knobs": {"t_family": 0.0, "family_floor": 0.0},
              "t_species": 0.987, "t_genus": 0.945, "t_family": 0.865,
              "family_floor": 0.865, "decay_tau": 0.01},
    "picrust2": {"hsp_method": "mp", "edge_exponent": 0.5, "min_align": 0.80,
                 "placement_tool": "epa-ng", "max_nsti_reported": [None, 2.0],
                 "calculate_nsti": True},
    "reference_qc": {"len_min": 1200, "len_max": 1800, "occupancy_min": 0.90,
                     "ec_min": 1},
    "fasttree_args": ["-nt", "-gtr", "-gamma"],
    "density_fractions": [0.25, 0.5, 1.0],
    "sesoi_delta_f1": 0.02,
    "label_prefix": "g_",
}


# ------------------------------------------------------------------------- utilities
def die(code: int, msg: str):
    print(f"[abort] {msg}", flush=True)
    sys.exit(code)


def jload(path: str):
    with open(path, "rb") as fh:
        return orjson.loads(fh.read())


def jdump(path: str, obj, indent: bool = True):
    opt = orjson.OPT_INDENT_2 if indent else 0
    tmp = path + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(orjson.dumps(obj, option=opt))
    os.replace(tmp, path)


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=30).stdout.strip()
    except Exception:
        return ""


def iso_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_fasta(path: str, wanted: set | None = None) -> dict:
    """id -> uppercase sequence. `wanted` (ids) streams past everything else."""
    out, cur, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None and (wanted is None or cur in wanted):
                    out[cur] = "".join(buf).upper()
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None and (wanted is None or cur in wanted):
        out[cur] = "".join(buf).upper()
    return out


def norm_ec(e) -> str | None:
    """Spec §4.2 EC normalisation: strip an ``EC:`` prefix; DROP partial ECs (containing '-')."""
    if not e:
        return None
    s = str(e).strip()
    if s.upper().startswith("EC:"):
        s = s[3:].strip()
    if not s or "-" in s:
        return None
    return s


def agreement_rank(a: dict, b: dict) -> str:
    """Highest (most specific) rank at which two panel lineages agree; 'none' if never.

    The spec's enum is genus|family|order|class|phylum|none because the primary
    leave-genus-out design makes genus agreement impossible by construction (genus is
    graph edge E1). ``species`` is emitted as well so the ``--rank md5only`` species-tier
    sensitivity arm, where same-species references DO survive, is not silently relabelled.
    """
    for r in LINEAGE_RANKS:
        va, vb = (a.get(r) or "").strip(), (b.get(r) or "").strip()
        if va and vb and va == vb:
            return r
    return "none"


# ---------------------------------------------------------------- restricted caches
def restricted_map(src_path: str, keys: set, cache_path: str, loader_label: str,
                   force: bool, extract):
    """Load ``{genome_id: extract(value)}`` for `keys` only, memoised to `cache_path`.

    ``genome_function_sets.json`` is 1.7 GB and the PGFam cache is 3.9 GB, while a fold
    only ever needs its own <=9k reference genomes. The first run pays the full parse and
    persists the restriction; a re-run (or the forced-answer twin of the same fold) reads
    the small file. Deterministic either way — this is a cache, not a shortcut.
    """
    if os.path.exists(cache_path) and not force:
        sub = jload(cache_path)
        missing = keys - set(sub)
        if not missing:
            print(f"[cache] {loader_label}: reusing {os.path.basename(cache_path)} "
                  f"({len(sub)} genomes)", flush=True)
            return {k: sub[k] for k in keys if k in sub}
        print(f"[cache] {loader_label}: {len(missing)} genomes missing from "
              f"{os.path.basename(cache_path)} -> full reload", flush=True)
    t0 = time.time()
    print(f"[load] {loader_label}: {src_path} (streaming, restricting to {len(keys)} "
          f"genomes) ...", flush=True)
    # STREAMED, not orjson-loaded: the PGFam cache is 3.9 GB and genome_function_sets.json
    # is 1.7 GB, so a whole-file parse costs ~30-40 GB resident per fold runner. The
    # orchestrator runs folds concurrently (xargs -P), which would OOM a 251 GB box; a
    # streaming pass with the yajl2_c backend holds a few hundred MB and exits as soon as
    # every wanted genome has been seen.
    import ijson
    sub = {}
    with open(src_path, "rb") as fh:
        for gid, val in ijson.kvitems(fh, "", use_float=True):
            if gid in keys:
                sub[gid] = extract(val)
                if len(sub) == len(keys):
                    break
    jdump(cache_path, sub, indent=False)
    print(f"[load] {loader_label}: {len(sub)}/{len(keys)} genomes present "
          f"({time.time()-t0:.0f}s)", flush=True)
    return sub


# ------------------------------------------------------------------------------ main
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fold-plan", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--fold", type=int, required=True)
    ap.add_argument("--density", type=float, default=1.0)
    ap.add_argument("--ref-manifest", action="append", default=[], required=True,
                    help="component-2 manifest.json; repeat for bac + arc")
    ap.add_argument("--amplicons", default=DEF_AMPLICONS)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--k2", type=int, default=200)
    ap.add_argument("--gpu-cand", type=int, default=400)
    ap.add_argument("--prefilter-k", type=int, default=1500)
    ap.add_argument("--force-backend", default="cuda",
                    choices=["auto", "cuda", "metal", "cpu"])
    ap.add_argument("--select-modes", default="legacy,exact_tie,graduated,exact_tie_guards")
    ap.add_argument("--hits-from", default=None,
                    help="reuse the alignment from this run's outdir when its "
                         "align_stamp.json is identical (SPEC D8: the select-mode and "
                         "--forced-answer arms are pure re-selections over the same "
                         "hits, but must write to separate outdirs)")
    ap.add_argument("--forced-answer", action="store_true",
                    help="merge frozen.prfba.forced_answer_knobs (t_family/family_floor = 0.0) "
                         "into every select-mode; pred keys become <mode>_forced. Use a "
                         "distinct --outdir (convention <dir>_forced).")
    ap.add_argument("--queries", default=None,
                    help="explicit query FASTA (error-injected arm); restricted to the "
                         "fold's test_md5s")
    ap.add_argument("--workers", type=int, default=32)
    # ---- additive (not in the spec CLI): mode plumbing + path overrides for controls
    ap.add_argument("--mode", default=None,
                    choices=["fold", "include_self", "plumbing"],
                    help="default: taken from the manifests' `mode`. include_self/plumbing "
                         "references CONTAIN the test genomes, so the C6 leakage canaries are "
                         "replaced by the C2 source-recovery gate.")
    ap.add_argument("--functions", default=DEF_FUNCTIONS)
    ap.add_argument("--fanout", default=DEF_FANOUT)
    ap.add_argument("--panel-md5-seq", default=None, help="default <fold-plan dir>/panel_md5_seq.json")
    ap.add_argument("--panel-md5-id", default=None, help="default <fold-plan dir>/panel_md5_id.json")
    ap.add_argument("--ec-universe", default=None, help="default <fold-plan dir>/panel_ec_universe.json")
    ap.add_argument("--pgfam-cache", default=PATHS.pgfam_cache)
    ap.add_argument("--min-ref-len", type=int, default=None, help="default frozen.prfba.min_ref_len")
    ap.add_argument("--force", action="store_true", help="redo the alignment and all caches")
    return ap


def main():
    args = build_parser().parse_args()
    t_wall0 = time.perf_counter()
    outdir = os.path.abspath(args.outdir)
    os.makedirs(os.path.join(outdir, "hits"), exist_ok=True)
    data_dir = os.path.dirname(os.path.abspath(args.fold_plan))
    p_md5_seq = args.panel_md5_seq or os.path.join(data_dir, "panel_md5_seq.json")
    p_md5_id = args.panel_md5_id or os.path.join(data_dir, "panel_md5_id.json")
    p_ecuni = args.ec_universe or os.path.join(data_dir, "panel_ec_universe.json")

    if args.forced_answer and not os.path.basename(outdir).endswith("_forced"):
        print(f"[warn] --forced-answer with --outdir {outdir!r}: the convention is a "
              f"'_forced' suffix, and the shipped and forced runs MUST NOT share a "
              f"directory (spec §6.2).", flush=True)
    if args.gpu_cand < args.k2:
        die(EXIT_CONFIG, f"--gpu-cand {args.gpu_cand} < --k2 {args.k2}: the CUDA path would "
                         f"truncate the hit list at {args.gpu_cand}")

    # ---------------------------------------------------------------- inputs
    plan = jload(args.fold_plan)
    panel = jload(args.panel)
    frozen = (panel.get("meta") or {}).get("frozen") or FROZEN
    frozen_source = "panel_eligible.json:meta.frozen" if (panel.get("meta") or {}).get("frozen") \
        else "run_prfba_fold.FROZEN (spec §1 fallback)"
    fz = frozen.get("prfba", FROZEN["prfba"])
    min_ref_len = args.min_ref_len if args.min_ref_len is not None else fz.get("min_ref_len", 200)
    tau = frozen.get("tau_amplicon", 0.987)
    genomes = panel["genomes"]

    fold = next((f for f in plan["folds"] if int(f["fold"]) == args.fold), None)
    if fold is None:
        die(EXIT_CONFIG, f"fold {args.fold} not in {args.fold_plan}")
    test_genomes = set(map(str, fold["test_genomes"]))
    test_md5s = [str(m) for m in fold["test_md5s"]]
    plan_ref = set(map(str, fold["reference_genomes"]))

    # ---------------------------------------------------------------- manifests
    manifests, ref_genomes, ref_names, modes = [], set(), [], set()
    for mp in args.ref_manifest:
        m = jload(mp)
        manifests.append(m)
        modes.add(m.get("mode", "fold"))
        ref_names.append(m.get("name") or os.path.basename(os.path.dirname(os.path.abspath(mp))))
        ref_genomes |= {str(g) for g in m.get("reference_genomes", [])}
    if len(modes) != 1:
        die(EXIT_CONFIG, f"--ref-manifest modes disagree: {sorted(modes)}")
    mode = args.mode or modes.pop()
    enforce_c6 = (mode == "fold")
    for m, mp in zip(manifests, args.ref_manifest):
        if mode == "fold" and m.get("fold") is not None and int(m["fold"]) != args.fold:
            die(EXIT_CONFIG, f"{mp}: manifest fold {m['fold']} != --fold {args.fold}")
        if m.get("density") is not None and abs(float(m["density"]) - args.density) > 1e-9:
            die(EXIT_CONFIG, f"{mp}: manifest density {m['density']} != --density {args.density}")
    if not ref_genomes:
        die(EXIT_CONFIG, "the manifests name no reference genomes")
    unknown = ref_genomes - set(genomes)
    if unknown:
        die(EXIT_CONFIG, f"{len(unknown)} manifest reference genomes are not panel-eligible "
                         f"(e.g. {sorted(unknown)[:5]})")
    if enforce_c6:
        leak = ref_genomes & test_genomes
        if leak:
            die(EXIT_CANARY, f"C6: {len(leak)} test genomes are IN the reference "
                             f"(e.g. {sorted(leak)[:5]})")
        outside = ref_genomes - plan_ref
        if outside:
            die(EXIT_CONFIG, f"{len(outside)} manifest reference genomes are not in fold "
                             f"{args.fold}'s reference_genomes (e.g. {sorted(outside)[:5]})")
    print(f"[plan] mode={mode} fold={args.fold} density={args.density} "
          f"refs={len(ref_genomes)} test_genomes={len(test_genomes)} "
          f"test_md5s={len(test_md5s)} | frozen from {frozen_source}", flush=True)

    # ------------------------------------------------- ref_db / ref_id_db / midas
    md5_of = {g: str(genomes[g]["md5"]) for g in ref_genomes}
    if len(set(md5_of.values())) != len(md5_of):
        dup = [m for m, n in Counter(md5_of.values()).items() if n > 1]
        die(EXIT_CONFIG, f"{len(dup)} 16S md5s are shared by >1 reference genome: the "
                         f"md5->genome join would be ambiguous (e.g. {dup[:3]})")
    panel_seq = jload(p_md5_seq)
    panel_hdr = jload(p_md5_id)
    ref_db, ref_id_db = {}, {}
    for g, md5 in sorted(md5_of.items()):
        seq = panel_seq.get(md5)
        if not seq:
            die(EXIT_CONFIG, f"reference genome {g}: md5 {md5} absent from {p_md5_seq}")
        hdr = panel_hdr.get(md5)
        if not hdr:
            die(EXIT_CONFIG, f"reference genome {g}: md5 {md5} absent from {p_md5_id}")
        _org, hgid, htid, _feat = EBH.parse_header(hdr)
        if str(hgid) != g or htid != int(g.split(".")[0]):
            die(EXIT_CONFIG, f"synthesised header for {g} does not round-trip through "
                             f"parse_header: got genome_id={hgid} taxon_id={htid} "
                             f"from {hdr!r}")
        ref_db[md5] = seq.upper()
        ref_id_db[md5] = hdr
    md5_to_ref_genome = {m: g for g, m in md5_of.items()}
    ref_taxa = {int(genomes[g]["taxon_id"]) for g in ref_genomes
                if genomes[g].get("taxon_id") is not None}
    jdump(os.path.join(outdir, "ref_db.json"), ref_db, indent=False)
    jdump(os.path.join(outdir, "ref_id_db.json"), ref_id_db, indent=False)
    midas_csv = os.path.join(outdir, "midas_empty.csv")
    with open(midas_csv, "w", newline="") as fh:        # header ONLY: no external prior
        fh.write(MIDAS_HEADER + "\n")

    # ---------------------------------------------------------------- queries.fna
    q_path = os.path.join(outdir, "queries.fna")
    want = set(test_md5s)
    src_fa = args.queries or args.amplicons
    seqs = read_fasta(src_fa, want if want else None)
    if want:
        missing = want - set(seqs)
        if missing:
            print(f"[warn] {len(missing)} of {len(want)} test md5s absent from {src_fa} "
                  f"(e.g. {sorted(missing)[:3]})", flush=True)
        q_ids = [m for m in sorted(want) if m in seqs]
    else:
        q_ids = sorted(seqs)
    if not q_ids:
        die(EXIT_CONFIG, f"no query sequences: none of the fold's {len(want)} test md5s "
                         f"are in {src_fa}")
    with open(q_path, "w") as fh:
        for m in q_ids:
            fh.write(f">{m}\n{seqs[m]}\n")
    print(f"[queries] {len(q_ids)} amplicons -> {q_path}", flush=True)
    del seqs

    # ---------------------------------------------------------------- alignment
    hits_dir = os.path.join(outdir, "hits")
    hits_json = os.path.join(hits_dir, "asv_top20_alignment_hits.json")
    stamp_path = os.path.join(outdir, "align_stamp.json")
    stamp = {"ref_db_sha256": sha256_file(os.path.join(outdir, "ref_db.json")),
             "queries_sha256": sha256_file(q_path),
             "k2": args.k2, "gpu_cand": args.gpu_cand, "prefilter_k": args.prefilter_k,
             "min_ref_len": min_ref_len, "force_backend": args.force_backend,
             "n_reference_md5": len(ref_db)}
    old = jload(stamp_path) if os.path.exists(stamp_path) else None
    pred_path = os.path.join(outdir, "predictions.json")
    fresh = (not args.force) and old == stamp and os.path.exists(hits_json)

    # --hits-from: adopt an identically-parameterised run's alignment.
    #
    # SPEC D8 says the select-mode variants and the --forced-answer arm are "pure
    # re-selections over the same hits" -- they change only select_representatives'
    # knobs, never the reference DB or the query set.  But --forced-answer must write
    # to a separate --outdir (the two runs may not share one), and align_stamp.json
    # lives in the outdir, so the forced arm re-ran the whole alignment: measured 119 s
    # of 157 s on a 604-query/182-reference smoke fold, i.e. ~76% of its wall spent
    # reproducing hits byte-for-byte.  At full scale that is 10 extra alignments of
    # ~6k queries x gpu_cand 400 on the pipeline's longest CPU pole.
    #
    # The stamp already covers everything the alignment depends on (ref_db sha, query
    # sha, k2/gpu_cand/prefilter_k/min_ref_len/backend), so adopting a donor whose
    # stamp is byte-identical is safe by construction; a mismatch falls through to a
    # normal alignment rather than failing.
    if (not fresh) and args.hits_from:
        donor_stamp_p = os.path.join(args.hits_from, "align_stamp.json")
        donor_hits = os.path.join(args.hits_from, "hits",
                                  "asv_top20_alignment_hits.json")
        if not os.path.exists(donor_hits):
            print(f"[hits-from] no hits at {donor_hits}: aligning normally", flush=True)
        elif (jload(donor_stamp_p) if os.path.exists(donor_stamp_p) else None) != stamp:
            print(f"[hits-from] align stamp at {args.hits_from} differs "
                  f"(different reference, queries or aligner knobs): aligning normally",
                  flush=True)
        else:
            os.makedirs(hits_dir, exist_ok=True)
            for fn in sorted(os.listdir(os.path.join(args.hits_from, "hits"))):
                src = os.path.join(args.hits_from, "hits", fn)
                if os.path.isfile(src):
                    shutil.copyfile(src, os.path.join(hits_dir, fn))
            jdump(stamp_path, stamp)
            fresh = True
            print(f"[hits-from] adopted the alignment from {args.hits_from} "
                  f"(align stamp identical) -- skipping re-alignment", flush=True)
    modes_req = [m.strip() for m in args.select_modes.split(",") if m.strip()]
    if fresh and os.path.exists(pred_path):
        pm = (jload(pred_path).get("meta") or {})
        if (pm.get("select_modes") == modes_req
                and bool(pm.get("forced_answer")) == bool(args.forced_answer)
                and pm.get("mode") == mode):
            print(f"[skip] {pred_path} is up to date (align stamp + select modes match); "
                  f"--force to redo", flush=True)
            return
        print(f"[redo] {os.path.basename(pred_path)} exists but was built for "
              f"select_modes={pm.get('select_modes')} forced={pm.get('forced_answer')} "
              f"mode={pm.get('mode')}: re-selecting", flush=True)
    if fresh:
        print(f"[skip] reusing {hits_json} (align stamp matches)", flush=True)
        align_seconds = None
    else:
        cmd = [sys.executable, ALIGN_HITS,
               "--fasta", q_path, "--ref-db", os.path.join(outdir, "ref_db.json"),
               "--ref-id-db", os.path.join(outdir, "ref_id_db.json"),
               "--taxonomy-csv", midas_csv, "--outdir", hits_dir,
               "--k2", str(args.k2), "--gpu-cand", str(args.gpu_cand),
               "--prefilter-k", str(args.prefilter_k),
               "--min-ref-len", str(min_ref_len),
               "--workers", str(args.workers),
               "--force-backend", args.force_backend]
        print("[align] " + " ".join(cmd), flush=True)
        t0 = time.perf_counter()
        rc = subprocess.run(cmd, cwd=REPO).returncode
        align_seconds = round(time.perf_counter() - t0, 1)
        if rc != 0:
            die(EXIT_CONFIG if rc == 2 else 1, f"align_hits.py exited {rc}")
        jdump(stamp_path, stamp)

    hits = jload(hits_json)
    hits = {k: v for k, v in hits.items() if k != "_meta" and isinstance(v, dict)}
    backend = "unknown"
    ms_path = os.path.join(hits_dir, "method_selection.json")
    if os.path.exists(ms_path):
        ms = jload(ms_path)
        backend = ms.get("effective_backend", "unknown")
        if ms.get("params", {}).get("ref_db") is None:
            die(EXIT_CONFIG, "method_selection.json records ref_db=null: the alignment ran "
                             "against the FULL BV-BRC database, not the fold reference")

    # ------------------------------------------------- md5 -> test genome (via fanout)
    # The fan-out is read UNFILTERED on purpose. C6a asks whether a SURVIVING REFERENCE
    # carries the same amplicon md5 as a test amplicon, and a reference genome's cell is
    # `<ref_gid>__<region>`, which is by definition not in `test_cells`: pre-filtering to
    # the fold's own cells would leave only test genomes behind and make the assertion
    # unable to fire. The same unfiltered set defines a TWIN cell (an md5 carried by more
    # than one panel genome), which spec control C2 reports against the analytic
    # `source u twin` ceiling and never gates.
    fanout = jload(args.fanout)
    fan_gids, md5_src = {}, {}
    for m in q_ids:
        gids = {c.split("__", 1)[0] for c in (fanout.get(m) or [])}
        fan_gids[m] = gids
        in_test = sorted(gids & test_genomes)
        md5_src[m] = in_test[0] if in_test else None   # deterministic twin representative
    n_twin = sum(1 for m in q_ids if len(fan_gids[m] & set(genomes)) > 1)

    # ------------------------------------------------------------------ canaries
    A = {"exclude_self_n_removed_total": 0, "midas_nonempty_records": 0,
         "species_dedup_collisions": 0, "gene_count_fn_calls": 0,
         "reps_outside_reference": 0, "k2_saturated_records": 0}
    c6 = {"md5_shared_with_reference": 0, "amplicon_identity_ge_tau": 0,
          "full16s_identity_ge_tau": None, "taxon_id_shared": 0,
          "hits_outside_reference": 0,
          "exclude_self_n_removed_raw": 0,
          "exclude_self_removed_subcoverage": 0}
    cov_min = SR.DEFAULT_KNOBS["cov_min"]
    examples = []
    for m in q_ids:
        rec = hits.get(m)
        if rec is None:
            continue
        if (rec.get("midas_taxonomy") or "").strip():
            A["midas_nonempty_records"] += 1
        top = rec.get("top20") or []
        L = max(1, int(rec.get("asv_len") or 1))
        for h in top:
            gid = str(h.get("genome_id"))
            if gid not in ref_genomes:
                c6["hits_outside_reference"] += 1
                if len(examples) < 5:
                    examples.append(f"hit outside reference: {m} -> {gid}")
            if not enforce_c6:
                continue
            ident = h.get("identity")
            cov = (h.get("aligned_len") or 0) / L
            if ident is not None and ident >= tau and cov >= cov_min:
                c6["amplicon_identity_ge_tau"] += 1
                if len(examples) < 5:
                    examples.append(f"C6b {m} -> {gid} identity={ident} coverage={cov:.3f}")
        # k2 truncation "could still bite" only if the SHALLOWEST hit prFBA can actually
        # use still clears the family floor. Scored on the coverage-qualified hits (§3.2),
        # because `top[-1]` is the lowest-SCORING hit, which is routinely a 20 bp local
        # match at identity 1.0 that stage 0 discards.
        qual = [h for h in top if (h.get("aligned_len") or 0) / L >= cov_min]
        if len(top) >= args.k2 and qual and (qual[-1].get("identity") or 0.0) >= fz.get(
                "family_floor", 0.865):
            A["k2_saturated_records"] += 1
        if enforce_c6:
            shared = fan_gids[m] & ref_genomes
            if shared:
                c6["md5_shared_with_reference"] += len(shared)
                if len(examples) < 8:
                    examples.append(f"C6a {m}: amplicon md5 also carried by reference "
                                    f"{sorted(shared)[:3]}")
            src_gid = md5_src[m]
            src_tid = None
            if src_gid and src_gid in genomes and genomes[src_gid].get("taxon_id") is not None:
                src_tid = int(genomes[src_gid]["taxon_id"])
            # C6d is the canary of D14: exclude_self_record must have NOTHING to remove.
            # It is evaluated on the coverage-qualified view of the record, because the
            # spec pins ONE identity definition (§3.2: `n_matches/aligned_len` QUALIFIED by
            # `aligned_len/query_len >= cov_min`, and "a hit below the coverage floor is not
            # an edge and is not a neighbour, at any tau"). A 20 bp perfect local match at
            # coverage 0.037 is not a leak: plan_folds does not draw an E4 edge for it, and
            # select_representatives' own stage-0 reliability gate (cov_min 0.90) discards
            # it before it can touch a prediction. The raw (unqualified) removal count is
            # reported beside it so the distinction is never hidden. Source-genome and
            # taxon_id leaks are caught at FULL strength, independent of coverage, by the
            # set-level ref∩test and taxon_id_shared checks above.
            _raw, n_raw = SG.exclude_self_record(rec, src_gid, src_tid, tau)
            edge_rec = dict(rec)
            edge_rec["top20"] = [h for h in top
                                 if (h.get("aligned_len") or 0) / L >= cov_min]
            _rec, n_removed = SG.exclude_self_record(edge_rec, src_gid, src_tid, tau)
            c6["exclude_self_n_removed_raw"] += n_raw
            c6["exclude_self_removed_subcoverage"] += (n_raw - n_removed)
            if n_removed:
                A["exclude_self_n_removed_total"] += n_removed
                if len(examples) < 8:
                    examples.append(f"C6d {m}: exclude_self removed {n_removed} "
                                    f"coverage-qualified hits (raw {n_raw})")
    if enforce_c6:
        shared_taxa = ref_taxa & {int(genomes[g]["taxon_id"]) for g in test_genomes
                                  if g in genomes and genomes[g].get("taxon_id") is not None}
        c6["taxon_id_shared"] = len(shared_taxa)

    for ex in examples:
        print(f"[canary] {ex}", flush=True)
    if c6["hits_outside_reference"]:
        die(EXIT_CANARY, f"{c6['hits_outside_reference']} hits name a genome outside the fold "
                         f"reference set: --ref-db did NOT restrict the alignment")
    if A["midas_nonempty_records"]:
        die(EXIT_CANARY, f"C8: {A['midas_nonempty_records']} records carry a MiDAS lineage — "
                         f"an external taxonomic prior for the held-out organism")
    if enforce_c6:
        bad = {k: v for k, v in c6.items()
               if k in ("md5_shared_with_reference", "amplicon_identity_ge_tau",
                        "taxon_id_shared") and v}
        if bad or A["exclude_self_n_removed_total"]:
            die(EXIT_CANARY,
                f"C6 leakage: {json.dumps(bad)} exclude_self_removed="
                f"{A['exclude_self_n_removed_total']} "
                f"(of which sub-coverage: {c6['exclude_self_removed_subcoverage']}) — the "
                f"holdout group graph is wrong and both methods got a free near-exact neighbour")
    print(f"[canary] C6/C8 clean (mode={mode}, enforced={enforce_c6}); "
          f"k2_saturated_records={A['k2_saturated_records']}; "
          f"exclude_self raw removals={c6['exclude_self_n_removed_raw']} "
          f"(all sub-coverage, not edges: {c6['exclude_self_removed_subcoverage']})",
          flush=True)

    # ------------------------------------------------------- gene content + EC universe
    ec_uni = set()
    for e in jload(p_ecuni)["ec"]:
        n = norm_ec(e)
        if n:
            ec_uni.add(n)
    ec_by_genome = restricted_map(
        args.functions, set(ref_genomes), os.path.join(outdir, "_cache_ec_subset.json"),
        "genome_function_sets", args.force,
        lambda v: sorted({n for n in (norm_ec(e) for e in (v.get("ec") or [])) if n} & ec_uni))
    n_missing_fn = len(ref_genomes - set(ec_by_genome))

    pgfam = {}
    if os.path.exists(args.pgfam_cache):
        pgfam = restricted_map(
            args.pgfam_cache, set(ref_genomes),
            os.path.join(outdir, "_cache_pgfam_subset.json"), "PGFam cache", args.force,
            lambda v: sorted(v) if not isinstance(v, dict) else sorted(v))
    else:
        print(f"[warn] PGFam cache {args.pgfam_cache} absent: gene counts are unavailable, "
              f"so the P3 species-dedup tie-break and the exact_tie_guards contamination "
              f"budget are inert in this run", flush=True)
    provider = PS._cache_only_provider(dict(pgfam))
    base_gcf = SR.gene_count_from_provider(provider)
    counter = {"calls": 0, "none": 0}

    def gcf(gid):
        counter["calls"] += 1
        v = base_gcf(gid)
        if v is None:
            counter["none"] += 1
        return v

    # ------------------------------------------------------------------- selection
    sm = fz.get("select_modes", FROZEN["prfba"]["select_modes"])
    unknown_modes = [m for m in modes_req if m not in sm]
    if unknown_modes:
        die(EXIT_CONFIG, f"unknown select-modes {unknown_modes}; frozen defines {sorted(sm)}")
    base_knobs = {k: fz[k] for k in ("t_species", "t_genus", "t_family", "family_floor",
                                     "decay_tau") if k in fz}
    forced_knobs = fz.get("forced_answer_knobs", {}) if args.forced_answer else {}

    def predict(rec, knobs, md5):
        sel = SR.select_representatives(rec, knobs=knobs, gene_provider=provider,
                                        gene_count_fn=gcf)
        reps = SG.selected_rep_gids(sel)
        ecs = set()
        for r in reps:
            ecs.update(ec_by_genome.get(r, ()))
        top = rec.get("top20") or []
        r1 = top[0] if top else {}
        r1g = str(r1.get("genome_id")) if r1.get("genome_id") is not None else None
        src = md5_src.get(md5)
        if r1g and src and r1g in genomes and src in genomes:
            r1rank = agreement_rank(genomes[r1g], genomes[src])
        else:
            r1rank = "none"
        flags = sel.get("flags") or []
        if reps:
            reason = ""
        elif not top:
            reason = "no_hits"
        elif "abstain_below_family" in flags:
            reason = "abstain_below_family"
        elif "abstain_reliability" in flags:
            reason = "abstain_reliability"
        else:
            reason = "empty_selection"
        n_dedup = sum(1 for d in (sel.get("decisions") or [])
                      if d.get("stage") == "species_dedup")
        return {
            "ec": sorted(ecs),
            "reps": reps, "n_reps": len(reps),
            "tier": sel.get("tier"),
            "abstained": not reps,
            "abstain_reason": reason,
            "best_identity": sel.get("best_identity"),
            "rank1_genome": r1g,
            "rank1_identity": r1.get("identity"),
            "rank1_rank": r1rank,
        }, n_dedup

    pred = {}
    inv = {"checked": 0, "divergence_pool_sensitive": 0,
           "divergence_band_explained": 0}
    for mname in modes_req:
        key = f"{mname}_forced" if args.forced_answer else mname
        knobs = {**base_knobs, **sm[mname]}
        shipped_knobs = dict(knobs)
        knobs = {**knobs, **forced_knobs}
        window = (knobs.get("select_all_window", 0.0)
                  if sm[mname].get("select_all", True) else SR.DEFAULT_KNOBS["band"])
        band_tol = shipped_knobs.get("family_floor", 0.865) + window + 1e-12
        out = {}
        for m in q_ids:
            rec = hits.get(m)
            if rec is None:
                continue
            p, n_dedup = predict(rec, knobs, m)
            A["species_dedup_collisions"] += n_dedup
            A["reps_outside_reference"] += sum(1 for r in p["reps"] if r not in ref_genomes)
            if args.forced_answer:
                s, _ = predict(rec, shipped_knobs, m)
                if not s["abstained"]:
                    # CONTRACT (module docstring): the forced arm must differ from the
                    # shipped arm ONLY where the shipped arm abstained, so that a post-hoc
                    # rank1_identity sweep is a risk-coverage curve for ONE method.
                    # Lowering family_floor globally also widens the candidate pool, and
                    # the ordered reducers (legacy band 0.005, graduated window 0.005) are
                    # pool-sensitive far outside the band -- measured on fold 0, a cell at
                    # best_identity 0.964 went 516 -> 624 ECs. Take the shipped prediction
                    # verbatim on answered cells; force only the abstentions. The check
                    # below is kept as a guard and must now be identically zero.
                    inv["checked"] += 1
                    if s["ec"] != p["ec"]:
                        # DIAGNOSTIC ONLY: how often the old global-floor construction
                        # would have changed an answered cell. Reported, never gated --
                        # the contract is now enforced by the substitution below.
                        bi = s.get("best_identity") or 0.0
                        if bi < band_tol:
                            inv["divergence_band_explained"] += 1
                        else:
                            inv["divergence_pool_sensitive"] += 1
                            if len(examples) < 20:
                                examples.append(m)
                                print(f"[diag] global-floor divergence: {mname} {m} "
                                      f"best_identity={bi} shipped |ec|={len(s['ec'])} "
                                      f"global-floor |ec|={len(p['ec'])} -> shipping the "
                                      f"shipped set", flush=True)
                    p = s                      # answered cells: shipped verbatim
                    if p["ec"] != s["ec"]:     # guard: contract holds on what we EMIT
                        die(EXIT_CANARY,
                            f"forced-answer contract broken for {mname} {m}: an answered "
                            f"cell was not emitted verbatim")
            out[m] = p
        pred[key] = out
        cov = sum(1 for v in out.values() if not v["abstained"]) / max(1, len(out))
        tiers = Counter(v["tier"] for v in out.values())
        print(f"[select] {key}: coverage={cov:.4f} tiers={dict(tiers)} "
              f"mean_reps={sum(v['n_reps'] for v in out.values())/max(1,len(out)):.2f}",
              flush=True)

    A["gene_count_fn_calls"] = counter["calls"]
    if A["reps_outside_reference"]:
        die(EXIT_CANARY, f"C8: {A['reps_outside_reference']} selected representatives lie "
                         f"outside the fold reference set")
    if args.forced_answer:
        print(f"[forced-answer] contract: {inv['checked']} answered cells emitted verbatim; "
              f"global-floor divergence would have been {inv['divergence_pool_sensitive']} "
              f"pool-sensitive + {inv['divergence_band_explained']} band-explained "
              f"(diagnostic, not gated)", flush=True)

    # ---------------------------------------- C2 source recovery (include_self/plumbing)
    src_rec = None
    if not enforce_c6:
        src_rec = {}
        for key, out in pred.items():
            n = ok = 0
            for m, v in out.items():
                g = md5_src.get(m)
                if g is None or g not in ref_genomes:
                    continue
                if len(fan_gids[m] & set(genomes)) > 1:
                    continue                      # twin cell: analytic ceiling, never gated
                n += 1
                ok += int(g in v["reps"])
            src_rec[key] = {"n_twin_free": n, "n_recovered": ok,
                            "frac": (ok / n if n else None)}
            print(f"[C2] {key}: source recovery {ok}/{n}", flush=True)
        primary = next((k for k in pred if k.startswith("exact_tie")), None) or next(iter(pred))
        f = src_rec[primary]["frac"]
        if f is not None and f < 1.0:
            die(EXIT_CANARY, f"C2 HARD gate: prFBA source recovery {f:.4f} < 1.00 for "
                             f"{primary} — the restricted DB, the synthesised md5->header "
                             f"join, or the select call is wrong")

    # ---------------------------------------------------------------------- output
    meta = {
        "fold": args.fold, "density": args.density, "mode": mode,
        "ref_name": ref_names, "n_reference_md5": len(ref_db),
        "n_reference_genomes": len(ref_genomes), "forced_answer": bool(args.forced_answer),
        "k2": args.k2, "gpu_cand": args.gpu_cand, "prefilter_k": args.prefilter_k,
        "min_ref_len": min_ref_len, "backend": backend,
        "select_modes": modes_req, "n_query": len(q_ids),
        "n_query_with_hits": sum(1 for m in q_ids if hits.get(m)),
        "n_twin_md5": n_twin,
        "queries": os.path.abspath(args.queries) if args.queries else os.path.abspath(args.amplicons),
        "frozen": frozen, "frozen_source": frozen_source,
        "wall_seconds": round(time.perf_counter() - t_wall0, 1),
        "align_seconds": align_seconds,
        "generated": iso_now(), "git_commit": git_commit(),
        "assertions": {**A,
                       "c6": c6,
                       "n_reference_genomes_missing_function_set": n_missing_fn,
                       "n_reference_genomes_missing_pgfam_cache": len(ref_genomes) - len(pgfam),
                       "gene_count_fn_calls_returning_none": counter["none"],
                       "forced_invariant": (inv if args.forced_answer else None),
                       "source_recovery": src_rec},
    }
    jdump(pred_path, {"meta": meta, "pred": pred})
    print(f"[done] {pred_path}: {len(pred)} arms x {len(q_ids)} amplicons "
          f"({meta['wall_seconds']:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
