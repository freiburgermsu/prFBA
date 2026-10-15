#!/usr/bin/env python
"""
plan_folds.py — panel QC, the pre-pass all-vs-all alignments, the similarity graph and
the GroupKFold fold plan for the SHARED-REFERENCE k-fold benchmark (SPEC.md §3, §4).

WHY THIS COMPONENT EXISTS
-------------------------
The predecessor head-to-head (picrust2_comparison/README.md) is confounded two ways: the
truth was BV-BRC's annotation, which is also the universe prFBA's predictions are drawn
from (scoring the SAME predictions against each system's own annotation flips the ranking:
prFBA +0.050 F1 under BV-BRC truth, -0.204 under PICRUSt2's), and the two reference sets
differ in size by 4x while only prFBA's is pruned. The shared-reference design removes both
by giving BOTH methods the SAME reference genome set and ONE annotation as both their gene
content and the answer key. That design is only honest if the held-out organisms are
genuinely held out — and "held out" is exactly what this script defines and then PROVES.

A naive leave-one-genome-out split is not a holdout: a genus-mate, a byte-identical amplicon
twin, or a 16S neighbour above prFBA's own species threshold is, for prediction purposes,
the organism itself. So the unit of holdout here is not a genome but a CONNECTED COMPONENT
of a similarity graph whose edges are the four ways one panel genome can stand in for
another (SPEC §3.2b):

  E1  same genus (or family/order, per --rank; the rank actually blocked)
  E2  a byte-identical amplicon md5 in any region
  E3  full-16S vs full-16S identity >= tau_full16s
  E4  an amplicon of u aligns to v's full 16S at identity >= tau_amplicon

E1/E2 are free (lineage strings, the amplicon fanout). E3/E4 require alignment, which is
why the pre-pass lives here: the graph cannot be built without it, and the fold plan cannot
be built without the graph. Three scripts would hand two intermediates to no other consumer,
hence one component with three phases (--phase panel|prepass|plan).

Per-amplicon identity is the criterion, not full-16S identity (SPEC D1): 27.4% of genome
pairs sharing a byte-identical amplicon fall BELOW 0.987 full-16S identity, so a full-16S
rule alone leaks exact twins. Group CV, not phylum stratification (SPEC D2): 133 of 206
phyla hold <10 genomes and 36 are singletons, so stratification is arithmetically
impossible; phylum/domain balance survives only as a bin-packing tiebreak.

ONE DEFINITION OF IDENTITY (SPEC §3.2). `identity` everywhere in this benchmark is the
field the prFBA aligner emits (edlib_biopython_hits._rescore -> n_matches / aligned_len,
local Smith-Waterman, BLASTN-like MATCH=2 MISMATCH=-3 GAP_OPEN=-5 GAP_EXTEND=-2), qualified
by a coverage floor aligned_len/query_len >= 0.90 (prFBA's shipped cov_min; for full-16S vs
full-16S the denominator is the shorter sequence). A hit below the floor is not an edge and
not a neighbour at any tau. That is precisely the quantity select_references._tier reads, so
the exclusion boundary is expressed in prFBA's native metric. glocal/global/MSA-derived
identities are not used anywhere.

WHAT THIS SCRIPT GUARANTEES
---------------------------
 1. Eligibility is applied ONCE and identically to a genome in its role as a reference and
    in its role as a test organism, so neither method ever sees a genome the other does not
    (SPEC §3.1). Dropped genomes leave the panel, both references, both methods' inputs and
    the truth; they are never scored.
 2. The >=tau neighbour sets are NOT truncated: the pre-pass asserts that the k2-th hit of
    every query falls below tau, and doubles --k2 (up to 3 times) if it does not. Without
    that assertion the graph silently loses edges exactly where it matters most.
 3. Every tip label carries --label-prefix (default "g_"). 3.61% of BV-BRC genome ids end in
    a trailing zero, and castor_nsti.R's read.table typing would delete those reference
    genomes from the NSTI known-set (SPEC C3). The prefix makes the label non-numeric.
 4. panel_md5_id.json is SYNTHESISED, not subset from 16S_md5_ID.json (which stores one
    header per md5 for the whole BV-BRC DB and may name a non-panel genome), and the
    parse_header round trip is asserted for every entry at write time.
 5. Every leakage_preflight counter is 0 for the rank actually blocked, or the script exits
    2. A non-zero counter means the graph is wrong and BOTH methods were handed a free
    near-exact neighbour, which is the one failure that would invalidate the whole benchmark.

Resumable: each phase skips work whose output already exists unless --force. The two
expensive loads (the 1.7 GB genome_function_sets.json and the 519 MB 16S_md5_seq.json) are
distilled to panel-sized caches under --outdir on first run.

Interpreter: /home/freiburger/Documents/py_venv/bin/python (the ONLY interpreter).
PICRUSt2 binaries (hmmalign) are reached by subprocess through
`bash -lc 'source ~/opt/picrust2/activate.sh && ...'` — never by importing picrust2 here.

    PY=~/Documents/py_venv/bin/python
    SR=/home/freiburger/Documents/prFBA/region_validation/shared_reference
    $PY $SR/scripts/plan_folds.py --phase panel  --limit-groups 40 --limit-genomes 300 --outdir $SR/data_smoke
    $PY $SR/scripts/plan_folds.py --phase prepass --outdir $SR/data_smoke --k2 200 --gpu-cand 300
    $PY $SR/scripts/plan_folds.py --phase plan    --outdir $SR/data_smoke --k 3

Exit codes: 0 ok | 2 design failure (leakage counter non-zero, or the largest component
exceeds 1.5*n_eligible/k) | 3 missing prerequisite artefact | 4 pre-pass truncation could
not be cleared by doubling --k2.
"""
from __future__ import annotations

import argparse
import collections
import gzip
import hashlib
import json
import math
import os
import random
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone

# --------------------------------------------------------------------------- paths
REPO = "/home/freiburger/Documents/prFBA"
SR_ROOT = os.path.join(REPO, "region_validation", "shared_reference")
DEFAULT_OUTDIR = os.path.join(SR_ROOT, "data")
RV_DATA = os.path.join(REPO, "region_validation", "data")
P2C_DATA = os.path.join(REPO, "region_validation", "picrust2_comparison", "data")

DEF_SELECTION = os.path.join(RV_DATA, "selection_10k.json")
DEF_FUNCTIONS = os.path.join(P2C_DATA, "genome_function_sets.json")
DEF_MD5_SEQ = "/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json"
DEF_AMPLICONS = os.path.join(RV_DATA, "combined_amplicons.unique.fasta")
DEF_FANOUT = os.path.join(RV_DATA, "amplicon_md5_fanout.json")
DEF_TRUTH = os.path.join(RV_DATA, "benchmark_truth.json")

P2_DEFAULTS = os.path.expanduser(
    "~/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files")
HMM = {"bac": os.path.join(P2_DEFAULTS, "bacteria", "bac_ref", "bac_ref.hmm"),
       "arc": os.path.join(P2_DEFAULTS, "archaea", "arc_ref", "arc_ref.hmm")}
PICRUST2_ACTIVATE = os.path.expanduser("~/opt/picrust2/activate.sh")
ALIGN_HITS = os.path.join(REPO, "align_hits.py")
PY = "/home/freiburger/Documents/py_venv/bin/python"

MIDAS_EMPTY_HEADER = "seq,Kingdom,Phylum,Class,Order,Family,Genus,Species,rel_ab"

# --------------------------------------------------- frozen configuration (SPEC §1)
# Recorded verbatim into every output JSON as meta.frozen. Nothing in this block may be
# changed after the first fold is scored.
FROZEN = {
    "tau_amplicon": 0.987,
    "tau_full16s": 0.987,
    "k_folds": 10,
    "seed": 1729,
    "prfba": {"k2": 200, "gpu_cand": 400, "prefilter_k": 1000, "min_ref_len": 200,
              "force_backend": "cuda",
              "select_modes": {
                  "legacy": {"select_all": False},
                  "exact_tie": {},
                  "graduated": {"select_all_window": 0.005,
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

# prFBA's shipped cov_min — the coverage floor of SPEC §3.2.
COV_MIN = 0.90

# Template is a pure function of the TRUE lineage (SPEC §9.2 / D10), computed here so every
# downstream component reads one value. Never auto_select_template: that would make the
# template depend on PREDICTED content.
GRAMPOS_PHYLA = {
    "bacillota", "firmicutes", "actinomycetota", "actinobacteria", "chloroflexota",
    "deinococcota", "thermotogota", "mycoplasmatota", "tenericutes",
}


def template_for(domain: str, phylum: str) -> str:
    """{Archaea -> archaea; GramPos phyla -> grampos; else gramneg} (SPEC §9.2)."""
    if (domain or "").strip().lower().startswith("arch"):
        return "archaea"
    p = (phylum or "").strip().lower()
    for stem in GRAMPOS_PHYLA:
        if p == stem or p.startswith(stem):
            return "grampos"
    return "gramneg"


def domain_key(domain: str) -> str:
    return "arc" if (domain or "").strip().lower().startswith("arch") else "bac"


# --------------------------------------------------------------------------- utilities
def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return "unknown"


def log(msg: str) -> None:
    print(f"[plan_folds] {msg}", flush=True)


def die(msg: str, code: int) -> None:
    print(f"[plan_folds][abort] {msg}", file=sys.stderr, flush=True)
    sys.exit(code)


def write_json(path: str, obj, indent=None) -> None:
    """Atomic write (tmp + os.replace) so a killed run never leaves a half-written artefact."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(os.path.abspath(path)), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(obj, fh, indent=indent, sort_keys=False)
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.unlink(tmp)


def load_json(path: str):
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        return json.load(fh)


def need(path: str, what: str) -> str:
    if not os.path.exists(path):
        die(f"missing {what}: {path} (run the earlier --phase first)", 3)
    return path


def fresh(path: str, force: bool) -> bool:
    """True when `path` already exists and may be reused (existence keyed; use fresh_fp for
    anything whose CONTENT depends on which genomes are in the panel)."""
    ok = (not force) and os.path.exists(path) and os.path.getsize(path) > 0
    if ok:
        log(f"reuse {os.path.basename(path)} (exists; --force to rebuild)")
    return ok


def panel_fingerprint(panel: dict) -> str:
    """sha256 over the eligible md5 set + label prefix — the identity of a panel.

    Existence-keyed resumability is a trap for this component: --outdir is shared between
    runs (and between components), so a smaller earlier run leaves behind a .sto, a prepass
    hits file and a neighbours table that a LARGER later run would happily reuse, silently
    dropping every genome the old artefact never saw. Measured in development: a 44-genome
    run's artefacts made a 300-genome panel drop ~260 genomes as 'no_alignment' while both
    phases reported success. Every panel-dependent artefact therefore carries this
    fingerprint and is rebuilt when it does not match.
    """
    md5s = sorted(r["md5"] for r in panel["genomes"].values())
    prefix = (panel.get("meta") or {}).get("label_prefix", "")
    return hashlib.sha256(("\n".join(md5s) + f"|{prefix}").encode()).hexdigest()[:16]


def fresh_fp(path: str, fp: str, force: bool) -> bool:
    """fresh(), but the artefact must also have been built for THIS panel fingerprint."""
    if force or not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    side = path + ".fp"
    got = None
    if os.path.exists(side):
        with open(side) as fh:
            got = fh.read().strip()
    if got != fp:
        log(f"rebuild {os.path.basename(path)}: panel fingerprint {got} != {fp}")
        return False
    log(f"reuse {os.path.basename(path)} (fingerprint {fp}; --force to rebuild)")
    return True


def stamp_fp(path: str, fp: str) -> None:
    with open(path + ".fp", "w") as fh:
        fh.write(fp + "\n")


def link_or_copy(src: str, dst: str) -> None:
    """Hardlink the aligner's hits file into place; copy only if the link fails. At full
    scale that file is GBs and a copy would double the disk for no reason."""
    if os.path.exists(dst):
        os.unlink(dst)
    try:
        os.link(src, dst)
    except OSError:
        shutil.copyfile(src, dst)


def iter_fasta(path: str):
    """(id, uppercase seq) for each record; id is the first whitespace-delimited field."""
    cur, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if cur is not None:
                    yield cur, "".join(buf).upper()
                cur, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if cur is not None:
        yield cur, "".join(buf).upper()


def write_fasta(path: str, items, width: int = 80) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in items:
            fh.write(f">{name}\n")
            if width:
                for i in range(0, len(seq), width):
                    fh.write(seq[i:i + width] + "\n")
            else:
                fh.write(seq + "\n")


def norm_ec(raw: str):
    """'EC:1.1.1.1' -> '1.1.1.1'; partial ECs containing '-' are DROPPED (SPEC §4.2).

    Mirrors compare.py's `c.replace("EC:", "")` so the vocabulary written here is the one
    the scorer and the PICRUSt2 trait table use.
    """
    e = (raw or "").replace("EC:", "").strip()
    if not e or "-" in e:
        return None
    return e


class UnionFind:
    def __init__(self, items):
        self.p = {x: x for x in items}
        self.r = {x: 0 for x in items}

    def find(self, x):
        p = self.p
        while p[x] != x:
            p[x] = p[p[x]]
            x = p[x]
        return x

    def union(self, a, b) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.r[ra] < self.r[rb]:
            ra, rb = rb, ra
        self.p[rb] = ra
        if self.r[ra] == self.r[rb]:
            self.r[ra] += 1
        return True

    def components(self):
        out = collections.defaultdict(list)
        for x in self.p:
            out[self.find(x)].append(x)
        return [sorted(v) for v in out.values()]


def bash_picrust2(cmd: str, log_path: str | None = None) -> None:
    """Run `cmd` inside the PICRUSt2 environment. Its python is NOT py_venv — nothing from
    this process is importable there, which is exactly why this is a subprocess."""
    full = f"source {shlex.quote(PICRUST2_ACTIVATE)} && {cmd}"
    log(f"picrust2$ {cmd}")
    proc = subprocess.run(["bash", "-lc", full], capture_output=True, text=True)
    if log_path:
        with open(log_path, "w") as fh:
            fh.write(f"$ {full}\n\n--- stdout ---\n{proc.stdout}\n--- stderr ---\n{proc.stderr}\n")
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "")[-2000:]
        die(f"picrust2 command failed (rc={proc.returncode}): {cmd}\n{tail}", 3)


# ======================================================================================
# PHASE panel — eligibility QC + the panel artefacts every later component reads
# ======================================================================================
def _cache_panel_ec(functions: str, panel_ids: set, cache: str, force: bool) -> dict:
    """{genome_id: [complete EC, ...]} for panel genomes only.

    genome_function_sets.json is 1.7 GB for 49,110 genomes; we need 10,000 of them. Streamed
    with ijson (C backend, ~12 s) and distilled to a panel-sized cache so every later run of
    this phase is instant. Genomes absent from the file get no entry (drop reason
    'no_function_set').
    """
    if fresh(cache, force):
        return load_json(cache)
    log(f"streaming EC sets for {len(panel_ids):,} panel genomes from {os.path.basename(functions)}")
    out = {}
    t0 = time.perf_counter()
    try:
        import ijson
        with open(functions, "rb") as fh:
            for gid, rec in ijson.kvitems(fh, "", use_float=True):
                if gid in panel_ids:
                    ecs = sorted({e for e in (norm_ec(x) for x in (rec.get("ec") or [])) if e})
                    out[gid] = {"ec": ecs, "n_ec_raw": len(rec.get("ec") or [])}
    except ImportError:                      # fallback: 1.7 GB in one gulp (works, slower)
        blob = load_json(functions)
        for gid in panel_ids:
            if gid in blob:
                rec = blob[gid]
                ecs = sorted({e for e in (norm_ec(x) for x in (rec.get("ec") or [])) if e})
                out[gid] = {"ec": ecs, "n_ec_raw": len(rec.get("ec") or [])}
    log(f"  {len(out):,}/{len(panel_ids):,} panel genomes have a function set "
        f"({time.perf_counter() - t0:.0f}s)")
    write_json(cache, out)
    return out


def _cache_panel_seq(md5_seq: str, md5s: set, cache: str, force: bool) -> dict:
    """{md5: UPPERCASE sequence} for the panel's representative 16S md5s (519 MB -> ~15 MB)."""
    if fresh(cache, force):
        return load_json(cache)
    log(f"streaming {len(md5s):,} panel 16S sequences from {os.path.basename(md5_seq)}")
    out = {}
    t0 = time.perf_counter()
    try:
        import ijson
        with open(md5_seq, "rb") as fh:
            for md5, seq in ijson.kvitems(fh, ""):
                if md5 in md5s:
                    out[md5] = seq.upper()
    except ImportError:
        blob = load_json(md5_seq)
        for md5 in md5s:
            if md5 in blob:
                out[md5] = blob[md5].upper()
    log(f"  {len(out):,}/{len(md5s):,} md5s found ({time.perf_counter() - t0:.0f}s)")
    write_json(cache, out)
    return out


def _free_edge_components(records, fanout_cells_by_genome):
    """Connected components of the FREE edge graph E1(genus) UNION E2(shared amplicon md5).

    Both are computable with no alignment, which is what makes the smoke panel sampling of
    SPEC §4.1 possible before the pre-pass has ever run.
    """
    ids = [r["genome_id"] for r in records]
    uf = UnionFind(ids)
    by_genus = collections.defaultdict(list)
    for r in records:
        g = (r.get("genus") or "").strip()
        if g:
            by_genus[g].append(r["genome_id"])
    for members in by_genus.values():
        for other in members[1:]:
            uf.union(members[0], other)
    by_md5 = collections.defaultdict(set)
    for gid, cells in fanout_cells_by_genome.items():
        if gid not in uf.p:
            continue
        for md5, _region in cells:
            by_md5[md5].add(gid)
    for owners in by_md5.values():
        owners = sorted(owners)
        for other in owners[1:]:
            uf.union(owners[0], other)
    return uf.components()


def _smoke_subsample(records, fanout_cells_by_genome, limit_groups, limit_genomes, seed):
    """SPEC §4.1 smoke panel: whole connected components of the free graph, not a flat
    subsample. A 'first N eligible' slice would contain almost no genus-mates and no md5
    twins, which would make the leakage assertions vacuous — they would pass because there
    is nothing to detect. Here the smoke panel CONTAINS leak structure by construction."""
    comps = _free_edge_components(records, fanout_cells_by_genome)
    by_id = {r["genome_id"]: r for r in records}
    multi = sorted([c for c in comps if len(c) > 1], key=lambda c: (-len(c), c[0]))
    single = sorted([c[0] for c in comps if len(c) == 1])
    rng = random.Random(seed)
    rng.shuffle(multi)
    rng.shuffle(single)
    budget = limit_genomes if limit_genomes > 0 else math.inf
    want_groups = limit_groups if limit_groups > 0 else math.inf
    chosen, n_groups = [], 0
    for comp in multi:
        if n_groups >= want_groups:
            break
        if len(chosen) + len(comp) > budget:
            continue                        # never truncate a component: that breaks the leak
        chosen.extend(comp)
        n_groups += 1
    for gid in single:                      # top up to --limit-genomes with singletons
        if len(chosen) >= budget:
            break
        chosen.append(gid)
    chosen = sorted(set(chosen))
    log(f"smoke panel: {n_groups} multi-genome free-graph components + singletons "
        f"-> {len(chosen)} genomes (of {len(records)} filter survivors)")
    return [by_id[g] for g in chosen]


def _parse_stockholm(path: str):
    """(name -> aligned seq, RF string) from an interleaved Stockholm alignment.

    Hand-rolled rather than via Bio.AlignIO because the ONLY thing needed is the `#=GC RF`
    match-column mask, and the HMM/MSA geometry assertion downstream depends on reading it
    exactly as hmmalign wrote it.
    """
    seqs = collections.defaultdict(list)
    rf = []
    order = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("# STOCKHOLM") or line == "//":
                continue
            if line.startswith("#=GC RF"):
                rf.append(line.split(None, 2)[2])
                continue
            if line.startswith("#"):
                continue                     # #=GS / #=GR (posterior probabilities) etc.
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            name, chunk = parts
            if name not in seqs:
                order.append(name)
            seqs[name].append(chunk.strip())
    aln = {n: "".join(seqs[n]) for n in order}
    return aln, "".join(rf)


def _hmmalign_occupancy(records, seq_by_md5, dom, outdir, force):
    """hmmalign the domain's panel 16S against the stock domain HMM, keep MATCH columns only
    (SPEC §5.1 step 3), and return {genome_id: (occupancy, match-column row)}.

    Insert columns are discarded: keeping them is the measured 69%-gap failure mode. Match
    columns are fixed by the HMM, so dropping a row after the fact (the occupancy filter)
    cannot change any surviving row — the alignment is never recomputed for that reason.
    """
    sub = [r for r in records if domain_key(r["domain"]) == dom]
    if not sub:
        return {}, {"n": 0, "width": 0, "mean_gap_frac": 0.0}
    work = os.path.join(outdir, "_hmmalign")
    os.makedirs(work, exist_ok=True)
    fna = os.path.join(work, f"panel_{dom}.prefilter.fna")
    write_fasta(fna, [(r["label"], seq_by_md5[r["best_16s_md5"]]) for r in sub])
    with open(fna, "rb") as fh:              # content-addressed: the alignment of THESE
        digest = hashlib.sha256(fh.read()).hexdigest()[:12]   # sequences, never another set
    sto = os.path.join(work, f"panel_{dom}.{digest}.sto")
    if not fresh(sto, force):
        bash_picrust2(
            f"hmmalign --trim --dna --outformat Stockholm -o {shlex.quote(sto)} "
            f"{shlex.quote(HMM[dom])} {shlex.quote(fna)}",
            log_path=os.path.join(work, f"hmmalign_{dom}.log"))
    aln, rf = _parse_stockholm(sto)
    if not rf:
        die(f"hmmalign produced no '#=GC RF' line for {dom}: cannot identify match columns", 3)
    mask = [i for i, c in enumerate(rf) if c != "."]
    label2gid = {r["label"]: r["genome_id"] for r in sub}
    rows, gaps = {}, 0
    for label, row in aln.items():
        gid = label2gid.get(label)
        if gid is None:
            continue
        if len(row) != len(rf):
            die(f"hmmalign row {label} is {len(row)} cols but RF is {len(rf)}", 3)
        kept = "".join(row[i] for i in mask).upper().replace(".", "-")
        occ = sum(1 for c in kept if c != "-") / max(1, len(mask))
        gaps += len(kept) - sum(1 for c in kept if c != "-")
        rows[gid] = (occ, kept)
    geom = {"n": len(rows), "width": len(mask),
            "mean_gap_frac": round(gaps / max(1, len(rows) * len(mask)), 6)}
    log(f"  {dom}: {geom['n']} rows x {geom['width']} match columns, "
        f"mean gap {geom['mean_gap_frac']:.4f}")
    if geom["mean_gap_frac"] > 0.10 or geom["width"] > 2000:
        log(f"  [warn] {dom} MSA geometry outside the SPEC §5.1 expectation "
            f"(gap<=0.10, width<=2000) — build_reference.py will abort on this")
    return rows, geom


def phase_panel(args) -> None:
    out = args.outdir
    os.makedirs(out, exist_ok=True)
    panel_path = os.path.join(out, "panel_eligible.json")

    sel = load_json(need(args.selection, "selection_10k.json"))
    panel = sel["genomes"]
    n_panel = len(panel)
    log(f"panel: {n_panel:,} genomes from {os.path.basename(args.selection)}")

    ec = _cache_panel_ec(args.functions, {g["genome_id"] for g in panel},
                         os.path.join(out, "_cache_panel_ec.json"), args.force)
    seq_by_md5 = _cache_panel_seq(args.md5_seq, {g["best_16s_md5"] for g in panel},
                                  os.path.join(out, "_cache_panel_md5_seq.json"), args.force)

    # ---- filters 1,2,3,5 (free); filter 4 (occupancy) needs the alignment, below -------
    drop_counts = collections.Counter()
    dropped: dict[str, str] = {}
    cand = []
    qc = FROZEN["reference_qc"]
    for g in panel:
        gid = g["genome_id"]
        rec = ec.get(gid)
        if rec is None:
            dropped[gid] = "no_function_set"; drop_counts["no_function_set"] += 1; continue
        if len(rec["ec"]) < qc["ec_min"]:
            dropped[gid] = "empty_ec"; drop_counts["empty_ec"] += 1; continue
        if not (qc["len_min"] <= g["best_16s_len"] <= qc["len_max"]):
            dropped[gid] = "length"; drop_counts["length"] += 1; continue
        seq = seq_by_md5.get(g["best_16s_md5"])
        if seq is None:
            dropped[gid] = "no_sequence"; drop_counts["no_sequence"] += 1; continue
        if set(seq) - set("ACGTN"):
            dropped[gid] = "bad_chars"; drop_counts["bad_chars"] += 1; continue
        r = dict(g)
        r["label"] = f"{args.label_prefix}{gid}"
        r["n_ec"] = len(rec["ec"])
        cand.append(r)
    log(f"free filters: {len(cand):,} survive | drops {dict(drop_counts)}")

    # ---- smoke panel (whole free-graph components; see _smoke_subsample) ---------------
    fanout = load_json(need(args.fanout, "amplicon_md5_fanout.json"))
    cells_by_genome = collections.defaultdict(list)
    for md5, cells in fanout.items():
        for cell in cells:
            gid, _, region = cell.partition("__")
            cells_by_genome[gid].append((md5, region))
    if args.limit_groups or args.limit_genomes:
        keep = _smoke_subsample(cand, cells_by_genome, args.limit_groups,
                                args.limit_genomes, FROZEN["seed"])
        kept_ids = {r["genome_id"] for r in keep}
        for r in cand:
            if r["genome_id"] not in kept_ids:
                dropped[r["genome_id"]] = "smoke_subsample"
                drop_counts["smoke_subsample"] += 1
        cand = keep

    # ---- filter 4: HMM match-column occupancy -----------------------------------------
    occ_rows, geom = {}, {}
    for dom in ("bac", "arc"):
        rows, g = _hmmalign_occupancy(cand, seq_by_md5, dom, out, args.force)
        occ_rows.update(rows)
        geom[dom] = g
    eligible, afa = [], {"bac": [], "arc": []}
    for r in cand:
        gid = r["genome_id"]
        if gid not in occ_rows:
            dropped[gid] = "no_alignment"; drop_counts["no_alignment"] += 1; continue
        occ, row = occ_rows[gid]
        if occ < qc["occupancy_min"]:
            dropped[gid] = "occupancy"; drop_counts["occupancy"] += 1; continue
        r["occupancy"] = round(occ, 4)
        eligible.append(r)
        afa[domain_key(r["domain"])].append((r["label"], row))
    n_elig = len(eligible)
    log(f"eligible: {n_elig:,} genomes (dropped {len(dropped):,})")
    if n_elig == 0:
        die("no eligible genomes survived QC", 3)

    # ---- artefacts ---------------------------------------------------------------------
    for dom in ("bac", "arc"):
        sub = [r for r in eligible if domain_key(r["domain"]) == dom]
        write_fasta(os.path.join(out, f"panel_{dom}.fna"),
                    [(r["label"], seq_by_md5[r["best_16s_md5"]]) for r in sub])
        write_fasta(os.path.join(out, f"panel_{dom}.afa"), afa[dom])

    panel_md5_seq = {r["best_16s_md5"]: seq_by_md5[r["best_16s_md5"]] for r in eligible}
    if len(panel_md5_seq) != n_elig:
        die(f"md5 collision: {len(panel_md5_seq)} distinct md5s for {n_elig} eligible genomes "
            f"— one md5 per eligible genome is required (SPEC §4.2)", 3)
    write_json(os.path.join(out, "panel_md5_seq.json"), panel_md5_seq)

    # SYNTHESISED md5 -> BV-BRC-style header, asserted through the real parser (SPEC §4.2).
    sys.path.insert(0, REPO)
    import edlib_biopython_hits as EBH        # noqa: E402  (reuse the real parse_header)
    panel_md5_id = {}
    for r in eligible:
        gid, md5 = r["genome_id"], r["best_16s_md5"]
        hdr = (f"fig|{gid}.rna.1|   16S ribosomal RNA   "
               f"[{r.get('species') or ''} | {gid}]")
        got = EBH.parse_header(hdr)
        # SPEC §4.2 states the 4th element is 'fig|<gid>.rna.1|'; the real parser's
        # _HDR_FEAT = r"^fig\|([^|]+)\|" captures the INNER feature id, so it is
        # '<gid>.rna.1'. The header string written here is exactly the specified one — only
        # this assertion's expectation is corrected to what parse_header actually returns.
        want = (r.get("species") or "", gid, int(gid.split(".")[0]), f"{gid}.rna.1")
        if got != want:
            die(f"panel_md5_id header round trip failed for {gid}: {got} != {want}", 3)
        panel_md5_id[md5] = hdr
    write_json(os.path.join(out, "panel_md5_id.json"), panel_md5_id)

    # EC universe: the FIXED trait/score vocabulary, identical in every fold and density.
    prev = collections.Counter()
    for r in eligible:
        for e in ec[r["genome_id"]]["ec"]:
            prev[e] += 1
    ec_sorted = sorted(prev)                 # plain lexicographic: stable and reproducible
    write_json(os.path.join(out, "panel_ec_universe.json"), {
        "meta": {"n_genomes": n_elig, "n_ec": len(ec_sorted)},
        "ec": ec_sorted,
        "prevalence": {e: round(prev[e] / n_elig, 6) for e in ec_sorted},
        "ec_index": {e: i for i, e in enumerate(ec_sorted)},
    })

    # queries for --phase prepass (headers = md5, so the amplicon-vs-reference join is the
    # same md5 join everywhere).
    write_fasta(os.path.join(out, "panel_16s_queries.fna"),
                [(r["best_16s_md5"], panel_md5_seq[r["best_16s_md5"]]) for r in eligible])
    with open(os.path.join(out, "midas_empty.csv"), "w") as fh:
        fh.write(MIDAS_EMPTY_HEADER + "\n")

    genomes_out = {}
    for r in eligible:
        genomes_out[r["genome_id"]] = {
            "label": r["label"], "taxon_id": r["taxon_id"], "md5": r["best_16s_md5"],
            "len": len(seq_by_md5[r["best_16s_md5"]]), "domain": r["domain"],
            "phylum": r.get("phylum") or "", "class": r.get("class") or "",
            "order": r.get("order") or "", "family": r.get("family") or "",
            "genus": r.get("genus") or "", "species": r.get("species") or "",
            "n_full_length_16s": r.get("n_full_length_16s") or 1,
            "n_ec": r["n_ec"], "occupancy": r["occupancy"],
            "template": template_for(r["domain"], r.get("phylum")),
        }
    fp = hashlib.sha256(("\n".join(sorted(r["best_16s_md5"] for r in eligible))
                         + f"|{args.label_prefix}").encode()).hexdigest()[:16]
    write_json(panel_path, {
        "meta": {"frozen": FROZEN, "generated": now_iso(), "git_commit": git_commit(),
                 "n_panel": n_panel, "n_eligible": n_elig, "script": "plan_folds.py",
                 "phase": "panel", "label_prefix": args.label_prefix,
                 "panel_fingerprint": fp,
                 "msa_geometry": geom,
                 "ec_filter_note": ("empty_ec counts genomes with zero COMPLETE EC numbers "
                                    "after norm_ec (partial 'x.x.x.-' ECs dropped), so it is "
                                    ">= the 161 raw-empty figure in SPEC §3.1"),
                 "drop_counts": {k: int(v) for k, v in sorted(drop_counts.items())}},
        "genomes": genomes_out,
        "dropped": dict(sorted(dropped.items())),
    }, indent=1)
    log(f"wrote {panel_path} (n_eligible={n_elig}, fingerprint={fp})")


# ======================================================================================
# PHASE prepass — the two all-vs-all alignments that E3/E4 and the novelty covariate need
# ======================================================================================
def _kvitems(path: str):
    """Stream a top-level JSON object as (key, value). The aligner's hits file is
    n_query x k2 rich records (GBs at full scale), so it is never json.load-ed."""
    try:
        import ijson
        with open(path, "rb") as fh:
            for k, v in ijson.kvitems(fh, "", use_float=True):
                yield k, v
    except ImportError:
        for k, v in load_json(path).items():
            yield k, v


def _run_aligner(queries: str, hits_dir: str, *, ref_db: str, ref_id_db: str, midas: str,
                 k2: int, gpu_cand: int, prefilter_k: int, workers: int,
                 force_backend: str) -> str:
    """align_hits.py against a RESTRICTED reference DB (component 0's --ref-db).

    The restriction is the whole point: a post-hoc filter of full-DB top-k hits reproduces
    the k2-truncation artifact (94.6% of self-exclusion abstentions came from emptying an
    already-truncated record), so the reference set is restricted at ALIGNMENT time. The
    same --ref-db mechanism is what run_prfba_fold.py uses per fold, so the pre-pass and the
    fold runs measure identity with one code path.
    """
    os.makedirs(hits_dir, exist_ok=True)
    cmd = [PY, ALIGN_HITS, "--fasta", queries, "--outdir", hits_dir,
           "--ref-db", ref_db, "--ref-id-db", ref_id_db, "--taxonomy-csv", midas,
           "--k2", str(k2), "--gpu-cand", str(max(gpu_cand, k2)),
           "--prefilter-k", str(max(prefilter_k, k2)),
           "--min-ref-len", str(FROZEN["prfba"]["min_ref_len"]),
           "--workers", str(workers)]
    if force_backend != "auto":
        cmd += ["--force-backend", force_backend]
    log("aligner$ " + " ".join(shlex.quote(c) for c in cmd))
    proc = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True)
    with open(os.path.join(hits_dir, "align_hits.log"), "w") as fh:
        fh.write(proc.stdout + "\n--- stderr ---\n" + proc.stderr)
    if proc.returncode != 0:
        die(f"align_hits.py failed (rc={proc.returncode}); tail:\n{proc.stderr[-2000:]}", 3)
    hits = os.path.join(hits_dir, "asv_top20_alignment_hits.json")
    if not os.path.exists(hits):
        die(f"align_hits.py produced no {hits}", 3)
    return hits


def _truncation_violations(hits: str, k2: int, n_ref: int, tau: float):
    """How many queries have their >=tau neighbour set possibly CUT by k2.

    A violation is a query whose k2-th (last) emitted hit is still at identity >= tau while
    the reference set is deeper than k2: there may be further >=tau neighbours the graph
    never sees, i.e. missing E3/E4 edges exactly where leakage is worst. If k2 >= n_ref the
    whole reference was emitted and truncation is impossible.
    """
    if n_ref <= k2:
        return 0, []
    bad = []
    for qid, entry in _kvitems(hits):
        top = entry.get("top20") or []
        if len(top) >= k2 and (top[-1].get("identity") or 0.0) >= tau:
            bad.append(qid)
    return len(bad), bad[:20]


def _aligner_with_truncation_guard(queries, hits_dir, *, ref_db, ref_id_db, midas, k2,
                                   gpu_cand, prefilter_k, workers, force_backend, n_ref,
                                   tau, label):
    """Run the aligner, assert the >=tau neighbour set is not truncated, double --k2 up to
    3 times if it is, then abort (SPEC §4.1). Silently accepting truncation here would
    quietly delete edges from the similarity graph."""
    for attempt in range(4):
        hits = _run_aligner(queries, hits_dir, ref_db=ref_db, ref_id_db=ref_id_db,
                            midas=midas, k2=k2, gpu_cand=gpu_cand, prefilter_k=prefilter_k,
                            workers=workers, force_backend=force_backend)
        viol, examples = _truncation_violations(hits, k2, n_ref, tau)
        log(f"{label}: k2={k2} truncation violations={viol}")
        if viol == 0:
            return hits, k2, 0
        if attempt == 3:
            die(f"{label}: {viol} queries still truncated at k2={k2} (e.g. {examples}); "
                f"the >=tau neighbour set cannot be trusted", 4)
        k2 *= 2
        gpu_cand = max(gpu_cand * 2, k2)
        prefilter_k = max(prefilter_k, k2)
        log(f"{label}: retrying with k2={k2} (doubling, attempt {attempt + 2}/4)")
    return hits, k2, viol                     # unreachable


def _inject_errors(src_fasta: str, dst: str, rate: float, indel_frac: float, seed: int):
    """Substitutions + indels on TEST AMPLICONS ONLY, same md5 headers, reference untouched.

    In-silico PCR queries are exact substrings of their reference 16S; real reads are not.
    This arm is the degradation slope, and it is the identical pipeline on mutated queries —
    which is why the headers (and therefore every md5 join downstream) are preserved.
    """
    rng = random.Random(seed)
    alphabet = "ACGT"
    out = []
    for name, seq in iter_fasta(src_fasta):
        s = list(seq)
        n_sub = int(round(rate * len(s)))
        n_indel = int(round(rate * indel_frac * len(s)))
        for pos in rng.sample(range(len(s)), min(n_sub, len(s))):
            cur = s[pos]
            s[pos] = rng.choice([b for b in alphabet if b != cur])
        for _ in range(n_indel):
            pos = rng.randrange(len(s))
            if rng.random() < 0.5 and len(s) > 1:
                del s[pos]
            else:
                s.insert(pos, rng.choice(alphabet))
        out.append((name, "".join(s)))
    write_fasta(dst, out)
    log(f"wrote {os.path.basename(dst)}: {len(out)} mutated queries "
        f"(rate={rate}, indel_frac={indel_frac}, seed={seed})")


def phase_prepass(args) -> None:
    out = args.outdir
    panel = load_json(need(os.path.join(out, "panel_eligible.json"), "panel_eligible.json"))
    genomes = panel["genomes"]
    ref_db = need(os.path.join(out, "panel_md5_seq.json"), "panel_md5_seq.json")
    ref_id_db = need(os.path.join(out, "panel_md5_id.json"), "panel_md5_id.json")
    midas = os.path.join(out, "midas_empty.csv")
    if not os.path.exists(midas):
        with open(midas, "w") as fh:
            fh.write(MIDAS_EMPTY_HEADER + "\n")
    md5_to_gid = {rec["md5"]: gid for gid, rec in genomes.items()}
    n_ref = len(md5_to_gid)
    tau_f, tau_a = FROZEN["tau_full16s"], FROZEN["tau_amplicon"]
    fp = panel_fingerprint(panel)
    log(f"panel fingerprint {fp} ({n_ref} reference md5s)")

    # ---- amplicon queries: every amplicon md5 that owns a cell of an eligible genome -----
    fanout = load_json(need(args.fanout, "amplicon_md5_fanout.json"))
    amp_owners: dict[str, list] = collections.defaultdict(list)
    for md5, cells in fanout.items():
        for cell in cells:
            gid, _, region = cell.partition("__")
            if gid in genomes:
                amp_owners[md5].append((gid, region))
    amp_q = os.path.join(out, "panel_amplicon_queries.fna")
    if not fresh_fp(amp_q, fp, args.force):
        want = set(amp_owners)
        recs = [(m, s) for m, s in iter_fasta(need(args.amplicons, "unique amplicon fasta"))
                if m in want]
        if len(recs) != len(want):
            log(f"[warn] {len(want) - len(recs)} amplicon md5s of eligible cells are absent "
                f"from {os.path.basename(args.amplicons)}")
        write_fasta(amp_q, recs)
        stamp_fp(amp_q, fp)
        log(f"wrote {os.path.basename(amp_q)}: {len(recs)} amplicon queries")
    full_q = need(os.path.join(out, "panel_16s_queries.fna"), "panel_16s_queries.fna")

    # The error-injected arm is the IDENTICAL pipeline on mutated queries, so its query file
    # is produced here, from the same amplicon FASTA, independently of everything below
    # (both fold runners take --queries PATH and read it instead of slicing the clean FASTA).
    if args.inject_error > 0:
        tag = f"{int(round(args.inject_error * 1000)):03d}"
        dst = os.path.join(out, f"queries_err{tag}.fna")
        if not fresh_fp(dst, fp, args.force):   # fingerprint-gated: it mirrors amp_q
            _inject_errors(amp_q, dst, args.inject_error, args.inject_indel_frac,
                           args.inject_seed)
            stamp_fp(dst, fp)

    # ---- the two alignments ---------------------------------------------------------------
    prepass_full = os.path.join(out, "prepass_full16s.json")
    prepass_amp = os.path.join(out, "prepass_amplicon.json")
    k2_used = {"full16s": args.k2, "amplicon": args.k2}
    viol_total = 0
    for label, qfa, dest in (("full16s", full_q, prepass_full), ("amplicon", amp_q, prepass_amp)):
        hits_dir = os.path.join(out, f"hits_{label}")
        if fresh_fp(dest, fp, args.force):
            # The reused run may have doubled k2 past --k2; read what it ACTUALLY used so
            # neighbors.meta.k2 never misreports the depth the edges were drawn from.
            ms = os.path.join(hits_dir, "method_selection.json")
            if os.path.exists(ms):
                try:
                    k2_used[label] = int(load_json(ms)["params"]["k2"])
                except Exception:
                    pass
            continue
        hits, k2u, viol = _aligner_with_truncation_guard(
            qfa, hits_dir, ref_db=ref_db, ref_id_db=ref_id_db,
            midas=midas, k2=args.k2, gpu_cand=args.gpu_cand, prefilter_k=args.prefilter_k,
            workers=args.threads, force_backend=args.force_backend, n_ref=n_ref,
            tau=(tau_f if label == "full16s" else tau_a), label=label)
        link_or_copy(hits, dest)
        stamp_fp(dest, fp)
        k2_used[label] = k2u
        viol_total += viol
        log(f"wrote {os.path.basename(dest)} (k2={k2u})")

    # ---- distil neighbors.json ------------------------------------------------------------
    nb_path = os.path.join(out, "neighbors.json")
    if fresh_fp(nb_path, fp, args.force):
        return                               # every other prepass artefact is already written
    full16s = collections.defaultdict(list)
    below_floor = collections.Counter()
    for qid, entry in _kvitems(prepass_full):
        gid = md5_to_gid.get(qid)
        if gid is None:
            continue
        qlen = entry.get("asv_len") or 0
        for h in entry.get("top20") or []:
            hg = md5_to_gid.get(h.get("md5"))
            if hg is None or hg == gid:
                continue
            denom = max(1, min(qlen, h.get("ref_seq_len") or qlen))
            if (h.get("aligned_len") or 0) / denom < COV_MIN:
                below_floor["full16s"] += 1
                continue
            ident = float(h.get("identity") or 0.0)
            if ident >= tau_f:
                full16s[gid].append([hg, round(ident, 4)])

    amplicon: dict = collections.defaultdict(lambda: collections.defaultdict(dict))
    rank1: dict = {}
    for qid, entry in _kvitems(prepass_amp):
        qlen = entry.get("asv_len") or 0
        ranked = []
        for h in entry.get("top20") or []:
            hg = md5_to_gid.get(h.get("md5"))
            if hg is None:
                continue
            if (h.get("aligned_len") or 0) / max(1, qlen) < COV_MIN:
                below_floor["amplicon"] += 1
                continue
            ranked.append((hg, round(float(h.get("identity") or 0.0), 4),
                           float(h.get("align_score") or 0.0)))
        # best-first by the ONE identity definition (then SW score, then id) so a consumer
        # taking the first surviving entry gets the highest-identity surviving reference.
        ranked.sort(key=lambda t: (-t[1], -t[2], t[0]))
        rank1[qid] = [[g, i] for g, i, _s in ranked]
        for gid, region in amp_owners.get(qid, []):
            for hg, ident, _s in ranked:
                if hg == gid or ident < tau_a:
                    continue
                prev = amplicon[gid][region].get(hg)
                if prev is None or ident > prev:
                    amplicon[gid][region][hg] = ident

    md5_share = collections.defaultdict(list)
    for md5, owners in amp_owners.items():
        if len({g for g, _r in owners}) < 2:
            continue
        for gid, region in owners:
            for other, _r2 in owners:
                if other != gid:
                    md5_share[gid].append([other, region])

    write_json(nb_path, {
        "meta": {"tau_full16s": tau_f, "tau_amplicon": tau_a,
                 "k2": k2_used["amplicon"], "k2_full16s": k2_used["full16s"],
                 "truncation_violations": viol_total,
                 "coverage_floor": COV_MIN,
                 "n_hits_below_coverage_floor": dict(below_floor),
                 "n_reference_md5": n_ref, "panel_fingerprint": fp,
                 "generated": now_iso(),
                 "git_commit": git_commit(), "frozen": FROZEN,
                 "rank1_note": ("full coverage-qualified hit list per amplicon md5, sorted by "
                                "identity desc; the source genome is NOT excluded (fold "
                                "membership excludes it downstream)")},
        "full16s": {g: sorted(v, key=lambda t: (-t[1], t[0])) for g, v in sorted(full16s.items())},
        "amplicon": {g: {r: sorted(([h, i] for h, i in d.items()), key=lambda t: (-t[1], t[0]))
                         for r, d in sorted(regions.items())}
                     for g, regions in sorted(amplicon.items())},
        "md5_share": {g: sorted(v) for g, v in sorted(md5_share.items())},
        "rank1": rank1,
    })
    stamp_fp(nb_path, fp)
    log(f"wrote {nb_path}: full16s edges for {len(full16s)} genomes, "
        f"amplicon edges for {len(amplicon)} genomes, md5_share for {len(md5_share)}, "
        f"rank1 for {len(rank1)} md5s")


# ======================================================================================
# PHASE plan — the similarity graph, the GroupKFold plan, and the leakage preflight
# ======================================================================================
RANK_EDGE = {"genus": "genus", "family": "family", "order": "order", "md5only": None}


def _composition(ids, genomes, phyla):
    """Domain + top-5-phylum fractions of a genome set (the soft balance target, SPEC §3.3)."""
    n = max(1, len(ids))
    dom = collections.Counter(genomes[g]["domain"] for g in ids)
    ph = collections.Counter(genomes[g]["phylum"] for g in ids)
    v = {f"dom:{d}": dom.get(d, 0) / n for d in ("Bacteria", "Archaea")}
    v.update({f"ph:{p}": ph.get(p, 0) / n for p in phyla})
    return v


def _comp_distance(a, b):
    return sum(abs(a[k] - b.get(k, 0.0)) for k in a)


def build_graph(genomes, neighbors, cells_by_genome, rank, disable):
    """Union-find over the four SPEC §3.2b edge types; returns (components, edge census).

    `disable` exists ONLY to prove the leakage assertions fire (smoke gate item 1): running
    --phase plan with --disable-edges e4 must exit 2, which is the evidence that a passing
    run means something. It is never used in a real run.
    """
    ids = sorted(genomes)
    uf = UnionFind(ids)
    census = collections.Counter()

    rank_key = RANK_EDGE[rank]
    if rank_key and "e1" not in disable:
        by_rank = collections.defaultdict(list)
        for g in ids:
            v = (genomes[g].get(rank_key) or "").strip()
            if v:
                by_rank[v].append(g)
        for members in by_rank.values():
            for other in members[1:]:
                census["e1_pairs"] += 1
                uf.union(members[0], other)
        # taxon_id is folded into E1 for every taxonomic rank: C6c asserts no reference
        # shares a test genome's taxon_id, and a genome with an EMPTY rank string would
        # otherwise escape E1 and make that assertion unsatisfiable. (On this panel all
        # 10,000 taxon_ids are distinct, so this adds no merge — it is a guard, not a policy.)
        by_taxon = collections.defaultdict(list)
        for g in ids:
            by_taxon[genomes[g]["taxon_id"]].append(g)
        for members in by_taxon.values():
            for other in members[1:]:
                census["e1_taxon_pairs"] += 1
                uf.union(members[0], other)

    if "e2" not in disable:
        by_md5 = collections.defaultdict(set)
        for gid, cells in cells_by_genome.items():
            for md5, _region in cells:
                by_md5[md5].add(gid)
        for owners in by_md5.values():
            owners = sorted(owners)
            for other in owners[1:]:
                census["e2_pairs"] += 1
                uf.union(owners[0], other)

    if "e3" not in disable:
        for gid, hits in (neighbors.get("full16s") or {}).items():
            if gid not in uf.p:
                continue
            for other, _ident in hits:
                if other in uf.p:
                    census["e3_pairs"] += 1
                    uf.union(gid, other)

    if "e4" not in disable:
        for gid, regions in (neighbors.get("amplicon") or {}).items():
            if gid not in uf.p:
                continue
            for _region, hits in regions.items():
                for other, _ident in hits:
                    if other in uf.p:
                        census["e4_pairs"] += 1
                        uf.union(gid, other)

    comps = sorted(uf.components(), key=lambda c: (-len(c), c[0]))
    return comps, census


def pack_kfold(comps, genomes, k, phyla):
    """Greedy size-descending bin packing into k bins: each component joins the currently
    SMALLEST bin; ties go to the bin whose domain/top-5-phylum composition is furthest from
    the panel marginal (soft balance only — group integrity is never traded for balance)."""
    marginal = _composition(sorted(genomes), genomes, phyla)
    bins, sizes = [[] for _ in range(k)], [0] * k
    for comp in comps:
        m = min(sizes)
        tied = [i for i in range(k) if sizes[i] == m]
        if len(tied) == 1:
            pick = tied[0]
        else:
            best, bestd = tied[0], -1.0
            for i in tied:
                members = [g for c in bins[i] for g in c]
                d = math.inf if not members else _comp_distance(
                    marginal, _composition(members, genomes, phyla))
                if d > bestd:
                    best, bestd = i, d
            pick = best
        bins[pick].append(comp)
        sizes[pick] += len(comp)
    return bins


def pack_blocks(comps, genomes, blocks, target_total, phyla, seed):
    """Secondary ladder (SPEC §3.4): the family graph has a 2,879-genome giant component
    (28.8% of the panel) and CANNOT be 10-folded, so leave-family-out runs a few HOLDOUT
    BLOCKS drawn from a phylum-proportional subsample instead of a partition. Genomes in no
    block are references in every block."""
    rng = random.Random(seed)
    panel_ph = collections.Counter(genomes[g]["phylum"] for g in genomes)
    n_panel = max(1, len(genomes))
    remaining = sorted(comps, key=lambda c: (-len(c), c[0]))
    rng.shuffle(remaining)
    per = max(1, int(math.ceil(target_total / max(1, blocks))))
    out = []
    chosen_ph = collections.Counter()
    for _b in range(blocks):
        bin_, size = [], 0
        while size < per and remaining:
            best_i, best_score = None, None
            for i, comp in enumerate(remaining):
                if size + len(comp) > per * 1.25:
                    continue
                ph = collections.Counter(genomes[g]["phylum"] for g in comp).most_common(1)[0][0]
                n_sel = max(1, sum(chosen_ph.values()))
                deficit = panel_ph[ph] / n_panel - chosen_ph[ph] / n_sel
                score = (deficit, len(comp), -i)
                if best_score is None or score > best_score:
                    best_i, best_score = i, score
            if best_i is None:
                break
            comp = remaining.pop(best_i)
            bin_.append(comp)
            size += len(comp)
            for g in comp:
                chosen_ph[genomes[g]["phylum"]] += 1
        out.append(bin_)
    return out


def leakage_preflight(test_ids, ref_ids, genomes, neighbors, cells_by_genome):
    """The five counters of SPEC §4.2 / control C6. Every one MUST be 0 for the rank blocked:
    a non-zero value means the graph is wrong and BOTH methods were handed a free near-exact
    neighbour of the held-out organism, which would invalidate the comparison outright."""
    ref = set(ref_ids)
    test = set(test_ids)
    ref_md5 = {md5 for g in ref for md5, _r in cells_by_genome.get(g, [])}
    test_md5 = {md5 for g in test for md5, _r in cells_by_genome.get(g, [])}
    amp = neighbors.get("amplicon") or {}
    full = neighbors.get("full16s") or {}
    n_amp = sum(1 for g in test for _r, hits in (amp.get(g) or {}).items()
                for h, _i in hits if h in ref)
    n_full = sum(1 for g in test for h, _i in (full.get(g) or []) if h in ref)
    ref_taxa = {genomes[g]["taxon_id"] for g in ref}
    ref_genera = {(genomes[g].get("genus") or "").strip() for g in ref} - {""}
    return {
        "md5_shared_with_reference": len(test_md5 & ref_md5),
        "amplicon_identity_ge_tau": n_amp,
        "full16s_identity_ge_tau": n_full,
        "taxon_id_shared": sum(1 for g in test if genomes[g]["taxon_id"] in ref_taxa),
        "genus_shared": sum(1 for g in test
                            if (genomes[g].get("genus") or "").strip() in ref_genera),
    }


def phase_plan(args) -> None:
    out = args.outdir
    panel = load_json(need(os.path.join(out, "panel_eligible.json"), "panel_eligible.json"))
    genomes = panel["genomes"]
    neighbors = load_json(need(os.path.join(out, "neighbors.json"),
                               "neighbors.json (run --phase prepass)"))
    # A neighbors.json built for a DIFFERENT panel is the one stale-artefact failure that
    # would pass silently: the E3/E4 lookups simply find nothing, the graph comes out all
    # singletons, and the leakage preflight reports 0 because it is consulting the same
    # empty tables. Measured in development when a parallel component overwrote
    # panel_eligible.json in a shared --outdir. Abort instead.
    nb_meta = neighbors.get("meta") or {}
    fp_panel, fp_nb = panel_fingerprint(panel), nb_meta.get("panel_fingerprint")
    if fp_nb is not None and fp_nb != fp_panel:
        die(f"neighbors.json was built for panel fingerprint {fp_nb} but panel_eligible.json "
            f"is {fp_panel}: the pre-pass is stale. Re-run --phase prepass "
            f"(fingerprint-gated artefacts rebuild themselves; nothing needs deleting).", 3)
    nb_n = nb_meta.get("n_reference_md5")
    if nb_n is not None and nb_n != len(genomes):
        die(f"neighbors.json was built for {nb_n} reference md5s but panel_eligible.json has "
            f"{len(genomes)} eligible genomes: the pre-pass is stale. Re-run "
            f"--phase prepass --force (a stale neighbours table yields an all-singleton graph "
            f"whose leakage preflight passes for the wrong reason).", 3)
    unknown = [g for g in (neighbors.get("full16s") or {}) if g not in genomes]
    unknown += [g for g in (neighbors.get("amplicon") or {}) if g not in genomes]
    if unknown:
        die(f"neighbors.json names {len(unknown)} genomes absent from panel_eligible.json "
            f"(e.g. {unknown[:5]}): mismatched artefacts", 3)
    disable = {s.strip().lower() for s in (args.disable_edges or "").split(",") if s.strip()}
    if disable:
        log(f"[warn] edge types DISABLED: {sorted(disable)} — this is a gate demonstration, "
            f"not a valid plan")

    fanout = load_json(need(args.fanout, "amplicon_md5_fanout.json"))
    truth = load_json(need(args.truth, "benchmark_truth.json"))
    cells_by_genome = collections.defaultdict(list)
    for md5, cells in fanout.items():
        for cell in cells:
            gid, _, region = cell.partition("__")
            if gid in genomes:
                cells_by_genome[gid].append((md5, region))
    cells_in_truth = collections.defaultdict(list)
    for cell, rec in truth.items():
        gid = rec.get("src_genome_id")
        if gid in genomes:
            cells_in_truth[gid].append(cell)

    n_elig = len(genomes)
    phyla = [p for p, _ in collections.Counter(
        g["phylum"] for g in genomes.values()).most_common(5)]
    comps, census = build_graph(genomes, neighbors, cells_by_genome, args.rank, disable)
    n_single = sum(1 for c in comps if len(c) == 1)
    largest = max(len(c) for c in comps)
    log(f"rank={args.rank}: {len(comps)} components ({n_single} singletons), "
        f"largest {largest} ({largest / n_elig:.1%}); edge census {dict(census)}")

    block_mode = args.rank in ("family", "order")
    if not block_mode:
        cap = 1.5 * n_elig / args.k
        if largest > cap:
            die(f"largest component {largest} exceeds 1.5*n_eligible/k = {cap:.1f}: the design "
                f"cannot be {args.k}-folded at rank '{args.rank}'. Use --rank with blocks, or "
                f"tell the orchestrator — do NOT silently re-fold.", 2)
        bins = pack_kfold(comps, genomes, args.k, phyla)
    else:
        target = args.block_target or (2000 if args.rank == "family" else 400)
        target = int(min(target, 0.5 * n_elig))
        bins = pack_blocks(comps, genomes, args.blocks, target, phyla, FROZEN["seed"])
        log(f"block mode: {args.blocks} holdout blocks, ~{target} test genomes total")

    group_id = {}
    groups = {}
    for i, comp in enumerate(comps):
        gid = f"grp{i:05d}"
        groups[gid] = comp
        for g in comp:
            group_id[g] = gid

    folds, preflight = [], {}
    gated = ["md5_shared_with_reference", "amplicon_identity_ge_tau", "full16s_identity_ge_tau"]
    if args.rank != "md5only":
        gated += ["taxon_id_shared", "genus_shared"]
    failures = []
    for i, bin_ in enumerate(bins):
        test_ids = sorted(g for comp in bin_ for g in comp)
        test_groups = sorted({group_id[g] for g in test_ids})
        ref_ids = sorted(set(genomes) - set(test_ids))
        test_cells = sorted(c for g in test_ids for c in cells_in_truth.get(g, []))
        test_md5s = sorted({md5 for g in test_ids for md5, _r in cells_by_genome.get(g, [])})
        dom = collections.Counter(genomes[g]["domain"] for g in test_ids)
        ph = collections.Counter(genomes[g]["phylum"] for g in test_ids)
        n_t = max(1, len(test_ids))
        folds.append({
            "fold": i, "test_groups": test_groups, "test_genomes": test_ids,
            "reference_genomes": ref_ids, "test_cells": test_cells, "test_md5s": test_md5s,
            "n_test_genomes": len(test_ids), "n_reference_genomes": len(ref_ids),
            "n_reference_bac": sum(1 for g in ref_ids if domain_key(genomes[g]["domain"]) == "bac"),
            "n_reference_arc": sum(1 for g in ref_ids if domain_key(genomes[g]["domain"]) == "arc"),
            "domain_frac": {d: round(dom.get(d, 0) / n_t, 4) for d in ("Bacteria", "Archaea")},
            "top5_phylum_frac": {p: round(ph.get(p, 0) / n_t, 4) for p in phyla},
        })
        counters = leakage_preflight(test_ids, ref_ids, genomes, neighbors, cells_by_genome)
        preflight[f"fold_{i}"] = counters
        bad = {k: v for k, v in counters.items() if k in gated and v}
        if bad:
            failures.append((i, bad))
        log(f"fold {i}: {len(test_ids)} test / {len(ref_ids)} reference genomes, "
            f"{len(test_cells)} cells, {len(test_md5s)} md5s | leakage {counters}")

    name = "fold_plan.json" if args.rank == "genus" else f"fold_plan_{args.rank}.json"
    path = os.path.join(out, name)
    write_json(path, {
        "meta": {"frozen": FROZEN, "rank_blocked": args.rank,
                 "k": (args.blocks if block_mode else args.k),
                 "mode": "blocks" if block_mode else "kfold",
                 "n_eligible": n_elig, "n_components": len(comps), "n_singletons": n_single,
                 "largest_component": largest, "edge_census": {k: int(v) for k, v in census.items()},
                 "disabled_edges": sorted(disable),
                 "gated_counters": gated, "generated": now_iso(), "git_commit": git_commit()},
        "groups": groups, "genome_group": group_id, "folds": folds,
        "leakage_preflight": preflight,
    }, indent=1)
    log(f"wrote {path}")

    if failures:
        die(f"leakage_preflight non-zero for {len(failures)} fold(s): {failures[:5]} — the "
            f"similarity graph is wrong at rank '{args.rank}'; both methods would see a free "
            f"near-exact neighbour of the held-out organism", 2)
    log("leakage_preflight: every gated counter is 0")


# ======================================================================================
def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", required=True, choices=["panel", "prepass", "plan"])
    ap.add_argument("--outdir", default=DEFAULT_OUTDIR)
    ap.add_argument("--force", action="store_true", help="rebuild artefacts that already exist")
    # --- panel ---
    ap.add_argument("--selection", default=DEF_SELECTION)
    ap.add_argument("--functions", default=DEF_FUNCTIONS)
    ap.add_argument("--md5-seq", default=DEF_MD5_SEQ)
    ap.add_argument("--limit-groups", type=int, default=0,
                    help="SMOKE: sample G connected components of the FREE edge graph "
                         "E1(genus) UNION E2(shared amplicon md5) and take ALL their members")
    ap.add_argument("--limit-genomes", type=int, default=0,
                    help="SMOKE: top up to N genomes with singleton components")
    ap.add_argument("--label-prefix", default=FROZEN["label_prefix"],
                    help='tip-label prefix; set to "" ONLY to prove control C3 fires')
    ap.add_argument("--threads", type=int, default=32)
    # --- prepass ---
    ap.add_argument("--amplicons", default=DEF_AMPLICONS)
    ap.add_argument("--fanout", default=DEF_FANOUT)
    ap.add_argument("--truth", default=DEF_TRUTH)
    ap.add_argument("--k2", type=int, default=400)
    ap.add_argument("--gpu-cand", type=int, default=600)
    ap.add_argument("--prefilter-k", type=int, default=1500)
    ap.add_argument("--inject-error", type=float, default=0.0)
    ap.add_argument("--inject-indel-frac", type=float, default=0.1)
    ap.add_argument("--inject-seed", type=int, default=FROZEN["seed"])
    ap.add_argument("--force-backend", default="auto",
                    choices=["auto", "cuda", "metal", "cpu"])
    # --- plan ---
    ap.add_argument("--rank", default="genus", choices=["genus", "family", "order", "md5only"])
    ap.add_argument("--k", type=int, default=FROZEN["k_folds"])
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--block-target", type=int, default=0,
                    help="total test genomes across blocks (0 = 2000 family / 400 order)")
    ap.add_argument("--disable-edges", default="",
                    help="comma list of e1,e2,e3,e4 to DISABLE — only to demonstrate that the "
                         "leakage gate fires; never use in a real plan")
    args = ap.parse_args()
    # Absolutise every path BEFORE use: align_hits.py is a subprocess run with cwd=REPO and
    # hmmalign runs through `bash -lc`, so a relative --outdir would resolve against the
    # wrong directory in those children (and silently in only some of them).
    for attr in ("outdir", "selection", "functions", "md5_seq", "amplicons", "fanout", "truth"):
        setattr(args, attr, os.path.abspath(os.path.expanduser(getattr(args, attr))))

    t0 = time.perf_counter()
    if args.phase == "panel":
        phase_panel(args)
    elif args.phase == "prepass":
        phase_prepass(args)
    else:
        phase_plan(args)
    log(f"--phase {args.phase} done in {time.perf_counter() - t0:.1f}s")


if __name__ == "__main__":
    main()
