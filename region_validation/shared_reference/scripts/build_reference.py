#!/usr/bin/env python
"""Build a PICRUSt2 custom reference *de novo* from exactly one fold's survivors.

WHY THIS COMPONENT EXISTS
-------------------------
The shared-reference k-fold benchmark asks prFBA and PICRUSt2 to predict the EC
content of held-out genomes *from the same reference genome set*, using one
annotation (the BV-BRC EC projection) as both methods' gene content and the
answer key.  prFBA's side of that is cheap: restrict its alignment database.
PICRUSt2's side is not: its reference is four coupled artefacts (aligned FASTA,
tree, HMM, substitution model) plus a trait table, and every shortcut for
producing them leaks the held-out organism back in.  This script is the honest
construction.  Four measured facts drive its design, all from
``picrust2_comparison/data/kfold_audit/audit_findings.json``:

  1. **No pruning (D3).**  ``place_seqs.py`` runs ``hmmalign --mapali <ref.fna>``,
     and hmmalign enforces an HMM/MSA checksum: feeding a pruned ``.fna`` with
     the old ``.hmm`` exits 1 ("isn't same as the one HMM came from").  Worse, a
     pruned tree still carries the held-out sequence's influence in the surviving
     branch lengths -- an information channel prFBA has no counterpart to.  So
     every fold reference is rebuilt from scratch from survivors only, and the
     only place this script ever prunes a tree is the recipe-calibration control,
     where pruning is the thing being measured.

  2. **Match columns only (§5.1).**  ``hmmalign`` against the stock 16S HMM
     returns a 4,934-column alignment at 69% gaps because it keeps insert
     columns.  Keeping only the ``#=GC RF`` match columns gives 1,529 columns at
     2.4% gaps -- the geometry the stock reference has.  The 0.10 gap-fraction
     and 2,000-column aborts below exist to catch a regression into that 69%
     failure mode.

  3. **The ``.model`` must be derived, never copied.**  epa-ng reports every
     parameter as ``(user)`` and does not validate the range suffix, so copying
     the stock ``.model`` onto a FastTree tree is a *silent* misspecification
     (FastTree's branch lengths run ~0.42x the RAxML tree's).  We reconstruct it
     from FastTree's own output: the final ``GTR_Opt freq ... rates ...`` line of
     ``-log`` (rates renormalised so gt = 1.0), and ``alpha`` from the
     ``Gamma(20) LogLk = ... alpha = ...`` line of **stderr**, which is not in
     the log.  Also measured: epa-ng accepts the ``, noname = 1-W`` suffix only
     when the model is handed over as a *file* (as a bare string it fails), and
     rejects a file *without* the suffix -- so the suffix is mandatory.

  4. **The ``g_`` prefix is load-bearing.**  ``castor_nsti.R`` reads tip labels
     with ``read.table(...)$V1``, which types all-numeric BV-BRC ids as numeric
     and silently truncates ``168471.10`` to ``168471.1``; 3.61% of the panel's
     ids end in a trailing zero, and 15 of 240 reference genomes became spurious
     *study* tips in the smoke run.  Labels therefore arrive here already
     prefixed (``panel_eligible.json:genomes.<id>.label``) and this script never
     strips them.

WHAT IT GUARANTEES
------------------
* ``NAME/NAME.{fna,tre,hmm,model}`` with directory basename == file basename,
  which is what ``picrust2.place_seqs.identify_ref_files`` requires.
* ``NAME.fna`` is an ALIGNED, match-column-only, uppercase FASTA with ``-`` gaps.
* ``set_equality_ok`` (control C4): ``set(fna ids) == set(tree tips) ==
  set(traits_EC rows) == set(traits_16S rows)``.  One differing id turns a
  reference genome into a study sequence, so this is asserted at build time and
  recorded in the manifest for the runner to re-assert.
* ``traits_EC.txt.gz`` survives ``hsp.py -m mp``.  castor's max-parsimony crashes
  with ``Error in matrix(0, nrow=Nclades, ncol=Nstates-common_state): invalid
  'ncol' value`` on all-zero columns and on 9 of 163 singleton columns of a
  binarised table.  §5.4's screen is mechanical: drop prevalence 0 and n, try;
  on exactly that castor error additionally drop prevalence 1 and n-1 and retry
  once; abort on a second failure.  ``dropped_reason`` records which rule fired
  (``zero`` = prevalence 0, ``full`` = prevalence n, ``singleton`` = prevalence 1,
  ``crash`` = prevalence n-1 -- the last two only ever on the retry).  Every
  dropped EC stays in the truth and is scored as predicted-0 by the scorer,
  never as unscoreable.
* Densities are *nested*: d025 is drawn from d050's survivors, which are drawn
  from d100's.  (A per-density independent seed, read literally from §5.2, would
  give three unrelated stratified draws and break the nesting the density sweep
  needs.)
* Resumable: a directory whose ``manifest.json`` exists is skipped unless
  ``--force``; within a build, FastTree is skipped if its tree, log and stderr
  are all present.

Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.  Every
PICRUSt2 / HMMER / epa-ng call goes through
``bash -lc 'source ~/opt/picrust2/activate.sh && ...'`` because that activates a
*different* Python; nothing from ``picrust2`` is ever imported here.  ete3 lives
only in the PICRUSt2 environment, so the four tree operations that need it (tip
list, ``prune(preserve_branch_length=True)``, mean pairwise tip distance,
genus/family monophyly) go through one JSON-in/JSON-out subprocess bridge.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import random
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.abspath(os.path.join(HERE, ".."))
RV = os.path.abspath(os.path.join(SR, ".."))
REPO = os.path.abspath(os.path.join(RV, ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))

DEFAULT_STOCK_REF = os.path.expanduser(
    "~/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files")
DEFAULT_FASTTREE = os.path.expanduser("~/opt/picrust2/bin/FastTree")
ACTIVATE = os.path.expanduser("~/opt/picrust2/activate.sh")
SELECTION_JSON = os.path.join(RV, "data", "selection_10k.json")
FUNCTIONS_JSON = os.path.join(
    RV, "picrust2_comparison", "data", "genome_function_sets.json")
KFOLD_AUDIT = os.path.join(RV, "picrust2_comparison", "data", "kfold_audit")

# ---------------------------------------------------------------------------
# §1 frozen configuration -- recorded verbatim into every output JSON.
# Nothing in this block may change after the first fold is scored.
# ---------------------------------------------------------------------------
FROZEN = {
    "tau_amplicon": 0.987,
    "tau_full16s": 0.987,
    "k_folds": 10,
    "seed": 1729,
    "prfba": {
        "k2": 200, "gpu_cand": 400, "prefilter_k": 1000, "min_ref_len": 200,
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

# Alignment geometry gates (§5.1 step 3).
MAX_MSA_GAP_FRAC = 0.10
MAX_MSA_WIDTH = 2000
# C5 / mean_pairwise_tip_distance are reported over a fixed seeded tip sample.
MPD_SAMPLE = 200
# check_monophyly is O(n) per group; cap the number of groups tested so a
# 8.6k-tip / 2k-genus build does not spend minutes on a descriptive.
MONOPHYLY_MAX_GROUPS = 500
# Stock-aligned V4 window: INCLUSIVE column pair of bac_ref.fna (width 1578),
# verbatim from kfold_audit/v4_window.json.  Slicing it as [lo:hi+1] reproduces
# the audit's q_stock.fna byte for byte, so the recipe control is a reproduction
# of the measured -0.0016 paired dF1, not a fresh experiment.
STOCK_V4_WINDOW = (529, 823)
PROBE_N = 3
PROBE_LEN = 250
CASTOR_NCOL_ERR = "invalid 'ncol' value"

DOMAIN_HMM = {
    "bac": ("bacteria", "bac_ref", "bac_ref.hmm"),
    "arc": ("archaea", "arc_ref", "arc_ref.hmm"),
}
DOMAIN_OF_NAME = {"Bacteria": "bac", "Archaea": "arc"}


class BuildError(RuntimeError):
    """Abort condition with exit code 2 (§5 gates)."""


# ---------------------------------------------------------------------------
# small utilities
# ---------------------------------------------------------------------------
def now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"],
                              capture_output=True, text=True, timeout=30
                              ).stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def log(msg: str) -> None:
    print(f"[build_reference {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def run_plain(cmd: list[str], stdout=None, stderr=None, cwd=None) -> None:
    """Run a command that needs no PICRUSt2 environment."""
    r = subprocess.run(cmd, stdout=stdout, stderr=stderr, cwd=cwd)
    if r.returncode != 0:
        raise BuildError(f"command failed ({r.returncode}): {' '.join(cmd)}")


def run_p2(cmd: str, logpath: str | None = None, check: bool = True):
    """Run *cmd* inside the PICRUSt2 environment.

    PICRUSt2's activate.sh puts a different Python on PATH; nothing from this
    process's interpreter is used by the child.  Returns (rc, combined_output).
    """
    full = f"source {ACTIVATE} && {cmd}"
    r = subprocess.run(["bash", "-lc", full], capture_output=True, text=True)
    out = (r.stdout or "") + (r.stderr or "")
    if logpath:
        with open(logpath, "w") as fh:
            fh.write(f"$ {cmd}\n\n{out}")
    if check and r.returncode != 0:
        tail = "\n".join(out.strip().splitlines()[-25:])
        raise BuildError(f"PICRUSt2-env command failed ({r.returncode}): "
                         f"{cmd}\n{tail}")
    return r.returncode, out


def read_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, str] = {}
    cur = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def write_fasta(seqs: dict[str, str], path: str, order=None) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        for sid in (order if order is not None else seqs):
            fh.write(f">{sid}\n{seqs[sid]}\n")
    os.replace(tmp, path)


def read_stockholm(path: str):
    """Return (ordered seq dict, GC annotation dict).  Handles wrapped blocks."""
    seq: dict[str, str] = {}
    gc: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("//"):
                continue
            if line.startswith("#=GC"):
                p = line.split()
                if len(p) >= 3:
                    gc[p[1]] = gc.get(p[1], "") + p[2]
                continue
            if line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 2:
                continue
            seq[p[0]] = seq.get(p[0], "") + p[1]
    return seq, gc


def sha256_labels(labels) -> str:
    joined = "\n".join(sorted(labels))
    return hashlib.sha256(joined.encode()).hexdigest()


def degap(s: str) -> str:
    return s.replace("-", "").replace(".", "").upper()


# ---------------------------------------------------------------------------
# ete3 bridge: the only tree code, and it lives in the PICRUSt2 environment
# ---------------------------------------------------------------------------
_TREEOPS = r'''
import json, random, sys, itertools
from ete3 import Tree
req = json.load(sys.stdin)
t = Tree(req["tree"], format=1, quoted_node_names=True)
out = {}

def mpd(tree, sample):
    """Mean pairwise patristic tip distance over *sample*.

    Pruning a copy to the sample with preserve_branch_length=True leaves every
    pairwise distance unchanged, so this is exact and avoids all-pairs traversal
    of a 27k-tip tree."""
    sub = tree.copy("newick")
    sub.prune(sample, preserve_branch_length=True)
    nodes = {n.name: n for n in sub.get_leaves()}
    tot = 0.0; n = 0
    names = [s for s in sample if s in nodes]
    for a, b in itertools.combinations(names, 2):
        tot += nodes[a].get_distance(nodes[b]); n += 1
    return (tot / n) if n else 0.0

for op in req["ops"]:
    if op == "tips":
        out["tips"] = [l.name for l in t.get_leaves()]
    elif op == "mpd":
        s = req["mpd_sample"]
        out["mpd"] = mpd(t, s) if len(s) >= 2 else 0.0
    elif op == "monophyly":
        res = {}
        for key, groups in req["groups"].items():
            ok = 0; tested = 0
            for gname, tips in groups.items():
                tested += 1
                try:
                    mono = t.check_monophyly(values=tips, target_attr="name",
                                             unrooted=True)[0]
                except Exception:
                    mono = False
                ok += 1 if mono else 0
            res[key] = {"frac": (ok / tested) if tested else 0.0,
                        "n_tested": tested}
        out["monophyly"] = res
    elif op == "prune":
        t.prune(req["prune_keep"], preserve_branch_length=True)
        t.write(outfile=req["prune_out"], format=1)
        out["pruned_tips"] = [l.name for l in t.get_leaves()]
json.dump(out, sys.stdout)
'''


def _run_treeops(tree: str, ops: list[str], **kw) -> dict:
    """JSON-in/JSON-out ete3 bridge into the PICRUSt2 environment.

    ete3 is not installed in py_venv (and installing it there is not this
    script's call to make), so the four tree operations that need it run as a
    throwaway script under ``~/opt/picrust2/activate.sh``.  The script is
    written beside the tree, so 60 concurrent builds never collide.
    """
    req = {"tree": tree, "ops": ops}
    req.update(kw)
    script = os.path.join(
        os.path.dirname(os.path.abspath(tree)),
        "_treeops_%s_%d.py" % (os.path.basename(tree), os.getpid()))
    with open(script, "w") as fh:
        fh.write(_TREEOPS)
    try:
        full = f"source {ACTIVATE} && python {script}"
        r = subprocess.run(["bash", "-lc", full], input=json.dumps(req),
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise BuildError(f"treeops {ops} failed: {r.stderr.strip()[-800:]}")
        try:
            return json.loads(r.stdout)
        except json.JSONDecodeError:
            raise BuildError(f"treeops {ops} produced no JSON: "
                             f"{r.stdout[-400:]} / {r.stderr[-400:]}")
    finally:
        if os.path.exists(script):
            os.remove(script)


# ---------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------
def load_panel(path: str) -> dict:
    with open(path) as fh:
        panel = json.load(fh)
    if "genomes" not in panel:
        raise BuildError(f"{path} has no 'genomes' block")
    frozen = panel.get("meta", {}).get("frozen")
    if frozen and frozen != FROZEN:
        diff = sorted(k for k in set(frozen) | set(FROZEN)
                      if frozen.get(k) != FROZEN.get(k))
        raise BuildError(
            "panel_eligible.json meta.frozen disagrees with this script's "
            f"frozen block on: {diff}.  §1 may not change after the first fold "
            "is scored; reconcile before building.")
    return panel


def load_ec_universe(path: str):
    with open(path) as fh:
        u = json.load(fh)
    if not isinstance(u, dict) or "ec" not in u or "ec_index" not in u:
        raise BuildError(
            f"{path} is not the §4.2 panel_ec_universe format "
            "({'meta','ec','prevalence','ec_index'}) -- component 1 writes it")
    order = sorted(u["ec"], key=lambda e: u["ec_index"][e])
    return order


def load_ec_sets(panel_path: str, genomes, explicit: str | None,
                 functions: str) -> dict[str, set]:
    """Per-genome EC sets, from the cheapest available source.

    ``panel_eligible.json`` carries only ``n_ec``, and
    ``genome_function_sets.json`` is 1.7 GB -- re-parsing it in each of 60
    concurrent builds is not viable.  Precedence: ``--ec-sets`` >
    ``plan_folds.py``'s own ``_cache_panel_ec.json`` (same shape, and its ECs are
    already ``norm_ec``-normalised, so the trait vocabulary matches
    ``panel_ec_universe.json`` exactly) > a cache beside the panel > the audit's
    panel subset > derive from ``--functions`` and cache.
    """
    need = set(genomes)
    panel_dir = os.path.dirname(os.path.abspath(panel_path))
    cache = os.path.join(panel_dir, "panel_ec_sets.json")
    for cand in [explicit,
                 os.path.join(panel_dir, "_cache_panel_ec.json"),
                 cache,
                 os.path.join(KFOLD_AUDIT, "panel_ec_sets.json")]:
        if cand and os.path.exists(cand):
            with open(cand) as fh:
                raw = json.load(fh)
            got = {g: set(raw[g]["ec"]) for g in need if g in raw}
            if len(got) == len(need):
                log(f"EC sets from {cand} ({len(got)} genomes)")
                return got
            if cand == explicit:
                raise BuildError(
                    f"--ec-sets {cand} is missing {len(need) - len(got)} of "
                    f"{len(need)} required genomes")
            log(f"{cand} covers only {len(got)}/{len(need)}; falling through")
    log(f"deriving EC sets from {functions} (large, cached to {cache})")
    if not os.path.exists(functions):
        raise BuildError(f"no EC-set source available; {functions} missing")
    with open(panel_path) as fh:
        allg = set(json.load(fh)["genomes"])
    sub = {}
    try:
        import ijson  # streams the 1.7 GB file instead of materialising it
        with open(functions, "rb") as fh:
            for gid, rec in ijson.kvitems(fh, "", use_float=True):
                if gid in allg:
                    sub[gid] = {"ec": sorted({e for e in (rec.get("ec") or [])
                                              if "-" not in e})}
    except ImportError:
        with open(functions) as fh:
            gfs = json.load(fh)
        sub = {g: {"ec": sorted({e for e in gfs[g]["ec"] if "-" not in e})}
               for g in allg if g in gfs}
    tmp = cache + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(sub, fh)
    os.replace(tmp, cache)
    missing = need - set(sub)
    if missing:
        raise BuildError(f"{len(missing)} reference genomes absent from "
                         f"{functions} (e.g. {sorted(missing)[:3]})")
    return {g: set(sub[g]["ec"]) for g in need}


def load_marker_counts(panel: dict, selection_path: str) -> dict[str, int]:
    """16S copy number per genome, floored at 1 (§5.5).

    ``plan_folds.py`` copies ``n_full_length_16s`` into ``panel_eligible.json``,
    so prefer that (it is the eligibility-filtered view both methods see);
    ``selection_10k.json`` is the fallback for a panel written without it.
    """
    out = {g: max(1, int(r.get("n_full_length_16s") or 1))
           for g, r in panel["genomes"].items()
           if r.get("n_full_length_16s") is not None}
    if len(out) == len(panel["genomes"]):
        return out
    with open(selection_path) as fh:
        sel = json.load(fh)
    for r in sel["genomes"]:
        out.setdefault(r["genome_id"],
                       max(1, int(r.get("n_full_length_16s") or 1)))
    return out


# ---------------------------------------------------------------------------
# density subsampling (D5, nested)
# ---------------------------------------------------------------------------
def density_seed(fold: int, frac: float) -> int:
    return FROZEN["seed"] * 1000 + fold * 10 + int(round(frac * 100))


def stratified_subsample(ids: list[str], stratum_of, target: int,
                         rng: random.Random) -> list[str]:
    """Proportional-by-stratum sample of exactly *target* ids (deterministic).

    Largest-remainder apportionment so the quotas sum to *target* exactly;
    strata and members are both visited in sorted order before shuffling, so the
    result depends only on the seed.
    """
    if target >= len(ids):
        return sorted(ids)
    buckets: dict[str, list[str]] = {}
    for g in sorted(ids):
        buckets.setdefault(stratum_of(g) or "_unknown", []).append(g)
    keys = sorted(buckets)
    exact = {k: len(buckets[k]) * target / len(ids) for k in keys}
    quota = {k: int(exact[k]) for k in keys}
    short = target - sum(quota.values())
    for k in sorted(keys, key=lambda k: (-(exact[k] - quota[k]), k))[:short]:
        quota[k] += 1
    out: list[str] = []
    for k in keys:
        pool = list(buckets[k])
        rng.shuffle(pool)
        out.extend(pool[:min(quota[k], len(pool))])
    # Rounding can still leave a deficit when a quota exceeded its stratum.
    if len(out) < target:
        rest = [g for g in sorted(set(ids) - set(out))]
        rng.shuffle(rest)
        out.extend(rest[:target - len(out)])
    return sorted(out)


def nested_density_set(full: list[str], stratum_of, fold: int,
                       frac: float) -> list[str]:
    """d100 -> d050 -> d025, each drawn from the previous (nesting is required).

    §5.2's seed formula is per-density, which read literally yields three
    independent stratified draws -- not the nested subsets the density sweep
    needs.  Chaining the draws keeps both the seed formula and the nesting.
    """
    cur = sorted(full)
    n_full = len(cur)
    for f in sorted([x for x in FROZEN["density_fractions"] if x < 1.0],
                    reverse=True):
        if f < frac - 1e-9:
            continue
        target = int(round(f * n_full))
        cur = stratified_subsample(cur, stratum_of, target,
                                   random.Random(density_seed(fold, f)))
        if abs(f - frac) < 1e-9:
            return cur
    if abs(frac - 1.0) < 1e-9:
        return cur
    raise BuildError(f"--density {frac} is not one of "
                     f"{FROZEN['density_fractions']}")


# ---------------------------------------------------------------------------
# §5.1 alignment recipe
# ---------------------------------------------------------------------------
def hmm_for(stock_ref: str, domain: str) -> str:
    sub, refdir, fn = DOMAIN_HMM[domain]
    p = os.path.join(stock_ref, sub, refdir, fn)
    if not os.path.exists(p):
        raise BuildError(f"stock HMM missing: {p}")
    return p


def align_match_columns(unaligned: str, hmm: str, sto: str, out_fna: str,
                        order: list[str]):
    """hmmalign -> keep only #=GC RF match columns -> aligned FASTA.

    Insert columns are what produce the measured 4,934-wide / 69%-gap alignment;
    the RF mask recovers the 1,529-wide / 2.4%-gap geometry of the stock
    reference.  A lowercase/dot heuristic is computed as an independent check
    that the RF mask really is the match mask.
    """
    run_p2(f"hmmalign --trim --dna --outformat Stockholm -o {sto} {hmm} "
           f"{unaligned}")
    seqs, gc = read_stockholm(sto)
    if "RF" not in gc:
        raise BuildError(f"{sto} has no '#=GC RF' line; cannot identify match "
                         "columns")
    rf = gc["RF"]
    widths = {len(s) for s in seqs.values()} | {len(rf)}
    if len(widths) != 1:
        raise BuildError(f"ragged Stockholm in {sto}: widths {sorted(widths)}")
    keep = [j for j, c in enumerate(rf) if c != "."]
    # Independent check: a match column is never all-dot/all-lowercase.
    cols = list(zip(*[seqs[i] for i in seqs]))
    heur = [j for j, col in enumerate(cols)
            if not all((c == ".") or c.islower() for c in col)]
    if heur != keep:
        raise BuildError(
            f"RF mask ({len(keep)} cols) disagrees with the insert-column "
            f"heuristic ({len(heur)} cols) in {sto}")
    aligned = {}
    gapfracs = []
    for sid, s in seqs.items():
        t = "".join(s[j] for j in keep).upper().replace(".", "-")
        aligned[sid] = t
        gapfracs.append(t.count("-") / len(t))
    width = len(keep)
    mean_gap = sum(gapfracs) / len(gapfracs)
    if width > MAX_MSA_WIDTH:
        raise BuildError(f"MSA width {width} > {MAX_MSA_WIDTH}: insert columns "
                         "were kept (§5.1 step 3)")
    if mean_gap > MAX_MSA_GAP_FRAC:
        raise BuildError(f"MSA mean gap fraction {mean_gap:.4f} > "
                         f"{MAX_MSA_GAP_FRAC}: insert columns were kept")
    missing = [g for g in order if g not in aligned]
    if missing:
        raise BuildError(f"hmmalign dropped {len(missing)} sequences, "
                         f"e.g. {missing[:3]}")
    write_fasta(aligned, out_fna, order=order)
    return width, mean_gap


def hmmbuild(fna: str, hmm_out: str) -> int:
    """Per-fold HMM from exactly this .fna (hmmalign --mapali checks the cksum)."""
    run_p2(f"hmmbuild --dna --informat afa {hmm_out} {fna}",
           logpath=hmm_out + ".log")
    leng = 0
    with open(hmm_out) as fh:
        for line in fh:
            if line.startswith("LENG"):
                leng = int(line.split()[1])
                break
            if line.startswith("HMM "):
                break
    return leng


def resolve_nonroot_polytomies(tre: str) -> int:
    """Binarise non-root internal nodes; return how many were resolved.

    epa-ng ABORTS on a tree containing multifurcations ("Input Tree contains
    multifurcations (polytomies)! Please resolve the tree fully."), which at full
    density kills every bac build after FastTree has already run.  FastTree emits a
    polytomy wherever it cannot separate near-identical 16S, and the stock PICRUSt2
    reference is fully binary, so this is our recipe's gap, not epa-ng's.

    A trifurcating ROOT is how an unrooted Newick is written and epa-ng accepts it
    (every reference that passed here has exactly one such node and no other), so
    the root is left alone.  Each remaining polytomy is split by balanced pairing
    under zero-length internal branches: the resolution is arbitrary but carries no
    length, so no pairwise tip distance changes and C5 is unaffected.

    Idempotent and byte-stable: a tree that needs nothing is never rewritten.
    """
    from Bio import Phylo
    from Bio.Phylo.BaseTree import Clade

    tree = Phylo.read(tre, "newick")
    resolved = 0
    for clade in tree.find_clades(order="level"):
        if clade is tree.root or len(clade.clades) <= 2:
            continue
        kids = list(clade.clades)
        while len(kids) > 2:                      # balanced -> depth log2(n)
            nxt = []
            for i in range(0, len(kids), 2):
                if i + 1 < len(kids):
                    nxt.append(Clade(branch_length=0.0,
                                     clades=[kids[i], kids[i + 1]]))
                else:
                    nxt.append(kids[i])
            kids = nxt
        clade.clades = kids
        resolved += 1
    if resolved:
        tmp = tre + ".tmp"
        Phylo.write(tree, tmp, "newick", format_branch_length="%1.6f")
        os.replace(tmp, tre)
        log(f"  resolved {resolved} non-root polytomy(ies) for epa-ng")
    return resolved


def fasttree(fna: str, tre: str, logf: str, errf: str, binary: str,
             force: bool) -> float:
    if (not force and os.path.exists(tre) and os.path.getsize(tre) > 0
            and os.path.exists(logf) and os.path.exists(errf)):
        log(f"  FastTree: reusing {os.path.basename(tre)}")
        resolve_nonroot_polytomies(tre)
        return 0.0
    t0 = time.time()
    cmd = [binary] + FROZEN["fasttree_args"] + ["-log", logf, fna]
    with open(tre + ".tmp", "w") as out, open(errf, "w") as err:
        run_plain(cmd, stdout=out, stderr=err)
    os.replace(tre + ".tmp", tre)
    resolve_nonroot_polytomies(tre)
    return time.time() - t0


_GTR_OPT = re.compile(
    r"^GTR_Opt\s+freq\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+rates"
    r"\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")
_ALPHA = re.compile(r"Gamma\(\d+\)\s+LogLk\s*=\s*\S+\s+alpha\s*=\s*([0-9.eE+-]+)")


def derive_model(logf: str, errf: str, width: int):
    """Reconstruct epa-ng's model string from FastTree's own output.

    The final ``GTR_Opt freq ... rates ...`` line of the ``-log`` carries the ML
    frequencies and the six rates in FastTree's own scale; epa-ng wants them with
    gt = 1.0, so they are divided by the sixth.  ``alpha`` appears only on
    stderr.  The ``, noname = 1-W`` suffix is mandatory for a model *file*
    (measured: epa-ng rejects a file without it) and is not range-validated.
    """
    last = None
    with open(logf) as fh:
        for line in fh:
            m = _GTR_OPT.match(line.rstrip("\n"))
            if m:
                last = m
    if last is None:
        raise BuildError(f"no 'GTR_Opt freq ... rates ...' line in {logf}; "
                         "was FastTree run with -gtr?")
    freqs = [last.group(i) for i in (1, 2, 3, 4)]
    rates = [float(last.group(i)) for i in (5, 6, 7, 8, 9, 10)]
    if rates[5] <= 0:
        raise BuildError(f"non-positive gt rate {rates[5]} in {logf}")
    norm = [r / rates[5] for r in rates]
    alpha = None
    with open(errf) as fh:
        for line in fh:
            m = _ALPHA.search(line)
            if m:
                alpha = m.group(1)
    if alpha is None:
        raise BuildError(f"no 'Gamma(N) LogLk = ... alpha = ...' line in "
                         f"{errf}; was FastTree run with -gamma?")
    model = ("GTR{" + "/".join(f"{r:.6f}" for r in norm) + "}"
             "+FU{" + "/".join(freqs) + "}"
             "+G4m{" + alpha + "}"
             f", noname = 1-{width}")
    return model, float(alpha)


# ---------------------------------------------------------------------------
# §5.3-5.5 trait and marker tables
# ---------------------------------------------------------------------------
def write_traits(path: str, header: list[str], rows: list[str],
                 values) -> None:
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt", compresslevel=6) as fh:
        fh.write("assembly\t" + "\t".join(header) + "\n")
        for label in rows:
            fh.write(label + "\t" + "\t".join(values(label)) + "\n")
    os.replace(tmp, path)


def write_marker(path: str, labels: list[str], counts: dict[str, int],
                 label_to_gid: dict[str, str]) -> int:
    def val(label):
        return [str(counts.get(label_to_gid[label], 1))]
    write_traits(path, ["16S_rRNA_Count"], labels, val)
    return len(labels)


def gz_header(path: str) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def gz_rows(path: str) -> list[str]:
    out = []
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            if line.strip():
                out.append(line.split("\t", 1)[0])
    return out


# ---------------------------------------------------------------------------
# §5.4 crash-column screen
# ---------------------------------------------------------------------------
def probe_queries(fna_aligned: str, out: str, window=None) -> int:
    """PROBE_N short slices of the first surviving reference sequences.

    ``hsp.py`` needs a tree that already contains *unknown* tips, so the screen
    cannot be run on the reference tree alone; these probes exist only to make
    ``place_seqs.py`` produce such a tree -- which doubles as the proof that
    epa-ng accepts the derived ``.model``.  Headers are single-field and
    deliberately *not* ``g_``-prefixed so they can never collide with a label.
    """
    seqs = read_fasta(fna_aligned)
    ids = list(seqs)[:PROBE_N]
    picked = {}
    for i, sid in enumerate(ids):
        s = seqs[sid]
        if window:
            lo, hi = window
            frag = degap(s[lo:hi + 1])
        else:
            lo = int(round(STOCK_V4_WINDOW[0] / 1578 * len(s)))
            frag = ""
            j = lo
            while j < len(s) and len(frag) < PROBE_LEN:
                c = s[j].upper()
                if c not in "-.":
                    frag += c
                j += 1
        frag = frag[:PROBE_LEN]
        if len(frag) >= 100:
            picked[f"probe_{i + 1}"] = frag
    if not picked:
        raise BuildError(f"could not cut probe queries from {fna_aligned}")
    write_fasta(picked, out)
    return len(picked)


def prevalence_screen(ec_order: list[str], labels: list[str],
                      ecsets: dict[str, set], label_to_gid, aggressive: bool):
    """Drop castor-crashing columns; return (kept, dropped_reason)."""
    n = len(labels)
    kept, dropped = [], {}
    for ec in ec_order:
        p = sum(1 for lab in labels if ec in ecsets[label_to_gid[lab]])
        if p == 0:
            dropped[ec] = "zero"
        elif p == n:
            dropped[ec] = "full"
        elif aggressive and p == 1:
            dropped[ec] = "singleton"
        elif aggressive and p == n - 1:
            dropped[ec] = "crash"
        else:
            kept.append(ec)
    return kept, dropped


def crash_screen(refdir: str, name: str, ec_order: list[str],
                 labels: list[str], ecsets, label_to_gid, probe_tre: str,
                 threads: int):
    """§5.4: write the trait table, prove ``hsp.py -m mp`` survives it.

    Only the exact castor ``invalid 'ncol' value`` error licenses the retry; any
    other failure aborts immediately with its log, because retrying a missing-R-
    package error as if it were a column problem would hide the real fault.
    """
    traits = os.path.join(refdir, "traits_EC.txt.gz")
    scratch = os.path.join(refdir, "hsp_screen")
    os.makedirs(scratch, exist_ok=True)
    retry_used = False
    for aggressive in (False, True):
        kept, dropped = prevalence_screen(ec_order, labels, ecsets,
                                          label_to_gid, aggressive)
        if not kept:
            raise BuildError(f"{name}: the prevalence screen emptied the trait "
                             f"vocabulary ({len(ec_order)} ECs, {len(labels)} "
                             "rows)")
        write_traits(traits, kept, labels,
                     lambda lab: ["1" if ec in ecsets[label_to_gid[lab]] else "0"
                                  for ec in kept])
        out = os.path.join(scratch, "probe_EC.tsv.gz")
        if os.path.exists(out):
            os.remove(out)
        rc, txt = run_p2(
            f"hsp.py -t {probe_tre} --observed_trait_table {traits} "
            f"-o {out} -m {FROZEN['picrust2']['hsp_method']} "
            f"-e {FROZEN['picrust2']['edge_exponent']} -n -p {threads} "
            f"--seed {FROZEN['seed']}",
            logpath=os.path.join(refdir, f"{name}.hsp_screen.log"), check=False)
        if rc == 0:
            return kept, dropped, retry_used
        if CASTOR_NCOL_ERR not in txt:
            tail = "\n".join(txt.strip().splitlines()[-25:])
            raise BuildError(
                f"{name}: hsp.py failed for a reason other than the castor "
                f"column crash -- not retrying.\n{tail}")
        if aggressive:
            break
        log(f"  crash screen: castor '{CASTOR_NCOL_ERR}' -> dropping "
            "prevalence-1 and prevalence-(n-1) columns and retrying once")
        retry_used = True
    raise BuildError(f"{name}: hsp.py -m mp still crashes after the §5.4 "
                     "retry; aborting this reference (exit 2)")


# ---------------------------------------------------------------------------
# descriptives (C13) and set equality (C4)
# ---------------------------------------------------------------------------
def tree_descriptives(tre: str, labels: list[str], panel_genomes: dict,
                      label_to_gid: dict[str, str]) -> dict:
    rng = random.Random(FROZEN["seed"])
    sample = sorted(labels)
    if len(sample) > MPD_SAMPLE:
        sample = sorted(rng.sample(sample, MPD_SAMPLE))
    groups = {}
    for rank in ("genus", "family"):
        by = {}
        for lab in labels:
            v = (panel_genomes.get(label_to_gid[lab], {}) or {}).get(rank) or ""
            if v:
                by.setdefault(v, []).append(lab)
        by = {k: v for k, v in by.items() if len(v) >= 2}
        keys = sorted(by)
        if len(keys) > MONOPHYLY_MAX_GROUPS:
            keys = sorted(random.Random(FROZEN["seed"] + 1
                                        ).sample(keys, MONOPHYLY_MAX_GROUPS))
        groups[rank] = {k: by[k] for k in keys}
    res = _run_treeops(tre, ["tips", "mpd", "monophyly"],
                       mpd_sample=sample, groups=groups)
    mono = res.get("monophyly", {})
    return {
        "tips": res["tips"],
        "mean_pairwise_tip_distance": res.get("mpd", 0.0),
        "mpd_sample_n": len(sample),
        "genus_monophyly_frac": mono.get("genus", {}).get("frac", 0.0),
        "genus_monophyly_n_tested": mono.get("genus", {}).get("n_tested", 0),
        "family_monophyly_frac": mono.get("family", {}).get("frac", 0.0),
        "family_monophyly_n_tested": mono.get("family", {}).get("n_tested", 0),
    }


def assert_set_equality(name: str, fna: str, tips, traits: str,
                        marker: str) -> bool:
    s_fna = set(read_fasta(fna))
    s_tre = set(t for t in tips if t)
    s_tr = set(gz_rows(traits))
    s_mk = set(gz_rows(marker))
    if not (s_fna == s_tre == s_tr == s_mk):
        def d(a, b):
            return sorted(a - b)[:5]
        raise BuildError(
            f"{name}: C4 set equality FAILED -- a reference genome would "
            f"silently become a study sequence. "
            f"fna\\tre={d(s_fna, s_tre)} tre\\fna={d(s_tre, s_fna)} "
            f"fna\\traits={d(s_fna, s_tr)} fna\\marker={d(s_fna, s_mk)}")
    return True


# ---------------------------------------------------------------------------
# fold / include_self / plumbing builds
# ---------------------------------------------------------------------------
def ref_name(mode: str, fold, domain: str, density: float) -> str:
    ddd = f"d{int(round(density * 100)):03d}"
    if mode == "fold":
        return f"f{int(fold):02d}_{domain}_{ddd}"
    if mode == "include_self":
        return f"incself_{domain}_{ddd}"
    if mode == "plumbing":
        return f"plumbing_{domain}_{ddd}"
    raise BuildError(f"no name rule for mode {mode}")


def build_one(args, mode: str, fold, domain: str, density: float,
              panel: dict, ec_order: list[str], ecsets, markers) -> str:
    name = ref_name(mode, fold, domain, density)
    refdir = os.path.join(args.outdir, name)
    manifest_p = os.path.join(refdir, "manifest.json")
    if os.path.exists(manifest_p) and not args.force:
        log(f"{name}: manifest.json exists -> skip (use --force to rebuild)")
        return manifest_p
    os.makedirs(refdir, exist_ok=True)
    genomes = panel["genomes"]

    # --- the surviving reference set -------------------------------------- #
    if mode == "fold":
        plan = args._plan
        try:
            fspec = next(f for f in plan["folds"] if int(f["fold"]) == int(fold))
        except StopIteration:
            raise BuildError(f"fold {fold} not in {args.fold_plan}")
        pool = list(fspec["reference_genomes"])
    else:
        # include_self / plumbing: the whole eligible panel is the reference,
        # which makes these references fold-independent by construction.
        pool = list(genomes)
    pool = [g for g in pool
            if DOMAIN_OF_NAME.get((genomes[g] or {}).get("domain")) == domain]
    if not pool:
        raise BuildError(f"{name}: no eligible {domain} genomes in the "
                         "reference set")
    seed_fold = int(fold) if (mode == "fold" and fold is not None) else 0
    refs = nested_density_set(pool, lambda g: genomes[g].get("phylum") or "",
                              seed_fold, density)
    if len(refs) < 4:
        raise BuildError(f"{name}: {len(refs)} tips is too few to build a tree "
                         "(FastTree needs >= 4); raise --density or widen the "
                         "panel")
    label_to_gid = {genomes[g]["label"]: g for g in refs}
    labels = sorted(label_to_gid)
    if len(labels) != len(refs):
        raise BuildError(f"{name}: duplicate labels among {len(refs)} genomes")

    log(f"{name}: {len(refs)} {domain} tips (density {density})")

    # --- §5.1 alignment ---------------------------------------------------- #
    panel_fna = os.path.join(args.paneldir, f"panel_{domain}.fna")
    if not os.path.exists(panel_fna):
        raise BuildError(f"missing {panel_fna} (plan_folds.py --phase panel)")
    allseq = read_fasta(panel_fna)
    unaligned = os.path.join(refdir, f"{name}.input.fna")
    missing = [lab for lab in labels if lab not in allseq]
    if missing:
        raise BuildError(f"{len(missing)} labels absent from {panel_fna}, "
                         f"e.g. {missing[:3]}")
    write_fasta({lab: allseq[lab].upper() for lab in labels}, unaligned,
                order=labels)
    fna = os.path.join(refdir, f"{name}.fna")
    sto = os.path.join(refdir, f"{name}.sto")
    if args.force or not os.path.exists(fna):
        width, mean_gap = align_match_columns(
            unaligned, hmm_for(args.stock_ref, domain), sto, fna, labels)
    else:
        aln = read_fasta(fna)
        width = len(next(iter(aln.values())))
        mean_gap = sum(s.count("-") / len(s) for s in aln.values()) / len(aln)
        log(f"  reusing {os.path.basename(fna)} (width {width})")

    hmm = os.path.join(refdir, f"{name}.hmm")
    hmm_leng = hmmbuild(fna, hmm) if (args.force or not os.path.exists(hmm)) \
        else _hmm_leng(hmm)
    tre = os.path.join(refdir, f"{name}.tre")
    logf = os.path.join(refdir, f"{name}.fasttree.log")
    errf = os.path.join(refdir, f"{name}.fasttree.stderr")
    secs = fasttree(fna, tre, logf, errf, args.fasttree_bin, args.force)
    model_string, alpha = derive_model(logf, errf, width)
    with open(os.path.join(refdir, f"{name}.model"), "w") as fh:
        fh.write(model_string + "\n")

    # --- §5.5 marker, §5.4 screened traits --------------------------------- #
    marker = os.path.join(refdir, "traits_16S.txt.gz")
    marker_rows = write_marker(marker, labels, markers, label_to_gid)

    probe = os.path.join(refdir, "probe_queries.fna")
    probe_queries(fna, probe)
    shared_probe = os.path.join(args.outdir, "probe_queries.fna")
    if not os.path.exists(shared_probe):
        try:
            shutil.copyfile(probe, shared_probe)
        except OSError:
            pass
    probe_tre = os.path.join(refdir, f"{name}.probe.tre")
    if args.force or not os.path.exists(probe_tre):
        if os.path.exists(probe_tre):
            os.remove(probe_tre)
        inter = os.path.join(refdir, "probe_int")
        shutil.rmtree(inter, ignore_errors=True)
        run_p2(f"place_seqs.py -s {probe} -r {refdir} -o {probe_tre} "
               f"-p {args.threads} "
               f"--min_align {FROZEN['picrust2']['min_align']} "
               f"--intermediate {inter} --verbose",
               logpath=os.path.join(refdir, f"{name}.probe_place.log"))
    kept_ec, dropped, retry_used = crash_screen(
        refdir, name, ec_order, labels, ecsets, label_to_gid, probe_tre,
        args.threads)
    traits = os.path.join(refdir, "traits_EC.txt.gz")

    desc = tree_descriptives(tre, labels, genomes, label_to_gid)
    set_ok = assert_set_equality(name, fna, desc["tips"], traits, marker)

    man = {
        "name": name, "mode": mode,
        # include_self / plumbing references are the whole eligible panel and
        # therefore fold-independent; --fold is recorded if supplied, else null.
        "fold": int(fold) if fold is not None else None,
        "domain": domain, "density": density,
        "n_tips": len(labels), "msa_width": width,
        "msa_mean_gap_frac": round(mean_gap, 6), "hmm_leng": hmm_leng,
        "model_string": model_string, "fasttree_alpha": alpha,
        "fasttree_seconds": round(secs, 3),
        "tip_labels_sha256": sha256_labels(labels),
        "trait_cols": len(kept_ec),
        "dropped_ec": sorted(dropped),
        "dropped_reason": dropped, "n_dropped": len(dropped),
        "retry_used": retry_used, "marker_rows": marker_rows,
        "set_equality_ok": set_ok,
        "mean_pairwise_tip_distance": desc["mean_pairwise_tip_distance"],
        "mpd_sample_n": desc["mpd_sample_n"],
        "genus_monophyly_frac": desc["genus_monophyly_frac"],
        "genus_monophyly_n_tested": desc["genus_monophyly_n_tested"],
        "family_monophyly_frac": desc["family_monophyly_frac"],
        "family_monophyly_n_tested": desc["family_monophyly_n_tested"],
        "ec_universe_n": len(ec_order),
        "reference_genomes": refs,
        "frozen": FROZEN, "generated": now_iso(), "git_commit": git_commit(),
    }
    _dump(manifest_p, man)
    log(f"{name}: OK  tips={len(labels)} width={width} "
        f"gap={mean_gap:.4f} traits={len(kept_ec)} dropped={len(dropped)} "
        f"retry={retry_used} mpd={desc['mean_pairwise_tip_distance']:.4f}")
    return manifest_p


def _hmm_leng(hmm: str) -> int:
    with open(hmm) as fh:
        for line in fh:
            if line.startswith("LENG"):
                return int(line.split()[1])
            if line.startswith("HMM "):
                break
    return 0


def _dump(path: str, obj) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=False)
        fh.write("\n")
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# §5.7 recipe-calibration control (the decisive control, C1 + C5)
# ---------------------------------------------------------------------------
def build_recipe_control(args, ec_order_unused=None) -> list[str]:
    """Three arms over the IDENTICAL retained stock genomes, traits and queries.

    A = stock RAxML tree pruned with ``preserve_branch_length=True`` plus the
    stock ``.model`` (the one place in this experiment where a reference is ever
    pruned, D15), on the stock MSA; B = our hmmbuild/FastTree/derived-model
    recipe on that same stock MSA (§5.1 steps 5-7 only -- ``bac_ref.fna`` is
    already aligned); C = the SAME retained sequences degapped and run through
    the FULL §5.1 recipe (steps 1-7: hmmalign against the stock domain HMM,
    ``#=GC RF`` match columns, per-reference hmmbuild, FastTree, derived model),
    i.e. exactly what every fold reference is built with.

    C1 is therefore a decomposition over one shared query set: dF1(C-A) is the
    whole build-recipe dialect (the gated number), dF1(B-A) the tree+model half
    with the MSA held fixed, dF1(C-B) the alignment half with the tree recipe
    held fixed.  Any |dF1(C-A)| >= 0.01 means our recipe handicaps PICRUSt2 and
    every PICRUSt2 number must carry that band.
    """
    bac = os.path.join(args.stock_ref, "bacteria")
    ref = os.path.join(bac, "bac_ref")
    stock_fna, stock_tre, stock_model = (
        os.path.join(ref, "bac_ref.fna"), os.path.join(ref, "bac_ref.tre"),
        os.path.join(ref, "bac_ref.model"))
    for p in (stock_fna, stock_tre, stock_model):
        if not os.path.exists(p):
            raise BuildError(f"stock reference file missing: {p}")
    seqs = read_fasta(stock_fna)
    widths = {len(s) for s in seqs.values()}
    if len(widths) != 1:
        raise BuildError("stock bac_ref.fna is not an alignment")
    width = widths.pop()

    kept, held = _recipe_split(args, sorted(seqs))
    log(f"recipe_control: {len(kept)} retained / {len(held)} held out "
        f"(stock alignment width {width})")

    # Shared artefacts: identical retained alignment, identical HMM.
    os.makedirs(args.outdir, exist_ok=True)
    shared_fna = os.path.join(args.outdir, "recipe_retained.fna")
    write_fasta({g: seqs[g] for g in kept}, shared_fna, order=kept)
    shared_hmm = os.path.join(args.outdir, "recipe_retained.hmm")
    if args.force or not os.path.exists(shared_hmm):
        # The stock .hmm cannot be reused: hmmalign --mapali enforces a checksum
        # against the exact .fna, and this .fna is a 960-row subset.
        hmmbuild(shared_fna, shared_hmm)
    hmm_leng = _hmm_leng(shared_hmm)

    # Queries + truth from the held-out stock rows.
    qpath = os.path.join(args.outdir, "recipe_control_queries.fna")
    lo, hi = STOCK_V4_WINDOW
    queries = {}
    for g in held:
        frag = degap(seqs[g][lo:hi + 1])
        if len(frag) >= 100:
            queries[_qlabel(g)] = frag
    if not queries:
        raise BuildError("no recipe-control queries could be cut")
    write_fasta(queries, qpath, order=sorted(queries))

    ec_hdr, ec_rows = _stock_ec(os.path.join(bac, "ec.txt.gz"),
                                set(kept) | set(held))
    nz = [j for j, _ in enumerate(ec_hdr)
          if any(ec_rows[g][j] not in ("0", "0.0") for g in kept)]
    cols = [ec_hdr[j] for j in nz]
    truth = {_qlabel(g): sorted(ec_hdr[j] for j in nz
                                if ec_rows[g][j] not in ("0", "0.0"))
             for g in held}
    _dump(os.path.join(args.outdir, "recipe_control_truth.json"), truth)
    log(f"recipe_control: {len(cols)} non-zero stock EC columns over the "
        f"retained set; {len(truth)} truth rows")

    marker_src = _stock_marker(os.path.join(bac, "16S.txt.gz"), set(kept))

    stock_gap = sum(seqs[g].count("-") / width for g in kept) / len(kept)

    out = []
    for arm in ("A", "B", "C"):
        name = f"recipe_{arm}"
        refdir = os.path.join(args.outdir, name)
        manifest_p = os.path.join(refdir, "manifest.json")
        if os.path.exists(manifest_p) and not args.force:
            log(f"{name}: manifest.json exists -> skip")
            out.append(manifest_p)
            continue
        os.makedirs(refdir, exist_ok=True)
        fna = os.path.join(refdir, f"{name}.fna")
        hmm_p = os.path.join(refdir, f"{name}.hmm")
        if arm == "C":
            # §5.1 steps 1-3 on the SAME retained residues: degap the stock
            # alignment and re-align it with our hmmalign + RF-match-column
            # recipe, so the alignment half of the dialect gets its own paired
            # comparison instead of being shared by construction.
            unaligned = os.path.join(refdir, f"{name}.input.fna")
            write_fasta({g: degap(seqs[g]).upper() for g in kept}, unaligned,
                        order=kept)
            sto = os.path.join(refdir, f"{name}.sto")
            arm_width, arm_gap = align_match_columns(
                unaligned, hmm_for(args.stock_ref, "bac"), sto, fna, kept)
            arm_hmm_leng = hmmbuild(fna, hmm_p)
            log(f"{name}: re-aligned {len(kept)} degapped stock sequences -> "
                f"width {arm_width} (stock {width}), mean gap {arm_gap:.4f}")
        else:
            shutil.copyfile(shared_fna, fna)
            shutil.copyfile(shared_hmm, hmm_p)
            arm_width, arm_gap, arm_hmm_leng = width, stock_gap, hmm_leng
        tre = os.path.join(refdir, f"{name}.tre")
        secs, alpha, model_string = 0.0, None, None
        mpd_before = mpd_after = None
        rng = random.Random(FROZEN["seed"])
        sample = sorted(kept) if len(kept) <= MPD_SAMPLE else \
            sorted(rng.sample(sorted(kept), MPD_SAMPLE))
        if arm == "A":
            # C5: the same seeded 200-tip sample before and after pruning.
            mpd_before = _run_treeops(stock_tre, ["mpd"],
                                      mpd_sample=sample)["mpd"]
            _run_treeops(stock_tre, ["prune"], prune_keep=kept, prune_out=tre)
            mpd_after = _run_treeops(tre, ["mpd"], mpd_sample=sample)["mpd"]
            ratio = (mpd_after / mpd_before) if mpd_before else 0.0
            if abs(ratio - 1.0) >= 1e-6:
                raise BuildError(
                    f"C5 FAILED: mean pairwise tip distance changed by "
                    f"{ratio:.6f}x under pruning ({mpd_before:.4f} -> "
                    f"{mpd_after:.4f}).  ete3's default (preserve_branch_length="
                    "False) yields 0.537x and would fabricate the NSTI/novelty "
                    "relationship.")
            with open(stock_model) as fh:
                model_string = fh.read().strip()
            with open(os.path.join(refdir, f"{name}.model"), "w") as fh:
                fh.write(model_string + "\n")
        else:
            logf = os.path.join(refdir, f"{name}.fasttree.log")
            errf = os.path.join(refdir, f"{name}.fasttree.stderr")
            secs = fasttree(fna, tre, logf, errf, args.fasttree_bin, args.force)
            model_string, alpha = derive_model(logf, errf, arm_width)
            with open(os.path.join(refdir, f"{name}.model"), "w") as fh:
                fh.write(model_string + "\n")
            mpd_before = mpd_after = _run_treeops(
                tre, ["mpd"], mpd_sample=sample)["mpd"]

        traits = os.path.join(refdir, "traits_EC.txt.gz")
        write_traits(traits, cols, kept,
                     lambda g: [ec_rows[g][j] for j in nz])   # counts (D12 n/a)
        marker = os.path.join(refdir, "traits_16S.txt.gz")
        write_traits(marker, ["16S_rRNA_Count"], kept,
                     lambda g: [str(marker_src.get(g, 1))])

        probe = os.path.join(refdir, "probe_queries.fna")
        # STOCK_V4_WINDOW is in stock-1578 coordinates; arm C's own MSA needs the
        # proportional cut used by every fold build.
        probe_queries(fna, probe,
                      window=(STOCK_V4_WINDOW if arm != "C" else None))
        probe_tre = os.path.join(refdir, f"{name}.probe.tre")
        if args.force or not os.path.exists(probe_tre):
            if os.path.exists(probe_tre):
                os.remove(probe_tre)
            inter = os.path.join(refdir, "probe_int")
            shutil.rmtree(inter, ignore_errors=True)
            run_p2(f"place_seqs.py -s {probe} -r {refdir} -o {probe_tre} "
                   f"-p {args.threads} "
                   f"--min_align {FROZEN['picrust2']['min_align']} "
                   f"--intermediate {inter} --verbose",
                   logpath=os.path.join(refdir, f"{name}.probe_place.log"))
        rc, txt = run_p2(
            f"hsp.py -t {probe_tre} --observed_trait_table {traits} "
            f"-o {os.path.join(refdir, 'hsp_screen_EC.tsv.gz')} "
            f"-m {FROZEN['picrust2']['hsp_method']} "
            f"-e {FROZEN['picrust2']['edge_exponent']} -n -p {args.threads} "
            f"--seed {FROZEN['seed']}",
            logpath=os.path.join(refdir, f"{name}.hsp_screen.log"), check=False)
        if rc != 0:
            tail = "\n".join(txt.strip().splitlines()[-20:])
            raise BuildError(f"{name}: hsp.py failed on the stock count table "
                             f"(measured OK for 960x3111) -- aborting\n{tail}")
        tips = _run_treeops(tre, ["tips"])["tips"]
        set_ok = assert_set_equality(name, fna, tips, traits, marker)
        man = {
            "name": name, "mode": "recipe_control", "fold": None,
            "domain": "bac", "density": 1.0, "n_tips": len(kept),
            "msa_width": arm_width, "msa_mean_gap_frac": round(arm_gap, 6),
            "hmm_leng": arm_hmm_leng, "model_string": model_string,
            "fasttree_alpha": alpha, "fasttree_seconds": round(secs, 3),
            "tip_labels_sha256": sha256_labels(kept),
            "trait_cols": len(cols), "dropped_ec": [], "dropped_reason": {},
            "n_dropped": 0, "retry_used": False, "marker_rows": len(kept),
            "set_equality_ok": set_ok,
            "mean_pairwise_tip_distance": mpd_after,
            "mean_pairwise_tip_distance_before_prune": mpd_before,
            "mean_pairwise_tip_distance_after_prune": mpd_after,
            "mpd_sample_n": len(sample),
            "genus_monophyly_frac": None, "family_monophyly_frac": None,
            "recipe_arm": arm,
            # which half of the build recipe this arm uses (C1 decomposition)
            "recipe_alignment": ("stock_ssu_align_msa" if arm != "C"
                                 else "ours_hmmalign_rf_match_columns"),
            "recipe_tree_model": ("stock_raxml_pruned_plus_stock_model"
                                  if arm == "A"
                                  else "ours_fasttree_plus_derived_model"),
            "n_holdout": len(held), "holdout": held,
            "queries": os.path.relpath(qpath, args.outdir),
            "truth": "recipe_control_truth.json",
            "reference_genomes": kept,
            "frozen": FROZEN, "generated": now_iso(),
            "git_commit": git_commit(),
        }
        _dump(manifest_p, man)
        log(f"{name}: OK  tips={len(kept)} width={arm_width} "
            f"traits={len(cols)} mpd={mpd_after:.4f}")
        out.append(manifest_p)
    return out


def _qlabel(acc: str) -> str:
    """Stock accession -> query label, reproducing kfold_audit/q_stock.fna."""
    return "qs_" + acc.replace(".", "_")


def _recipe_split(args, all_ids: list[str]):
    """Retained / held-out stock genomes.

    Reuses ``kfold_audit/stockctl_split.json`` when its sizes match the request,
    so smoke gate 5 reproduces the measured -0.0016 paired dF1 instead of
    sampling a fresh split.
    """
    split = args.stock_split
    if split and os.path.exists(split):
        with open(split) as fh:
            s = json.load(fh)
        kept, held = list(s.get("kept", [])), list(s.get("held", []))
        if (len(kept) == args.n_stock and len(held) == args.n_holdout
                and set(kept + held) <= set(all_ids)):
            log(f"recipe_control: reusing the audit split {split}")
            return kept, held
        log(f"recipe_control: {split} is {len(kept)}/{len(held)}, requested "
            f"{args.n_stock}/{args.n_holdout} -> drawing a fresh seeded split")
    rng = random.Random(FROZEN["seed"])
    pool = sorted(all_ids)
    take = args.n_stock + args.n_holdout
    if take > len(pool):
        raise BuildError(f"--n-stock + --n-holdout = {take} exceeds the "
                         f"{len(pool)} stock reference sequences")
    sub = rng.sample(pool, take)
    return sorted(sub[:args.n_stock]), sorted(sub[args.n_stock:])


def _stock_ec(path: str, want: set):
    hdr, rows = None, {}
    with gzip.open(path, "rt") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] in want:
                rows[p[0]] = p[1:]
    missing = want - set(rows)
    if missing:
        raise BuildError(f"{len(missing)} stock genomes absent from {path}, "
                         f"e.g. {sorted(missing)[:3]}")
    return hdr, rows


def _stock_marker(path: str, want: set) -> dict[str, int]:
    out = {}
    with gzip.open(path, "rt") as fh:
        fh.readline()
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if p[0] in want:
                try:
                    out[p[0]] = max(1, int(float(p[1])))
                except ValueError:
                    out[p[0]] = 1
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Build PICRUSt2 custom references for the shared-reference "
                    "k-fold benchmark (one per fold x domain x density) and the "
                    "control references.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--mode", required=True,
                    choices=["fold", "recipe_control", "include_self",
                             "plumbing"])
    ap.add_argument("--fold-plan",
                    default=os.path.join(SR, "data", "fold_plan.json"))
    ap.add_argument("--panel",
                    default=os.path.join(SR, "data", "panel_eligible.json"))
    ap.add_argument("--fold", type=int, default=None,
                    help="required for --mode fold")
    ap.add_argument("--domain", default="both", choices=["bac", "arc", "both"])
    ap.add_argument("--density", type=float, default=1.0,
                    help="one of frozen.density_fractions")
    ap.add_argument("--outdir", default=os.path.join(SR, "data", "refs"))
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--fasttree-bin", default=DEFAULT_FASTTREE)
    ap.add_argument("--stock-ref", default=DEFAULT_STOCK_REF)
    ap.add_argument("--n-stock", type=int, default=960,
                    help="recipe_control: retained stock reference tips")
    ap.add_argument("--n-holdout", type=int, default=500,
                    help="recipe_control: held-out stock genomes (queries)")
    ap.add_argument("--force", action="store_true")
    # --- additive conveniences (see module docstring / deviations) --------- #
    ap.add_argument("--ec-universe", default=None,
                    help="default <panel dir>/panel_ec_universe.json")
    ap.add_argument("--ec-sets", default=None,
                    help="per-genome EC sets cache; default <panel dir>/"
                         "panel_ec_sets.json, else derived from --functions")
    ap.add_argument("--functions", default=FUNCTIONS_JSON,
                    help="genome_function_sets.json (only read when no EC-set "
                         "cache exists)")
    ap.add_argument("--selection", default=SELECTION_JSON,
                    help="n_full_length_16s for the marker table")
    ap.add_argument("--stock-split",
                    default=os.path.join(KFOLD_AUDIT, "stockctl_split.json"),
                    help="recipe_control: reuse this retained/held split when "
                         "its sizes match --n-stock/--n-holdout")
    args = ap.parse_args(argv)

    args.outdir = os.path.abspath(args.outdir)
    os.makedirs(args.outdir, exist_ok=True)
    args.paneldir = os.path.dirname(os.path.abspath(args.panel))

    if not os.path.exists(args.fasttree_bin):
        log(f"WARNING: {args.fasttree_bin} not found; FastTree is required for "
            "every arm except recipe arm A")

    try:
        if args.mode == "recipe_control":
            build_recipe_control(args)
            return 0

        if args.mode == "fold" and args.fold is None:
            raise BuildError("--mode fold requires --fold N")
        panel = load_panel(args.panel)
        args._plan = None
        if args.mode == "fold":
            with open(args.fold_plan) as fh:
                args._plan = json.load(fh)
        ec_path = args.ec_universe or os.path.join(args.paneldir,
                                                   "panel_ec_universe.json")
        ec_order = load_ec_universe(ec_path)
        domains = ["bac", "arc"] if args.domain == "both" else [args.domain]

        # Which genomes will we need EC sets for?  Union over the domains built.
        need = set()
        for dom in domains:
            if args.mode == "fold":
                fspec = next((f for f in args._plan["folds"]
                              if int(f["fold"]) == int(args.fold)), None)
                if fspec is None:
                    raise BuildError(f"fold {args.fold} not in {args.fold_plan}")
                pool = list(fspec["reference_genomes"])
            else:
                pool = list(panel["genomes"])
            need |= {g for g in pool
                     if DOMAIN_OF_NAME.get(
                         (panel["genomes"][g] or {}).get("domain")) == dom}
        ecsets = load_ec_sets(args.panel, need, args.ec_sets, args.functions)
        markers = load_marker_counts(panel, args.selection)

        for dom in domains:
            build_one(args, args.mode, args.fold, dom, args.density, panel,
                      ec_order, ecsets, markers)
        return 0
    except BuildError as e:
        print(f"\nABORT (exit 2): {e}\n", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
