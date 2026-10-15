#!/usr/bin/env python
"""Run PICRUSt2 for ONE (fold, domain, density) against a custom shared reference.

WHY THIS EXISTS
---------------
The shared-reference k-fold design (``../SPEC.md``) requires PICRUSt2 to predict the
EC content of held-out genomes from *exactly* the reference set prFBA is restricted
to, with one annotation (the BV-BRC EC projection) serving as both methods' gene
content and the answer key.  ``picrust2_pipeline.py`` cannot be used for that:

  * it owns the domain-routing decision (D9 makes routing an explicit arm, not a
    hidden default), and it insists on an abundance table and a metagenome step that
  * would drag a synthetic 16S copy-number table into a comparison that must turn on
    gene content alone (D17);
  * it writes its outputs under names that encode its own two-reference wrapper.

So this component calls ``place_seqs.py`` and ``hsp.py`` directly and then does the
one piece of glue the downstream scorer needs.

WHAT IT GUARANTEES
------------------
1. ``hsp.py -n`` appends TWO non-trait columns to the EC table -- ``metadata_NSTI``
   (documented) and ``closest_reference_genome`` (undocumented; produced by
   ``castor_nsti.R`` and concatenated wholesale by ``picrust2.wrap_hsp``).  Left in
   place, the first makes ``compare.load_picrust2``'s ``df.to_numpy() > 0`` score NSTI
   as a predicted EC for every query with NSTI > 0, and the second makes the frame an
   object array so ``> 0`` raises TypeError outright.  Both are stripped into
   ``combined_EC_predicted.tsv.gz``; the metadata is re-emitted as
   ``combined_marker_predicted_and_nsti.tsv.gz`` with the ``best_domain`` column that
   ``compare.load_picrust2`` reads with a plain ``[]`` (and which only
   ``picrust2_pipeline.py`` would otherwise produce).  With those two files present
   ``compare.load_picrust2(outdir, ec_index)`` is reused verbatim.
2. C3 -- the tip-label round trip.  ``castor_nsti.R`` reads the known-tip list with
   ``read.table(...)$V1``, which types all-numeric BV-BRC ids as numeric and silently
   drops reference genomes whose id ends in a trailing zero out of the NSTI known-set.
   This component asserts that the SET of unknown-tip rows in ``hsp.py``'s output
   equals the SET of study sequences that actually survived ``--min_align`` -- set
   equality, not count equality, because the failure signature is one reference tip
   turning into one extra study row.  A fail exits 3.
3. C4 -- artefact set equality: ``NAME.fna`` ids == ``NAME.tre`` tips == trait-table
   rows, checked again at run time (the builder checks it at build time), so a
   reference genome can never silently become a study sequence.  A fail exits 3.
4. C7 -- placement-loss accounting.  Every md5 dropped by ``--min_align`` is listed
   with its raw and aligned length and its regions, so the scorer can remove it from
   BOTH methods' denominators rather than charging PICRUSt2 for a shared data limit.
   Reported, never gated.

DELIBERATE DEVIATIONS FROM THE LITERAL CONTRACT (all recorded in ``meta``)
-------------------------------------------------------------------------
* C3 is gated as ``set(unknown rows) == set(placed md5s)`` i.e.
  ``n_unknown_rows == n_query - n_min_align_dropped``.  The literal
  ``n_unknown_rows == n_query`` is unsatisfiable whenever ``--min_align`` drops a
  query, and C7 already removes those from both denominators.  Both counts are in
  ``meta``; ``n_min_align_dropped == 0`` makes the two definitions coincide.
* ``closest_reference_genome`` is taken from castor's own column rather than
  recomputed by a patristic walk over ``placed.tre``.  castor's value IS the nearest
  known tip by patristic distance (``find_nearest_tips`` on the same tree, same
  known-tip list that produced NSTI), so recomputing it could only introduce
  disagreement.  If the column is absent the run exits 3 rather than guessing.
* ``--mode recipe_control`` subsets the given ``--traits`` table to the reference
  tree's tips before calling ``hsp.py``.  ``castor_hsp.R`` aborts when
  ``num_known / nrow(trait_values) < 0.1``, which the stock 26,868-row ``ec.txt.gz``
  against a ~920-tip control reference always trips.  Subsetting is prediction-neutral
  (``castor_hsp.R`` reindexes the table to the tree's tips anyway, and trait rows that
  are not tips cannot be anyone's nearest tip) -- it only removes the guard.
* ``--ref-dir`` is repeatable so ``--routing nsti`` can be given both domain
  references; with one reference the sibling is derived by swapping ``_bac_``/``_arc_``.
  Under ``--routing nsti`` the per-domain placements live in ``<outdir>/{bac,arc}/``,
  so ``<outdir>/placed.tre`` exists only under ``--routing forced`` (one placement per
  tree is inherent).  ``--routing forced`` -- the primary arm -- matches the contract
  byte for byte.
* ``--queries`` in a fold mode replaces the *sequence source* (so the error-injected
  arm runs the identical pipeline on mutated queries), but the md5 selection still
  comes from the fold plan.  In ``--mode recipe_control`` the file is used verbatim.

Every PICRUSt2 call goes through ``bash -lc 'source ~/opt/picrust2/activate.sh && ...'``
because PICRUSt2 lives in a different interpreter; nothing from ``picrust2`` is ever
imported into py_venv (``/home/freiburger/Documents/py_venv/bin/python``, the only
interpreter for this file).

Resumable: a stage whose output exists is skipped unless ``--force``.
Exit codes: 0 ok, 2 usage/config error, 3 failed assertion (C3/C4).
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)                      # region_validation/shared_reference
RV = os.path.dirname(SR)                        # region_validation
REPO = os.path.dirname(RV)                      # prFBA
sys.path.insert(0, os.path.join(RV, "scripts"))

ACTIVATE = os.path.expanduser("~/opt/picrust2/activate.sh")
DOMAIN_OF = {"Bacteria": "bac", "Archaea": "arc"}
DOMAIN_LONG = {"bac": "Bacteria", "arc": "Archaea"}

# ``frozen`` (SPEC Section 1) verbatim.  Used only as the fallback: the fold plan's
# ``meta.frozen`` (or the reference manifest's) wins, so one authority -- plan_folds --
# fixes the configuration of a run and this component cannot drift from it.
FROZEN_DEFAULT = {
    "tau_amplicon": 0.987,
    "tau_full16s": 0.987,
    "k_folds": 10,
    "seed": 1729,
    "picrust2": {"hsp_method": "mp", "edge_exponent": 0.5, "min_align": 0.80,
                 "placement_tool": "epa-ng", "max_nsti_reported": [None, 2.0],
                 "calculate_nsti": True},
    "density_fractions": [0.25, 0.5, 1.0],
    "label_prefix": "g_",
}

EXIT_USAGE, EXIT_ASSERT = 2, 3


# --------------------------------------------------------------------------- #
# small utilities
# --------------------------------------------------------------------------- #
def log(msg: str) -> None:
    print(f"[run_picrust2_fold] {msg}", flush=True)


def die(code: int, msg: str) -> "None":
    print(f"[run_picrust2_fold] FATAL: {msg}", file=sys.stderr, flush=True)
    sys.exit(code)


def now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


def git_commit() -> str:
    try:
        return subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"],
                              capture_output=True, text=True,
                              timeout=30).stdout.strip()
    except Exception:
        return ""


def load_json(path: str):
    try:
        import orjson
        with open(path, "rb") as fh:
            return orjson.loads(fh.read())
    except ImportError:
        with open(path) as fh:
            return json.load(fh)


def dump_json(path: str, obj) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=False, default=str)
        fh.write("\n")
    os.replace(tmp, path)


def picrust2_call(cmd: str, log_path: str, label: str, check: bool = True) -> int:
    """Run one PICRUSt2 command in the PICRUSt2 environment, teeing to log_path.

    ``check=False`` returns the exit code instead of aborting, so the caller can
    inspect the log (the §5.4 castor crash-column retry needs exactly that).
    """
    full = f"source {ACTIVATE} && {cmd}"
    log(f"{label}: {cmd}")
    t0 = time.time()
    with open(log_path, "w") as fh:
        fh.write(f"# {now_iso()}\n# {cmd}\n\n")
        fh.flush()
        proc = subprocess.run(["bash", "-lc", full], stdout=fh,
                              stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        if not check:
            log(f"{label}: exit {proc.returncode} (see {log_path})")
            return proc.returncode
        tail = ""
        try:
            with open(log_path) as fh:
                tail = "".join(fh.readlines()[-25:])
        except OSError:
            pass
        die(EXIT_USAGE, f"{label} failed (exit {proc.returncode}); see {log_path}\n{tail}")
    log(f"{label}: ok ({time.time() - t0:.1f}s)")
    return 0


# --------------------------------------------------------------------------- #
# FASTA / Stockholm / Newick / trait-table readers (py_venv only, no picrust2)
# --------------------------------------------------------------------------- #
def read_fasta(path: str) -> dict:
    """{first header field: uppercase sequence}; rejects duplicate ids."""
    seqs, name, buf = {}, None, []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line[0] == ">":
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                if name in seqs:
                    die(EXIT_USAGE, f"duplicate FASTA id {name!r} in {path}")
                buf = []
            else:
                buf.append(line.strip().upper())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def fasta_ids(path: str) -> set:
    ids = set()
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def write_fasta(seqs: dict, path: str, order=None) -> None:
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        for name in (order if order is not None else sorted(seqs)):
            fh.write(f">{name}\n{seqs[name]}\n")
    os.replace(tmp, path)


def stockholm_aligned_lengths(path: str, wanted: set) -> dict:
    """{id: len(degapped aligned seq)} for ids in ``wanted``.

    Mirrors ``picrust2.util.read_stockholm(clean_char=True)`` followed by the
    ``check_alignments`` degapping (``replace('-','')`` then ``replace('.','')``),
    which is the exact quantity ``--min_align`` thresholds.
    """
    out = {k: 0 for k in wanted}
    if not wanted or not os.path.exists(path):
        return out
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line or line[0] == "#" or line[:2] == "//":
                continue
            parts = line.split(None, 1)
            if len(parts) != 2 or parts[0] not in out:
                continue
            seq = parts[1].strip().upper().replace(".", "-").replace("-", "")
            out[parts[0]] += len(seq)
    return out


_NEWICK_TOKEN = re.compile(r"""('(?:[^']|'')*'|[^(),:;\[\]]+)""")


def newick_tip_labels(path: str) -> set:
    """Tip labels of a newick tree: a label immediately following '(' or ','.

    Comments ``[...]`` are stripped first; single-quoted labels are unquoted.  This
    is deliberately a tiny local parser -- ete3 is not installed in py_venv, and the
    only thing needed from the tree here is its tip-label set (C4/C3).
    """
    with open(path) as fh:
        text = fh.read()
    text = re.sub(r"\[[^\]]*\]", "", text)
    tips, i, prev = set(), 0, ""
    n = len(text)
    while i < n:
        ch = text[i]
        if ch in "(),:;":
            prev = ch
            i += 1
            continue
        if ch.isspace():
            i += 1
            continue
        m = _NEWICK_TOKEN.match(text, i)
        if not m:
            i += 1
            continue
        label = m.group(1).strip()
        if label.startswith("'") and label.endswith("'") and len(label) >= 2:
            label = label[1:-1].replace("''", "'")
        if prev in ("(", ",") and label:
            tips.add(label)
        i = m.end()
        prev = "LABEL"
    return tips


def topen(path: str):
    return gzip.open(path, "rt", newline="") if path.endswith(".gz") else open(
        path, "rt", newline="")


def trait_header(path: str) -> list:
    with topen(path) as fh:
        return next(csv.reader(fh, delimiter="\t"))


def trait_rows(path: str) -> set:
    rows = set()
    with topen(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd, None)
        for row in rd:
            if row:
                rows.add(row[0])
    return rows


def subset_trait_table(src: str, keep: set, dest: str, strict: bool = False) -> dict:
    """``dest`` = ``src`` restricted to rows in ``keep``, minus degenerate columns.

    Two things happen here and both are forced by downstream tools, not chosen:

    * **Rows.** ``castor_hsp.R`` aborts when ``num_known / nrow(trait_values) < 0.1``,
      which the stock 26,868-row ``ec.txt.gz`` against a ~920-tip control reference
      always trips.  Restricting the table to the tree's tips is prediction-neutral
      (``castor_hsp.R`` reindexes to the tree's tips anyway, and a trait row that is
      not a tip cannot be anyone's nearest tip, so NSTI is unchanged too).
    * **Columns.** Once rows are restricted, a column with only ONE distinct value
      across the survivors gives castor's ``mp`` a single state and crashes it with
      ``invalid 'ncol' value``.  This is SPEC §5.4's mechanical screen; for a fold
      reference ``build_reference.py`` has already applied it (rows == tips, so this
      function is never called), and it is applied here only on the path where the
      table had to be subset -- i.e. ``--mode recipe_control``.  Both recipe arms
      share the same retained row set, so they get an identical screen.

      The test is ``len(set(values)) == 1``, NOT "nonzero prevalence is 0 or n".
      On the binary fold tables the two are identical, but the stock ``ec.txt.gz``
      holds COPY NUMBERS: a column present in all 160 tips with copy numbers
      {1, 2, 3} has nonzero prevalence n yet three states, so castor runs on it
      fine.  Screening on prevalence would throw away exactly the universally
      present core ECs and depress both recipe arms' absolute recall.
      ``strict=True`` additionally drops columns in which some state occurs exactly
      once (the retry of §5.4 step 3, generalised from "prevalence 1 or n-1").
      A dropped column is scored as predicted-0 by the scorer, never as
      unscoreable.
    """
    rows, order = {}, []
    with topen(src) as fin:
        rd = csv.reader(fin, delimiter="\t")
        head = next(rd)
        for row in rd:
            if row and row[0] in keep:
                rows[row[0]] = row
                order.append(row[0])
    n = len(order)
    cols = head[1:]
    states = [{} for _ in cols]
    for rid in order:
        vals = rows[rid][1:]
        for j, v in enumerate(vals):
            c = states[j]
            c[v] = c.get(v, 0) + 1
    keep_j, dropped = [], {}
    for j in range(len(cols)):
        c = states[j]
        if len(c) < 2:
            dropped[cols[j]] = f"constant_{next(iter(c), '')}"
        elif strict and min(c.values()) == 1:
            dropped[cols[j]] = ("singleton_state_"
                                + min(c, key=lambda k: c[k]))
        else:
            keep_j.append(j)
    tmp = dest + ".tmp"
    with gzip.open(tmp, "wt", newline="") as fout:
        wr = csv.writer(fout, delimiter="\t", lineterminator="\n")
        wr.writerow([head[0]] + [cols[j] for j in keep_j])
        for rid in order:
            vals = rows[rid][1:]
            wr.writerow([rid] + [vals[j] for j in keep_j])
    os.replace(tmp, dest)
    return {"n_rows": n, "n_cols_in": len(cols), "n_cols_used": len(keep_j),
            "n_dropped": len(dropped), "dropped_reason": dropped, "strict": strict}


def identify_ref_files(ref_dir: str) -> dict:
    """Mirror ``picrust2.place_seqs.identify_ref_files`` for the epa-ng pipeline."""
    ref_dir = ref_dir.rstrip("/")
    base = os.path.basename(ref_dir)
    fastas = [os.path.join(ref_dir, base + ext)
              for ext in (".fna.gz", ".fasta.gz", ".fna", ".fasta")]
    found = [p for p in fastas if os.path.isfile(p)]
    if len(found) != 1:
        die(EXIT_USAGE, f"expected exactly one of {fastas} in {ref_dir}; found {found}")
    out = {"name": base, "fna": found[0]}
    for key, ext in (("tre", ".tre"), ("hmm", ".hmm"), ("model", ".model")):
        p = os.path.join(ref_dir, base + ext)
        if not os.path.isfile(p):
            die(EXIT_USAGE, f"missing reference file {p}")
        out[key] = p
    return out


# --------------------------------------------------------------------------- #
# configuration resolution
# --------------------------------------------------------------------------- #
def resolve_frozen(fold_plan, ref_dirs) -> dict:
    """``frozen`` from the fold plan, else the first reference manifest, else §1."""
    if fold_plan and isinstance(fold_plan.get("meta"), dict) \
            and fold_plan["meta"].get("frozen"):
        return fold_plan["meta"]["frozen"]
    for d in ref_dirs:
        man = os.path.join(d, "manifest.json")
        if os.path.exists(man):
            m = load_json(man)
            if m.get("frozen"):
                return m["frozen"]
    return FROZEN_DEFAULT


def frozen_picrust2(frozen: dict) -> dict:
    p = dict(FROZEN_DEFAULT["picrust2"])
    p.update(frozen.get("picrust2") or {})
    return p


def domain_from_refname(name: str):
    for dom in ("bac", "arc"):
        if f"_{dom}_" in name or name.endswith(f"_{dom}") or name.startswith(f"{dom}_"):
            return dom
    return None


def sibling_ref_dir(ref_dir: str):
    """f00_bac_d100 -> f00_arc_d100 (and back), in the same parent directory."""
    ref_dir = ref_dir.rstrip("/")
    base = os.path.basename(ref_dir)
    for a, b in (("_bac_", "_arc_"), ("_arc_", "_bac_")):
        if a in base:
            cand = os.path.join(os.path.dirname(ref_dir), base.replace(a, b))
            return cand if os.path.isdir(cand) else None
    return None


# --------------------------------------------------------------------------- #
# query construction
# --------------------------------------------------------------------------- #
def fold_record(fold_plan: dict, fold: int) -> dict:
    for f in fold_plan.get("folds", []):
        if int(f.get("fold", -1)) == fold:
            return f
    die(EXIT_USAGE, f"fold {fold} not present in the fold plan "
                    f"(folds: {[f.get('fold') for f in fold_plan.get('folds', [])]})")


def build_queries(args, frozen):
    """-> (ordered md5 list, {md5: seq}, {md5: [region,...]}, {md5: domain|None}, info)"""
    source = args.queries or args.amplicons
    if not source or not os.path.exists(source):
        die(EXIT_USAGE, f"query sequence source not found: {source!r}")

    if args.mode == "recipe_control":
        seqs = read_fasta(source)
        if not seqs:
            die(EXIT_USAGE, f"no sequences in {source}")
        order = sorted(seqs)
        return (order, seqs, {q: [] for q in order}, {q: None for q in order},
                {"query_source": source, "selection": "verbatim", "n_source": len(seqs),
                 "n_mixed_domain_md5": 0, "md5_missing_from_source": []})

    # ---- fold / include_self / plumbing: md5 selection comes from the fold plan ----
    if not args.fold_plan:
        die(EXIT_USAGE, f"--fold-plan is required for --mode {args.mode}")
    plan = load_json(args.fold_plan)
    rec = fold_record(plan, args.fold)
    test_cells = list(rec.get("test_cells") or [])
    if not test_cells:
        die(EXIT_USAGE, f"fold {args.fold} has no test_cells in {args.fold_plan}")

    truth = load_json(args.truth or os.path.join(RV, "data", "benchmark_truth.json"))
    fanout_path = os.path.join(RV, "data", "amplicon_md5_fanout.json")

    # cell -> md5 via the md5 fanout (md5 -> cells), restricted to this fold's cells.
    if not os.path.exists(fanout_path):
        die(EXIT_USAGE, f"missing {fanout_path}")
    fanout = load_json(fanout_path)
    cellset = set(test_cells)
    md5_cells, mixed = {}, 0
    for md5, cells in fanout.items():
        hit = [c for c in cells if c in cellset]
        if hit:
            md5_cells[md5] = hit
    plan_md5s = set(rec.get("test_md5s") or [])
    if plan_md5s and plan_md5s != set(md5_cells):
        log(f"WARNING: fold_plan.test_md5s ({len(plan_md5s)}) differs from the fanout "
            f"of test_cells ({len(md5_cells)}); using the intersection "
            f"({len(plan_md5s & set(md5_cells))})")
        md5_cells = {m: c for m, c in md5_cells.items() if m in plan_md5s}

    md5_dom, md5_regions = {}, {}
    for md5, cells in md5_cells.items():
        doms, regions = set(), []
        for c in cells:
            t = truth.get(c)
            if t is None:
                continue
            d = DOMAIN_OF.get(t.get("domain"))
            if d:
                doms.add(d)
            regions.append(t.get("region", c.split("__", 1)[-1]))
        if len(doms) > 1:
            mixed += 1
        md5_dom[md5] = sorted(doms)
        md5_regions[md5] = sorted(set(regions))

    wanted = [m for m, doms in md5_dom.items()
              if (args.routing == "nsti") or (args.domain in doms)]
    if not wanted:
        die(EXIT_USAGE, f"no fold-{args.fold} test md5 has truth domain "
                        f"{DOMAIN_LONG.get(args.domain)}")

    src = read_fasta(source)
    missing = sorted(m for m in wanted if m not in src)
    if missing:
        die(EXIT_USAGE, f"{len(missing)} selected md5s absent from {source} "
                        f"(first: {missing[:3]})")
    order = sorted(wanted)
    seqs = {m: src[m] for m in order}
    by_dom = {}
    for m in wanted:
        for d in (md5_dom[m] or ["unknown"]):
            by_dom[d] = by_dom.get(d, 0) + 1
    info = {"query_source": source, "n_query_by_truth_domain": by_dom,
            "selection": f"fold_plan fold {args.fold} test_cells -> md5 "
                         f"({'all domains' if args.routing == 'nsti' else args.domain})",
            "n_source": len(src), "n_mixed_domain_md5": mixed,
            "md5_missing_from_source": []}
    return order, seqs, md5_regions, {m: (md5_dom[m][0] if md5_dom[m] else None)
                                      for m in order}, info


# --------------------------------------------------------------------------- #
# hsp.py output handling
# --------------------------------------------------------------------------- #
META_COLS = ("metadata_NSTI", "closest_reference_genome")


def hsp_scan(path: str):
    """One pass: -> (ec column names, {md5: (nsti, closest_ref)}, row order)."""
    with topen(path) as fh:
        rd = csv.reader(fh, delimiter="\t")
        head = next(rd)
        if head[0] != "sequence":
            die(EXIT_ASSERT, f"{path}: first column header is {head[0]!r}, expected "
                             f"'sequence' (hsp.py writes index_label='sequence')")
        idx = {c: i for i, c in enumerate(head)}
        for c in META_COLS:
            if c not in idx:
                die(EXIT_ASSERT,
                    f"{path}: hsp.py output is missing the {c!r} column. "
                    f"'metadata_NSTI' requires `hsp.py -n`; "
                    f"'closest_reference_genome' is emitted by castor_nsti.R in "
                    f"PICRUSt2 2.6.3 -- its absence means a different PICRUSt2 "
                    f"version is on PATH and the nearest-tip label cannot be "
                    f"recovered without recomputing it.")
        ec_cols = [c for c in head[1:] if c not in META_COLS]
        meta, order = {}, []
        ni, ci = idx["metadata_NSTI"], idx["closest_reference_genome"]
        for row in rd:
            if not row:
                continue
            md5 = row[0]
            try:
                nsti = float(row[ni])
            except (ValueError, IndexError):
                nsti = float("nan")
            meta[md5] = (nsti, row[ci])
            order.append(md5)
    return ec_cols, meta, order


def _positive(value: str) -> bool:
    try:
        return float(value) > 0
    except ValueError:
        return False


def write_combined(dom_tables, winner, union_cols, out_ec_gz):
    """Stream the winning row of each md5 into the combined EC table.

    ``dom_tables`` = [(domain, hsp_path)] in a deterministic order.  ``winner`` maps
    md5 -> domain.  Columns are the union over domains (the §5.4 crash-column screen
    can drop different ECs in the two domain references), 0-filled where a domain
    lacks a column -- which is the honest value: a dropped column is scored as
    predicted-0, never as unscoreable (SPEC §5.4).
    """
    col_index = {c: i for i, c in enumerate(union_cols)}
    pred, n_rows = {}, 0
    tmp = out_ec_gz + ".tmp"
    with gzip.open(tmp, "wt", newline="") as fout:
        wr = csv.writer(fout, delimiter="\t", lineterminator="\n")
        wr.writerow(["sequence"] + union_cols)
        for dom, path in dom_tables:
            with topen(path) as fh:
                rd = csv.reader(fh, delimiter="\t")
                head = next(rd)
                keep = []
                for i, c in enumerate(head):
                    if i == 0 or c in META_COLS:
                        continue
                    c = c[3:] if c.startswith("EC:") else c
                    keep.append((i, col_index[c]))
                if len({d for _, d in keep}) != len(keep):
                    die(EXIT_ASSERT, f"{path}: two EC columns collapse to the same "
                                     f"name after stripping an 'EC:' prefix")
                for row in rd:
                    if not row or winner.get(row[0]) != dom:
                        continue
                    out = ["0"] * len(union_cols)
                    ecs = []
                    for src_i, dst_i in keep:
                        v = row[src_i]
                        out[dst_i] = v
                        if _positive(v):
                            ecs.append(union_cols[dst_i])
                    wr.writerow([row[0]] + out)
                    pred[row[0]] = sorted(ecs)
                    n_rows += 1
    os.replace(tmp, out_ec_gz)
    return pred, n_rows


def write_marker(path: str, rows) -> None:
    """``sequence, metadata_NSTI, best_domain, closest_reference_genome`` -- exactly.

    ``compare.load_picrust2`` reads ``r['best_domain']`` with a plain ``[]``, so the
    column must exist; only ``picrust2_pipeline.py``'s two-reference wrapper would
    otherwise produce it.
    """
    tmp = path + ".tmp"
    with gzip.open(tmp, "wt", newline="") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n")
        wr.writerow(["sequence", "metadata_NSTI", "best_domain",
                     "closest_reference_genome"])
        for md5, nsti, dom, closest in rows:
            wr.writerow([md5, nsti, dom, closest])
    os.replace(tmp, path)


# --------------------------------------------------------------------------- #
# controls C3 / C4
# --------------------------------------------------------------------------- #
def check_c4(ref: dict, traits_path: str, strict: bool, controls: dict) -> dict:
    """C4 -- artefact set equality.  One differing id turns a reference genome into a
    study sequence (or deletes it from the HSP known-set), so this is checked again at
    run time even though build_reference.py checks it at build time."""
    fna = fasta_ids(ref["fna"])
    tips = newick_tip_labels(ref["tre"])
    rows = trait_rows(traits_path)
    measured = {"ref_name": ref["name"], "n_fna": len(fna), "n_tips": len(tips),
                "n_trait_rows": len(rows), "strict": strict,
                "fna_not_in_tree": sorted(fna - tips)[:10],
                "tree_not_in_fna": sorted(tips - fna)[:10],
                "tips_not_in_traits": sorted(tips - rows)[:10],
                "traits_not_in_tips": sorted(rows - tips)[:10]}
    ok = (fna == tips) and (tips <= rows) and (strict is False or rows == tips)
    controls.setdefault("C4", {})[ref["name"]] = {
        "status": "pass" if ok else "fail", "measured": measured,
        "gate": ("set(fna) == set(tre tips) == set(trait rows)" if strict
                 else "set(fna) == set(tre tips) <= set(trait rows)"),
        "falsifies": "a reference genome silently became a study sequence, or was "
                     "deleted from the HSP known-set"}
    if not ok:
        die(EXIT_ASSERT, f"C4 failed for {ref['name']}: {json.dumps(measured)}")
    return measured, tips, rows


def _r_numeric_like(label: str) -> bool:
    """Would R's ``read.table`` parse this single field as a number?

    R types a column numeric only if EVERY field in it parses, so this is the
    per-field half of the test in :func:`check_c3_preflight`.  Bare BV-BRC genome
    ids ("1144682.3") parse; the ``g_``-prefixed labels do not.
    """
    try:
        float(label.strip().strip("'\""))
    except ValueError:
        return False
    return True


def check_c3_preflight(dom: str, tips: set, controls: dict) -> dict:
    """C3, checked BEFORE place_seqs/hsp instead of only after (see below).

    ``castor_nsti.R`` reads the known-tip list with ``read.table(...)$V1``.  R types
    that column numeric only when *every* label parses as a number -- which is exactly
    what happens to an unprefixed BV-BRC panel, where "168471.10" comes back as
    168471.1 and the tip vanishes from the known-set.

    :func:`check_c3` catches the consequence, but only after ``hsp.py`` has run.  In
    practice the coerced known-set also perturbs per-column state counts, so castor
    dies first inside ``hsp_max_parsimony`` with ``invalid 'ncol' value`` -- an error
    whose documented cause is a degenerate *trait column*, pointing an operator at
    build_reference.py's §5.4 screen instead of at the labels.  Measured on the real
    143-tip fold-0 reference with the ``g_`` prefix stripped: exit 2, two opaque R
    tracebacks, no mention of labels.  Hence this preflight.
    """
    coercible = sorted(t for t in tips if _r_numeric_like(t))
    all_numeric = bool(tips) and len(coercible) == len(tips)
    at_risk = sorted(t for t in coercible
                     if t.strip().strip("'\"").split(".")[-1].endswith("0"))
    measured = {"domain": dom, "n_tips": len(tips),
                "n_r_numeric_like": len(coercible),
                "all_labels_numeric_like": all_numeric,
                "n_would_be_coerced": len(at_risk) if all_numeric else 0,
                "would_be_coerced_sample": at_risk[:10] if all_numeric else []}
    controls.setdefault("C3_preflight", {})[dom] = {
        "status": "fail" if all_numeric else "pass", "measured": measured,
        "gate": "NOT every reference tip label parses as an R number",
        "falsifies": "castor_nsti.R's read.table would type the known-tip column "
                     "numeric and silently drop every id with a trailing zero; this "
                     "is the bug the 'g_' label prefix exists to prevent"}
    if all_numeric:
        die(EXIT_ASSERT,
            f"C3 preflight failed for domain {dom}: all {len(tips)} reference tip "
            f"labels parse as R numbers, so castor_nsti.R's read.table will type the "
            f"known-tip column numeric and drop the {len(at_risk)} id(s) ending in a "
            f"trailing zero (e.g. {at_risk[:3]}). Rebuild the reference with a "
            f"non-numeric label prefix (frozen.label_prefix, default 'g_'). "
            f"NOTE: without this preflight the symptom is castor's "
            f"\"invalid 'ncol' value\" crash, which misattributes the fault to a "
            f"degenerate trait column.")
    return measured


def check_c3(dom: str, hsp_rows: set, placed: set, n_query: int, controls: dict) -> dict:
    """C3 -- tip-label round trip (set equality, not count equality).

    ``castor_nsti.R`` reads the known-tip list with ``read.table(...)$V1``; a label
    set that R types as numeric loses every id with a trailing zero, which turns
    reference tips into extra unknown-tip rows.  Comparing SETS catches the
    one-tip-in / one-query-out case that count equality would mask.
    """
    measured = {"domain": dom, "n_unknown_rows": len(hsp_rows), "n_placed": len(placed),
                "n_query": n_query,
                "rows_not_placed": sorted(hsp_rows - placed)[:10],
                "placed_not_in_rows": sorted(placed - hsp_rows)[:10]}
    ok = hsp_rows == placed
    controls.setdefault("C3", {})[dom] = {
        "status": "pass" if ok else "fail", "measured": measured,
        "gate": "set(hsp unknown-tip rows) == set(study seqs surviving --min_align)",
        "falsifies": "castor_nsti.R coerced the known-tip labels to numeric and "
                     "deleted reference genomes from the NSTI known-set (this is why "
                     "every reference label carries the 'g_' prefix)"}
    if not ok:
        die(EXIT_ASSERT, f"C3 failed for domain {dom}: {json.dumps(measured)}")
    return measured


# --------------------------------------------------------------------------- #
# one domain: place_seqs -> hsp
# --------------------------------------------------------------------------- #
def run_domain(dom, ref_dir, workdir, queries_fna, query_seqs, p2, args,
               controls, strict_c4):
    os.makedirs(workdir, exist_ok=True)
    ref = identify_ref_files(ref_dir)

    traits_in = args.traits or os.path.join(ref_dir.rstrip("/"), "traits_EC.txt.gz")
    if not os.path.exists(traits_in):
        die(EXIT_USAGE, f"trait table not found: {traits_in} (pass --traits)")

    c4, tips, rows = check_c4(ref, traits_in, strict_c4, controls)
    check_c3_preflight(dom, tips, controls)

    # --- stage 1: placement (resumable) ------------------------------------- #
    placed_tre = os.path.join(workdir, "placed.tre")
    inter = os.path.join(workdir, "intermediate")
    ps_log = os.path.join(workdir, "place_seqs.log")
    study_aln = os.path.join(inter, "study_seqs_hmmalign.fasta")
    stamp_path = os.path.join(workdir, "stage_place.json")
    # The stamp hashes the reference ARTEFACTS, not just the directory path: the
    # orchestrator rebuilds refs/<NAME>/ in place, and a path-only stamp would make
    # this runner silently reuse a placement computed against the PREVIOUS tree.
    stamp = {"ref": os.path.abspath(ref_dir), "n_query": len(query_seqs),
             "min_align": p2["min_align"], "placement_tool": p2["placement_tool"],
             "queries_sha256": sha256_file(queries_fna),
             "ref_sha256": {k: sha256_file(ref[k]) for k in ("fna", "tre", "model")},
             "traits_sha256": sha256_file(traits_in)}
    prev = load_json(stamp_path) if os.path.exists(stamp_path) else None
    reuse = (not args.force and prev == stamp
             and os.path.exists(placed_tre) and os.path.getsize(placed_tre) > 0
             and os.path.exists(study_aln))
    if reuse:
        log(f"[{dom}] reusing existing placement in {workdir}")
    else:
        shutil.rmtree(inter, ignore_errors=True)
        for p in (placed_tre, os.path.join(workdir, "EC_predicted.tsv.gz")):
            if os.path.exists(p):
                os.remove(p)
        picrust2_call(
            f"place_seqs.py -s {queries_fna} -r {ref_dir.rstrip('/')} "
            f"-o {placed_tre} -p {args.processes} "
            f"-t {p2['placement_tool']} --min_align {p2['min_align']} "
            f"--chunk_size {args.chunk_size} --intermediate {inter} --verbose",
            ps_log, f"[{dom}] place_seqs")
        dump_json(stamp_path, stamp)

    if not os.path.exists(study_aln):
        die(EXIT_USAGE, f"[{dom}] place_seqs produced no {study_aln}")
    placed = fasta_ids(study_aln)
    dropped = set(query_seqs) - placed
    extra = placed - set(query_seqs)
    if extra:
        die(EXIT_ASSERT, f"[{dom}] {len(extra)} aligned study sequences are not "
                         f"queries (first: {sorted(extra)[:3]})")
    aligned_len = stockholm_aligned_lengths(
        os.path.join(inter, "query_align.stockholm"), dropped)
    log(f"[{dom}] placed {len(placed)}/{len(query_seqs)}; "
        f"{len(dropped)} dropped by --min_align {p2['min_align']}")

    # --- stage 2: trait table for hsp (recipe_control needs the subset) ------ #
    # tips / rows come from check_c4 above -- re-reading a 26k x 3.7k gzip table
    # twice per domain is the one measurable cost in this script.
    screen = {"applied": False, "n_rows": len(rows), "n_cols_in": None,
              "n_cols_used": None, "n_dropped": 0, "strict": False}
    if rows == tips:
        # A fold reference: build_reference.py already applied the §5.4 screen and
        # recorded dropped_ec in its manifest. Touch nothing.
        traits_used, n_used = traits_in, len(rows)
    else:
        traits_used = os.path.join(workdir, "traits_used.tsv.gz")
        screen_path = traits_used + ".screen.json"
        if (not args.force) and os.path.exists(traits_used) \
                and os.path.exists(screen_path):
            screen = load_json(screen_path)
        else:
            screen = subset_trait_table(traits_in, tips, traits_used)
            dump_json(screen_path, screen)
        screen["applied"] = True
        n_used = screen["n_rows"]
        log(f"[{dom}] trait table {len(rows)} -> {n_used} rows (tree tips), "
            f"{screen['n_cols_in']} -> {screen['n_cols_used']} columns "
            f"({screen['n_dropped']} constant across the survivors)")

    # --- stage 3: hidden-state prediction (resumable, §5.4 retry) ----------- #
    hsp_out = os.path.join(workdir, "EC_predicted.tsv.gz")
    hsp_log = os.path.join(workdir, "hsp.log")
    hsp_stamp_path = os.path.join(workdir, "stage_hsp.json")
    hsp_stamp = {"placed_tre_sha256": sha256_file(placed_tre),
                 "traits_used_sha256": sha256_file(traits_used),
                 "hsp_method": p2["hsp_method"],
                 "edge_exponent": p2["edge_exponent"], "seed": args.seed,
                 "calculate_nsti": bool(p2.get("calculate_nsti", True))}
    prev_hsp = load_json(hsp_stamp_path) if os.path.exists(hsp_stamp_path) else None
    hsp_fresh = (not args.force and prev_hsp == hsp_stamp
                 and os.path.exists(hsp_out) and os.path.getsize(hsp_out) > 0)
    if not hsp_fresh:
        if os.path.exists(hsp_out):
            os.remove(hsp_out)
        nsti_flag = "-n " if p2.get("calculate_nsti", True) else ""
        # hsp.py --chunk_size is COLUMNS per job, so the 500 default caps the job
        # count at ceil(n_cols/500) -- 5 jobs for a 2,236-EC table, however many
        # processes are asked for.  Chunks are independent columns, so sizing them
        # to the process count is result-identical and actually uses the cores.
        n_cols = max(1, len(trait_header(traits_used)) - 1)
        hsp_chunk = min(500, max(25, -(-n_cols // max(1, args.processes))))

        def _hsp(tt, check):
            return picrust2_call(
                f"hsp.py -t {placed_tre} --observed_trait_table {tt} "
                f"-o {hsp_out} -m {p2['hsp_method']} -e {p2['edge_exponent']} "
                f"{nsti_flag}-p {args.processes} --chunk_size {hsp_chunk} "
                f"--seed {args.seed}",
                hsp_log, f"[{dom}] hsp", check=check)

        if not screen["applied"]:
            _hsp(traits_used, True)
        elif _hsp(traits_used, False) != 0:
            # SPEC §5.4 step 3: on castor's `invalid 'ncol' value`, additionally drop
            # prevalence-1 and prevalence-(n-1) columns and retry ONCE; a second
            # failure aborts.
            with open(hsp_log) as fh:
                logtext = fh.read()
            if "invalid 'ncol' value" not in logtext:
                die(EXIT_USAGE, f"[{dom}] hsp.py failed for a reason other than "
                                f"castor's degenerate-column crash; see {hsp_log}")
            log(f"[{dom}] castor 'invalid ncol' -- retrying with the strict "
                f"column screen (§5.4 step 3)")
            screen = subset_trait_table(traits_in, tips, traits_used, strict=True)
            dump_json(traits_used + ".screen.json", screen)
            screen["applied"] = True
            screen["retry_used"] = True
            n_used = screen["n_rows"]
            hsp_chunk = min(500, max(25, -(-max(1, len(
                trait_header(traits_used)) - 1) // max(1, args.processes))))
            _hsp(traits_used, True)
            hsp_stamp["traits_used_sha256"] = sha256_file(traits_used)
        dump_json(hsp_stamp_path, hsp_stamp)
    else:
        log(f"[{dom}] reusing existing {hsp_out}")

    ec_cols, meta, order = hsp_scan(hsp_out)
    check_c3(dom, set(order), placed, len(query_seqs), controls)
    if set(ec_cols) != (set(trait_header(traits_used)[1:])):
        die(EXIT_ASSERT, f"[{dom}] hsp.py EC columns ({len(ec_cols)}) differ from the "
                         f"trait table's columns "
                         f"({len(trait_header(traits_used)) - 1})")
    return {"domain": dom, "ref": ref, "ref_dir": os.path.abspath(ref_dir),
            "hsp": hsp_out, "traits_in": traits_in, "traits_used": traits_used,
            "n_trait_rows_used": n_used, "trait_screen": screen,
            "ec_cols": ec_cols, "meta": meta,
            "order": order, "placed": placed, "dropped": dropped,
            "aligned_len": aligned_len, "c4": c4, "placed_tre": placed_tre}


def compact_screen(screen: dict) -> dict:
    """predictions.json records the screen's shape plus a sample, not 3k EC names."""
    out = {k: v for k, v in screen.items() if k != "dropped_reason"}
    dropped = sorted(screen.get("dropped_reason") or {})
    out["dropped_ec_sample"] = dropped[:20]
    out["dropped_ec_full"] = "see <workdir>/traits_used.tsv.gz.screen.json"
    return out


def sha256_file(path: str) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_versions(path: str) -> str:
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path) as fh:
            head = fh.read()
    else:
        cmd = ("python -c \"from importlib.metadata import version; "
               "print('picrust2', version('PICRUSt2'))\" ; "
               "epa-ng --version 2>&1 | head -1 ; "
               "echo -n 'gappa ' ; gappa --version 2>&1 | head -1 ; "
               "hmmalign -h 2>&1 | sed -n 2p ; "
               "Rscript -e 'cat(\"castor\", as.character(packageVersion(\"castor\")), "
               "\"\\n\")' 2>&1 | tail -1")
        out = subprocess.run(["bash", "-lc", f"source {ACTIVATE} && {cmd}"],
                             capture_output=True, text=True)
        head = out.stdout + out.stderr
        with open(path, "w") as fh:
            fh.write(head)
    m = re.search(r"picrust2\s+(\S+)", head)
    return m.group(1) if m else "unknown"


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True,
                    choices=["fold", "recipe_control", "include_self", "plumbing"])
    ap.add_argument("--ref-dir", required=True, action="append", metavar="DIR",
                    help="custom reference directory NAME/ holding NAME.{fna,tre,hmm,"
                         "model} (+ traits_EC.txt.gz). Repeatable: give both domain "
                         "references for --routing nsti.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--fold-plan", default=None,
                    help="required for --mode fold|include_self|plumbing")
    ap.add_argument("--panel", default=None,
                    help="panel_eligible.json; used to check that every "
                         "closest_reference_genome is a known panel label")
    ap.add_argument("--fold", type=int, default=0)
    ap.add_argument("--domain", default="bac", choices=["bac", "arc"])
    ap.add_argument("--density", type=float, default=1.0)
    ap.add_argument("--amplicons",
                    default=os.path.join(RV, "data", "combined_amplicons.unique.fasta"))
    ap.add_argument("--queries", default=None,
                    help="explicit query FASTA. In a fold mode it replaces the "
                         "SEQUENCE SOURCE (error-injected arm) while the md5 "
                         "selection still comes from the fold plan; in "
                         "--mode recipe_control it is used verbatim and is REQUIRED.")
    ap.add_argument("--truth", default=None,
                    help='{"<query_id>": ["<ec>", ...]} for --mode recipe_control '
                         "(REQUIRED there); otherwise benchmark_truth.json, which "
                         "supplies the forced-routing domain of each test cell.")
    ap.add_argument("--traits", default=None,
                    help="override <ref-dir>/traits_EC.txt.gz")
    ap.add_argument("--processes", type=int, default=30)
    ap.add_argument("--chunk-size", type=int, default=5000,
                    help="place_seqs --chunk_size (query seqs per EPA-ng chunk)")
    ap.add_argument("--routing", default="forced", choices=["forced", "nsti"])
    ap.add_argument("--force", action="store_true",
                    help="recompute every stage instead of reusing existing output")
    args = ap.parse_args()

    t0 = time.time()
    if not os.path.exists(ACTIVATE):
        die(EXIT_USAGE, f"PICRUSt2 activation script not found: {ACTIVATE}")
    if args.mode == "recipe_control":
        if not args.queries or not args.truth:
            die(EXIT_USAGE, "--mode recipe_control requires both --queries and --truth")
        if args.routing != "forced":
            die(EXIT_USAGE, "--mode recipe_control supports --routing forced only")
    elif not args.fold_plan:
        die(EXIT_USAGE, f"--fold-plan is required for --mode {args.mode}")
    if args.routing == "nsti" and abs(args.density - 1.0) > 1e-9:
        die(EXIT_USAGE, "--routing nsti is restricted to --density 1.0 (SPEC D9): it "
                        "doubles the PICRUSt2 run count and is a secondary arm only")

    # ---- reference directories -------------------------------------------- #
    ref_dirs = [d.rstrip("/") for d in args.ref_dir]
    for d in ref_dirs:
        if not os.path.isdir(d):
            die(EXIT_USAGE, f"--ref-dir {d} is not a directory")
    if args.routing == "nsti" and len(ref_dirs) == 1:
        sib = sibling_ref_dir(ref_dirs[0])
        if not sib:
            die(EXIT_USAGE, "--routing nsti needs both domain references: pass "
                            "--ref-dir twice (the _bac_/_arc_ sibling of "
                            f"{ref_dirs[0]} was not found)")
        ref_dirs.append(sib)
        log(f"derived sibling reference {sib}")
    if args.routing == "forced" and len(ref_dirs) != 1:
        die(EXIT_USAGE, "--routing forced takes exactly one --ref-dir")

    dom_of_ref = {}
    for d in ref_dirs:
        dom = domain_from_refname(os.path.basename(d))
        if args.routing == "forced":
            if dom and dom != args.domain:
                log(f"WARNING: --domain {args.domain} but reference name "
                    f"{os.path.basename(d)} looks like domain {dom}")
            dom = args.domain
        if dom is None:
            die(EXIT_USAGE, f"cannot infer bac/arc from reference name "
                            f"{os.path.basename(d)}; needed for --routing nsti")
        if dom in dom_of_ref:
            die(EXIT_USAGE, f"two references map to domain {dom}: "
                            f"{dom_of_ref[dom]} and {d}")
        dom_of_ref[dom] = d

    plan = load_json(args.fold_plan) if args.fold_plan else None
    frozen = resolve_frozen(plan, ref_dirs)
    p2 = frozen_picrust2(frozen)
    args.seed = int(frozen.get("seed", FROZEN_DEFAULT["seed"]))
    label_prefix = frozen.get("label_prefix", "g_")

    os.makedirs(args.outdir, exist_ok=True)
    controls = {}

    # ---- queries ----------------------------------------------------------- #
    order, seqs, md5_regions, _md5_dom, qinfo = build_queries(args, frozen)
    queries_fna = os.path.join(args.outdir, "study_seqs.fna")
    write_fasta(seqs, queries_fna, order=order)
    log(f"{len(order)} query sequences -> {queries_fna} "
        f"(source {qinfo['query_source']})")

    if args.mode == "recipe_control":
        # --truth is the query-id registry here (there is no benchmark_truth cell for
        # a stock reference genome); the scoring itself is score_folds.py's job.
        truth_ids = set(load_json(args.truth))
        missing = [q for q in order if q not in truth_ids]
        if len(missing) == len(order):
            die(EXIT_USAGE, f"none of the {len(order)} queries appear in {args.truth}; "
                            f"the query/truth pair does not match")
        if missing:
            log(f"WARNING: {len(missing)} queries have no entry in {args.truth} "
                f"(first: {missing[:3]}); the scorer will drop them")
        qinfo["n_truth_queries"] = len(truth_ids)
        qinfo["n_queries_without_truth"] = len(missing)

    picrust2_version = write_versions(os.path.join(args.outdir, "versions.txt"))

    # ---- per-domain placement + HSP ---------------------------------------- #
    strict_c4 = args.mode != "recipe_control"
    results = {}
    for dom in sorted(dom_of_ref):
        workdir = (args.outdir if args.routing == "forced"
                   else os.path.join(args.outdir, dom))
        results[dom] = run_domain(dom, dom_of_ref[dom], workdir, queries_fna, seqs,
                                  p2, args, controls, strict_c4)

    # ---- routing: pick the winning domain per md5 -------------------------- #
    # Under --routing forced there is one domain and every placed row wins trivially.
    # Under --routing nsti the lower-NSTI domain wins (D9); NaN never beats a number,
    # and exact ties go to the alphabetically first domain, so the result is
    # deterministic.
    winner, nsti_of, closest_of = {}, {}, {}
    for dom in sorted(results):
        for md5, (nsti, closest) in results[dom]["meta"].items():
            if md5 in winner:
                old = nsti_of[md5]
                is_better = (nsti == nsti) and (old != old or nsti < old)
            else:
                is_better = True
            if is_better:
                winner[md5], nsti_of[md5], closest_of[md5] = dom, nsti, closest

    # ---- combined outputs the scorer reads --------------------------------- #
    union_cols = []
    seen = set()
    for dom in sorted(results):
        for c in results[dom]["ec_cols"]:
            c = c[3:] if c.startswith("EC:") else c
            if c not in seen:
                seen.add(c)
                union_cols.append(c)
    ec_gz = os.path.join(args.outdir, "combined_EC_predicted.tsv.gz")
    pred_ec, n_rows = write_combined(
        [(dom, results[dom]["hsp"]) for dom in sorted(results)], winner, union_cols,
        ec_gz)
    write_marker(os.path.join(args.outdir, "combined_marker_predicted_and_nsti.tsv.gz"),
                 [(md5, nsti_of[md5], winner[md5], closest_of[md5])
                  for md5 in sorted(winner)])
    if n_rows != len(winner):
        die(EXIT_ASSERT, f"combined EC table has {n_rows} rows for {len(winner)} "
                         f"routed queries")

    # ---- placement loss (C7: reported, never gated) ------------------------ #
    dropped_all = set(order)
    for dom in results:
        dropped_all &= results[dom]["dropped"]
    drop_rec = {}
    for md5 in sorted(dropped_all):
        alen = max(results[d]["aligned_len"].get(md5, 0) for d in results)
        drop_rec[md5] = {"len": len(seqs[md5]), "aligned_len": int(alen),
                         "regions": md5_regions.get(md5, [])}
    dump_json(os.path.join(args.outdir, "min_align_dropped.json"), drop_rec)
    by_region = {}
    for rec in drop_rec.values():
        for r in rec["regions"] or ["unknown"]:
            by_region[r] = by_region.get(r, 0) + 1
    controls["C7"] = {
        "status": "reported", "gate": "reported, not gated",
        "measured": {"n_query": len(order), "n_min_align_dropped": len(drop_rec),
                     "min_align": p2["min_align"], "by_region": by_region},
        "falsifies": "nothing -- scoring these as PICRUSt2 abstentions would charge "
                     "PICRUSt2 for a data limitation shared with prFBA; the scorer "
                     "must drop them from BOTH methods' denominators"}

    # ---- label hygiene (the reason C3 exists) ------------------------------ #
    label_check = {"label_prefix": label_prefix}
    ref_labels = set()
    for dom in results:
        ref_labels |= fasta_ids(results[dom]["ref"]["fna"])
    bad = sorted({c for c in closest_of.values()
                  if c.strip("'\"") not in ref_labels})
    label_check["closest_ref_not_a_reference_tip"] = bad[:10]
    label_check["n_closest_ref_not_a_reference_tip"] = len(bad)
    if bad:
        log(f"WARNING: {len(bad)} closest_reference_genome labels are not reference "
            f"tips (first: {bad[:3]}); C3 passed, so this is a label-quoting "
            f"difference rather than known-set corruption")
    if label_prefix:
        n_unprefixed = sum(1 for c in ref_labels if not c.startswith(label_prefix))
        label_check["n_reference_labels_without_prefix"] = n_unprefixed
    if args.panel and os.path.exists(args.panel):
        panel = load_json(args.panel)
        plabels = {g.get("label", f"{label_prefix}{gid}")
                   for gid, g in (panel.get("genomes") or {}).items()}
        label_check["n_closest_ref_outside_panel"] = len(
            {c for c in closest_of.values() if c not in plabels})

    # ---- predictions.json -------------------------------------------------- #
    n_unknown = len(winner)
    pred = {md5: {"ec": pred_ec.get(md5, []),
                  "nsti": (None if nsti_of[md5] != nsti_of[md5] else nsti_of[md5]),
                  "best_domain": winner[md5],
                  "closest_reference_genome": closest_of[md5]}
            for md5 in sorted(winner)}
    ref_name = (results[args.domain]["ref"]["name"] if args.routing == "forced"
                else [results[d]["ref"]["name"] for d in sorted(results)])
    meta = {
        "mode": args.mode,
        "fold": (args.fold if args.mode != "recipe_control" else None),
        "domain": (args.domain if args.routing == "forced" else "both"),
        "density": args.density,
        "ref_name": ref_name,
        "routing": args.routing,
        "frozen": frozen,
        "n_query": len(order),
        "n_placed": len(winner),
        "n_min_align_dropped": len(drop_rec),
        "min_align_dropped_md5": sorted(drop_rec),
        "n_unknown_rows": n_unknown,
        # C3 gate: unknown-tip rows must match the queries that SURVIVED --min_align.
        # The literal `== n_query` is unsatisfiable once --min_align drops anything,
        # and C7 removes those cells from both methods' denominators.
        "unknown_rows_match_queries": n_unknown == len(order) - len(drop_rec),
        "n_unknown_rows_equals_n_query": n_unknown == len(order),
        "n_unknown_rows_per_domain": {d: len(results[d]["order"])
                                      for d in sorted(results)},
        "n_reference_tips": {d: results[d]["c4"]["n_tips"] for d in sorted(results)},
        "n_trait_columns": len(union_cols),
        "trait_table": {d: {"input": results[d]["traits_in"],
                            "used": results[d]["traits_used"],
                            "n_rows_used": results[d]["n_trait_rows_used"],
                            "screen": compact_screen(results[d]["trait_screen"])}
                        for d in sorted(results)},
        "query_info": qinfo,
        "label_check": label_check,
        "controls": controls,
        "picrust2_version": picrust2_version,
        "wall_seconds": round(time.time() - t0, 1),
        "generated": now_iso(),
        "git_commit": git_commit(),
        "command": " ".join(sys.argv),
    }
    if not meta["unknown_rows_match_queries"]:
        die(EXIT_ASSERT, "C3 aggregate failed: n_unknown_rows "
                         f"{n_unknown} != n_query - n_min_align_dropped "
                         f"{len(order) - len(drop_rec)}")
    dump_json(os.path.join(args.outdir, "predictions.json"),
              {"meta": meta, "pred": pred})
    log(f"wrote {args.outdir}/predictions.json: {len(pred)} predictions, "
        f"{len(union_cols)} EC columns, {len(drop_rec)} min_align drops "
        f"({meta['wall_seconds']}s)")


if __name__ == "__main__":
    main()
