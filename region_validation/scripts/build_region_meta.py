#!/usr/bin/env python
"""S11 — region metadata + amplification tables (DESIGN §8 "Precomputed tables").

Emits the two flat tables the renderers (and the cost/survivorship audits) read:

  * ``PATHS.region_meta_csv`` — one row per region with the amplicon-length
    summary, amplicon count, overall and per-domain extract rate, the primer
    pair / runtime reverse-complement / insert bounds, the region's panel
    ``priority`` + ``domain`` scope, and the reproducibility version stamps
    (python / numpy / cupy / CUDA-driver) demanded by DESIGN §6 (verifier A3) so
    every figure caption can cite the exact toolchain.
  * ``PATHS.amplification_by_phylum_csv`` — long ``(region, domain, phylum)`` rows
    of ``extract_rate = n_amplified / n_source`` (the F13 selection-bias / "which
    taxa fail to amplify" answer), with explicit numerator/denominator columns and
    a Bacteria/Archaea ``domain`` split.

The region-level counts, length summaries, and the per-domain extract rates all
come from ``PATHS.combined_stats`` (written by ``insilico_pcr.py --combined-fasta``):
that stats file already carries a ``regions[r].by_domain[dom]`` block with the exact
per-domain ``extract_rate`` / ``n_amplicons`` (the same Bacteria/Archaea split
DESIGN §4 asks for), so the ``region_meta`` domain columns are read straight off it
and stay byte-consistent with the upstream run.  The finer ``(region, phylum)``
stratification (which the stats do NOT carry) is the only thing computed here by a
join: the amplified ``(genome, region)`` set — enumerated from ``PATHS.amplicon_fanout``
(S4), or the combined FASTA itself as a fallback — against the per-source-genome
domain/phylum labels read from the S1 selection manifest (``PATHS.selection_json``,
shape ``{"meta":..., "genomes":[record,...]}``), with ``PATHS.truth_json`` (S7) as a
backfill.  If the upstream stats predate the ``by_domain`` split, the per-domain
region_meta columns gracefully fall back to the same join.

Resumable: if both output CSVs already exist (and ``--force`` is not given) the
stage short-circuits.  Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.

I/O contract
------------
reads : PATHS.combined_stats     (data/combined_amplicons.stats.json)   [required]
        PATHS.selection_json     (data/selection_10k.json)              [optional]
        PATHS.truth_json         (data/benchmark_truth.json)            [optional]
        PATHS.amplicon_fanout    (data/amplicon_md5_fanout.json)        [optional]
        PATHS.combined_fasta     (data/combined_amplicons.fasta)        [fallback]
writes: PATHS.region_meta_csv             (data/region_meta.csv)
        PATHS.amplification_by_phylum_csv (data/amplification_by_phylum.csv)
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import sys

# Import the shared contract (paths/constants/helpers) — never hardcode a path.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C  # noqa: E402


# --------------------------------------------------------------------------- #
# Version stamps (DESIGN §6 / §8 — verifier A3 reproducibility pins)
# --------------------------------------------------------------------------- #
def version_stamps() -> dict:
    """Best-effort python / numpy / cupy / CUDA-driver version strings.

    This stage is a pure CPU post-step that may run in an environment without
    CuPy or a visible GPU, so every probe is wrapped: a missing library or no
    device yields the literal ``"unavailable"`` rather than crashing the table
    build.  The stamps are recorded verbatim on every ``region_meta`` row so a
    figure caption can cite the exact toolchain that produced the alignment.
    """
    stamps = {
        "python": platform.python_version(),
        "numpy": "unavailable",
        "cupy": "unavailable",
        "cuda_driver": "unavailable",
        "cuda_runtime": "unavailable",
        "gpu": "unavailable",
    }
    try:
        import numpy as np

        stamps["numpy"] = np.__version__
    except Exception:
        pass
    try:
        import cupy as cp

        stamps["cupy"] = cp.__version__
        try:
            drv = cp.cuda.runtime.driverGetVersion()
            stamps["cuda_driver"] = f"{drv // 1000}.{(drv % 1000) // 10}"
        except Exception:
            pass
        try:
            rt = cp.cuda.runtime.runtimeGetVersion()
            stamps["cuda_runtime"] = f"{rt // 1000}.{(rt % 1000) // 10}"
        except Exception:
            pass
        try:
            stamps["gpu"] = cp.cuda.runtime.getDeviceProperties(0)["name"].decode()
        except Exception:
            pass
    except Exception:
        pass
    return stamps


# --------------------------------------------------------------------------- #
# Loaders (tolerant — sibling producers S1/S7 may format slightly differently)
# --------------------------------------------------------------------------- #
def load_stats(path: str) -> dict:
    """Load the combined-amplicons panel stats JSON (``insilico_pcr.py`` output)."""
    with open(path) as fh:
        stats = json.load(fh)
    if "regions" not in stats:
        raise SystemExit(
            f"[region_meta] {path} is not panel-mode stats "
            f"(no 'regions' block); rerun insilico_pcr.py --combined-fasta"
        )
    return stats


def _genome_id_of(rec_id: str) -> str:
    """``genome_id`` from a ``{genome_id}__{region}`` amplicon key (region may hold no '__')."""
    return rec_id.rsplit("__", 1)[0]


def _region_of(rec_id: str) -> str:
    """``region`` from a ``{genome_id}__{region}`` amplicon key."""
    return rec_id.rsplit("__", 1)[-1]


def _norm_label(v) -> str:
    """Normalize a domain/phylum label to a non-empty display string."""
    if v is None:
        return "Unknown"
    s = str(v).strip()
    return s if s else "Unknown"


def _manifest_records(man):
    """Yield per-genome record dicts from any S1 selection-manifest shape.

    ``select_10k.py`` writes ``{"meta": {...}, "genomes": [record, ...]}`` (each
    record carries lower-case ``domain``/``phylum``/... top-level), so the canonical
    case is "dict with a ``genomes`` list".  For robustness we also accept a bare
    ``[record, ...]`` list and a ``{genome_id: record}`` mapping (the latter only if
    its values look like per-genome records, never the ``meta`` sub-dict).
    """
    if isinstance(man, dict) and isinstance(man.get("genomes"), list):
        yield from man["genomes"]
        return
    if isinstance(man, list):
        yield from man
        return
    if isinstance(man, dict):
        for v in man.values():
            if isinstance(v, dict) and "genome_id" in v:
                yield v


def load_genome_labels(manifest_path: str, truth_path: str) -> dict:
    """Build ``{genome_id: {"domain": str, "phylum": str}}`` for every source genome.

    The S1 selection manifest (``PATHS.selection_json``) is authoritative (it carries
    the full diverse pool's domain/phylum, including MAG/unscoreable genomes).  Its
    canonical shape is ``{"meta": {...}, "genomes": [record, ...]}``; a bare list or a
    ``{genome_id: record}`` mapping are also tolerated.  Rank keys may be lower-case
    (``domain``/``phylum``, the ``_common.RANKS`` convention) or capitalized.  The S7
    truth file backfills any genome the manifest missed (truth stores ``domain``
    top-level and the phylum under ``truth_lineage["Phylum"]``).  Either source alone
    is sufficient; both being absent is a hard error (no way to stratify).
    """
    labels: dict[str, dict[str, str]] = {}

    def _rec_label(rec: dict, key_lower: str, key_cap: str) -> str:
        # accept either casing; tolerate the truth nested-lineage layout too
        if key_lower in rec:
            return _norm_label(rec[key_lower])
        if key_cap in rec:
            return _norm_label(rec[key_cap])
        return "Unknown"

    if manifest_path and C.done(manifest_path):
        with open(manifest_path) as fh:
            man = json.load(fh)
        for rec in _manifest_records(man):
            if not isinstance(rec, dict):
                continue
            gid = rec.get("genome_id")
            if not gid:
                continue
            labels[str(gid)] = {
                "domain": _rec_label(rec, "domain", "Domain"),
                "phylum": _rec_label(rec, "phylum", "Phylum"),
            }

    if truth_path and C.done(truth_path):
        with open(truth_path) as fh:
            truth = json.load(fh)
        for key, rec in truth.items():
            if not isinstance(rec, dict):
                continue
            gid = str(rec.get("src_genome_id") or _genome_id_of(key))
            if gid in labels:
                continue  # manifest already authoritative for this genome
            lin = rec.get("truth_lineage") or {}
            labels[gid] = {
                "domain": _norm_label(rec.get("domain")),
                "phylum": _norm_label(lin.get("Phylum") or lin.get("phylum")),
            }

    if not labels:
        raise SystemExit(
            "[region_meta] no domain/phylum labels available — need at least one of\n"
            f"        {manifest_path}  (S1 selection_10k.json)\n"
            f"        {truth_path}  (S7 benchmark_truth.json)"
        )
    return labels


def amplified_pairs(fanout_path: str, combined_fasta: str):
    """Yield distinct amplified ``(genome_id, region)`` pairs.

    Prefers the S4 fan-out (``{insert_md5: [genome_id__region, ...]}``) since it is
    compact and already the canonical "what amplified" record; falls back to a
    single streaming pass over the combined FASTA headers if the fan-out is absent.
    Each ``(genome, region)`` is yielded once (the fan-out membership is already a
    set of distinct amplicon keys; the FASTA writer emits one record per amplifying
    pair).
    """
    if fanout_path and C.done(fanout_path):
        with open(fanout_path) as fh:
            fanout = json.load(fh)
        for members in fanout.values():
            for key in members:
                yield _genome_id_of(key), _region_of(key)
        return

    if combined_fasta and C.done(combined_fasta):
        with open(combined_fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    rid = line[1:].strip().split(None, 1)[0]
                    yield _genome_id_of(rid), _region_of(rid)
        return

    raise SystemExit(
        "[region_meta] cannot enumerate amplified amplicons — need one of\n"
        f"        {fanout_path}  (S4 amplicon_md5_fanout.json)\n"
        f"        {combined_fasta}  (S3 combined_amplicons.fasta)"
    )


# --------------------------------------------------------------------------- #
# Core aggregation
# --------------------------------------------------------------------------- #
def panel_domain(panel_by_name: dict, region: str) -> str:
    """Panel ``domain`` scope for a region (``Bacteria`` / ``Bact+Arch``)."""
    r = panel_by_name.get(region)
    return r.get("domain", "") if r else ""


def build_tables(stats: dict, labels: dict, pairs_iter, panel) -> tuple:
    """Return ``(region_meta_rows, amplification_rows, domains_seen)``.

    * ``region_meta_rows`` — one dict per region (panel order), region-level
      length/count from ``stats`` plus per-domain extract rates derived from the
      amplified-pairs × genome-labels join.
    * ``amplification_rows`` — long ``(region, domain, phylum)`` extract-rate rows.
    """
    panel_by_name = {r["name"]: r for r in panel}

    # ---- denominators: source-genome counts per domain and per (domain, phylum)
    n_source_total = len(labels)
    n_source_by_domain: dict[str, int] = {}
    n_source_by_dp: dict[tuple, int] = {}  # (domain, phylum) -> n source genomes
    for lab in labels.values():
        d, p = lab["domain"], lab["phylum"]
        n_source_by_domain[d] = n_source_by_domain.get(d, 0) + 1
        n_source_by_dp[(d, p)] = n_source_by_dp.get((d, p), 0) + 1

    # ---- numerators: amplified (genome, region) joined to labels
    # n_amp_by_region_domain[(region, domain)] and ...[(region, domain, phylum)]
    n_amp_rd: dict[tuple, int] = {}
    n_amp_rdp: dict[tuple, int] = {}
    n_amp_unlabeled = 0  # amplified genome not in the label set (provenance gap)
    for gid, region in pairs_iter:
        lab = labels.get(gid)
        if lab is None:
            n_amp_unlabeled += 1
            continue
        d, p = lab["domain"], lab["phylum"]
        n_amp_rd[(region, d)] = n_amp_rd.get((region, d), 0) + 1
        n_amp_rdp[(region, d, p)] = n_amp_rdp.get((region, d, p), 0) + 1

    domains_seen = sorted(n_source_by_domain)
    regions_stats = stats.get("regions", {})
    # Panel order is the canonical region order; append any stats-only regions.
    region_order = [r["name"] for r in panel if r["name"] in regions_stats]
    region_order += [r for r in regions_stats if r not in region_order]

    stamps = version_stamps()

    region_meta_rows = []
    for region in region_order:
        rs = regions_stats.get(region, {})
        amp_len = rs.get("amplicon_len") or {}
        # n_amplicons / extract_rate: prefer the upstream stats (authoritative,
        # md5-byte-consistent); fall back to the join sum if stats omitted them.
        joined_total = sum(v for (r, _d), v in n_amp_rd.items() if r == region)
        n_amplicons = rs.get("n_amplicons")
        if n_amplicons is None:
            n_amplicons = joined_total
        n_input = stats.get("n_input_records") or n_source_total
        extract_rate = rs.get("extract_rate")
        if extract_rate is None:
            extract_rate = round(n_amplicons / n_input, 4) if n_input else 0.0

        # primer pair + runtime revcomp + insert bounds: prefer panel, fall back
        # to whatever the stats recorded.
        pr = panel_by_name.get(region, {})
        if pr:
            fwd = f"{pr.get('fwd_name', '')} {pr.get('fwd', '')}".strip()
            rev = f"{pr.get('rev_name', '')} {pr.get('rev', '')}".strip()
            rev_revcomp = C.revcomp(pr["rev"]) if pr.get("rev") else ""
            insert_lo = pr.get("insert_lo")
            insert_hi = pr.get("insert_hi")
            ecoli = pr.get("ecoli", "")
            priority = pr.get("priority", "")
            panel_dom = pr.get("domain", "")
        else:
            fwd = rs.get("fwd", "")
            rev = rs.get("rev", "")
            rev_revcomp = rs.get("rev_revcomp", "")
            bounds = rs.get("insert_bounds") or [None, None]
            insert_lo, insert_hi = (bounds + [None, None])[:2]
            ecoli = rs.get("ecoli_coords", "")
            priority = rs.get("priority", "")
            panel_dom = ""

        row = {
            "region": region,
            "priority": priority,
            "panel_domain_scope": panel_dom,
            "ecoli_coords": ecoli,
            "fwd_primer": fwd,
            "rev_primer": rev,
            "rev_revcomp_on_sense": rev_revcomp,
            "insert_lo": insert_lo,
            "insert_hi": insert_hi,
            "amp_len_median": amp_len.get("median"),
            "amp_len_mean": amp_len.get("mean"),
            "amp_len_p5": amp_len.get("p5"),
            "amp_len_p95": amp_len.get("p95"),
            "amp_len_min": amp_len.get("min"),
            "amp_len_max": amp_len.get("max"),
            "n_amplicons": n_amplicons,
            "n_source_genomes": n_source_total,
            "extract_rate": extract_rate,
        }
        # per-domain extract_rate columns (one pair of cols per domain seen).
        # PREFER the upstream stats' regions[r].by_domain[d] block — it already holds
        # the exact per-domain (n_amplicons, extract_rate) with the correct per-domain
        # input-record denominator, so these columns stay byte-consistent with the run.
        # Fall back to the amplified-pairs x labels join only when the stats predate the
        # by_domain split (older insilico_pcr.py) or omit a domain.
        stats_by_dom = rs.get("by_domain") or {}
        for d in domains_seen:
            sd = stats_by_dom.get(d)
            if sd is not None and sd.get("extract_rate") is not None:
                n_amp = sd.get("n_amplicons", 0)
                rate = sd.get("extract_rate", 0.0)
            else:
                n_amp = n_amp_rd.get((region, d), 0)
                n_src = n_source_by_domain.get(d, 0)
                rate = round(n_amp / n_src, 4) if n_src else 0.0
            row[f"n_amplicons_{d}"] = n_amp
            row[f"extract_rate_{d}"] = rate
        # version stamps
        row["py_version"] = stamps["python"]
        row["numpy_version"] = stamps["numpy"]
        row["cupy_version"] = stamps["cupy"]
        row["cuda_driver"] = stamps["cuda_driver"]
        row["cuda_runtime"] = stamps["cuda_runtime"]
        row["gpu"] = stamps["gpu"]
        region_meta_rows.append(row)

    # ---- amplification_by_phylum: long (region, domain, phylum) rows
    amplification_rows = []
    for region in region_order:
        for (d, p), n_src in sorted(n_source_by_dp.items()):
            n_amp = n_amp_rdp.get((region, d, p), 0)
            rate = round(n_amp / n_src, 4) if n_src else 0.0
            amplification_rows.append({
                "region": region,
                "domain": d,
                "phylum": p,
                "n_amplified": n_amp,
                "n_source": n_src,
                "extract_rate": rate,
            })

    return region_meta_rows, amplification_rows, domains_seen, n_amp_unlabeled


# --------------------------------------------------------------------------- #
# CSV writers (atomic)
# --------------------------------------------------------------------------- #
def write_csv(path: str, rows: list, fieldnames: list) -> None:
    """Write ``rows`` as CSV with ``fieldnames`` header (atomic via .tmp)."""
    tmp = path + ".tmp"
    with open(tmp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    os.replace(tmp, path)


def region_meta_fieldnames(domains_seen: list) -> list:
    """Stable column order for region_meta.csv (per-domain cols interleaved)."""
    base = [
        "region", "priority", "panel_domain_scope", "ecoli_coords",
        "fwd_primer", "rev_primer", "rev_revcomp_on_sense",
        "insert_lo", "insert_hi",
        "amp_len_median", "amp_len_mean", "amp_len_p5", "amp_len_p95",
        "amp_len_min", "amp_len_max",
        "n_amplicons", "n_source_genomes", "extract_rate",
    ]
    for d in domains_seen:
        base += [f"n_amplicons_{d}", f"extract_rate_{d}"]
    base += ["py_version", "numpy_version", "cupy_version",
             "cuda_driver", "cuda_runtime", "gpu"]
    return base


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="rebuild even if both output CSVs already exist")
    args = ap.parse_args()

    out_meta = C.PATHS.region_meta_csv
    out_amp = C.PATHS.amplification_by_phylum_csv

    # Resumability: short-circuit if both checkpoints already exist.
    if not args.force and C.done(out_meta) and C.done(out_amp):
        print("[region_meta] outputs present — skipping (use --force to rebuild)")
        print(f"[region_meta]   {out_meta}")
        print(f"[region_meta]   {out_amp}")
        return 0

    stats_path = C.PATHS.combined_stats
    if not C.done(stats_path):
        raise SystemExit(
            f"[region_meta] missing input {stats_path}\n"
            f"        run S3 first: insilico_pcr.py --combined-fasta {C.PATHS.combined_fasta}"
        )

    stats = load_stats(stats_path)
    panel = C.load_panel()
    labels = load_genome_labels(C.PATHS.selection_json, C.PATHS.truth_json)
    pairs = amplified_pairs(C.PATHS.amplicon_fanout, C.PATHS.combined_fasta)

    region_meta_rows, amplification_rows, domains_seen, n_unlabeled = build_tables(
        stats, labels, pairs, panel
    )

    write_csv(out_meta, region_meta_rows, region_meta_fieldnames(domains_seen))
    write_csv(out_amp, amplification_rows,
              ["region", "domain", "phylum", "n_amplified", "n_source", "extract_rate"])

    print(f"[region_meta] {len(region_meta_rows)} regions x {len(domains_seen)} domain(s) "
          f"({', '.join(domains_seen)})")
    if n_unlabeled:
        print(f"[region_meta] WARNING: {n_unlabeled:,} amplified (genome,region) pairs had "
              f"no domain/phylum label (source genome absent from manifest+truth)")
    print(f"[region_meta] wrote {out_meta}  ({len(region_meta_rows)} rows)")
    print(f"[region_meta] wrote {out_amp}  ({len(amplification_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
