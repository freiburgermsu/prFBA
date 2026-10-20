#!/usr/bin/env python
"""S7 — build the per-amplicon benchmark truth table (DESIGN.md §6).

For every amplicon key ``{genome_id}__{region}`` emitted by the in-silico PCR
(panel mode of ``insilico_pcr.py``), resolve the SOURCE organism's true lineage
and record it as ground truth for taxonomic-accuracy scoring.

    src_taxon_id  = int(genome_id.split('.')[0])
    truth_lineage = _common.lineage_for(src_taxon_id)        # 7 capitalized ranks
    domain        = the taxopy ``domain`` rank (Bacteria / Archaea)

Output (``PATHS.truth_json``):

    {
      "{genome_id}__{region}": {
        "src_genome_id": "<genome_id>",
        "src_taxon_id":  <int|null>,
        "region":        "<region>",
        "domain":        "Bacteria" | "Archaea" | null,
        "truth_lineage": {Kingdom, Phylum, Class, Order, Family, Genus, Species}
      },
      ...
    }

Per DESIGN §6: non-numeric / MAG ``genome_id``s (no leading integer taxon id)
get an all-None lineage + null taxon id + null domain.  They are written (so the
table is complete) but logged, and are excluded later from the benchmark as
unscoreable — that is NOT a pipeline failure.

Amplicon keys come from ``PATHS.combined_fasta`` headers by default; if that is
absent the script falls back to ``PATHS.selection_json`` (the 10k manifest) and
expands each source genome over all 8 panel regions.  Both yield the identical
``{genome_id}__{region}`` key space; the FASTA is preferred because it is the
exact set of amplicons that actually amplified.

Resumable: if ``PATHS.truth_json`` already exists and is non-empty, the script
short-circuits unless ``--force`` is given.

Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.

I/O contract
------------
Reads (one of, in preference order):
    PATHS.combined_fasta   = .../data/combined_amplicons.fasta
    PATHS.selection_json   = .../data/selection_10k.json   (fallback)
  + the prFBA taxdump (TAXDUMP_NODES / TAXDUMP_NAMES) via _common.lineage_for.
Writes:
    PATHS.truth_json       = .../data/benchmark_truth.json
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys

# Make `import _common` work regardless of CWD (agent threads reset cwd).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import _common as C  # noqa: E402
from _common import PATHS, RANKS, done, lineage_for, load_panel  # noqa: E402

# A BV-BRC genome_id is "<taxon_id>.<assembly_index>", both integers.
_GENOME_ID_RE = re.compile(r"^\d+\.\d+$")


def _src_taxon_id(genome_id: str):
    """Leading integer of a ``<taxon>.<asm>`` genome id, or None for MAG/non-numeric ids."""
    if not _GENOME_ID_RE.match(genome_id):
        return None
    try:
        return int(genome_id.split(".")[0])
    except (ValueError, IndexError):
        return None


def _domain_for(taxon_id):
    """Resolve the taxopy ``domain`` (Bacteria/Archaea) for a taxon id.

    ``_common.lineage_for`` puts the *kingdom* (e.g. 'Pseudomonadati') in the
    ``Kingdom`` slot, which is NOT the domain — so the Bacteria/Archaea
    stratification key (DESIGN §1, §6) is read straight off the ``domain`` rank
    (``superkingdom`` fallback, mirroring select_10k §2).  Reuses the cached
    TaxDb singleton; any failure / None id -> None.
    """
    if taxon_id is None:
        return None
    try:
        import taxopy

        rd = taxopy.Taxon(int(taxon_id), C._get_taxdb()).rank_name_dictionary
        return rd.get("domain") or rd.get("superkingdom")
    except Exception:
        return None


def _amplicon_keys_from_fasta(path: str):
    """Yield ``{genome_id}__{region}`` keys from a combined-amplicons FASTA's headers."""
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].split(None, 1)[0].strip()


def _amplicon_keys_from_selection(path: str):
    """Yield keys by crossing the 10k source genomes with all panel region names.

    Fallback when the combined FASTA has not been built yet.  The select_10k.py
    manifest is ``{"meta": {...}, "genomes": [{"genome_id": ...}, ...]}``; we also
    tolerate a bare list of records and a ``{genome_id: record}`` dict.  In every
    case we recover the genome ids and expand over the panel.  This is a superset
    of the FASTA key space (it includes non-amplifying (genome,region) pairs),
    which is acceptable: those rows still carry a valid truth lineage and are
    simply never referenced by a real amplicon during scoring.
    """
    with open(path) as fh:
        sel = json.load(fh)

    def _gid(rec):
        return rec.get("genome_id") if isinstance(rec, dict) else rec

    if isinstance(sel, dict) and isinstance(sel.get("genomes"), list):
        # select_10k.py manifest shape: {"meta": ..., "genomes": [record, ...]}.
        genome_ids = [_gid(rec) for rec in sel["genomes"]]
    elif isinstance(sel, list):
        genome_ids = [_gid(rec) for rec in sel]
    elif isinstance(sel, dict):
        # {genome_id: record} mapping.
        genome_ids = list(sel.keys())
    else:
        raise SystemExit(f"unrecognized selection JSON shape in {path}: {type(sel)}")

    region_names = [r["name"] for r in load_panel()]
    for gid in genome_ids:
        if not gid:
            continue
        for region in region_names:
            yield f"{gid}__{region}"


def _iter_amplicon_keys():
    """Choose the amplicon key source (FASTA preferred, selection fallback)."""
    if done(PATHS.combined_fasta):
        return PATHS.combined_fasta, _amplicon_keys_from_fasta(PATHS.combined_fasta)
    if done(PATHS.selection_json):
        return PATHS.selection_json, _amplicon_keys_from_selection(PATHS.selection_json)
    raise SystemExit(
        "no amplicon key source found — expected one of:\n"
        f"  {PATHS.combined_fasta}  (combined-amplicons FASTA, preferred)\n"
        f"  {PATHS.selection_json}  (10k selection manifest, fallback)\n"
        "Run insilico_pcr.py --combined-fasta (S3) or select_10k.py (S1) first."
    )


def build_truth(force: bool = False) -> dict:
    """Build (or load) the benchmark truth table; resumable."""
    if done(PATHS.truth_json) and not force:
        with open(PATHS.truth_json) as fh:
            truth = json.load(fh)
        print(
            f"[truth] {PATHS.truth_json} already exists "
            f"({len(truth):,} amplicons) — skipping (use --force to rebuild)",
            flush=True,
        )
        return truth

    src_path, keys = _iter_amplicon_keys()
    print(f"[truth] reading amplicon keys from {src_path}", flush=True)

    # Per-genome caches so each source genome resolves its lineage/domain once,
    # even though it owns up to 8 amplicon keys (one per region).
    lineage_cache: dict = {}
    domain_cache: dict = {}

    truth: dict = {}
    n_amp = 0
    n_mag = 0          # non-numeric / MAG genome ids (all-None lineage)
    n_no_lineage = 0   # numeric id whose lineage_for returned all-None
    mag_examples: list = []

    for key in keys:
        n_amp += 1
        # Region name may contain '-' / '_' but never '__'; split on first '__'.
        if "__" not in key:
            print(f"[truth] WARNING: malformed amplicon key (no '__'): {key!r}", flush=True)
            continue
        genome_id, region = key.split("__", 1)

        taxon_id = _src_taxon_id(genome_id)

        if genome_id not in lineage_cache:
            lineage_cache[genome_id] = lineage_for(taxon_id)
            domain_cache[genome_id] = _domain_for(taxon_id)
        lineage = lineage_cache[genome_id]
        domain = domain_cache[genome_id]

        if taxon_id is None:
            n_mag += 1
            if len(mag_examples) < 10:
                mag_examples.append(genome_id)
        elif all(v is None for v in lineage.values()):
            n_no_lineage += 1

        truth[key] = {
            "src_genome_id": genome_id,
            "src_taxon_id": taxon_id,
            "region": region,
            "domain": domain,
            "truth_lineage": lineage,
        }

    # Domain tally over scoreable (lineage-resolved) amplicons.
    dom_counts: dict = {}
    for rec in truth.values():
        dom_counts[rec["domain"]] = dom_counts.get(rec["domain"], 0) + 1

    tmp = PATHS.truth_json + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(truth, fh)
    os.replace(tmp, PATHS.truth_json)

    print(
        f"[truth] wrote {len(truth):,} amplicon truth rows to {PATHS.truth_json}\n"
        f"[truth]   distinct source genomes: {len(lineage_cache):,}\n"
        f"[truth]   domain split: "
        + ", ".join(f"{k}={v:,}" for k, v in sorted(dom_counts.items(), key=lambda x: (x[0] is None, x[0])))
        + f"\n[truth]   MAG/non-numeric ids -> all-None lineage (excluded later): {n_mag:,}"
        + (f" e.g. {mag_examples}" if mag_examples else "")
        + f"\n[truth]   numeric ids that failed to resolve any rank: {n_no_lineage:,}",
        flush=True,
    )
    assert set(next(iter(truth.values()))["truth_lineage"]) == set(C._LINEAGE_KEYS), (
        "truth_lineage keys must be the 7 capitalized lineage ranks"
    )
    return truth


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--force",
        action="store_true",
        help="rebuild even if PATHS.truth_json already exists",
    )
    args = ap.parse_args()
    build_truth(force=args.force)


if __name__ == "__main__":
    main()
