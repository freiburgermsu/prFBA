#!/usr/bin/env python
"""S4 — collapse byte-identical amplicon inserts before alignment (DESIGN §5.0).

The single highest-leverage alignment speedup: conserved primer flanks make many
``(genome, region)`` amplicons byte-identical across source genomes, so the
aligner only needs to see each distinct insert once.  This stage reads the
combined panel FASTA (one record per amplifying ``(genome, region)`` written by
``insilico_pcr.py --combined-fasta``), groups records by the md5 of their insert
sequence, and emits:

  * ``PATHS.unique_fasta``    — one record per unique insert, id = the insert md5,
                                so the aligner keys on the md5 (deterministic,
                                provenance-free, collision-safe).
  * ``PATHS.amplicon_fanout`` — JSON ``{insert_md5: [genome_id__region, ...]}``
                                so alignment results fan back out to every source
                                ``(genome, region)`` that produced that insert.

Resumable: if both outputs already exist (and ``--force`` is not given) the stage
short-circuits.  Sole interpreter: ``/home/freiburger/Documents/py_venv/bin/python``.

I/O contract
------------
reads : PATHS.combined_fasta   (data/combined_amplicons.fasta)
writes: PATHS.unique_fasta     (data/combined_amplicons.unique.fasta)
        PATHS.amplicon_fanout  (data/amplicon_md5_fanout.json)
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys

# Import the shared contract (paths/constants/helpers) — never hardcode a path.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C  # noqa: E402


def iter_fasta(path):
    """Yield ``(record_id, sequence)`` from a FASTA file.

    ``record_id`` is the header token up to the first whitespace (the panel
    writer emits ``>{genome_id}__{region}`` with no description).  The sequence
    is the concatenation of all wrapped lines, uppercased and whitespace-free, so
    the md5 is identical whether or not ``insilico_pcr.py`` was run with ``--wrap``.
    """
    rec_id = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if rec_id is not None:
                    yield rec_id, "".join(chunks)
                rec_id = line[1:].strip().split(None, 1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if rec_id is not None:
        yield rec_id, "".join(chunks)


def dedup(combined_fasta: str):
    """Group combined-FASTA records by insert md5.

    Returns ``(fanout, md5_to_seq, n_records)`` where ``fanout`` is
    ``{md5: [genome_id__region, ...]}`` and ``md5_to_seq`` is ``{md5: seq}``.
    Membership lists and md5 iteration order are kept sorted for byte-identical
    (reproducible) outputs regardless of FASTA record order.
    """
    fanout: dict[str, list[str]] = {}
    md5_to_seq: dict[str, str] = {}
    n_records = 0
    for rec_id, seq in iter_fasta(combined_fasta):
        if not seq:
            continue  # skip empty/malformed record defensively
        seq = seq.upper()
        n_records += 1
        md5 = hashlib.md5(seq.encode("ascii")).hexdigest()
        fanout.setdefault(md5, []).append(rec_id)
        md5_to_seq.setdefault(md5, seq)
    for md5 in fanout:
        fanout[md5].sort()
    return fanout, md5_to_seq, n_records


def write_unique_fasta(path: str, md5_to_seq: dict[str, str]) -> None:
    """Write one record per unique insert: ``>{md5}\\n{seq}\\n`` (md5-sorted)."""
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        for md5 in sorted(md5_to_seq):
            fh.write(f">{md5}\n{md5_to_seq[md5]}\n")
    os.replace(tmp, path)


def write_fanout(path: str, fanout: dict[str, list[str]]) -> None:
    """Write the ``{insert_md5: [genome_id__region, ...]}`` JSON (md5-sorted, atomic)."""
    tmp = path + ".tmp"
    ordered = {md5: fanout[md5] for md5 in sorted(fanout)}
    with open(tmp, "w") as fh:
        json.dump(ordered, fh, indent=2)
    os.replace(tmp, path)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--force", action="store_true",
                    help="rebuild even if both outputs already exist")
    args = ap.parse_args()

    in_fasta = C.PATHS.combined_fasta
    out_fasta = C.PATHS.unique_fasta
    out_fanout = C.PATHS.amplicon_fanout

    # Resumability: short-circuit if both checkpoints already exist.
    if not args.force and C.done(out_fasta) and C.done(out_fanout):
        with open(out_fanout) as fh:
            n_unique = len(json.load(fh))
        print(f"[dedup] outputs present — skipping (use --force to rebuild)")
        print(f"[dedup]   {out_fasta}")
        print(f"[dedup]   {out_fanout}  ({n_unique:,} unique inserts)")
        return 0

    if not C.done(in_fasta):
        raise SystemExit(
            f"[dedup] missing input {in_fasta}\n"
            f"        run S3 first: insilico_pcr.py --combined-fasta {in_fasta}"
        )

    fanout, md5_to_seq, n_records = dedup(in_fasta)
    n_unique = len(md5_to_seq)

    write_unique_fasta(out_fasta, md5_to_seq)
    write_fanout(out_fanout, fanout)

    reduction = (n_records / n_unique) if n_unique else 0.0
    print(f"[dedup] {n_records:,} amplicon records -> {n_unique:,} unique inserts "
          f"({reduction:.2f}x reduction)")
    print(f"[dedup] wrote {out_fasta}")
    print(f"[dedup] wrote {out_fanout}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
