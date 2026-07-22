"""Build the hermetic fixtures the embedding validators read.

`validate_analysis.py` and `validate_query.py` were written against an
out-of-band store at ``/tmp/prFBA_test`` plus ``/tmp/sample_fasta.json`` that no
committed script produced, so they could not run from a clean checkout.  This
script regenerates those fixtures deterministically:

  * a small synthetic 16S-like sequence set (random ACGT) with a few
    near-duplicate clusters (so kNN / density have real structure) and both
    full-length (1400-1600 bp) and short (200-350 bp) members, which
    ``validate_query`` requires;
  * embedded with the SAME encoder the pipeline uses (``nt_embed``), so a query
    re-embedded by ``query.py`` matches its stored vector;
  * written as a ``query.py``-compatible store:
        /tmp/prFBA_test/embeddings.f16.npy   (N, 1024) L2-normalized float16
        /tmp/prFBA_test/index.parquet        row-aligned metadata
        /tmp/prFBA_test/members.parquet      feature_id -> row_id
        /tmp/sample_fasta.json               {header: sequence}

The sequences are synthetic, not real 16S: the validators check numerical
correctness (GPU kNN == CPU brute force, PCA == numpy, cosine == brute force,
exact-md5 semantics, novelty flag), none of which needs biological realism.

Requires the encoder + GPU (the validators need them too).  Deterministic
(seed 0); safe to re-run.

    python tests/make_fixtures.py            # (re)build the fixtures
    python tests/make_fixtures.py --force    # rebuild even if present

Interpreter: /home/freiburger/Documents/py_venv/bin/python
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import nt_embed  # noqa: E402

STORE = "/tmp/prFBA_test"
SAMPLE_FASTA = "/tmp/sample_fasta.json"
SEED = 0
BASES = np.array(list("ACGT"))


def _rand_seq(rng, n):
    return "".join(rng.choice(BASES, size=n))


def _mutate(rng, seq, n_sub):
    """A near-duplicate: n_sub random substitutions (gives tight kNN clusters)."""
    s = list(seq)
    for pos in rng.choice(len(s), size=min(n_sub, len(s)), replace=False):
        s[pos] = rng.choice(BASES)
    return "".join(s)


def build_sequences():
    """Deterministic sequence set with clusters + the length bands the tests need.

    The FIRST full-length (1400-1600 bp) and FIRST short (200-350 bp) sequences
    are ISOLATED probes with no near-duplicates: ``validate_query`` picks these as
    its exact / short queries and needs each to be uniquely its own rank-1 hit
    (near-twins would let fp16 embedding noise displace rank-1, breaking the
    exact-md5 and batch==single checks).  The clustered sequences that follow give
    ``validate_analysis`` its dense kNN / radius structure.
    """
    rng = np.random.default_rng(SEED)
    seqs = {}
    # Probe FIRST (uniquely its own top-1): a full-length sequence + a graded
    # ladder of relatives at 2/4/6/9% divergence.  validate_query picks this as
    # its 'exact' query and compares the WHOLE top-5 between batch and single
    # invocations, so ranks 2-5 must be a deterministic ladder (well-separated in
    # cosine), not a background near-tie that fp16 matmul noise would reorder.
    probe_full = _rand_seq(rng, 1500)
    seqs["probeFull_iso"] = probe_full             # -> validate_query 'exact' (rank-1 = itself)
    # widely-spaced divergences so the cosine gaps between rungs exceed the ~few-e-3
    # fp16 batch-shape noise (the exact query is embedded in a 4-seq mixed-length
    # batch vs alone), keeping the top-5 order identical batch-vs-single.
    for pct in (5, 12, 20, 30):
        seqs[f"probeFull_r{pct:02d}"] = _mutate(rng, probe_full, int(1500 * pct / 100))
    seqs["probeShort_iso"] = _rand_seq(rng, 300)   # -> validate_query 'short'
    # 12 near-duplicate clusters spanning the length bands, each a base sequence
    # plus 15-25 low-divergence variants -> dense NN for the geometry checks.
    lengths = [1500, 1550, 1450, 1420, 1580,      # full-length (1400-1600)
               300, 250, 340, 220,                # short (200-350)
               900, 1100, 700]                    # mid
    for ci, L in enumerate(lengths):
        base = _rand_seq(rng, L)
        seqs[f"clust{ci}_base"] = base
        for j in range(int(rng.integers(15, 26))):
            n_sub = int(rng.integers(1, max(2, L // 60)))   # ~1.5% divergence
            seqs[f"clust{ci}_v{j}"] = _mutate(rng, base, n_sub)
    # a handful of isolated random sequences (sparse frontier)
    for k in range(12):
        seqs[f"iso{k}"] = _rand_seq(rng, int(rng.integers(250, 1550)))
    return seqs


def main():
    force = "--force" in sys.argv
    if not force and os.path.exists(f"{STORE}/embeddings.f16.npy") and os.path.exists(SAMPLE_FASTA):
        print(f"[make_fixtures] fixtures already present ({STORE}, {SAMPLE_FASTA}); "
              f"use --force to rebuild.")
        return

    os.makedirs(STORE, exist_ok=True)
    seqs = build_sequences()
    headers = list(seqs)
    seq_list = [seqs[h] for h in headers]
    print(f"[make_fixtures] {len(seq_list)} synthetic sequences "
          f"({sum(1400 <= len(s) <= 1600 for s in seq_list)} full-length, "
          f"{sum(200 <= len(s) <= 350 for s in seq_list)} short); embedding ...")

    tok, model = nt_embed.load_model()
    # embed in modest batches (full-length ~260 tokens); L2-normalized already
    vecs = []
    B = 64
    for i in range(0, len(seq_list), B):
        v = nt_embed.embed_batch(seq_list[i:i + B], tok, model).float().cpu().numpy()
        vecs.append(v)
    E = np.concatenate(vecs, 0).astype(np.float16)

    np.save(f"{STORE}/embeddings.f16.npy", E)
    rows = []
    members = []
    for r, (h, s) in enumerate(zip(headers, seq_list)):
        md5 = hashlib.md5(s.encode()).hexdigest()
        L = len(s)
        rows.append(dict(
            row_id=r, md5=md5, seq_len=L, multiplicity=1, n_genomes=1,
            feature_id=f"fig|{r}.rna.1", organism=h.split("_")[0],
            genome_id=f"{r}.1", taxon_id=r, genome_name=h,
            atypical_len=bool(L > 2000 or L < 200)))
        members.append(dict(feature_id=f"fig|{r}.rna.1", row_id=r, seq_len=L))
    pd.DataFrame(rows).to_parquet(f"{STORE}/index.parquet")
    pd.DataFrame(members).to_parquet(f"{STORE}/members.parquet")
    json.dump(seqs, open(SAMPLE_FASTA, "w"))

    print(f"[make_fixtures] wrote:\n  {STORE}/embeddings.f16.npy  {E.shape} {E.dtype}\n"
          f"  {STORE}/index.parquet  ({len(rows)} rows)\n"
          f"  {STORE}/members.parquet\n  {SAMPLE_FASTA}")
    print("[make_fixtures] now run:  python tests/validate_analysis.py "
          "&& python tests/validate_query.py")


if __name__ == "__main__":
    main()
