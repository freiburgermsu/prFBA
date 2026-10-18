#!/usr/bin/env python
"""Select a 10,000-genome panel by UNIFORM RANDOM SAMPLING of admissible genomes.

The third arm of the panel-composition comparison. The other two deliberately reshape
the database's taxonomic distribution in opposite directions:

  * select_10k.py          taxonomy-stratified, sqrt-flattening availability, which
                           over-weights rare phyla and archaea relative to the database.
  * select_kmeans.py       embedding-uniform: k-means over 16S embedding space, so the
                           panel is uniform over SEQUENCE space rather than over taxa.
  * select_random.py (this) uniform over admissible GENOMES, so the panel inherits
                           BV-BRC's own taxonomic skew unchanged.

That last property is the point: a purely random panel is the null composition against
which the other two schemes' corrections are measured. It is expected to be dominated by
the heavily-sequenced genera and to contain few archaea, and the comparison tables should
show exactly that.

Admissibility is IDENTICAL to select_kmeans.py (and therefore to select_10k.py's gate):
full-length (>= 1400 bp, non-atypical) 16S, domain in {Archaea, Bacteria}, named phylum.
Sampling is over genomes, not sequences, so a genome with many operons is no likelier to
be drawn than one with a single operon -- the opposite of what clustering over sequences
does.

The representative 16S for each drawn genome follows select_10k.py's DESIGN rule verbatim
(longest sequence; tie -> lexicographically-smallest md5) rather than k-means' centroid-
nearest rule, because that rule is a property of the clustering and has no analogue here.

Emits selection_random.json and src16s_random.fasta in exactly the schema
run_validation.sh consumes, so the rest of the 16S pipeline runs unchanged.
"""
import argparse, hashlib, json, os, sys, time
import numpy as np
import pandas as pd

ROOT = os.path.expanduser("~/Documents/prFBA")
sys.path.insert(0, os.path.join(ROOT, "region_validation", "scripts"))
import _common                                             # noqa: F401  (taxdb bootstrap)

SEED = 1729
FULL_LEN = 1400
ALLOWED_DOMAINS = ("Archaea", "Bacteria")
MD5_SEQ = "/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json"
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]


def log(m):
    print("[random] %s" % m, flush=True)


_LIN = {}


def lineage(taxon_id):
    """select_kmeans.py's lineage rule, verbatim, so the admissible pool is identical
    across the two alternative panels and any difference between them is the sampling
    scheme alone."""
    if taxon_id in _LIN:
        return _LIN[taxon_id]
    import taxopy
    try:
        rd = taxopy.Taxon(int(taxon_id), _common._get_taxdb()).rank_name_dictionary
    except Exception:
        _LIN[taxon_id] = None
        return None
    domain = rd.get("domain") or rd.get("superkingdom")
    phylum = rd.get("phylum")
    if domain not in ALLOWED_DOMAINS or not phylum:
        _LIN[taxon_id] = None
        return None
    lin = {"domain": domain, "phylum": phylum}
    pr, pn = "phylum", phylum
    for rank in ("class", "order", "family", "genus"):
        nm = rd.get(rank)
        if nm:
            lin[rank] = nm; pr, pn = rank, nm
        else:
            lin[rank] = "unclassified_%s:%s" % (pr, pn)
    lin["species"] = rd.get("species") or "sp_taxon_%s" % taxon_id
    _LIN[taxon_id] = lin
    return lin


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--outdir", default=os.path.join(ROOT, "random_panel"))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    t0 = time.time()

    df = pd.read_parquet(os.path.join(ROOT, "index.parquet"))
    mask = (df.seq_len >= FULL_LEN) & (~df.atypical_len)
    sub = df[mask].reset_index(drop=True)
    log("full-length: %d sequences over %d genomes" % (len(sub), sub.genome_id.nunique()))

    tids = sorted({int(x) for x in sub.taxon_id.dropna().unique()})
    log("resolving %d distinct taxa ..." % len(tids))
    ok = {t for t in tids if lineage(t) is not None}
    log("admissible taxa (domain in {Archaea,Bacteria} and named phylum): %d/%d"
        % (len(ok), len(tids)))
    sub = sub[sub.taxon_id.notna() & sub.taxon_id.astype("Int64").isin(ok)].reset_index(drop=True)
    log("admissible: %d sequences over %d genomes"
        % (len(sub), sub.genome_id.nunique()))

    # ---- uniform draw over GENOMES (not sequences) ----------------------------------
    pool = np.array(sorted(sub.genome_id.astype(str).unique()))   # sorted => reproducible
    if a.n > len(pool):
        raise SystemExit("asked for %d genomes but only %d are admissible" % (a.n, len(pool)))
    rng = np.random.default_rng(a.seed)
    drawn = rng.choice(len(pool), size=a.n, replace=False)
    drawn.sort()                                      # stable output order
    chosen = pool[drawn]
    log("drew %d genomes uniformly at random from %d admissible (seed %d)"
        % (len(chosen), len(pool), a.seed))

    # ---- representative 16S per genome: longest; tie -> smallest md5 ------------------
    chosen_set = set(chosen)
    g = sub[sub.genome_id.astype(str).isin(chosen_set)].copy()
    g["genome_id"] = g.genome_id.astype(str)
    g["md5"] = g.md5.astype(str)
    g = g.sort_values(["genome_id", "seq_len", "md5"], ascending=[True, False, True])
    best = g.groupby("genome_id", sort=True).first().reset_index()
    nfull = g.groupby("genome_id", sort=True).size().rename("n_full").reset_index()
    best = best.merge(nfull, on="genome_id")
    log("representatives resolved for %d/%d genomes" % (len(best), len(chosen)))

    # index.parquet keys on md5(sequence.upper()); 16S_md5_seq.json keys on md5 of the
    # raw sequence, so bridge by re-hashing (same fix as select_kmeans.py).
    raw = json.load(open(MD5_SEQ))
    seqs = {hashlib.md5(v.upper().encode()).hexdigest(): v for v in raw.values()}
    log("sequence map: %d raw md5 -> %d upper-md5 keys" % (len(raw), len(seqs)))
    del raw

    genomes, fasta, missing = [], [], []
    for r in best.itertuples(index=False):
        seq = seqs.get(str(r.md5))
        if seq is None:
            missing.append(str(r.genome_id))
            continue
        tid = int(r.taxon_id)
        lin = lineage(tid)
        genomes.append(dict(
            genome_id=str(r.genome_id), taxon_id=tid,
            best_16s_md5=str(r.md5), best_16s_len=int(r.seq_len),
            n_full_length_16s=int(r.n_full),
            **{k: lin[k] for k in RANKS},
            seed=a.seed, full_len=FULL_LEN))
        fasta.append(">%s\n%s\n" % (r.genome_id, seq))
    if missing:
        log("WARNING: %d drawn genomes had no retrievable sequence for their best md5: %s"
            % (len(missing), ", ".join(missing[:10])))

    dom = {}
    for e in genomes:
        dom[e["domain"]] = dom.get(e["domain"], 0) + 1
    log("domains: %s" % dom)

    sel = os.path.join(a.outdir, "selection_random.json")
    json.dump(dict(meta=dict(method="uniform_random_over_genomes", n=a.n, seed=a.seed,
                             full_len=FULL_LEN, selected_n=len(genomes),
                             n_admissible_genomes=int(len(pool)),
                             n_admissible_sequences=int(len(sub)),
                             n_missing_sequence=len(missing),
                             rep_rule="longest full-length 16S; tie -> smallest md5",
                             index=os.path.join(ROOT, "index.parquet")),
                   genomes=genomes), open(sel, "w"))
    fa = os.path.join(a.outdir, "src16s_random.fasta")
    open(fa, "w").write("".join(fasta))
    log("wrote %s (%d genomes) and %s  [%.0fs]"
        % (sel, len(genomes), fa, time.time() - t0))


if __name__ == "__main__":
    main()
