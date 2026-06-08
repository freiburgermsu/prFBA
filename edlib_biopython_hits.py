#!/usr/bin/env python
"""
edlib_biopython_hits.py — two-stage alignment mapping of the EmilyKin EBPR 16S
V4-V5 ASVs against the BV-BRC 16S reference database.

Pipeline (per ASV, fully parallelized across all cores):
  Stage 1 — edlib (infix / "HW" mode) edit-distance scan over EVERY unique BV-BRC
            16S reference sequence; keep the top-K1 (default 500) lowest-distance
            references. An adaptive k-band (edlib `k` = current K1-th best distance)
            makes the >99% of clearly-dissimilar references reject in microseconds.
  Stage 2 — Biopython PairwiseAligner local Smith-Waterman rescoring of those K1
            candidates; keep the top-K2 (default 20) by alignment score, then
            re-align each to record % identity, aligned length, and reference coords.

References are the DEDUPLICATED unique sequences from 16S_md5_seq.json (md5 -> seq);
each unique sequence is mapped back to a representative genome via 16S_md5_ID.json,
with taxon_id parsed from the BV-BRC genome_id and a full NCBI lineage attached via
taxopy (prFBA taxdump). Identical reference sequences (same md5) collapse to one hit.

Outputs (in --outdir):
  asv_top20_alignment_hits.json   per-ASV mapping: asv_len, midas_taxonomy, rel_ab,
                                  best_align_score, best_identity, n_edlib_candidates,
                                  and top20[] each with align_score, identity,
                                  aligned_len, edlib_distance, edlib_identity,
                                  organism, genome_id, taxon_id, feature_id, md5,
                                  ref_seq_len, ref_aln_start/end, lineage{Kingdom..Species}
  asv_alignment_summary.csv       one row per ASV (best hit + key fields)
  run_stats.json                  parameters, counts, per-stage CPU time, wall time

Interpreter: ~/Documents/py_venv/bin/python   (per project convention)
"""
from __future__ import annotations
import argparse, csv, heapq, json, os, re, time
from multiprocessing import Pool

import numpy as np

# ----------------------------------------------------------------------------- paths
CODIF = "/home/freiburger/Documents/codiffusion_bioreactor/model_inputs"
DB_MD5_SEQ = os.path.join(CODIF, "16S_md5_seq.json")      # md5 -> unique sequence
DB_MD5_ID  = os.path.join(CODIF, "16S_md5_ID.json")       # md5 -> representative BV-BRC header
EMILYKIN   = "/home/freiburger/Documents/EmilyKin"
ASVS_FASTA = os.path.join(EMILYKIN, "asvs.fasta")
TAXONOMY_CSV = os.path.join(EMILYKIN, "bvbrc_embedding_hits/inputs/taxonomy.csv")  # MiDAS lineage + rel_ab
TAXDUMP_NODES = "/home/freiburger/Documents/prFBA/nodes.dmp"
TAXDUMP_NAMES = "/home/freiburger/Documents/prFBA/names.dmp"

RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
# biopython local (Smith-Waterman) scoring — BLASTN-like nucleotide scheme
MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2.0, -3.0, -5.0, -2.0

# --------------------------------------------------------------- shared (fork) globals
REF_BLOB: bytes = b""          # all uppercased reference sequences concatenated
REF_OFF: np.ndarray = None     # int64[N+1] offsets into REF_BLOB
REF_LEN: np.ndarray = None     # int32[N] lengths
ALIGNER = None                 # per-worker Bio.Align.PairwiseAligner


# ----------------------------------------------------------------------------- loaders
def load_references(min_ref_len: int):
    """Parse 16S_md5_seq.json -> (md5_list, blob bytes, offsets, lengths). Filters short seqs."""
    with open(DB_MD5_SEQ) as fh:
        md5_seq = json.load(fh)
    md5s, offsets, ba = [], [0], bytearray()
    for md5, seq in md5_seq.items():
        if len(seq) < min_ref_len:
            continue
        ba += seq.upper().encode("ascii", "ignore")
        md5s.append(md5)
        offsets.append(len(ba))
    blob = bytes(ba)
    off = np.asarray(offsets, dtype=np.int64)
    lens = np.diff(off).astype(np.int32)
    return md5s, blob, off, lens


def load_asvs(path: str):
    ids, seqs, cur, buf = [], [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                if cur is not None:
                    seqs.append("".join(buf).upper())
                cur = line[1:].split()[0]
                ids.append(cur); buf = []
            else:
                buf.append(line)
    if cur is not None:
        seqs.append("".join(buf).upper())
    return ids, seqs


def load_midas(path: str):
    """seq-hash -> (midas_taxonomy string, rel_ab)."""
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            lineage = " ".join(row[r] for r in RANKS if row.get(r))
            try:
                rel = float(row.get("rel_ab", "") or 0.0)
            except ValueError:
                rel = 0.0
            out[row["seq"]] = (lineage, rel)
    return out


# --------------------------------------------------------------------- BV-BRC metadata
_HDR_TAIL = re.compile(r"\[(.+?)\s*\|\s*([0-9]+\.[0-9]+)\]\s*$")
_HDR_FEAT = re.compile(r"^fig\|([^|]+)\|")


def parse_header(hdr: str):
    """BV-BRC header -> (organism, genome_id, taxon_id, feature_id)."""
    organism, genome_id, taxon_id = None, None, None
    m = _HDR_TAIL.search(hdr)
    if m:
        organism, genome_id = m.group(1).strip(), m.group(2)
        taxon_id = int(genome_id.split(".")[0])
    fm = _HDR_FEAT.search(hdr)
    feature_id = fm.group(1) if fm else None
    return organism, genome_id, taxon_id, feature_id


# --------------------------------------------------------------------------- alignment
def make_aligner():
    from Bio import Align
    a = Align.PairwiseAligner()
    a.mode = "local"
    a.match_score = MATCH
    a.mismatch_score = MISMATCH
    a.open_gap_score = GAP_OPEN
    a.extend_gap_score = GAP_EXTEND
    return a


def init_worker():
    global ALIGNER
    ALIGNER = make_aligner()


def _identity_and_coords(aln):
    """(% identity over aligned columns, n_matches, aligned_len, ref_start, ref_end)."""
    rows = [str(s) for s in aln]
    q, t = rows[0], rows[1]
    matches = sum(1 for a, b in zip(q, t) if a == b and a != "-")
    alen = len(q)
    ref_start = ref_end = None
    try:
        tblocks = aln.aligned[1]
        ref_start = int(tblocks[0][0]); ref_end = int(tblocks[-1][1])
    except Exception:
        pass
    return (matches / alen if alen else 0.0), matches, alen, ref_start, ref_end


def process_asv(task):
    """One ASV through both stages. Returns (asv_index, top20 records, t_edlib, t_bio)."""
    import edlib
    asv_idx, asv_seq, k1, k2 = task
    q = asv_seq.encode("ascii", "ignore")
    n = REF_LEN.shape[0]

    # ---- Stage 1: edlib infix scan, adaptive-k top-k1 (min-distance) ----
    t0 = time.perf_counter()
    heap = []          # max-heap of (-dist, ref_idx); size <= k1
    kthr = -1          # edlib k threshold (-1 = unbounded until heap fills)
    worst = None       # current k1-th best distance (= -heap[0][0]) once full
    for i in range(n):
        s = REF_OFF[i]
        ref = REF_BLOB[s:REF_OFF[i + 1]]
        d = edlib.align(q, ref, mode="HW", task="distance", k=kthr)["editDistance"]
        if d == -1:                       # exceeded band -> not a candidate
            continue
        if len(heap) < k1:
            heapq.heappush(heap, (-d, i))
            if len(heap) == k1:
                worst = -heap[0][0]; kthr = worst
        elif d < worst:
            heapq.heapreplace(heap, (-d, i))
            worst = -heap[0][0]; kthr = worst
    cand = sorted(((-nd, idx) for nd, idx in heap))   # (dist, ref_idx) ascending
    t_edlib = time.perf_counter() - t0

    # ---- Stage 2: biopython local rescoring of the k1 candidates ----
    t1 = time.perf_counter()
    asv_str = asv_seq
    scored = []
    for dist, idx in cand:
        ref = REF_BLOB[REF_OFF[idx]:REF_OFF[idx + 1]].decode("ascii")
        sc = ALIGNER.score(asv_str, ref)
        scored.append((sc, dist, idx, ref))
    scored.sort(key=lambda x: (-x[0], x[1], x[2]))     # score desc, edlib dist asc, idx
    out = []
    for sc, dist, idx, ref in scored[:k2]:
        aln = ALIGNER.align(asv_str, ref)[0]
        ident, matches, alen, rs, re_ = _identity_and_coords(aln)
        out.append({
            "ref_idx": int(idx), "align_score": float(sc),
            "identity": round(ident, 4), "n_matches": matches, "aligned_len": alen,
            "ref_aln_start": rs, "ref_aln_end": re_,
            "edlib_distance": int(dist), "ref_seq_len": int(REF_LEN[idx]),
        })
    t_bio = time.perf_counter() - t1
    return asv_idx, out, t_edlib, t_bio


# --------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(EMILYKIN, "bvbrc_alignment_hits"))
    ap.add_argument("--workers", type=int, default=max(1, os.cpu_count() - 4))
    ap.add_argument("--k1", type=int, default=500, help="edlib candidates per ASV")
    ap.add_argument("--k2", type=int, default=20, help="final biopython hits per ASV")
    ap.add_argument("--min-ref-len", type=int, default=200)
    ap.add_argument("--limit", type=int, default=0, help="benchmark on first N ASVs (0=all)")
    ap.add_argument("--chunksize", type=int, default=2)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    global REF_BLOB, REF_OFF, REF_LEN

    wall0 = time.perf_counter()
    print(f"[load] references from {os.path.basename(DB_MD5_SEQ)} ...", flush=True)
    md5s, REF_BLOB, REF_OFF, REF_LEN = load_references(args.min_ref_len)
    n_ref = len(md5s)
    print(f"[load] {n_ref:,} unique refs kept (>= {args.min_ref_len} bp), "
          f"blob {len(REF_BLOB)/1e6:.0f} MB", flush=True)

    asv_ids, asv_seqs = load_asvs(ASVS_FASTA)
    if args.limit:
        asv_ids, asv_seqs = asv_ids[:args.limit], asv_seqs[:args.limit]
    n_asv = len(asv_ids)
    midas = load_midas(TAXONOMY_CSV)
    print(f"[load] {n_asv:,} ASVs | workers={args.workers} | k1={args.k1} k2={args.k2}", flush=True)

    tasks = [(i, asv_seqs[i], args.k1, args.k2) for i in range(n_asv)]
    results = {}
    cpu_edlib = cpu_bio = 0.0
    t_scan = time.perf_counter()
    with Pool(args.workers, initializer=init_worker) as pool:
        done = 0
        for asv_idx, recs, te, tb in pool.imap_unordered(process_asv, tasks, chunksize=args.chunksize):
            results[asv_idx] = recs
            cpu_edlib += te; cpu_bio += tb
            done += 1
            if done % 100 == 0 or done == n_asv:
                el = time.perf_counter() - t_scan
                eta = el / done * (n_asv - done)
                print(f"[scan] {done}/{n_asv}  elapsed={el:6.0f}s  eta={eta:6.0f}s  "
                      f"(cpu edlib={cpu_edlib:,.0f}s bio={cpu_bio:,.0f}s)", flush=True)
    scan_wall = time.perf_counter() - t_scan

    # --------- enrichment: ref_idx -> md5 -> genome/taxonomy/lineage ----------
    print("[enrich] loading md5->header + taxdump ...", flush=True)
    with open(DB_MD5_ID) as fh:
        md5_hdr = json.load(fh)
    import taxopy
    taxdb = taxopy.TaxDb(nodes_dmp=TAXDUMP_NODES, names_dmp=TAXDUMP_NAMES)
    lin_cache: dict[int, dict] = {}

    def lineage_for(taxon_id):
        if taxon_id is None:
            return {r: None for r in RANKS}
        if taxon_id in lin_cache:
            return lin_cache[taxon_id]
        try:
            rd = taxopy.Taxon(taxon_id, taxdb).rank_name_dictionary
            lin = {
                "Kingdom": rd.get("superkingdom") or rd.get("kingdom") or rd.get("domain"),
                "Phylum": rd.get("phylum"), "Class": rd.get("class"), "Order": rd.get("order"),
                "Family": rd.get("family"), "Genus": rd.get("genus"), "Species": rd.get("species"),
            }
        except Exception:
            lin = {r: None for r in RANKS}
        lin_cache[taxon_id] = lin
        return lin

    mapping, summary_rows = {}, []
    for i in range(n_asv):
        asv = asv_ids[i]
        recs = results.get(i, [])
        midas_tax, rel_ab = midas.get(asv, ("", 0.0))
        top = []
        for rank, r in enumerate(recs, 1):
            md5 = md5s[r["ref_idx"]]
            org, gid, taxid, feat = parse_header(md5_hdr.get(md5, ""))
            top.append({
                "rank": rank,
                "align_score": r["align_score"],
                "identity": r["identity"],
                "n_matches": r["n_matches"],
                "aligned_len": r["aligned_len"],
                "edlib_distance": r["edlib_distance"],
                "edlib_identity": round(1 - r["edlib_distance"] / max(1, len(asv_seqs[i])), 4),
                "organism": org,
                "genome_id": gid,
                "taxon_id": taxid,
                "feature_id": feat,
                "md5": md5,
                "ref_seq_len": r["ref_seq_len"],
                "ref_aln_start": r["ref_aln_start"],
                "ref_aln_end": r["ref_aln_end"],
                "lineage": lineage_for(taxid),
            })
        best = top[0] if top else {}
        mapping[asv] = {
            "asv_len": len(asv_seqs[i]),
            "midas_taxonomy": midas_tax,
            "rel_ab": rel_ab,
            "best_align_score": best.get("align_score"),
            "best_identity": best.get("identity"),
            "n_edlib_candidates": min(args.k1, n_ref),
            "top20": top,
        }
        summary_rows.append({
            "asv": asv, "asv_len": len(asv_seqs[i]), "rel_ab": rel_ab,
            "midas_taxonomy": midas_tax,
            "best_organism": best.get("organism"), "best_genome_id": best.get("genome_id"),
            "best_align_score": best.get("align_score"), "best_identity": best.get("identity"),
            "best_edlib_distance": best.get("edlib_distance"),
            "best_genus": best.get("lineage", {}).get("Genus"),
            "best_family": best.get("lineage", {}).get("Family"),
        })

    # --------------------------------- write outputs ---------------------------------
    out_json = os.path.join(args.outdir, "asv_top20_alignment_hits.json")
    with open(out_json, "w") as fh:
        json.dump(mapping, fh)
    out_csv = os.path.join(args.outdir, "asv_alignment_summary.csv")
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary_rows[0].keys()))
        w.writeheader(); w.writerows(summary_rows)

    wall = time.perf_counter() - wall0
    stats = {
        "n_asvs": n_asv, "n_unique_refs": n_ref, "min_ref_len": args.min_ref_len,
        "k1_edlib_candidates": args.k1, "k2_final_hits": args.k2, "workers": args.workers,
        "edlib_mode": "HW (infix) adaptive-k", "edlib_task": "distance",
        "biopython": {"mode": "local", "match": MATCH, "mismatch": MISMATCH,
                      "gap_open": GAP_OPEN, "gap_extend": GAP_EXTEND},
        "cpu_seconds_edlib_stage1": round(cpu_edlib, 1),
        "cpu_seconds_biopython_stage2": round(cpu_bio, 1),
        "wall_seconds_scan": round(scan_wall, 1),
        "wall_seconds_total": round(wall, 1),
        "parallel_speedup_scan": round((cpu_edlib + cpu_bio) / max(scan_wall, 1e-9), 1),
    }
    with open(os.path.join(args.outdir, "run_stats.json"), "w") as fh:
        json.dump(stats, fh, indent=2)
    print(f"[done] wrote {out_json}\n        {out_csv}\n        run_stats.json", flush=True)
    print(f"[time] scan wall={scan_wall:.0f}s  total wall={wall:.0f}s  "
          f"cpu edlib={cpu_edlib:,.0f}s bio={cpu_bio:,.0f}s", flush=True)


if __name__ == "__main__":
    main()
