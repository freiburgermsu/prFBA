#!/usr/bin/env python
"""
insilico_pcr.py — extract a 16S sub-region from reference sequences by in-silico PCR.

The experimental ASVs in this study are V4–V5 amplicons (E. coli ~533–906, ~373 bp insert)
generated with the universal 515F / 926R primer pair (covers bacteria *and* archaea). To compare
amplicons against the BV-BRC reference space *like-with-like* (region-matched), we excise the same
V4–V5 insert from every reference 16S gene: locate the forward primer and the reverse-complement of
the reverse primer (IUPAC-degenerate, ≤3-mismatch fuzzy fallback for divergent taxa) and keep the
sequence between them.

Default primers (Parada/Quince, EMP V4–V5):
  515F  GTGYCAGCMGCCGCGGTAA
  926R  CCGYCAATTYMTTTRAGTTT      (its rev-comp is matched on the sense strand)

Input  : {header: full_16S_sequence} JSON (e.g. model_inputs/BV_BRC_16S.json)
Output : {header: region_insert} JSON for references where the region was found, + a stats JSON.
The output plugs directly into embed_16s.py (same header format → same metadata join).
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import regex

# IUPAC-expanded regexes. FWD = 515F on sense strand; REVc = reverse-complement of 926R on sense strand.
FWD_515F = r'GTG[CT]CAGC[AC]GCCGCGGTAA'
REVC_926R = r'AAACT[CT]AAA[GT][AG]AATTG[AG]CGG'


def make_extractor(fwd=FWD_515F, revc=REVC_926R, emax=3, lo=300, hi=450):
    fx, rx = regex.compile(fwd), regex.compile(revc)
    ff, rf = regex.compile(f'(?:{fwd}){{e<={emax}}}'), regex.compile(f'(?:{revc}){{e<={emax}}}')

    def extract(seq: str):
        f = fx.search(seq) or ff.search(seq)            # exact-degenerate first, fuzzy fallback
        if not f:
            return None
        r = rx.search(seq, f.end()) or rf.search(seq, f.end())
        if not r or not (lo <= r.start() - f.end() <= hi):
            return None
        return seq[f.end():r.start()]
    return extract


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fasta", type=Path, required=True, help="{header: full 16S sequence} JSON")
    ap.add_argument("--out", type=Path, required=True, help="{header: region insert} JSON to write")
    ap.add_argument("--emax", type=int, default=3, help="max primer mismatches (fuzzy fallback)")
    ap.add_argument("--lo", type=int, default=300); ap.add_argument("--hi", type=int, default=450)
    args = ap.parse_args()

    t0 = time.time()
    print(f"[pcr] loading {args.fasta}", flush=True)
    fasta = json.load(open(args.fasta))
    extract = make_extractor(emax=args.emax, lo=args.lo, hi=args.hi)
    out, n_in, n_full, n_full_hit = {}, 0, 0, 0
    lens = []
    for header, raw in fasta.items():
        n_in += 1
        seq = raw.strip().upper()
        is_full = 1400 <= len(seq) <= 1600
        n_full += is_full
        ins = extract(seq)
        if ins:
            out[header] = ins
            lens.append(len(ins))
            n_full_hit += is_full
        if n_in % 200000 == 0:
            print(f"[pcr] {n_in:,} scanned, {len(out):,} extracted ({time.time()-t0:.0f}s)", flush=True)

    json.dump(out, open(args.out, "w"))
    import numpy as np
    L = np.array(lens)
    stats = {
        "primers": {"515F": "GTGYCAGCMGCCGCGGTAA", "926R": "CCGYCAATTYMTTTRAGTTT"},
        "region": "V4-V5 (E. coli ~533-906)",
        "n_input": n_in, "n_extracted": len(out), "extract_rate": round(len(out) / max(n_in, 1), 4),
        "n_full_length_input": n_full, "n_full_length_extracted": n_full_hit,
        "full_length_extract_rate": round(n_full_hit / max(n_full, 1), 4),
        "insert_len": {"median": int(np.median(L)), "mean": round(float(L.mean()), 1),
                       "p5": int(np.percentile(L, 5)), "p95": int(np.percentile(L, 95)),
                       "min": int(L.min()), "max": int(L.max())},
        "wall_seconds": round(time.time() - t0, 1),
    }
    json.dump(stats, open(args.out.with_suffix(".stats.json"), "w"), indent=2)
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
