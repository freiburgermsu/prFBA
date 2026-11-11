#!/usr/bin/env python
"""Which synthetic MAGs still carry a 16S gene? -- PICRUSt2's attainable MAG coverage.

A MAG reaches PICRUSt2 only through a 16S gene on one of its contigs. The MAG arm's
synthetic MAGs were deleted after the skani search, but ``synth_mags.py`` is
deterministic, so this replays it on COORDINATES instead of sequence: the
synthesizer's own ``fragment()`` and ``take_to_fraction()`` are called unchanged on
``Span`` objects (``len`` + slicing, all those functions use), with the synthesizer's
per-genome / per-tier seeding copied verbatim. Each kept span is then intersected
with the 16S loci from ``mag_16s_loci.py``.

Replay fidelity is checked against every manifest the MAG arm wrote: ``src_bp`` must
match for every MAG, and ``n_contigs`` for every uncontaminated one. (The contaminated
manifests predate the ``allow_empty`` graft fix, so their contamination side is not
comparable; contaminated tiers here use the CURRENT synthesizer and are new data.)

Per MAG, ``piece_max`` is the longest contiguous retained stretch of any 16S locus
(a fragment boundary can split a gene across two contigs). Status:
  full     -- a near-full-length locus (>= FULL_MIN bp) is wholly inside one kept contig
  partial  -- otherwise, the longest retained piece is >= PARTIAL_MIN bp (amplicon-sized;
              PICRUSt2 places 250 bp V4 ASVs routinely)
  none     -- nothing PICRUSt2 could take (short 16S fragments at draft-contig ends,
              wholly retained or not, do not count)

    python mag_16s_retention.py            # writes data/mag_arm/retention.csv.gz + summary
"""
import argparse
import collections
import csv
import glob
import gzip
import json
import os
import random
import sys
import zlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mag_common as M  # noqa: E402

sys.path.insert(0, M.MV)
import synth_mags as S  # noqa: E402

PARTIAL_MIN = 250
FULL_MIN = 1200
CLEAN = [1.0, 0.9, 0.7, 0.5, 0.3, 0.25, 0.2, 0.15, 0.1, 0.07, 0.05, 0.03, 0.02, 0.01]
DIRTY = [1.0, 0.9, 0.7, 0.5, 0.3]


class Span:
    """A contig interval that slices like the sequence string it stands in for."""
    __slots__ = ("ci", "a", "b")

    def __init__(self, ci, a, b):
        self.ci, self.a, self.b = ci, a, b

    def __len__(self):
        return self.b - self.a

    def __getitem__(self, sl):
        start, stop, _ = sl.indices(self.b - self.a)
        return Span(self.ci, self.a + start, self.a + max(start, stop))


def retained(loci, spans):
    """(n near-full-length loci wholly kept, longest contiguous kept piece of any locus,
    indices of loci >= PARTIAL_MIN wholly kept, index of the locus giving that piece)."""
    by_ci = collections.defaultdict(list)
    for s in spans:
        by_ci[s.ci].append(s)
    n_full, piece, whole, piece_i = 0, 0, [], None
    for i, L in enumerate(loci):
        best, full = 0, False
        for s in by_ci.get(L["ci"], ()):
            ov = min(L["hi"], s.b) - max(L["lo"], s.a)
            if ov > best:
                best = ov
            if s.a <= L["lo"] and L["hi"] <= s.b:
                full = True
        n_full += full and L["len"] >= FULL_MIN
        if full and L["len"] >= PARTIAL_MIN:
            whole.append(i)
        if best > piece:
            piece, piece_i = best, i
    return n_full, piece, whole, piece_i


def status(n_full, piece):
    return "full" if n_full else ("partial" if piece >= PARTIAL_MIN else "none")


def replay(gid, rec, donor_rec, grid):
    """Mirror synth_mags.synth_one for one genome over the grid."""
    base_rng = random.Random(S.SEED + zlib.crc32(gid.encode()))
    src = S.fragment([Span(i, 0, n) for i, n in enumerate(rec["contig_lens"])], base_rng)
    donor_frags, rows = None, []
    pool_full, pool_piece, _, _ = retained(rec["loci"], src)    # what fragmenting leaves
    for comp, contam in grid:
        rng = random.Random(S.SEED + zlib.crc32(("%s|%g|%g" % (gid, comp, contam)).encode()))
        kept, acc, total = S.take_to_fraction(src, comp, rng)
        cfr = []
        if contam > 0:
            if donor_frags is None:
                donor_frags = S.fragment([Span(i, 0, n) for i, n in
                                          enumerate(donor_rec["contig_lens"])], random.Random(S.SEED))
            want = contam * acc
            dtot = sum(map(len, donor_frags)) or 1
            cfr, _, _ = S.take_to_fraction(donor_frags, min(1.0, want / dtot), rng, allow_empty=True)
        n_full, piece, whole, piece_i = retained(rec["loci"], kept)
        row = dict(mag="%s_c%d_x%d" % (gid, round(comp * 100), round(contam * 100)),
                   genome_id=gid, completeness=comp, contamination=contam,
                   src_bp=acc, n_contigs=len(kept) + len(cfr),
                   genome_loci=len(rec["loci"]),
                   genome_loci_ge1200=sum(1 for L in rec["loci"] if L["len"] >= 1200),
                   pool_loci_full=pool_full, pool_piece_max=pool_piece,
                   n_full=n_full, piece_max=piece, status=status(n_full, piece),
                   kept_loci=";".join(map(str, whole)),
                   piece_locus="" if piece_i is None or piece < PARTIAL_MIN else piece_i)
        if contam > 0:
            d_full, d_piece, _, _ = retained(donor_rec["loci"], cfr)
            row.update(contam_bp=sum(map(len, cfr)), donor_n_full=d_full,
                       donor_piece_max=d_piece, donor_status=status(d_full, d_piece))
        rows.append(row)
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--loci", default=os.path.join(M.DATA, "loci.jsonl"))
    ap.add_argument("--out", default=os.path.join(M.DATA, "retention.csv.gz"))
    ap.add_argument("--summary", default=os.path.join(M.P2C, "data", "mag_retention_summary.json"))
    a = ap.parse_args()

    loci = {}
    for line in open(a.loci):
        r = json.loads(line)
        loci[r["id"]] = r
    donors = json.load(open(os.path.join(M.MV, "donor_assignment.json")))
    grid = [(c, 0) for c in CLEAN] + [(c, 0.1) for c in DIRTY]
    rows = []
    for gid in M.panel_ids():
        rec = loci[gid]
        d = donors.get(gid)
        dg = [g for g in grid if g[1] == 0 or (d and d in loci)]
        rows += replay(gid, rec, loci.get(d), dg)
    by_mag = {r["mag"]: r for r in rows}

    # ---- fidelity: the replay must reproduce what the MAG arm actually searched ----
    check = collections.Counter()
    for p in sorted(glob.glob(os.path.join(M.MV, "manifests", "*.json")) +
                    glob.glob(os.path.join(M.MV, "lowc_manifests", "*.json"))):
        for m in json.load(open(p))["rows"]:
            r = by_mag.get(m["mag"])
            if r is None:
                check["manifest_mag_not_replayed"] += 1
                continue
            check["src_bp_match" if r["src_bp"] == m["src_bp"] else "src_bp_MISMATCH"] += 1
            if m["contamination_target"] == 0:
                check["n_contigs_match" if r["n_contigs"] == m["n_contigs"]
                      else "n_contigs_MISMATCH"] += 1
    print("[retention] replay vs MAG-arm manifests:", dict(check), flush=True)

    fields = list(dict.fromkeys(k for r in rows for k in r))
    with gzip.open(a.out, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)

    # ---- summary ----
    summ = {"partial_min_bp": PARTIAL_MIN, "full_min_bp": FULL_MIN, "replay_check": dict(check),
            "clean": [], "contaminated": []}
    for comp in CLEAN:
        rs = [r for r in rows if r["completeness"] == comp and r["contamination"] == 0]
        n = len(rs)
        c = collections.Counter(r["status"] for r in rs)
        summ["clean"].append(dict(completeness=comp, n=n,
                                  full=round(c["full"] / n, 4), partial=round(c["partial"] / n, 4),
                                  none=round(c["none"] / n, 4),
                                  usable=round((c["full"] + c["partial"]) / n, 4)))
    for comp in DIRTY:
        rs = [r for r in rows if r["completeness"] == comp and r["contamination"] == 0.1]
        n = len(rs)
        c = collections.Counter((r["status"] != "none", r["donor_status"] != "none") for r in rs)
        summ["contaminated"].append(dict(
            completeness=comp, n=n,
            source_only=round(c[(True, False)] / n, 4), donor_only=round(c[(False, True)] / n, 4),
            both=round(c[(True, True)] / n, 4), neither=round(c[(False, False)] / n, 4)))
    c100 = [r for r in rows if r["completeness"] == 1.0 and r["contamination"] == 0]
    why = collections.Counter(
        "genome_has_no_16s" if r["genome_loci"] == 0 else
        "only_sub_250bp_fragments_in_genome" if max((L["len"] for L in loci[r["genome_id"]]["loci"]), default=0) < PARTIAL_MIN else
        "16s_only_on_contigs_dropped_by_min_frag" if r["pool_piece_max"] < PARTIAL_MIN else
        "kept" if r["status"] != "none" else "other"
        for r in c100)
    summ["c100_reasons"] = dict(why)
    json.dump(summ, open(a.summary, "w"), indent=1)
    for s in summ["clean"]:
        print("  c%-5g usable %.3f (full %.3f, partial %.3f)" % (s["completeness"] * 100, s["usable"],
                                                                 s["full"], s["partial"]))
    for s in summ["contaminated"]:
        print("  c%-5g x10  src-only %.3f donor-only %.3f both %.3f neither %.3f" % (
            s["completeness"] * 100, s["source_only"], s["donor_only"], s["both"], s["neither"]))
    print("  c100 reasons:", dict(why))


if __name__ == "__main__":
    sys.exit(main())
