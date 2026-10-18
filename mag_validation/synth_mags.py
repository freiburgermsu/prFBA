#!/usr/bin/env python
"""Synthesize MAGs from complete BV-BRC assemblies for the 10k MAG validation.

Two-step degradation, not one:

  1. FRAGMENT each deposited contig into synthetic contigs with a target mean
     length. This is required, not cosmetic: 19.8% of panel genomes are closed
     (1 contig) and 40.4% have <10 contigs, so dropping *deposited* contigs
     cannot reach a 10%-granular completeness grid, and would make degradation
     severity depend on each reference's assembly quality -- an uncontrolled
     confound across a panel spanning 1 to 5,265 contigs. Fragmenting to a fixed
     contiguity first removes it and is closer to what binning actually sees.
     Deposited contig boundaries are never joined (they are not adjacent in
     reality); each is cut independently.

  2. DROP whole fragments until the retained fraction meets the target
     completeness. Bases are never mutated -- skani compares k-mer sketches, and
     incompleteness, not divergence, is what a MAG actually suffers.

Contamination splices in fragments from a donor genome of a different phylum.
Every emitted contig header is prefixed with the MAG id, because skani keys
Query_name on the FIRST CONTIG HEADER and silently collapses queries whose first
header collides (verified: 8 variants sharing '>src_0' returned 5 of 8 queries).

Determinism: the RNG is seeded from a stable CRC32 of the genome id (Python's
str hash is salted per process), so a given (genome, tier) is byte-identical
across runs and machines.
"""
import argparse, glob, json, os, random, sys, zlib

MEAN_FRAG = 50_000          # target mean synthetic contig length (short-read-MAG scale)
MIN_FRAG  = 2_500           # MetaBAT2's default minimum contig length
SEED      = 1729


def read_fasta(path):
    recs, h, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            if line[0] == '>':
                if h is not None:
                    recs.append(''.join(buf))
                h, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if h is not None:
        recs.append(''.join(buf))
    return recs


def fragment(contigs, rng, mean_len=MEAN_FRAG, min_len=MIN_FRAG):
    """Cut each deposited contig into synthetic contigs of ~exponential length."""
    out = []
    for seq in contigs:
        i, n = 0, len(seq)
        while i < n:
            L = max(min_len, int(rng.expovariate(1.0 / mean_len)))
            piece = seq[i:i + L]
            # absorb a runt tail rather than emit a sub-min_len contig
            if n - (i + len(piece)) < min_len:
                piece = seq[i:]
            if len(piece) >= min_len or not out:
                out.append(piece)
            i += len(piece)
    return out


def take_to_fraction(frags, frac, rng, allow_empty=False):
    """Shuffle, then accumulate whole fragments up to frac of total bases.

    BUGFIX (audit): the original guard was ``if acc + len(f) > target and kept``,
    which could never skip the FIRST fragment in shuffled order. For the source
    genome that is harmless, but for the contamination graft -- where the target is
    only ~10% of the retained source -- it floored the dose at one whole fragment
    however small the target, so realized contamination reached 0.97 on small or
    heavily-degraded MAGs and the over-dose rate rose from 1.7% to 13.5% as
    completeness fell. Pass allow_empty=True for the graft: the first fragment is
    now skippable, and a fragment that would overshoot is truncated rather than
    taken whole, so the realized dose tracks the target.
    """
    total = sum(map(len, frags))
    target = frac * total
    if target <= 0:
        return [], 0, total
    order = list(frags)
    rng.shuffle(order)
    kept, acc = [], 0
    for f in order:
        if acc >= target:
            break
        room = target - acc
        if len(f) <= room:
            kept.append(f); acc += len(f)
        elif allow_empty:
            # truncate to fit rather than overshoot by a whole fragment
            if room >= MIN_FRAG:
                kept.append(f[:int(room)]); acc += int(room)
            break
        elif not kept:
            kept.append(f); acc += len(f)   # source side: never return an empty MAG
    return kept, acc, total


def write(path, tag, src_frags, contam_frags):
    with open(path + ".tmp", 'w') as fh:
        for i, s in enumerate(src_frags):
            fh.write(">%s|s%d\n" % (tag, i))
            for j in range(0, len(s), 80):
                fh.write(s[j:j + 80] + "\n")
        for i, s in enumerate(contam_frags):
            fh.write(">%s|x%d\n" % (tag, i))
            for j in range(0, len(s), 80):
                fh.write(s[j:j + 80] + "\n")
    os.replace(path + ".tmp", path)


def synth_one(gid, src_path, donor_path, grid, outdir):
    """Emit every (completeness, contamination) variant for one source genome."""
    base_rng = random.Random(SEED + zlib.crc32(gid.encode()))
    src_contigs = read_fasta(src_path)
    deposited_bp = sum(map(len, src_contigs))      # before the MIN_FRAG filter
    src_frags = fragment(src_contigs, base_rng)
    donor_frags = None
    rows = []
    for comp, contam in grid:
        rng = random.Random(SEED + zlib.crc32(("%s|%g|%g" % (gid, comp, contam)).encode()))
        kept, acc, total = take_to_fraction(src_frags, comp, rng)
        cfr = []
        if contam > 0:
            if donor_frags is None:
                donor_frags = fragment(read_fasta(donor_path), random.Random(SEED))
            want = contam * acc
            dtot = sum(map(len, donor_frags)) or 1
            cfr, cacc, _ = take_to_fraction(donor_frags, min(1.0, want / dtot), rng,
                                            allow_empty=True)
        tag = "%s_c%d_x%d" % (gid, round(comp * 100), round(contam * 100))
        write(os.path.join(outdir, tag + ".fna"), tag, kept, cfr)
        rows.append(dict(mag=tag, genome_id=gid, completeness_target=comp,
                         contamination_target=contam,
                         src_bp=acc, contam_bp=sum(map(len, cfr)),
                         total_src_bp=total, deposited_bp=deposited_bp,
                         n_contigs=len(kept) + len(cfr),
                         # against the DEPOSITED genome, not the post-MIN_FRAG pool
                         completeness_actual=round(acc / deposited_bp, 4) if deposited_bp else 0,
                         completeness_vs_pool=round(acc / total, 4) if total else 0,
                         contamination_actual=round(sum(map(len, cfr)) / (acc + sum(map(len, cfr))), 4)
                         if (acc + sum(map(len, cfr))) else 0))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--donors", required=True, help="json: genome_id -> donor genome_id")
    ap.add_argument("--donor-contigs", required=True,
                    help="dir holding donor FASTAs (kept OUT of the reference sketch so "
                         "splicing cannot create an exact ANI tie against the donor)")
    ap.add_argument("--grid", required=True, help="json list of [completeness, contamination]")
    ap.add_argument("--ids", required=True, help="file with one genome_id per line")
    ap.add_argument("--manifest", required=True)
    a = ap.parse_args()

    grid = [tuple(x) for x in json.loads(a.grid)]
    donors = json.load(open(a.donors))
    ids = [l.strip() for l in open(a.ids) if l.strip()]
    os.makedirs(a.outdir, exist_ok=True)
    rows, missing = [], []
    for gid in ids:
        sp = os.path.join(a.contigs, gid + ".fna")
        if not os.path.exists(sp):
            missing.append(gid); continue
        dp = os.path.join(a.donor_contigs, (donors.get(gid) or "") + ".fna")
        if not os.path.exists(dp):
            dp = None
            if any(c > 0 for _, c in grid):
                missing.append(gid + ":donor"); continue
        rows += synth_one(gid, sp, dp, grid, a.outdir)
    json.dump(dict(n_mags=len(rows), n_missing=len(missing), missing=missing[:50],
                   mean_frag=MEAN_FRAG, min_frag=MIN_FRAG, seed=SEED, rows=rows),
              open(a.manifest, 'w'))
    print("synthesized %d MAGs from %d genomes (%d skipped)" % (len(rows), len(ids) - len(missing), len(missing)))


if __name__ == "__main__":
    main()
