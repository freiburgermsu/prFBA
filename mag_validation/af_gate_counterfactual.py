#!/usr/bin/env python
"""The 2% aligned-fraction-gate counterfactual behind the manuscript's AF-gate paragraph
(sec:magval) and the deck's abstention slide: re-select every 100%-complete, uncontaminated
synthetic MAG from its stored skani hits under the shipped 15% gate and a 2% gate (self-
exclusion at ANI 95), score each answer's conditional gene-capture F1, and count the newly
admitted MAGs that beat the full-length 16S level. Writes out/af_gate_counterfactual.json.
Moved from a session scratch directory on 2026-10-27; first run 2026-10-18.

The full-length 16S levels are read at run time from the paired MAG-vs-16S comparisons
(key ``exclude|c1_x0|FullLength16S|f1``, conditional ``amp``): out/paired_comparison_d100_legacy.json
(legacy reducer, depth 100), out/paired_comparison_d100.json (exact-tie union, depth 100) and
out/paired_comparison.json (depth 20), so run compare_arms.py first. Until 2026-11 they were
hard-coded (0.6912 / 0.71 / 0.7092) and went stale when those comparisons were re-scored.

    python mag_validation/af_gate_counterfactual.py
"""
import sys, os, csv, json, collections, pickle
sys.path.insert(0, '/home/freiburger/Documents/prFBA/mag_validation')
import score_mag as SM
from mag_skani_hits import ANI_KNOBS
import select_references as SR

MV = os.path.dirname(os.path.abspath(__file__))
LEVEL_KEY = "exclude|c1_x0|FullLength16S|f1"
LEVEL_SOURCES = (("legacy reducer, depth 100", "paired_comparison_d100_legacy.json"),
                 ("exact-tie union, depth 100", "paired_comparison_d100.json"),
                 ("depth 20", "paired_comparison.json"))


def fl16s_levels(out_dir=os.path.join(MV, "out")):
    """[(label, level, file)]: the conditional full-length 16S gene-capture F1 of the paired,
    self-excluded, 100%-complete uncontaminated comparison in each paired_comparison file."""
    return [(label, json.load(open(os.path.join(out_dir, name)))[LEVEL_KEY]["conditional"]["amp"], name)
            for label, name in LEVEL_SOURCES]


def main():
    levels = fl16s_levels()
    print('16S levels', levels)
    truth, pgf = SM.load_truth(), SM.load_pgfams()
    meta = {}
    for i in range(20):
        for r in json.load(open('/home/freiburger/Documents/prFBA/mag_validation/manifests/b%03d.json' % i))['rows']:
            if r['completeness_target'] == 1.0 and r['contamination_target'] == 0:
                meta[r['mag']] = r
    hits = collections.defaultdict(list)
    for i in range(20):
        with open('/home/freiburger/Documents/prFBA/mag_validation/hits/b%03d.tsv' % i) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                q = os.path.basename(row['Query_file'])[:-4]
                if q not in meta: continue
                hits[q].append(dict(ref=os.path.basename(row['Ref_file'])[:-4], ani=float(row['ANI']),
                                    af_query=float(row['Align_fraction_query']), af_ref=float(row['Align_fraction_ref'])))
    def run(gate):
        ANI_KNOBS['cov_min'] = gate; SM.AF_GATE = gate * 100.0
        res = {}
        for mag, m in sorted(meta.items()):
            gid = m['genome_id']; tg = pgf.get(gid) or set()
            hh = SM.exclude_hits(hits.get(mag, []), gid, SM.src_taxon(gid))
            rec, b = SM.build_record(mag, hh, m['src_bp'] + m['contam_bp'])
            sel = SR.select_representatives(rec, knobs=ANI_KNOBS, gene_provider=lambda g: set(pgf.get(g) or ()))
            picks = [s['genome_id'] for s in (sel.get('selected') or [])]
            if not tg: res[mag] = ('truth_empty', None); continue
            if not picks: res[mag] = ('abstained', None); continue
            pg = set()
            for g in picks: pg |= (pgf.get(g) or frozenset())
            if not pg: res[mag] = ('pred_empty', None); continue
            tp = len(tg & pg); rc = tp/len(tg); pr = tp/len(pg)
            res[mag] = ('ok', 2*rc*pr/(rc+pr) if rc+pr else 0.0)
        return res
    a = run(0.15); b = run(0.02)
    for n, r in (('15', a), ('2', b)): print(n, collections.Counter(v[0] for v in r.values()))
    newly = [m for m in meta if a[m][0] == 'abstained' and b[m][0] == 'ok']
    already = [m for m in meta if a[m][0] == 'ok']
    fn = [b[m][1] for m in newly]; fo = [a[m][1] for m in already]
    print('newly', len(newly), 'mean', sum(fn)/len(fn), 'already', len(already), 'mean', sum(fo)/len(fo))
    print('transitions', collections.Counter((a[m][0], b[m][0]) for m in meta))
    for label, thr, _src in levels:
        k = sum(1 for f in fn if f > thr); print(thr, label, k, '%.6f' % (k/len(newly)))
    # alternative: newly = not ok at 15 (abstained or pred_empty) and ok at 2
    alt = [m for m in meta if a[m][0] != 'ok' and b[m][0] == 'ok']
    fa = [b[m][1] for m in alt]
    print('alt newly', len(alt), sum(fa)/len(fa), {thr: sum(1 for f in fa if f > thr) for _l, thr, _s in levels})
    pool = fo + fn; print('pooled mean', sum(pool)/len(pool), 'coverage rel gain', len(newly)/len(already))
    out = {"_about": "Per-MAG outcome (status, conditional gene-capture F1) of the 100%-complete, "
                     "uncontaminated synthetic MAGs under the shipped 15% aligned-fraction gate and a 2% "
                     "gate, self-exclusion at ANI 95. Newly admitted = abstained at 15%, answered at 2%.",
           "summary": {"n_mags": len(a), "n_newly_admitted": len(newly),
                       "mean_f1_newly": round(sum(fn) / len(fn), 6),
                       "n_already_answered": len(already), "mean_f1_already": round(sum(fo) / len(fo), 6),
                       "newly_above_full_length_16S_level": {
                           "%s (%s)" % (thr, label): sum(1 for f in fn if f > thr) for label, thr, _s in levels},
                       "full_length_16S_levels": {
                           label: {"level": thr, "source": "mag_validation/out/" + src, "key": LEVEL_KEY,
                                   "scope": "conditional", "field": "amp"} for label, thr, src in levels}},
           "mags": {m: {"gate_0.15": list(a[m]), "gate_0.02": list(b[m])} for m in sorted(a)}}
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out",
                                     "af_gate_counterfactual.json"), "w"), indent=0, sort_keys=True)


if __name__ == "__main__":
    main()
