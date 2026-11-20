import json, collections, sys, statistics
sys.path.insert(0, '/home/freiburger/Documents/prFBA/manuscript_numbers')
import derived_numbers as DN
R = '/home/freiburger/Documents/prFBA/'
res = DN.Resolver('post', R.rstrip('/'), '/home/freiburger/Documents/prFBA_rescore_2026-11')
meta, hits = DN._load_c100_hits(res)
gid_by_mag = {m: r['genome_id'] for m, r in meta.items()}
hits_by_gid = {gid_by_mag[m]: hits.get(m, []) for m in meta}
truth = {g: v['truth_lineage'] for g, v in DN.load_truth(R + 'region_validation/data/benchmark_truth.json').items()}
mr = DN.marker_reach(truth, hits_by_gid)
print('W6 reproduce:', mr['n_congener'], mr['n_congener_with_mag'], round(mr['median_best_af'], 3), round(100*mr['below_gate'], 2), round(100*mr['nohit'], 2))
sel = {g['genome_id']: g for g in json.load(open(R + 'region_validation/data/selection_10k.json'))['genomes']}
sk = set(open(R + 'mag_validation/all_ids.txt').read().split())
sl2 = lambda g: (None if (not sel[g].get('genus') or sel[g]['genus'].startswith('unclassified')) else sel[g]['genus'])
gl = lambda g: (truth[g].get('Genus') if g in truth else sl2(g))
cnt = collections.Counter(gl(g) for g in sk if gl(g))
pop = [g for g in truth if g in sk and g in hits_by_gid and gl(g) and cnt[gl(g)] >= 2]
best, nohit = [], 0
for g in pop:
    hh = DN.exclude_hits(hits_by_gid[g], g)
    if not hh: nohit += 1; continue
    best.append(max(h['af_query'] for h in hh))
print('strict (genus-mate among the 9,997 sketched):', len(pop), 'median AF %.3f' % statistics.median(best),
      'below gate %.2f%%' % (100*sum(a < 15 for a in best)/len(best)), 'no hit %.2f%%' % (100*nohit/len(pop)))
