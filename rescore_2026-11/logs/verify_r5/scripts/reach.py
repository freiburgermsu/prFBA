import json, collections, glob, sys
sys.path.insert(0, '/home/freiburger/Documents/prFBA/manuscript_numbers')
R = '/home/freiburger/Documents/prFBA/'
import derived_numbers as DN
truth = {g: v['truth_lineage'] for g, v in DN.load_truth(R + 'region_validation/data/benchmark_truth.json').items()}
sel = {g['genome_id']: g for g in json.load(open(R + 'region_validation/data/selection_10k.json'))['genomes']}
sk = set(open(R + 'mag_validation/all_ids.txt').read().split())
print('truth', len(truth), 'panel', len(sel), 'sketched', len(sk), 'truth&sketched', len(set(truth) & sk),
      'sketched-not-truth', len(sk - set(truth)), 'truth-not-sketched', sorted(set(truth) - sk))
# agreement of genus names between truth_lineage and selection_10k on the overlap
dis = [g for g in truth if (truth[g].get('Genus') or None) != (sel[g].get('genus') or None)]
print('genus name disagreements truth vs selection_10k:', len(dis), dis[:5], [(truth[g].get('Genus'), sel[g].get('genus')) for g in dis[:5]])
def n_cong(universe, pool, lin):
    cnt = collections.Counter(lin(g) for g in pool if lin(g))
    return [g for g in universe if lin(g) and cnt[lin(g)] >= 2]
tl = lambda g: truth[g].get('Genus') if g in truth else None
sl = lambda g: sel[g].get('genus')
a = n_cong(truth, truth, tl)
print('W6: truth organisms with a genus-mate in truth:', len(a), '; with MAG:', sum(g in sk for g in a))
c = n_cong([g for g in truth if g in sk], sk, sl)
print('truth organisms with a MAG and a genus-mate among the 9,997 sketched (selection_10k genus):', len(c))
c2 = n_cong(truth, sel, sl)
print('truth organisms with a genus-mate in the 10,000 panel (selection_10k genus):', len(c2))
print('--- with unclassified_* placeholders treated as no genus')
sl2 = lambda g: (None if (not sel[g].get('genus') or sel[g]['genus'].startswith('unclassified')) else sel[g]['genus'])
dis2 = [g for g in truth if (truth[g].get('Genus') or None) != sl2(g)]
print('remaining disagreements:', len(dis2), [(g, truth[g].get('Genus'), sl2(g)) for g in dis2[:8]])
# pool = sketched genomes; organism's own genus from truth, mates' genus from truth where available else selection_10k
gl = lambda g: tl(g) if g in truth else sl2(g)
c3 = n_cong([g for g in truth if g in sk], sk, gl)
print('truth organisms with a MAG and a genus-mate among the 9,997 sketched genomes:', len(c3))
extra = sorted(set(c3) - set(a)); lost = sorted(set(a) - set(c3))
print('  gained (mate only among the 254 sketched non-truth genomes):', len(extra), extra[:10])
print('  lost (no MAG, or only mate unsketched):', lost)
