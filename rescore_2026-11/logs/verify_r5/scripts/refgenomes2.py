import json, re, collections
M = '/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/'
seq = json.load(open(M + '16S_md5_seq.json'))
ids = json.load(open(M + '16S_md5_ID.json'))
T = re.compile(r"\[(.+?)\s*\|\s*([0-9]+\.[0-9]+)\]\s*$")
def g_tail(h):
    m = T.search(h); return m.group(2) if m else None
def g_fig(h): return h.split('|')[1].rsplit('.rna.', 1)[0]
a = {g_tail(h) for h in ids.values()}
b = {g_tail(ids[k]) for k, s in seq.items() if len(s) >= 200}
fullen = {g_tail(ids[k]) for k, s in seq.items() if len(s) >= 1400}
dis = sum(1 for h in ids.values() if g_tail(h) != g_fig(h))
print('pipeline parser: all md5 %d, >=200 bp %d, >=1400 bp %d; None in set: %s; tail!=fig prefix: %d' % (len(a), len(b), len(fullen), None in a, dis))
keys = set(open('/home/freiburger/Documents/prFBA_rescore_2026-11/stamps/R1.6.cache_keys.txt').read().split())
print('cache keys %d; >=200 bp reference genomes with a cache key %d; cache keys outside the reference %d' % (len(keys), len(b & keys), len(keys - a)))
