import json
M = '/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/'
seq = json.load(open(M + '16S_md5_seq.json'))
ids = json.load(open(M + '16S_md5_ID.json'))
def gid(h):
    # header: fig|<genome_id>.rna.N|...
    return h.split('|')[1].rsplit('.rna.', 1)[0]
all_g = {gid(h) for h in ids.values()}
ge200 = {gid(ids[k]) for k, s in seq.items() if len(s) >= 200}
print('distinct representative genome ids: all md5 %d; >=200 bp md5 %d' % (len(all_g), len(ge200)))
keys = set(open('/home/freiburger/Documents/prFBA_rescore_2026-11/stamps/R1.6.cache_keys.txt').read().split())
print('cache keys %d; reference genomes (>=200 bp) with a cache entry %d' % (len(keys), len(ge200 & keys)))
