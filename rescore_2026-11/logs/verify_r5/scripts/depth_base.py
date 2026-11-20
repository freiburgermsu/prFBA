import json, subprocess, os
R = '/home/freiburger/Documents/prFBA/'
W = '/home/freiburger/Documents/prFBA_rescore_2026-11/'
def load(rel, state):
    if state == 'post': return json.load(open(R + rel))
    for base in (W + 'pre_repair/',):
        if os.path.exists(base + rel): return json.load(open(base + rel))
    try:
        return json.loads(subprocess.check_output(['git', '-C', R, 'show', 'rescore-2026-11-start:' + rel], stderr=subprocess.DEVNULL))
    except subprocess.CalledProcessError:
        return None
for tree in ('region_validation', 'region_validation_d100_legacy', 'region_validation_d100'):
    rel = tree + '/data/taxacc_region_summary.json'
    for state in ('pre', 'post'):
        d = load(rel, state)
        if d is None: print(tree, state, 'missing'); continue
        n = nc = ncor = 0
        for reg, Rg in d['regions'].items():
            g = Rg['anchor']['exclude_0.987']['all']['all']['per_rank']['Genus']
            n += g['n']; nc += g['n_called']; ncor += g['n_correct']
        print('%-32s %-4s select_mode %-8s genus pooled: coverage %.4f  pop acc %.4f  cond acc %.4f  (n %d)' % (
            tree, state, d['_meta'].get('select_mode'), nc / n, ncor / n, ncor / nc, n))
