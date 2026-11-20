import json, subprocess
R = '/home/freiburger/Documents/prFBA/'
W = '/home/freiburger/Documents/prFBA_rescore_2026-11/'
F = {'genus': 'region_validation/shared_reference/data/summary_genus.json',
     'family': 'region_validation/shared_reference/data_family/summary_family.json'}
def load_pre(rel):
    try:
        return json.loads(subprocess.check_output(['git', '-C', R, 'show', 'rescore-2026-11-start:' + rel], stderr=subprocess.DEVNULL))
    except subprocess.CalledProcessError:
        import os
        for base in (W + 'pre_repair/', W + 'start/'):
            if os.path.exists(base + rel): return json.load(open(base + rel))
        return None
for state in ('pre', 'post'):
    for rank, rel in F.items():
        d = load_pre(rel) if state == 'pre' else json.load(open(R + rel))
        if d is None: print(state, rank, 'no copy'); continue
        nc = d['nested_cv']; ms = []
        for k in sorted(x for x in nc if x.startswith('fold_')):
            f = nc[k]['inner_f1']; w = nc[k]['inner_winner']
            srt = sorted(f.items(), key=lambda kv: -kv[1])
            ms.append((k, w, srt[0][0], round(srt[0][1] - srt[1][1], 4), srt[1][0], round(f['legacy'] - f['exact_tie'], 4)))
        print(state, rank, 'winner_instability', nc.get('winner_instability'))
        for m in ms: print('   %s winner %s (top %s) margin over runner-up %.4f (%s); legacy-exact_tie %.4f' % m)
