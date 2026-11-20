import json
W = '/home/freiburger/Documents/prFBA_rescore_2026-11/pre_repair/'
R = '/home/freiburger/Documents/prFBA/'
f = 'region_validation_d100/picrust2_comparison/data/comparison_summary.json'
METHODS = ('picrust2', 'prfba_exclude')
def cf1(cell): 
    c = cell.get('conditional') or {}
    return (c.get('f1') or {}).get('value'), cell.get('n_called')
def run(path, tag):
    d = json.load(open(path))
    print('==', tag)
    br = d['by_region']
    for m in METHODS:
        e, _ = cf1(br['ec']['all']['all'][m]); r, _ = cf1(br['rxn']['all']['all'][m])
        print('  overall %-14s ec %.4f rxn %.4f diff %+.4f' % (m, e, r, r - e))
    # region x domain
    worst = {}
    for reg in br['ec']:
        for dom in br['ec'][reg]:
            for m in METHODS:
                ce, ne = cf1(br['ec'][reg][dom][m]); cr, nr = cf1(br['rxn'][reg][dom][m])
                if ce is None or cr is None: continue
                dd = cr - ce
                key = (dom, m)
                if dom != 'all' or reg != 'all':
                    if abs(dd) > 0.01:
                        print('   >0.01 region %-14s %-8s %-14s n_called %5s diff %+.4f' % (reg, dom, m, ne, dd))
                    worst[key] = max(worst.get(key, 0), abs(dd))
    print('  max |rxn-ec| by (domain, method):', {k: round(v, 4) for k, v in worst.items()})
    bs = d['by_stratum']
    for s in bs['ec']:
        for b in bs['ec'][s]:
            for m in METHODS:
                ce, ne = cf1(bs['ec'][s][b][m]); cr, nr = cf1(bs['rxn'][s][b][m])
                if ce is None or cr is None: continue
                dd = cr - ce
                flag = ' <--' if abs(dd) > 0.01 else ''
                print('   stratum %-20s %-12s %-14s n_called %6s diff %+.4f%s' % (s, b, m, ne, dd, flag))
run(W + f, 'PRE (pre_repair)')
run(R + f, 'POST (live)')
