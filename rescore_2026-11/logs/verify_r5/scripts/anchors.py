import json, sys, gzip, csv, collections
sys.path.insert(0, '/home/freiburger/Documents/prFBA/region_validation/scripts')
import score_taxacc as ST
from _common import lineage_for
W = '/home/freiburger/Documents/prFBA_rescore_2026-11/pre_repair/'
R = '/home/freiburger/Documents/prFBA/'
truth = json.load(open(R + 'region_validation/data/benchmark_truth.json'))
f = 'region_validation/picrust2_comparison/data/prfba_selection.json'
cf = 'region_validation/picrust2_comparison/data/comparison_cells.csv.gz'
def ec_called(path):
    # amplicons where prFBA exclude answered and the cell is EC-scoreable
    out = set(); hdr = None
    with gzip.open(path, 'rt') as fh:
        rd = csv.DictReader(fh)
        hdr = rd.fieldnames
        for r in rd:
            out.add(r)
    return hdr
def run(path, cells_path, tag):
    sel = json.load(open(path))
    # cell table: which amplicons are EC-scoreable
    with gzip.open(cells_path, 'rt') as fh:
        rd = csv.DictReader(fh); cols = rd.fieldnames; rows = list(rd)
    space_col = 'space' if 'space' in cols else None
    method_col = 'method' if 'method' in cols else None
    ec_amp = set()
    for r in rows:
        if (space_col is None or r[space_col] == 'ec') and (method_col is None or r[method_col] == 'prfba_exclude'):
            if r.get('called', r.get('n_called', '1')) in ('1', 'True', 'true'):
                ec_amp.add(r.get('amplicon') or r.get('key'))
    res = {}
    for variant, keep in (('all answered', lambda k: True), ('EC-scoreable answered', lambda k: k in ec_amp)):
        acc = {}
        for rank in ('Family', 'Genus'):
            n = c = 0
            for k, v in sel.items():
                if k == '_meta' or not isinstance(v, dict): continue
                e = v.get('exclude_0.987') or {}
                if not (e.get('n_selected') or 0) > 0 or not keep(k): continue
                t = truth.get(k)
                if t is None: continue
                a = e['reps'][0]
                pl = lineage_for(int(a.split('.')[0]))
                rc = ST.rank_correct(t['truth_lineage'].get(rank), pl.get(rank))
                if rc is None: continue
                n += 1; c += bool(rc)
            acc[rank] = (c / n if n else None, n)
        res[variant] = acc
    print(tag, 'cells cols:', cols[:14], 'ec_amp', len(ec_amp))
    for v, acc in res.items():
        print('   %-24s family %.6f (n %d)  genus %.6f (n %d)' % (v, acc['Family'][0], acc['Family'][1], acc['Genus'][0], acc['Genus'][1]))
run(W + f, W + cf, 'PRE')
run(R + f, R + cf, 'POST')
