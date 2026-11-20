import json
d = json.load(open('/home/freiburger/Documents/prFBA/region_validation/data/taxacc_region_summary.json'))
printed = {'FullLength16S': ('0.987','0.858','0.730','+0.128'), 'V1-V2': ('0.974','0.830','0.720','+0.110'), 'V1-V3': ('0.977','0.840','0.717','+0.123'),
 'V3-V4': ('0.969','0.807','0.639','+0.168'), 'V4': ('0.964','0.790','0.599','+0.191'), 'V4-V5': ('0.968','0.815','0.626','+0.189'),
 'V4-V5_944R': ('0.968','0.820','0.631','+0.189'), 'V7-V9': ('0.971','0.819','0.653','+0.193')}
for reg, R in d['regions'].items():
    a = R['anchor']; pr = lambda m, rk: a[m]['all']['all']['per_rank'][rk]
    fi, fe, ge = pr('include', 'Family'), pr('exclude_0.987', 'Family'), pr('exclude_0.987', 'Genus')
    vals = (fi['n_correct']/fi['n_called'], fe['n_correct']/fe['n_called'], ge['n_correct']/ge['n_called'])
    cur = ('%.3f' % vals[0], '%.3f' % vals[1], '%.3f' % vals[2], '%+.3f' % (vals[1]-vals[2]))
    flag = ['' if c == p else '<-- %s vs printed %s' % (c, p) for c, p in zip(cur, printed[reg])]
    print('%-14s %.5f %.5f %.5f %+.5f %s' % (reg, *vals, vals[1]-vals[2], ' '.join(f for f in flag if f)))
