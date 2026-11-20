import json, subprocess, sys
def table(d):
    out = {}
    for reg, R in d['regions'].items():
        a = R['anchor']
        fi = a['include']['all']['all']['per_rank']['Family']['correct']
        fe = a['exclude_0.987']['all']['all']['per_rank']['Family']['correct']
        ge = a['exclude_0.987']['all']['all']['per_rank']['Genus']['correct']
        out[reg] = (fi, fe, ge, fe - ge)
    return out
cur = json.load(open('/home/freiburger/Documents/prFBA/region_validation/data/taxacc_region_summary.json'))
old = json.loads(subprocess.check_output(['git', '-C', '/home/freiburger/Documents/prFBA', 'show', 'rescore-2026-11-start:region_validation/data/taxacc_region_summary.json']))
T, O = table(cur), table(old)
print('%-14s %s' % ('region', 'current: Fam incl, Fam excl, Gen excl, delta(unrounded) | delta(rounded) | start-tag same 4'))
for r in T:
    fi, fe, ge, dl = T[r]
    print('%-14s %.4f %.4f %.4f %+.4f | %+.3f | %s' % (r, fi, fe, ge, dl, round(fe,3)-round(ge,3), ' '.join('%.4f' % x for x in O[r])))
