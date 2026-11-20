import json, statistics as st
P = {'pre': '/home/freiburger/Documents/prFBA_rescore_2026-11/pre_repair/region_validation/shared_reference/data/weighted_objective/per_amplicon.jsonl',
     'post': '/home/freiburger/Documents/prFBA/region_validation/shared_reference/data/weighted_objective/per_amplicon.jsonl'}
for state, p in P.items():
    rows = [json.loads(l) for l in open(p)]
    arms0 = sorted(rows[0]['arms'])
    plt1, exch, internal, same = [], [], [], 0
    ex100 = None
    for r in rows:
        a = r['arms']['real']
        plt1.append(a['internal_term_split']['gene_backed_plt1'] / a['objective'])
        exch.append(a['exchange_term'] / a['objective'])
        internal.append(a['internal_term'] / a['objective'])
        # arms at exchange weight 100 = all arms whose name does not end in _ex1
        E = [v['exchange_term'] for k, v in r['arms'].items() if not k.endswith('_ex1') and v.get('status') == 'optimal']
        if E and (max(E) - min(E)) / max(abs(max(E)), 1e-12) <= 1e-6: same += 1
    n100 = [k for k in arms0 if not k.endswith('_ex1')]
    print(state, 'rows', len(rows), 'arms', arms0)
    print('   real arm: gene-backed p<1 term / objective median %.4f%% (max %.2f%%); exchange share median %.2f%%; internal share median %.2f%%'
          % (100*st.median(plt1), 100*max(plt1), 100*st.median(exch), 100*st.median(internal)))
    print('   exchange term identical (rel spread <=1e-6) across %d ex-weight-100 arms %s on %d of %d rows' % (len(n100), n100, same, len(rows)))
