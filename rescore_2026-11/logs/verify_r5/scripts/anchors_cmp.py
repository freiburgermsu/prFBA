import json, collections
W = '/home/freiburger/Documents/prFBA_rescore_2026-11/pre_repair/'
R = '/home/freiburger/Documents/prFBA/'
f = 'region_validation/picrust2_comparison/data/prfba_selection.json'
pre = json.load(open(W + f)); post = json.load(open(R + f))
print('keys pre %d post %d same keyset %s' % (len(pre), len(post), set(pre) == set(post)))
c = collections.Counter()
for k in pre:
    a, b = pre[k].get('exclude_0.987') or {}, post[k].get('exclude_0.987') or {}
    ans_a, ans_b = (a.get('n_selected') or 0) > 0, (b.get('n_selected') or 0) > 0
    c['answered pre'] += ans_a; c['answered post'] += ans_b
    if ans_a != ans_b: c['answered status differs'] += 1
    if ans_a and ans_b:
        if a['reps'][0] != b['reps'][0]: c['anchor differs'] += 1
        if a['reps'] != b['reps']: c['rep set/order differs'] += 1
    if a != b: c['exclude record differs (any field)'] += 1
    if pre[k].get('include') != post[k].get('include'): c['include record differs'] += 1
print(dict(c))
# which fields differ?
diff_fields = collections.Counter()
for k in pre:
    a, b = pre[k].get('exclude_0.987') or {}, post[k].get('exclude_0.987') or {}
    for x in set(a) | set(b):
        if a.get(x) != b.get(x): diff_fields[x] += 1
print('differing exclude fields:', dict(diff_fields))
