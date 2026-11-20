import json, hashlib, collections, sys
d = json.load(open('/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json'))
print('keys', len(d))
lens = collections.Counter()
short = []
bad = 0
for k, s in d.items():
    L = len(s)
    if L < 20: short.append((k, L, s[:30]))
    if L < 200: lens['<200'] += 1
    if hashlib.md5(s.upper().encode()).hexdigest() != k: bad += 1
print('n<20', len(short), short[:10])
print('n<200', lens['<200'], 'n>=200', len(d)-lens['<200'])
print('md5 mismatch (upper)', bad)
