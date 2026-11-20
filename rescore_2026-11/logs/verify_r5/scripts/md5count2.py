import json, hashlib
d = json.load(open('/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_seq.json'))
bad = sum(1 for k, s in d.items() if hashlib.md5(s.encode()).hexdigest() != k)
up = len({s.upper() for s in d.values()})
print('md5 mismatch (as stored)', bad, 'distinct upper-cased seqs', up, 'any uppercase chars', sum(1 for s in d.values() if s != s.lower()))
ids = json.load(open('/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_ID.json'))
print('ID map keys', len(ids), 'same keyset', set(ids) == set(d), 'empty-seq header:', ids.get('d41d8cd98f00b204e9800998ecf8427e'))
import pyarrow.parquet as pq
t = pq.read_table('/home/freiburger/Documents/prFBA/index.parquet', columns=['md5','seq_len'])
print('index.parquet rows', t.num_rows, 'min seq_len', min(t.column('seq_len').to_pylist()))
