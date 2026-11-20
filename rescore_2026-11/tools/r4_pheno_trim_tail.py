#!/usr/bin/env python
"""Before resuming a stopped phenotype_fold pass: the parent appends one genome's rows at a time (buffered, flushed
per genome), so a stop can leave a truncated last line and/or an incomplete last genome block. Drop any unparseable
tail and the whole last genome block (it is rebuilt on resume; rows are per-genome deterministic). Every other genome
block must be contiguous. A copy of the file as stopped is kept beside the report.

usage: r4_pheno_trim_tail.py PER_GENOME_JSONL REPORT_JSON
"""
import json
import shutil
import sys

path, rep = sys.argv[1], sys.argv[2]
shutil.copy2(path, rep.replace(".json", ".per_genome.at_stop.jsonl"))
lines = open(path).read().split("\n")
good, bad_tail = [], 0
for i, line in enumerate(lines):
    if not line.strip():
        continue
    try:
        good.append(json.loads(line))
    except json.JSONDecodeError:
        if any(l.strip() for l in lines[i + 1:]):
            raise SystemExit(f"unparseable line {i} is not the tail")
        bad_tail += 1
order, seen = [], set()
for r in good:
    g = r["src_genome_id"]
    if g not in seen:
        seen.add(g)
        order.append(g)
    elif order[-1] != g:
        raise SystemExit(f"genome {g} rows are not contiguous")
last = order[-1] if order else None
keep = [r for r in good if r["src_genome_id"] != last]
with open(path, "w") as fh:
    for r in keep:
        fh.write(json.dumps(r, sort_keys=True) + "\n")
out = {"rows_at_stop_parsed": len(good), "unparseable_tail_lines": bad_tail, "genomes_at_stop": len(order),
       "dropped_last_genome": last, "dropped_rows": len(good) - len(keep), "rows_kept": len(keep),
       "genomes_kept": len(order) - (1 if last else 0)}
json.dump(out, open(rep, "w"), indent=1)
print(json.dumps(out))
