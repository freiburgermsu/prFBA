#!/usr/bin/env python
"""Per-tree and per-stage memory of the R1.3 POST unit: align the unit's 30 s cgroup samples
(logs/rv-rescore-post.mem.log: memory.current incl. page cache) with the driver log's stage stamps."""
import re, sys, json, datetime as dt
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
R = "/home/freiburger/Documents/prFBA"
P = lambda s: dt.datetime.fromisoformat(s)
samples = []
for ln in open(f"{W}/logs/rv-rescore-post.mem.log"):
    m = re.match(r"(\S+) cur=(\d+) peak=(\d+)", ln)
    if m: samples.append((P(m.group(1)), int(m.group(2))))
events, tree = [], None
for ln in open(f"{R}/region_validation_d100/logs/post_all_rescore.log"):
    m = re.match(r"\[post-all\] (\S+) POST (start|done) (\S+)", ln)
    if m:
        tree = m.group(1)
        events.append((P(m.group(3)), tree, "TREE_" + m.group(2)))
        continue
    m = re.match(r"\[run \] (.+?)\s+\((\d{4}-\d\d-\d\dT[^)]+)\)", ln)
    if m and tree:
        events.append((P(m.group(2)), tree, m.group(1).split(" (")[0].split(" ->")[0].strip()))
out = {"trees": {}, "stages": []}
for i, (t, tree, what) in enumerate(events):
    if what == "TREE_done":
        continue
    nxt = next((e[0] for e in events[i + 1:]), samples[-1][0])
    win = [c for (s, c) in samples if t <= s < nxt] or [0]
    if what == "TREE_start":
        end = next(e[0] for e in events[i + 1:] if e[1] == tree and e[2] == "TREE_done")
        tw = [c for (s, c) in samples if t <= s <= end]
        out["trees"][tree] = {"start": t.isoformat(), "end": end.isoformat(), "minutes": round((end - t).seconds / 60, 1),
                              "max_cgroup_GiB": round(max(tw) / 2**30, 1)}
    else:
        out["stages"].append({"tree": tree, "stage": what, "start": t.isoformat(),
                              "minutes": round((nxt - t).seconds / 60, 1), "max_cgroup_GiB": round(max(win) / 2**30, 1)})
json.dump(out, open(f"{W}/control/r13/post_stage_mem.json", "w"), indent=1)
for k, v in out["trees"].items(): print(f"{k:40s} {v['minutes']:6.1f} min  max {v['max_cgroup_GiB']:6.1f} GiB")
for s in out["stages"]: print(f"  {s['tree']:38s} {s['stage'][:45]:45s} {s['minutes']:5.1f} min {s['max_cgroup_GiB']:6.1f} GiB")
