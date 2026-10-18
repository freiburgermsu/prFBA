#!/usr/bin/env python
"""Fetch BV-BRC contig FASTAs for the 10,000-genome validation panel.

Resumable: skips any genome whose .fna already exists and is non-empty.
Polite: bounded concurrency, retry with backoff, explicit User-Agent
(the bare python UA is what triggers BV-BRC's 406).
"""
import json, os, sys, time, threading, queue, urllib.request, urllib.error

ROOT = os.path.expanduser("~/Documents/prFBA")
OUT  = os.path.join(ROOT, "mag_validation_rnd", "contigs")
PANEL= os.path.join(ROOT, "mag_validation_rnd", "panel_manifest.json")
API  = "https://www.bv-brc.org/api/genome_sequence/"
UA   = "prFBA-validation/1.0 (academic; metabolic-model benchmarking)"
NWORK= int(os.environ.get("NWORK", "16"))

def fetch(gid, tries=4):
    url = (API + "?eq(genome_id,%s)&select(accession,sequence)"
           "&limit(50000)&http_accept=application/json" % gid)
    for a in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.load(r)
        except Exception as e:
            if a == tries - 1:
                raise
            time.sleep(2 ** a)

def write(gid, rows):
    tmp = os.path.join(OUT, gid + ".fna.tmp")
    with open(tmp, "w") as fh:
        for i, row in enumerate(rows):
            seq = (row.get("sequence") or "").strip()
            if not seq:
                continue
            # header must be globally unique per synthetic MAG downstream; keep
            # accession (unique within genome) and prefix the genome id
            fh.write(">%s|%s\n" % (gid, row.get("accession") or ("ctg%d" % i)))
            for j in range(0, len(seq), 80):
                fh.write(seq[j:j+80] + "\n")
    os.replace(tmp, os.path.join(OUT, gid + ".fna"))

def main():
    panel = json.load(open(PANEL))["genomes"]
    ids = [g["genome_id"] for g in panel]
    todo = [g for g in ids
            if not (os.path.exists(os.path.join(OUT, g + ".fna"))
                    and os.path.getsize(os.path.join(OUT, g + ".fna")) > 0)]
    print("panel %d, already have %d, fetching %d with %d workers"
          % (len(ids), len(ids) - len(todo), len(todo), NWORK), flush=True)
    q = queue.Queue()
    for g in todo:
        q.put(g)
    done = {"ok": 0, "empty": 0, "fail": 0}
    failed = []
    lock = threading.Lock()
    t0 = time.time()

    def worker():
        while True:
            try:
                gid = q.get_nowait()
            except queue.Empty:
                return
            try:
                rows = fetch(gid)
                if not rows:
                    with lock:
                        done["empty"] += 1; failed.append((gid, "no rows"))
                else:
                    write(gid, rows)
                    with lock:
                        done["ok"] += 1
            except Exception as e:
                with lock:
                    done["fail"] += 1; failed.append((gid, repr(e)[:120]))
            finally:
                n = done["ok"] + done["empty"] + done["fail"]
                if n % 250 == 0:
                    el = time.time() - t0
                    print("  %d/%d  ok=%d empty=%d fail=%d  %.0fs  (%.1f/s)"
                          % (n, len(todo), done["ok"], done["empty"], done["fail"],
                             el, n / el if el else 0), flush=True)
                q.task_done()

    ws = [threading.Thread(target=worker, daemon=True) for _ in range(NWORK)]
    [w.start() for w in ws]
    [w.join() for w in ws]
    print("DONE in %.0fs: ok=%d empty=%d fail=%d"
          % (time.time() - t0, done["ok"], done["empty"], done["fail"]), flush=True)
    json.dump(failed, open(os.path.join(ROOT, "mag_validation", "logs", "fetch_failures.json"), "w"), indent=1)

if __name__ == "__main__":
    main()
