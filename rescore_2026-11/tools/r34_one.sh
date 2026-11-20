#!/usr/bin/env bash
# One R3.4 prFBA re-select (FM or TF-inside-prfba-ns): r34_one.sh <root> <name> <command...>
#  1. pre-check: the run dir's align_stamp.json matches its own ref_db.json/queries.fna (the hits are reusable);
#  2. MOVE predictions.json + _cache_ec_subset.json + _cache_pgfam_subset.json to $W/pre_repair/... (COPY the
#     small re-written inputs there too);
#  3. hits/ read-only during the run (guard: a stamp mismatch must fail, never re-align);
#  4. run the ledger argv (no --force); require "[skip] reusing" and no "[align]"; check the re-written
#     ref_db/ref_id_db/queries/midas/align_stamp are byte-identical to the copies.
set -u
root=$1; name=$2; shift 2
PY=/home/freiburger/Documents/py_venv/bin/python
W=/home/freiburger/Documents/prFBA_rescore_2026-11; R=/home/freiburger/Documents/prFBA
SR=$R/region_validation/shared_reference
d=$SR/$root/prfba/$name; b=$W/pre_repair/region_validation/shared_reference/$root/prfba/$name
mkdir -p $W/logs/r3; log=$W/logs/r3/prfba_${root}_${name}.log
t0=$(date +%s)
{ echo "# $(date -Is) $(hostname) $root/$name"; echo "# $*"; } > $log
$PY - $d >> $log 2>&1 <<'PYEOF' || { echo "$(date -Is) $root/$name PRECHECK FAILED"; exit 4; }
import hashlib, json, os, sys
d = sys.argv[1]
s = json.load(open(os.path.join(d, "align_stamp.json")))
h = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
ok = (h(os.path.join(d, "ref_db.json")) == s["ref_db_sha256"] and h(os.path.join(d, "queries.fna")) == s["queries_sha256"]
      and os.path.getsize(os.path.join(d, "hits", "asv_top20_alignment_hits.json")) > 0)
print("[precheck] align stamp consistent with ref_db.json/queries.fna and hits present:", ok)
sys.exit(0 if ok else 1)
PYEOF
mkdir -p $b
for f in predictions.json _cache_ec_subset.json _cache_pgfam_subset.json; do
  if [ -e $d/$f ]; then
    [ -e $b/$f ] && { echo "$(date -Is) $root/$name: both $d/$f and $b/$f exist" | tee -a $log; exit 5; }
    mv $d/$f $b/$f
  fi
done
for f in ref_db.json ref_id_db.json queries.fna midas_empty.csv align_stamp.json; do
  [ -e $b/$f ] || cp -p $d/$f $b/$f
done
chmod a-w $d/hits; trap 'chmod u+w '$d'/hits' EXIT
/usr/bin/time -v -o $log.time "$@" >> $log 2>&1; rc=$?
chmod u+w $d/hits
if grep -q '^\[align\]' $log || ! grep -q '^\[skip\] reusing' $log; then echo "[guard] alignment path taken or no hit reuse" >> $log; [ $rc = 0 ] && rc=9; fi
same=1; for f in ref_db.json ref_id_db.json queries.fna midas_empty.csv align_stamp.json; do cmp -s $d/$f $b/$f || { same=0; echo "[check] $f differs from the pre-run copy" >> $log; }; done
rss=$(awk -F': ' '/Maximum resident/{print $2}' $log.time)
echo "EXIT=$rc wall=$(( $(date +%s) - t0 ))s maxrss_kB=$rss inputs_identical=$same" >> $log
echo "$(date -Is) $(hostname) $root/$name rc=$rc wall=$(( $(date +%s) - t0 ))s maxrss_kB=$rss inputs_identical=$same"
exit $rc
