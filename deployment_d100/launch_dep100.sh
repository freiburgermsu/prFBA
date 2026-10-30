#!/usr/bin/env bash
# launch_dep100.sh -- start the depth-100 deployment chain as ONE transient systemd user unit,
# dep100-chain: deployment A (codiffusion, 3,276 ASVs) then deployment B (EmilyKin, 3,950 ASVs),
# each align -> tie-set gene cache -> prefetch -> build -> stats (deployment_d100/dep100_chain.sh).
#
#   CPU : <= 12 workers at any moment (align pool 12; prefetch/warm 8 threads; merge pool 8),
#         hard-capped by CPUQuota=1200%.          GPU: yes (align stage 2; ~0.9 GB measured).
#   RAM : MemoryMax=26G = 1.3 x the 20 GB full-run peak estimate (see README.md "Resources").
#         The chain also waits before each heavy stage until MemAvailable >= 32 GB.
#   Logs: deployment_d100/logs/dep100-chain.log (+ codif_d100.log, emilykin_d100.log)
#   Time: estimated 2.5-5.5 h in total (README.md).
#
# Preflight (refuses to launch unless FORCE=1): unit not already present, MemAvailable >= 40 GB,
# GPU memory free >= 4 GB, >= 150 GB free disk. Writes only NEW *_d100 outputs (README.md).
set -uo pipefail
PRFBA=/home/freiburger/Documents/prFBA
D=$PRFBA/deployment_d100
UNIT=dep100-chain
MEMMAX=26G
LOG=$D/logs/dep100-chain.log
FORCE=${FORCE:-0}
mkdir -p $D/logs

fail() { echo "[launch] $*"; [ "$FORCE" = 1 ] && echo "[launch] FORCE=1 -> continuing" || exit 1; }
if systemctl --user list-units --all --no-legend "$UNIT.service" | grep -q "$UNIT"; then
  echo "[launch] unit $UNIT already exists:"; systemctl --user status "$UNIT" --no-pager | head -5; exit 1
fi
grep -E '^(MemTotal|MemAvailable):' /proc/meminfo; free -g | head -2
avail=$(awk '/^MemAvailable:/{print int($2/1048576)}' /proc/meminfo)
[ "$avail" -ge 40 ] || fail "MemAvailable=${avail} GB < 40 GB (MemoryMax $MEMMAX + headroom)"
gfree=$(nvidia-smi --query-gpu=memory.free --format=csv,noheader,nounits 2>/dev/null | head -1)
[ -n "$gfree" ] && [ "$gfree" -ge 4096 ] || fail "GPU memory free=${gfree:-?} MiB < 4096 MiB (or no GPU)"
dfree=$(df -BG --output=avail /home/freiburger/Documents | tail -1 | tr -dc 0-9)
[ "$dfree" -ge 150 ] || fail "free disk ${dfree} GB < 150 GB (synthetic genomes: est. 30-60 GB per deployment)"
systemctl --user list-units --no-legend 'codif-*' 'gtdbtk*' 2>/dev/null | sed 's/^/[launch] other job: /'

set -x
systemd-run --user --unit="$UNIT" --collect \
  -p MemoryMax=$MEMMAX -p CPUQuota=1200% -p WorkingDirectory=$PRFBA \
  -E OMP_NUM_THREADS=1 -E OPENBLAS_NUM_THREADS=1 -E MKL_NUM_THREADS=1 -E PYTHONUNBUFFERED=1 \
  /bin/bash -c "/bin/bash $D/dep100_chain.sh > $LOG 2>&1; echo EXIT=\$? >> $LOG"
set +x
echo "[launch] started $UNIT; follow: tail -f $LOG ; status: systemctl --user status $UNIT"
