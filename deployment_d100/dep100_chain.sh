#!/usr/bin/env bash
# dep100_chain.sh -- the body of the dep100-chain unit: deployment A (codiffusion) then
# deployment B (EmilyKin), each align -> tieset -> prefetch -> build -> stats, one stage at a
# time (<= 12 CPU workers at any moment). B runs even if A fails; exit is non-zero if either failed.
D=/home/freiburger/Documents/prFBA/deployment_d100
mkdir -p $D/logs
echo "[$(date -Is)] dep100-chain start"; free -g | head -2
bash $D/run_codif_d100.sh > $D/logs/codif_d100.log 2>&1; A=$?
echo "[$(date -Is)] deployment A (codif_d100) exit=$A -> $D/logs/codif_d100.log"
bash $D/run_emilykin_d100.sh > $D/logs/emilykin_d100.log 2>&1; B=$?
echo "[$(date -Is)] deployment B (emilykin_d100) exit=$B -> $D/logs/emilykin_d100.log"
cg=$(cut -d: -f3 /proc/self/cgroup); [ -r /sys/fs/cgroup$cg/memory.peak ] && \
  echo "[$(date -Is)] unit memory.peak=$(( $(cat /sys/fs/cgroup$cg/memory.peak) / 1048576 )) MiB"
[ $A -eq 0 ] && [ $B -eq 0 ]
