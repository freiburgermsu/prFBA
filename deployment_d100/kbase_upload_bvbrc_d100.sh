#!/bin/bash
# Upload both BV-BRC-only depth-100 builds (PATRIC annotation only; E3 feature probabilities for
# multi-genome ASVs, k/n for single-genome ones) to KBase. User request 2026-11-15: "upload all of
# the new builds to KBase".
#   A  codiffusion -> workspace 261948 ("codiffusion bioreactor"; it also holds the 20-hit build under
#      other object names, which this run never touches)
#   B  EmilyKin    -> workspace 275568 ("EmilyKin prFBA synthetic genomes (depth 100)")
# Mode --new-version: each genome is saved as a NEW VERSION of the same-named (ASV md5) object that the
# 2026-11-04 mixed-annotation upload (kbase_upload_d100.sh) created as version 1. Checked read-only on
# 2026-11-15: the earlier GenomeSets 261948/6467/1 and 275568/3757/1 store every genome ref WITH its
# version (all '/1'), so they keep resolving to the earlier genomes. --prior-upload-log restricts the new
# versions to the objects that upload recorded. New GenomeSets (exact ws/obj/ver refs):
#   codiffusion_synthetic_genomes_prFBA_d100_bvbrc, emilykin_synthetic_genomes_prFBA_d100_bvbrc
# Resumable: re-run after an interruption; an ASV is done when its object's LATEST version carries this
# build's stamp (metadata genome_json_sha256 = sha256 of {asv}.json), and the GenomeSet is written only
# once every built genome is verified. Each build keeps kbase_upload_log.json in its own directory.
# DRY_RUN=1 adds --dry-run (reads the token and the workspaces, saves nothing); LOG=<path> redirects
# the log. Launch (no dry run):
#   systemd-run --user --unit=kbase-bvbrc-d100-upload --collect -p MemoryMax=16G \
#     -p WorkingDirectory=/home/freiburger/Documents/prFBA bash deployment_d100/kbase_upload_bvbrc_d100.sh
set -u
cd /home/freiburger/Documents/prFBA
PY=/home/freiburger/Documents/py_venv/bin/python
LOG=${LOG:-/home/freiburger/Documents/prFBA/deployment_d100/logs/kbase_upload_bvbrc_d100.log}
DRY=${DRY_RUN:+--dry-run}
exec >>"$LOG" 2>&1
CG=/home/freiburger/Documents/codiffusion_bioreactor/synthetic_genomes_bvbrc_d100
EG=/home/freiburger/Documents/EmilyKin/synthetic_genomes_bvbrc_d100
CG_PRIOR=/home/freiburger/Documents/codiffusion_bioreactor/synthetic_genomes_d100/kbase_upload_log.json
EG_PRIOR=/home/freiburger/Documents/EmilyKin/synthetic_genomes_d100/kbase_upload_log.json
HEAD=$(git rev-parse --short HEAD)

facts() {  # built|n_asvs|abstained|E3-applied|k/n|build commit|estimator id, from the build's own files
  $PY - "$1" <<'EOF'
import json, os, sys
d = sys.argv[1]
m = json.load(open(os.path.join(d, "synthetic_genome_manifest.json")))
b = [e for a, e in m.items() if not a.startswith("_") and e.get("status") == "built"]
e3 = sum(bool((e.get("pfunc") or {}).get("applied")) for e in b)
g = json.load(open(os.path.join(d, "mapping_provenance.json")))["pipeline"]["git"]
built = g["commit_short"] + ("" if g.get("clean") else " plus the uncommitted working-tree changes listed in mapping_provenance.json")
print(f"{len(b):,}|{m['_meta']['n_asvs']:,}|{m['_meta']['abstain']:,}|{e3:,}|{len(b) - e3:,}|{built}|{m['_meta']['pfunc']['estimator_id']}")
EOF
}
desc() {  # dataset, build dir, facts, earlier set
  local ds=$1 dir=$2 f=$3 old=$4
  local built n abst e3 kn bc est
  IFS='|' read -r built n abst e3 kn bc est <<<"$f"
  echo "prFBA synthetic genomes for the $ds 16S ASVs at the depth-100 standard (100-hit alignment records), \
built from BV-BRC's own (PATRIC, RASTtk) annotation only: $built of $n ASVs ($abst abstained: no BV-BRC reference \
passed the selection gates). Feature probabilities: E3 ($est) for the $e3 ASVs with >= 2 built source genomes; \
k/n for the $kn single-genome ASVs. Each genome is a new version of the same-named object, and this set references \
those exact versions (ws/obj/ver). The earlier mixed-annotation (PATRIC + RefSeq) k/n build remains: GenomeSet $old \
still references version 1 of the same objects. Uploaded with prFBA $HEAD; genomes built with prFBA $bc. \
Selection = exact-tie de-duplicated union + P3 gene tie-break; see mapping_provenance.json / abstained_asvs.md in $dir."
}

echo "[kbase] start $(date -Is) prFBA $HEAD ${DRY:-}"
CF=$(facts $CG) || { echo "[kbase] cannot read $CG"; echo "EXIT=1 $(date -Is)"; exit 1; }
EF=$(facts $EG) || { echo "[kbase] cannot read $EG"; echo "EXIT=1 $(date -Is)"; exit 1; }
C_DESC="$(desc "codiffusion bioreactor (V4-V5)" $CG "$CF" "codiffusion_synthetic_genomes_prFBA_d100 (261948/6467/1)") \
The 20-hit build in this narrative is GenomeSet codiffusion_synthetic_genomes_prFBA_54d1150."
E_DESC="$(desc "EmilyKin" $EG "$EF" "emilykin_synthetic_genomes_prFBA_d100 (275568/3757/1)")"

echo "[kbase] A start $(date -Is)"
$PY upload_kbase_genomes.py --genome-dir $CG --ws-id 261948 --new-version \
    --prior-upload-log $CG_PRIOR --dataset-label "codiffusion bioreactor, V4-V5" \
    --genomeset-name codiffusion_synthetic_genomes_prFBA_d100_bvbrc \
    --genomeset-description "$C_DESC" $DRY
rc_a=$?
echo "[kbase] A exit $rc_a $(date -Is)"
echo "[kbase] B start $(date -Is)"
$PY upload_kbase_genomes.py --genome-dir $EG --ws-id 275568 --new-version \
    --prior-upload-log $EG_PRIOR --dataset-label "EmilyKin" \
    --genomeset-name emilykin_synthetic_genomes_prFBA_d100_bvbrc \
    --genomeset-description "$E_DESC" $DRY
rc_b=$?
echo "[kbase] B exit $rc_b $(date -Is)"
rc=0
{ [ $rc_a -eq 0 ] && [ $rc_b -eq 0 ]; } || rc=1
echo "EXIT=$rc $(date -Is)"
exit $rc
