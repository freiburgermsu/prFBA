#!/bin/bash
# Upload both depth-100 synthetic-genome builds to KBase (unit kbase-d100-upload), user-approved 2026-11-04.
#   A  codiffusion -> workspace 261948 ("codiffusion bioreactor", which also holds the 20-hit build under
#      different object names; checked 2026-11-04: 0 md5-named objects there, 0 of the 3,276 d100 names present)
#   B  EmilyKin    -> workspace 275568 ("EmilyKin prFBA synthetic genomes (depth 100)", created 2026-11-04)
# Resumable: the uploader skips objects already saved (by name) and keeps kbase_upload_log.json per build.
set -u
cd /home/freiburger/Documents/prFBA
PY=/home/freiburger/Documents/py_venv/bin/python
CG=/home/freiburger/Documents/codiffusion_bioreactor/synthetic_genomes_d100
EG=/home/freiburger/Documents/EmilyKin/synthetic_genomes_d100
echo "[kbase] A start $(date -Is)"
$PY upload_kbase_genomes.py --genome-dir $CG --ws-id 261948 \
    --genomeset-name codiffusion_synthetic_genomes_prFBA_d100 \
    --genomeset-description "prFBA synthetic genomes for the codiffusion bioreactor 16S V4-V5 ASVs at the depth-100 standard (100-hit alignment records): 3,225 of 3,276 ASVs (51 abstained: no BV-BRC reference passed the selection gates). prFBA commit 315c6d7; selection = exact-tie de-duplicated union + P3 gene tie-break; see mapping_provenance.json / abstained_asvs.md in $CG. The 20-hit build in this narrative is GenomeSet codiffusion_synthetic_genomes_prFBA_54d1150."
rc_a=$?
echo "[kbase] A exit $rc_a $(date -Is)"
echo "[kbase] B start $(date -Is)"
$PY upload_kbase_genomes.py --genome-dir $EG --ws-id 275568 \
    --genomeset-name emilykin_synthetic_genomes_prFBA_d100 \
    --genomeset-description "prFBA synthetic genomes for the EmilyKin 16S ASVs at the depth-100 standard (100-hit alignment records): 3,755 of 3,950 ASVs (195 abstained: no BV-BRC reference passed the selection gates). prFBA commit 315c6d7; selection = exact-tie de-duplicated union + P3 gene tie-break; see mapping_provenance.json / abstained_asvs.md in $EG."
rc_b=$?
echo "[kbase] B exit $rc_b $(date -Is)"
[ $rc_a -eq 0 ] && [ $rc_b -eq 0 ]
