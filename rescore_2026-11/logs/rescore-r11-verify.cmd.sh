set -e
sha256sum /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair.json /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair_d20.json /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair_d20june.json
cp -p /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair.json /home/freiburger/Documents/prFBA_rescore_2026-11/control/cl_snap_copy.json
mkdir -p /home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl
cp -p /home/freiburger/Documents/prFBA_rescore_2026-11/control/genome_gene_families.pre_repair_d20.json /home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl/snap20_copy.json
sha256sum /home/freiburger/Documents/prFBA_rescore_2026-11/control/cl_snap_copy.json /home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl/snap20_copy.json
stat -c '%n %s %Y' /home/freiburger/Documents/prFBA_rescore_2026-11/control/cl_snap_copy.json /home/freiburger/Documents/prFBA_rescore_2026-11/control/d20ctl/snap20_copy.json
