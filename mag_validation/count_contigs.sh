#!/bin/bash
# Per-genome contig counts of the original panel's deposited assemblies (mag_validation/contigs/*.fna,
# the 9,997 genomes fetch_contigs.py fetched; three came back empty from BV-BRC and have no file).
# Output (stdout, with a header line): genome_id <TAB> n_contigs <TAB> n_contigs_ge2500 <TAB> total_bp, sorted by file name.
# n_contigs_ge2500 counts contigs of at least synth_mags.py's MIN_FRAG (2,500 bp).
# It reads all 32 GB of FASTA once; run it niced and memory-capped, e.g.
#   systemd-run --user --scope -p MemoryMax=2G -p IOWeight=10 \
#     ionice -c3 nice -n19 mag_validation/count_contigs.sh > mag_validation/data/contig_counts.tsv
# then summarise with mag_validation/contig_counts_summary.py.
set -euo pipefail
cd "$(dirname "$0")/contigs"
printf "genome_id\tn_contigs\tn_contigs_ge2500\ttotal_bp\n"
for f in $(ls *.fna | sort); do
  awk -v g="${f%.fna}" 'BEGIN{n=0;n25=0;bp=0;L=-1}
    /^>/{ if(L>=0){ if(L>=2500)n25++; bp+=L } n++; L=0; next }
    { L+=length($0) }
    END{ if(L>=0){ if(L>=2500)n25++; bp+=L } printf "%s\t%d\t%d\t%d\n", g, n, n25, bp }' "$f"
done
