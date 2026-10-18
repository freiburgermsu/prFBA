# MAG-matching validation on the 10,000-genome panel — design v1 (FOR REVIEW)

Goal: measure the whole-genome (MAG) arm of prFBA on the SAME 10,000 BV-BRC genomes
used for the 16S amplicon arm, so the two markers are compared on one organism set
under one truth and one metric definition.

## Established facts (verified, do not re-litigate)
- Panel: `prFBA/region_validation/data/selection_10k.json` -> 10,000 genome_ids with
  taxon_id and a lineage. 8,711 Bacteria / 1,289 Archaea, 206 phyla, seed 1729.
- Truth: `prFBA/region_validation/data/benchmark_truth.json`, keyed
  `<genome_id>__<region>`, carries `truth_lineage` {Kingdom..Species} with **null at
  unnamed ranks**. 9,746 unique source genomes; 2,457 (25.2%) have null Genus.
  The amplicon arm treats a null truth rank as UNSCOREABLE (excluded from the
  denominator), not as a prediction failure.
- PGFam truth: `prFBA/bvbrc_cache/genome_gene_families.json`, 111,789 genomes,
  covers 10,000/10,000 panel genomes.
- Contigs: fetchable at `/api/genome_sequence/?eq(genome_id,X)&select(accession,sequence)`
  (needs a User-Agent; ftp.bvbrc.org is unreachable from here). ~0.43 s, ~4 MB per genome.
- `assembly_accession` is on `/api/genome/`, queryable in BOTH directions and in bulk,
  but partial: 261/300 panel genomes have one, and 2/3 sampled GTDB accessions
  reverse-resolved to a BV-BRC genome_id.
- skani v0.3.2 at `prFBA/bin/skani`. GTDB r226 sketch present (57 GB, 732,475 accessions).
- Machine: 64 cores, 1.2 TB free. Target ~61 threads (95%).
- Selection code: `prFBA/select_references.py` + `prFBA/mag_skani_hits.py` (`ANI_KNOBS`,
  `best_hit`, `classify`). ANI knobs: species 0.95, genus 0.85, floor 0.83,
  cov_min 0.15 (= AF_query gate), band 0.005.
- Amplicon-arm scorer rule (`region_validation/scripts/score_taxacc.py`), to be matched
  EXACTLY: per-rank cell is None (no-call: pred absent/ambiguous OR truth undefined),
  True, or False; `correct` = correct/called (conditional); `correct_population` =
  correct x coverage over all scoreable, where the coverage denominator is the FULL
  include-self population so exclusion abstentions cannot shrink it (anti-inflation).

## Pilot findings that shaped this design (verified, see mag_validation_pilot/FINDINGS.md)
1. skani keys `Query_name` on the FIRST CONTIG HEADER. Eight variants of one genome all
   starting `>src_0` returned 5 of 8 queries, silently. Unique per-MAG contig headers
   fixed it (8/8). Every synthetic MAG MUST have uniquely-prefixed contig headers.
2. ANI is invariant to incompleteness: 100.00 (self) and 99.45-99.50 (to a GTDB
   same-species relative) from 100% down to 5% completeness. AF_query stays ~95.
   AF_ref ~= 0.95 x completeness is the only completeness-tracking quantity, and the
   selector ignores it.
3. Therefore completeness is a NULL axis under self-inclusion and a SHARP axis under
   self-exclusion: after removing ANI>=95 hits, the surviving relatives sit at ~92.5 ANI
   with AF_query 16-18%, so incompleteness pushes them under the 15% gate. One pilot
   organism flipped from resolving to abstaining between 70% and 50% completeness.
4. Contamination spliced from a donor that is itself in the reference DB produces ANI
   100.00 to BOTH source and donor (an exact ANI tie). The AF-aware rule resolves it
   correctly, but it is an artifact of synthetic splicing.

## Proposed design
**Reference DB.** GTDB r226 (already sketched, dense enough that panel organisms have
real same-species and same-genus relatives). Gene capture requires mapping a selected
GCA/GCF accession back to a BV-BRC genome_id to look up PGFams; build that map in bulk
only for accessions that actually appear as selected hits.
Rejected alternative: sketching the 10,000 panel genomes as the reference universe. The
panel is species-stratified (one genome per species cluster, sqrt-flattened), so it
contains few same-species pairs; self-exclusion would remove little and the next-best
hit would be distant and low-AF, making the self-exclusion arm degenerate.

**Synthetic MAGs.** Drop whole contigs (never mutate bases) to completeness
{100, 90, 70, 50, 30}% x contamination {0, 10}% = 10 variants per genome = 100,000
queries. Contig headers prefixed `<genome_id>_c<comp>_x<contam>|`. Contamination contigs
spliced from a donor of a DIFFERENT phylum chosen deterministically from the panel by
seed 1729. Generate -> search -> delete in batches to bound peak disk.

**Search.** `skani search -d <gtdb> --ql <batch> -n 100 --min-af 0 -t 61`.
Persistence depth 100 from the outset (the amplicon arm's top-20 truncation is what
makes its self-exclusion figures a lower bound).

**Regimes.** self-inclusion: hits as returned. self-exclusion (tau_ANI = 95): drop any
hit whose accession maps to the source genome, whose taxon matches the source taxon, or
whose ANI >= 95; then re-apply the AF gate and re-select.

**Selection.** `mag_skani_hits.best_hit` for the AF-aware best hit and
`select_references.select_representatives` with `ANI_KNOBS` for the reference set.

**Metrics.** (a) per-rank taxonomic accuracy, conditional and population, Wilson 95% CI,
truth-undefined excluded, anti-inflation denominator; (b) PGFam gene capture: truth =
source genome's PGFam set, prediction = union of selected genomes' PGFams, macro
recall/precision/F1; (c) abstention rate and its cause (low_AF_unresolved / below floor /
no_hit); (d) paired MAG-vs-amplicon comparison on the 9,746-genome intersection.

## Questions for review
1. Is GTDB r226 the right reference DB, given gene capture needs a BV-BRC genome_id and
   the accession map is partial? What is the fallback for selected hits that do not map?
2. Does the degradation model (random whole-contig dropping) bias the result, and how?
   Real binning loses contigs non-randomly (low coverage, short, atypical composition).
3. Are contamination donors drawn from the panel acceptable given finding 4, or must they
   come from outside the reference DB? What does each choice cost?
4. Is the self-exclusion rule the correct analogue of the amplicon arm's tau=0.987, given
   GTDB accessions do not carry BV-BRC taxon_ids directly?
5. What confound would most plausibly make the headline result (the completeness at which
   a MAG beats an amplicon) wrong or uninterpretable?
6. Is anything in the metric definitions not exactly comparable to the amplicon arm?
