# Feasibility pilot: in-silico binning for the 10k MAG-matching validation

Date 2026-10-08. Scratch pilot, not part of any committed analysis.

## Infrastructure present
| dependency | state |
|---|---|
| cores / free disk | 64 cores, 1.2 TB free |
| skani | v0.3.2 at `prFBA/bin/skani` |
| GTDB r226 sketch | present, 57 GB (+37 GB removable tarball) |
| BV-BRC contigs | fetchable via `/api/genome_sequence/`, 0.43 s per 4 Mb genome (set a User-Agent) |
| PGFam truth | cached for 10,000/10,000 panel genomes (`bvbrc_cache/genome_gene_families.json`, 111,789 genomes) |
| `assembly_accession` | exposed by `/api/genome/`, bulk-queryable; present for only 261/300 (87%) of a panel sample |
| scorer parameterization | TAU / knobs / ident-bin edges are module constants, not CLI flags (only `--force`, `--legacy`) |

## Operational hazard: duplicate contig headers silently drop queries
skani reports `Query_name` as the first contig header. Eight synthetic MAGs from one
source genome, all whose first contig was named `src_0`, returned **5 of 8 queries**
(1 hit each instead of `-n 5`) — reproducibly, with no warning. Eight *distinct*
genomes returned 8/8. Re-emitting the same eight with per-MAG unique contig headers
returned **8/8 at 5 hits each**. The proposed benchmark generates 12 near-identical
variants per genome, i.e. exactly this case: the synthesizer must prefix contig
headers per MAG.

## Completeness: a null axis under self-inclusion, a sharp one under self-exclusion
Source `1219030.9` (2,979 contigs, 4.22 Mb) degraded by dropping whole contigs.

Against `skani_db/bvbrc_sludge`, where the source itself is in the database, the
AF-aware best hit is the source at **ANI 100.00 and AF_query ~95 at every
completeness from 100% down to 5%** — a clean contig subset still aligns almost
entirely, so the selector's AF_query coverage gate never sees completeness.
**AF_ref ~= 0.95 x completeness** (95.3, 85.9, 66.9, 47.6, 28.6, 19.0, 9.5, 4.8) is
the quantity that tracks it, and the selector ignores AF_ref — a free completeness
estimator the pipeline currently discards.

Against `skani_gtdb_r226-v0.3` (16 queries, 18 s, Learned-ANI), the same invariance
holds for ANI to a *relative*: the top hit is GCF_001895025.1 at **99.45-99.50 ANI
across all eight completeness tiers**, a 0.05 spread. So incompleteness does not
change *which* reference is selected.

What it does change is the size of the candidate pool, from the bottom. Simulating
self-exclusion (drop every hit at ANI >= 95, then apply the AF >= 15 gate):

| completeness | self-incl. best ANI / n_AFqual | self-excl. best ANI / n_AFqual | outcome |
|---|---|---|---|
| 100% | 99.45 / 3 | 92.60 / 2 | resolves |
| 90%  | 99.45 / 3 | 92.51 / 2 | resolves |
| 70%  | 99.45 / 2 | 92.56 / 1 | resolves |
| 50%  | 99.49 / 1 | -- / 0 | **abstains** |
| 30% - 5% | 99.45-99.50 / 1 | -- / 0 | **abstains** |

The surviving relatives after exclusion sit at ~92.5 ANI with AF_query of only
16-18%, so incompleteness pushes them under the 15% gate and the pipeline flips
from resolving to abstaining between 70% and 50% completeness for this organism.
**The completeness axis is therefore informative, but only crossed with
self-exclusion, and what it moves is coverage/abstention rather than taxonomic
accuracy.** n=1 organism here; where the cliff sits depends on how marginal a given
organism's surviving relatives are, which is exactly what a 10,000-organism panel
would measure.

## Verbatim-spliced contamination manufactures exact ANI ties
10% contamination spliced from a donor genome that is itself in the database yields
ANI 100.00 to *both* the source and the donor. Naive max-ANI picks the donor at 20%
and 10% completeness; the AF-aware rule picks the source correctly (donor AF_query
8.79 < 15 gate) — a live demonstration of the guard. But it means "exact ANI ties do
not occur at genome scale" holds for the 291 real MAGs and not for synthetic
splice-contaminated data. Only 5% completeness + contamination fails, correctly, as
`low_AF_unresolved`.

Reproduce: `synth.py` builds the grid; `grid.tsv` is the skani output.
