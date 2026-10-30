# Depth-100 deployments: results (2026-10-30)

The chain `dep100-chain` ran on 2026-10-29 from 20:41 to 22:33 (CDT) and exited 0
(`logs/dep100-chain.log`). Both real-amplicon deployments were rebuilt at the 100-hit standard
(`align_hits.py --force-backend cuda --k2 100 --prefilter-k 500 --gpu-cand 200`; shipped
exact-tie union with the P3 gene tie-break on; synthetic-genome build through `build_d100.py`).
**The depth-100 build is local only.** Nothing was uploaded. KBase workspace 261948 still
holds the 20-hit codiffusion build (3,224 genomes). The upload command is in README section 8.

## Summary

* **The outputs are complete and clean.** All 3,276 + 3,950 records hold exactly 100 hits.
  Every ASV has a selection. Every answered ASV has one synthetic genome, one provenance file,
  one genome object and one taxonomy file. No member genome failed, no stage returned non-zero,
  and the tie-set guard never fired (0 of 21,042 and 0 of 21,364 queries fell outside the
  tie set). The unit's `memory.peak` of 26,624 MiB equals its 26 GiB `MemoryMax`, but it was
  page cache, not a kill (section 1.4).
* **What depth alone changed (EmilyKin).** The 20-hit hits re-selected with P3 on reproduce
  every n_genomes figure. Cutting the 100-hit records to 20 hits gives nearly the same numbers.
  So the move from 20 to 100 hits changed, in effect:
  * the mean selection size, 2.034 to 2.289;
  * the largest selection, 19 to 77 genomes;
  * the cross-phylum share of larger merges, 59.5% (464/780) to 62.2% (488/784);
  * answered ASVs, 3,753 to 3,755.

  The single-genome share stayed at 62.5% (2,344 to 2,345). The two-genome cross-phylum
  share stayed at 31.9% (145/455 to 145/454).
* **What P3 alone changed (EmilyKin).** P3 changed which conspecific represents a species in
  478 ASVs. It changed no selection size and no species set. The published cross-phylum
  figures (143/451, 462/783) come from the P3-off selection. With P3 on and the current cache
  the 20-hit figures are 145/456 = 31.8% and 463/779 = 59.4%. The depth-100 build uses P3 on,
  so adopting it also removes the mismatch between the paper's Methods (P3 on) and its
  deployment numbers (P3 off).
* **The aligner path is not a confound.** The 20-hit records came from exhaustive GPU
  Smith-Waterman. The 100-hit records came from an edlib top-500 prefilter plus CUDA
  rescoring. The first 20 hits of the two differ in order in 790 of 3,950 EmilyKin records,
  yet the selections differ in 1 ASV. For codiffusion they differ in 0 ASVs.
* **Selected sets that changed with depth** (P3-on 20-hit cut versus depth 100, same code and
  cache):
  * EmilyKin: 87 of 3,753 ASVs answered at both depths. In 79 the set grew, by 965 genomes
    in all and by up to 59 in one ASV. All 79 are among the 86 records whose 20th hit tied
    the best hit. 2 sets shrank, 6 kept their size but changed members, and 2 ASVs moved
    from abstaining to answered.
  * Codiffusion: 136 of 3,224. In 130 the set grew (129 of them among the 144 tie-cut
    records), 4 shrank, 2 kept their size but changed members, and 1 ASV moved from
    abstaining to answered.
* **Some tie clouds are still cut at 100 hits.** In 25 EmilyKin records and 50 codiffusion
  records the 100th hit still ties the best. Those selections are lower bounds on the
  co-optimal set. They hold 7 to 77 genomes for EmilyKin and 0 to 79 for codiffusion, where
  one of the 50 abstains.

## 1. Validation

### 1.1 Hits (`compare/compare_*.json` → `validation`; logs `[check]`)

| | codiffusion (A) | EmilyKin (B) |
|---|---|---|
| FASTA ASVs (unique) | 3,276 (3,276) | 3,950 (3,950) |
| hits records / = FASTA set | 3,276 / yes | 3,950 / yes |
| hits per record | 100 in all 3,276 | 100 in all 3,950 |
| ranks contiguous 1..100 | yes | yes |
| GPU re-score fallbacks | 0 | 0 |
| shortlist pairs scored | 1,638,000 | 1,975,000 |

No record has fewer than 100 hits. Every ASV had at least 100 references. The file name
`asv_top20_alignment_hits.json` and its key `top20` are the pipeline's fixed schema.

### 1.2 Selection and build outputs

| | A | B |
|---|---|---|
| selection records / = hits set | 3,276 / yes | 3,950 / yes |
| answered / abstaining | 3,225 / 51 | 3,755 / 195 |
| manifest `_meta` built / abstain / failed | 3,225 / 51 / 0 | 3,755 / 195 / 0 |
| manifest member genomes `failed` | 0 | 0 |
| manifest genomes ≠ selection | 0 | 0 |
| top-level `*.json` | 3,226 = 3,225 ASV genomes + `synthetic_genome_manifest.json` | 3,756 = 3,755 + manifest |
| `provenance/`, `genome_objects/`, `taxonomies/` | 3,225 each, = answered set | 3,755 each, = answered set |
| `ASVset_taxonomies/` | empty (as in the 20-hit builds) | empty |
| zero-byte genome files | 0 | 0 |
| `provenance_check` (deployment_stats.py) | 0 mismatched | 0 mismatched |
| disk | 26 GB | 32 GB |

The only extra top-level entries are the manifest and the four subdirectories. The 20-hit
`synthetic_genomes_54d1150/` also has upload, rename and CSV side files. The local d100
directory has none of them because it was never uploaded.

### 1.3 Logs

* Every stage logged `rc=0`. The build logged `built=3225 abstain=51 failed=0` (A) and
  `built=3755 abstain=195 failed=0` (B).
* The genome-cache warm step fetched 984 genomes for A and 480 for B, with `err=0` both times.
* The gene-family prefetch needed 0 fetches for A. For B it fetched 19, with 1 empty and
  `err=0`.
* `modelseedpy imported: False`.
* Every log line matching "traceback" is the progress line "Biopython traceback of the GPU
  top-100". No line matches error, exception, Killed or OOM.
* No ASV was skipped: the hits, selection and manifest each cover the full ASV set.
* Tie-set guard: A made 21,042 `gene_count_fn` queries and B made 21,364; 0 fell outside the
  8,125- and 7,346-genome tie sets. The shared cache held 7,418/8,125 (A) and 6,709/7,346 (B)
  of those tie-set genomes.

### 1.4 memory.peak = MemoryMax (26,624 MiB): page cache, not a kill

* The cgroup's `memory.peak` was 5,604 MiB after alignment and 7,500 MiB after the 5.09 GB
  gene-cache stream. It reached 26,624 MiB during the codiffusion build, which wrote 26 GB of
  genome files.
* Between stages no pipeline process is alive, yet `memory.current` read 16,210 to
  22,583 MiB. That residue is file-backed page cache charged to the unit.
* `MemoryMax` reclaims page cache before it invokes the OOM killer. The user journal for
  `dep100-chain` holds only "Started" and "Consumed 19h 27min CPU". systemd logs an OOM kill
  explicitly, and no such message appears. The kernel journal was not readable here, so it is
  not cited.
* The completeness counts above (failed=0, 0 provenance mismatches, every record at 100
  hits) rule out a silently killed worker.

## 2. Statistics, 20 hits vs 100 hits

All rows come from `deployment_stats.py`, the counting method that reproduces the published
numbers (cross-phylum via `--taxon-source index`).

| Column | Selection | Source file |
|---|---|---|
| d20 pub | the published 20-hit selection | `stats_d20/emilykin_gpu.json` |
| d20 P3 June | `.p3.json` | `stats_d20/emilykin_gpu.p3.json` |
| d20 P3 | control: the 20-hit hits re-selected now with P3 on, against the d100 tie-set cache | `compare/stats_emilykin_ctl_d20hits_p3.json` |
| d100@20 P3 | control: the depth-100 hits cut to 20 per record, P3 on, same cache | `compare/stats_emilykin_ctl_d100top20_p3.json` |
| **d100** | the depth-100 build | `stats_d100_emilykin.json` (= `EmilyKin/bvbrc_alignment_hits_d100/deployment_stats_d100.json`) |

### 2.1 EmilyKin (B), 3,950 ASVs (the deployment the paper quotes)

| | d20 pub (P3 off) | d20 P3 June | d20 P3 | d100@20 P3 | **d100** |
|---|---|---|---|---|---|
| answered | 3,753 | 3,753 | 3,753 | 3,753 | **3,755** |
| abstaining (below family / reliability) | 197 (139 / 58) | 197 (139 / 58) | 197 (139 / 58) | 197 (145 / 52) | **195 (145 / 50)** |
| n = 1 (% of answered) | 2,344 (62.5) | 2,344 | 2,344 | 2,344 (62.5) | **2,345 (62.5)** |
| n = 2 | 599 (16.0) | 599 | 599 | 598 (15.9) | **597 (15.9)** |
| n ≥ 3 | 810 (21.6) | 810 | 810 | 811 (21.6) | **813 (21.7)** |
| mean n over answered | 2.034 | 2.034 | 2.034 | 2.034 | **2.289** (8,597 / 3,755) |
| mean n over all ASVs | 1.933 | 1.933 | 1.933 | 1.933 | 2.177 |
| largest selection | 19 | 19 | 19 | 19 | **77** |
| selections ≥ 10 / ≥ 20 genomes | 62 / 0 | | | 62 / 0 | 98 / 33 |
| cross-phylum, two-genome | 143/451 = 31.7% | 145/452 = 32.1% | 145/456 = 31.8% | 145/455 = 31.9% | **145/454 = 31.9%** |
| cross-phylum, ≥ 3 genomes | 462/783 = 59.0% | 462/777 = 59.5% | 463/779 = 59.4% | 464/780 = 59.5% | **488/784 = 62.2%** |
| cross-phylum, all multi-genome | 605/1,234 = 49.0% | 607/1,229 = 49.4% | 608/1,235 = 49.2% | 609/1,235 = 49.3% | **633/1,238 = 51.1%** |
| resolved / multi-genome merges | 1,234 / 1,409 | 1,229 / 1,409 | 1,235 / 1,409 | 1,235 / 1,409 | **1,238 / 1,410** |
| Bacteria + Archaea merges | 18 | 18 | 18 | 18 | **20** |
| None-counted variant (2 / ≥ 3) | 40.6 / 66.2% | 41.2 / 66.7% | 40.8 / 66.6% | 40.9 / 66.7% | 41.0 / 69.6% |
| rule (inferred) | union, P3 off | union, P3 on | union, P3 on | union, P3 on | union, P3 on |

### 2.2 Codiffusion (A), 3,276 ASVs (the KBase build; not quoted in the paper, used in the deck)

The first column is the KBase build from `stats_d20/codif_gpu.54d1150.json`. The other
columns follow section 2.1 with the `codif` files.

| | d20 KBase 54d1150 (P3 on) | d20 P3 off | d20 P3 | d100@20 P3 | **d100** |
|---|---|---|---|---|---|
| answered | 3,224 | 3,224 | 3,224 | 3,224 | **3,225** |
| abstaining (below family / reliability) | 52 (38 / 14) | 52 (38 / 14) | 52 (38 / 14) | 52 (41 / 11) | **51 (41 / 10)** |
| n = 1 (% of answered) | 2,188 (67.9) | 2,188 | 2,188 | 2,188 | **2,191 (67.9)** |
| n = 2 | 525 (16.3) | 525 | 525 | 525 | **520 (16.1)** |
| n ≥ 3 | 511 (15.8) | 511 | 511 | 511 | **514 (15.9)** |
| mean n over answered | 1.895 | 1.895 | 1.895 | 1.895 | **2.314** (7,462 / 3,225) |
| largest selection | 20 | 20 | 20 | 20 | **79** |
| selections ≥ 10 / ≥ 20 genomes | 71 / 1 | | | | 129 / 47 |
| cross-phylum, two-genome | 116/397 = 29.2% | 115/393 = 29.3% | 116/399 = 29.1% | 116/399 = 29.1% | **116/398 = 29.1%** |
| cross-phylum, ≥ 3 genomes | 191/484 = 39.5% | 200/494 = 40.5% | 191/484 = 39.5% | 191/484 = 39.5% | **237/492 = 48.2%** |
| cross-phylum, all multi-genome | 307/881 = 34.8% | 315/887 = 35.5% | 307/883 = 34.8% | 307/883 = 34.8% | **353/890 = 39.7%** |
| resolved / multi-genome merges | 881 / 1,036 | 887 / 1,036 | 883 / 1,036 | 883 / 1,036 | **890 / 1,034** |
| Bacteria + Archaea merges | 33 | 32 | 33 | 33 | **32** |

### 2.3 Separating depth from P3 and from the aligner path

Each comparison holds everything else fixed:

* **P3 off → on.** Same 20-hit hits; compare `d20 pub` with `d20 P3`.
  * EmilyKin: no selection size moved. 478 ASVs swapped a conspecific representative and
    none changed species set. Because the index-based lineage lookup resolves some
    representatives and not others, cross-phylum moved from 143/451 to 145/456 and from
    462/783 to 463/779.
  * Codiffusion: 432 swaps.
* **Gene-cache drift, June to now.** Same rule; compare `d20 P3 June` (or the codiffusion
  KBase build) with `d20 P3`. This gave 71 (B) and 20 (A) representative swaps, all of the
  same size and the same species set (README caveat on P3 and cache contents).
* **Aligner path.** Exhaustive GPU 20-hit records against the depth-100 records cut to 20;
  compare `d20 P3` with `d100@20 P3`.
  * The selections differ in 1 ASV for B (one set grew by 1 genome) and 0 for A.
  * The top-20 lists themselves differ more. Identical ordered lists: 3,160/3,950 (B) and
    2,666/3,276 (A). Identical sets: 3,590 and 2,974. Identical rank-1 hits: 3,916 and 3,261.
  * The differences sit almost entirely below the exact-tie pool.
* **Depth 20 → 100.** Same hits, rule and cache; compare `d100@20 P3` with `d100`. This is
  the column pair the paper can attribute to depth.

Control-selection integrity:
* The controls ran `exact_select`'s body (`S.bvbrc_gene_provider` +
  `S.gene_count_from_provider` + `S.build_audit`) with network access disabled (0 calls).
* They ran against the d100 run's tie-set cache, plus 36 (A) and 1 (B) entries streamed
  read-only from the shared cache for control tie-set genomes it lacked. One of A's 37
  missing genomes is absent from the shared cache as well.
* A tie-set guard ran: 0 queries fell outside, of 14,769 (A) and 17,651 / 17,654 (B).
* The control selections are scratch files and are not kept. Their stats are in `compare/`.

## 3. Selected-set changes per ASV (`compare/compare_*.json` → `set_changes`)

| | B: published → d100 | B: d100@20 P3 → d100 (depth only) | A: KBase 54d1150 → d100 | A: d100@20 P3 → d100 (depth only) |
|---|---|---|---|---|
| answered at both | 3,753 | 3,753 | 3,224 | 3,224 |
| identical genome set | 3,260 | 3,666 | 3,075 | 3,088 |
| changed genome set | 493 | 87 | 149 | 136 |
| grew | 80 | 79 | 130 | 130 |
| grew as a superset of the old set | 4 | 22 | 25 | 29 |
| same size, other members (same species set) | 411 (408) | 6 (3) | 15 (15) | 2 (2) |
| shrank | 2 | 2 | 4 | 4 |
| species set changed | 85 | 84 | 134 | 134 |
| abstain → answered / answered → abstain | 2 / 0 | 2 / 0 | 1 / 0 | 1 / 0 |
| growth by 1 / 2–4 / 5–9 / ≥ 10 genomes | 9 / 9 / 25 / 37 | 8 / 9 / 25 / 37 | 14 / 19 / 41 / 56 | 14 / 19 / 41 / 56 |
| genomes added (sum) / largest growth | 966 / 59 | 965 / 59 | 1,379 / 61 | 1,379 / 61 |
| records with the 20th hit tied to the best | 86 | 86 | 144 | 144 |
| … of which the set grew | 79 | 79 | 129 | 129 |

Size-bucket transitions, depth only:

| | B | A |
|---|---|---|
| 0 → 1 | 1 | 1 |
| 0 → ≥ 3 | 1 | 0 |
| 1 → ≥ 3 | 2 | 0 |
| 2 → ≥ 3 | 1 | 6 |
| ≥ 3 → 1 | 2 | 2 |
| ≥ 3 → 2 | 0 | 1 |
| unchanged bucket | everything else | everything else |

Reading:
* **Growth.** Depth mostly lets tie clouds that the 20th hit had cut grow into the
  already-multi-genome bucket. Few ASVs change bucket, so the n = 1 / 2 / ≥ 3 counts barely
  move while the mean and the large-merge cross-phylum share rise.
* **Shrinks.** The few shrinks are ASVs whose extra hits changed the consensus family or the
  reliability gate. For example, codiffusion `a8a3fd00` changed consensus family from
  Mycobacteriaceae to Bacillaceae, 5 → 2 genomes.
* **Abstain to answered.** The ASVs that moved from abstaining to answered all passed the
  reliability gate once more hits were seen.
* **Tie clouds cut at 20.** "Tie cloud cut" means the k-th hit has the rank-1 identity *and*
  align_score. This stricter definition gives 86 and 144, against the README's identity-only
  89 and 148.

**Still cut at 100.** The 100th hit ties the best in 25 EmilyKin records and 50 codiffusion
records. Their depth-100 selections hold:
* EmilyKin: 7–77 genomes;
* codiffusion: 0–79 genomes, one of them abstaining, 13 of them at exactly 29.

These selections are lower bounds on the full co-optimal set.

## 4. Files

| file | content |
|---|---|
| `stats_d100_emilykin.json`, `stats_d100_codif.json` | `deployment_stats.py` on the d100 builds (README section 7). The codiffusion file is identical to the chain's own `deployment_stats_d100.json` apart from wall time. |
| `compare_d20_d100.py` | validation, hits concordance, control selections and set changes |
| `compare/compare_{codif,emilykin}.json` | its output |
| `compare/stats_*_ctl_*.json` | `deployment_stats.py` on the control selections |

Reproduce (about 70 s and 3.3 GB peak RSS per deployment):

```
~/Documents/py_venv/bin/python deployment_d100/compare_d20_d100.py --out-dir deployment_d100/compare --scratch <dir>
```
