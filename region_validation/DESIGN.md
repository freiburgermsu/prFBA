# 16S Hypervariable-Region Validation — Implementation Specification

> **USER DECISIONS (authoritative — override anything below that conflicts):**
> 1. **Scope:** taxonomic accuracy AND % gene capture (full study).
> 2. **Archaea:** kept, stratified by domain (only V4 / V4-V5 amplify them; other regions are Bacteria-scoped).
> 3. **Execution:** validate the whole pipeline on a ~500-genome pilot first, then launch the full 10k.
> 4. **Pipeline UNALTERED; report BOTH self scenarios for every examined sequence.** Each amplicon passes
>    through the unaltered pipeline against the full DB (source organism included) — the realistic experimental run.
>    From that single alignment we report, per amplicon, TWO numbers:
>    - **self-inclusion (primary, realistic):** selection/scoring on the unfiltered hits (source in the DB).
>    - **self-exclusion (generalization):** post-hoc only — drop hits that are the source organism or a
>      near-identical neighbor (`hit.genome_id == src_genome_id` OR `hit.taxon_id == src_taxon_id` OR
>      `hit.identity ≥ τ`, τ = 0.987 species threshold; also report τ = 0.995), then re-run the *pure*
>      `select_references.select_representatives` on the filtered record and re-score — performance as if the
>      source organism were NOT in the reference set (novel-organism generalization).
>    The pipeline RUN is never altered (exclude-self is analysis on the same hits). No `16S_metadata.json` parse /
>    region-collision index needed — the identity-τ filter subsumes near-identical-16S relatives. KEEP figure F11
>    (self-inclusion vs self-exclusion delta per region = the leakage/generalization gap). Both scenarios appear
>    on every per-record output and every region aggregate.
> 5. **Persistence/resumability:** ALL intermediate files live under `/home/freiburger/Documents/prFBA/region_validation/data/`;
>    every stage is checkpointed and scripts skip stages whose outputs already exist, so the effort resumes across sessions.

A single, coherent, implementation-ready scheme. Every verifier correction is folded in: the V6–V8/V7–V9 mislabel, the corrected primer revcomps (K not RR; single-R 1492R), the self-hit redefinition against per-operon ground truth + region-level identity + a near-identical tier, the gene-ID-set persistence blocker, archaea scoping, the corrected cost model (~0.6 GB GPU, hyphenated BV-BRC host, ~4.5k gated PGFam fetches, dedup-first), and reproducibility pins.

All paths absolute. Sole interpreter: `/home/freiburger/Documents/py_venv/bin/python` (`$PY` below). Long jobs run in detached `tmux`.

---

## 1. Overview and goal

**Question.** For 16S amplicon surveys, how much does the choice of hypervariable region trade off **taxonomic accuracy** (can the pipeline place an amplicon on the correct lineage?) against **functional gene capture** (does the resolved reference set reconstruct the source organism's gene content?) — quantified **per region**, **honestly** (excluding the source organism from the reference set), with calibrated uncertainty.

**Design.** A closed-loop in-silico benchmark:

1. **Select** 10,000 taxonomically diverse, full-length-16S BV-BRC genomes (deterministic).
2. **Extract** each region's amplicon from every source genome's 16S by in-silico PCR (one front-loaded combined FASTA).
3. **Align** all amplicons against the 459k-reference 16S DB in one GPU-prefiltered pass.
4. **Select** representative genomes per amplicon (`select_references.py`), in two branches: **include-self** (upper bound) and **exclude-self** (honest generalization).
5. **Score** taxonomic accuracy (per-rank concordance vs. the source's true lineage) and gene capture (PGFam set recall/precision/F1/Jaccard vs. the source genome).
6. **Render** the full figure slew, every accuracy/gene-capture number stratified by region, tier, and domain, with Wilson CIs.

**Two populations, reconciled (verifier C3).** This study evaluates the **10,000 selected source genomes**, each treated as a synthetic "query organism." This is distinct from the existing 3,950-ASV community study. The 10k genomes are the queries; the **full 459,301-sequence BV-BRC 16S DB is the reference store** (region-matched per §4). The join is: source genome → its 16S → region amplicon → aligned against the *whole DB* → self removed per §5 → selection scored. The 3,950-ASV artifacts are reused only as **schema references** (audit format, tier thresholds, render conventions), not as the evaluation set.

**Scope decision (verifier C2 — archaea).** The core region panel uses bacterial-domain primers. Only **515F/926R (V4–V5)** and **515F/806R (V4)** reliably amplify archaea; 27F, 341F/805R, and 1100F have documented archaeal mismatches. Therefore: **report a `domain` column on every record and stratify every metric Bacteria/Archaea.** Archaeal numbers in the 6 bacterial-only regions are expected to be near-zero *extract-rate* (not "inaccuracy") and must be labeled as amplification failure, not pipeline failure. The selection's √-flattening over-weights the 4,964-archaea pool; we keep them so V4/V4–V5 archaeal performance is measurable, but document that 6/8 regions are Bacteria-scoped.

---

## 2. Genome selection (final algorithm — 10,000 diverse)

**Inputs:** `/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_md5_ID.json` (`{md5: representative_header}`), `16S_md5_seq.json` (`{md5: seq}`); taxonomy `/home/freiburger/Documents/prFBA/{nodes,names}.dmp` via `taxopy 0.14.0`.

**Measured pool:** 460,117 unique 16S; 147,286 genomes with ≥1 full-length (≥1400 bp) 16S; 140,912 resolve to a phylum (135,948 Bacteria / 4,964 Archaea); 201 phyla / 1,447 families / 36,419 species. The 36,419-species ceiling ⇒ 10k draws ~27% of species with one genome per species — comfortably satisfiable.

**Constants:** `TARGET=10_000`, `FULL_LEN=1400`, `FALLBACK_LEN=1200`, `SEED=1729`, `α=0.5` (√-flattening), `CAP_FRAC=0.35`, `RANKS=["domain","phylum","class","order","family","genus","species"]`. **Determinism rule:** iterate only over `sorted(...)`; one `random.Random(SEED)` consumed in sorted order ⇒ byte-identical output.

**Algorithm.**
1. **Candidate table (one pass).** Parse `genome_id` from each representative header via `header[rfind("[")+1:rfind("]")]` then `rsplit("|",1)`; validate `^\d+\.\d+$`. Per `genome_id` track `best_len`, `best_md5` (longest seq, tie → lexicographically-smallest md5), `n_16s`, `n_full` (count ≥1400).
2. **Lineage with rank fallback.** `taxon_id = int(gid.split('.')[0])`. Read `domain` with `superkingdom` fallback; require `domain ∈ {Bacteria,Archaea}` and a non-null phylum, else drop. For each missing intermediate rank synthesize a stable placeholder `unclassified_<parentrank>:<parentname>` so genus-less environmental clades form their **own** bin (genus missing ~32%, family ~24%) instead of collapsing. `species` falls back to `sp_taxon_<tid>`.
3. **Species dedup → one representative.** Keep only `best_len ≥ FULL_LEN`. Per `(domain,phylum,species)`, pick the rep by the deterministic ladder: `(-n_full, -best_len, -int(assembly_index), genome_id)`.
4. **Stratified allocation (recursive, capped, deterministic).** Build the nested tree domain→…→species. Allocate top-down with `weight(child)=avail(child)**α`, `share=budget·w/Σw`, `alloc=min(round(share), avail, ceil(CAP_FRAC·budget))`. Fill to the exact budget with the **largest-remainder (Hamilton)** method, ties broken on sorted child key; relax caps only when a clade is genuinely exhausted. At leaves, take all species if `budget≥|species|`, else `rng.shuffle(sorted(species))` and take a re-sorted prefix.
5. **Guarantees.** Full-length is enforced *before* allocation (every leaf owns a concrete `best_md5≥1400`). Largest-remainder makes the root sum exactly `TARGET`; assert `len(chosen)==TARGET`. `FALLBACK_LEN=1200` is an optional second pass admitting `1200≤best_len<1400` **only** for species with no ≥1400 member (backfill for 16S-poor clades).
6. **Output manifest** `bvbrc_16s_validation_10k.json`: per genome `{genome_id, taxon_id, best_16s_md5, best_16s_len, n_full_length_16s, domain..species, seed/alpha/cap_frac}`.

Expected coverage: all 201 phyla and effectively all 1,447 families present; top phyla flattened from raw 37/20/14% toward ~13–16% each.

---

## 3. Region panel (FINAL corrected table)

IUPAC: R=A/G, Y=C/T, S=G/C, W=A/T, K=G/T, M=A/C, B=C/G/T, D=A/G/T, H=A/C/T, V=A/C/G, N=ACGT. Sequences in published 5'→3'. The extractor matches the **forward primer on the sense strand** and the **reverse-complement of the reverse primer on the sense strand**, both computed at runtime by `revcomp()` (so the `revc` column is informational — the code never reads it).

**Corrections applied (verifier PRIMERS):** (1) the old "V6–V8" entry was a **mislabel** — 1100F's 3' end (~1114) sits downstream of V6 (ends 1043), so 1100F/1492R amplifies **V7–V9**, identical to the "V7–V9" entry. The two are collapsed to one **V7–V9** row. (2) The two prose revcomp typos are corrected here: 926R sense-revcomp uses **K** (`...AAACTYAAAKRAATTGRCGG`), and 1492R has a single R (`AAGTCGTAACAAGGTARCCGTA`). (3) An **optional true V6–V8** row (926F/1378R) is listed separately, included only if V6 must be covered.

| Region | Fwd (name, 5'→3') | Rev (name, 5'→3') | Rev-comp on sense (corrected) | E.coli coords | Amplicon bp | Insert bounds [lo,hi] | Priority | Domain |
|---|---|---|---|---|---|---|---|---|
| FullLength16S (pos. ctrl) | 27F `AGAGTTTGATCMTGGCTCAG` | 1492R `TACGGYTACCTTGTTACGACTT` | `AAGTCGTAACAAGGTARCCGTA` | 8–1507 | ~1500 | 1380–1520 | core | Bacteria |
| V1–V2 | 27F `AGAGTTTGATCMTGGCTCAG` | 338R `TGCTGCCTCCCGTAGGAGT` | `ACTCCTACGGGAGGCAGCA` | 8–355 | ~330 | 270–340 | core | Bacteria |
| V1–V3 | 27F `AGAGTTTGATCMTGGCTCAG` | 534R `ATTACCGCGGCTGCTGG` | `CCAGCAGCCGCGGTAAT` | 8–534 | ~525 | 420–520 | core | Bacteria |
| V3–V4 | 341F `CCTACGGGNGGCWGCAG` | 805R `GACTACHVGGGTATCTAATCC` | `GGATTAGATACCCBDGTAGTC` | 341–805 | ~465 | 390–450 | core | Bacteria |
| V4 | 515F-Parada `GTGYCAGCMGCCGCGGTAA` | 806R-Apprill `GGACTACNVGGGTWTCTAAT` | `ATTAGAWACCCBNGTAGTCC` | 515–806 | ~290 | 240–300 | core | **Bact+Arch** |
| **V4–V5 (anchor)** | 515F-Parada `GTGYCAGCMGCCGCGGTAA` | 926R-Quince `CCGYCAATTYMTTTRAGTTT` | `AAACTYAAAKRAATTGRCGG` | 515–926 | ~411 | 300–450 | core | **Bact+Arch** |
| V4–V5_944R (variant) | 515F-Parada `GTGYCAGCMGCCGCGGTAA` | 944R `GAATTAAACCACATGCTC` | `GAGCATGTGGTTTAATTC` | 515–944 | ~430 | 320–470 | variant | Bact+Arch |
| **V7–V9** (was "V6–V8" + "V7–V9") | 1100F `YAACGAGCGCAACCC` | 1492R `TACGGYTACCTTGTTACGACTT` | `AAGTCGTAACAAGGTARCCGTA` | 1100–1507 | ~410 | 350–430 | core | Bacteria |
| *(opt) true V6–V8* | 926F `AAACTYAAAKGAATTGRCGG` | 1378R-class | — | ~907–1100 | — | tune | optional | Bacteria |

**Code change required** in `/home/freiburger/Documents/prFBA/insilico_pcr.py` `REGION_PANEL`: collapse the two byte-identical 1100F/1492R entries (lines 128–141) to one named `V7-V9`. No extraction-logic change — `revcomp()` already computes the correct sense-strand patterns at runtime (verified: `REVC_926R = AAACT[CT]AAA[GT][AG]AATTG[AG]CGG`, i.e. the K form). Re-tune `lo/hi` on a phylum-stratified ≥5k sample (verifier R2) and record observed p1/p99.

---

## 4. In-silico amplicon extraction (algorithm + front-loaded combined FASTA)

**Module:** `/home/freiburger/Documents/prFBA/insilico_pcr.py`, **panel mode** (`--combined-fasta`). Reuse it directly — do **not** write bespoke `build_combined_fasta.py` (verifier F8: panel mode already emits provenance-bearing keys).

**Per-region extractor** (`make_extractor`): forward regex matched exact-degenerate first, then fuzzy `{e<=emax}` (default `emax=3`) fallback; reverse primer's revcomp matched downstream; insert length gated by `lo ≤ r.start()-f.end() ≤ hi`. `include_primers=True` (panel default) emits the primer-to-primer amplicon `seq[f.start():r.end()]`; non-amplification (any primer absent or insert out of band) → record skipped.

**Required extraction upgrades (verifiers G1, G2):**
- **Both-strand pass.** "Both strands" is currently a misnomer — only the given sense is scanned. Add: if `extract(seq)` is `None`, retry `extract(revcomp(seq))`. Catches minus-strand-deposited 16S genes (otherwise a silent selection bias).
- **Mismatch accounting.** Record per-amplicon `fwd_edits`, `rev_edits`; provide a **strict mode** (`--emax 1` or a total-edit cap) so a 3+3-mismatch "amplification" doesn't count identically to a perfect match. Verify fuzzy-hit amplicon boundaries (indels can frame-shift the insert).

**Front-loaded combined output.** One pass over all 10k source genomes' full-length 16S against the whole panel → **ONE** `combined_amplicons.fasta`, one record per amplifying `(genome, region)`, id = `>{genome_id}__{region}` (unique, provenance-native — no header-collision risk from insert-hashing). Sibling `combined_amplicons.stats.json` carries per-region primer pair, runtime `rev_revcomp`, insert bounds, `n_amplicons`, `extract_rate`, and amplicon-length p5/median/p95 — **per domain** (add a Bacteria/Archaea split to the existing `lstats`).

**Provenance + dedup map (verifiers F8, cost).** Build `amplicon_provenance.json` keyed by `{genome_id}__{region}` → `{src_genome_id, region, src_amplicon_md5, src_fulllen_md5_set}` and an inverse `amplicon_md5_fanout.json` mapping `region_amplicon_md5 → [genome_id__region, …]`. The fan-out enables (a) **dedup** of byte-identical inserts before alignment (the single highest-leverage speedup) and (b) **region-collision self-detection** (§5). `src_fulllen_md5_set` and per-operon md5s come from `16S_metadata.json` (see §5), built once and cached.

Invocation:
```
$PY /home/freiburger/Documents/prFBA/insilico_pcr.py \
  --fasta /home/freiburger/Documents/EmilyKin/16S_validation/src16s_10k.json \
  --combined-fasta /home/freiburger/Documents/EmilyKin/16S_validation/combined_amplicons.fasta \
  --emax 3
```
(`src16s_10k.json` = `{header: full_16S}` for the 10k selected genomes, one record per operon, header carries `genome_id`.)

---

## 5. Pipeline run (commands, parallel/GPU settings, self-hit handling)

**Hardware (confirmed):** RTX 5070 Ti, 15.7 GB free; 64 cores; 188 GB RAM. CuPy GPU SW path auto-selected. `MAXQ=600` in `gpu_align.py` safely covers all panel inserts ≤525 bp (the FullLength 1500 bp control **exceeds** MAXQ → either bump `MAXQ` to 1600 before running, or run the full-length control in a separate pass; this is a required pre-flight check, not optional).

**5.0 — Dedup first (cost-critical).** Before alignment, collapse identical inserts via `amplicon_md5_fanout.json`; align only the unique set (often a 2–5× reduction across conserved-flank V-regions), then fan results back out. The aligner keys on header, cost is per unique insert.

**5.1 — Build the per-operon self-truth maps (verifier FINDINGS 1–2, the central fix).** The reference DB is md5-deduplicated to **one representative `genome_id` per sequence** (460,117 md5 → 315,050 reps). Therefore the naïve "drop other genome_ids sharing the source md5" is structurally inert, and the real leak — a *different* genome being the DB representative for the source's sequence, and *region-level* collisions where genomes differ across full 16S but are identical within the amplicon — is missed. Fix: parse `/home/freiburger/Documents/codiffusion_bioreactor/model_inputs/16S_metadata.json` **once** (~25 min, 3.4 GB RSS) into slim cached maps and never re-parse:
- `md5_to_genomes.json`: `{operon_md5: [genome_id,…]}`
- `genome_to_md5s.json`: `{genome_id: [operon_md5,…]}`
- From §4's in-silico PCR over the **whole DB** (already run to build the combined FASTA): `region_amplicon_md5 → {ref genome_ids}` per region (free region-level collision index).

**5.2 — Single fully-parallel GPU alignment pass** (detached tmux):
```
$PY /home/freiburger/Documents/prFBA/align_hits.py \
  --fasta   /home/freiburger/Documents/EmilyKin/16S_validation/combined_amplicons.fasta \
  --outdir  /home/freiburger/Documents/EmilyKin/16S_validation/hits_withself \
  --prefilter-k 500 --k2 20 --gpu-cand 100 --workers 60 \
  --min-ref-len 200 --force-backend cuda \
  2>&1 | tee .../hits_withself/align.log
```
Settings + rationale: `--workers 60` (measured 59.5× speedup at 60/64 cores; Stage-1 edlib prefilter and Stage-3 biopython rescore are the CPU-parallel stages; GPU SW is a single stream). `--prefilter-k 500` — **do not lower** (`gpu_align.py` checks for hits outside top-500; tie clusters reach ~200). `--gpu-cand 100` — rescore margin above `k2=20`. `--force-backend cuda` — pin + audit. GPU memory is a **non-issue** (measured ~0.6 GB peak; `pairs_batch=4_000_000` copies int32 scores back per batch). Note Stages 1→2→3 run **strictly sequentially** (GPU idle during edlib and vice-versa) — do not expect overlap.

**Output** `asv_top20_alignment_hits.json`: per amplicon, top-20 hits each with `genome_id`, `md5`, `taxon_id`, `organism`, `identity`, `align_score`, `aligned_len`, `edlib_identity`, `lineage{Kingdom..Species}`. Both `genome_id` and `md5` present ⇒ self-filtering is pure post-processing.

**Top-20 depth guard (verifier F8).** If self + sisters fill the top-20, filtering can empty an amplicon's hit list. Log `n_removed_self` and the **post-filter empty rate**; if high, base the exclude-self branch on the existing `asv_top500_alignment_hits_gpu.json` (1.1 GB) instead of top-20.

**5.3 — Include-self selection branch (upper bound):**
```
$PY /home/freiburger/Documents/prFBA/select_references.py \
  --hits .../hits_withself/asv_top20_alignment_hits.json \
  --out  .../selection_withself.json --gene-provider auto
```

**5.4 — Exclude-self via in-memory re-selection (NOT a JSON-rewriter; verifier FINDING 4).** Drop the standalone `filter_self_hits.py` entirely (it used `best_identity=None`, which crashes `_tier(None,k)` on the first all-self amplicon, and mutated `best_align_score` which `select_representatives` never reads). Keep **one** code path: an in-memory `exclude_self(asv_record, …) → select_representatives(…)`. `select_representatives` is a **pure function of one record** (reads `best_identity` + `top20`), so re-selection is exact, cheap, no re-alignment.

**Self-exclusion is a configurable cascade; report the whole leakage ladder, headline = strictest (verifiers FINDINGS 1–3,5,6):**
1. `exact_self` — hit `genome_id == src_genome_id`, **or** hit `md5 ∈ genome_to_md5s[src_genome_id]` (drops *all* operons of the source, multi-copy-safe).
2. `region_identical` — hit's reference is in `region_amplicon_md5_fanout[src_amplicon_md5]` (genomes byte-identical to the source **within this region's amplicon**, even if their full 16S differs — the leak unique to a hypervariable-region study).
3. `same_taxon` — hit `taxon_id == src_taxon_id` (other assemblies of the same organism).
4. `near_identical` — hit `identity ≥ τ` (τ = species threshold 0.987, also report 0.995) regardless of id/md5.

For each level, drop matching hits, **re-rank 1..n**, recompute `best_identity = kept[0].identity if kept else 0.0` (0.0 → `below` tier → clean abstain; never `None`), and re-run `select_representatives`. Persist `n_removed_self` and the removed `(genome_id, md5, identity, level)` per amplicon for audit. Emit `selection_noself_<level>.json` for each level; the **headline honest number uses `near_identical`** (the strictest), with `exact_self`/`region_identical`/`same_taxon` as the leakage gradient.

```
$PY /home/freiburger/Documents/EmilyKin/16S_validation/scripts/select_excludeself.py \
  --hits .../hits_withself/asv_top20_alignment_hits.json \
  --prov .../amplicon_provenance.json \
  --md5maps .../{md5_to_genomes.json,genome_to_md5s.json} \
  --region-fanout .../amplicon_md5_fanout.json \
  --levels exact_self,region_identical,same_taxon,near_identical \
  --tau 0.987 --outdir .../selection_noself/
```

**(Optional) PGFam pre-warm** before selection so novelty uses exact gene-set Jaccard (`prefetch_gene_families.py`, 12 workers; §7).

---

## 6. Accuracy metrics (definitions; include/exclude-self)

**The one comparability rule.** Truth and predicted lineages are **both** resolved through the *same* `edlib_biopython_hits.lineage_for(taxon_id)` (taxopy + prFBA taxdump), which builds `Kingdom = superkingdom or kingdom or domain`. Under modern NCBI, `Kingdom→'Pseudomonadati'` for bacteria — fine, because both sides use the identical function. Never compare against hand-written name strings. Normalize both sides with `select_references._norm` (strips `Candidatus`/`Ca.`/`midas_`, casefolds) before `==`.

**Truth.** Per amplicon, `src_taxon_id = int(src_genome_id.split('.')[0])`, `truth = lineage_for(src_taxon_id)`. Persist `benchmark_truth.json`: `{ampliconKey: {src_genome_id, src_taxon_id, region, domain, truth_lineage{7 ranks}}}`. Non-numeric/MAG source ids → all-None lineage → **excluded from the benchmark** (unscoreable, logged, not a failure).

**Ranks.** `SCORED = [Phylum,Class,Order,Family,Genus,Species]`; Kingdom is a sanity gate.

**Two predictors per amplicon, both re-resolved via `lineage_for(taxon_id)`:**
- **Anchor** = `selected[0]` (`role=="anchor"`); `n_selected==0` → `ABSTAIN`.
- **Consensus** = per-rank majority over `selected[]`; ties → `AMBIG`; record `consensus_strength[rank]=top/n`.

**Three-state per-rank cell** (the backbone of every aggregate):
```
rank_correct(truth_v, pred_v):
  if pred_v is None or AMBIG: return None      # no-call (abstention) — excluded from %correct denom, counted in coverage
  if not truth_v:             return None      # truth undefined at this rank — don't penalize
  return _norm(pred_v) == _norm(truth_v)
```
`None ≠ False`: `None` is a no-call; `False` is an active misassignment.

**Deepest correct rank** = deepest `SCORED` rank with `True` **and a contiguous correct prefix** (a deep match over a wrong Family does not count as resolution to that depth).

**Misassignment characterization** (only where a cell is `False`): `first_wrong_rank` (shallowest False), `shared_to` (deepest agreeing rank), `off_by_ranks`, and `kind ∈ {sibling, cross_lineage, overconfident, self_displaced}`. `self_displaced` (self was in top-20 but anchor is a wrong neighbor) is the selector-quality KPI. `overconfident` (pred deeper than truth defines) must **never** count as `True`.

**Per-region aggregates** — every one computed in **both** include-self and exclude-self, and stratified by **domain** and **tier**:
1. **(correct, coverage) pair** per rank, per predictor: `correct = #{True}/#{not None}`, `coverage = #{not None}/N`. Always report the pair; never `correct` alone (verifier FINDING 7).
2. **Self-recovery** (include-self only): `#{self_selected}/#{self_in_top20}`, plus `self_top20_rate`, `self_anchor_rate`.
3. **Rank-resolution profile**: histogram/CDF of `deepest_correct_rank` (the headline stacked bar).
4. **Abstention rate** split by reason (`below_family`, `reliability`), cross-tabbed against `best_identity` bins (Yarza floors 0.987/0.945/0.865).
5. **Misassignment spectrum** (`kind` counts, `first_wrong_rank` histogram, mean `off_by_ranks`).
6. **Anchor-vs-consensus delta** per rank.
7. **Tier-conditioned accuracy** (a family-tier amplicon shouldn't be graded on Species).

**Include vs exclude-self (verifier FINDING 7 — the anti-inflation rule).** Exclude-self selectively pushes amplicons below threshold (their ~1.0 self hit vanished), shrinking the denominator toward easy cases. So **compute exclude-self `correct` against the include-self (full) denominator** — abstentions count as no-call against the full set — so the include→exclude *drop* reflects reality, not denominator drift. Report all four exclusion levels (§5.4); headline = `near_identical`.

**Uncertainty (verifier A1).** Every proportion carries a **Wilson 95% CI** (`statsmodels.stats.proportion.proportion_confint(..., method='wilson')`). `n_evaluable` differs across regions (amplification differs), so cross-region accuracy differences are tested with a **paired McNemar test on the common-amplified-in-all-regions subset**, not eyeballed heatmap deltas. Also report accuracy on that common subset as a survivorship-bias sensitivity check (verifier R2).

**Reproducibility (verifier A3).** Selection is already deterministic. For the alignment/embedding numbers, pin `torch.use_deterministic_algorithms(True)`, fixed `PYTHONHASHSEED`, deterministic top-k tie order, and version-stamp GPU/driver/lib in `region_meta.csv`.

**Tie/ambiguity (verifier R3).** Short amplicons produce many references at ~identical identity. Record per amplicon the **tie-cluster size** (refs within ε of best identity); make the top-k cutoff deterministic; add an **ambiguity-adjusted concordance** (correct if the true name is anywhere in the full tie cluster, not just the truncated top-20).

**Outputs:** `taxacc_per_record.jsonl` (one record, all self-levels nested), `taxacc_region_summary.json` (`region→predictor→self_level→domain→tier→metrics+CI`), `taxacc_region_summary.csv` (flat, for plotting).

---

## 7. Gene-capture metrics (definitions, PGFam fetch, exclude-self)

**BLOCKER first (verifier C1).** The synthetic-genome manifest's `gene_union` stores **counts only** (`n_functions_union`, `n_core`, …) — no gene-ID sets — so set recall/precision are *uncomputable* as-is. **Required:** modify the gene-union builder (`build_synthetic_genomes.py`) to persist `gene_union.pgfam_ids: [sorted pgfam_id]` per genome. Until this lands, F4/F5b/F8/F9/F12 cannot be produced. We define the metric over **PGFam sets** fetched per `genome_id`.

**Definitions.** `PGF(g)` = set of distinct `pgfam_id` over genome g's CDS. `T = PGF(src_genome)` (truth). `P = ⋃_{r∈R_eff} PGF(r)` (union over selected reps). With `t=|T|`, `p=|P|`, `i=|T∩P|`:

| Score | Formula | Reads as |
|---|---|---|
| **recall** (= "% gene capture") | `i/t` | fraction of true genes recovered |
| **precision** | `i/p` | fraction of predicted genes that are real |
| **F1** | `2i/(t+p)` | balanced |
| **Jaccard** | `i/(t+p−i)` | set agreement |

Recall is reported as `100·i/t`; the other three guard against recall being gamed by selecting more genomes.

**Two modes (verifier FINDING 6 — exclude at organism level, not just genome_id).**
- **`exclude_source` (PRIMARY, honest):** `R_eff = R \ {reps with taxon_id == src_taxon_id}` (and optionally the `near_identical` set), **not** merely `genome_id == g*`. Removing only the exact genome_id leaves near-identical sisters in, inflating recall to ~1.0 and falsely labeling it generalization.
- **`include_source` (diagnostic ceiling):** `R_eff = R`.

**Edge outcomes (every cell classified, never silently dropped):** `truth_empty` (T=∅, 16S-only genome) → all NaN, **excluded from means**; `no_amplicon` (region absent in this copy) → tracked separately; `abstained` (n_selected=0) → recall/F1/Jaccard=0, precision NaN; `source_only` (only rep was self) → recall 0; `pred_empty` (all reps gene-less) → recall 0; **`pgfam_missing` (any rep's fetch failed) → poison the whole cell** (flag, exclude from denominators — do not silently shrink P, which would inflate precision, verifier FINDING 6).

**Multi-copy 16S aggregation.** Run per operon copy that amplifies; aggregate to one `(genome,region)` cell by **union of selected reps** (primary) and report the per-copy recall median/IQR (robustness). No amplifying copy → `no_amplicon`.

**Reporting granularity.** Per `(genome,region,mode)` cell table; per-region **macro-mean** (genome-equal) *and* **micro/pooled** (`Σi/Σt`). Stratify by **tier** (species/genus/family — gene-union size is tier-dominated, so cross-region gene-capture deltas can be tier-mix artifacts; faceting by tier is mandatory, verifier C2/figures) and by **best-non-self-hit identity bins** (the recall-vs-identity curve binned by *distance to nearest genuinely-different reference*, verifier FINDING 6 — not include-self best_identity which is ~1.0 everywhere). Report explicit per-outcome denominators.

**Circularity caveat (verifier completeness).** "Gold = source genome's PGFam set" is the true source content (good). If instead a full-16S-resolved gold is ever used, label it **relative** capture (vs. longer amplicon), not absolute.

**PGFam fetch plan (verifier cost corrections).**
- **Host is hyphenated:** `https://www.bv-brc.org/api/genome_feature/` (the design's `www.bvbrc.org` 404s).
- Per genome: `?eq(genome_id,<gid>)&eq(feature_type,CDS)&select(patric_id,pgfam_id,plfam_id,product)&limit(25000)&http_accept=application/json`. `limit(25000)` truncates >25k-CDS genomes — **flag** any rep hitting the cap.
- **Single cache** `/home/freiburger/Documents/prFBA/bvbrc_cache/genome_gene_families.json` (`.json.gz` supported), `gid→set(pgfam)`, loaded via `select_references.bvbrc_gene_provider` — one namespace for truth and predicted ⇒ apples-to-apples.
- **Threaded, resumable prefetch** (`prefetch_gene_families.py`, `workers=12`, 3× retry/backoff, save-every-200 — atomic, unlike the lazy provider). Two passes into the same cache: **truth** = the 10k source genome_ids (~35–45 min, ~300 MB); **predicted reps** = the *distinct selected* reps (gated to stage-7 candidates, ~79% already cached → realistically <5k new fetches, ~15–25 min). **Not** 10k blind calls; ~4.5k gated. Total cache ~600–700 MB uncompressed / ~100 MB gzipped. Do **not** raise workers past ~12–16 (BV-BRC Solr is the limiter, not your cores).

**Outputs:** `gene_capture_cells.csv`, `gene_capture_by_region.json` (macro/micro recall/precision/F1/Jaccard + recall-vs-identity curve + per-outcome denominators, both modes, with CIs). **Headline = macro-mean recall under `exclude_source`, per region, with CI and explicit denominators.**

---

## 8. Figure suite (the full slew)

All long-format CSVs consumed by standalone `scripts/render_region_<name>.py` (Agg backend, `REPO`-hop boilerplate, `~/Documents/py_venv/bin/python`, `region_validation/figures/`, `.png@300dpi` + `.pdf`). Shared module constants: `RANKS`, `REGION_ORDER` (by `amplicon_len_median`), `REGION_PALETTE` (consistent color per region everywhere), `DEPTH_RAMP`. **Every accuracy/gene-capture figure: exclude-self default, include-self as overlay/audit; Wilson CIs; Bacteria/Archaea split; tier facet where gene-capture.**

| # | Title | Type | x / y / encoding | Message |
|---|---|---|---|---|
| F1 | Accuracy heatmap region×rank | annotated heatmap (1×3: best/any20/abw) | x=rank (genus→phylum), y=region (short→long), cell=concordance, viridis 0–1 | headline grid of which (region,rank) transfer names; best-vs-any20 gap |
| F2 | Per-rank concordance by region | grouped lines + markers | x=rank, y=any20 concordance, hue=region | whole-curve vertical shift + slope per region |
| F3 | Deepest-correct-rank distribution | 100%-stacked horizontal bar (+ violin variant) | x=fraction, y=region, segments=deepest rank | most intuitive "how deep can this region resolve" |
| F4 | Gene-capture recall/precision/F1 | violins (`inner='quartile'`) + median strip, **tier-faceted** | x=region, y=metric, 3 panels | functional-capture distributions (bimodal cultured/MAG) |
| F4b | **% gene capture per region** | bar + Wilson CI | x=region, y=mean recall % | the literal "% gene capture by region" answer |
| F5 / F5b | Amplicon length vs accuracy / vs gene-capture | scatter+regression, point per region + faint per-amplicon underlay | x=median amplicon bp; y=genus concordance / mean F1 | does more sequence buy accuracy and/or gene content; where they diverge |
| F6 | Self-recovery rate by region | bar + Wilson CI, top-1 vs any-top-20 | x=region, y=fraction recovered | in-silico round-trip ceiling / region fidelity |
| F7 | Genus-failure confusion | confusion heatmap (truth×pred family, LogNorm), per region | x=pred family, y=true family, cell=count | random vs structured misassignment hotspots |
| F8 | Recall–precision scatter | scatter, point per amplicon, hue=region, size=tier, F1 iso-contours, per-region centroid diamond | x=precision, y=recall | each region's PR operating point + spread |
| F9 | Accuracy vs gene-capture | scatter + per-region regression, Spearman ρ+CI | x=**best_cosine/identity (continuous)**, y=gene-capture F1, hue=region | are taxonomically-resolvable amplicons also clean gene-capture (test decoupling) |
| F10 | Per-phylum accuracy | heatmap phylum×region (+ phylum color strip, shortened Proteobacteria legend), cells n<10 masked | x=region, y=phylum, cell=genus concordance | which lineages each region handles |
| F11 | Include-self vs exclude-self deltas | diverging/paired bar or slope chart | x=Δconcordance (incl−excl), y=rank faceted by region | leakage audit: robust generalization vs memorizing |
| F12 | Accuracy×gene-capture dashboard (capstone) | bubble scatter | x=genus concordance, y=gene-capture F1, size=self-recovery, color/label=region+len | one-figure decision aid; Pareto-dominant region |
| **F13** | **Amplification-rate / selection-bias** | grouped bar, per region × per phylum | x=region, y=extract_rate, facet/hue=phylum, **Bact vs Arch** | the "which taxa fail to amplify" answer; survivorship-bias guard |
| **F14** | **Tie/ambiguity cluster size** | violin or ECDF | x=region, y=#refs within ε of best identity | short-region ambiguity that inflates any-top-k |

**Precomputed tables (renderers read ONLY these):** `concordance_long.csv`, `asv_region_concordance.csv`, `gene_capture_long.csv`, `self_recovery.csv`, `confusion_genus.csv`, `amplification_by_phylum.csv` (new, F13), `ambiguity.csv` (new, F14), `region_meta.csv` (+ `domain`, version stamps). Each carries `domain` and (where gene-capture) `tier`.

---

## 9. Script architecture (ordered; front-loaded vs single parallel launch)

Existing modules reused **as-is** (no rewrite): `insilico_pcr.py` (panel mode + the §4 both-strand/edits upgrades), `align_hits.py`, `edlib_biopython_hits.py`, `gpu_align.py`, `select_references.py`, `prefetch_gene_families.py`, `build_synthetic_genomes.py` (+ §7 `pgfam_ids` upgrade), `render_diff_heatmaps.py` (convention donor).

| # | Script (to write, absolute under `/home/freiburger/Documents/EmilyKin/16S_validation/scripts/`) | Responsibility | In → Out | Phase |
|---|---|---|---|---|
| S1 | `select_10k.py` | §2 deterministic selection | `16S_md5_{ID,seq}.json`, taxdump → `bvbrc_16s_validation_10k.json`, `src16s_10k.json` | front-load |
| S2 | `build_self_maps.py` | parse `16S_metadata.json` once → per-operon maps | `16S_metadata.json` → `md5_to_genomes.json`, `genome_to_md5s.json` | front-load (once, ~25 min) |
| S3 | *(run `insilico_pcr.py --combined-fasta`)* + `build_provenance.py` | panel extraction (both-strand, edits) + provenance/fan-out | `src16s_10k.json`, panel → `combined_amplicons.fasta`, `.stats.json`, `amplicon_provenance.json`, `amplicon_md5_fanout.json` | front-load |
| S4 | `dedup_amplicons.py` | collapse identical inserts via fan-out | `combined_amplicons.fasta` → `combined_amplicons.unique.fasta` | front-load |
| **S5** | *(run `align_hits.py`)* | **THE single parallel/GPU launch** | `*.unique.fasta` → `hits_withself/asv_top20_alignment_hits.json` | **parallel (tmux)** |
| S6 | `select_excludeself.py` | §5.4 cascade: in-memory `exclude_self → select_representatives` for all 4 levels; audit | hits + prov + maps + fan-out → `selection_withself.json`, `selection_noself_<level>.json`, `.removed.json` | post (in-mem, fast) |
| S7 | `build_truth.py` | §6 `benchmark_truth.json` via `lineage_for` | prov + taxdump → `benchmark_truth.json` | front-load (can run anytime) |
| S8 | `prefetch_pgfams.py` (thin wrapper over `prefetch_gene_families.py`) | two-pass PGFam fetch (truth 10k + distinct reps), hyphenated host | selections → `bvbrc_cache/genome_gene_families.json` | post (tmux, ~1 hr) |
| S9 | `score_taxacc.py` | §6 metrics, both self-modes, domain/tier strata, Wilson CI, McNemar | hits + selections + truth → `taxacc_*` + `concordance_long.csv`, `asv_region_concordance.csv`, `self_recovery.csv`, `confusion_genus.csv`, `ambiguity.csv` | post |
| S10 | `score_genecap.py` | §7 `score_cell(T,P)` + edge logic, both modes, tier strata | selections + PGFam cache → `gene_capture_cells.csv`, `gene_capture_by_region.json`, `gene_capture_long.csv` | post |
| S11 | `build_region_meta.py` | `region_meta.csv` (+ domain, extract_rate, version stamps), `amplification_by_phylum.csv` | `*.stats.json` + selection + truth → tables | post |
| S12 | `render_region_*.py` ×14 | F1–F14, standalone, read only the flat tables | tables → `figures/*.png/.pdf` | render (seconds) |

**Front-loaded:** S1–S4, S7 (everything before alignment, plus truth). **Single parallel launch:** S5 (`align_hits.py`, the only GPU/60-worker job; S8 PGFam prefetch is a second, lighter tmux job). **Post (pure, fast, no GPU):** S6, S9–S12.

---

## 10. Cost / wall-time per stage + total; tractability verdict

Anchored to measured `run_stats.json` (12.27 edlib cpu-s/query, 59.5× at 60 workers) and `gpu_align_stats.json` (57.3 GCUPS). N = unique amplicons after dedup; raw is up to 8 regions × ~10k ≈ 80k, dedup typically → ~20–40k.

| Stage | Model | N=80k (no dedup) | N=40k | N≈20k (dedup) |
|---|---|---|---|---|
| S1 select_10k | one pass + taxopy | ~2–5 min | — | — |
| S2 build_self_maps | one-time `16S_metadata.json` parse | ~25 min (once) | — | — |
| S3 extraction + provenance | scan 10k×panel | ~10–20 min | — | — |
| S4 dedup | md5 of inserts | < 5 min | — | — |
| **S5 edlib prefilter** (CPU, dominant) | 12.27 cpu-s/q ÷ 60 | **~4.6 wall-hr** | ~2.3 hr | **~1.1 hr** |
| S5 GPU SW | N·500 pairs ÷ 57.3 GCUPS | ~7 min | ~3.5 min | ~2 min |
| S5 biopython rescore (gpu_cand=100) | 2.47 cpu-s/q ×100/500 ÷ 60 | ~11 min | ~5.5 min | ~3 min |
| S5 enrichment (taxopy) | cached | few min | few min | few min |
| **S5 align total** | | **≈ 5.0 wall-hr** | **≈ 2.6 hr** | **≈ 1.3 hr** |
| S6 exclude-self ×4 levels | in-memory re-selection | < 5 min total | — | — |
| S8 PGFam prefetch | 12 workers, ~4.5k gated, API-bound | ~50–70 min (once) | — | — |
| S9/S10 scoring | pure post | ~10–30 min | — | — |
| S12 render ×14 | flat tables | seconds each | — | — |

**Totals:** end-to-end **~6–7 wall-hr** at no-dedup N=80k, **~3–4 hr** after dedup (~20k), plus the one-time ~25-min metadata parse and ~1-hr PGFam prefetch (parallelizable with alignment). **Edlib-prefilter-bound, not GPU-bound** (GPU SW is ~2% of wall; GPU memory ~0.6 GB, a non-issue). The exhaustive no-prefilter GPU alternative is ~73 hr at N=80k — the top-500 prefilter is correctly ~15× better and stays.

**Verdict: TRACTABLE in a single run, entirely local** (RTX 5070 Ti + 64-core, no cloud/$). One `align_hits.py` process handles N=80k with no chunking (memory verified). **Dedup-first is the single highest-leverage move** (2–5× fewer unique inserts → 1.3 hr alignment). Chunk into 4×20k shards only if a single ~5-hr process is unacceptable for fault tolerance (re-pays ~30 s ref load per shard).

---

## 11. Risks and open decisions for the user

1. **[BLOCKER] Gene-ID sets not persisted.** Gene-capture recall/precision/F1/Jaccard and figures F4/F4b/F5b/F8/F9/F12 are **uncomputable** until `build_synthetic_genomes.py`'s gene-union builder persists `gene_union.pgfam_ids`. **Decide:** patch the builder now (recommended), or restrict v1 to taxonomic accuracy only.
2. **Archaea scope.** 6/8 core regions cannot amplify archaea. **Decide:** (a) formally scope the study to Bacteria (drop the 4,964 archaea from selection, simpler), or (b) keep archaea and add archaeal primer variants (515F-Arch/Arch806R) so V4/V4–V5 archaeal performance is reportable. Default: keep, stratify, document.
3. **FullLength16S exceeds `MAXQ=600`.** The ~1500 bp positive control will trip `assert maxq<=MAXQ` in `gpu_align.py`. **Decide:** bump `MAXQ` to 1600 (more GPU buffer, still <1 GB) or run the full-length control in a separate CPU/edlib pass. Pre-flight check required before S5.
4. **Self-exclusion strictness for the headline.** Four levels offered (`exact_self`→`near_identical`). **Decide:** confirm `near_identical` (τ=0.987) as the headline honest number, or pick τ=0.995. This materially changes the reported generalization accuracy.
5. **Top-20 vs top-500 substrate for exclude-self.** If self+sisters routinely fill the top-20, exclude-self empties hit lists. **Decide a threshold** on post-filter empty rate above which S6 switches to `asv_top500_alignment_hits_gpu.json` (1.1 GB).
6. **Insert-bound re-tuning.** Current `lo/hi` came from a 200-genome convenience sample. **Decide:** re-tune on a phylum-stratified ≥5k sample before the production run (cheap; prevents per-region clipping bias).
7. **`emax=3` leniency.** A 3+3-mismatch "amplification" counts as a hit. **Decide** the strict-mode policy (total-edit cap) — affects whether borderline regions look artificially universal.
8. **V6–V8 coverage.** The mislabeled entry is collapsed to V7–V9. **Decide** whether V6 genuinely needs coverage (add the optional 926F/1378R row) or whether 3'-coverage via V7–V9 suffices.
9. **Determinism of the embedding/GPU numbers.** Selection is bit-reproducible; the alignment numbers are not unless `torch.use_deterministic_algorithms(True)`, `PYTHONHASHSEED`, and deterministic tie-breaks are pinned. **Decide** whether full reproducibility of figures (not just selection) is required.
10. **Dedup vs provenance.** Dedup is the key speedup but requires the `amplicon_md5_fanout.json` fan-out to restore per-source provenance. Confirm the fan-out is built before dedup (S3 before S4).

**Headline deliverables:** per-region taxonomic accuracy (per-rank, exclude-self/`near_identical`, Wilson CI, domain/tier-stratified) and per-region % gene capture (macro-mean recall, exclude_source at organism level, with precision/F1/Jaccard and explicit denominators) — the two numbers that answer "accuracy and % gene capture by 16S region," with F11 quantifying how much each region's apparent accuracy depends on the source genome being in the reference set.