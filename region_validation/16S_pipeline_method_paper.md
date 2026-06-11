# Hardware-Adaptive 16S-to-Reference-Genome Mapping and Probabilistic Synthetic-Genome Reconstruction for Amplicon-Resolved Metabolic Modeling

## Abstract

16S rRNA amplicon surveys enumerate the members of a microbial community but reveal neither their genomes nor their metabolism, leaving a gap between taxonomic census and mechanistic, genome-scale metabolic modeling. We present a pipeline that closes this gap by mapping each amplicon sequence variant (ASV) to representative reference genome(s) in the BV-BRC 16S database (1,428,909 raw headers; 460,117 md5-deduplicated unique sequences; 459,301 retained at min_ref_len ≥ 200 bp) and reconstructing a probability-weighted synthetic genome suitable for protein-constrained flux balance analysis (prFBA). We first explored an embedding plus vector-database approach for approximate nearest-neighbor retrieval but abandoned it in favor of exact alignment, after a classical k-mer retriever beat the embedding on every closest-hit metric and the dense cosine space failed to resolve the equal-score tie clusters that determine the correct representative. The pipeline makes four contributions: (1) a hardware-adaptive aligner combining a hardware-independent edlib edit-distance prefilter with a runtime-compiled CuPy/NVRTC GPU Smith-Waterman engine (57.3 GCUPS, 7.8× over a 60-core CPU baseline, bit-exact to Biopython on 400/400 pairs); (2) decision-tracked, tie-aware representative selection over a top-500 shortlist; (3) probability-weighted synthetic-genome reconstruction; and (4) a deterministic 10,000-organism, 8-region, 67,991-amplicon in-silico validation reporting both self-inclusion (realistic) and self-exclusion (novel-organism generalization) taxonomic accuracy and PGFam gene-capture, exposing the leakage gap between the two. Across an eight-region panel and 67,991 amplicons from 10,000 source organisms, include-self genus accuracy of 0.94–0.99 collapses to 0.60–0.73 under self-exclusion, and gene capture (macro recall) falls from 95–99% to 68–72% — a near-uniform ~27-point leakage gap. Family-level assignment is far more robust to novelty (0.79–0.86 under self-exclusion), so family-level metabolism transfers even for organisms absent from the database, while genus-level transfer is where short amplicons pay the penalty — quantifying the cost of amplicon brevity for inferring the metabolism of novel organisms.

## 1. Introduction

High-throughput 16S rRNA gene amplicon sequencing is the dominant, cost-effective survey for microbial community composition. By clustering reads into amplicon sequence variants (ASVs) and assigning taxonomy, a single run can census the membership of an environment such as an enhanced biological phosphorus removal (EBPR) bioreactor. What such a survey cannot deliver is the metabolism of those members: a short hypervariable marker carries no information about the thousands of protein-coding genes, pathways, and stoichiometric constraints that a genome-scale metabolic model requires. To build a genome-scale, protein-constrained flux balance (prFBA) model directly from an amplicon survey, each ASV must therefore be connected to one or more whole genomes whose gene content can stand in for the un-sequenced organism. This connection is the central problem this paper addresses, and it decomposes into three coupled subproblems.

The first is **scalable, accurate alignment**. Each ASV must be compared against the BV-BRC 16S reference set — 1,428,909 raw header-to-sequence records, collapsed by md5 to 460,117 unique sequences, of which 459,301 survive a 200 bp minimum-length filter. Exhaustive Smith-Waterman scoring of every ASV against all 459,301 references is accurate but expensive, while fast heuristics risk missing the true representative. We initially pursued a sequence-embedding plus vector-database strategy, hoping approximate nearest-neighbor search would sidestep the cost of exact alignment; we abandoned it because nucleotide 16S embeddings could not reliably separate the correct genome from near-identical decoys, and because the biologically decisive signal lives in large clusters of references that tie at the top alignment score (empirically up to ~200 genomes). Recovering those clusters demands exact local alignment, not approximate retrieval.

The second is **principled multi-genome representative selection**. Because a marker as short as the V4 or V4–V5 region is frequently identical across many genomes, the correct answer is rarely a single best hit but a tie cluster that must be resolved by an auditable, deterministic rule rather than an arbitrary first hit. The third is **functional reconstruction**: converting the selected representative(s) into a single synthetic genome whose gene complement reflects the uncertainty in the mapping, so that downstream metabolic reconstruction inherits a calibrated, probability-weighted gene set rather than a hard guess.

This paper contributes a complete, hardware-adaptive solution to all three. (1) An aligner pairs an always-on edlib infix edit-distance **prefilter** — which discards ~99.9% of references before any costly scoring — with a Smith-Waterman scoring stage whose backend adapts to available hardware (CUDA GPU > Apple Metal > CPU). The GPU engine is a runtime-compiled CuPy/NVRTC kernel needing no CUDA toolkit, reaching 57.3 GCUPS, a 7.8× speedup over a 60-core CPU baseline, with results bit-exact to Biopython on 400/400 control pairs and 100% top-5 overlap with CPU-exhaustive search. (2) A decision-tracked selection procedure resolves tie clusters reproducibly over a **top-500 shortlist** and records why each representative was chosen. (3) A probability-weighted reconstruction propagates mapping uncertainty into the synthetic genome. (4) A deterministic 10,000-organism, 8-region, 67,991-amplicon in-silico benchmark evaluates the whole pipeline under both **self-inclusion** (the realistic experimental run, with the source organism present in the database) and **self-exclusion** (the source and its near-identical neighbors removed, measuring generalization to novel organisms), reporting per-rank taxonomic concordance and PGFam gene-capture for each.

## 2. Methods

### 2.1 Reference database

All amplicon sequence variants (ASVs) were mapped against a 16S rRNA reference set derived from the Bacterial and Viral Bioinformatics Resource Center (BV-BRC) [9]. The raw pull (`BV_BRC_16S.json`, 2.0 GB) comprised 1,428,909 header→sequence entries. Because BV-BRC stores one 16S record per annotated feature, the raw set is heavily redundant — multiple genomes, and multiple rRNA operons within a genome, contribute byte-identical sequences. We therefore deduplicated by MD5 digest of the uppercased sequence (`16S_md5_seq.json`, a `{md5: sequence}` map, 496 MB), collapsing the raw set to 460,117 unique sequences. Identical sequences (same MD5) are represented exactly once in the search space and resolved to a representative genome only at the reporting stage. Applying a minimum-length filter of ≥200 bp (`min_ref_len`, which removes fragmentary 16S entries) leaves **459,301 unique references**, the value used throughout.

Each unique sequence is linked back to its source via `16S_md5_ID.json` (81 MB, MD5→BV-BRC header). A header parser extracts the organism name, BV-BRC genome ID, taxon ID (the genome ID truncated at the first "."), and feature ID. Full NCBI lineages are then attached with **taxopy**, reading the `nodes.dmp`/`names.dmp` taxdump shipped with the prFBA repository, so every reference carries a complete species-through-domain lineage for the taxonomy-aware selection logic described in §2.3.

### 2.2 Alignment

Mapping each EBPR V4–V5 ASV to 459,301 references by full Smith-Waterman (SW) is wasteful, because the overwhelming majority of references are obviously irrelevant. We therefore used a three-stage strategy in which a fast, hardware-independent edit-distance prefilter eliminates ~99.9% of candidates before any SW work, an accelerated SW engine scores the surviving shortlist, and an exact CPU re-alignment recovers the per-hit statistics that a score-only kernel cannot provide.

**Stage 1 — edlib edit-distance prefilter (top-500).** For every ASV we scan all unique references with `edlib.align(query, ref, mode="HW", task="distance", k=kthr)` [5]. The "HW" (infix) mode computes the edit distance of the best placement of the ASV *within* each reference, which is the correct semantics for a short amplicon against a longer 16S gene. To avoid paying full edit-distance cost on hopeless candidates, an adaptive *k*-band is used: a size-`prefilter_k` max-heap of `(−distance, ref_idx)` is maintained, and once full, `kthr` is set to the current worst-kept distance, so edlib returns −1 (early reject) for any reference exceeding the band. The scan returns the top-500 references ascending by edit distance. This stage runs identically on every machine.

The prefilter depth of `prefilter_k = 500` is the central tuning choice. Empirically, the named representative of an ASV is frequently embedded in a large cluster of references that tie at the top SW score (up to ~200 genomes observed at a single top score, e.g. multi-strain species with identical 16S). A shallow prefilter would truncate such tie clusters arbitrarily and could exclude the correct representative; 500 leaves a wide margin above any observed cluster while still scoring <0.11% of the database.

**Stage 2 — Smith-Waterman scoring of the 500-candidate shortlist.** Only the 500 prefiltered candidates per ASV are scored by local, affine-gap SW. The scoring scheme is a BLASTN-like nucleotide scheme — **match +2, mismatch −3, gap-open −5, gap-extend −2** — defined once and shared by both backends so scores are directly comparable. On a CUDA host this stage runs on the GPU (CuPy SW kernel); otherwise it runs on the CPU (Biopython `PairwiseAligner.score`, `mode="local"`) [6].

**Stage 3 — Biopython re-alignment (top-20).** A score-only engine yields the maximum SW score but *not* percent identity, coverage, or reference coordinates — none of these are derivable from the score alone. We therefore re-align the best-scoring candidates with Biopython, counting matching non-gap columns over aligned columns for identity and reading `aln.aligned[1]` for the reference alignment start/end. On the GPU path the top `gpu_cand = 100` scored candidates are re-aligned (the keep count is `max(gpu_cand, k2)`; the margin above `k2` absorbs equal-score tie clusters), then trimmed to the final **top-20** (`k2`) hits per ASV with full identity and coordinates. The top-20 hit list (`asv_top20_alignment_hits.json`) is the input to reference selection (§2.3). The remaining parameters are summarized in Table 1.

**Table 1. Alignment-pipeline parameters.**

| Parameter | Default | Rationale |
|---|---|---|
| `prefilter-k` (K1) | 500 | Captures large equal-SW-score tie clusters (up to ~200) without scoring all references |
| `k2` | 20 | Final hits emitted per ASV with full identity/coordinates |
| `gpu-cand` | 100 | GPU candidates Biopython-rescored; margin above `k2` to absorb ties |
| `min-ref-len` | 200 bp | Drops fragmentary 16S references |
| `workers` | `max(1, ncpu − 4)` | Leaves 4 logical cores free (reference host: 60 workers) |
| MAXQ | 1600 | Max query length for per-thread int16 DP buffers (fits the ~1500 bp full-length positive control; longest hypervariable ASV = 561 bp) |
| SW pairs batch | 4,000,000 | GPU pair-scoring batch size |

**GPU kernel design.** Off-the-shelf GPU aligners were infeasible here for two compounding reasons: the only widely available GPU SW engines (MMseqs2-GPU, CUDASW++4.0) operate on protein/PSSM data and cannot accelerate a nucleotide task (16S rRNA is non-coding, so protein translation is invalid), and the alternatives (GASAL2, ADEPT, WFA-GPU, NVBIO, Parabricks) require `nvcc` or a container that the no-toolkit target environment does not provide. We therefore wrote a custom kernel whose C++ source is compiled **at runtime via CuPy `RawModule`/NVRTC** (`--std=c++14`), with no `nvcc` and no CUDA toolkit; NVRTC ships inside the `cupy-cuda13x` wheel and targets the live GPU's `sm_120` (Blackwell; NVIDIA RTX 5070 Ti, 16 GB, driver 595.71.05 / CUDA 13.2). This "no-nvcc plus brand-new sm_120" combination is precisely what disqualified the off-the-shelf tools.

The kernel assigns **one thread per (query, reference) pair** and runs affine-gap (Gotoh) local SW; two `__global__` entry points exist — `sw_cross` (full cross-product, exhaustive benchmarking) and `sw_pairs` (arbitrary pairs, used for shortlist rescoring) — with a thread-block size of 128. Each thread holds two int16 rolling DP buffers (`short Hp[MAXQ+1]; short F[MAXQ+1]`, ≈6.4 KB/thread); the short ASV is the inner DP dimension and the long reference the outer loop. Because scores stay well below the int16 range for sequences ≤561 bp, int16 halves local-memory traffic for roughly 2× throughput. Base equality uses a 256-entry lookup table (A/C/G/T→0..3) so the GPU reproduces Biopython's character-equality on uppercased sequences exactly. The kernel is **score-only** (it returns only the maximum H value and stores no traceback), which is why the Biopython re-alignment of Stage 3 is required to recover identity, coverage, and coordinates.

**Performance.** Against a Biopython reference, the kernel was bit-exact on 400/400 random real ASV×reference pairs (`max_abs_diff = 0`). On an exhaustive benchmark of 20 validation ASVs × 459,301 references (3.78 trillion DP cells), the GPU completed in **65.92 s at 57.3 GCUPS**, versus 511 s for the 60-core CPU reference (**7.8× speedup**); GCUPS is computed as `(Σ ASV lengths × Σ reference lengths) / GPU_wall / 1e9`. Concordance with the CPU-exhaustive baseline was 100% top-5 MD5 overlap (100/100) and 20/20 top-1 score matches. We also ran the full prefilter-free search — 3,950 ASVs × 459,301 references, retaining the top 500 per ASV — to completion in 12,985 s (≈3.6 h) at 57.3 GCUPS, versus an extrapolated 28.0 h on CPU-exhaustive. In production, the prefiltered GPU path scores only the 500-candidate shortlist (~2M pairs across ~3,950 ASVs) rather than the full cross-product, so the operational cost is far below even the 65.92 s exhaustive figure.

### 2.3 Reference selection

The top-20 hits per ASV are reduced to a representative handful of genomes by `select_references.py`, which emits a JSON audit justifying every inclusion and exclusion. The selected `genome_id`s are exactly what feed the synthetic-genome merge (§2.4). Selection is a joint function of alignment quality *and* contamination control, executed as an **ordered 8-stage pipeline whose order is the load-bearing design choice**: each stage removes a class of artifact that, if left in place, would corrupt a later stage (Table 2).

**Table 2. The ordered selection pipeline.**

| Stage | Function |
|---|---|
| 0 | Reliability gate: drop hits with coverage `< cov_min` or `\|SW−edlib identity\| > swed_max` |
| 1 | Family floor + tier: classify the ASV by best identity; abstain if below family |
| 2 | Identity band: keep only hits within `band` identity of the best (on identity, not score) |
| 3 | Family-consensus + blank-family rule: drop taxonomically inconsistent off-targets |
| 4 | genome_id dedup: collapse multi-copy 16S operons to the best copy |
| 5 | Species dedup: ≤1 representative per species |
| 6 | Identity-decay representativeness rank |
| 7 | Marginal-gain greedy inclusion + contamination budget |
| 8 | Provenance tagging (anchor/congener, core/accessory, confidence) |

In code, the ASV-level "below family" abstain (stage 1) is tested before the per-hit reliability gate (stage 0), short-circuiting hopeless ASVs; the remaining stages proceed in numeric order. The ordering rationale is explicit and verifiable: the **reliability gate runs first** so that no partial or method-discordant (possibly chimeric) alignment can ever cast a consensus vote or become a representative; the **family-consensus filter is placed before dedup** so a "silent off-target" cannot survive dedup and be chosen as a representative; the **identity band is computed on identity, not score**, so the tie cloud reflects genuine sequence similarity rather than score artifacts; and **representativeness ranking precedes greedy inclusion** so the anchor is the single best representative and each added congener is the next-best genome that still contributes novelty.

**Identity tiers.** Each ASV is classified by its best-hit identity against thresholds `t_species = 0.987`, `t_genus = 0.945`, `t_family = 0.865` (Yarza/Stackebrandt bands [1,2]). Below the family floor (0.865) the ASV abstains entirely — "too divergent for a genome-level model (would inject false metabolism)" — yielding a taxonomy-only record and no genome.

**Reliability and consensus knobs.** Stage 0 uses `cov_min = 0.90` (minimum aligned-length / ASV-length) and `swed_max = 0.02`. Stage 2 uses `band = 0.005`. Stage 3 enforces consensus at the **Family** rank (genus is blank for ~59% of hits, family is stable): an off-consensus hit may override only if it is ≥`override_margin = 0.03` more identical *and* present in MiDAS, and a blank/unclassified-family hit is retained only if within `blank_family_tol = 0.001` identity of a named consensus-consistent hit (otherwise excluded as a probable silent off-target).

**Representativeness (stage 6).** Survivors are ranked by `rep_score = 0.55·identity + 0.30·(align_score/(2·asv_len)) + 0.15·coverage`, multiplied by `uncult_factor = 0.98` for uncultured organisms (a gentle nudge, not a gate), and then weighted by an identity decay `exp(−(best_id − identity)/decay_tau)` with `decay_tau = 0.01`.

**Marginal-gain inclusion and contamination budget (stage 7).** Candidates are added greedily in rank order while each new genome contributes ≥`novel_min = 0.15` *new* genes relative to the running union, and while the accumulated off-target gene mass stays within a tier-scaled contamination budget. The marginal-gain stop is the primary, data-driven cap; confidence-scaled per-tier ceilings act as secondary bounds (`cap_species = 2`, `cap_genus = 5`, `cap_family = 2`; absolute `hard_cap = 6`), so a less-certain (e.g. family-level) ASV cannot accrue many genomes' worth of speculative metabolism. At most `max_low_conf = 1` uncultured/MAG genome is admitted. The contamination test excludes a candidate when `contam_genes + (1−identity)·novel_genes > tier_budget · anchor_genes`, with tier budgets `contam_species = 0.25`, `contam_genus = 0.40`, `contam_family = 0.20`.

**Exact vs estimated novelty.** When a gene provider is supplied, "novelty" is computed from **exact PGFam gene sets**: for each candidate genome the running union of PGFam IDs is accumulated and `novel_frac = |gset − union| / |gset|` — a one-directional, Jaccard-style novelty of the candidate's own gene set against the accumulated union — and records are tagged `gene_data = "bvbrc"`. PGFam (global protein families) makes "same gene" comparable across genomes; the provider fetches each genome's CDS PGFam IDs from the BV-BRC API (`eq(genome_id,GID)&eq(feature_type,CDS)&select(pgfam_id)`) and caches them to disk so the ~21.7k-genome corpus is fetched at most once. A companion prefetch utility (`prefetch_gene_families.py`) populates this cache for *exactly* the ~4.5k genomes that reach the gene-dependent stage-7 logic (stages 0–6 are gene-independent), a ~79% saving over fetching all top-20 genomes. When no provider is available, a taxonomy-based estimator is used (shared-fraction by rank: `ov_species = 0.90`, `ov_genus = 0.75`, `ov_family = 0.55`, `ov_other = 0.45`, with accumulation decay `accum_decay = 0.12` and a `default_gene_count = 3000`); every such record is tagged `gene_data = "estimated"` so provenance stays honest, and downstream genome building refuses estimated selections unless explicitly permitted.

**Provenance (stage 8).** The first included genome is tagged `anchor`, subsequent ones `congener`. Confidence is `high` for a genuine binomial species name with consensus/MiDAS-consistent family, `medium` for a family-consistent but unnamed species, and `low` for uncultured or blank-family genomes. The output `_meta` records the SW scheme, the full knob set, tier counts, and the handful-size distribution.

### 2.4 Synthetic-genome reconstruction

For each ASV with at least one selected genome, the selected references are merged into a single **synthetic KBase Genome object** [8] that serves as the input to prFBA. ASVs with no selected genome receive a taxonomy-only manifest entry and no genome (abstain).

**Function-union merge.** The merge keys on gene **function**, not sequence. Iterating over each selected genome's features, the first occurrence of a function string creates a new merged feature; later genomes carrying the same function only increment a counter, with a per-genome guard so each genome contributes at most +1 per function. The synthetic genome therefore holds exactly one feature per distinct function across the selected set — the gene union by function.

**Gene→reaction ontology hook.** Each merged feature copies its source `aliases` and `ontology_terms`, and crucially ensures `ontology_terms['RAST']` contains the function string, because the downstream MSBuilder uses `ontology_terms['RAST']` to map genes to reactions. Each feature also receives a parallel CDS entry and a synthetic single-feature contig (`location = [[contig, 1, "+", dna_length]]`), with contig counts, `dna_size`, and the genome MD5 accumulated across features.

**Per-function probability (the prFBA weight).** The counter that accumulates how many selected genomes carry each function is, after the merge, normalized over all source genomes:

`probability = (number of selected genomes carrying the function) / (number of selected genomes)`

This value is written onto both the merged feature and its CDS, and is the per-function (and hence per-reaction) weight consumed by prFBA — functions present in all selected genomes carry weight 1, accessory functions carry a low fractional weight.

**KBase scaffold.** The template scaffold sets `source = 'Synthetic'`, `domain = 'Bacteria'` (overridden from consensus taxonomy), `genetic_code = 11`, `molecule_type = 'DNA'`, `genome_tiers = ['User']`, and `gc_content` = mean of the source GC contents, with feature counts (CDS, gene, protein-encoding) equal to the union size.

**Core/accessory provenance.** A per-function prevalence count across the ASV's selected genomes tags each feature as `core` (present in all selected genomes), `accessory` (unique to one genome, when more than one was selected), `shared` (intermediate), or `unknown`; accessory genes are retained as gap-fill candidates for the FBA reconstruction rather than as hard reactions. A provenance sidecar records the core functions (majority prevalence), accessory functions, and the full prevalence dictionary.

**Fast vs full builder.** The default `fast` builder issues two BV-BRC API calls per genome (one genome-core query, one `genome_feature` CDS query selecting `patric_id, product, figfam_id, pgfam_id, plfam_id, gene`), emitting per-CDS functions and aliases with empty DNA/protein sequences. Because the merge keys only on function, this reproduces the full gene union at ~0.5 s/genome (cached) and matches the shape of the prior sequence-free `genome_objects` run. The `full` builder additionally fetches contig and per-feature DNA/protein sequences for downstream KBase modeling (including a patch for BV-BRC now returning the `go` field as a list). Built genomes are cached to `bvbrc_cache/kbase_genome_cache/{gid}.json`. By default the pipeline is exact: it prefetches candidate gene sets, runs selection with exact PGFam (Jaccard) novelty, and writes `gene_data = "bvbrc"`; building from an estimated selection is refused unless explicitly allowed.

Both production runs used the fast builder with exact (`bvbrc`) gene data and zero failures, selecting from the GPU-alignment-based selection (`asv_reference_selection_gpu.json`): the EBPR (EmilyKin) dataset built **3,724 of 3,950** ASV genomes (226 abstained), and the codiffusion-bioreactor dataset built **3,219 of 3,276** (57 abstained).

### 2.5 Hardware-adaptive scoring

The pipeline detects its compute backend at runtime and adapts only the SW-scoring stage; the edlib prefilter is always identical. A hardware probe queries CUDA availability (CuPy device count > 0), Apple Metal/MPS availability (`torch.backends.mps.is_available()` on Darwin), logical core count, and free disk. A pure mapping then selects the backend in strict priority **CUDA > Metal > CPU**: a CUDA host scores the 500-candidate shortlist on the NVRTC SW kernel (§2.2); otherwise the CPU path scores the shortlist with Biopython's local affine-gap aligner using the identical match +2 / mismatch −3 / gap-open −5 / gap-extend −2 scheme, so results are backend-independent. Metal detection and strategy are in place for a future port, but because no Metal SW kernel is yet implemented, a detected Metal backend currently falls back to CPU scoring. In all cases Stage 1 (edlib top-500) and Stage 3 (Biopython identity/coordinate recovery) are unchanged; only the Stage-2 scoring engine differs, which is what makes the GPU and CPU pipelines produce concordant top hits (100% top-5 MD5 overlap, 20/20 top-1 score matches) while differing only in wall-clock time.

## 3. Validation and Development

The mapping pipeline that assigns each ASV to a representative BV-BRC genome was not designed in one pass. Each component — the alignment scorer, the prefilter that feeds it, the retrieval model itself, and the downstream representative-selection logic — was adopted only after a head-to-head validation against an exact ground truth retired a simpler or more fashionable alternative. This section reports those experiments in the order they shaped the production architecture: what was tested, what the result was, and how that result guided the design. It closes with the in-silico validation framework built to characterize the finished pipeline across 16S regions, for which a 500-genome pilot is reported as preliminary and a 10,000-organism run is pending.

### 3.1 GPU vs CPU Smith-Waterman scoring

The scoring stage computes a local SW alignment (match +2, mismatch −3, gap-open −5, gap-extend −2) between each ASV and its candidate references, emitting an exact percent-identity. Because percent-identity is the load-bearing quantity for every downstream taxonomic and gene-capture decision, the scorer had to be both correct and fast enough to consider a large candidate pool. We validated the CUDA kernel against the Biopython CPU oracle.

**Result — bit-exact equivalence.** Across 400 random real ASV×reference pairs the GPU and CPU scores were identical at every position: 400/400 exact, maximum absolute difference 0. The packaged test suite extends this to integer-level equality (`np.array_equal`) across nine scoring schemes — DNA local/global/semiglobal variants, EDNAFULL, BLOSUM62 — plus edge cases and int32-promotion cases. The GPU kernel is therefore not an approximation of the CPU scorer; it is the same scorer on different hardware.

**Result — throughput and timing.** An exhaustive 20-ASV × 459,301-reference sweep completed in 65.9 s at 57.3 GCUPS (3.78 trillion cells) on an RTX 5070 Ti, a 7.8× speedup over the 60-core CPU-exhaustive run (511 s). Critically, the full prefilter-free run was actually carried out: 3,950 ASVs × 459,301 references, retaining the top 500 per ASV, finished in 12,985 s (≈3.6 h) at 57.3 GCUPS, versus an extrapolated 28.0 h for CPU-exhaustive.

**How it guided development.** Bit-exactness meant the GPU backend could replace the CPU scorer with zero change in scientific output, removing the usual "is the fast path also the correct path?" caveat. The 3.6-hour exhaustive figure proved that even a brute-force whole-database scoring of every ASV was now tractable, which in turn meant the GPU could comfortably rescore a 500-candidate shortlist for the entire ASV set in a small fraction of that budget. The GPU SW kernel was therefore made the default scorer of the top-500 shortlist, with the bit-identical CPU path retained as a fallback.

### 3.2 edlib edit-distance prefilter vs full alignment scoring

Running SW against all 459,301 references for every ASV is exact but wasteful: the overwhelming majority of references are unrelated. We therefore prefilter with edlib — a fast edit-distance computation — keeping the top-500 candidates, then rescore those with SW. The validation question is whether the edit-distance prefilter ever discards a candidate that SW would have ranked highly, i.e. whether the prefilter loses recall.

**Result — prefilter recall.** Against the full 3,950-ASV GPU-exhaustive ground truth, the edlib-prefilter pipeline recovered the true top-1 hit for 99.47% of ASVs (the prefilter top-1 equaled the exhaustive top-1 for 99.14%), with a mean top-20 overlap of 97.6% and 90.9% of ASVs identical 20-for-20. The complementary failure mode is real but small: 219 of 3,950 ASVs (5.54%) lost at least one genuine top-20 hit, totaling 1,696 missed hits. An older 1,000-ASV benchmark against the Biopython oracle corroborates the picture — Spearman(−edit distance, %identity) = 0.994, top-1 agreement 91.1%, the Biopython #1 inside the edlib top-20 for 99.9% of ASVs, and the true #1 never falling outside the edlib top-50.

**How it guided development.** Two extremes were rejected by these numbers. Edit-distance alone is insufficient: a top-1 agreement near 91% and the existence of deep-hit losses mean edit distance is a *near*-surrogate for identity, not a replacement for it. Exhaustive SW is unnecessary: at ~99% top-1 recall the prefilter almost never costs the right answer. The design that survives is the middle path actually shipped — a generous **top-500 edlib prefilter feeding exact SW rescoring**. The depth of 500 (rather than, say, 20 or 50) is what keeps the ~5.5% deep-hit loss from touching the representative selection, as the concordance analysis in §3.4 confirms; for the 219 affected ASVs the computable best-score gap has a median of 0, i.e. the missed hits are score-tied with hits already retained.

### 3.3 The abandoned vector-database / embedding approach

Before the alignment pipeline, we built a dense-retrieval system: every BV-BRC 16S sequence (460,116 unique) was embedded with `nucleotide-transformer-v2-500m-multi-species` (NT-v2-500m, 498 M parameters, non-overlapping 6-mer tokens) [7] into an L2-normalized 1024-dimensional vector, stored as a (460,116 × 1024) float16 matrix, with cosine similarity as the retriever. Because the experimental amplicons are ~373 bp V4–V5 fragments rather than full 16S genes, we added an in-silico-PCR *region-matching* stage that excised the V4–V5 insert from every reference (1,089,954 inserts → 97,624 unique) and re-embedded them, lifting mean cosine from 0.9912 to 0.9974 and reaching a true cosine of 1.0 for 803 ASVs.

**Result — classical retrievers match or beat the embedding.** Benchmarked against the Biopython %identity ground truth over the full 97,624-insert database (n=200 ASVs), a no-ML, no-GPU **k-mer-8 TF-IDF** retriever decisively outperformed the NT-v2 cosine retriever on every closest-hit metric: exact top-1 75.5% vs 50.5%, recall@100 of the true #1 99.0% vs 84.5%, and all-true-top-5-present@100 96% vs 37%. The mechanism is interpretable — one substitution flips up to *k* k-mers, so k-mer overlap tracks nucleotide identity, whereas mean-pooling averages that signal away. Reciprocal-rank fusion of cosine and k-mer added nothing (it dropped top-1 to 57%), and a hierarchical embedding tree was catastrophically worse than flat search (greedy 2-level reached the true cluster only 8.7% of the time). Two structural problems compounded this: the 16S space is so concentrated that random pairs average cosine 0.90 — a hit at 0.92 "is essentially background" — and cosine correlated with true %identity at only Spearman 0.55 within the candidate cloud.

**How it guided development.** Alignment, unlike cosine, *is* the definition of percent-identity, so it reproduces the ground truth by construction and yields interpretable literature identity bands (≈98.7% species, 94.5% genus, 86.5% family [1,2]) rather than a space-specific cosine scale that must be hand-calibrated against a 0.90 background. Production therefore reverted to exact alignment — edlib prefilter → SW rescoring — for fine identity mapping, at ~279× lower per-pair cost than the Biopython oracle it reproduces. The embedding was not discarded everywhere: it remained a competent coarse *locator* of the community (median best cosine 0.9985; 99.88% of ASVs ≥0.97; only 3 ASVs novel at <0.90), but it was retired from the fine-ranking role that determines representative genomes.

### 3.4 edlib-vs-GPU representative concordance

The two preceding sections validate alignment *scores*; this section validates the *decision* those scores feed. We ran the identical representative-selection algorithm and identical BV-BRC PGFam provider on two different upstream alignments — the GPU-exhaustive top-20 versus the edlib-top-500→Biopython pipeline — and asked whether the chosen representative genome set agreed.

**Result — near-total concordance.** On the EmilyKin set (3,950 ASVs, 3,724 co-built) the representative genome set was identical for 99.97% of ASVs (1 differing), the anchor genome matched for 99.97%, and the genome-set Jaccard had mean and median 1.0; all 5,549 GPU representative genomes were present in the edlib top-20. On the new_data set (3,276 ASVs, 3,219 co-built) the genome set was identical for 99.94% (2 differing), with the anchor genome matching for 100.0%. Every divergence resolved to an equal-SW-score tie, not a recall failure: both new_data "misses" were already inside the K1=500 prefilter and fell only at the k2=20 cut inside a cluster of SW-score-identical genomes (one had 26 tied members at ranks 4–29; the other ≥199 tied members filling ranks 2–200+), and edlib's alternative choices carried the same SW score. A top-100 recovery analysis confirmed the first would be recovered at depth 100, while the second — sitting in a ~200-member tie — needs depth ~1000 for a guaranteed tie-break.

**How it guided development.** This was the decisive evidence that the prefilter-plus-selection stack is robust: switching the entire upstream scorer changes the final answer for roughly 1 ASV in 1,000, and those cases are interchangeable tie members rather than errors. It also justified retaining the **top-500 prefilter depth specifically to span tie clusters** — a shallower cut would let large equal-score tie groups (which can exceed 100 members) shuffle the selected representative, whereas depth 500 keeps the selection deterministic and reproducible across backends.

### 3.5 The 10,000-organism in-silico validation

The preceding experiments validate the pipeline's internal consistency; the final framework validates its *biological* accuracy. We treat 10,000 BV-BRC genomes as synthetic query organisms, perform in-silico PCR to generate amplicons for eight 16S regions, run each amplicon through the unaltered production pipeline, and score the result against the organism's true NCBI lineage and PGFam gene set.

**Why these 10,000 organisms.** The goal is to probe the whole tree of life without redundancy bias. Selection is deterministic and taxonomy-stratified (`SEED = 1729`, byte-identical output), drawing ≤1 genome per species down a `domain→phylum→…→species` tree. Over-redundant clades are down-weighted by sqrt-flattening (α = 0.5), which pulls the top phyla from raw 37/20/14% toward ~13–16% each, and no clade may exceed 35% of its parent's budget (`CAP_FRAC = 0.35`), so dominant taxa cannot swamp the draw. A full-length ≥1400 bp 16S is required *before* allocation, and archaea are kept and stratified by domain so that archaeal performance in the two archaea-competent regions is measurable.

**What it examines.** Each amplicon is scored two ways from a single unaltered pipeline run. **Self-inclusion** (primary, realistic) leaves the source organism in the database — the experimental upper bound. **Self-exclusion** (generalization) post-hoc drops the source genome, its taxon, region-identical neighbors, and near-identical hits (headline τ = 0.987 species threshold), re-runs pure representative selection on the filtered record, and re-scores — performance as if the organism were novel. An anti-inflation rule scores exclude-self correctness against the include-self denominator. The framework reports per-region, per-rank taxonomic accuracy (Wilson 95% CIs; paired McNemar across regions) and PGFam gene-capture recall/precision/F1/Jaccard, across the eight-region panel (V1-V2, V1-V3, V3-V4, V4, V4-V5, V4-V5_944R, V7-V9, and FullLength16S as positive control). In-silico PCR (`insilico_pcr.py`, `emax=3`, primers included) generated 67,991 amplicons across the 10,000 sources.

**How it guided development.** Building this framework hardened several design choices: the **eight-region panel** (with V4–V5 as the headline anchor and FullLength16S as positive control) was fixed here, including the relabeling of a mislabeled "V6–V8" region into V7–V9; the **self-inclusion vs self-exclusion framing** was adopted to separate the realistic experimental number from honest novel-organism generalization; and a **prefetch-before-select** hardening was added so that representative selection runs against a complete candidate record.

#### Pilot (500-genome) results — preliminary

A 500-genome pilot extracted 2,501 amplicons. Values are anchor-predictor, summed across Bacteria+Archaea, for include-self vs exclude-self (τ=0.987). **These results are preliminary and will be superseded by the full run.**

**Table 3. Pilot genus accuracy by region (anchor).**

| Region | Incl. acc | Excl_0.987 acc |
|---|---|---|
| FullLength16S | 1.0000 | 0.8605 |
| V1-V2 | 0.9786 | 0.7647 |
| V1-V3 | 0.9930 | 0.7959 |
| V3-V4 | 0.9867 | 0.7215 |
| V4 | 0.9593 | 0.6612 |
| V4-V5 | 0.9713 | 0.6860 |
| V4-V5_944R | 0.9786 | 0.7294 |
| V7-V9 | 0.9827 | 0.7405 |

**Table 4. Pilot gene-capture macro-mean recall by region (= % gene capture).**

| Region | Incl. recall | Excl_0.987 recall |
|---|---|---|
| FullLength16S | 0.9982 (99.82%) | 0.6819 (68.19%) |
| V1-V2 | 0.9821 (98.21%) | 0.6680 (66.80%) |
| V1-V3 | 0.9890 (98.90%) | 0.6699 (66.99%) |
| V3-V4 | 0.9705 (97.05%) | 0.6330 (63.30%) |
| V4 | 0.9631 (96.31%) | 0.6335 (63.35%) |
| V4-V5 | 0.9709 (97.09%) | 0.6377 (63.77%) |
| V4-V5_944R | 0.9684 (96.84%) | 0.6487 (64.87%) |
| V7-V9 | 0.9720 (97.20%) | 0.6734 (67.34%) |

In the pilot, include-self genus accuracy is ~0.96–1.00 across all regions and collapses to ~0.66–0.86 under exclude-self, with longer regions generalizing better (FullLength16S 0.86, V7-V9 0.74) than short V4 (0.66). Gene capture mirrors this: ~96–100% include falls to ~63–68% exclude, FullLength16S highest and V3-V4/V4 lowest — the expected leakage/generalization gap.

#### Full 10,000-organism validation — results

The full validation ran 10,000 deterministically stratified source genomes through in-silico PCR to **67,991 amplicons** across the eight regions (8 regions × ~8,500 amplicons each), every amplicon passed through the unaltered production pipeline and scored against its source organism's true NCBI lineage and PGFam gene set. All 67,991 amplicons were scored (0 unscoreable, 0 without a selection record). Taxonomy was resolved for every source; the include-self gene-capture pass scored 66,978 cells (958 sources had an empty BV-BRC PGFam set and were marked `truth_empty`), and the self-exclusion pass scored 50,495 cells, with 16,538 amplicons (24.3%) abstaining after their self / near-identical hits were removed — the honest novel-organism denominator. Values below are the anchor predictor; accuracy is pooled across Bacteria and Archaea, and gene capture is macro-mean recall (= % of the source's PGFam genes recovered by the selected representatives).

**Table 5. Full 10,000-organism validation results (anchor; accuracy pooled Bacteria+Archaea; gene capture = macro recall).**

| Region | Genus acc (incl) | Genus acc (excl) | Gene capture (incl) | Gene capture (excl) |
|---|---|---|---|---|
| FullLength16S | 0.986 | 0.730 | 99.4% | 72.2% |
| V1-V2 | 0.967 | 0.720 | 97.5% | 70.7% |
| V1-V3 | 0.973 | 0.717 | 98.1% | 71.2% |
| V3-V4 | 0.953 | 0.639 | 96.5% | 69.0% |
| V4 | 0.936 | 0.599 | 95.1% | 67.9% |
| V4-V5 | 0.950 | 0.626 | 96.0% | 68.6% |
| V4-V5_944R | 0.951 | 0.631 | 96.0% | 69.1% |
| V7-V9 | 0.959 | 0.653 | 96.8% | 69.6% |

The full run confirms and sharpens the pilot. Under realistic **self-inclusion**, genus accuracy is 0.94–0.99 and gene capture 95–99% across all eight regions — when the source organism is present in the database, even a short amplicon recovers nearly all of its metabolism. Under honest **self-exclusion**, genus accuracy falls to 0.60–0.73 and gene capture to 68–72%, with the now-familiar **region-length gradient**: full-length 16S and the long V1–V3 generalize best (genus 0.73/0.72, capture 72.2%/71.2%), while the short V4 generalizes worst (genus 0.60, capture 67.9%); V4–V5 and V7–V9 sit between. Notably, the gene-capture **leakage gap is near-uniform at ~27 points** across every region (0.268–0.275), whereas the genus-accuracy gap widens for the shorter amplicons — i.e. read length costs taxonomic resolution more than it costs raw functional recovery.

A per-rank view exposes a second, practically important effect: **family-level assignment is far more robust to novelty than genus-level**. Pooling Bacteria+Archaea at the Family rank, exclude-self accuracy holds at 0.79–0.86 — 12–19 points above the genus figure for the same amplicon (Table 6).

**Table 6. Family- vs genus-rank accuracy under self-exclusion (anchor; pooled Bacteria+Archaea).**

| Region | Family acc (incl) | Family acc (excl) | Genus acc (excl) | Δ (family − genus, excl) |
|---|---|---|---|---|
| FullLength16S | 0.987 | 0.858 | 0.730 | +0.128 |
| V1-V2 | 0.974 | 0.830 | 0.720 | +0.110 |
| V1-V3 | 0.977 | 0.840 | 0.717 | +0.123 |
| V3-V4 | 0.969 | 0.807 | 0.639 | +0.168 |
| V4 | 0.964 | 0.790 | 0.599 | +0.191 |
| V4-V5 | 0.968 | 0.815 | 0.626 | +0.189 |
| V4-V5_944R | 0.968 | 0.820 | 0.631 | +0.189 |
| V7-V9 | 0.971 | 0.819 | 0.653 | +0.193 |

For the short regions the family→genus gap is ~0.19 — roughly **one in five novel short-amplicon organisms lands in the correct family but is assigned a wrong-genus relative** — shrinking to ~0.12 for full-length/V1–V3, where longer reads more often resolve the genus itself. The practical reading is that **family-level metabolism transfers robustly across all regions even for organisms absent from the database**, while genus-level transfer is where short amplicons pay their penalty.

## 4. Discussion

Three design choices distinguish this pipeline and its validation from embedding-based alternatives. First, **exact alignment with decision-tracked reference selection** yields interpretable, threshold-grounded genome mapping that opaque dense embeddings cannot. A nucleotide-transformer encoder [7] embeds every 16S sequence into a 1024-dim cosine space, but that space is globally so concentrated — random sequence pairs average cosine 0.90 — that absolute similarities are uninformative, and cosine ordering within a near-identical V4–V5 cloud correlates only ~0.55 (Spearman) with true percent identity. Exact alignment instead reproduces the percent-identity definition by construction: an edlib edit-distance prefilter [5] reranked by a Biopython Smith-Waterman pass [6] recovers the alignment ground truth almost exactly (Spearman 0.994, 91.1% top-1 agreement) at ~279× lower per-pair cost. Critically, alignment maps cleanly onto the literature identity bands — 98.7% species, 94.5% genus, 86.5% family [1,2] — so every selection is auditable against a published threshold rather than a hand-calibrated, space-specific cosine cutoff. A classical k-mer-8 TF-IDF retriever even beat the embedding on every closest-hit metric (75.5% vs 50.5% exact top-1), confirming that dense pooling averages away the substitution-level signal that fine 16S mapping requires. The concordance analysis of §3.4 then showed that, with this aligner, switching the entire upstream scorer from CPU to GPU changes the final representative for roughly 1 ASV in 1,000, and only among interchangeable equal-score tie members.

Second, **probability-weighted synthetic genomes preserve core/accessory uncertainty** for downstream flux-balance analysis [8]. Rather than collapsing each amplicon to a single best reference, the pipeline retains a weighted set of representatives whose PGFam content reflects genuine ambiguity, so core (high-confidence) and accessory (uncertain) gene content propagate into prFBA as per-function probabilities rather than spurious certainties; accessory genes enter the reconstruction as gap-fill candidates rather than hard reactions.

Third, the full 10,000-organism validation exposes a clear **region-length versus generalization tradeoff**, and resolves it at two taxonomic ranks. Include-self genus accuracy is near-ceiling (0.94–0.99) for all eight regions, but under honest self-exclusion (τ = 0.987 [1,2]) it collapses to 0.60–0.73, and gene capture falls from 95–99% to 68–72% — a near-uniform ~27-point leakage gap whose genus-accuracy counterpart, by contrast, widens for shorter reads. Short V4 (515F/806R [3,4]) generalizes worst (0.60 genus, 67.9% capture); full-length 16S best (0.73, 72.2%); V4–V5 and V7–V9 sit between — quantifying the cost of amplicon brevity for novel-organism inference. Crucially, the penalty is rank-dependent: family-level assignment stays robust under self-exclusion (0.79–0.86 across all regions), 12–19 points above genus, so for the short markers roughly one novel organism in five lands in the right family but a wrong-genus relative. The practical reading is that for a community dominated by organisms with close database relatives, even a short marker recovers nearly all of the metabolism; for genuinely novel organisms a short amplicon still places them confidently at the family level but leaves roughly a third of the gene complement and a third of genus calls unresolved — a gap that longer reads (full-length 16S, V1–V3) measurably narrow.

## Limitations and Future Work

These conclusions rest on the full 10,000-genome, 8-region, 67,991-amplicon run (the 500-genome pilot of §3.5 is retained only as the developmental checkpoint); the reported accuracies carry Wilson 95% intervals and the per-region comparisons are resolved at full statistical power. Six of eight regions are Bacteria-scoped, so only V4 and V4–V5 measure archaeal performance despite the archaeal sources in the draw. In-silico PCR (emax = 3) cannot capture wet-lab amplification bias, chimeras, or sequencing error, so reported accuracies are optimistic upper bounds. On the retrieval side, the GPU Smith-Waterman backend is currently CUDA-only; a detected Apple Metal device falls back to CPU scoring until a Metal kernel is implemented. Future work should propagate probability-weighted genome assignments through end-to-end community prFBA, validate against paired metagenomes, extend the panel with archaea-competent primers and long-read full-length 16S to test whether the generalization gap narrows with read length, and complete the Metal scoring backend to broaden hardware portability.

## References

1. Yarza P, et al. (2014) Uniting the classification of cultured and uncultured bacteria and archaea using 16S rRNA gene sequences. *Nature Reviews Microbiology*.
2. Stackebrandt E, Goebel BM (1994) Taxonomic note: a place for DNA-DNA reassociation and 16S rRNA sequence analysis in the present species definition in bacteriology. *International Journal of Systematic Bacteriology*.
3. Caporaso JG, et al. (2011) Global patterns of 16S rRNA diversity at a depth of millions of sequences per sample (Earth Microbiome 515F/806R primers). *PNAS*.
4. Parada AE, Needham DM, Fuhrman JA (2016) Every base matters: assessing small subunit rRNA primers for marine microbiomes; Apprill A, et al. (2015) Minor revision to V4 region 806R primer. *Environmental Microbiology / Aquatic Microbial Ecology*.
5. Šošić M, Šikić M (2017) edlib: a C/C++ library for fast, exact sequence alignment using edit distance. *Bioinformatics*.
6. Cock PJA, et al. (2009) Biopython: freely available Python tools for computational molecular biology and bioinformatics (PairwiseAligner). *Bioinformatics*.
7. Dalla-Torre H, et al. (2024) Nucleotide Transformer: building and evaluating robust foundation models for human genomics. *Nature Methods*.
8. Arkin AP, et al. (2018) KBase: The United States Department of Energy Systems Biology Knowledgebase; Henry CS, et al. (2010) ModelSEED high-throughput generation of genome-scale metabolic models. *Nature Biotechnology*.
9. Olson RD, et al. (2023) Introducing the Bacterial and Viral Bioinformatics Resource Center (BV-BRC/PATRIC). *Nucleic Acids Research*.
10. Okuta R, et al. (2017) CuPy: a NumPy-compatible library for NVIDIA GPU calculations. *Proceedings of Workshop on ML Systems (NIPS)*.
