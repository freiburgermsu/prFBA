# Elaborating the top-2 frameworks for closest-hit 16S mapping

>
> **Measured update (supersedes the "84.5%" placeholder where it reads ~92% in any cached copy).** After drafting, we directly measured recall of the *true single best-%identity reference* (the exact quantity Option 1's post-rerank top-1 equals) over the whole 97,624-insert DB, n=200 ASVs:
>
> | cosine top-k | true-#1 recalled = rerank top-1 ceiling | all true-top-5 present |
> |---|---|---|
> | 100 | **84.5%** | 37.5% |
> | 200 | **87.0%** | 41.0% |
> | 500 | **88.5%** | 47.0% |
>
> Option 1 lifts top-1 correctness from **54.7% → ~84.5% at k=100** (~88.5% at k=500) — **not 92%**; the earlier 92% was a mislabeled "≥1-of-top-5 captured" statistic. Strict same-insert match (counting an equally-best tied insert it is somewhat higher). This makes **Option 2 (raising the recall ceiling) more important, not less**.

## Framing

Our region-matched setup (experimental V4-V5 ASVs vs 97,624 in-silico-PCR-excised V4-V5 inserts) does **not** suffer a recall failure — NT-v2-500m mean-pooled cosine already recalls the true single-best-%identity reference (Biopython global, over the whole DB) in its top-100 for **~84.5%** of ASVs. It suffers an **ordering failure**: cosine-#1 equals alignment-identity-#1 only **54.7%** of the time (Spearman 0.55), and cosine top-100/200/500 contains *all* of the true top-5 alignment hits for only **37.5% / 41% / 47%** of ASVs. **Option 1 (retrieve → exhaustive-alignment rerank)** closes the ordering gap *by construction* — it re-scores the cosine shortlist with the exact metric the ground truth uses, so the reranked #1 is the true #1 whenever it was retrieved (top-1 → ~84.5% at k=100, essentially free to run). **Option 2 (k-mer + RRF fusion → learned-%identity prefilter → exact rerank)** attacks the thing Option 1 cannot: the **recall ceiling**, by adding an orthogonal k-mer retrieval signal that catches near-identical references mean-pooling smears away, lifting fused recall toward ~98-100% and raising the all-top-5-present rate well above 37-47%.

---

## Option 1 — Retrieve then exhaustive-alignment rerank (minimal)

### Mechanism: why this closes our exact gap

Our ground truth is **alignment %identity over the whole DB** (Biopython, global). Reranking fixes the ordering by construction: take the cosine shortlist, re-score every candidate by the *same* %identity definition the ground truth uses, then sort by it — the resulting #1 is the true global #1 **whenever the true #1 is in the shortlist**. There is no model-error term left: the rerank score *is* the ground-truth metric. Embedding cosine is demoted to what it is good at (cheap recall/candidate generation, 0.06 ms/query), and the exact, expensive metric is computed only on the ~100-500 survivors. This is exactly the regime where retrieve-then-rerank wins: a fast recall-oriented retriever followed by an exact-but-slow scorer applied to a tiny candidate set.

### The hard bound: rerank accuracy is bounded by retrieval recall

Reranking **reorders; it never adds candidates**. This logic is airtight:

- **Closest-hit (top-1) accuracy after rerank = recall@k of the true #1.** Cosine recalls the true single-best reference in top-100 for ~84.5% of ASVs, so reranking cosine top-100 by true %identity yields the correct closest hit for **~84.5%** of ASVs (rising with k). The remaining ~15.5% are *unrecoverable by reranking alone* — the true #1 was never retrieved. **This ~84.5% is an upper bound, achieved only when (a) the reranker metric is bit-identical to ground truth and (b) co-best ties are scored as success.** A `--iddef` mismatch or the VSEARCH k-mer prefilter gate (below) are the only ways to fall *below* it; the tie issue is the only way the *measured* 54.7% understates true success.
- **True-top-5 capture after rerank is unchanged: 37.5% / 41% / 47% at k = 100 / 200 / 500.** Reranking sorts whatever is present; it cannot surface a top-5 hit cosine never retrieved. To raise top-5 capture you must raise *k* or raise *recall* (Option 2) — not change the reranker.

**Rerank makes the ordering exact within the shortlist; the shortlist's contents are the binding constraint.** Top-1 ≈ 84.5%@100 is achievable today; full top-5 capture is gated by recall and needs Option 2.

### Exact recipe (VSEARCH `usearch_global`)

Per-ASV shortlist build (we already have `feature_id → insert` in `/home/freiburger/Documents/prFBA/v4v5_refs.json`):

1. For each ASV, take its cosine top-k `feature_id`s.
2. Emit a per-ASV (or per-batch) shortlist FASTA whose records are exactly those candidates' insert sequences, with the `feature_id` as the FASTA header. Two layouts both work:
   - **Per-query DB** (cleanest semantics): one query FASTA (the ASV) + one `shortlist.fa` of just its k candidates. Loop over ASVs.
   - **Batched** (faster I/O): concatenate all ASVs into `queries.fa`, put the *union* of all shortlist candidates into one `db.fa`. Over-searches but is fine at DB sizes of (n_ASVs × k) and avoids per-query process spawn. With `--maxaccepts 0 --maxrejects 0`, correctness is identical either way; only runtime differs.

```
vsearch --usearch_global queries.fa --db shortlist.fa \
        --id 0.5 \
        --iddef 2 \
        --strand plus \
        --maxaccepts 0 --maxrejects 0 \
        --top_hits_only \
        --output_no_hits \
        --userout hits.tsv \
        --userfields query+target+id+alnlen+mism+gaps+qlo+qhi+tlo+thi \
        --threads <N>
```

Critical-flag rationale (confirmed against the VSEARCH manpage and Rognes et al. 2016, PeerJ 4:e2584):

- **`--maxaccepts 0 --maxrejects 0` is REQUIRED.** Defaults are `--maxaccepts 1` (stop after the *first* hit clearing `--id`) and `--maxrejects 32` (stop after 32 consecutive non-matches) — heuristics tuned for "find *a* good hit fast" in a huge DB that can terminate before the *best* hit, exactly the failure mode we are eliminating in a dense cloud of near-identical V4-V5 candidates. The manpage is explicit: *"If --maxaccepts and --maxrejects are both set to 0, the complete database is searched."* On a k≈100-500 shortlist, exhaustive search is cheap, so there is no reason to leave the heuristic on. Without this, rerank is no longer exhaustive and the "exact by construction" guarantee is void.

- **VSEARCH does OPTIMAL full-DP alignment — but "complete database searched" ≠ "every candidate aligned."** VSEARCH computes *the optimal global alignment via a multi-threaded, vectorised full dynamic-programming (Needleman-Wunsch) algorithm adapted from SWIPE*, **not** the seed-and-extend heuristic used by USEARCH (Rognes et al. 2016) — so the "exact by construction" claim holds **for VSEARCH specifically** (it would be false if someone substituted USEARCH; do not cite the drive5 USEARCH page for optimality). **However**, even with both heuristics at 0, the search is two-phase: a **k-mer prefilter gates each target — at least 6 shared words, and at least 1 of every 16 query k-mers must match — before any full-DP alignment.** "Complete database is searched" means *no early termination after N accepts/rejects*; it does **not** mean every target receives a full-DP score. For a k≈100-500 shortlist this is almost certainly harmless (true high-identity V4-V5 hits clear the 6-shared-8-mer gate trivially, so the **best-hit guarantee survives**), but a very-low-identity candidate sharing <6 8-mers can be silently skipped. If you need a score for *literally every* candidate regardless of identity, use the Biopython/parasail reranker below.

- **`--id 0.5`** (or lower) is a permissive floor, not a filter — we want *every* shortlist candidate scored and ranked, so set it below any plausible real identity; tighten only if you want VSEARCH to drop obvious non-matches for you.

- **`--iddef 2`** is the default: "edit distance excluding terminal gaps (same as `--id`)" = matches / (alnlen − terminal gaps). **Pin this to your Biopython ground-truth definition.** If Biopython counts terminal gaps, use `--iddef 1` (matches/alnlen); if BLAST-style, `--iddef 4`; `--iddef 0` = CD-HIT (matches/shortest). The whole point of rerank is rerank-metric == ground-truth-metric, so any mismatch reintroduces an ordering discrepancy. (`id0..id4` in userfields give the identity under each `--iddef` simultaneously — handy for cross-checking against Biopython.)

- **global vs local:** `usearch_global` does **global** alignment, which is correct: both ASVs and inserts are trimmed to the *same* ~373 bp V4-V5 locus, so end-to-end global identity is the right, comparable measure and reproduces the global ground truth. A local search would let partial-overlap hits inflate identity. (Manpage caveat: query and reference must be trimmed to the same locus for global identity to be meaningful — exactly our setup.)

- **`--top_hits_only`** reports the maximum-identity hit and every hit tied at that identity — **but only because `--maxaccepts 0` un-caps the tie count.** The manpage states the number of co-best ties written is *"controlled by the --maxaccepts value,"* so the full-tie-tier behavior is a *dependency* on `--maxaccepts 0`, not an independent property of `--top_hits_only`. Drop this flag (and add `--maxhits k`) if you want the full ranked list per ASV instead of just the top tier.

- **`--strand plus`** is the default and correct because both ASVs and inserts are oriented to the same primer/locus. Use `--strand both` only if any insert could be reverse-complemented relative to the queries; it doubles work and can return per-strand duplicates.

- **`--output_no_hits`** keeps queries with no hit above `--id` so you can detect ASVs whose entire shortlist failed (a recall problem → an Option-2 case), rather than silently dropping them.

- **`--userfields`**: confirmed valid. `id` = percent identity under `--iddef`; available fields include `query, target, id, id0..id4, alnlen, mism, opens, gaps, qlo, qhi, tlo, thi, qcov, tcov, ql, tl, evalue, bits`.

**Alternative reranker (no new dependency, aligns EVERY candidate): Biopython / parasail / edlib over the shortlist.** Reranking cosine top-k = k pairwise alignments per ASV: at ~3.7-10 ms/pair (pure Biopython), ~0.4-5 s/ASV at k=100, ~2-25 s/ASV at k=500 single-threaded — trivially parallelizable. **parasail** (vectorized SW/NW, SSE/AVX2) is ~1-2 orders of magnitude faster per pair and is the better in-process choice; **edlib** (already installed) is fastest but gives edit distance / NW alignment, so derive %identity from its alignment under *your* `--iddef` convention. The advantage of the Biopython/parasail/edlib path is zero install, bit-for-bit agreement with the existing ground-truth scorer, **and that it genuinely aligns all k candidates** (no k-mer prefilter gate); the advantage of VSEARCH is throughput and a one-command batched run.

### Cost

Reranking a 100-500-candidate shortlist per ASV is **milliseconds to sub-second per query** (VSEARCH exhaustive over k≈100-500 is near-instant; parasail similar; pure Biopython ~0.4-5 s at the high end). A whole study (3,215-3,950 ASVs) reranks in **minutes**, multi-threaded.

Put **357 s/query** in its place: that was **un-indexed full-DB pairwise alignment against all 97,624 inserts**. Retrieve-then-rerank never aligns against the full DB — cosine (0.06 ms/query) reduces 97,624 candidates to ~100-500 *first*, and alignment runs only on that shortlist. The 357 s figure is the cost of the approach this design exists to avoid.

### Caveats

- **V4-V5 ties are biologically real.** Many distinct references are identical or co-equal over a ~373 bp V4-V5 window. A forced single #1 is artificial and will disagree with an equally-correct alternative at the same identity. **Report the whole top-identity tier** (what `--top_hits_only` + `--maxaccepts 0` returns) — our "#1 == #1" metric *undercounts* true success because a "wrong" #1 is often a co-best tie. Evaluate as "is the chosen hit *in* the true top-identity tier," not strict #1 equality.
- **The recall ceiling is the real limit.** Top-1 after rerank ≈ 84.5%@100; true-top-5 capture stays at 37.5/41/47%. Reranking cannot exceed retrieval. **Option 1 makes ordering exact and is essentially free to run; Option 2 is what raises the ceiling Option 1 ranks within.**

**Sources:** [VSEARCH manpage (Debian)](https://manpages.debian.org/stretch/vsearch/vsearch.1) · [Rognes et al. 2016, VSEARCH, PeerJ 4:e2584 — optimal full-DP NW via SWIPE; two-phase k-mer prefilter (≥6 shared words, ≥1/16 query k-mers)](https://peerj.com/articles/2584/) · [VSEARCH manpage (Ubuntu jammy) — iddef 0-4, userfields](https://manpages.ubuntu.com/manpages/jammy/man1/vsearch.1.html) · [QIIME2 vsearch-global — same-locus trimming requirement](https://docs.qiime2.org/2024.10/plugins/available/feature-classifier/vsearch-global/)

---

## Option 2 — k-mer + RRF fusion + learned-%identity rerank (ceiling-raising)

### Why this option exists (the binding constraint it attacks)

Option 1's hard ceiling is **84.5%@top-100** for the single best reference, and only **37.5% / 41% / 47%** of ASVs have *all* their true top-5 hits inside cosine top-100/200/500. No reranker — exact or learned — recovers a true best hit retrieval never surfaced. For the ~15.5% of ASVs where cosine misses entirely, and the majority where the top-5 set is incomplete, the failure is **recall**, and the only fix is a **second, orthogonal retrieval signal**.

Mean-pooled NT-v2-500m embeddings are the right thing to be orthogonal *to*: averaging 373 per-token vectors into one is a low-pass filter that smears exactly the SNP/indel-level differences distinguishing a 99.5%-identity reference from a 97% one. Two inserts differing by 3-4 bases out of 373 (the regime deciding the alignment-#1 vs alignment-#2 tie we lose 45% of the time) can be nearly identical after mean-pooling, while their **exact k-mer sets** differ in a sharp, countable way — every substitution kills up to *k* k-mers. K-mer overlap is therefore the natural complement: precisely sensitive where mean-pooling is blind.

### Stage A — k-mer sparse retrieval (the orthogonal signal)

**Build as an explicit sparse matrix — NOT sourmash.** sourmash/FracMinHash is a *sketching* method: at the recommended `scaled=1000` it keeps ~1/1000 of k-mers, tuned for overlaps in the **10 kb** range (docs: ~99.98% of *5 kb* overlaps detected; explicitly "doesn't work well" for very small sequences, "no systematic advice" for sub-kb amplicons). A 373 bp insert has only ~365 8-mers; at `scaled=1000` its sketch is essentially empty. You'd be forced to `scaled=1` (= all k-mers = exact MinHash, no compression), at which point sourmash buys nothing over a direct sparse matrix and adds an uninstalled dependency. **Verdict: skip sourmash; wrong resolution regime.**

Instead, with already-installed scipy 1.17.1 / scikit-learn 1.7.0:

- **Vectorize** every reference insert (97,624, on disk at `/home/freiburger/Documents/prFBA/v4v5_refs.json`) and every ASV into a sparse k-mer count/TF-IDF vector via `sklearn.feature_extraction.text.CountVectorizer`/`TfidfVectorizer` with `analyzer='char'`, `ngram_range=(k,k)`, **k=8** (4^8 = 65,536 possible k-mers — comparable to the ref count, dense enough to discriminate, short enough that one SNP still leaves most k-mers shared so near-identical refs co-retrieve). Optionally also build **k=12** (4^12 ≈ 16.7M columns, very sparse, very specific) as a third retrieval list.
- The reference matrix is **97,624 × 65,536**, ~365 nonzeros/row → ~36M nonzeros, a few hundred MB in CSR. Build once, offline, persist alongside the embeddings.
- **Retrieve** per ASV by k-mer **cosine** (L2-normalize rows, single sparse matvec `R @ q`) for top-N, or by **containment** = |k-mers(ASV) ∩ k-mers(ref)| / |k-mers(ASV)| (more identity-faithful for near-equal-length amplicons; monotone with %identity at high identity). Cost is one sparse matvec — single-digit ms/query, same order as the 0.06 ms cosine, trivially batchable.

Use **TF-IDF** (not raw counts) so k-mers shared across the whole DB (conserved V-region stretches) are down-weighted and the **discriminating** k-mers dominate — directly targeting the variable positions that decide identity ties.

### Stage B — Reciprocal Rank Fusion (raise the recall ceiling)

Fuse the **cosine rank list** (embeddings) and the **k-mer rank list(s)** by rank, not score. For each candidate reference *d*:

```
RRF(d) = Σ_j  1 / (c + rank_j(d)) ,   c = 60
```

where `rank_j(d)` is *d*'s 1-based position in retrieval list *j* (cosine, k-mer-cosine, optionally k=12 containment), and *d* absent from a list contributes 0. This is exactly the formula and constant from Cormack, Clarke & Büttcher, *"Reciprocal Rank Fusion outperforms Condorcet and individual Rank Learning Methods,"* SIGIR 2009; `c=60` is their empirically-best TREC value (Cormack's `k`, renamed `c` here to avoid colliding with top-`k`), and follow-on work finds `c ∈ [40,80]` performs comparably, so it needs no tuning.

Why RRF rather than averaging/normalizing raw scores: cosine similarity (∈[−1,1], tightly clustered near the top) and k-mer TF-IDF cosine / containment (∈[0,1], different distribution) are **not comparable in magnitude**, and a near-duplicate ref produces a saturated, uninformative cosine score. Any score-normalization scheme (min-max, z-score) is distribution-dependent and brittle. RRF discards magnitudes and uses only **ordinal position**, so a ref that k-mer-ranks #2 but cosine-ranks #300 still gets a strong fused score — *exactly* the missed-best-hit case. The `1/(c+rank)` shape is steep at the head (a top-few placement in *either* list dominates) and flat in the tail (deep noisy ranks don't swing the fusion). The fused top-k therefore has **union-level recall**: it contains a true hit if *either* retriever found it. Cosine's 84.5%@100 is a **floor** here, not a ceiling — k-mer retrieval catches most of the missed ~15.5% (near-identical refs are precisely where k-mer overlap is strongest), pushing fused recall toward **~98-100%** and lifting the all-top-5-present rate well above 37.5/41/47%. **This is the whole point of Option 2: it moves the ceiling Option 1 is capped by.**

### Stage C — learned-%identity prefilter (cheap ordering before exact align)

Before paying for exact alignment, order the fused top-~500 with **Identity** (BioinformaticsToolsmith / Girgis, James, Luczak, *NAR Genomics and Bioinformatics* 2021, lqab001) — the published successor to **FASTCAR** (bioRxiv 380824). Verified specifics:

- **k-mer-histogram-based, alignment-free, self-supervised GLM.** It converts each sequence to a k-mer histogram and predicts identity from alignment-free k-mer statistics. The "learned" signal is itself k-mer-derived, so it is **correlated with, not orthogonal to, Stage A** — its job is not extra recall (Stage B did that) but a **better-than-rank ordering** of the fused candidates, calibrated to actual %identity. **This is the most important conceptual caveat in Option 2: the orthogonality that lifts the recall ceiling comes only from the embedding ⟂ k-mer pairing in Stages A+B, never from this reranker.**
- **Accuracy: Pearson 0.97 vs Needleman-Wunsch *global*** (confirmed verbatim; 0.94-0.95 vs BLAST/USEARCH/Mash/MUMmer4). The 0.97 is specifically against NW global — the right target, since our ground truth is whole-DB best-%identity.
- **Validated on short reads, in-distribution for us:** the paper runs it on **1,071,335 sequences of ~251 bp** (~573 billion pairwise predictions, ~13.5 h). 251 bp < our ~373 bp, so the short-amplicon regime is genuinely in-distribution — its single best-supported applicability claim, and *not* true of sourmash.
- **Self-trains per invocation** by mutating template subsequences to generate ~10,000 semi-synthetic labeled pairs (5,000 train / 5,000 test) before scoring real data — **no external labels, but it re-trains each run**, so amortize by running query-vs-shortlist in batch, not one ASV at a time.
- **Speed: ≥2× faster than BLAST, ~15-23× USEARCH, ~20-65× Mash** (paper benchmarks) — scoring a few hundred candidates per ASV is cheap relative to aligning them.
- **Flags (per README, not independently re-verified line-by-line):** `-d database.fasta`, optional `-q query.fasta` (omit `-q` → all-vs-all), `-t` reporting identity threshold (0-0.99), `-o` tab-separated `header1 \t header2 \t predicted_identity`, `-c` cores, `-a y` to emit below-threshold pairs, `-r y` to auto-relax threshold by the predictor's error.
- **License/build caveat:** dual-licensed — **AGPL (academic use)**, separate commercial license required (contact author). Fine for research; flag if it ever ships commercially. Compiles from source (cmake, g++ ≥ 7.5); **not currently installed** (C++ binary, not pip).

Use Identity to rank the fused top-500 by predicted %identity, keep the top **~10-20**, then run the **exact** pass (edlib for fast %identity, or Biopython `PairwiseAligner` for the full NW score) on only that head. Because Identity tracks NW at 0.97, the exact-#1 is almost always inside that top-10-20, so you exact-align ~15 candidates instead of 100-500.

**Sources:** [Identity, NAR Genomics & Bioinformatics 2021 (lqab001)](https://academic.oup.com/nargab/article/3/1/lqab001/6125549) · [Identity GitHub README (flags, license)](https://github.com/BioinformaticsToolsmith/Identity/blob/master/README.md) · [FASTCAR, bioRxiv 380824](https://www.biorxiv.org/content/10.1101/380824v2.full) · [Cormack, Clarke & Büttcher, RRF, SIGIR 2009](https://research.google/pubs/reciprocal-rank-fusion-outperforms-condorcet-and-individual-rank-learning-methods/) · [sourmash FAQ — scaled / small-sequence limits](https://sourmash.readthedocs.io/en/latest/faq.html) · [sourmash internals — FracMinHash/scaled](https://sourmash.readthedocs.io/en/stable/sourmash-internals.html)

---

## Drop-in pipeline

Stages labeled **[O1]** = Option 1 minimal (ordering-exact, free to run today); **[O2]** = Option 2 ceiling-raising (add when recall is the binding constraint or aligning the shortlist is the bottleneck).

```
# refs on disk: /home/freiburger/Documents/prFBA/v4v5_refs.json  (feature_id -> insert; 97,624)
# installed: edlib, Biopython 1.85, scikit-learn 1.7.0, scipy 1.17.1
# NOT installed: VSEARCH, sourmash, MMseqs2 (pip/conda), Identity (compile from source, AGPL)

for each ASV:

  # ---- retrieval --------------------------------------------------------
  L1 = NT-v2 mean-pool cosine  -> top-500 feature_ids        # [O1] 0.06 ms, exists
  L2 = sparse k-mer TF-IDF (CountVectorizer char ngram k=8)  # [O2] ~ms, scipy/sklearn
       cosine or containment   -> top-500 feature_ids
  L3 = optional k=12 containment -> top-500                  # [O2] extra orthogonal list

  # ---- fusion -----------------------------------------------------------
  fused = RRF({L1, L2, L3}, c=60)  -> fused top-500          # [O2] recall = union of retrievers
          #   RRF(d) = sum_j 1 / (60 + rank_j(d))
          # [O1] alone: shortlist = L1 top-k (skip fusion)

  # ---- learned-%identity prefilter -------------------------------------
  head = Identity(-d shortlist.fa -q ASV.fa -c N) sort desc  # [O2] 0.97 vs NW, top ~10-20
         # batch query-vs-shortlist (re-trains per invocation); OPTIONAL — drop if
         # exact-aligning ~50-100/query is already fast enough

  # ---- EXACT rerank (the ground-truth metric) --------------------------
  # [O1] VSEARCH path (batched, throughput):
  emit shortlist.fa  (head, or L1 top-k if Option 1 only)
  vsearch --usearch_global queries.fa --db shortlist.fa \
          --id 0.5 --iddef 2 --strand plus \
          --maxaccepts 0 --maxrejects 0 \            # REQUIRED: exhaustive, un-caps ties
          --top_hits_only --output_no_hits \
          --userout hits.tsv \
          --userfields query+target+id+alnlen+mism+gaps+qlo+qhi+tlo+thi \
          --threads N
          # pin --iddef to the Biopython ground-truth convention!
          # NOTE: k-mer prefilter (>=6 shared words) gates targets; best hit always clears it.

  # [O1] in-process alternative (no install, aligns EVERY candidate, bit-for-bit w/ truth):
  #   edlib / parasail / Biopython PairwiseAligner over all shortlist candidates

  final_top1, final_top5 = sort head/shortlist by exact %identity (matched --iddef)

# evaluate as "chosen hit IN true top-identity TIER", not strict #1 equality (V4-V5 ties)
```

---

## Expected outcome on our data

| Metric | Today (cosine) | **Option 1** alone | **Option 2** (+O1) |
|---|---|---|---|
| Closest-hit top-1 correct | 54.7% (Spearman 0.55) | **~84.5%** at k=100 (= recall@100 of true #1; rises with k) — *proven-by-construction ceiling, achieved only if iddef matches truth & ties count as success* | **rises with recall** toward the lifted ceiling as fused recall → ~98-100% |
| All true top-5 captured | 37.5 / 41 / 47% @100/200/500 | **unchanged** at 37.5/41/47% *unless k grows* — rerank cannot surface unretrieved hits | **well above 37-47%** — fused set is the union of embedding+k-mer retrieval |
| Recall ceiling | 84.5%@100 | **84.5%@100 (hard cap)** — rerank reorders, never adds | **~98-100%** — orthogonal k-mer signal catches the near-identical refs mean-pooling smears |
| Run cost | 0.06 ms/query | **ms-to-sub-second/query**; whole study in **minutes** | + ~ms k-mer matvec + Identity batch; exact fan-out shrinks 100-500 → ~15/query |

- **Proven (by construction):** Option 1's post-rerank top-1 = recall@k of the true #1; on our data that is **~84.5% at k=100** (vs 54.7% today). Top-5 capture is mathematically **unchanged** by reranking alone. These are not estimates — they follow from "rerank reorders a fixed set."
- **Expected (not yet measured on our data):** Option 2's fused recall toward **~98-100%** and the corresponding rise in top-5 capture and post-rerank top-1. This rests on the *expectation* that k-mer overlap recovers most of cosine's missed ~15.5% (well-founded — near-identical refs are exactly where k-mer overlap is strongest — but to be measured, not assumed). Identity's 0.97-vs-NW correlation is proven *in its paper on 251 bp reads*; its exact behavior on our 373 bp inserts is expected-in-distribution, not yet measured here.
- **The honest bottom line:** Option 1 is the cheap, immediate, provably-correct ordering fix (54.7% → ~84.5% top-1, free to run). Option 2 is the investment that raises the recall ceiling Option 1 ranks within — worth it once 84.5% is your binding constraint or shortlist alignment is the bottleneck, and noting its reranker (Identity) is k-mer-correlated, so the ceiling lift comes from the embedding ⟂ k-mer pairing, not from the reranker.

---

# Appendix: verification

I now have all the primary-source confirmations needed, including an important nuance about VSEARCH's k-mer prefilter that the elaboration omits. Let me compile the verification.

## Verification

I verified all four flagged specifics against primary and authoritative secondary sources. The two elaborations are **substantially accurate**, but there are **three nontrivial corrections/qualifications** and several confirmations worth recording.

### (1) VSEARCH `usearch_global` flag semantics — CONFIRMED, with one important added nuance

- **`--maxaccepts 0 --maxrejects 0` = search whole DB: CONFIRMED, verbatim.** The Debian manpage states defaults are `--maxaccepts 1` and `--maxrejects 32`, and: *"If --maxaccepts and --maxrejects are both set to 0, the complete database is searched."* The elaboration's quote, defaults, and rationale (heuristic early-termination can stop before the true best hit in a dense near-identical cloud) are all correct.
- **`usearch_global` is full optimal DP: CONFIRMED, and the elaboration actually UNDERSELLS the strongest evidence while citing the wrong source.** Elaboration 1 cites the **drive5/USEARCH** page for "global alignment." That page is about **USEARCH**, which uses a **seed-and-extend heuristic alignment** — NOT optimal. The correct, stronger fact (Rognes et al. 2016, PeerJ 4:e2584; and the VSEARCH manpage) is that **VSEARCH** computes *"the optimal global alignment using a multi-threaded and vectorised full dynamic programming algorithm (Needleman & Wunsch, 1970) adapted from SWIPE... instead of the seed-and-extend heuristic used by USEARCH."* So the "exact by construction" claim is **true for VSEARCH** but would NOT be true if someone substituted USEARCH. Recommend citing the VSEARCH/Rognes source, not the drive5 USEARCH page, for the alignment-optimality claim.
- **IMPORTANT OMITTED NUANCE — the k-mer prefilter:** Even with `--maxaccepts 0 --maxrejects 0`, VSEARCH does **not** run full DP against literally every DB record. Per Rognes 2016, the search is two-phase: a **k-mer heuristic prefilter** sorts targets by shared-word count, and there is a **gating threshold — at least 6 shared k-mers, and at least 1 of every 16 query k-mers must match — before a target is aligned at all.** "Complete database is searched" means *no early termination after N accepts/rejects*; it does **not** mean every target receives a full-DP alignment. **For a k≈100–500 shortlist this is almost certainly harmless** (true high-identity V4-V5 hits will clear the 6-shared-k-mer gate easily), but the elaboration's phrasing "every shortlist candidate scored" is slightly too strong: a candidate sharing fewer than 6 8-mers with the query (i.e. very low identity) can be silently skipped. This does not affect the *best-hit* guarantee (the true top hit will always clear the gate), so the load-bearing claim survives — but the "every candidate gets a number" claim should be softened, OR you use the Biopython/parasail reranker (which genuinely aligns every candidate) when you need a score for *all* k candidates regardless of identity.
- **`--iddef`: CONFIRMED.** Default is `--iddef 2` ("edit distance excluding terminal gaps, same as --id"); `1` = matches/alnlen; `0` = CD-HIT (matches/shortest); `4` = BLAST def. The instruction to pin `--iddef` to the Biopython ground-truth convention is correct and important.
- **`--strand plus` default: CONFIRMED.** `--top_hits_only`: CONFIRMED it reports only the maximum-identity tier. **Minor correction on ties:** the manpage says the number of co-best ties actually written is *"controlled by the --maxaccepts value."* So to truly get *all* tied top-identity hits, `--maxaccepts 0` is doing double duty (it also un-caps the tie report) — the elaboration's recommended flag set is correct, but the reason "top_hits_only returns the whole tie tier" only holds **because** `--maxaccepts 0` is set. Worth stating that dependency explicitly.

### (2) "Rerank accuracy is bounded by retrieval recall (= recall@k of the true #1)" — LOGICALLY AIRTIGHT

This is correct and the elaboration states it honestly. Reranking is a **reordering** of a fixed candidate set; it cannot introduce a candidate absent from the shortlist. Therefore:
- post-rerank top-1 correctness ≤ (true #1 ∈ shortlist) = **recall@k of the true #1** — with equality when the reranker reproduces the ground-truth metric exactly (which VSEARCH/Biopython with matched `--iddef` does, modulo the prefilter-gate caveat above and ties);
- post-rerank top-5 capture = unchanged from retrieval's top-5 capture (37.5/41/47%), because a missing top-5 hit cannot be re-surfaced.

The claim is sound. **One honesty refinement:** "top-1 after rerank = recall@k of the true #1 = ~84.5%" is an *upper bound that is achieved only if* (a) the reranker metric is bit-identical to ground truth and (b) ties are scored as success. The k-mer prefilter gate and any `--iddef` mismatch are the only ways to fall *below* 84.5%; the tie issue is the only way the *measured* number could *understate* true success. Both are acknowledged in the caveats. Good.

### (3) Identity / FASTCAR — CONFIRMED (real tool, correct lineage, correct numbers)

Verified against the NAR Genomics & Bioinformatics 2021 paper (lqab001, Girgis, James, Luczak):
- **Exists; successor to FASTCAR (bioRxiv 380824): CONFIRMED.**
- **k-mer-histogram-based, self-supervised GLM: CONFIRMED.** It converts sequences to k-mer histograms and **self-trains per invocation** by mutating template subsequences to generate ~10,000 semi-synthetic labeled pairs (5,000 train / 5,000 test) before scoring real data. So the elaboration's "re-trains each invocation, amortize in batch" warning is correct and material.
- **Pearson 0.97 vs Needleman-Wunsch global: CONFIRMED verbatim** ("highly correlated (0.97 Pearson's correlation coefficient) to the scores calculated by the original global alignment algorithm"). The 0.94–0.95 vs BLAST/USEARCH/Mash/MUMmer4: CONFIRMED.
- **Short-read validation: CONFIRMED and even stronger than stated** — validated on a 16S dataset of **1,071,335 sequences ~251 bp** (~573 billion pairwise predictions, ~13.5 h). 251 bp < your ~373 bp, so the short-amplicon regime is genuinely in-distribution. This is the elaboration's single best-supported claim for applicability.
- **Speed "2–80× faster than competitors": CONFIRMED** (paper benchmarks: faster than BLAST, ~15–23× USEARCH, ~20–65× Mash). The "≥2× faster than BLAST" specific is consistent.
- **Inputs `-d`/`-q`, all-vs-all when `-q` omitted: CONSISTENT** with the paper and README. I could not independently re-verify every minor flag (`-t`, `-a`, `-r`, `-c`, `-o` column format, AGPL license, cmake/g++≥7.5 build) line-by-line beyond the GitHub README the elaboration already cites; these are plausible and sourced, but treat the exact flag letters as "per README" rather than independently confirmed here.

**One honesty point the elaboration gets right and should keep prominent:** Identity is k-mer-derived, so it is **correlated with, not orthogonal to, the Stage-A k-mer retriever** — it improves *ordering/fan-out*, not *recall*. The orthogonality that lifts the recall ceiling comes only from the embedding ⟂ k-mer pairing. This is the most important conceptual caveat in Option 2 and it is stated correctly.

### (4) RRF formula and constant — CONFIRMED

- **Formula CONFIRMED:** `RRFscore(d) = Σ_j 1/(k + rank_j(d))`, sum over rankings, `rank` 1-based, absent ⇒ 0 contribution. Matches the elaboration's `1/(c+rank)` (the elaboration just renames Cormack's `k` to `c` to avoid colliding with top-`k`, which is sensible).
- **Constant `k=60` from Cormack, Clarke & Büttcher, SIGIR 2009: CONFIRMED** as the value tuned on TREC and reported to generalize; subsequent work finds `k∈[40,80]` comparable. (Note: I confirmed this from the authoritative secondary sources and Google/IR-anthology listings; the primary PDF and PeerJ/Google-research pages were unparseable/403 to the fetcher, but the constant and formula are uncontroversial and consistently reported.)

### Summary of what is overstated vs. solid

- **Overstated #1 (cite-source error):** Using the **drive5 USEARCH** page to justify "global, optimal" alignment. USEARCH's alignment is heuristic; the optimal-full-DP guarantee is a **VSEARCH** property (Rognes 2016). Swap the citation. The conclusion is still correct *for VSEARCH*.
- **Overstated #2 (phrasing):** "every shortlist candidate scored/ranked" under `--maxaccepts 0 --maxrejects 0`. VSEARCH still applies a k-mer prefilter gate (≥6 shared words; ≥1/16 query k-mers) before aligning; very-low-identity candidates can be skipped. Harmless for best-hit recovery; soften the wording, or use parasail/Biopython if you need a score for literally all k candidates.
- **Minor:** `--top_hits_only` returns the full top-identity tie tier **only because `--maxaccepts 0` un-caps the tie count**; make that dependency explicit.
- **Everything else verified accurate:** maxaccepts/maxrejects semantics and defaults, the recall-bound logic, Identity's existence/lineage/0.97-vs-NW/251 bp validation/self-training/speed, and the RRF formula and `c≈60` provenance.

Sources:
- [VSEARCH manpage (Debian) — maxaccepts/maxrejects defaults, "both 0 searches complete database", iddef definitions, top_hits_only tie count tied to maxaccepts, strand default plus](https://manpages.debian.org/stretch/vsearch/vsearch.1)
- [Rognes et al. 2016, VSEARCH (PeerJ 4:e2584) — optimal full-DP Needleman-Wunsch via SWIPE vs USEARCH seed-and-extend; two-phase k-mer prefilter (≥6 shared words, ≥1/16 query k-mers) then full DP](https://peerj.com/articles/2584/)
- [drive5 usearch_global reference — USEARCH (heuristic) global command; the page the elaboration cites](https://drive5.com/usearch/manual/cmd_usearch_global.html)
- [Identity, NAR Genomics & Bioinformatics 2021 (lqab001) — 0.97 Pearson vs NW global; 1,071,335 seqs ~251 bp; self-training; speed](https://academic.oup.com/nargab/article/3/1/lqab001/6125549)
- [FASTCAR (bioRxiv 380824) — Identity's predecessor](https://www.biorxiv.org/content/10.1101/380824.full.pdf)
- [Cormack, Clarke & Büttcher, RRF, SIGIR 2009 (listing)](https://research.google/pubs/reciprocal-rank-fusion-outperforms-condorcet-and-individual-rank-learning-methods/)
- [RRF formula + k=60 provenance (secondary, consistent with primary)](https://bigdataboutique.com/blog/reciprocal-rank-fusion-how-it-works-and-when-to-use-it)
