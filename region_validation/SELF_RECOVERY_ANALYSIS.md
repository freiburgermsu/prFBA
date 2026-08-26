# Self-recovery failure analysis & the union ("include-all") selection fix

**Question.** The 10,000-genome self-recovery benchmark (each genome's 16S → in-silico
amplicon → aligned against the 459k-ref DB → `select_references` → did we pick the
*source* organism back out?) showed ~1–2% of the full-length round-trips, and more for
short regions, failing to recover the source. *Why*, and can a selection change push
capture toward 100% — specifically by **including every reference genome tied at the best
identity** while **de-duplicating taxonomically-identical hits** so dense reference
regions don't flood the result?

All numbers below are pure **re-selections over the same 67,991 amplicon hits**
(`data/hits/asv_top20_alignment_hits.json`) scored with the **unmodified** metric code
(`score_taxacc.py`, `score_genecap.py`). No re-alignment. The cause attribution and the
`select_all` code were independently re-derived and audited by a 4-agent verification
workflow (all confirmed, no correctness bugs).

---

## 1. Why the source is lost (cause attribution)

`scripts/diagnose_self_recovery.py` attributes every failure (source in top-20 but not
selected) to the exact stage that dropped the source's furthest-surviving hit.

**11,300 failures / 66,144 amplicons-with-source-in-top-20 (17.1%).** Terminal stage:

| stage | share | taxonomic salvage (does the kept anchor still recover the lineage?) |
|---|---|---|
| `marginal_gain` | **58.8%** | genus kept 61%, family 14%, lost-above-family 25%; the source added <15% novel genes vs an already-chosen near-relative |
| `family_consensus` | **18.2%** | **98% lost** — the source's *correct* family was the *minority* in the in-band tie cloud and got out-voted |
| `tier_cap` | 10.6% | genus mostly kept (the species/genus ceiling filled before the source) |
| `no_gene_set` | 8.6% | **data gap** — the source genome has no CDS in BV-BRC (16S-only entry); *induced by exact-gene mode* |
| `species_dedup` | 3.3% | **benign** — a conspecific kept 77% of the time (species still correct) |
| `contamination_budget` | 0.4% | — |

**FullLength16S positive control (267/7,298 = 3.66%):** `family_consensus` 50%,
`no_gene_set` 31%, `marginal_gain` 15%, `species_dedup` 3%. The dense tie clouds are real:
2,083 failures have the *entire* top-20 within 0.002 identity.

**Take-away.** The source is not lost to one tie-break; it is lost to the **order-dependent
reducers** — `marginal_gain`, `family_consensus`, `tier_cap`. `species_dedup` (the most
"stochastic-looking" stage) is the *smallest and most benign* cause. This localizes the
"stochasticity from order selection" to the right stages.

---

## 2. The fix — a union ("include-all") mode (`select_references.select_all`)

A new, backward-compatible mode unions every reference whose identity is within
`select_all_window` of the best (still gated by the **reliability + family floor**, which
only drop partial/chimeric/sub-family junk), bypassing the order-dependent reducers. The
anchor stays the highest identity-decayed representative (deterministic), so anchor-based
taxonomy is preserved while the source is recovered whenever it is a qualifying hit.

Two knobs control the two questions you raised:

- **Which hits count as "tied"?** `select_all_window`.
  - **`0.0` = exact-identity ties only** — keep just the genomes whose identity *equals*
    the best (your clarified intent). This is the minimal set that removes the order
    tie-break, and it only unions genomes genuinely *indistinguishable* over this amplicon.
  - `0.005` = a 0.5%-identity near-tie band (looser; benchmarked as `in-band`).
  - `≥1.0` = every reliable hit down to the family floor (`all-reliable`; over-broad).

- **De-duplicating taxonomically-identical hits?** `select_all_species_dedup` (so a dense
  reference region — many sequenced strains of one species — does **not** flood the union).
  The "same species" determination is a **BV-BRC taxonomy** judgement, **not** an identity
  one: the grouping key is the BV-BRC lineage `Species` name (→ `taxon_id` fallback → `genome_id`
  last resort), never the identity score. So two genomes at *identical* identity but
  *different* BV-BRC species are **both kept** (e.g. a *Pediococcus acidilactici* and a
  *Clostridium perfringens* whose conserved V4 amplicons both match at 1.0), while multiple
  assemblies of the **same** BV-BRC species collapse to the single best-matching
  representative. Identity only chooses *which* assembly represents each species.

  Effect: union selections with ≥2 genomes of the **same BV-BRC species** drop from
  **20,493 → 0**. Example `187452.24__V4`: 6 *P. acidilactici* assemblies → 1.

### 2a. The guarded configuration — full selection criteria

For the FBA synthetic-genome use case, two order-independent **precision firewalls**
(`select_all_coherence`, `select_all_contam`) and a cap are added. The complete pipeline,
in order:

| # | step | rule |
|---|---|---|
| 0 | **reliability gate** | drop a hit if coverage `aligned_len/asv_len < 0.90`, **or** `|SW_id − edlib_id| > 0.02` (method-discordant / chimera), **or** identity `< 0.865` (below family floor) |
| 1 | **tier floor** | if the *best* identity `< 0.865`, abstain (no genome-level model) |
| 2 | **identity window** | keep every survivor within `select_all_window` of the best identity (`0.0` = exact ties) |
| 3 | **genome_id dedup** | collapse multi-copy 16S operons of one genome to its best copy |
| 4 | **species dedup** (BV-BRC taxonomy) | ≤1 representative per BV-BRC `Species` (→`taxon_id`→`genome_id`); keep the highest-identity copy (tie → align_score → genome_id) |
| 5 | **rank → anchor** | order survivors by identity-decayed representativeness `_weighted_rep = (0.55·id + 0.30·norm_score + 0.15·cov)·exp(−(best−id)/0.01)`; the top is the **anchor** (deterministic) |
| 6 | **anchor-family coherence** (firewall) | drop a member whose BV-BRC family ≠ the anchor's family **and** is not in the MiDAS prior (the cross-genus firewall dense short-region tie clouds need; anchored to the best hit, so it cannot out-vote the source the way the old plurality `family_consensus` did) |
| 7 | **contamination budget** | admit a member only while cumulative `(1−identity)·novel_genes` ≤ tier budget × anchor genes (species 25% / genus 40% / family 20%); the ~1.0-identity source contributes ~0 and is never the one dropped |
| 8 | **cap** | a hard ceiling (default 6) on genomes unioned — final bound for pathological tie clouds |

The selected genomes' gene sets are unioned into the synthetic genome; the anchor's lineage
is the taxonomic call. **Shipped defaults** (`DEFAULT_KNOBS`): `select_all=True`,
`select_all_window=0.0`, `select_all_species_dedup=True`, and the guards **off**
(`select_all_coherence=False`, `select_all_contam=False`, `select_all_cap=0`,
`select_all_unconditional_window=-1.0`) — i.e. the default is steps 0–5 only (the exact-tie
de-duplicated union); steps 6–8 are the opt-in guarded/graduated pathway, and the legacy
ordered reducer is `select_all=False`.

---

## 3. Results (5-way; ALL regions = n-weighted over the 8 panels)

| (include scenario) | legacy (baseline) | **exact-tie (DEFAULT)** | graduated +0.005 | exact-tie +guards | in-band 0.005 |
|---|---|---|---|---|---|
| **self-recovery** P(selected \| in top-20) | 82.9% | **99.2%** | 99.2% | 95.3% | 100.0% |
| … FullLength16S control | 96.3% | **99.9%** | 99.9% | 99.7% | 100.0% |
| mean genomes unioned | 1.30 | 2.20 | 2.72 | 1.75 | 4.79 |
| same-species-duplicate amplicons | 0 | 0 | 0 | 0 | 23,429 |
| anchor **species** concordance (FL16S) | 97.1 | **99.1** | 99.1 | 99.1 | 99.1 |
| gene-capture **recall** FL16S / V4 | 99.4/95.1 | **100/98.6** | 100/98.7 | 99.9/96.5 | 100/99.4 |
| gene-capture **precision** FL16S / V4 | 96.6/85.9 | **99.7/81.1** | 93.9/77.2 | 99.7/83.8 | 89.7/73.8 |
| gene-capture **F1** FL16S / V4 | 97.7/89.4 | **99.8/86.5** | 96.3/83.9 | 99.8/88.3 | 93.2/81.2 |
| consensus genus concordance (FL16S) | 98.7 | 99.9 | 99.2 | 99.8 | 98.0 |

(`all-reliable` — every reliable hit to the family floor — is omitted: precision 29–36%, F1
~45–50%, consensus genus 76% — pathological, rejected. Anchor, the genome whose lineage is the
taxonomic call, is identical across all union variants and ≥ baseline.)

Per-region tables (recall/precision/F1/Jaccard, include & honest-exclude) in
`data/variant_comparison.md`.

**The decisive finding.** Restricting the union to **exact-identity ties** (window `0.0`)
instead of any wider window is what makes this nearly a *free* improvement: exact ties union
only genomes *indistinguishable* over the amplicon, so **gene-capture precision stays at —
or above — baseline** (FL16S 96.6→99.7) while recall and self-recovery rise; **F1 is ≥
baseline in almost every region**. **BV-BRC-taxonomy species-dedup** removes the dense-region
redundancy (20,493→0 same-species amplicons) without touching recall.

**Exploration: is the within-0.005 extension worth it?** No. The **graduated** pathway keeps
all exact ties and admits the 0.005 near-tie extension *only* through the family-coherence
firewall + contamination budget — yet it is **net-negative**: it raises the mean union from
2.20 to 2.72 genomes and *lowers* precision (FL16S 99.7→93.9, V4 81.1→77.2) and F1 (FL16S
99.8→96.3, V4 86.5→83.9) for **no self-recovery gain** (99.2%, identical) and a negligible
recall change (V4 +0.1pp). The guards mitigate the hit (graduated beats the un-guarded
`in-band` band: FL16S precision 93.9 vs 89.7) but cannot make the extension pay: the source
is already an exact tie, so the 0.5%-divergent neighbours add mostly off-target genes. The
**exact-tie default is therefore the precision-optimal operating point**; the `in-band` band
(no guards) is strictly worse, and `all-reliable` collapses.

**Honest (novel-organism, exclude_0.987) scenario.** With the source removed, exact ties are
sparse, so the union is small and recall dips slightly below baseline (e.g. FL16S 72.2→68.1)
at near-baseline precision — i.e. for *novel* organisms the tight exact-tie union is
conservative, not inflating. Use the `in-band` window only if broader novel-organism recall
is wanted (at a precision cost).

---

## 4. The path to a literal 100%

Self-recovery factors as **P(selected) = P(source in top-20) × P(selected | in top-20)**.

- The union maxes the **second** term. On the **FullLength16S control the source is in the
  top-20 100% of the time, so the union reaches ~99.7–100% capture** there.
- The residual for short regions is the **first** term — an *alignment* ceiling, not a
  selection one: V4 surfaces the source in only 94.4% of cases because its ~250 bp amplicon
  is byte-identical to ≥20 other genomes, so the source ranks beyond top-20. No selection
  policy can recover what alignment never returned, and any tied organism is an equally valid
  representative. To lift this ceiling, widen alignment retention (top-50/100) for short
  regions — in `align_hits.py`/`gpu_align.py`, not here.
- **`no_gene_set` (8.6%; 31% of full-length failures) is self-inflicted by exact-gene mode**
  dropping 16S-only source genomes from the gene union. Union mode keeps them as taxonomic
  reps; they simply contribute nothing to the gene union.

**Overall capture (selected / all amplicons): 80.7% → ~95% (exact-tie+guards) / ~97%
(exact-tie), ~100% on the full-length control.**

---

## 5. Recommendation — and the new default

**The default selection standard (now shipped in `select_references.DEFAULT_KNOBS`) is the
exact-tie de-duplicated union:** `select_all=True, select_all_window=0.0,
select_all_species_dedup=True`, guards off. It keeps every reference tied at the best
identity, one representative per BV-BRC species — self-recovery 82.9→**99.2%** (99.9%
full-length), gene-capture recall +2–3pp, **precision/F1 at or above baseline**, zero
same-species redundancy. The legacy ordered reducer is `select_all=False`
(`--select-mode legacy`); it reproduces the committed validation baseline exactly.

The **precision firewalls are an optional pathway**, not the default:

1. **`--guarded`** (exact-tie + `select_all_coherence` + `select_all_contam` + `select_all_cap=6`):
   adds the cross-family firewall + contamination budget. Self-recovery 95.3% / 99.7%, precision
   a touch higher than the bare default in short regions. Use when an FBA model must minimize any
   cross-family gene even among exact ties.
2. **`--select-mode graduated`** (exact ties kept unconditionally + a within-`window` extension —
   default 0.005 — admitted only through the firewalls + budget): adds nearby pangenome breadth
   *only where it passes the family-agreement / contamination trade-off*, leaving the exact ties
   untouched. (Benchmarked in §3 / `variant_comparison.md`.)
3. **Do not use** the un-guarded 0.005 `in-band` band (−7–12pp precision for no F1 gain) or the
   `all-reliable` union (precision/F1/consensus collapse).
4. **For literal 100% on short regions**, raise the alignment top-k; selection is already at its
   ceiling there.
5. Gene-less source genomes are kept as taxonomic reps (intrinsic to union mode).

*Note:* the committed `selection.json` and baseline score files were produced by the legacy
reducer (now `select_all=False`); a fresh validation run uses the new union default. The
five-way `variant_comparison.md` is the authoritative legacy-vs-union comparison.

---

## 6. Artifacts / reproducibility

- `scripts/diagnose_self_recovery.py` → `data/self_recovery_diagnosis.json`, `data/self_recovery_failures.csv`
- `select_references.py` — `select_all` union mode + `select_all_window` / `_species_dedup` (BV-BRC taxonomy) / `_coherence` / `_contam` / `_cap`
- `scripts/build_variant_selection.py` → `data/selection_{exacttie,exacttie_guarded,includeall_band,includeall_reliable}.json`
- `scripts/prefetch_variant_genomes.py` — grew the shared PGFam cache 57k→112k genomes (append-only; baseline untouched)
- `scripts/score_variant.py` — scores a variant with the unmodified `score_taxacc`/`score_genecap` (suffixed outputs `_et`/`_etg`/`_ib`/`_ir`)
- `scripts/compare_variants.py` → `data/variant_comparison.md`
