# Exclusion stringency: tau = 99 ANI on the existing 10,000-genome panel

> **Re-score 2026-11-16.** The numbers below were re-derived from the artefacts re-scored on the repaired
> BV-BRC PGFam cache (2026-11-13 repair; re-score runs R2.2/R2.3; the depth-20 amplicon cells the pairing reads
> were re-scored in R1.4). Headline changes: resolve rate at 100% completeness tau 95 -> 99 is now
> 0.313 -> 0.372 (was 0.307 -> 0.366); conditional F1 0.8066 -> 0.8220 (was 0.8068 -> 0.8230); paired against
> full-length 16S at tau 99, n 1,815 -> 1,857 and mean delta +0.130 -> +0.131; MAG population F1 at tau 99
> 0.3008 -> 0.3040. Unchanged: genus and family accuracy, species 0.000 and the 1,285 references retained in the
> 95-99 ANI band (selections did not change). Nothing is left historical.

A pure re-score of the existing hits (no new alignment or synthesis): self-exclusion
drops the source genome, its taxon, and every hit at >= 99 ANI instead of >= 95.
All 20 batches, 199,940 cells, 0% contamination.

**tau = 99 is a WEAKER exclusion, not a stricter one.** 95 ANI is the species
boundary, so tau = 99 retains references that are conspecific with the source by the
ANI definition. It is therefore NOT matched to the amplicon arm's tau = 0.987; the
stringency-matched amplicon comparator would be exclude_0.995, which the shipped
gene_capture_cells.csv does not carry. Every MAG-vs-16S delta below is measured
against a STRICTER amplicon exclusion than the MAG one.

## Effect on the MAG arm

| completeness | resolve 95 -> 99 | cond. F1 95 -> 99 | genus 95 -> 99 | family 95 -> 99 |
|---|---|---|---|---|
| 100% | 0.313 -> **0.372** | 0.8066 -> **0.8220** | 0.910 -> 0.895 | 0.960 -> 0.951 |
| 90%  | 0.314 -> 0.374 | 0.8067 -> 0.8224 | 0.910 -> 0.893 | 0.957 -> 0.949 |
| 70%  | 0.311 -> 0.372 | 0.8070 -> 0.8221 | 0.910 -> 0.894 | 0.959 -> 0.951 |
| 50%  | 0.312 -> 0.372 | 0.8063 -> 0.8232 | 0.911 -> 0.898 | 0.959 -> 0.950 |
| 30%  | 0.313 -> 0.375 | 0.8049 -> 0.8230 | 0.912 -> 0.896 | 0.958 -> 0.951 |

Coverage rises +0.059 (+19% relative) and conditional gene capture +0.015, because
1,285 references in the 95-99 ANI band are retained that tau = 95 removed. Genus
accuracy falls slightly (0.910 -> 0.895): the newly admitted neighbours sit a little
further off taxonomically than those that already qualified.

**Completeness remains a null axis at tau = 99**, flat across all five tiers exactly
as at tau = 95.

## Paired against full-length 16S

| completeness | tau=95 (n / MAG / 16S / delta) | tau=99 (n / MAG / 16S / delta) |
|---|---|---|
| 100% | 1541 / 0.8142 / 0.7097 / **+0.104** | 1857 / 0.8288 / 0.6983 / **+0.131** |
| 50%  | 1545 / 0.8127 / 0.7092 / +0.104 | 1863 / 0.8289 / 0.6983 / +0.131 |
| 30%  | 1541 / 0.8114 / 0.7096 / +0.102 | 1882 / 0.8283 / 0.6973 / +0.131 |

Population F1 (abstention scored 0) at 100% completeness: MAG 0.2502 -> 0.3040
against the amplicon's unchanged 0.5388, so the population gap narrows from -0.289
to -0.235 but does not close. Part of the widened conditional delta is the
stringency mismatch noted above, not a genuine gain.

## Species accuracy stays at 0.000

Even at tau = 99, which by construction retains ANI-conspecific references, not one
selected reference is the source's named NCBI species. The panel is
species-dereplicated -- no two panel genomes share a species name -- so every
retained 95-99 ANI neighbour is conspecific by ANI and a different named species.
ANI-95 and NCBI species nomenclature simply disagree for these pairs, which is worth
stating directly rather than reporting as a 0.
