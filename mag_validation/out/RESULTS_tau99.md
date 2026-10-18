# Exclusion stringency: tau = 99 ANI on the existing 10,000-genome panel

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
| 100% | 0.307 -> **0.366** | 0.8068 -> **0.8230** | 0.910 -> 0.895 | 0.960 -> 0.951 |
| 90%  | 0.308 -> 0.367 | 0.8070 -> 0.8233 | 0.910 -> 0.893 | 0.957 -> 0.949 |
| 70%  | 0.306 -> 0.366 | 0.8071 -> 0.8228 | 0.910 -> 0.894 | 0.959 -> 0.951 |
| 50%  | 0.307 -> 0.366 | 0.8067 -> 0.8240 | 0.911 -> 0.898 | 0.959 -> 0.950 |
| 30%  | 0.307 -> 0.368 | 0.8052 -> 0.8239 | 0.912 -> 0.896 | 0.958 -> 0.951 |

Coverage rises +0.058 (+19% relative) and conditional gene capture +0.016, because
1,285 references in the 95-99 ANI band are retained that tau = 95 removed. Genus
accuracy falls slightly (0.910 -> 0.895): the newly admitted neighbours sit a little
further off taxonomically than those that already qualified.

**Completeness remains a null axis at tau = 99**, flat across all five tiers exactly
as at tau = 95.

## Paired against full-length 16S

| completeness | tau=95 (n / MAG / 16S / delta) | tau=99 (n / MAG / 16S / delta) |
|---|---|---|
| 100% | 1510 / 0.8136 / 0.7092 / **+0.104** | 1815 / 0.8295 / 0.6991 / **+0.130** |
| 50%  | 1514 / 0.8121 / 0.7087 / +0.103 | 1821 / 0.8291 / 0.6990 / +0.130 |
| 30%  | 1510 / 0.8107 / 0.7092 / +0.102 | 1840 / 0.8287 / 0.6980 / +0.131 |

Population F1 (abstention scored 0) at 100% completeness: MAG 0.2477 -> 0.3008
against the amplicon's unchanged 0.5346, so the population gap narrows from -0.287
to -0.234 but does not close. Part of the widened conditional delta is the
stringency mismatch noted above, not a genuine gain.

## Species accuracy stays at 0.000

Even at tau = 99, which by construction retains ANI-conspecific references, not one
selected reference is the source's named NCBI species. The panel is
species-dereplicated -- no two panel genomes share a species name -- so every
retained 95-99 ANI neighbour is conspecific by ANI and a different named species.
ANI-95 and NCBI species nomenclature simply disagree for these pairs, which is worth
stating directly rather than reporting as a 0.
