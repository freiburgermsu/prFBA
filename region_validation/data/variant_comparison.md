# Include-all (“union”) selection vs baseline — 1e4 self-recovery benchmark

Three selection policies, same alignment hits, same metric code:
- **baseline** = ordered reducer pipeline (family-consensus → dedup → marginal-gain → tier-cap).
- **in-band union** = keep every reliable hit within 0.005 identity of the best (remove order stochasticity among co-optimal hits).
- **all-reliable union** = keep every reliable hit down to the family floor.

## 1. Self-recovery — P(source selected | source in top-20), include scenario

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 96.3 | 99.9 | 99.9 | 99.7 | 100.0 |
| V1-V2 | 86.0 | 99.4 | 99.4 | 96.3 | 100.0 |
| V1-V3 | 88.6 | 99.5 | 99.5 | 97.6 | 100.0 |
| V3-V4 | 81.1 | 98.9 | 98.9 | 94.6 | 100.0 |
| V4 | 75.7 | 98.8 | 98.8 | 92.2 | 100.0 |
| V4-V5 | 78.6 | 99.0 | 99.0 | 93.8 | 100.0 |
| V4-V5_944R | 79.0 | 99.0 | 99.0 | 93.9 | 100.0 |
| V7-V9 | 81.9 | 99.1 | 99.1 | 95.3 | 100.0 |
| ALL | 82.9 | 99.2 | 99.2 | 95.3 | 100.0 |

_self-as-anchor | source in top-20 (include):_

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 96.1 | 99.1 | 99.1 | 99.1 | 98.9 |
| V1-V2 | 82.2 | 82.1 | 82.1 | 82.1 | 83.5 |
| V1-V3 | 86.0 | 87.4 | 87.4 | 87.4 | 88.2 |
| V3-V4 | 76.1 | 76.9 | 76.9 | 76.9 | 78.1 |
| V4 | 67.8 | 68.1 | 68.1 | 68.1 | 69.6 |
| V4-V5 | 72.3 | 73.3 | 73.3 | 73.3 | 74.6 |
| V4-V5_944R | 72.7 | 73.7 | 73.7 | 73.7 | 74.9 |
| V7-V9 | 77.0 | 78.5 | 78.5 | 78.5 | 79.6 |

## 2. Taxonomic concordance (correct | called)

### Anchor — Genus, include

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 98.6 | 99.6 | 99.6 | 99.6 | 99.7 |
| V1-V2 | 96.7 | 95.1 | 95.1 | 95.1 | 96.5 |
| V1-V3 | 97.2 | 97.0 | 97.0 | 97.0 | 97.7 |
| V3-V4 | 95.3 | 93.2 | 93.2 | 93.2 | 95.3 |
| V4 | 93.6 | 90.8 | 90.8 | 90.8 | 93.1 |
| V4-V5 | 95.0 | 93.0 | 93.0 | 93.0 | 94.6 |
| V4-V5_944R | 95.0 | 93.0 | 93.0 | 93.0 | 94.5 |
| V7-V9 | 95.9 | 94.6 | 94.6 | 94.6 | 96.0 |

### Anchor — Species, include

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 97.1 | 99.2 | 99.2 | 99.2 | 99.1 |
| V1-V2 | 83.1 | 82.2 | 82.2 | 82.2 | 83.7 |
| V1-V3 | 87.2 | 87.7 | 87.7 | 87.7 | 88.6 |
| V3-V4 | 75.3 | 75.8 | 75.8 | 75.8 | 76.9 |
| V4 | 65.7 | 65.5 | 65.5 | 65.5 | 66.8 |
| V4-V5 | 71.1 | 71.8 | 71.8 | 71.8 | 72.7 |
| V4-V5_944R | 71.5 | 72.0 | 72.0 | 72.0 | 72.9 |
| V7-V9 | 76.8 | 77.8 | 77.8 | 77.8 | 78.6 |

### Consensus — Genus, include

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 98.7 | 99.9 | 99.2 | 99.8 | 98.0 |
| V1-V2 | 97.2 | 98.5 | 98.2 | 96.4 | 97.2 |
| V1-V3 | 97.7 | 99.2 | 98.7 | 98.0 | 97.5 |
| V3-V4 | 96.0 | 97.8 | 96.8 | 94.8 | 95.3 |
| V4 | 95.1 | 96.2 | 95.6 | 92.9 | 94.2 |
| V4-V5 | 96.1 | 97.3 | 96.2 | 94.8 | 94.9 |
| V4-V5_944R | 96.3 | 97.3 | 96.2 | 94.8 | 94.8 |
| V7-V9 | 96.8 | 98.2 | 97.0 | 96.1 | 95.4 |

### Consensus — Genus, exclude_0.987 (honest)

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 75.1 | 72.6 | 73.6 | 72.3 | 74.7 |
| V1-V2 | 73.2 | 72.5 | 72.8 | 71.2 | 73.1 |
| V1-V3 | 73.6 | 72.1 | 72.7 | 71.1 | 73.2 |
| V3-V4 | 66.1 | 63.9 | 64.9 | 62.8 | 65.7 |
| V4 | 62.1 | 61.7 | 61.8 | 59.7 | 61.6 |
| V4-V5 | 65.0 | 64.2 | 64.0 | 62.4 | 64.4 |
| V4-V5_944R | 65.5 | 64.9 | 65.1 | 63.2 | 64.7 |
| V7-V9 | 66.8 | 66.6 | 66.2 | 64.7 | 66.7 |

### Anchor — Genus, exclude_0.987 (honest)

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 73.0 | 72.1 | 72.1 | 72.1 | 72.0 |
| V1-V2 | 72.0 | 70.6 | 70.6 | 70.6 | 71.1 |
| V1-V3 | 71.7 | 70.7 | 70.7 | 70.7 | 70.7 |
| V3-V4 | 63.9 | 62.0 | 61.9 | 62.0 | 62.5 |
| V4 | 59.9 | 58.4 | 58.3 | 58.4 | 59.1 |
| V4-V5 | 62.6 | 61.3 | 61.2 | 61.3 | 61.1 |
| V4-V5_944R | 63.1 | 62.0 | 61.9 | 62.0 | 62.3 |
| V7-V9 | 65.3 | 64.0 | 64.0 | 64.0 | 64.2 |

## 3. Gene capture (macro, PGFam sets)

### Recall — include (% gene capture)

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 99.4 | 100.0 | 100.0 | 99.9 | 100.0 |
| V1-V2 | 97.5 | 99.4 | 99.4 | 98.1 | 99.8 |
| V1-V3 | 98.1 | 99.6 | 99.6 | 98.9 | 99.9 |
| V3-V4 | 96.5 | 99.1 | 99.1 | 97.3 | 99.7 |
| V4 | 95.1 | 98.6 | 98.7 | 96.5 | 99.4 |
| V4-V5 | 96.0 | 98.9 | 99.0 | 97.5 | 99.6 |
| V4-V5_944R | 96.0 | 98.9 | 99.0 | 97.4 | 99.5 |
| V7-V9 | 96.8 | 99.2 | 99.2 | 97.8 | 99.7 |

### Precision — include

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 96.6 | 99.7 | 93.9 | 99.7 | 89.7 |
| V1-V2 | 93.3 | 91.9 | 89.0 | 93.0 | 86.2 |
| V1-V3 | 94.3 | 94.6 | 90.5 | 95.4 | 87.0 |
| V3-V4 | 89.4 | 87.5 | 82.3 | 89.3 | 78.2 |
| V4 | 85.9 | 81.1 | 77.2 | 83.8 | 73.8 |
| V4-V5 | 87.4 | 85.2 | 79.3 | 87.4 | 75.3 |
| V4-V5_944R | 87.4 | 85.6 | 79.7 | 87.5 | 75.6 |
| V7-V9 | 89.2 | 88.7 | 82.3 | 90.4 | 78.0 |

### F1 — include

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 97.7 | 99.8 | 96.3 | 99.8 | 93.2 |
| V1-V2 | 94.9 | 94.4 | 92.5 | 94.9 | 90.6 |
| V1-V3 | 95.7 | 96.3 | 93.7 | 96.7 | 91.2 |
| V3-V4 | 92.2 | 91.2 | 87.8 | 92.2 | 84.7 |
| V4 | 89.4 | 86.5 | 83.9 | 88.3 | 81.2 |
| V4-V5 | 90.8 | 89.6 | 85.7 | 91.0 | 82.5 |
| V4-V5_944R | 90.8 | 89.9 | 86.1 | 91.2 | 82.9 |
| V7-V9 | 92.2 | 92.1 | 88.0 | 93.1 | 84.6 |

### Jaccard — include

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 96.4 | 99.7 | 93.9 | 99.7 | 89.7 |
| V1-V2 | 92.1 | 91.6 | 88.7 | 92.4 | 86.2 |
| V1-V3 | 93.4 | 94.4 | 90.3 | 94.9 | 87.0 |
| V3-V4 | 87.7 | 87.1 | 81.9 | 88.5 | 78.2 |
| V4 | 83.7 | 80.6 | 76.8 | 82.8 | 73.7 |
| V4-V5 | 85.6 | 84.8 | 79.0 | 86.6 | 75.2 |
| V4-V5_944R | 85.6 | 85.2 | 79.4 | 86.8 | 75.6 |
| V7-V9 | 87.7 | 88.3 | 82.0 | 89.6 | 78.0 |

### Recall — exclude_0.987 (honest generalization)

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 72.2 | 68.3 | 73.4 | 68.1 | 75.4 |
| V1-V2 | 70.7 | 69.1 | 72.0 | 68.6 | 73.6 |
| V1-V3 | 71.2 | 68.8 | 72.7 | 68.4 | 74.4 |
| V3-V4 | 69.0 | 66.0 | 70.7 | 65.4 | 72.9 |
| V4 | 67.9 | 66.2 | 69.9 | 65.2 | 71.8 |
| V4-V5 | 68.6 | 65.5 | 70.4 | 64.8 | 72.7 |
| V4-V5_944R | 69.1 | 65.7 | 70.8 | 65.0 | 73.0 |
| V7-V9 | 69.6 | 66.1 | 71.2 | 65.4 | 73.4 |

### Precision — exclude_0.987

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 61.8 | 66.3 | 59.5 | 66.4 | 56.4 |
| V1-V2 | 63.6 | 64.5 | 61.1 | 65.2 | 58.9 |
| V1-V3 | 63.2 | 65.3 | 60.5 | 65.7 | 57.8 |
| V3-V4 | 56.8 | 59.3 | 53.9 | 60.2 | 50.7 |
| V4 | 55.7 | 56.5 | 52.4 | 57.7 | 49.6 |
| V4-V5 | 55.9 | 58.3 | 52.5 | 59.2 | 49.3 |
| V4-V5_944R | 56.2 | 58.8 | 52.9 | 59.6 | 49.7 |
| V7-V9 | 56.5 | 59.5 | 53.6 | 60.1 | 50.5 |

### F1 — exclude_0.987

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 64.5 | 65.6 | 63.4 | 65.7 | 61.6 |
| V1-V2 | 65.0 | 64.7 | 63.8 | 65.0 | 62.7 |
| V1-V3 | 64.9 | 65.1 | 63.6 | 65.3 | 62.1 |
| V3-V4 | 59.8 | 60.0 | 58.1 | 60.5 | 56.3 |
| V4 | 58.5 | 58.0 | 56.7 | 58.7 | 55.0 |
| V4-V5 | 59.0 | 59.1 | 57.1 | 59.6 | 55.1 |
| V4-V5_944R | 59.4 | 59.5 | 57.5 | 59.9 | 55.6 |
| V7-V9 | 59.9 | 60.2 | 58.3 | 60.4 | 56.4 |

_OK (scorable) cells per region, include:_

### n OK cells — include

| region | legacy (baseline) | exact-tie (DEFAULT) | graduated | exact-tie+guards | in-band 0.005 |
|---|---|---|---|---|---|
| FullLength16S | 721400.0 | 721500.0 | 721500.0 | 721500.0 | 721500.0 |
| V1-V2 | 744200.0 | 744500.0 | 744500.0 | 744500.0 | 744900.0 |
| V1-V3 | 750200.0 | 750400.0 | 750400.0 | 750400.0 | 750700.0 |
| V3-V4 | 824900.0 | 824700.0 | 825000.0 | 824600.0 | 825200.0 |
| V4 | 943600.0 | 944000.0 | 944200.0 | 943800.0 | 944900.0 |
| V4-V5 | 939900.0 | 940300.0 | 940400.0 | 940300.0 | 941000.0 |
| V4-V5_944R | 899400.0 | 899900.0 | 899900.0 | 899700.0 | 900400.0 |
| V7-V9 | 874200.0 | 874200.0 | 874500.0 | 874200.0 | 874700.0 |
