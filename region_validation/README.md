# 16S region validation

In-silico validation of the ASV → reference-genome mapping + synthetic-genome pipeline.
10,000 deterministically stratified BV-BRC genomes are run through in-silico PCR for eight
16S regions, each amplicon passed through the **unaltered** production pipeline, and scored
against the source organism's true NCBI lineage and PGFam gene set under two regimes:

- **self-inclusion** — source organism left in the database (realistic experimental upper bound)
- **self-exclusion** (τ = 0.987) — source + near-identical neighbors removed (novel-organism generalization)

The method, design rationale, and full results are written up in
[`16S_pipeline_method_paper.docx`](16S_pipeline_method_paper.docx) (source: `.md`); the design
spec is [`DESIGN.md`](DESIGN.md).

## Layout

| path | contents |
|---|---|
| `scripts/` | pipeline driver (`run_validation.sh`) + every stage (select → PCR → align → select_references → score → render) |
| `data/` | intermediates + result-summary JSONs (`taxacc_region_summary.json`, `gene_capture_by_region.json`) + figure-input CSVs |
| `figures/` | 22 rendered figures (PNG + PDF) |
| `data_pilot_500/`, `figures_pilot_500/` | 500-genome pilot (developmental checkpoint) |

## Run / resume

```
cd scripts && bash run_validation.sh full      # fully resumable; add --force to rebuild from scratch
```

## Large files (> 100 MB, git-ignored)

Regenerable by re-running `scripts/run_validation.sh`; absent from a clone (also listed in `../README.md`):

| file | size |
|---|---|
| `data/hits/asv_top20_alignment_hits.json` | 774 MB |
| `data/hits/raw/asv_top20_alignment_hits.json` | 690 MB |
| `data/selection.json` | 666 MB |
| `data/taxacc_per_record.jsonl` | 270 MB |

Everything else here (result summaries, CSVs, figures, scripts, manuscript, pilot archive) is < 100 MB and tracked.
