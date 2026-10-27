# Manuscript figures rendered from the depth-100 exact-tie union trees

Scripts that render the prFBA_paper figures whose data depend on the selection rule. They
read the union trees (`region_validation_d100`, `mag_validation/out/*_d100.json`), the
rule the paper uses throughout since 2026-10-27, and write to `out/` (git-ignored). Copy
the outputs to `prFBA_paper/images/prfba/` under these names:

| Script | Output in `out/` | Paper file |
|---|---|---|
| `rv/render_accuracy.py` | `F1_accuracy_heatmap.pdf` | `F1_accuracy_heatmap.pdf` (fig:accuracy) |
| `rv/render_accuracy.py` | `F11_self_inclusion_exclusion_delta.pdf` | `F11_self_inclusion_exclusion_delta.pdf` (fig:leakage) |
| `rv/render_genecap.py` | `F4b_percent_genecapture_by_region.pdf` | `F4b_percent_genecapture_by_region.pdf` (fig:genecap) |
| `rv/run_exclude_alt.py` | `F2_perrank_concordance_exclude_self_ALT.pdf` | `F2_perrank_concordance.pdf` (fig:perrank) |
| `mag/plot_results_d100.py` | `F_mag_completeness.pdf` | `F14_mag_completeness.pdf` (fig:magcomp) |
| `mag/plot_crossover_d100.py` | `F_mag_crossover.pdf` | `F15_mag_crossover.pdf` (fig:magcross) |

Run each from its own directory, for example
`cd rv && ~/Documents/py_venv/bin/python render_accuracy.py`.

The `rv/` scripts are copies of `region_validation_d100/scripts/` with `_common.py` pointed
at the union tree and at `out/`. `render_genecap.py` has two layout patches: y headroom and
the legend placement. The `mag/` scripts are copies of `mag_validation/plot_results.py` and
`plot_crossover.py` with the union inputs and the union rule text.

The change list behind the text is `region_validation_d100/UNION_THROUGHOUT_CHANGES.md`.
