# Re-score of 2026-11: run records

Records of the full re-score (2026-11-15 to 2026-11-17), which recomputed every manuscript and repo value with the
current code on the repaired, BV-BRC-only gene-family data. Copied here from the work directory
`~/Documents/prFBA_rescore_2026-11` so they travel with the repo.

- `plan/`: the plan (`PLAN.md`), decisions D1-D11 and the claim inventory.
- `logs/`: one report per run (`RUN1_REPORT.md` ... `RUN5B_REPORT.md`, `RUN5C_FIXES.md`), the independent
  verifiers' reports (`verify_r5/`), commit messages and run logs.
- `tools/`: helper scripts the runs used outside the repo's own code.
- `TF_RUNBOOK.md`: how the second machine ran prFBA code in a path-isolating namespace.

What was not kept: the multi-GB pre-repair controls and snapshots (`control/`, `pre_repair/`, deleted 2026-11-19).
The pre-repair state stays recoverable: committed results at tag `rescore-2026-11-start`, and the old cache values in
`bvbrc_cache/CORRECTIONS_2026-11-13.json` (`backup_old_values`). The work directory keeps the small hash stamps
(`stamps/`).
