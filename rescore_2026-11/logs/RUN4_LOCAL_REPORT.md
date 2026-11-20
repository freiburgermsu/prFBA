# Run 4 report, FM (local) part: wobj, phenotype passes, score pass 2 (R4.1, R4.2, R4.3 wobj), 2026-11-16

Executor: the Run 4 local agent. The TF agent ran calibrate, Stage E, the Stage O weights, Stage O, FVA, the wrapper
and verification on TF; none of that is reported here.

Plan sections used: `plan/PLAN.md` §0–4 (Run 4), §5 (R4 claims) and §8; `plan/pipelines_shared_ref.md` §3D/3F/6;
`RUN3_REPORT.md` "What Run 4 needs"; `TF_RUNBOOK.md` §4 (D11) and the wobj barrier note.

Where things are:
- progress log with times: `logs/RUN4_LOCAL_PROGRESS.md`;
- every memory decision with its arithmetic: `logs/r4/r4_mem_decisions.log`;
- controls and checks: `control/r4/`;
- unit logs: `logs/r4/`, plus the repo logs `SR/logs/{wobj_run,pheno_c13_atp,pheno_c13_gf,score_genus_c13,rescore_compare_*}.log`;
- unit wall times and cgroup peaks: `stamps/units.jsonl`;
- moved-aside outputs: `pre_repair/region_validation/shared_reference/{data/weighted_objective,data/phenotype_{atp,gf}_c13,data/phenotype_summary.json,logs}`.

**Outcome.**
- All three local stages ran on FM, and every control passed. Pass 2 is valid (exit 0; C1-C15 all pass; C11/C12
  bound to the new gapfilled summary). The genus summary is no longer stamped INVALID, and its arms equal pass 1's,
  so **Run 3's genus numbers are final**.
- The cache was never written: `stamp_state.py --compare R3.final` reads CACHE UNCHANGED at the start
  (`stamps/R4L.0.json`) and at the end (`stamps/R4L.final.json`).
- Commits are path-limited and not pushed: **2520f4f** (wobj), **5fac436** (phenotype + default summary),
  **bb41947** (genus pass 2).
- **Headline for Run 5.** The phenotype numbers move for a reason other than the repaired data itself.
  - The new data change the ΔF1 of 200 phenotype candidates. That reshuffles the stratified samples: 558 of the
    1,000 ATP genomes and 138 of the 200 gapfilled genomes are new.
  - On the new sample, the **prFBA-vs-nn1 ATP advantage (−0.146 media, p = 0.017, q = 0.028; L635-637, L1186) is
    no longer significant (−0.030, p = 0.41, q = 0.51)**.
  - McNemar now favours PICRUSt2 on **4** of 22 media, not 2.
  - The PICRUSt2 null stays null (0.13 → 0.069 media, p 0.49 → 0.79).
  - Every row of a genome kept in both samples is identical unless its own inputs changed (0 unattributed), so this
    is sample composition, not code or nondeterminism. The text in §2 and §5 must be rewritten, not just
    refreshed.
- wobj moves only in the 4th decimal: 2 values change at print (0.0120 → 0.0121, 27.9 → 27.8).
- **Side effect on the user's other work.** Our phenotype passes and pass 2 held the salterns rollout's prep queue
  for about 5.75 h in total (§0, §4).

---------------------------------------------------------------------------------------------------------------

## 0. Start state

- `stamp_state.py --compare stamps/R3.final.json`: **CACHE UNCHANGED**, catalogue current (`stamps/R4L.0.json`,
  `logs/r4/stamp_R4L0.log`).
- The d20 baseline fails on exactly the 17 files plan R5.5 expects. Re-stamping it is Run 5's job.
- Code: weighted_objective_test.py, phenotype_fold.py and phenotype_arm.py are unchanged since c91cf19. Their
  sha256s equal the A-C13 binding (`defa403b`, `0d3342ee`).
- Pins verified clean:
  - ModelSEEDpy-c2e31a0, ModelSEEDpy-d884a50;
  - ModelSEEDTemplates-c12run 21a8084, ModelSEEDTemplates-main2294299;
  - the ModelSEEDDatabase-e8446e3e worktree.

  No main checkout was touched, and no new worktree was needed.

**A competing workload on FM.** The salterns R2A assembly rollout (`~/Documents/salterns_rollout_local`, approved
by Andrew on 11-16, started 11:10) runs in a scope capped at 220G.
- Its prep tasks (merge/tadpole) reserve 150 GB each.
- It launches a task only when MemAvailable − its outstanding reserve − the new reserve ≥ 40 GB, so a prep task
  needs MemAvailable ≥ 190.
- On this machine, idle MemAvailable is about 219 GB. So **any job of ours holding more than about 29 GB stops
  salterns from starting its next prep task**, whatever its MemoryMax.
- There is no safe "window" between salterns tasks: it relaunches within one 20 s poll, faster than a phenotype
  parent ramps up.

The memory rule was applied as follows: budget = MemAvailable − salterns' outstanding reserve − 15 GB margin −
the remaining growth of our running units up to their caps. Caps follow the rule:
- gapfilled pass: 15 + 10 GB per worker;
- ATP pass: 15 + 5 GB per worker.

Sequence:
1. The ATP pass started first at the size the strict budget allowed during a salterns merge. That holding blocked
   salterns' next prep task.
2. The gapfilled pass started at 12:22, once salterns was idle and held for memory.
3. Consequence: **salterns' prep was held from 12:09:31 to 16:32 (4.4 h)**, and again while pass 2 runs
   (§3).
4. Nothing was frozen or killed. Salterns' mem_guard never fired, and MemAvailable never went below 88 GB.

---------------------------------------------------------------------------------------------------------------

## 1. R4.3 weighted-objective test (FM; D11 pin ModelSEEDDatabase e8446e3e)

**Move-aside.** These went to `pre_repair` (sha256 in `control/r4/wobj_moveaside.sha256`): `population.json`,
`per_amplicon.jsonl`, `summary.json`, `recomputed_contrasts.json`, the 4 figures and `logs/wobj_run.log`.
`msprob_investigation/` and `patches/` stayed in place: they are tracked records, not outputs of this run.

**Commands.**
- `launch_wobj.sh` was read, not run, because it calls systemd-run itself. Its two python lines were run under
  `U.sh` with the pins explicit:
  - `PYTHONPATH=ModelSEEDpy-c2e31a0`, `--expect-modelseedpy c2e31a0`;
  - `--msd ~/Documents/ModelSEEDDatabase-e8446e3e`;
  - templates: the script default, c12run 21a8084.
- Units:
  - `rescore-r4-wobj-pop`: 12G, `--population-only`;
  - `rescore-r4-wobj-run`: 20G, `--workers 8 --maxtasksperchild 2`;
  - then `wobj_results.py` with `SOURCE_DATE_EPOCH=1794700800`.

| step | wall | cgroup peak | result |
|---|---|---|---|
| population | 83 s | 4.0 GiB | `population.json` **5d13cf1e…**; marker `logs/R4_WOBJ_POPULATION_DONE` at 11:56 (TF's Stage E commit bb5dc5b cites it) |
| control (HEAD, pre-repair inputs, 3 amplicons) | 296 s | 7.8 GiB | PASS (below) |
| main run | 7,094 s (11:56–13:54) | 12.3 GiB | 360/360 ok |
| wobj_results.py | 2 s | 0.1 GiB | agreement 645 cells, 0 differ > 1e-12 |

**Population checks** (`control/r4/wobj_population_check.json`):
- `spec.inputs` carry the new cells **d57709e3** and panel EC **8eb296b1**; panel, EC universe and fold plan are
  unchanged.
- The settings equal the 10-29 spec, and the sample is 100/100/100/60.
- Same 54,243 candidates and the same 360 keys as 10-29. 358 jobs are identical; 2 changed, both through repaired
  genomes:
  - 2975960.3: truth ECs 810 → 828;
  - members 1971901.6 and 888845.10 of 3069120.67: 608 → 640 and 578 → 638 ECs.

**Control.** On the pre-repair inputs (cells ffc9d3a7, panel EC 185ca60f, the 10-14 predictions), today's code with
the same pins reproduces:
- the committed population (jobs, spec without paths, and every count identical);
- the 3 subset rows (n2, n3, n4+; 8 arms each), bitwise in every field except timing and memory, **barrier
  re-solves included**.

Evidence: `control/r4/wobj_ctl_head_vs_committed.json`. The 11-15 FM rerun is also bitwise
(`wobj_ctl_fm1115_vs_committed.json`).

**Attribution** (`tools/r4_wobj_rowcheck.py` → `control/r4/wobj_rowcheck.json`):
- the 358 rows of unchanged jobs are identical to the 10-29 rows in every field except timing and memory;
- only the 2 changed jobs' rows differ;
- 0 unattributed.

**Old → new.**
- Verdicts unchanged:
  - H1 "detectable but below SESOI";
  - H1s specific to p;
  - H2 reduces degeneracy, not specific to p (0.2139; 1e-7 face 0.2165).
- `collect.py --only` the three wobj items (`control/r4/collect_r4_wobj.{json,md}`): 30 values, every old value
  equals the print, and 2 change at print precision:
  - L307-311 "beat shuffled" **0.0120 → 0.0121** (−0.012046 → −0.012079);
  - L307-311 "shuffled free reactions" **27.9 → 27.8**.
  - released_check is still 360/360.

Commit **2520f4f**.

---------------------------------------------------------------------------------------------------------------

## 2. R4.1 phenotype passes, A-C13 configuration

**Configuration.** The exact argv of the 10-29 `c13-atp` and `c13-gf` units, recovered from the user journal (the
10-29 run logged no argv anywhere else):
- `-W ignore::SyntaxWarning`, cwd SR;
- `--templates ModelSEEDTemplates-c12run/templates/v7.0`, `--integrate-reversals`;
- `--expect-modelseedpy c2e31a0` with PYTHONPATH c2e31a0;
- `--summary` into the pass directory;
- ATP: `--n-genomes 1000 --readouts atp,structure --maxtasksperchild 8`;
- gapfilled: `--n-genomes 200 --readouts atp,growth,essentiality,structure --maxtasksperchild 3`.

Two changes from that argv:
- **`--msd ~/Documents/ModelSEEDDatabase-e8446e3e`.** The 10-29 meta records only the default path, not a commit.
  The main checkout's reflog puts it at e8446e3e from 10-14 18:53 to 10-30 16:42, which covers both 10-29 passes,
  so the commit is inferred from that window (and `C12_SCOPE.md`). The existing clean worktree was used.
- **`--workers`.** ATP ran at 3, then 8 after a stop and resume; gapfilled ran at 12. On 10-29 they were 10 and 16.
  The counts were set by the memory budget (§0). Rows are per genome and seeded (worker-independence record
  84/84, 42/42; `--maxtasksperchild` unchanged).

The A-C13 directories and logs were moved to `pre_repair` first, so both passes wrote fresh canonical directories.

**The ATP stop and resume** (`control/r4/atp_restart/`).
- When wobj ended (13:54), the ATP unit was stopped at 225 genomes. The whole last genome block (1896223.3, 44 rows)
  was dropped as a precaution against a partial write; no line was truncated.
- The unit was relaunched with 8 workers (`rescore-r4-pheno-atp2`, 55G). phenotype_fold's own resume read
  "1000 genomes selected, 224 already done, 776 to build", and its final pass rewrites `per_genome.jsonl` sorted.
- The identity control re-ran on the resumed run's first genome. `identity_control.json` is byte-identical to the
  first run's and to 10-29's.
- This saved about 4 h on the ATP pass and on salterns' hold. The bit-for-bit row checks below cover it.

| pass | unit(s) | wall | cgroup peak / MemoryMax | rows | gated | identity |
|---|---|---|---|---|---|---|
| ATP | rescore-r4-pheno-atp (3 workers, 12:05–13:54) + rescore-r4-pheno-atp2 (8 workers) | 1 h 49 m + 2 h 17 m | 28.9 / 30G; 49.8 / 55G | 43,174 (1,000 genomes) | 0 | pass |
| gapfilled | rescore-r4-pheno-gf (12 workers) | 4 h 10 m (12:22–16:32) | 88.5 / 135G | 8,660 (200 genomes) | **0 (0 truth)** | pass (atp and gapfill identical) |

**Checks.**
- `--expect-modelseedpy` passed.
- `meta.binding`: cells d57709e3, functions 9463ed94 (new), fold plan, panel and EC universe unchanged, code =
  10-29 code.
- Every one of the 2,372 gapfilled models grows; 10-29 had 2,362.

**Control: c12_compare and row attribution.**
- `c12_compare.py` old → new: `data/rescore_compare_{atp,gf}.json` and `SR/logs/rescore_compare_*.log`.
  No gated rows before or after.
- The plan's identity claim was made strict with `tools/r4_pheno_attrib.py`:
  - It rebuilds phenotype_fold's own `build_jobs` (sample, every arm's input EC set, donors) from the pre-repair
    inputs and from the new inputs.
  - The pre-repair inputs are the pre_repair cells and GFS, the 10-14 predictions, and an old picrust2 root with 28
    moved-aside runs. The rebuilt samples equal the old and new `meta.binding.genomes` exactly.
  - It then compares rows on (genome, arm, pipeline), ignoring wall_seconds.

| pass | genomes in both samples | rows compared | identical | changed (all attributed) | unattributed |
|---|---|---|---|---|---|
| ATP | 442 of 1,000 | 19,072 | 18,966 | 106: truth 80 (2 genomes, 813.301 and 1503470.4), PICRUSt2 prediction 20 (5 genomes), prFBA prediction 4 (1), order donor 2 (1) | **0** |
| gapfilled | 62 of 200 | 2,684 | **2,684** | 0 | **0** |

Every genome whose inputs did not change has identical rows (`genomes_unchanged_inputs_all_rows_identical`).
Reports: `control/r4/pheno_attrib_{atp,gf}.json`.

**Why the samples changed** (`tools/r4_pheno_sample_diff.py` → `control/r4/pheno_sample_diff.json`).
- The candidate set is the same 7,855 genomes, with the same amplicons.
- On the repaired data, 200 candidates have a different Δ(prFBA, PICRUSt2) EC F1. That moves the ΔF1 quintile
  boundaries.
- 13 of the 15 domain × template × quintile strata then change membership. The per-stratum shuffle is seeded by the
  stratum key, so a changed member list draws a different set, and so does the top-up.
- Result: ATP keeps 442 of 1,000 genomes, gapfilled 62 of 200. The plan foresaw this ("the sample can change").

**Old → new (controlled pipeline).** `control/r4/pheno_headlines.txt` and
`control/r4/pheno_manual_numbers_atp.json`. The manual recipe reproduces every old printed value: 755, 747,
−0.028, 22/54, ρ −0.14 with p 1.7e-4, and 2 McNemar media.

| quantity (tex) | old | new |
|---|---|---|
| ATP-informative paired genomes (L79, 632, 1193, 1222) | 755 (747 V4) | **759 (750 V4)** |
| ΔATP Hamming prFBA − PICRUSt2 (L79) | 0.13 (p 0.49, q 0.49) | **0.069 (p 0.79, q 0.79)** |
| all-genome sensitivity | 0.221 (p 0.11) | 0.284 (p 0.057) |
| EC F1 prFBA − PICRUSt2 on those genomes (L79, 1193) | −0.028 | **−0.031** |
| Spearman (EC-F1 Δ, phenotype Δ) (L634, 1190) | −0.14 (p 2e-4) | **−0.19 (p 1.5e-7)** |
| media with ATP signal (L79, 638, 1222) | 22 of 54 | 22 of 54 |
| McNemar media favouring PICRUSt2 at q ≤ 0.05 (L79, 631, 639, 1185) | 2 (Ac.NO, Etho.O2; q 0.002) | **4 (Ac.NO q 0.029, Ac.O2 0.040, Etho.O2 0.046, For.O2 0.046)** |
| **prFBA vs nn1 ΔATP Hamming (L635-637, 1186)** | **−0.146 (p 0.017, q 0.028)** | **−0.030 (p 0.41, q 0.51): no longer significant** |
| prFBA vs nn1, all genomes | −0.06 (p 0.28) | +0.06 (p 0.46) |
| draft rxn Jaccard median PICRUSt2 / prFBA / nn1 (L641-646) | 0.737 / 0.707 / 0.705 | 0.735 / 0.704 / 0.700 |
| essentiality disagreement PICRUSt2 / prFBA / nn1 (gapfilled) | 0.065 / 0.068 / 0.069 | **0.059 / 0.063 / 0.063** |
| growth log-ratio prFBA / nn1 / PICRUSt2 | 0.008 / 0.008 / −0.006 | **0.020 / 0.011 / −0.003** |
| gapfill burden median truth / prFBA / PICRUSt2 / nn1 | 82.5 / 80 / 84 / 82.5 | 79.5 / 78 / 81 / 81.5 |

Commit **5fac436**. It holds:
- both pass summaries;
- the default `data/phenotype_summary.json` (step 4 below);
- `data/rescore_compare_{atp,gf}.json`, force-added because `data/*` is ignored;
- identity files byte-unchanged.

---------------------------------------------------------------------------------------------------------------

## 3. Step 4 and R4.2

**The default summary.** `data/phenotype_summary.json`, the file score_folds reads beside `--cells` when no summary
is named, is now a **byte copy** of the new `phenotype_gf_c13/phenotype_summary.json` (sha256 8682122e…,
`control/r4/default_phenotype_summary.sha256`). Two reasons:
- that is the file phenotype_fold writes by default (`dirname(--outdir)`) for the gapfilled pass;
- it is the summary SPEC §15 binds the genus rescoring to, and C11/C12 need the gapfilled pass.

The deck reads the per-pass files, not this one. The replaced file is the 10-16 pre-registered summary (957f284b),
now in pre_repair.

**Check that the default binds when read by the default path**
(`control/r4/default_summary_binding_check.txt`). score_folds' own `phenotype_binding` and `phenotype_controls`
were called with the new genus inventory:
- read as the auto-discovered summary beside `--cells`, with `--require-phenotype-binding`, it is bound with no
  defects, and its `per_genome.jsonl` check resolves to `phenotype_gf_c13/per_genome.jsonl`;
- C11 and C12 pass, the same result as the explicitly named file.

**R4.2, pass 2.** Only the `score_folds` line of `c13_finish.sh` was run (`tools/r42_genus_pass2.cmd.sh`):
- `--require-phenotype-binding --phenotype-summary data/phenotype_gf_c13/phenotype_summary.json`;
- pins as in c13_finish and Run 3 pass 1: ModelSEEDpy-d884a50 and templates main2294299;
- no `--force`.

Unit `rescore-r42-genus`: 17:04-18:29, 5,125 s, cgroup peak 45.7 GiB (anon about 31 GiB; MemoryMax 50G); log
`SR/logs/score_genus_c13.log`. The 10-29 log is in pre_repair.

Checks (`control/r4/pass2_vs_pass1.txt`; the pass-1 copies are in `control/r4/pass1/`):
- exit 0. The phenotype arm reads "bound=True", with no binding defects, against cells d57709e3 and functions
  9463ed94.
- **C11 pass** (atp and gapfill identical, Jaccard 1.0).
- **C12 pass**: 0 of 8,660 models gated, 0 truth models gated, limit 0.10.
- Every control C1-C15 reads pass. C3-C10 and C13-C15 are identical to pass 1's apart from their provenance
  timestamp and commit; C1 and C2 are the other modes' entries, unchanged.
- The cells table was **reused by fingerprint**: inventory a32df47c, answer key c2b3c809cf33. `cells.csv.gz`
  d57709e3 and `.fp` 224f90fb are byte-unchanged.
- **Arms identical to pass 1** apart from `run_valid`/`run_invalid_reasons`; contrasts, nested_cv and
  annotation_ceiling identical.
- `meta.run_valid` False → True. The headline sentence no longer starts with INVALID, and it equals pass 1's
  provisional sentence: prfba_nestedcv 0.6867 vs PICRUSt2 0.7787, Δ −0.092 [−0.1073, −0.0812], skill 0.282.

**What else needed pass 2.** The plan lists only genus. The family scorer has no phenotype summary beside its cells
table, and Run 3 already passed its controls (valid). The C2 and C1 modes do not evaluate C11/C12.

**W8 attribution where it applies.** The cells table is unchanged since R3, where W8 ran: 0 unattributed. The W8
analogues for the new row-level outputs are:
- wobj `rowcheck` (§1);
- phenotype `attrib` (§2).

Both have 0 unattributed.

**collect.py over the 18 local R4 items** (R4P, R4S2, and the wobj/abstract part of R4H;
`control/r4/collect_r4_local.{json,md}`; the repo's `number_diff.*` is untouched):
- 77 values; every old value equals the printed number;
- 0 errors, 0 pending;
- 38 change at print precision.

Every one of the 38 is a phenotype or wobj value listed in §1-§2:
- PH-01/abs: 755 → 759; 0.13 → 0.07; p 0.49 → 0.79;
- PH-02: Ac.NO q 0.002 → 0.029, Etho.O2 q 0.002 → 0.046, and two more media now significant;
- PH-03/refusal_credit: ρ −0.14 → −0.19;
- PH-04/atp_vs_nn1: 0.146 → 0.030, p 0.017 → 0.407, q 0.028 → 0.509, p(all) 0.28 → 0.46;
- PH-06: Jaccard 0.737/0.707/0.705 → 0.735/0.704/0.700; essentiality 0.065/0.068/0.069 → 0.059/0.063/0.063;
  growth 0.008/0.008/−0.006 → 0.020/0.011/−0.003;
- wobj: 0.0120 → 0.0121, 27.9 → 27.8.

The manual values that no summary key stores were recomputed with one recipe that reproduces every old print
(`tools/r4_pheno_manual_numbers.py` → `control/r4/pheno_manual_numbers_atp.json`):
- 747 → 750 V4;
- ΔEC-F1 −0.028 → −0.031;
- 22 of 54 media (unchanged);
- p 2e-4 → 1.5e-7;
- McNemar 2 → 4 media.

Not recomputed: "EC F1 0.782 vs 0.781 (p = 0.64)" (PH-04, disc.atp_vs_nn1) needs nn1's per-genome EC F1 from the
cells table. It is left for Run 5.

---------------------------------------------------------------------------------------------------------------

## 4. Timings and memory (cgroup `memory.peak` read inside the unit; `stamps/units.jsonl`)

| unit | stage | wall | cgroup peak | MemoryMax | rc |
|---|---|---|---|---|---|
| rescore-r4-wobj-pop | R4.3 | 83 s | 4.0 GiB | 12G | 0 |
| rescore-r4-wobj-ctl | R4.3 control | 296 s | 7.8 GiB | 14G | 0 |
| rescore-r4-wobj-run | R4.3 | 7,094 s | 12.3 GiB | 20G | 0 |
| rescore-r4-pheno-atp | R4.1 (3 workers) | 6,563 s (stopped) | 28.9 GiB | 30G | stopped |
| rescore-r4-pheno-atp2 | R4.1 (8 workers) | 8,248 s | 49.8 GiB | 55G | 0 |
| rescore-r4-pheno-gf | R4.1 (12 workers) | 15,018 s | 88.5 GiB | 135G | 0 |
| rescore-r4-atp-post | R4.1 checks | 342 s | 26.7 GiB | 30G | 0 |
| rescore-r4-gf-post | R4.1 checks | 579 s | 16.7 GiB | 30G | 0 |
| rescore-r42-genus | R4.2 | 5,125 s | 45.7 GiB | 50G | 0 |

- Measured per-worker footprints: about 6.1 GiB per gapfilled worker at the 88.5 GiB peak (15 GiB parent) and
  about 4.3 GiB per ATP worker at the 49.8 GiB peak. The rule's 10 and 5 GB held, with margin.
- Lowest MemAvailable seen: 88 GB (15:09).
- No memory-max events in any unit. No kill or freeze anywhere: salterns' guard never fired.

**Salterns hold, with the arithmetic.**
- Idle MemAvailable is about 219 GB. Salterns' prep gate is 150 GB reserve + 40 GB floor = 190 GB. So its next
  merge or tadpole waits whenever our units hold more than about 29 GB.
- Held 12:09:31 → 16:34:04, the two phenotype passes: 4 h 25 m.
- Held 17:10:08 → 18:30:18, pass 2: 1 h 20 m.
- Total about 5.75 h. Its driver log shows every prep task skipped for "memory", then resuming normally.
- The alternative was to wait for salterns' prep phase, about 3.5 h of tasks, then its SPAdes phase (2 × up to
  160 GB). That would have pushed R4.1 by about a day.

---------------------------------------------------------------------------------------------------------------

## 5. What Run 5 needs

1. **Manuscript rewrites, not refreshes (R4P), at L79, 631-639, 1185-1186, 1190-1193 and 1222.** On the
   re-drawn sample:
   - the nn1 sentence (L635-637, 1186: "against nn1 the phenotype favours prFBA … p = 0.017, q = 0.028") is no
     longer supported: −0.030 media, p = 0.41, q = 0.51, all-genome p = 0.46;
   - the "two of 22 media" sentence becomes **four** (Ac.NO, Ac.O2, Etho.O2, For.O2; q 0.029-0.046; all favour
     PICRUSt2);
   - 755/747 → 759/750; 0.13 (p 0.49) → 0.07 (p 0.79); ρ −0.14 (p 2e-4) → −0.19 (p 1.5e-7); ΔEC-F1 −0.028 →
     −0.031;
   - PH-06 medians as in §2.

   Disclose why the sample changed: the stratified draw reacts to ΔF1 changes of 200 candidates
   (`control/r4/pheno_sample_diff.json`). 442 of 1,000 and 62 of 200 genomes are kept, and their rows are
   identical unless their own inputs changed.
2. **Still to compute:** "EC F1 0.782 vs 0.781 (p = 0.64)" (PH-04, disc.atp_vs_nn1), from per_genome.jsonl joined
   to nn1 cell F1. The other manual phenotype numbers are in `control/r4/pheno_manual_numbers_atp.json`, which
   reproduces every old print; it could become a W6 recipe.
3. **SPEC §15 / C12_SCOPE text** that quotes A-C13 counts:
   - 43,140 → 43,174 rows;
   - 8,614 → 8,660 rows;
   - 2,362 → 2,372 growing models.

   The 10-29 A-C13 record (`c13_compare_*`, `c13_vs_prereg_*`) stays as history. The new
   `rescore_compare_{atp,gf}.json` are the old → new records.
4. **Docs.** Re-derive `WEIGHTED_OBJECTIVE_RESULTS.md` from the new `summary.json` / `recomputed_contrasts.json`:
   only 4th-decimal moves, and "beat shuffled" 0.0121.
5. **Deck.** It reads `phenotype_{atp,gf}_c13/phenotype_summary.json` (A-C13 tag) and
   `weighted_objective/summary.json`; both are regenerated. `data/phenotype_summary.json` is now the gapfilled
   pass's summary (it was the 10-16 pre-registered one).
6. **TF pull hazard.** `tf/lists/pull_R4.3_posthoc.txt` lists `data/weighted_objective/` and
   `logs/wobj_run.log`. TF holds the 11-15 pre-repair copies, and installing that pull would overwrite FM's new
   files. The marker says so.
   - Checked at the end of this run: `population.json` 5d13cf1e…, `per_amplicon.jsonl` d44a27e7…,
     `summary.json` cca76fcc… are intact, and `logs/wobj_run.log` is this run's.
   - The TF Stage O commit 06251c8 touches no `weighted_objective` or `logs` path.
   - No `from_tf/pull_R4.3_posthoc*` directory exists yet.

   Exclude those paths when the R4.3 pull is installed.
7. **Disk.** `pre_repair/region_validation/shared_reference/data/phenotype_{atp,gf}_c13` (119 MB) and
   `control/r4/atp_restart/` (20 MB) can go after review. The temporary symlink roots
   `control/r4/r4_old_p2_*` are small.

Tools written for this run, in `$W/tools`:
- `r4_wobj_pop.cmd.sh`, `r4_wobj_ctl.cmd.sh`, `r4_wobj_run.cmd.sh`, `r4_wobj_popcheck.py`, `r4_wobj_rowcheck.py`;
- `r4_pheno_{atp,gf}.cmd.tmpl`, `r4_pheno_trim_tail.py`, `r4_pheno_attrib.py`, `r4_pheno_sample_diff.py`,
  `r4_pheno_manual_numbers.py`;
- `r4_atp_post.cmd.sh`, `r4_gf_post.cmd.sh`, `r42_genus_pass2.cmd.sh`.

No prFBA code changed.

Untouched:
- the ModelSEEDpy, ModelSEEDTemplates and ModelSEEDDatabase main checkouts;
- `presentations/` (including the user-edits deck);
- codiffusion, EmilyKin, KBase, the paper repo;
- the salterns rollout (only read).
