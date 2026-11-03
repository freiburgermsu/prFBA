# prFBA jobs: final status (audit of 2026-11-02)

Checked 2026-11-02 22:00 CST to 2026-11-03 09:01 CST. Read-only throughout: nothing was
stopped, started, fetched, staged or committed.

**Sources:**
- `systemctl --user list-units --all`;
- the user journal (`journalctl --user -u <unit>`) and the kernel journal
  (`journalctl -k`);
- each job's log, for its `EXIT=` line. Most units wrap their job as
  `...; echo EXIT=$? >> log`, so the unit itself always "succeeds" and the log's `EXIT=` is
  the job's real exit status;
- the outputs themselves.

**Times** are local: CDT until 11-01, CST after.

**Every prFBA unit was launched with `--collect`.** None is still loaded except
`otune-run`, which is running. The rest are known only from the journal.

## 1. Summary: this project's recent jobs

| Job | Unit(s) | Start → end | Exit | Outputs | Left undone |
|---|---|---|---|---|---|
| Template-edit validation, smoke | `tmpl7-smoke` ×2 | 10-29 18:42:50 → 18:44:03 (stopped: 15G cap too small); 18:46:24 → 19:03:10 | 2nd run `EXIT=0` | 72/72 rows ok (`validation/smoke_2genomes/`) | — |
| Template-edit validation, gate | `tmpl7-gate` | 10-29 20:41:59 → launched `tmpl7-chain` at 22:12:00 | process gone; no deactivation line in the journal | — | its own MemAvailable check printed empty (quoting); `launch_tmpl7.sh` re-checked: 216 GB |
| Template-edit validation, stage a + stage b part 1 | `tmpl7-chain` (MemoryMax 77G, swap 0) | 10-29 22:12:00 → 22:54:18 | **`oom-kill`** (memcg, at the unit's 77 GiB cap) | stage a complete (2,362 rows, exit 0); stage b 20/200 genomes (854 rows) | resumed, below |
| Template-edit validation, resume attempt 1 | `tmpl7-resume` (120G) | 10-30 10:33:48 → 10:36:20 | stopped by hand (log misread) | 0 genomes; rows unchanged (854) | — |
| Template-edit validation, stage b part 2 | `tmpl7-resume2` (120G, swap 0) | 10-30 10:39:10 → 14:34:21 | **`EXIT=0`**; stage b 0; post-steps 0 | 200/200 genomes, 8,614 rows, all `ok`; both comparisons written; peak 72.2 GiB | outputs untracked; branch not pushed; PR not opened (section 5) |
| Weighted-objective test | `wobj-run` (MemoryMax 20G) | 10-29 20:40:42 → 22:25:16 | **`EXIT=0`** (`logs/wobj_run.log`) | 360 rows, all `ok`, 360 distinct keys (n1 60, n2 100, n3 100, n4+ 100); `summary.json` written 22:25:13 with verdict; analysis figures and `recomputed_contrasts.json` (10-30); committed in prFBA `c91cf19` | the single-model fix it led to, ModelSEEDpy-prfba `b6857a0` (+ `09856f5`), is on no remote |
| Depth-100 deployments (A codiffusion, B EmilyKin) | `dep100-chain` | 10-29 20:41:02 → 22:33:15 | **`EXIT=0`**; A exit 0 (21:33:05); B exit 0 (22:33:15) | complete, below | depth-100 build not uploaded to KBase (RESULTS.md) |
| A-C13 ATP pass | `c13-atp` | 10-29 10:15:06 → 12:27:01 | `EXIT=0` | `phenotype_atp_c13`: 43,140 rows ok, 1,000 genomes, 0 duplicates, identity pass | — |
| A-C13 gapfilled pass | `c13-gf` | 10-29 10:15:06 → 13:09:53 | `EXIT=0` | `phenotype_gf_c13`: 8,614 rows ok, 200 genomes, 0 gated, identity pass | — |
| A-C13 finish (compare + genus rescoring) | `c13-finish` | 10-29 10:17:51 → 14:34:19 | `EXIT=0`, `FINISH_DONE` | `summary_genus.json`, `controls.json` (14:34); controls pass (`score_genus_c13.log`); A-C12 copies kept as `*.a-c12` | — |
| Reversal-fix validation | `revcheck-fixed`, `revcheck-unfixed`, `harness-on-fix` | 10-29 10:21:51 → 10:35:35; 10:37:12 → 10:49:06 | all `EXIT=0` | scratchpad `revcheck/{fixed,unfixed}.json`, `harness_on_fix/out` | — |
| A-C12 ATP pass (pinned) | `c12pin-atp` | 10-27 18:29:10 → 21:12:58 | `EXIT=0` | `phenotype_atp_c12fix`: 43,140 rows ok, ModelSEEDpy d884a50, identity pass | — |
| A-C12 gapfilled pass (pinned) | `c12pin-gf` | 10-27 18:42:31 → 22:29:33 | `EXIT=0` | `phenotype_gf_c12fix`: 8,614 rows ok, 0 gated, identity pass | — |
| A-C12 finish | `c12fix-finish` | 10-27 22:29:40 → 23:46:13 | `EXIT=0`, `FINISH_DONE` | compare files and genus rescoring (`*.a-c12`) | — |
| A-C12 first attempts (unpinned) | `c12fix-atp`, `c12fix-gf`; `c12fix-atp2/gf2`; `c12fix-atp3/gf3` | 10-27 17:02 → 18:23 | `oom-kill` ×2 (gf: 125.0G peak, 2.0G swap); stopped ×4 | partial rows kept as `phenotype_{atp,gf}_c12fix.msp-c2e31a0-partial` (9,104 rows / 212 genomes; 1,954 / 46) | superseded by the pinned runs |
| A-C12 checks | `c12-smoke`, `c12chk-atp/gf`, `pin-repro-atp/gf`, `pin-smoke-gf`, `draftdiff-pin/moved` | 10-27 16:44 → 18:49 | all exit 0 | scratchpad `c12_smoke`, `worker_history_check`, `pin_control`, `draft_diff` | — |
| EGC draft scans | `c12-scope-scan`, `c12-draft-loops`, `c12-draft-loops-fix` (unpinned); `egc-scope-pin`, `egc-scope-var`, `egc-loops-var` (pinned) | 10-27 17:05 → 19:19 | pinned ones `EXIT=0`; unpinned ones superseded | `data/egc_analysis/scope_scan_{200,fixonly,variants}.jsonl`, `draft_loops_*`; unpinned copies kept as `*.msp-c2e31a0.jsonl` | — |
| GPU end-to-end tests | `gpucand200-e2e`, `traceback-e2e` | 10-27 10:32:24 → 11:46:20; 15:34:58 → 16:32:29 | `exit=0` both | scratchpad `gpucand200_e2e` (csv + method_selection; its hits JSON is no longer there), `traceback_e2e` (holds a 3.56 GB hits JSON) | optional cleanup of the 3.56 GB scratchpad file |
| Self-recovery checks (three panels) | `msv-selfrec-region-validation-{d100,kmeans-d100,random-d100}` | 10-26 13:02:47 → 13:03:35 | unit ok | scratchpad `ms_verify/twopanel/selfrec_*.json`, `n_missing_record` 0 in each | — |
| Three-panel bootstraps | none (ran outside systemd) | — | — | `three_panel/*.json` regenerated and committed in prFBA `5d2a204` (10-29 20:43); working tree clean | the paper's new intervals are in prFBA_paper `31e37a8`, unpushed |
| Depth-100 rerun of the amplicon arm | `rv-d100-align`, `rv-d100-post`, `rv-d100-picrust2`, `rv-d100-downstream`, `rv-d100-checkall` | 10-20 10:56 → 10-21 14:39 | final runs ok: "all three panels aligned", "post-all finished", "[p2] done", "[down] done", "[check-all] exit 0" | all six trees hold `post.done` and `close_exclude_rep_gap.done`; committed in prFBA `d4d2f76` | — (first `post` run `oom-kill`, first `downstream` exit 1 and first `align` stopped, each rerun) |
| Function-probability estimators | `pfunc-est` ×3 (same command) | 11-02 22:38:55 → 23:06:06 | last run `EXIT=0` (the log is overwritten each time) | `data/pfunc_estimators/` (scores, selected, weights E1/E3, assertions) | — |
| Objective tuning, stub + smoke | `otune-stub`, `otune-smoke` | 11-02 22:25:26 → 22:25:43; 22:42:35 → 22:47:50 | `EXIT=0` both | `data/objective_tuning/smoke/` | — |
| **Objective tuning, full run** | **`otune-run`** (MemoryMax 34 GiB) | **11-03 08:50:22 → running** | — | 45 of 360 jobs ok at 09:01; 0 errors in the log | **in progress; not touched** |

## 2. Details

### 2.1 Template-edit validation (`ModelSEEDTemplates-draftloops`)

**The OOM kill.** The kernel journal, 10-29 22:54:15, records a memory-cgroup OOM in
`tmpl7-chain.service`.
- Usage was 80,740,352 kB against a limit of 80,740,352 kB (77.0 GiB).
- Anonymous memory was 76.8 GiB.
- The kernel killed the stage b parent (anon-rss 14.2 GiB).
- Its ten workers held 5.7–6.8 GiB each, against the README's estimate of 4.4 GB.

**The stopped attempt.** `tmpl7-resume` was stopped because its log was misread. Its own
log line said "20 already done, 180 to build"; the "0 already done" line belonged to the
first run. The resume script's `rc_a=skipped` was changed to `rc_a=0` before `tmpl7-resume2`.

**Outputs:**
- `validation/scan/` (stage a);
- `validation/gf_pass/`: `per_genome.jsonl` with 8,614 rows, `phenotype_summary.json`,
  `identity_control.json` (pass), `compare_vs_c13.json`, `c12_compare_vs_c13.json`, and the
  854-row backup `per_genome.oom-2254.jsonl`, byte-identical to 854 of the final rows;
- `validation/logs/`.

**Results:**
- 0 of 8,614 rows gated;
- all 2,362 gapfilled models grow;
- ATP vectors identical in every row;
- burden change mean +0.013 reactions.

The write-up is in `validation/PROVENANCE.md` (updated) and `validation/REPORT.md`.
`REPORT.md`'s text is in the auditing agent's final message, because the write tool
refused to create that file from a subagent. Save it from there.

**Left undone:**
- everything under `validation/` is untracked except the scripts committed in `1f84188`:
  `scan/`, `gf_pass/`, `logs/`, `smoke_2genomes/`, `chain_tmpl7_resume.sh`, and the new
  `PR_TEXT.md`, `report_rows.py` and `gf_pass/report_rows.json` (`REPORT.md` is still to
  be saved);
- the branch is not pushed to the fork, and the PR is not opened (draft text in
  `validation/PR_TEXT.md`);
- the proposed pointer paragraph for prFBA's `C12_SCOPE.md` is not applied.

### 2.2 Weighted-objective test

`scripts/launch_wobj.sh` ran two stages in one unit: the population stage, then 360
amplicons with 8 workers and `--maxtasksperchild 2`.

**The log ends:**
- "[wobj] wrote 360 rows -> …/summary.json";
- a verdict: H1 "detectable but below SESOI", and H2 "reduces degeneracy … NOT specific
  to p";
- peak RSS, parent 0.743 GB and largest worker 1.804 GB;
- `EXIT=0`.

**Outputs:** `per_amplicon.jsonl` holds 360 rows, all `ok`, with 360 distinct keys. The
summary's `n_ok_by_stratum` is n1 60, n2 100, n3 100, n4+ 100. The run used templates
`21a8084` and ModelSEEDpy `c2e31a0`. The data directory is tracked: 20 files, committed in
`c91cf19`.

**Left undone:** the single-model-path fix (ModelSEEDpy-prfba `b6857a0`, plus `09856f5`
added on 11-03) is not pushed;
its PR text is prFBA `region_validation/shared_reference/PR_ModelSEEDpy_prfba_single_model.md`.

### 2.3 Depth-100 deployments

`logs/dep100-chain.log`: A exit 0 at 21:33:05, B exit 0 at 22:33:15, unit
`memory.peak=26624 MiB`, `EXIT=0`.

**Codiffusion (A):** `compare/compare_codif.json` and the files, read-only.
- 3,276 ASVs, with 100 hits in every one of the 3,276 records;
- answered 3,225, abstaining 51, failed 0;
- `synthetic_genomes_d100`: 3,226 top-level JSON (3,225 genomes plus the manifest), with
  3,225 each in `provenance/`, `genome_objects/` and `taxonomies/`, and no zero-byte genome
  file;
- newest file 10-29 21:32, so nothing has changed since.

**EmilyKin (B):** `compare_emilykin.json`.
- 3,950 ASVs, with 100 hits in all 3,950 records;
- answered 3,755, abstaining 195, failed 0;
- 3,756 top-level JSON, with 3,755 each in the three subdirectories, and no zero-byte
  file;
- newest file 10-29 22:32.

**No OOM.** The kernel journal for 10-29 20:41–22:34 holds no OOM or killed-process line.
This settles what `RESULTS.md` section 1.4 left open, since it could not read the kernel
journal: the 26,624 MiB `memory.peak` was page cache, not a kill.

**Left undone:**
- the depth-100 build is local only. KBase workspace 261948 still holds the 20-hit
  codiffusion build, and the upload command is in `deployment_d100/README.md` section 8;
- `deployment_d100/smoke/` is untracked; the rest is tracked (31 files).

### 2.4 A-C12 and A-C13 phenotype passes

| Output (`region_validation/shared_reference/data/`) | Rows | Genomes | Status | Duplicates | ModelSEEDpy | Templates | Identity control |
|---|---|---|---|---|---|---|---|
| `phenotype_atp_c12fix` | 43,140 | 1,000 | all ok | 0 | d884a50 | 21a8084 | pass |
| `phenotype_gf_c12fix` | 8,614 | 200 | all ok | 0 | d884a50 | 21a8084 | pass |
| `phenotype_atp_c13` | 43,140 | 1,000 | all ok | 0 | c2e31a0 | 21a8084 | pass |
| `phenotype_gf_c13` | 8,614 | 200 | all ok | 0 | c2e31a0 | 21a8084 | pass |

Each summary's `n_models` equals its row count, and `n_gated_out` is empty. All six logs
(`pheno_c12fix_{atp,gf}`, `pheno_c13_{atp,gf}`, `c12fix_finish`, `c13_finish`) end in
`EXIT=0`. Nothing is left undone.

### 2.5 Three-panel bootstraps

No systemd unit ran them. Commit `5d2a204` (10-29 20:43) made them reproducible: the strata
are sorted before weighting. It regenerated `three_panel_comparison{,_d100,_d100_legacy}.json`,
`three_panel_mag_gtdb.json` and `three_panel_mag_panelarm.json`. Its message says they are
byte-identical under PYTHONHASHSEED 0, 777 and 12345, with interval bounds moving by up to
0.0074.

`git status` is clean for `three_panel/`. The manuscript's intervals were updated in
prFBA_paper `31e37a8`, which is unpushed (section 5).

### 2.6 Earlier October units (10-15 to 10-21): journal results only

These predate the A-C12/A-C13 amendments and the depth-100 standard, and later work
superseded them. Their outputs were not re-verified here.

**Failed or stopped, then rerun:**

| Failed or stopped | Result | Rerun |
|---|---|---|
| `prfba-pheno-gf`, `prfba-pheno-gf2` (10-15) | oom-kill | `prfba-pheno-gf3` (10-16 11:36 → 13:25), "pass 1 COMPLETED cleanly", 63,922 rows |
| `p2mag-loci` (10-16) | exit-code | `p2mag-loci2`, "[loci] done" |
| `prfba-rv-random`, `prfba-rv-random2` | exit-code | `prfba-rv-random3` ok |
| `prfba-tmplval-B`, `prfba-tmplval-BC_full` (10-21) | oom-kill | `prfba-tmplval-rerun` ok |
| `prfba-tmplval2-BC_full_rev` | oom-kill | `prfba-tmplval2-BC_full_rev2` ok |
| `rv-d100-post` (1st) | oom-kill | 2nd run ok (section 1) |
| `rv-d100-downstream` (1st) | exit 1 | 2nd run ok (section 1) |
| `rv-d100-align` (1st) | stopped | relaunched (section 1) |

**`p2mag-sample` (10-16, exit-code) is an artefact, not a failure.** Its only batch
finished ("batch_00 done 12:09:15", with a `DONE` marker). The wrapper then failed with
"line 33: f: command not found". `run_picrust2.sh` was modified at 12:02, while the unit
was running, and bash reads a script as it runs.

**Unit ok:** `prfba-campaign-family`, `prfba-pheno-atp`, `prfba-score-family`, `p2mag-prfba`,
`p2mag-panel`, `prfba-mag-rnd-{fetch,fetch-full,panel,score}`, `prfba-gtdb-arms`,
`prfba-pheno-fixed`, `prfba-score-genus-pheno`, `prfba-egc-{rebuild,analyze,cf-*,tmplfix-*}`,
`prfba-c2`, `prfba-cache-snap`, `prfba-calib-pfunc{,2}`, `calib-fix-full`, and the other
`prfba-tmplval{,2,3}-*` units.

**Trivial:** `run-u1683`, `run-u1735` and `run-u1736` were `ulimit` checks.

## 3. Not this project's (current state only; not touched)

**codiffusion session, `codif-*` units:**
- 72 distinct `codif-*` unit names appear in the journal since 10-24;
- none is loaded or active now, and no codif timer remains;
- the most recent, `codif-orient-uncap-dc-174`, ran 11-02 04:12:54 → 04:18:45 and ended
  with exit-code;
- the last batch (11-01 23:00 → 11-02 04:18) ended with a mix of unit-ok and stopped.

**EmilyKin session:** no systemd units; it uses tmux.
- **27 `ek_*` tmux sessions** exist (ek_pilot_atp3–9, ek_r11_*, ek_r12_*, ek_smoke*,
  ek_verify_*), every pane alive at a bash prompt.
- **On 11-02 at about 22:10:**
  - `ek_r12_pilot` was running `run_r12_pilot.sh`, with `simulate_comm.py --days 199`
    (started 21:17);
  - a separate `simulate_comm.py --days 560` (started 20:42) and `diag2.py` in
    `ek_r12_diag2` (started 21:58) were also running.
- **On 11-03 at 09:01** none of these is running. The only EmilyKin-related process is a
  python script from that session's own scratchpad (`…-EmilyKin/…/scratchpad/r12…`),
  12 minutes old.

**Wetland / GTDB-Tk, another project:**
- `gtdbtk-wetland.scope` ended 10-29 22:08:57 (146.2G memory peak);
- `gtdbtk-wetland-run2.scope` ran 10-30 17:11:41 → 10-31 00:33:52;
- `cgtest-wetland.scope` (10-29 18:30) was a cgroup test.

## 4. Repositories: unpushed commits and uncommitted state

Read-only git (`git status -sb`, `git log`); nothing was fetched. The remote-tracking refs
are as of each repo's last fetch.

| Repo (branch) | Unpushed | Working tree | Last fetch |
|---|---|---|---|
| `prFBA_paper` (`depth100-standard`) | ahead 2 of `origin/depth100-standard`: `62fbaf1` "Report the weighted-objective test and move the deployments to depth 100"; `31e37a8` "Update the three-panel intervals to the reproducible bootstrap" | 4 untracked backups (`k_md5.txt`, `prFBA.tex.bak-audit`, `prFBA_andrew_voice.tex.bak-audit`, `prfba_refs.bib.bak-recon`) | 10-29 18:19 |
| `ModelSEEDpy-prfba` (`fix-prfba-single-model`) | no upstream configured. On no remote-tracking ref: `09856f5` "MSProbability.prFBA: refuse a model that cannot grow instead of returning all-zero fluxes" (committed during this audit, on 11-03 before 09:08) and `b6857a0` "MSProbability.prFBA: make the single-model path run, and make its two policies explicit" | clean | 11-01 11:38 |
| `ModelSEEDpy-reversals` (`fix-gapfill-reversals`) | none: HEAD `e1254ba` = `origin/fix-gapfill-reversals` (freiburgermsu/ModelSEEDpy PR #19, per project memory); no upstream configured, but pushed | clean | 11-01 11:38 |
| `ModelSEEDTemplates-draftloops` (`fix-draft-atp-loops`) | no upstream; `1f84188` and `7fd1ab1` are on no remote-tracking ref. The parent `21a8084` is on `fork/fix-irreversible-reaction-directions` | untracked `validation/` outputs and new files (section 2.1) | 11-01 11:45 |
| `prFBA` (`hardware-aware-alignment-dispatch`), for reference | even with origin at HEAD `2abcd1a` (09:08 on 11-03; it was `b46598b` when the audit started, so another session committed and pushed meanwhile) | 2 modified files from the 10-16 MAG/PICRUSt2 work (`picrust2_comparison/README.md`, +260 lines; `scripts/run_picrust2.sh`, a `TRAITS` env var), and many untracked files | — |

## 5. Left undone, consolidated

1. **Template-edit validation:**
   - commit `validation/`. Suggested: summaries and logs only, keeping the 30.7 MB
     `per_genome.jsonl` out;
   - push `fix-draft-atp-loops` to the fork, and open the PR from `validation/PR_TEXT.md`.
     It stacks on the five-reaction PR;
   - optionally, add the `C12_SCOPE.md` pointer (proposed at the end of `REPORT.md`).
2. **ModelSEEDpy-prfba `b6857a0` and `09856f5`:** push and open the PR.
3. **prFBA_paper:** push `62fbaf1` and `31e37a8`.
4. **Depth-100 deployment:** upload to KBase if wanted (`deployment_d100/README.md`
   section 8).
5. **prFBA:** decide on the two 10-16 modifications (`picrust2_comparison/README.md`,
   `run_picrust2.sh`) and the untracked MAG-arm files.
6. **`otune-run`:** still running (45 of 360 at 09:01). Check `logs/otune_run.log` for
   `EXIT=` when it ends.
7. **Optional scratchpad cleanup:** `traceback_e2e/asv_top20_alignment_hits.json` is
   3.56 GB.

## 6. Addendum (2026-11-03, after the audit)

- **`otune-run`** (Stage O, 360 amplicons × 90 objectives) finished 10:27, `EXIT=0`: 360 / 360
  rows ok, 32,400 / 32,400 solves optimal, released check 360 / 360. **`otune-fva`** (the
  confirm-half degeneracy pass, 193 amplicons × 4 finalists) finished `EXIT=0`. Results:
  `OBJECTIVE_TUNING_RESULTS.md` (prFBA b52984f), independently re-derived.
- **Stage E** (`pfunc-est`) results committed with their independent re-derivation (f2ec19a).
- **Template-edit validation** (item 5.1): `validation/` committed on `fix-draft-atp-loops`
  as `627ffaa` (summaries, logs, the smoke, row files gzipped); the `C12_SCOPE.md` pointer is
  prFBA bcd8c49. Push and PR remain.
- **No unit of this project is running.**
- **Still unpushed:** `prFBA_paper` `depth100-standard` (31e37a8, 62fbaf1, and the manuscript
  update that follows these results); `ModelSEEDpy-prfba` `fix-prfba-single-model` (b6857a0,
  09856f5); `ModelSEEDTemplates-draftloops` `fix-draft-atp-loops` (7fd1ab1, 1f84188,
  627ffaa). `fix-gapfill-reversals` (e1254ba) is already on origin.
