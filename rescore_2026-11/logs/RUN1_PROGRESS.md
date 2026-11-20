# Run 1 progress log (durable notes; the final report is RUN1_REPORT.md)

- 16:2x memory-peak method: systemd 255's end-of-unit "Memory peak" reads an emptied cgroup (256K for a 1.5 GB
  allocation). tools/U.sh + tools/unit_run.sh read memory.peak from inside the unit; self-test 1.48 GB OK.
  Every unit appends a line to stamps/units.jsonl (cgroup peak incl. page cache; GNU time max RSS of the largest process).
- Incident: `close_exclude_rep_gap.py --help` is not argparse; it started a gap-close on region_validation_kmeans (d20)
  and was killed during "loading selection + hits" (before any cache access; cache stamp unchanged, no .tmp, no tree writes).
- R1.0: stamps/R1.0.json (cache 5118271515 B, 2026-11-13T10:38:50, sha256 05f03bcd..., keys 4ef0950c..., 159,753 entries;
  catalogue current; CORRECTIONS 429c7b7a...; d20 baseline 32/32 OK; prFBA 23773d3, paper 0214a14); key list and
  1,006 repaired values saved; stamps/ignored_before.sha256 (442 files). Step 5 (uploader commit) skipped: D5.
- R1.1: SNAP f0a84bf1 (= logged pre-repair), SNAP20 86cfdbe1, SNAP20J e57a129c re-verified; copies cl_snap_copy.json,
  d20ctl/snap20_copy.json byte-identical. tools/test_tools.py 9/9 pass.
- R1.2 C-U PASS: score_taxacc --force + mcnemar --force in region_validation_d100: git diff clean on all 9 tracked outputs,
  taxacc_per_record.jsonl sha256 091cc592... == R1.0 stamp. 263 s, cgroup peak 22.3 GB.
- R1.2 C-L1 PASS (gate for R1.3): select_references legacy on a SNAP copy reproduces region_validation_d100_legacy
  selection.json record for record (67,991/67,991 identical; _meta identical with generated_from masked; 0 lazy misses:
  copy stat unchanged; live cache stat unchanged). touched records 20,644 (= plan exposure table). control/cl1_*.json.
- reps old: region_validation_kmeans_d100 dump on SNAP reproduces committed cells/summary/long byte-identically (456 s).
- R1.2 C-L2 PASS: kmeans_d100_legacy 74,704/74,704 identical, touched 13,265; C-L3 PASS: random_d100_legacy 75,328/75,328
  identical, touched 12,994. 0 lazy misses (cl_snap_copy stat unchanged throughout). select_references did not drift
  after 10-20 -> the legacy re-score is data-only.
- reps old: random_d100 and d100_legacy (--legacy) dumps reproduce committed cells/summary/long byte-identically.
- R1.2 C-L4 on SNAP20: 17 records change selection, all `no_gene_set -> selected` on the 4 SNAP20J genomes
  (2072416.12 x6, 294671.3 x5, 33953.3 x5, 176233.20 x1; + one tier_cap knock-on 2874285.4) -> per handoff item 3
  switched to SNAP20J (kept as control/d20_code_delta*_snap20.json). C-L4 on SNAP20J: 0 substantive differences,
  n_selected identical in 67,991/67,991, _meta summary identical; 67,427 records differ only by schema fields added
  since June (selected[].equivalent_genomes, selected[].gene_set_empty, provenance.n_equivalent_genomes); knobs gain
  select_all_* keys and select_mode. 0 lazy misses. => d20 selection code delta = output schema only.
  control/d20_code_delta.json (+ _substantive.json). C-D20 run anyway on SNAP20J (gap-close postdates the July genecap).
- reps old: kmeans_d100_legacy reproduces committed cells byte-identically.
- C-D20 control (SNAP20J copy, today's code, control tree $W/control/region_validation_ctl, all steps via ctl_run):
  gap close 0 missing; snapshot copy unchanged (0 lazy misses); 699 s, cgroup peak 42.6 GB. Committed -> control:
  taxacc tables + taxacc_per_record.jsonl byte-identical (summary differs only in _meta paths); mcnemar pairs CSV,
  self_recovery_diagnosis/failures identical (mcnemar_regions.json: paths + correct_definition text); gene capture
  17/135,982 cells differ (exclude_0.987 only; reps absent from the July cache, filled by the October gap close).
  Summaries written to region_validation/data/control_currentcode_precache/ (to commit).
- reps old: d20 original on SNAP20 + June selection does NOT reproduce the July cells (35 cells; truth_empty 958 vs 953):
  the July cache state != SNAP20. d20 kmeans reproduces on SNAP byte-identically.
- R1.3 launched 17:38 (unit rv-rescore-post, MemoryMax=150G, MemAvailable 209-220 GiB at launch; KBase unit no longer
  running at that time; only codif units beside). Move-aside: 66 items -> pre_repair (logs/r13_moveaside.log).
- R1.3 original union done 17:54:57 (16.7 min): P1/P3/P4 to_fetch=0; check 5 all PASS; gene_capture_{cells,by_region,long}
  BYTE-IDENTICAL to expected/repro_0987 (11-13 prediction). include ok 67010->67700, truth_empty 953->278, pred_empty 28->13;
  exclude_0.987 ok 59571->61047, abstained 6407->6423, pred_empty 1060->243.
- R1.3 original legacy done 18:26:46 (31.7 min): 0 fetches (P1/P3/P4 to_fetch=0; gap close 0 missing); check 5 PASS.
  Selection old->new (control/r13/sel_region_validation_d100_legacy.json): 1,506 records changed, ALL within the 20,644
  touched (0 untouched); n_selected changed in 674; mean n_selected 1.307129 -> 1.317630 (non-zero 1.317769 -> 1.320875);
  abstained 549 -> 167; no_gene_set exclusions 5,098 -> 1,053; stage moves no_gene_set->marginal_gain 1,859,
  ->included 1,148, ->tier_cap 1,023. Genecap include ok 67,706 / abstained 7 / truth_empty 278; exclude ok 61,137 / abst 6,576.
- R1.3 kmeans union done 18:43:52 (17.1 min), 0 fetches.
- R1.3 kmeans legacy done 19:16:10 (32.2 min), 0 fetches. Selection: 1,103 changed, all within 13,265 touched; n_selected
  changed 170; mean 1.426122 -> 1.428424; abstained 211 -> 129; no_gene_set 4,718 -> 1,728.
- R1.3 random union done 19:33:00 (16.8 min), 0 fetches. POST cgroup peak so far 101.6 GiB (page cache incl.).
- R1.3 random legacy done 20:05:41 (32.6 min). POST unit total: rc 0, 8,844 s (2 h 27 min), cgroup peak 111.97 GiB
  (page cache incl.; grows tree by tree), max RSS of the largest process 62.5 GiB. 0 fetches in every P1/P3/P4/P5;
  post-R1.3 cache sha256 == R1.0 (05f03bcd...). Per-stage memory: control/r13/post_stage_mem.json.
- R1.3 checks (logs/r13_checks.log): check 5 PASS x6; original union cells/summary/long byte-identical to the 11-13
  prediction; union trees: selection.json + taxacc_per_record.jsonl sha256 == R1.0, git shows only the 2 genecap files
  changed; legacy selections: changed records all touched (1,506/20,644; 1,103/13,265; 1,822/12,994).
  Headlines in control/r13/genecap_headline_d100.json.
- R1.4 (d20, flock held non-blocking, unit rescore-r14-d20): 1,721 s, cgroup 48.3 GiB, max RSS 43.3 GiB. Added a legacy
  candidate prefetch of the d20 hits before select_references (as the d100 driver does): 28,639 cached, to_fetch=0.
  select 0 lazy misses; gap close 0 missing; cache stat unchanged throughout. d20 union trees: selection.json and
  taxacc_per_record.jsonl sha256 == R1.0; only the 2 genecap files changed.
  d20 selection control -> new (--d20): 1,417 changed, 1,411 touched by repaired/poisoned genomes, the other 6 all on
  the SNAP20J reset genomes (33953.3 x5, 176233.20 x1: no_gene_set -> selected) = pre-11-11 cache history (these 4
  are neither repaired nor poisoned; they hold PGFams in SNAP). mean n_selected 1.298275 -> 1.309217.
  d20 original genecap committed -> new: include ok 66,978->67,706, abst 55->7, truth_empty 958->278; exclude_0.987
  population 52.968 -> 53.707, ok 50,495->51,733, abst 16,538->15,980. self-recovery failures now 10,821 (15.92%).
- R1.5 E1 (unit e1-rescore): smoke (--limit 300) 600/600 rows match; full 615 s (unit 703 s), max RSS 9.3 GiB; 0 provider
  misses, supplement empty (not written) -> no --merge-supplement, no cache write. Built-in 0.987 reproduction
  byte-identical vs the NEW P7 files (3/3). gene_capture_tau_compare.json == genecap_repaired_impact copy in content
  (only _meta.repro_0987_byte_identical differs: False there [pre-P7], True now). tau0995 cells/summary/long
  byte-identical to the 11-13 prediction (expected/genecap_repaired/*_tau0995*).
- R1.6 (20:47): cache UNCHANGED vs R1.0 (size/mtime/sha256/keys); 1,006 repaired values unchanged; catalogue current;
  no rebuild. FREEZE 2026-11-15T20:47:47. stamps/R1.6.json.
- d20 W8: reps dumps 4/4 reproduce; k-means PASS; original control->new 14 + random 1 non-repair cells, all explained
  (SNAP20J genomes / rep cached 10-21): control/w8/d20_unattributed_explained.json.
- Commits: 81e9b34 (R1.3), 000fcaa (R1.5), 3ffc97d (R1.4 + C-D20 summaries). Not pushed. Report logs/RUN1_REPORT.md.
