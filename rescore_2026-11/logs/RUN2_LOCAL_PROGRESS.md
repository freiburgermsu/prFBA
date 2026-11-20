# Run 2 (local part) progress log

- 2026-11-15T22:08:56-06:00 start. stamp_state --compare R1.6: CACHE UNCHANGED (stamps/R2L.0.json). Writer pgrep: none. MemAvailable 210 GiB.
- 2026-11-15T22:14:57-06:00 code: 772832a (fetcher UA/rate limit/both selections, D10; README note; tests 18 pass). Launched unit rescore-r2-selctl (prfba_selection controls: d100 on SNAP, d20 on SNAP20J).
- 22:14-22:33 prfba_selection controls (control/r2/cmp_*.json):
  - d100 legacy on SNAP: RECORDS IDENTICAL (67,991).
  - d100 default on SNAP: 305 records differ (reps only); every one swaps a genome cached before the 10-20 run for one cached after it
    (insertion index >= 148,300 = the cache size the 10-20 log reports). Emulation SNAP[:148,300] (control/r2/snap_d100_1020.json): RECORDS IDENTICAL.
  - d20 on SNAP20J: default 3 / legacy 20 records differ (late-cached genomes + 33953.3, a June-reset key that had PGFams by 09-17).
    Emulation SNAP20[:111,789] (09-17 log size): default and legacy RECORDS IDENTICAL (meta n_reps_missing 0->23: 20 poisoned keys
    + 3 re-inserted later sit past the cut; absent == [] for the cache-only provider). SNAP20J[:111,789]: default identical, legacy 5 differ (33953.3 only).
- 22:24 R2.4 move-aside (logs/r24_moveaside.log): mv 8 outputs, cp products/source genomes/figures.
- 2026-11-15T22:27:45-06:00 launched rescore-r24-p2d (run_picrust2_d100.sh, FMAPS e8446e3e, SEL_WORKERS=12; then dialect).
- 22:29 rescore-r24-p2d OOM-killed (cgroup 90G; anon 95.7 GB in the 12-worker fork pool of prfba_selection; nothing written). 2026-11-15T22:32:57-06:00 relaunched as rescore-r24-p2d-b, SEL_WORKERS=8, MemoryMax=150G.
- 22:52 R2.4 done (unit rescore-r24-p2d-b: 1,155 s, cgroup peak 117.1 GiB, max RSS 48.7 GiB). Fetch: genomes to_fetch=0; pgfams 39,229 (default) + 2,183 (legacy) new ids. EC n 67,038 -> 67,713. 2026-11-15T22:52:31-06:00 launched rescore-r24-checks.
- 22:52-22:55 R2.4 checks (control/r2/d100/): source genomes unchanged (to_fetch 0, cmp equal); domain_routing.csv unchanged; selection meta missing 0/0;
  legacy old->new 2,503 changed, all touched (PASS); default 755 changed, 200 untouched = all in the 305 cache-history records
  (committed->SNAP control); SNAP control->new 473 changed, 0 untouched (PASS). GFS: 729 changed, all repaired, 0 unexplained;
  maps: +39 pgfam_ec, +40 pgfam_rxn, 0 changed, ec_rxn identical; EC vocab unchanged (2,277 / shared 2,132).
  W8 comparison cells: 24,827 changed = 11,270 source_repaired + 13,176 rep_repaired + 381 cache_history_1020; 0 unattributed (PASS).
- F13 rendered to manuscript_figures/out (n 67,038 -> 67,713; median dF1 0.0716-0.0951 -> 0.0727-0.0952, label "+0.07 to +0.10").
- 2026-11-15T22:57:16-06:00 stamp R2L.1: cache unchanged. R2.5 move-aside (logs/r25_moveaside.log). Launched rescore-r25-p2c (160G, 12 workers).
- 2026-11-15T23:00:55-06:00 R2.2_DONE seen (22:37; commit abdfca6). R2.5 selections done (peak 66.4 GiB). Launched rv-rescore-downstream (R2.3 part A: run_downstream_d100.sh + hash-seed proof), 60G.
- 23:15 R2.5 done (rescore-r25-p2c: 1,071 s, cgroup peak 66.5 GiB). Fetch: 40,895 (default) + 350,476 (legacy, D10) new ids.
  depth-probe sample keys IDENTICAL to committed -> k500 hits reused; eval answered 283/320.
- 23:21 R2.3 part A done (rv-rescore-downstream: 1,241 s, 36.5 GiB). check_d100 all: only check 6 fails, on the 10 R1.4 d20 files
  + 5 P2C as-shipped files R2.5 was rewriting (re-stamp in R5.5). compare_three_panels _d100 and _d100_legacy byte-identical under
  PYTHONHASHSEED 0 and 777.
- 2026-11-15T23:21:55-06:00 launched rescore-r25-checks (SNAP20 control + diffs + W8) and rescore-r23b (MAG panel arm, pairings, d20 cross-panel).
- 23:24 rescore-r25-magarm OOM (60G) in mag_prfba_ec.py (Pool(40) workers x ~3.3 GB); mag_direct_feed.py had completed. 2026-11-15T23:31:49-06:00 relaunched mag_prfba_ec + mag_score as rescore-r25-magarm-b (180G).
- 2026-11-16T07:43:30-06:00 all W8 PASS (comparison cells d100/d20, rows: selfrec, ceiling, depth probe, AF gate, MAG arm). Commits c21b029 (R2.4), d354551 (R2.5), 7b14dc3 (AF-gate levels), ec1d0c3 (R2.3). collect: control/collect_run2_local.*. Report logs/RUN2_LOCAL_REPORT.md. Stamp R2L.final unchanged.
- 2026-11-16T07:46:04-06:00 report corrected per review: fresh d20 baseline check (17/32 failing, all re-score outputs incl. R2.3's 2 three_panel files), R2.2 items collected (control/collect_run2_r22.*), mag_panel_comparison 4dp diff, depth-probe wording.
