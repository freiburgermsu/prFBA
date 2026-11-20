# Run 4 local (FM) progress log

- 11:4x stamp_state --compare R3.final: CACHE UNCHANGED, catalogue current (stamps/R4L.0.json; logs/r4/stamp_R4L0.log).
- 11:46 salterns R2A rollout (user-approved 2026-11-16) found running on FM: scope cap 220G, prep tasks reserve 150 GB,
  launch gate MemAvailable - outstanding - reserve >= 40. Memory decisions: logs/r4/r4_mem_decisions.log.
- 11:54 wobj outputs moved aside (logs/r4/r4_wobj_moveaside.log; control/r4/wobj_moveaside.sha256).
- 11:54-11:55 rescore-r4-wobj-pop (83 s, 4.0 GiB) -> population.json 5d13cf1e; checks pass; marker logs/R4_WOBJ_POPULATION_DONE.
- 11:54-11:59 rescore-r4-wobj-ctl: HEAD on pre-repair inputs reproduces committed population jobs + 3 rows bitwise
  (control/r4/wobj_ctl_head_vs_committed.json); 11-15 FM rerun also bitwise (control/r4/wobj_ctl_fm1115_vs_committed.json).
- 11:56 wobj main run launched (rescore-r4-wobj-run, 20G, 8 workers, MSD e8446e3e).
- 11:56 A-C13 phenotype dirs + logs moved aside; data/phenotype_summary.json copied to pre_repair (logs/r4/r4_pheno_moveaside.log).
- 12:05 ATP pass launched, --workers 3, 30G (strict budget 31.5 GB beside salterns merge).
- 12:09 salterns merge:R2A_A_D2 done; every salterns prep task then skipped for 'memory'.
- 12:22 GF pass launched, --workers 12, 135G (budget 145 GB).
- 13:54 wobj main run done (7,094 s, 12.3 GiB peak); rowcheck: 358 unchanged-job rows bitwise = 10-29 rows; wobj_results agree.
- 13:55 ATP pass stopped at 225 genomes and resumed (rescore-r4-pheno-atp2, 8 workers, 55G); last genome block dropped
  before resume (control/r4/atp_restart/trim_report.json); resume read 224 done, 776 to build.
- 13:58 commit 2520f4f (wobj results).
- 16:12 ATP pass done (resumed run 8,248 s; 43,174 rows; peak 49.8 GiB). Post unit: c12_compare, attribution PASS
  (442 genomes in both samples; 106 changed rows all attributed), manual numbers.
- 16:32 GF pass done (15,018 s; 8,660 rows; 0 gated; peak 88.5 GiB). salterns resumed (tadpole:R2A_A_D2) at 16:32:
  held from 12:09:31 to 16:32 (4.4 h).
- 16:54-17:04 GF post unit: attribution PASS (62 in both, 2,684/2,684 rows identical); sample-diff analysis.
- 17:04 data/phenotype_summary.json := byte copy of phenotype_gf_c13/phenotype_summary.json (8682122e).
- 17:05 R4.2 pass 2 launched (rescore-r42-genus, 50G).
- 17:1x commit 5fac436 (R4.1).
- 18:29 R4.2 pass 2 done: exit 0, C1-C15 pass, bound, arms = pass 1, headline valid. salterns resumed 18:30.
- 19:1x commit bb41947 (R4.2); stamp_state --compare R3.final: CACHE UNCHANGED (stamps/R4L.final.json);
  collect.py over 18 local R4 items -> control/r4/collect_r4_local.{json,md}; report logs/RUN4_LOCAL_REPORT.md.
