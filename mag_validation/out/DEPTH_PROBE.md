# Does post-hoc self-exclusion match a true leave-one-out?

Self-exclusion in both arms is a POST-HOC filter on returned hits, not a re-search
against a database with the source removed. skani ranks with the source present,
reports it at rank 1, and `exclude_hits()` drops it afterwards. A true leave-one-out
would need one sketch per genome (9,997 of them), which is not tractable; the
amplicon arm does the same thing (`exclude_self_record` is documented as "a pure
post-hoc transform of the SAME hits" on a record already truncated to top-20).

## Why the mechanism is benign

skani's `-s` screen is an ABSOLUTE identity threshold (~80%), not relative to the
best hit, and the search ran `--min-af 0`. A 100-ANI self hit therefore cannot make
skani withhold more distant references -- it only CONSUMES A RESULT SLOT. The entire
distortion reduces to: how often did excluded hits fill the top-N?

## Measured on the production run (99,970 queries, -n 100)

| | |
|---|---|
| median references at >=95 ANI per query | 1 (just the source) |
| p95 / p99 | 7 / 45 |
| queries hitting the -n 100 cap | 3,687 (3.69%) |
| capped AND emptied by exclusion | **156 (0.156%)** |

## Direct test: depth 100 vs depth 2000

447 capped queries + 200 random non-capped controls re-searched at `-n 2000` against
the identical sketch, then put through the same exclusion and selection:

| | n | identical selection | survivor count differs |
|---|---|---|---|
| controls (<100 hits returned) | 200 | **200 (100.0%)** | 0 |
| capped at n=100 | 447 | 436 (97.5%) | 437 |

* For any query that did not hit the cap, post-hoc exclusion is **provably identical**
  to a true leave-one-out -- exactly as the absolute-screen argument predicts.
* Capped queries gain a median of 39 extra post-exclusion survivors at depth 2000
  (max 423), but these are more distant than the reference already chosen, so the
  selection changes for only 11 of 447.
* **All 11 changes are abstention -> resolution.** No selection flips to a different
  reference. Depth converts "no answer" into "an answer", never one answer into another.

Scaled to the full run, roughly 90 queries (0.09%) abstain when a deeper search would
have resolved them. The residual bias is one-directional and makes MAG coverage look
slightly WORSE than it is.

## Correction

RESULTS_v2 says "Nothing here is truncated: persistence depth was 100 from the
outset." That is wrong for the 156 capped-and-emptied queries. The accurate statement
is that truncation censors 0.16% of queries here against the amplicon arm's 94.6% of
abstentions -- still a qualitative difference, but not zero. The planned `-n 500`
re-run on BOTH panels addresses it, and the depth-100 vs depth-500 delta will be
reported rather than the depth silently changed.

## How deep is deep enough?

Re-selecting the same 647 probe queries from the depth-2000 hit lists, truncated to
successive depths:

| depth | selection matches depth-2000 | queries resolving |
|---|---|---|
| 100 | 636/647 (98.3%) | 340 |
| **200** | **647/647 (100.0%)** | **351** |
| 500 | 647/647 | 351 |
| 1000 | 647/647 | 351 |
| 2000 | -- | 351 |

The candidate pool exhausts naturally: no probe query returns more than 529 hits even
at `-n 2000`, and none hits the 2000 cap (p50 119, p95 205, p99 292). **Depth 200 is
already fully converged**, so the review's `-n 500` recommendation is sufficient with
margin rather than marginal -- its stated basis ("largest class 304") was not the
binding constraint. The re-runs still use `-n 500`, both for the margin and so both
panels are reported under one setting, but the `-n 100` -> `-n 500` delta is now known
to be the COMPLETE correction, not a partial one.
