# Submitting the ModelSEEDpy reversal fix

For `freiburgermsu/ModelSEEDpy`. Branch `fix-gapfill-reversals` (commit e1254ba), one
commit on top of `atp-safe-and-drain-fixes` (c2e31a0). The branch lives in
`~/Documents/ModelSEEDpy`, checked out at `~/Documents/ModelSEEDpy-reversals`.

Upstream (`cshenry/ModelSEEDpy` main and dev) has no `MSBuilder.build_and_gapfill`. The
function exists only in this fork, restored in cb6115f, so the bug and the fix are
fork-specific.

**Submit** (from a shell with push access to the fork):

```bash
cd ~/Documents/ModelSEEDpy
git push origin fix-gapfill-reversals
gh pr create --repo freiburgermsu/ModelSEEDpy --base atp-safe-and-drain-fixes \
  --head fix-gapfill-reversals \
  --title "build_and_gapfill: apply the gapfiller's reversals, not only its new reactions" \
  --body-file ~/Documents/prFBA/region_validation/shared_reference/PR_ModelSEEDpy_gapfill_reversals.md
```

The commit touches only `build_and_gapfill`'s integration step and adds a static helper,
so it should also cherry-pick onto `remove-community-module` (d884a50) if that line needs
it.

The PR body is `PR_ModelSEEDpy_gapfill_reversals.md`.
