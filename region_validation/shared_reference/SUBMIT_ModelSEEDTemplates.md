# Submitting the ModelSEEDTemplates changes

All branches are local, in `~/Documents/ModelSEEDTemplates`, whose `origin` is
`ModelSEED/ModelSEEDTemplates`. That clone's working tree is on `remove-rxn39860`, which
another session uses. Pushing branches from it does not change the checkout.
`phenotype_arm.MST` reads that working tree.

| Branch | Commit | Parent | Submit as | Description file |
|---|---|---|---|---|
| `fix-irreversible-reaction-directions` | `7d104ac` | `main` (`2294299`) | PR | `PR_ModelSEEDTemplates_directions.md` |
| `fix-opposite-template-directions` | `1f45108` | `7d104ac` | PR | `PR_ModelSEEDTemplates_tierB.md` |
| `fix-reverse-penalty-sign` | `4a39553` | `7d104ac` | PR | `PR_ModelSEEDTemplates_penalty.md` |
| `negate-all-negative-reverse-penalties` | `a82ced5` | `7d104ac` | push only; the penalty PR cites it | — |
| `db-align-opposite-directions` | `eb891a5` | `7d104ac` | push only; reference | — |
| `align-reversible-template-directions` | `75f972c` | `eb891a5` | push only; reference | — |
| (issue) | — | — | issue | `ISSUE_ModelSEED_template_direction_audit.md` |

The three PRs stack on the six-reaction fix. From a fork they all target `main`, so the
second and third show both commits until the first merges. Their descriptions say so.

## Before pushing

1. **Check upstream `main`.** The branches were cut from `main` at `2294299`. If upstream
   has moved, rebase each branch onto it and re-run
   `scripts/template_direction_audit.py apply` to confirm the edits still apply cleanly.
2. **Make the cited analysis visible.** The PR and issue texts cite freiburgermsu/prFBA
   (branch `hardware-aware-alignment-dispatch`). If that repository is private, the
   maintainers cannot follow those references. Either make the referenced files visible,
   or attach `data/template_direction_audit/` to the issue.
3. **Point the issue at the reference branches.** The issue names them by branch. Once
   they are pushed to your fork, add the fork's URL.

## Commands

```bash
cd ~/Documents/ModelSEEDTemplates
git fetch origin main && git log --oneline -1 origin/main     # expect 2294299
git remote add fork https://github.com/<you>/ModelSEEDTemplates.git   # once
git push -u fork fix-irreversible-reaction-directions fix-opposite-template-directions \
    fix-reverse-penalty-sign negate-all-negative-reverse-penalties \
    db-align-opposite-directions align-reversible-template-directions

D=~/Documents/prFBA/region_validation/shared_reference
pr() {  # pr <branch> <description file> [extra gh flags]
  gh pr create --repo ModelSEED/ModelSEEDTemplates --base main --head "<you>:$1" \
    --title "$(head -1 "$D/$2" | sed 's/^# //')" --body "$(tail -n +3 "$D/$2")" "${@:3}"
}
pr fix-irreversible-reaction-directions PR_ModelSEEDTemplates_directions.md
pr fix-opposite-template-directions     PR_ModelSEEDTemplates_tierB.md
pr fix-reverse-penalty-sign             PR_ModelSEEDTemplates_penalty.md
gh issue create --repo ModelSEED/ModelSEEDTemplates \
  --title "$(head -1 $D/ISSUE_ModelSEED_template_direction_audit.md | sed 's/^# //')" \
  --body "$(tail -n +3 $D/ISSUE_ModelSEED_template_direction_audit.md)"
```

The database corrections in the issue (22 reversed directions, 7 to `=`, 2
stoichiometry or mass-balance fixes) belong to ModelSEED/ModelSEEDDatabase. File them
there, or cross-post the issue.
