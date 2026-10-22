# Submitting the ModelSEEDTemplates changes

All branches are local, in `~/Documents/ModelSEEDTemplates`, whose `origin` is
`ModelSEED/ModelSEEDTemplates`. That clone's working tree is on `remove-rxn39860`, which
another session uses. Pushing branches from it does not change the checkout.
`phenotype_arm.MST` reads that working tree.

| Branch | Commit | Parent | Submit as | Description file |
|---|---|---|---|---|
| `fix-irreversible-reaction-directions` | `21a8084` | `main` (`2294299`) | PR | `PR_ModelSEEDTemplates_directions.md` |
| `fix-opposite-template-directions` | `5516611` | `21a8084` | PR | `PR_ModelSEEDTemplates_tierB.md` |
| `fix-reverse-penalty-sign` | `09f6bab` | `21a8084` | PR (optional) | `PR_ModelSEEDTemplates_penalty.md` |
| `negate-all-negative-reverse-penalties` | `6430128` | `21a8084` | push only; the penalty PR cites it | — |
| `db-align-opposite-directions` | `50ed6b4` | `21a8084` | push only; reference | — |
| `align-reversible-template-directions` | `62fcc52` | `50ed6b4` | push only; reference | — |
| `six-reaction-fix-nadx-variant` | `7d104ac` | `main` | push only; the first PR cites it as the NadX alternative | — |
| (issue) | — | — | issue | `ISSUE_ModelSEED_template_direction_audit.md` |

The PRs stack on the five-reaction fix. From a fork they all target `main`, so the
second and third show both commits until the first merges. Their descriptions say so.

## Before pushing

1. **Push prFBA first.** The PR and issue texts cite scripts and data in freiburgermsu/prFBA,
   on branch `hardware-aware-alignment-dispatch`. Push that branch, and make sure the
   maintainers can see it. If the repository is private, attach
   `data/template_direction_audit/` to the issue instead.
2. **Check upstream `main`.** The branches were cut from `main` at `2294299`. If upstream
   has moved, recreate each branch with `scripts/template_direction_audit.py apply` on the
   new `main` rather than rebasing. The Gram templates are single-line JSON, so any
   two edits conflict in git. The apply commands are in the PR texts' Reproduce sections.
3. **Point the issue at the reference branches.** The issue names them by branch. Once they
   are pushed to your fork, add the fork's URL.

## Commands

```bash
cd ~/Documents/ModelSEEDTemplates
git fetch origin main && git log --oneline -1 origin/main     # expect 2294299
git remote add fork https://github.com/<you>/ModelSEEDTemplates.git   # once
git push -u fork fix-irreversible-reaction-directions fix-opposite-template-directions \
    fix-reverse-penalty-sign negate-all-negative-reverse-penalties \
    db-align-opposite-directions align-reversible-template-directions \
    six-reaction-fix-nadx-variant

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

The database corrections in the issue belong to ModelSEED/ModelSEEDDatabase: 24
reversed directions, 8 set to `=`, and the stoichiometry, mass-balance, compartment and
name fixes. File them there, or cross-post the issue.
