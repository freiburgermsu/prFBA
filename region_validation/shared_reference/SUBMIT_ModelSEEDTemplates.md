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

**Status (2026-10-27):** the issue is filed as ModelSEED/ModelSEEDTemplates#10. A first
attempt at the PRs failed because this file used a `<you>` placeholder, which bash reads
as input redirection: no fork remote was added, and the PR heads were blank. Every command
below is literal and needs no editing. Do not re-run `gh issue create`.

```bash
cd ~/Documents/ModelSEEDTemplates
git fetch origin main && git log --oneline -1 origin/main     # expect 2294299

# 1. the fork (prints a notice and changes nothing if it already exists)
gh repo fork ModelSEED/ModelSEEDTemplates --clone=false

# 2. the fork as a remote; push the branches (no -u, so the shared clone's tracking
#    config is left alone; pushing never changes its checkout)
git remote get-url fork 2>/dev/null || \
  git remote add fork https://github.com/freiburgermsu/ModelSEEDTemplates.git
git push fork fix-irreversible-reaction-directions fix-opposite-template-directions \
    fix-reverse-penalty-sign negate-all-negative-reverse-penalties \
    db-align-opposite-directions align-reversible-template-directions \
    six-reaction-fix-nadx-variant

# 3. the PRs (the third is optional)
D=~/Documents/prFBA/region_validation/shared_reference
pr() {
  gh pr create --repo ModelSEED/ModelSEEDTemplates --base main \
    --head "freiburgermsu:$1" \
    --title "$(head -1 "$D/$2" | sed 's/^# //')" --body "$(tail -n +3 "$D/$2")"
}
pr fix-irreversible-reaction-directions PR_ModelSEEDTemplates_directions.md
pr fix-opposite-template-directions     PR_ModelSEEDTemplates_tierB.md
pr fix-reverse-penalty-sign             PR_ModelSEEDTemplates_penalty.md

# 4. point issue #10 at the pushed reference branches
gh issue comment 10 --repo ModelSEED/ModelSEEDTemplates --body "The branches this issue \
names are on https://github.com/freiburgermsu/ModelSEEDTemplates: \
fix-irreversible-reaction-directions (21a8084, the five-reaction fix), \
fix-opposite-template-directions (5516611), fix-reverse-penalty-sign (09f6bab), and for \
reference negate-all-negative-reverse-penalties (6430128), db-align-opposite-directions \
(50ed6b4), align-reversible-template-directions (62fcc52) and \
six-reaction-fix-nadx-variant (7d104ac)."
```

The database corrections in the issue belong to ModelSEED/ModelSEEDDatabase: 24
reversed directions, 8 set to `=`, and the stoichiometry, mass-balance, compartment and
name fixes. File them there, or cross-post the issue.
