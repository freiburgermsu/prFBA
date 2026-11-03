"""Proposed corrections to the uncommitted Stage E / Stage O manuscript edits.
Each old string must occur exactly once in each variant. Usage:
    python corrections.py <dir_with_both_tex_files> [--apply]
--apply writes the corrected files into that directory (use a scratch copy only)."""
import sys

# Required: factual errors, contradictions, or misattributions.
REQUIRED = [
 ("C1 Stage E union Brier is E1's value, not the estimate's (E3 0.1267)",
  r"and over the union's functions from 0.209 to 0.130, and it raises",
  r"and over the union's functions from 0.209 to 0.127, and it raises"),
 ("C2 the n>=2 within-union F1 gain exceeds the 0.02 SESOI; the parenthetical reads as the SESOI",
  r"below the smallest effect of interest ($+0.025$ at selections of two or more).",
  r"below the smallest effect of interest of 0.02, although at selections of two or more it is $+0.0249$ [$+0.0228$, $+0.0265$], above it."),
 ("C3 Discussion still says 'only in one post-hoc in-silico test' before 'A second post-hoc test'",
  r"and \Cref{eq:objective} only in one post-hoc in-silico test on selections of the controlled benchmark (the \nameref{sec:prfba} section): the weighting cuts",
  r"and \Cref{eq:objective} only in two post-hoc in-silico tests on selections of the controlled benchmark (the \nameref{sec:prfba} section). In the first, the weighting cuts"),
 ("C4 community kinetics: on c2e31a0 the community path raises before any solve (FP3b); the private-copy result is on the fix branch, which leaves it unfixed",
  r"In the released \texttt{prFBA} function (ModelSEEDpy c2e31a0), however, the community package installs this row on a private copy of the community model rather than on the model that is solved, so neither kinCoef nor the \texttt{commkinetics} option affects the solution; no result here uses the community path.",
  r"In the released \texttt{prFBA} function (ModelSEEDpy c2e31a0), however, the community path, too, fails at the carbon-uptake constraint built twice (below), and where that is corrected the community package installs this row on a private copy of the community model rather than on the model that is solved, so neither $\text{kinCoef}$ nor the \texttt{commkinetics} option affects the solution; the fix described below does not change this, and no result here uses the community path."),
 ("C5 abstract: the SESOI-passing configuration includes exchange weight 1; truth is the source genome's draft, not the organism",
  r"whereas with probabilities estimated from each query's 100-hit record and its nearest available family it moves flux off reactions the source organism lacks by more than the declared smallest effect, in multi-genome selections only).",
  r"whereas with probabilities estimated from each query's 100-hit record and its nearest available family, and a lower exchange weight, it moves flux off reactions the source genome's draft does not support by more than the declared smallest effect of interest, in multi-genome selections only)."),
 ("C6 G is a probability (price 1 - p), and the measured effect is up to 0.0002, not nothing",
  r"and pricing gapfilled reactions at 0.084 or 0.2 instead of 0 changed nothing, because",
  r"and giving gapfilled reactions a probability of 0.084 or 0.2 instead of 0 changed truth-absent flux by at most 0.0002, because"),
]

# Recommended: overstatement by omission, imprecision, or stale wording.
RECOMMENDED = [
 ("R1 the face-minimum gain is below the SESOI (source: 'below the SESOI')",
  r"with both face endpoints lowered (minimum 0.0076 [0.0047, 0.0133]);",
  r"with both face endpoints lowered (the minimum by 0.0076 [0.0047, 0.0133], below the smallest effect of interest);"),
 ("R2 the shuffled comparator is itself worse than the released objective (+0.0110 [+0.0031, +0.0222])",
  r"it beat the same costs on shuffled estimates by 0.0282 [0.0195, 0.0382],",
  r"it beat the same costs on shuffled estimates by 0.0282 [0.0195, 0.0382], partly because shuffled estimates raise truth-absent flux (by 0.0110 [0.0031, 0.0222] against the released objective),"),
 ("R3 omitted caution: exchange weight 1 raises total exchange flux, composition not recorded",
  r"so they were held out from the selection but not from the design.",
  r"so they were held out from the selection but not from the design. The lower exchange weight also raises total exchange flux (by 2.5 [1.2, 3.4] units, against a median of about 165), and which exchanges carry the increase was not recorded."),
 ("R4a the filter is the reliability gate plus the family floor, one hit per genome, not coverage alone",
  r"the identity-weighted share of the coverage-qualified hits among a query's first 100 that carry it",
  r"the identity-weighted share of the hits among a query's first 100 that pass the reliability gate and the family floor (one per genome) and carry it"),
 ("R4b the clade rule has a third case: no prior when the order also fails (7,591 of 57,021 cells)",
  r"(pseudo-count 2; the order when the family is unnamed or has fewer than three reference genomes);",
  r"(pseudo-count 2; the order when the family is unnamed or has fewer than three reference genomes, and no prior, for 13\% of cells, when the order fails too);"),
 ("R5 Stage E numbers are over all 57,021 scored cells, not the multi-reference cells of the preceding sentence",
  r"Against $p_{\text{func}}$ it lowers the Brier score over the whole EC vocabulary",
  r"Against $p_{\text{func}}$, over all 57{,}021 scored cells (49{,}048 of them single-genome selections), it lowers the Brier score over the whole EC vocabulary"),
 ("R6 p_func was also the released arm of the second test",
  r"used as the weight of \Cref{eq:objective} in one post-hoc in-silico test, it lowers",
  r"used as the weight of \Cref{eq:objective} in the first post-hoc in-silico test, it lowers"),
 ("R7a 167 of 360 is not half the amplicons; groups split 121/80 by hash",
  r"A configuration was selected on the holdout groups of one half of the amplicons (167) and tested on the other (193).",
  r"The holdout groups were split in two by a hash of their identifiers, and a configuration was selected on one half (121 groups, 167 amplicons) and tested on the other (80 groups, 193 amplicons)."),
 ("R7b Discussion: give the realised split",
  r"chosen on half of the holdout groups, truth-absent flux fell by 0.0172 [0.0118, 0.0221] against the released objective on the other half,",
  r"chosen on 121 of the 201 holdout groups, truth-absent flux fell by 0.0172 [0.0118, 0.0221] against the released objective on the other 80,"),
 ("R7c Limitations: give the realised split",
  r"and the second's configuration chosen on half of the holdout groups,",
  r"and the second's configuration chosen on 121 of the 201 holdout groups,"),
 ("R8a 'weight 0' for gapfilled reactions reads as free flux; it is p = 0, cost 1",
  r"the weight of gapfilled reactions (0 as released, 0.084,",
  r"the probability given to gapfilled reactions (0 as released, so full cost; 0.084,"),
 ("R8b same",
  r"gapfilled weight 0 and exchange weight 1",
  r"gapfilled probability 0 and exchange weight 1"),
 ("R9 the estimate uses the whole 100-hit record and a clade prior, not only the reference set",
  r"the objective can carry more of the reference set's signal:",
  r"the objective can carry more of the alignment record's signal:"),
 ("R10 burden 'at most one reaction' is the controlled pipeline; as shipped it is 1.5",
  r"and moved the method arms' medians by at most 0.004 and their gapfill burden by at most one reaction.",
  r"and moved the method arms' medians by at most 0.004 and, in the controlled pipeline, their gapfill burden by at most one reaction (1.5 as shipped)."),
 ("R11 sign convention: every other reduction in the edits is a positive magnitude",
  r"($-0.0311$ [$-0.0322$, $-0.0297$])",
  r"(by 0.0311 [0.0297, 0.0322])"),
 ("R13 the first test's description lists five weight sets; the exchange-weight-1 arm is never introduced",
  r"whose exploratory exchange-weight-1 arm had been seen when the factors were chosen,",
  r"whose exploratory exchange-weight-1 arm, not reported above, had been seen when the factors were chosen,"),
 ("R12 pre-existing, now sharper: the new passage below compares the weighting with an alternative",
  r"and no number below should be read as evidence for or against the weighting.",
  r"and no selector-arm number below should be read as evidence for or against the weighting."),
]

if __name__ == "__main__":
    d = sys.argv[1].rstrip("/") + "/"
    files = [d + "prFBA.tex", d + "prFBA_andrew_voice.tex"]
    text = {f: open(f).read() for f in files}
    ok = True
    for tag, old, new in REQUIRED + RECOMMENDED:
        for f in files:
            n = text[f].count(old)
            if n != 1:
                ok = False
                print("COUNT", n, f.split("/")[-1], tag)
    # apply in order, re-checking uniqueness after earlier replacements
    for f in files:
        t = text[f]
        for tag, old, new in REQUIRED + RECOMMENDED:
            n = t.count(old)
            if n != 1:
                ok = False
                print("SEQUENTIAL COUNT", n, f.split("/")[-1], tag)
                continue
            t = t.replace(old, new)
        text[f] = t
    print("all old strings unique in both files, also sequentially:", ok,
          "|", len(REQUIRED), "required +", len(RECOMMENDED), "recommended")
    if "--apply" in sys.argv and ok:
        for f, t in text.items():
            open(f, "w").write(t)
        print("applied to", d)
