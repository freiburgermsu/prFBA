# Template directions against the ModelSEED database: audit results (for reference, not for merge)

We compared every reaction direction in the v7.0 templates (GramNeg, GramPos, Archaea,
Core) with the ModelSEED database's reversibility. The template base was the
five-reaction fix (`fix-irreversible-reaction-directions`, `21a8084`). The database was
ModelSEEDDatabase dev (`078a395f`).

The committed audit file was produced from a fork at `e8446e3e`. That is dev at
`2ebc29ed` plus one unmerged commit correcting 33 stoichiometries and 63 acyl-ACP
charges. Its Biochemistry is otherwise identical to dev `078a395f`'s; the later dev
commits touch only Papers/. Against dev, every reaction is classified the same way.

The results bear on both repositories:
- some template entries are wrong;
- many database directions are wrong;
- applying the database directions wholesale costs growth.

Two branches apply the database directions in full, for reference only:
- `db-align-opposite-directions` (`50ed6b4`): the 41 opposite-direction reactions;
- `align-reversible-template-directions` (`62fcc52`): the complete alignment, stacked on
  the 41.

The curated subsets are in PRs:
- `fix-opposite-template-directions`: the five reactions whose template direction is
  wrong;
- `fix-reverse-penalty-sign`: the related reverse-penalty bias.

## Method

For every template reaction the database marks irreversible (`>` or `<`), we first
checked that the template writes it in the database's orientation. The test compares the
signs of the shared compounds, protons ignored. Every reaction resolved to the same
orientation except two, rxn08355 and rxn09314 (the periplasmic menaquinol DMSO and TMAO
reductases), which we left out. They share no compound in the same compartment with
their database entries: the template puts DMSO/TMAO outside and menaquinol in the
cytoplasm, and the database has the reverse. Their draft directions already match, so
only a reversible `GapfillDirection` remains.

We then classified each of the two direction fields:
- **opposite**: the template allows only the database-forbidden direction;
- **reversible**: the template allows both directions.

No affected reaction has non-default bounds.

| Template | Opposite (reactions; `direction` / `GapfillDirection` entries) | Reversible (reactions; `direction` / `GapfillDirection` entries) | Distinct reactions |
|---|---|---|---|
| GramNeg, GramPos | 41; 39 / 27 | 1,166; 874 / 1,164 | 1,192 |
| Archaea | 41; 39 / 27 | 1,170; 873 / 1,168 | 1,196 |
| Core | 5; 4 / 3 | 41; 8 / 41 | 44 |

15 reactions per genome-scale template (2 in Core) fall in both columns, with one field
opposite and the other reversible.

## The 41 opposite-direction reactions, read against physiology

<details>
<summary>Per-reaction readings (dG'0 in the database's orientation, kcal/mol)</summary>

| Reaction | Name | Database | dG'0 forward, grade | Reading | Confidence | Note |
|---|---|---|---|---|---|---|
| rxn00206 | superoxide:superoxide oxidoreductase | `>` | -16.47, BRONZE | template error | high | superoxide dismutase; the templates allow only O2 + H2O2 -> 2 superoxide |
| rxn00380 | ADP:sulfate adenylyltransferase | `<` | 13.29, GOLD | template error | medium-high | ADP-sulfurylase runs APS + Pi -> ADP + sulfate (sulfite oxidation); database '<' is GOLD and all three estimators agree |
| rxn09566 | pyridoxal oxidase | `<` | 90.17, SILVER | template error | high | pyridoxamine oxidase direction; the database stoichiometry is also wrong (consumes 2 H2O2, releases 0.5 O2; should be PM + O2 + H2O -> PL + NH3 + H2O2) |
| rxn11937 | Nitrogen:(NAD) oxidoreductase (N2O-forming) | `<` | 72.33, SILVER | template error | high | N2O reductase (NAD): the draft direction already allows only N2O reduction, but GapfillDirection '>' points the other way, which leaves the gapfill candidate with no direction (bounds 0/0 in ModelSEEDpy); GapfillDirection '<' makes it addable in the physiological direction |
| rxn16828 | R09982 | `>` | -19.1, BRONZE | template error | medium | aminoacrylate peracid reduction (annotated to RutC; the enzyme is uncertain); the database reaction is mass-imbalanced (H, O) |
| rxn00721 | ATP:anthranilate N-adenylyltransferase | `<` | 40.34, BRONZE | database error | medium | anthranilate adenylation; releases PPi, which pyrophosphatase removes; '=' would also be defensible |
| rxn00837 | inosine-5'-phosphate:NADP+ oxidoreductase (aminating) | `>` | -46.76, BRONZE | database error | high | GMP reductase: GMP + NADPH -> IMP + NH3 + NADP+ |
| rxn00986 | ATP:propanoate adenylyltransferase | `<` | 29.28, SILVER | database error | medium | propanoate adenylation; releases PPi, which pyrophosphatase removes; the acyl-adenylate step itself is endergonic, so '=' would also be defensible |
| rxn01225 | Guanine aminohydrolase | `<` | 52.85, SILVER | database error | high | guanine deaminase: guanine + H2O -> xanthine + NH3 |
| rxn01538 | ATP:thiamine phosphotransferase | `<` | 10.51, BRONZE | database error | medium-high | thiamine kinase (ThiK); the database estimators disagree in sign (eQuilibrator +10.5, group contribution -23.5) |
| rxn01802 | Glutaryl-CoA:(acceptor) 2,3-oxidoreductase (decarboxylating) | `<` | 16.71, SILVER | database error | high | glutaryl-CoA dehydrogenase (decarboxylating), degradation direction |
| rxn01885 | CTP:choline-phosphate cytidylyltransferase | `<` | 29.81, BRONZE | database error | medium | CTP:2-trimethylaminoethylphosphonate cytidylyltransferase (KEGG R02590; the database name says choline phosphate); releases PPi, which pyrophosphatase removes; the only estimate is an unpaired BRONZE dGPredictor value; '=' would also be defensible |
| rxn02569 | ATP:selenide, water phosphotransferase | `<` | 77.15, BRONZE | database error | medium-high | selenophosphate synthetase (SelD); dGPredictor value only |
| rxn03481 | Arbutin 6-phosphate glucohydrolase | `<` | 13.89, SILVER | database error | high | arbutin 6-phosphate hydrolysis |
| rxn03482 | Salicin 6-phosphate glucohydrolase | `<` | 13.88, SILVER | database error | high | salicin 6-phosphate hydrolysis |
| rxn03501 | 1-Carboxyvinyl carboxyphosphonate phosphorylmutase (decarboxylating) | `<` | 50.54, BRONZE | database error | medium-high | decarboxylating carboxyphosphonoenolpyruvate phosphonomutase; releases CO2 |
| rxn03545 | R05230 | `<` | 19.37, BRONZE | database error | high | carbamoyl-phosphate-driven O-carbamoyltransferase (cephamycin C biosynthesis); carbamoyl phosphate (hydrolysis about -12 kcal/mol) makes the transfer downhill, so the positive eQuilibrator and dGPredictor values are implausible; group contribution gives the favourable sign |
| rxn03755 | R05455 | `<` | 61.77, BRONZE | database error | medium | spontaneous lactonisation of OSB-CoA |
| rxn04588 | R06772 | `<` | 25.08, BRONZE | database error | high | carbamoyl-phosphate-driven O-carbamoyltransferase (novobiocin biosynthesis); as rxn03545 |
| rxn05117 | L-aspartate:NADP+ oxidoreductase (deaminating) | `<` | 13.26, SILVER | database error | high | NADP form of rxn05119 (NadX) |
| rxn05119 | L-aspartate:NAD+ oxidoreductase (deaminating) | `<` | 13.27, SILVER | database error | high | L-aspartate dehydrogenase (NadX, EC 1.4.1.21): L-aspartate + NAD+ -> iminoaspartate + NADH is the first step of NAD synthesis in organisms without NadB, pulled by quinolinate synthase |
| rxn07307 | (2R)-O-phospho-3-sulfolactate sulfo-lyase | `>` | -22.36, BRONZE | database error | high | phosphosulfolactate synthase (ComA), coenzyme M biosynthesis: PEP + sulfite -> phosphosulfolactate |
| rxn08167 | L-aspartate oxidase | `<` | 12.62, BRONZE | database error | medium-high | L-aspartate oxidase (NadB) with menaquinone-8: aspartate oxidation is the physiological direction, pulled by quinolinate synthase (NadA); the database '<' rests on a BRONZE eQuilibrator value (+12.6) that is inconsistent with its own NadX value (+13.3), since menaquinone is about 11 kcal/mol the better acceptor; a consistent estimate is near zero (about +2), and group contribution and dGPredictor make no call (-7.1 +/- 13.1, -17.1 +/- 51.5, both '?'); the five-reaction fix sets it '>' |
| rxn08656 | Glycolate oxidase | `<` | 9.87, BRONZE | database error | high | glycolate oxidation by menaquinone (GlcDEF) |
| rxn08679 | GTP amine hydrolysis (spontaneous) | `<` | 65.55, BRONZE | database error | high | spontaneous GTP deamination |
| rxn09992 | 6-phospho-beta-glucosidase | `<` | 12.33, SILVER | database error | high | 6-phospho-beta-glucoside (salicin 6-phosphate) hydrolysis |
| rxn10002 | Heme D synthase | `<` | 36.02, BRONZE | database error | high | heme d synthesis from heme b; the dev database equation (3 O2 + heme <= 2 H+ + heme D) is mass-imbalanced; the fork commit e8446e3e corrects it to O2 + heme + 2 H+ <= heme D |
| rxn12147 | 3D-(3,5/4)-Trihydroxycyclohexane-1,2-dione acylhydrolase (decyclizing) | `<` | 6.24, BRONZE | database error | high | IolD hydrolytic ring opening (myo-inositol degradation) |
| rxn16371 | sulfite:O-phospho-L-serine sulfotransferase (phosphate-hydrolysing, L-cysteate-forming) | `<` | 47.86, BRONZE | database error | high | cysteate synthase, coenzyme M biosynthesis II |
| rxn00288 | Succinate:(acceptor) oxidoreductase | `<` | 20.27, SILVER | reversible or near equilibrium | high | the flavoprotein half-reaction shared by succinate dehydrogenase and fumarate reductase; the database should be '='; succinate oxidation also runs through rxn09272 (succinate + UQ8 -> fumarate + UQ8H2), which is '>' in the draft |
| rxn00379 | ATP:sulfate adenylyltransferase | `<` | 11.34, GOLD | reversible or near equilibrium | high | ATP sulfurylase (Sat) runs forward in sulfate reduction and backwards in sulfur oxidation; the database '<' is a measured GOLD value; '=' is already the genome-scale draft direction |
| rxn01030 | ATP:agmatine N4-phosphotransferase | `<` | 68.0, BRONZE | reversible or near equilibrium | medium-high | agmatine kinase is a phosphagen kinase, near equilibrium |
| rxn01572 | Feruloyl-CoA:quinate O-(hydroxycinnamoyl)transferase | `<` | 9.36, SILVER | reversible or near equilibrium | low | hydroxycinnamoyl-CoA:quinate transferase (HQT); forward in chlorogenic-acid synthesis and reversible in vitro, so '=' rather than either one-way direction |
| rxn03108 | ATP:4-amino-2-methyl-5-phosphomethylpyrimidine phosphotransferase | `<` | 15.59, SILVER | reversible or near equilibrium | medium | HMP-P kinase (ThiD): a phosphoanhydride-to-phosphoanhydride transfer near equilibrium (group contribution 0.0, dGPredictor +2.2); runs forward in cells, pulled by ThiE; the database should be '=' |
| rxn05136 | p-coumaroyl-CoA:quinate hydroxycinnamoyltransferase | `<` | 9.4, SILVER | reversible or near equilibrium | low | as rxn01572 |
| rxn05887 | Hydrogen:NAD+ oxidoreductase | `<` | 14.37, GOLD | reversible or near equilibrium | high | bidirectional NAD-reducing hydrogenase; the database '<' agrees only with the dGPredictor estimate (+14.4) and is graded GOLD despite eQuilibrator (-8.9) and group contribution (-17.9) disagreeing |
| rxn07200 | hydrogen:NADP+ oxidoreductase | `<` | 14.37, SILVER | reversible or near equilibrium | high | bidirectional NADP-reducing hydrogenase; as for rxn05887, the database '<' agrees only with dGPredictor (+14.4; graded SILVER, cross-source disputed), against eQuilibrator (-8.9) and group contribution (-17.9) |
| rxn11944 | Methyltransferase | `>` | -35.82, BRONZE | mis-specified | high | a copy of MetE (EC 2.1.1.14) with SAH/SAM in place of homocysteine/methionine; the database direction would regenerate SAM from SAH without the MetK ATP cost; rxn03052 already covers MetE; remove or keep blocked |
| rxn13150 | 5-formyltetrahydrofolate cyclo-ligase | `<` | 100.91, BRONZE | mis-specified | high | the equation is not 5-formyltetrahydrofolate cyclo-ligase and neither direction is physiological; the database direction would allow O2-driven ATP synthesis; remove or block |
| rxn01349 | acetyl-CoA:hydrogen-sulfide S-acetyltransferase | `<` | 24.35, SILVER | not resolved | - | acetyl-CoA:hydrogen-sulfide S-acetyltransferase |
| rxn09219 | quinol monooxygenase (Ubiquinol-8) | `<` | 18.9, SILVER | not resolved | - | superoxide formation from ubiquinol-8 and O2 |

</details>

**Only five of the 41 are template errors.** For 34 of the other 36, applying the database
direction would block a physiological reaction or open a mis-specified one; 2 are
unresolved. These 41 include three that the five-reaction fix leaves against the
database on purpose:
- **NadX (rxn05119, rxn05117):** its forward reaction is the first step of NAD synthesis in
  organisms without NadB.
- **L-aspartate oxidase (rxn08167):** the fix sets it forward.

**Several database calls rest on one estimator that the others contradict:**
- dGPredictor alone, e.g. the hydrogenases;
- eQuilibrator alone, with the other estimators contradicting it or making no call:
  rxn01538, rxn08656, rxn12147, rxn08167;
- for the phospho-glycoside hydrolases rxn03481, rxn03482 and rxn09992, dGPredictor agrees
  in sign but is undetermined (`?`), and group contribution (about -3) gives the
  opposite sign.

The PPi-releasing transfers are pulled forward by pyrophosphate hydrolysis, which a
standard-state dG does not see.

**Database entries we suggest correcting** (ModelSEED/ModelSEEDDatabase):
- **The 24 "database error" rows:** reverse the irreversible direction, or set `=` where
  the note says so.
- **The 8 "reversible or near equilibrium" rows:** set to `=`.
- **rxn09566:** fix the stoichiometry (it consumes 2 H2O2 and releases 0.5 O2).
- **rxn16828:** fix the mass balance.
- **rxn10002:** the `dev` equation (3 O2 + heme <= 2 H+ + heme D) is mass-imbalanced.
- **rxn08355, rxn09314:** the compartments of DMSO/TMAO and menaquinol look inverted
  (medium-high confidence).
- **rxn01885:** the name says choline phosphate, but the equation is the
  2-trimethylaminoethylphosphonate reaction (KEGG R02590).

**Template entries we suggest changing** (beyond the curated PR):
- **rxn11944 and rxn13150:** remove from the templates, or block them.
  - rxn11944 is a mis-written copy of MetE. Its database direction would regenerate SAM
    from SAH without ATP. In the complete alignment it is open in that direction in 70
    validation drafts, and the gate does not catch it.
  - rxn13150's database direction would allow O2-driven ATP synthesis.
- **Pyruvate kinase, pyruvate formate-lyase and ADP-forming acetate:CoA ligase (rxn00148,
  rxn00157, rxn00172):** `direction` `<` conflicts with `GapfillDirection` `>` in all four
  templates, so gapfilling cannot add them in either direction. (rxn11937 had the same
  conflict; the curated PR fixes it.)

## The 1,166 reversible-in-template reactions

Restricting these to the database direction is the larger part of the alignment.

**The database's evidence for the GramNeg template's 2,038 entries (874 + 1,164):**
- grades: 320 GOLD, 1,002 SILVER, 716 BRONZE;
- median |dG'0|: 14.1 kcal/mol;
- 8% have |dG'0| below 5 kcal/mol.

The database direction agrees in sign with its best dG for 2,030 of the 2,038 entries.
The irreversibility is, in effect, the sign of a standard-state estimate.

We have not curated these 1,151 reactions (those not also in the opposite set) one by
one. The opposite set is a caution, not an error rate: it was selected because the
template curators had already disagreed with the database, and its evidence is weaker
(35 of its 66 GramNeg entries BRONZE, against 716 of 2,038).

**Validation of the complete alignment.** We rebuilt 218 gapfilled models from 20 genomes:
MSATPCorrection, then MSGapfill on complete medium, with the gapfiller's reversals
integrated in every row.

| Templates | Leaking models | New leaks | Models that do not grow | Lost growth (grew under the original templates) | Passing models' growth, median ratio (still growing) | Passing models moved > 10% (of 174) |
|---|---|---|---|---|---|---|
| original v7.0 | 44 | — | 0 | — | — | — |
| five-reaction fix | 0 | 0 | 0 | 0 | 1.000 | 18 |
| **complete alignment** (`62fcc52`) | **0** | **0** | **11** | **11** | **0.81** | **103** |

- **One archaeon loses all 11 of its models.** All come from genome 2837493.4, and
  gapfilling has no solution for them. Two of the 44 leaking models are among them, so
  they count as fixed only because gapfilling failed.
- **Three reactions each stop it growing.** We probed one of its models built with the
  five-reaction fix and reversals integrated: restricting any one of three reactions to
  the database direction stops it growing. The other 113 restrictions do not, one at a
  time.
  - **rxn42902** (methanofuran b from glutamate; database `<`): eQuilibrator +2.6 BRONZE
    and dGPredictor +3.3 support the database, against group contribution -18.1. It is
    probably mis-specified, since it adds two glutamates without ATP (low-medium
    confidence).
  - **rxn03923** (acetyl-phosphate:sulfite S-acetyltransferase, Xsc): the model runs it as
    sulfite + acetyl phosphate -> sulfoacetaldehyde. The database direction is the
    enzyme's physiological, degradative one (high confidence), so this loss removes a
    route that only gapfilling created.
  - **rxn10598** (sulfoethylcysteine synthase): the model needs the forward reduction.
    The database's estimators favour `<` (dGPredictor +32.6 SILVER, group contribution
    +14.3). Whether the forward reduction is the physiological step is unclear to us
    (low confidence).
- **The models that still grow lose about a fifth of their growth.** The median ratio is
  0.81, or 0.75 counting the non-growers as zero. 103 of the 174 passing models move by
  more than 10%: 94 of the 165 that still grow, almost all of them downwards. None of the gapfill solutions matches the five-reaction fix's (median
  Jaccard 0.84).

**Aligning `direction` also restricts gapfilling.** `GapfillingPkg` bounds each candidate
by the template's `lower_bound` and `upper_bound` fields, which these branches set with
`direction`, before applying `GapfillDirection`. One variant was built on the earlier
six-reaction base and rebuilt without reversal integration: it left the reversible set's
`GapfillDirection` entries at `=` and aligned only their `direction`. It gave exactly the
same models as the complete alignment.

## Reproduce

In freiburgermsu/prFBA (branch `hardware-aware-alignment-dispatch`), under
`region_validation/shared_reference/`:
- `scripts/template_direction_audit.py`: `audit [--db <Biochemistry dir>]`, `apply`,
  `table`;
- `scripts/egc_rebuild.py --template-dir <checkout>/templates/v7.0 --variant reversed_only`;
- `scripts/alignment_growth_probe.py`;
- `scripts/template_validation_summary.py`: entries `BC_full+rev`, `PR1v2+rev`,
  `base+rev`;
- `data/template_direction_audit/`: `direction_audit.tsv`, `direction_audit_summary.json`,
  `tierB_reactions.tsv`, `tierB_readings.tsv`, `alignment_growth_probe.json`,
  `validation_summary.json`, `validation_per_model.csv`.

The rebuild inputs (the benchmark's EC sets and the model files) are not in the
repository. ModelSEEDpy was a fork of `main`; the gapfilling bound logic used here is
the same as upstream's.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
