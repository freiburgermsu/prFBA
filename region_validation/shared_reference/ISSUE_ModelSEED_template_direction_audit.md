# Template directions against the ModelSEED database: audit results (for reference, not for merge)

We compared every reaction direction in the v7.0 templates (GramNeg, GramPos, Archaea,
Core) with the ModelSEED database's reversibility. The template base was the six-reaction
fix (`fix-irreversible-reaction-directions`, `7d104ac`), checked against ModelSEEDDatabase
`dev` (`078a395f`). The results bear on both repositories:
- some template entries are wrong;
- many database directions are wrong;
- applying the database directions wholesale costs growth.

Two branches apply the database directions in full, for reference only:
- `db-align-opposite-directions` (`eb891a5`): the 38 opposite-direction reactions;
- `align-reversible-template-directions` (`75f972c`): the complete alignment, stacked on
  the 38.

The curated subsets are in PRs:
- `fix-opposite-template-directions`: the five reactions whose template direction is wrong;
- `fix-reverse-penalty-sign`: the related reverse-penalty bias.

## Method

For every template reaction the database marks irreversible (`>` or `<`), we first
checked that the template writes it in the database's orientation. The test compares the
signs of the shared compounds, protons ignored. Every reaction resolved to the same
orientation except two, rxn08355 and rxn09314 (the periplasmic menaquinol DMSO and TMAO
reductases), which we left out. We then classified each of the two direction fields:
- **opposite**: the template allows only the database-forbidden direction;
- **reversible**: the template allows both directions.

No affected reaction has non-default bounds.

| Template | Opposite (reactions; `direction` / `GapfillDirection` entries) | Reversible (reactions; `direction` / `GapfillDirection` entries) |
|---|---|---|
| GramNeg, GramPos | 38; 36 / 26 | 1,165; 875 / 1,163 |
| Archaea | 38; 36 / 26 | 1,169; 874 / 1,167 |
| Core | 5; 4 / 3 | 41; 8 / 41 |

## The 38 opposite-direction reactions, read against physiology

<details>
<summary>Per-reaction readings (dG'0 in the database's orientation, kcal/mol)</summary>

| Reaction | Name | Database | dG'0 forward, grade | Reading | Confidence | Note |
|---|---|---|---|---|---|---|
| rxn00206 | superoxide:superoxide oxidoreductase | `>` | -16.47, BRONZE | template error | high | superoxide dismutase; the templates allow only O2 + H2O2 -> 2 superoxide |
| rxn00380 | ADP:sulfate adenylyltransferase | `<` | 13.29, GOLD | template error | medium-high | ADP-sulfurylase runs APS + Pi -> ADP + sulfate (sulfite oxidation); database '<' is GOLD and all three estimators agree |
| rxn09566 | pyridoxal oxidase | `<` | 90.17, SILVER | template error | high | pyridoxamine oxidase direction; the database stoichiometry is also wrong (consumes 2 H2O2, releases 0.5 O2; should be PM + O2 + H2O -> PL + NH3 + H2O2) |
| rxn11937 | Nitrogen:(NAD) oxidoreductase (N2O-forming) | `<` | 72.33, SILVER | template error | high | N2O reductase (NAD); the template GapfillDirection allows N2 + NAD+ -> N2O + NADH |
| rxn16828 | R09982 | `>` | -19.1, BRONZE | template error | medium | RutC aminoacrylate peracid reduction; the database reaction is mass-imbalanced (H, O) |
| rxn00721 | ATP:anthranilate N-adenylyltransferase | `<` | 40.34, BRONZE | database error | high | anthranilate adenylation; releases PPi |
| rxn00837 | inosine-5'-phosphate:NADP+ oxidoreductase (aminating) | `>` | -46.76, BRONZE | database error | high | GMP reductase: GMP + NADPH -> IMP + NH3 + NADP+ |
| rxn00986 | ATP:propanoate adenylyltransferase | `<` | 29.28, SILVER | database error | high | propanoate adenylation; releases PPi |
| rxn01225 | Guanine aminohydrolase | `<` | 52.85, SILVER | database error | high | guanine deaminase: guanine + H2O -> xanthine + NH3 |
| rxn01538 | ATP:thiamine phosphotransferase | `<` | 10.51, BRONZE | database error | medium-high | thiamine kinase (ThiK); the database estimators disagree in sign (eQuilibrator +10.5, group contribution -23.5) |
| rxn01802 | Glutaryl-CoA:(acceptor) 2,3-oxidoreductase (decarboxylating) | `<` | 16.71, SILVER | database error | high | glutaryl-CoA dehydrogenase (decarboxylating), degradation direction |
| rxn01885 | CTP:choline-phosphate cytidylyltransferase | `<` | 29.81, BRONZE | database error | high | CTP:2-trimethylaminoethylphosphonate cytidylyltransferase (KEGG R02590; the database name says choline phosphate); releases PPi |
| rxn02569 | ATP:selenide, water phosphotransferase | `<` | 77.15, BRONZE | database error | medium-high | selenophosphate synthetase (SelD); dGPredictor value only |
| rxn03108 | ATP:4-amino-2-methyl-5-phosphomethylpyrimidine phosphotransferase | `<` | 15.59, SILVER | database error | medium-high | HMP-P kinase (ThiD) |
| rxn03481 | Arbutin 6-phosphate glucohydrolase | `<` | 13.89, SILVER | database error | high | arbutin 6-phosphate hydrolysis |
| rxn03482 | Salicin 6-phosphate glucohydrolase | `<` | 13.88, SILVER | database error | high | salicin 6-phosphate hydrolysis |
| rxn03501 | 1-Carboxyvinyl carboxyphosphonate phosphorylmutase (decarboxylating) | `<` | 50.54, BRONZE | database error | medium-high | decarboxylating carboxyphosphonoenolpyruvate phosphonomutase; releases CO2 |
| rxn03545 | R05230 | `<` | 19.37, BRONZE | database error | high | carbamoyl-phosphate-driven O-carbamoyltransferase (cephamycin C biosynthesis); dGPredictor and group contribution disagree in sign |
| rxn03755 | R05455 | `<` | 61.77, BRONZE | database error | medium | spontaneous lactonisation of OSB-CoA |
| rxn04588 | R06772 | `<` | 25.08, BRONZE | database error | high | carbamoyl-phosphate-driven O-carbamoyltransferase (novobiocin biosynthesis); dGPredictor and group contribution disagree in sign |
| rxn07307 | (2R)-O-phospho-3-sulfolactate sulfo-lyase | `>` | -22.36, BRONZE | database error | high | phosphosulfolactate synthase (ComA), coenzyme M biosynthesis: PEP + sulfite -> phosphosulfolactate |
| rxn08656 | Glycolate oxidase | `<` | 9.87, BRONZE | database error | high | glycolate oxidation by menaquinone (GlcDEF) |
| rxn08679 | GTP amine hydrolysis (spontaneous) | `<` | 65.55, BRONZE | database error | high | spontaneous GTP deamination |
| rxn09992 | 6-phospho-beta-glucosidase | `<` | 12.33, SILVER | database error | high | 6-phospho-beta-glucoside (salicin 6-phosphate) hydrolysis |
| rxn10002 | Heme D synthase | `<` | 36.02, BRONZE | database error | high | heme d synthesis from heme b; the template stoichiometry (3 O2) also differs from the database's (1 O2, protons on the other side) |
| rxn12147 | 3D-(3,5/4)-Trihydroxycyclohexane-1,2-dione acylhydrolase (decyclizing) | `<` | 6.24, BRONZE | database error | high | IolD hydrolytic ring opening (myo-inositol degradation) |
| rxn16371 | sulfite:O-phospho-L-serine sulfotransferase (phosphate-hydrolysing, L-cysteate-forming) | `<` | 47.86, BRONZE | database error | high | cysteate synthase, coenzyme M biosynthesis II |
| rxn00288 | Succinate:(acceptor) oxidoreductase | `<` | 20.27, SILVER | reversible in vivo | high | the flavoprotein half-reaction shared by succinate dehydrogenase and fumarate reductase; the database should be '='; succinate oxidation also runs through rxn09272, which this audit does not touch |
| rxn00379 | ATP:sulfate adenylyltransferase | `<` | 11.34, GOLD | reversible in vivo | high | ATP sulfurylase (Sat) runs forward in sulfate reduction and backwards in sulfur oxidation; the database '<' is a measured GOLD value; '=' is already the genome-scale draft direction |
| rxn01030 | ATP:agmatine N4-phosphotransferase | `<` | 68.0, BRONZE | reversible in vivo | medium-high | agmatine kinase is a phosphagen kinase, near equilibrium |
| rxn01572 | Feruloyl-CoA:quinate O-(hydroxycinnamoyl)transferase | `<` | 9.36, SILVER | reversible in vivo | low | hydroxycinnamoyl-CoA:quinate transferase (HQT); forward in chlorogenic-acid synthesis, reversible in vitro |
| rxn05136 | p-coumaroyl-CoA:quinate hydroxycinnamoyltransferase | `<` | 9.4, SILVER | reversible in vivo | low | as rxn01572 |
| rxn05887 | Hydrogen:NAD+ oxidoreductase | `<` | 14.37, GOLD | reversible in vivo | high | bidirectional NAD-reducing hydrogenase; the database '<' rests on dGPredictor alone, contradicted by eQuilibrator (-8.9) and group contribution (-17.9) |
| rxn07200 | hydrogen:NADP+ oxidoreductase | `<` | 14.37, SILVER | reversible in vivo | high | bidirectional NADP-reducing hydrogenase; as rxn05887 |
| rxn11944 | Methyltransferase | `>` | -35.82, BRONZE | mis-specified | high | a copy of MetE (EC 2.1.1.14) with SAH/SAM in place of homocysteine/methionine; the database direction would regenerate SAM from SAH without the MetK ATP cost; rxn03052 already covers MetE; remove or keep blocked |
| rxn13150 | 5-formyltetrahydrofolate cyclo-ligase | `<` | 100.91, BRONZE | mis-specified | high | the equation is not 5-formyltetrahydrofolate cyclo-ligase and neither direction is physiological; the database direction would allow O2-driven ATP synthesis; remove or block |
| rxn01349 | acetyl-CoA:hydrogen-sulfide S-acetyltransferase | `<` | 24.35, SILVER | not resolved | - | acetyl-CoA:hydrogen-sulfide S-acetyltransferase |
| rxn09219 | quinol monooxygenase (Ubiquinol-8) | `<` | 18.9, SILVER | not resolved | - | superoxide formation from ubiquinol-8 and O2 |

</details>

Only five of the 38 are template errors. In the other 33, applying the database
direction would block a physiological reaction or open a mis-specified one.

**Database entries we suggest correcting:**
- **The 22 "database error" rows:** reverse the irreversible direction.
  - Several rest on a single dGPredictor value that eQuilibrator or group contribution
    contradicts in sign.
  - The PPi-releasing transfers are pulled forward by pyrophosphate hydrolysis, which a
    standard-state dG does not see.
- **The 7 "reversible in vivo" rows:** set to `=`.
- **rxn09566:** fix the stoichiometry (it consumes 2 H2O2 and releases 0.5 O2).
- **rxn16828:** fix the mass balance.
- **rxn11944 and rxn13150:** remove from the templates, or keep blocked.

## The 1,165 reversible-in-template reactions

Restricting these to the database direction is the larger part of the alignment.

**The database's evidence for the GramNeg template's 2,038 entries (875 + 1,163):**
- grades: 320 GOLD, 1,000 SILVER, 718 BRONZE;
- median |dG'0|: 14.1 kcal/mol;
- 8% have |dG'0| below 5 kcal/mol.

The opposite-direction set shows how often such calls are wrong, so we have not tried to
curate these 1,165 one by one.

**Validation of the complete alignment.** We rebuilt 218 gapfilled models from 20
genomes, with MSATPCorrection, then MSGapfill on complete medium, and the gapfiller's
reversals integrated:

| Templates | Leaks | New leaks | Models that do not grow | Lost growth (grew before) | Passing models' growth, median ratio (still growing) | Passing models moved > 10% (of 164) |
|---|---|---|---|---|---|---|
| original v7.0 | 44 | — | 10 | — | — | — |
| six-reaction fix | 0 | 0 | 0 | 0 | 1.000 | 13 |
| **complete alignment** (`75f972c`) | **0** | **0** | **11** | **9** | **0.81** | **96** |

- **One archaeon loses every model.** The 9 lost models all come from genome 2837493.4,
  and gapfilling has no solution for them. We probed one of its models built with the
  six-reaction fix: restricting any one of three reactions alone stops it growing.
  - rxn42902 (methanofuran b from glutamate; +2.6 kcal/mol, BRONZE, within the uncertain
    range);
  - rxn03923 (sulfoacetaldehyde from acetyl phosphate and sulfite);
  - rxn10598 (sulfoethylcysteine synthesis).

  Two of the 44 leaking models, both from this archaeon, count as fixed only because
  gapfilling failed.
- **The models that still grow lose about a fifth of their growth** (median ratio 0.81;
  0.75 counting the non-growers as zero), and 96 of the 164 passing models move by more
  than 10%. None of the gapfill solutions matches the six-reaction fix's (median Jaccard
  0.85).
- **Dropping the gapfiller's reversals adds three losses.** It stops three more models
  of a GramNeg genome; with reversals integrated they grow. Those losses are an
  integration artefact, not a template problem.

**Aligning `direction` also restricts gapfilling.** `GapfillingPkg` bounds each candidate
by the template's lower and upper bounds, which follow `direction`, before applying
`GapfillDirection`. Aligning `direction` alone therefore gave exactly the same models as
the complete alignment.

## Reproduce

In freiburgermsu/prFBA (branch `hardware-aware-alignment-dispatch`), under
`region_validation/shared_reference/`:
- `scripts/template_direction_audit.py` (`audit`, `apply`, `table`);
- `scripts/tierc_growth_probe.py`;
- `scripts/template_validation_summary.py`;
- `data/template_direction_audit/`: `direction_audit.tsv`, `direction_audit_summary.json`,
  `tierB_reactions.tsv`, `tierB_readings.tsv`, `tierC_growth_probe.json`,
  `validation_summary.json`, `validation_per_model.csv`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
