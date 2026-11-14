# Hand fixes to the Decagon mapping

Research for issue #55, `docs/spec-burden.md` §4 and §11.1. Run on 2026-11-09. The verdicts are in `data/decagon_overrides.csv`. This note records the evidence behind them and the verdicts on name-check hits that need no override row.

**Status:** reviewed by the user on 2026-11-09 (section 3). All 11 side-effect-identified drugs are joined, and the four judgement calls on dropped rows are approved.

Markers: **[verified: source]** = read from a primary source or a live API call; **[derived]** = computed here from the files and API responses.

## 1. How the sheet was built

`python -m sief_druco.burden.handfix` writes the evidence sheet to `local/data/handfix/` (`evidence.csv`, `name_hits.csv`, `uncovered_units.csv`). Lookups are cached in `local/data/cache/` (PubChem PUG REST and RxNav, live 2026-11-09). [verified: PubChem, RxNav, SIDER 4.1 `drug_names.tsv`]

- **Drugs covered:** the 34 unmatched Decagon drugs, the 17 that reached more than one `IN`, and the Decagon drugs hit by the name check.
- **Name check:** run on the source units outside the 645 in the NHANES main pool (`I` + `P_`) and in MEPS 2015–2017, matched by name (split at `-`). Exact means the name equals a Decagon name, ignoring case. Loose means an `IN` name appears as a whole word inside a Decagon name. The Decagon names are those of `research/nhanes-coverage-pilot.md` §3.3 plus the SIDER name of the same flat CID. Result: 14 exact and 68 loose hits on 48 Decagon drugs. [derived]
- **Evidence per candidate `IN`:** the SIDER name; the Decagon CID's formula and InChIKey next to those of the PubChem record the `IN` name resolves to (same InChIKey, same connectivity meaning the same first InChIKey block, or same formula); the RxNav type and status of the `IN`.
- **Extra evidence for drugs with no structural match:** the Decagon single-drug side effects (`bio-decagon-mono.csv`), searched for terms typical of the candidate drug. This evidence is indirect. It shows what the FAERS reports behind the STITCH ID were about, not what the compound is.

## 2. Result

54 override rows: 41 matched, 10 dropped_in and 3 unmatched. With the overrides, 642 of the 645 drugs reach an `IN`, no `unchecked` alias is left, and 654 `IN`s join a Decagon drug. The only collision is still hydrocortisone (`IN` 5492 → `CID000003640`, 396 labelled pairs, not the butyrate ester `CID000003642`, 26). [derived]

- **Aliases:** 9 racemate-enantiomer, 2 acid/salt (risedronate → risedronic acid; mycophenolate → mycophenolic acid) and 2 prodrug (tenofovir → tenofovir disoproxil; mycophenolate mofetil → mycophenolic acid).
- **Footnote, tenofovir:** NHANES and MEPS "tenofovir" (`IN` 117466) joins Decagon's tenofovir disoproxil, because the sources don't separate TDF from TAF.
- **No epimers or diastereomers are joined:** quinine and quinidine, dexamethasone and betamethasone, doxorubicin and epirubicin, estradiol and alfatradiol, adenosine and vidarabine, mannitol and sorbitol, ursodiol and chenodiol, and the three retinoic acids each keep one `IN`.

## 3. For the user's review

**Dropped `IN`s where the choice of main `IN` is a judgement (all four approved by the user, 2026-11-09):**
- `CID000001065`: SIDER says quinidine, while PubChem's title stereo is quinine's. Kept quinidine.
- `CID000003003`: SIDER says dexamethasone, while PubChem's title is Celestone, a betamethasone brand. Kept dexamethasone.
- `CID000000772`: SIDER says "LMWH" and a synonym is Clexane. Kept enoxaparin and dropped heparin, dalteparin and semuloparin.
- `CID000005538`: SIDER says "retinoic" (acid). Kept tretinoin and dropped isotretinoin and alitretinoin.

**Identified by side effects, joined (user decision, 2026-11-09):** 11 STITCH IDs whose PubChem structure is a fragment, isomer or impurity of a drug, and whose single-drug side effects point to that drug. None of the 11 structures is a marketed drug. The structure is a slip in how STITCH matched the FAERS drug name to a compound, and the labels behind the ID are reports on the named drug.

- **Why join:** the model's only drug input is the mono side-effect profile (`docs/spec-model.md`, no per-drug embedding), so the profile, not the structure, is what the model knows about the drug. Leaving the IDs unmatched would report montelukast, rosuvastatin, allopurinol and lithium as *not in the model*, though the model has their labels.
- **In the spec:** §4 now names these *structure slips* as a rule of their own: the drug the side effects identify is joined when the compound's name or structure points the same way. Each row's `reason` starts with "identity from the Decagon single-drug side effects", and §12 lists the joins as a limitation.
- **Weakest join:** sucralfate. An aluminium sulfate cluster with ulcer side effects would also fit aluminium hydroxide antacids.
- **Mycophenolic acid:** joined with mycophenolate mofetil (68149, prodrug) and mycophenolate (265323, the anion of mycophenolate sodium, acid/salt) as aliases.

| STITCH ID | Structure | Joined `IN` | Side-effect evidence | Source mentions (NHANES pool / MEPS) |
|---|---|---|---|---|
| `CID005311297` | montelukast minus its thio side chain | montelukast 88249 | allergic granulomatous angiitis | 457 / 8,109 |
| `CID004475485` | "Rosuvastatin EP Impurity L" (C22H30 vs C22H28) | rosuvastatin 301542 | CPK increased, muscle injury | 320 / 5,346 |
| `CID000083786` | ring-fusion isomer of allopurinol | allopurinol 519 | Stevens-Johnson syndrome | 228 / 4,266 |
| `CID000028112` | lithium aluminum hydride | lithium 6448 | bipolar I disorder, mania, hypomania | 20 / 854 |
| `CID006398525` | aluminium sulfate cluster | sucralfate 10156 | peptic ulcer, ulcer haemorrhage | 41 / 683 |
| `CID000087177` | maltulose, an isomer of lactulose | lactulose 6218 | ammonia increased, gastric varices | 17 / 466 |
| `CID000508230` | C12H12O6 fragment of mycophenolic acid | mycophenolic acid 7145 | kidney transplant rejection, CMV | 4 + 23 mofetil / 429 mofetil |
| `CID000443871` | 3-isopropyl isomer of nateglinide | nateglinide 274332 | diabetes terms | 5 / – |
| `CID005311048` | cisplatin without platinum | cisplatin 2555 | ototoxicity, febrile neutropenia | – |
| `CID005311039` | "Caspofungin Impurity F" | caspofungin 140108 | candida sepsis, aspergillosis | – |
| `CID005311082` | peptide close to eptifibatide | eptifibatide 75635 | haemorrhages, thrombocytopenia | – |

Three more stay unmatched with no candidate: `CID000027686` (a positional isomer of cromolyn, side effects not telling), `CID000504578` (an unnamed aminoglycoside) and `CID003086258` (an unnamed porphyrin, 1 labelled pair).

## 4. Name-check hits that need no override row

All 82 hits have a verdict. Those on the 34 + 17 drugs are covered by their override rows, along with epinephrine, risedronate and tenofovir. The rest leave the automatic mapping as it is: [derived]

- **Substrings of chemical names:** "acetic acid" (in 19 IUPAC names), "urea" (7), "salicylic acid" (aspirin, diflunisal, mesalamine), "glycerin" (nitroglycerin), "lysine" (glatiramer), "selenium" (selenium sulfide), "everolimus" ("everolimus EP impurity A", a sirolimus synonym). The source drug is a different compound.
- **PENICILLIN** (failed source name, NHANES 9 and MEPS 2 mentions) matches penicillin G (`CID000002349`) and penicillin V (`CID000004730`). It is a source name that reached no `IN`, so it stays *mapping failed* on the source side (§4: no source-side fixes).
- **calcium** (`IN` 1895, MEPS 1,218 mentions) matches `CID000000271` "Calcium Cation", which reached `IN` 1316581 calcium ion. Not joined: the MEPS records are "CALCIUM-VITAMIN D" supplements of an unnamed salt. Calcium carbonate, acetate and citrate are separate `IN`s, so joining one of them would open all of them.
- **vitamin D** (`IN` 11253, MEPS 1,218) matches ergocalciferol (`CID000003249`, SIDER name truncated to "vitamin"). Not joined: "vitamin D" can be D2 or D3.

## 5. Lithium and hydrocodone (§11.1)

- **Lithium:** no Decagon drug has lithium's structure. The only name hit is lithium aluminum hydride (`CID000028112`), whose side effects point to lithium (bipolar I disorder, mania). It is joined to lithium (`IN` 6448, section 3), so lithium is in the model. [derived] Lithium carbonate (`IN` 42351, its salt) is an alias (#59): DDInter names lithium only as "Lithium carbonate", so without it every lithium pair was *not checkable*. NHANES reaches only 6448, so coverage doesn't move. [verified: RxNav properties, DDInter files, 2026-11-14]
- **Hydrocodone** (`IN` 5489; NHANES pool 353 mentions, MEPS 15,503): no name hit, and none of the same-connectivity CIDs of PubChem's hydrocodone record (CID 5284569) is among the 645. Hydrocodone is *not in the model*, so gabapentin + hydrocodone can't be a positive control. [verified: PubChem, live 2026-11-09]
