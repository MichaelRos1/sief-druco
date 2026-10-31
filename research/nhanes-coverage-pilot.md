# Coverage pilot: NHANES 2017–2018 → Decagon's 645 drugs

Research for issue #16 (map #1), points `4.2`–`4.3` of `PROJECT_CONTEXT.md`. Run on 2026-10-29 on NHANES 2017–2018 (cycle J). It uses the reverse route decided on #12 (`research/nhanes-stitch-mapping.md`). It records numbers and decides nothing.

**Status: done for #16.** Every number below was computed. The name-check hits and the two RxNorm questions were checked by hand on 2026-10-29 (section 7.1). Two hand checks are left: the 34 unmatched Decagon drugs and the 17 Decagon drugs that reached more than one ingredient. They belong to "Handling failed drug-name mappings" (#26, decided 2026-10-31: Decagon-side hand fixes only, by an alias rule, done when the mapping pipeline is built). Section 7 bounds how much they could move the numbers. The numbers in this note are **not** reclassified after the hand checks.

Markers: **[verified: source]** = read from a primary source or a live API call; **[derived]** = computed here from the files and API responses; **[unverified]** = not checked.

---

## 1. Short answer

Units are RxNorm ingredients (`IN`). Combination products are split into their ingredients. Weights are `WTINT2YR`. "Mapping failed" covers both sides: the NHANES drug never reached an ingredient, or a Decagon drug never reached RxNorm (found by the name check). [derived]

| Measure | Covered (both / all in the 645) | Not in the model | Mapping failed |
|---|---|---|---|
| Decagon drugs reaching an RxNorm `IN` | **611 of 645** (94.7%) | — | 34 |
| Distinct NHANES ingredients (569) | 353 (62.0%) | 190 (33.4%) | 26 (25 NHANES side, 1 Decagon side) |
| Ingredient mentions (15,560), unweighted | **84.0%** | 14.8% | 1.2% |
| Ingredient mentions, weighted | **82.6%** | 16.2% | 1.2% |
| Participants with ≥ 2 ingredients (2,941): **all** mapped, unweighted / weighted | **50.7% / 47.6%** | 46.8% / 50.0% | 2.6% / 2.4% |
| Same 2,941: **at least 2** mapped, unweighted / weighted | **88.8% / 86.1%** | — | — |
| Co-prescribed pairs (42,246 person-pairs), unweighted | **71.0%** | 28.3% | 0.7% |
| Co-prescribed pairs, weighted | **69.2%** | 30.1% | 0.7% |

- Coverage is limited by drugs that are **not in the model**, not by mapping failures. Mapping failures are about 1% at every level. [derived]
- The two checks carried over from #11: `RXDDRGID` codes keep the same name between 2015–16 and 2017–18 (547 shared codes, no conflict), and `SDMVSTRA` values do **not** overlap between the two cycles (119–133 vs 134–148). [derived]

## 2. Inputs

| File | Source | SHA-256 |
|---|---|---|
| `RXQ_RX_J.xpt` (9,627,680 B) | already in `local/data/raw/nhanes/` (downloaded 2026-10-28) | `6196e6753f9bf23b47456f4fcc2d32de05f89388b512594d02b0735949165545` |
| `DEMO_J.xpt` (3,412,720 B) | `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/DEMO_J.xpt`, 2026-10-29 | `c0b46e0345ea19404928656277c8b0d10b0cca348a9b2fe4fc3c67e8b7ee73ec` |
| `RXQ_RX_I.xpt` (9,629,600 B) | `…/Public/2015/DataFiles/RXQ_RX_I.xpt`, 2026-10-29 | `a148b3c38e6823f3750b5835385a85ac1d65aeeefd8d8b351ba8ad9c9f6095aa` |
| `DEMO_I.xpt` (3,756,480 B) | `…/Public/2015/DataFiles/DEMO_I.xpt`, 2026-10-29 | `c9297c6c37ae8f78f29be9568fa2a03cf3b112616a39afee04030fc775a66a0d` |
| `bio-decagon-combo.csv` | `local/data/raw/decagon/` (see `research/decagon-files.md`) | — |
| SIDER 4.1 `meddra_all_se.tsv.gz` | `local/data/scratch/severity-weights/data/` (used for the stereo-only list only) | — |

- APIs: RxNav (data version `05-Oct-2026`, API `3.1.356`, read from `/REST/version.json` on 2026-10-29) and PubChem PUG-REST / PUG-View, 2026-10-29. [verified: live]
- `.xpt` files were read with `pandas.read_sas` (no `pyreadstat` needed). [derived]
- NHANES 2017–18 counts: 9,254 participants in `DEMO_J`, all with positive `WTINT2YR`. `RXQ_RX_J` has 19,643 rows, of which 14,434 have `RXDUSE` = 1 (a drug was reported), from 4,045 participants. Of the 14,434, 134 rows have no drug code (names `55555` / `77777` / `99999`). [derived]

## 3. Method

Scripts are in `local/data/scratch/coverage-pilot/` (git-ignored). Every API response is cached on disk; failures go to `api_failures.log`.

### 3.1 Decagon side (`decagon_map.py` → `decagon_rxnorm.csv`)

For each of the 645 STITCH IDs, strip `CID` and leading zeros to get a PubChem CID. Then try, in order, stopping at the first step that reaches RxNorm:

1. **Exact CID:** UNIIs from the CID's PubChem synonyms → RxNav `rxcui.json?idtype=UNII_CODE` → `related.json?tty=IN`.
2. **Parent CID** (`cids_type=parent`) → same.
3. **Same-connectivity CIDs** of the parent (`cids_type=same_connectivity`) → same. Every `IN` reached is kept. More than one `IN` is flagged for a hand check.
3b. **New step.** Some CIDs have no title, no synonyms, no parent and no same-connectivity list in PubChem, but they do have a structure (e.g. CID 143, 2022, 2182). For these, the same-connectivity CIDs were found from the CID's SMILES with `compound/fastidentity/smiles/cids?identity_type=same_connectivity`. [verified: live, 2026-10-29]
4. **Name match:** PubChem title → RxNav `rxcui.json?name=…&search=2` → `IN`.
5. **`approximateTerm`** on the title (rank 1). **New rule:** a hit is accepted only if the `IN` name also appears among the CID's own PubChem title and synonyms. Unconfirmed hits are rejected and logged. The rule was added because rank-1 hits on IUPAC titles were often wrong. For example, the cyanocobalamin records (CID 2891, 4479097) gave `mavacamten`. CID 2676, whose title describes a methoxymethyl-pyridine statin, gave `rosuvastatin`. [derived]

**PubChem synonym format has changed. [verified: live, 2026-10-29]** The #12 notes found UNIIs as `UNII-<code>` synonyms. Some records now list the UNII only as a bare code. Aztreonam (CID 5742832) has `G2B4VE5GH8` among its 126 synonyms and no `UNII-` entry. PUG-View `heading=UNII` lists `G2B4VE5GH8` for it. Metoprolol (4171) and ampicillin (6249) still carry both forms. So bare 10-character codes (at least one letter and one digit) were also tried in RxNav. A bare code that reached RxNorm was kept only if PUG-View lists it as that CID's UNII. 40 bare-code RxNorm hits were rejected by this check. Examples: CID 2153 → `C137DTR5RG` theophylline, and CID 33613 → `804826J2HU` amoxicillin, where the code is in the synonyms but not in the CID's UNII section. [derived]

### 3.2 NHANES side (`nhanes_map.py` → `nhanes_codes_J.csv`)

For each of the 666 distinct `RXDDRGID` codes in `RXQ_RX_J`:

- Codes starting `a`, `c` or `h` are logged as unmappable by code.
- `GN` + code → `rxcui.json?idtype=MMSL_CODE` (reusing `rxnav_nhanes_2017.json` unchanged). If that returns nothing, try `allsrc=1`. If that returns nothing too, split `RXDDRUG` at `;` and name-match each part (`search=2`).
- Every RxCUI → `related.json?tty=IN`. If that is empty → `historystatus.json` → `remappedConcept` → `IN`.

### 3.3 Join and units (`coverage.py`)

- **Unit:** an RxNorm `IN`. A combination row counts once for each of its ingredients. An NHANES code that never reached an `IN` stays one unit, and each uncoded row is its own unit. Within a participant, units are de-duplicated: two metformin products count as one.
- **Unit status:**
  - `in645`: the `IN` is in the Decagon table.
  - `not_in_model`: an `IN` was reached but is not in the table.
  - `failed_nhanes`: no code; an `a` / `c` / `h` code; or no `IN` reached.
  - `failed_decagon`: an `IN` was reached and is not in the table, but its name exactly matches (case-insensitive) a Decagon drug's name. Decagon names are the PubChem title, the first 30 synonyms of the exact and parent CID, and the RxNorm `IN` name. For unmatched Decagon drugs, the first 10 synonyms of each same-connectivity CID are added.
- **Pairs:** all unordered pairs of distinct units within a participant. A pair is "covered" if both units are `in645`. It is "mapping failed" if at least one unit failed. Otherwise it is "not in the model".
- **Weights:** each mention, participant or person-pair carries the participant's `WTINT2YR`. Weighted shares are ratios of weight sums. No variance estimates were computed (section 10).

## 4. Decagon side: 611 of 645 reach an RxNorm ingredient

| Step | Drugs |
|---|---|
| 1 exact CID UNII | 321 |
| 2 parent CID UNII | 7 |
| 3 same-connectivity UNII | 234 |
| 3b same-connectivity by SMILES | 30 |
| 4 name match (title) | 18 |
| 5 `approximateTerm`, name-confirmed | 1 (CID 24011 → selenium sulfide) |
| **Matched** | **611** (594 to one `IN`, 17 to more than one) |
| Unmatched: only an unconfirmed `approximateTerm` hit | 28 |
| Unmatched: UNII found, no RxNorm concept | 2 (CID 9052, 83786) |
| Unmatched: no UNII, title not found | 4 (CID 60754 "Omniscan", 5311048, 5487301, 6398525) |

All numbers are [derived] from `decagon_rxnorm.csv`. As #12 expected, the fallbacks did much of the work: 290 of 611 matches (47%) came from steps 2–5. [derived]

**The 28 unconfirmed `approximateTerm` hits** (title → rejected rank-1 `IN`), left for a hand check. Some look plausible from the title (`Oxaliplatin (Eloxatin)` → oxaliplatin, `Ridauran` → auranofin, `Calcium chloride 6-hydrate` → calcium chloride). Others are clearly wrong (cobalamin → mavacamten). None of these is counted as matched. [derived]

CID 85 → levocarnitine; 2283 → phentolamine; 2622 → cefepime; 2650 → *Stellaria media* extract; 2676 → rosuvastatin; 2891 → mavacamten; 4609 → oxaliplatin; 4675 → 1-octen-3-yl acetate; 4812 → vitamin K1; 24486 → calcium chloride; 27686 → cromolyn; 28112 → aluminum hydride; 38904 → carboplatin; 41774 → streptozocin; 87177 → lactulose; 443871 → cinnamic aldehyde; 504578 → migalastat; 508230 → mycophenolic acid; 3085017 → colestipol; 3086258 → proxigermanium; 3468412 → deserpidine; 4475485 → rosuvastatin; 4479097 → mavacamten; 5311039 → aniline; 5311082 → sulfadiazine / sulfapyridine / sulfathiazole; 5311297 → isopropyl alcohol; 5473385 → thioguanine; 6918453 → auranofin.

**17 drugs reached more than one `IN`.** Fifteen came from step 3 and one (CID 143) from step 3b; there stereoisomers or epimers collapse. CID 772 came from step 1, where several UNIIs on the one record lead to several heparins. All are kept in the join (the #12 route says "keep the hits"); section 7 gives the effect. [derived]

| STITCH ID | `IN`s reached |
|---|---|
| CID000000143 | leucovorin; levoleucovorin |
| CID000000191 | adenosine; vidarabine |
| CID000000206 | glucose; alpha-D-glucopyranose |
| CID000000450 | estradiol; alfatradiol |
| CID000000453 | mannitol; sorbitol |
| CID000000772 | heparin; enoxaparin; dalteparin; semuloparin |
| CID000001065 | quinidine; quinine |
| CID000001206 | methamphetamine; levo-methamphetamine |
| CID000001690 | doxorubicin; epirubicin |
| CID000003003 | betamethasone; dexamethasone |
| CID000003410 | formoterol; arformoterol |
| CID000003661 | atropine; hyoscyamine |
| CID000004158 | methylphenidate; dexmethylphenidate |
| CID000004542 | levonorgestrel; norgestrel |
| CID000005538 | isotretinoin; tretinoin; alitretinoin |
| CID000005645 | chenodeoxycholate; ursodeoxycholate |
| CID000060612 | dexmedetomidine; medetomidine |

- Two Decagon drugs, CID000003640 and CID000003642, both map to `IN` 5492 hydrocortisone. The 611 matched drugs give 630 distinct `IN`s. [derived]

**The 10 stereo-only IDs handed over from #12.** These are numbers that occur in SIDER only as stereo IDs. None of their SIDER flat partners is also among the 645. [derived]

| STITCH ID | PubChem title | Result |
|---|---|---|
| CID003002190 | Telithromycin | step 1 → telithromycin |
| CID005281104 | Paricalcitol | step 1 → paricalcitol |
| CID005282226 | Travoprost | step 1 → travoprost |
| CID005311027 | Bimatoprost | step 1 → bimatoprost |
| CID005311181 | Iloprost | step 1 → iloprost |
| CID005329102 | Sunitinib | step 1 → sunitinib |
| CID005481350 | Tenofovir Disoproxil | step 4 → `tenofovir disoproxil` (not the `tenofovir` `IN` that NHANES reaches; see section 7.1) |
| CID005493444 | Aliskiren | step 1 → aliskiren |
| CID006436173 | Rifaximin | step 1 → rifaximin |
| CID006918453 | Ridauran | **unmatched** (`approximateTerm` → auranofin, not confirmed by name) |

## 5. NHANES side: 641 of 666 codes reach an ingredient

| Outcome | Codes | Mentions (rows) |
|---|---|---|
| Step 1 `MMSL_CODE` | 632 | 14,214 |
| Step 2 `allsrc=1` + `historystatus` remap (d03457 aspirin / butalbital / caffeine) | 1 | 1 |
| Step 3 name match on `RXDDRUG` (pre-mixed insulins d04374, d04510, d04839; pancrelipase d01002; three tenofovir combinations; d08128) | 8 | 34 |
| **Reached an `IN`** | **641** | **14,249** |
| RxCUI found (`allsrc=1`) but `NotCurrent` with no remap → no `IN` | 8 | 25 |
| `c` codes | 13 | 19 |
| `h` codes (`h00015` ethinyl estradiol; gestodene, `h00024` antihypertensives unspecified, `h00039` etoricoxib) | 3 | 6 |
| `a` code (`a54115` oxygen) | 1 | 1 |
| No code (`55555` / `77777` / `99999`) | — | 134 |

The eight `NotCurrent` codes are penicillin (d00116, 6 mentions), fluoride topical (9), ammonium lactate topical (3), cannabis (2), conjugated estrogens topical (2), terconazole topical, azelaic acid topical and loteprednol ophthalmic (1 each). [derived] The #12 note counted 20 codes with no RxCUI. The `allsrc=1` and name steps fixed 9 of them; the other 11 are the 8 `NotCurrent` codes plus the 3 `h` codes (#12 counted `h` codes under "no RxCUI"). [derived]

## 6. Coverage in detail

**Distinct drugs. [derived]**

| Level | Total | All in 645 | Some in 645 (combinations) | Not in the model | Failed, NHANES side | Failed, Decagon side |
|---|---|---|---|---|---|---|
| Ingredients | 569 | 353 | — | 190 | 25 | 1 (risedronate) |
| Multum codes (products) | 666 | 421 | 48 | 171 | 25 | 1 |

**Mentions. [derived]**

| Status | Product rows (14,434) | share | weighted | Ingredient mentions (15,560) | share | weighted |
|---|---|---|---|---|---|---|
| all `IN`s in 645 | 12,084 | 83.7% | 82.3% | 13,064 | 84.0% | 82.6% |
| some `IN`s in 645 | 460 | 3.2% | 3.8% | — | — | — |
| not in the model | 1,702 | 11.8% | 12.6% | 2,308 | 14.8% | 16.2% |
| failed, NHANES side, coded | 51 | 0.35% | 0.46% | 185 (incl. 134 uncoded) | 1.2% | 1.2% |
| failed, NHANES side, no code | 134 | 0.93% | 0.81% | (in the row above) | | |
| failed, Decagon side | 3 | 0.02% | 0.01% | 3 | 0.02% | 0.01% |

**Participants. [derived]**
- 4,045 participants reported at least one drug: 48.96% of the weighted population.
- 2,941 have at least 2 distinct ingredient units. Counted by distinct Multum codes instead, the number is 2,719.

| Among the 2,941 | n | share | weighted |
|---|---|---|---|
| all units in the 645 | 1,490 | 50.7% | 47.6% |
| not all: at least one mapping failure | 75 | 2.6% | 2.4% |
| not all: only not-in-model units | 1,376 | 46.8% | 50.0% |
| at least 2 units in the 645 | 2,612 | 88.8% | 86.1% |

**Co-prescribed pairs. [derived]**

| Pair status | Person-pairs | share | weighted | Distinct pairs* |
|---|---|---|---|---|
| both in the 645 | 29,977 | 71.0% | 69.2% | 8,884 |
| mapping failed (≥ 1 side) | 306 | 0.7% | 0.7% | 139 |
| not in the model (no failed side) | 11,963 | 28.3% | 30.1% | 5,789 |

\* Distinct pairs leave out pairs that include an uncoded row (14,812 distinct pairs in total).

- **Extra:** of the 8,884 distinct covered pairs, 8,009 are labelled pairs in `bio-decagon-combo.csv` (for at least one STITCH ID behind each `IN`). They make up 95.6% of the weighted covered person-pairs. [derived]

## 7. Name check and what could still move

**Exact name check** (rule in section 3.3). Two hits among all unmapped NHANES units. [derived]
- `risedronate` (`IN` 73056, 3 mentions): Decagon CID 5245 "Risedronate" reached RxNorm `IN` 55685 `risedronic acid`, a different `IN`. It is counted as **failed, Decagon side**.
- `PENICILLIN` (d00116, 6 mentions): matches Decagon CID 2349 "Penicillins", which reached `penicillin G`. The NHANES code never reached an `IN`, so it stays **failed, NHANES side**.

**Loose name check:** the NHANES `IN` name occurs as a whole word inside a Decagon name. These hits are listed for a hand check and **not** reclassified. [derived]

| NHANES `IN` | Mentions | Weighted prevalence | Decagon hit |
|---|---|---|---|
| rosuvastatin | 104 | 1.31% | CID004475485, unmatched, title "Rosuvastatin EP Impurity L" |
| nystatin | 10 | 0.11% | CID003468412, unmatched, synonym "Nystatin (JP18/USP/INN)" |
| tenofovir | 9 | 0.12% | CID005481350, matched to `tenofovir disoproxil`, a different `IN` |
| sevelamer | 8 | 0.05% | CID003085017, unmatched, synonym "Sevelamer (INN)" |
| lactulose | 3 | 0.05% | CID000087177, unmatched, synonym "LACTULOSE [INN]" |
| acarbose | 3 | 0.03% | CID000041774, unmatched, synonym "Acarbose (JAN/USAN/INN)" |
| epinephrine | 3 | 0.04% | CID000000838 "(+-)-Adrenaline", matched to `racepinephrine` |
| ceftazidime | 1 | 0.05% | CID000002650, unmatched, title "Ceftazidime, Antibiotic for Culture Media Use Only" |
| lithium | 8 | 0.07% | only inside "Lithium aluminum hydride" (CID 28112) |
| urea | 1 | 0.00% | only inside chemical names of sulfonylureas, carmustine and others |

**Bounds** (ingredient mentions, of 15,560): [derived]
- **Toward "failed, Decagon side":** the first eight rows of the table add up to 141 mentions (0.9%). These would move from "not in the model" to "mapping failed, Decagon side". The lithium and urea hits are substrings of other compounds' names.
- **Away from "in the 645":** 231 mentions (1.5%) count as in the 645 only because a Decagon drug reached more than one `IN`. The largest are formoterol 57, methylphenidate 53, estradiol 30, dexmethylphenidate 19, levonorgestrel 18 and dexamethasone 15. If those extra `IN`s were dropped, these mentions would become "not in the model".
- **Limit of the name check:** it compares against a limited name list per Decagon drug (section 3.3). A common NHANES drug with no hit, such as allopurinol or montelukast, is therefore "not found in these names", not "proven absent from Decagon". [unverified: absence]

### 7.1 Hand check of the name hits (2026-10-29)

Evidence used:
- PubChem properties (title, InChIKey, formula) of the Decagon CID, compared with the PubChem record the drug name resolves to;
- the drug name SIDER 4.1 gives the same STITCH flat ID (`drug_names.tsv`);
- RxNav properties and history status of the RxNorm `IN`s.

Script: `local/data/scratch/coverage-pilot/handcheck/check.py`. Same connectivity means the first 14 characters of the InChIKey are equal (same molecule apart from stereochemistry). [verified: PubChem PUG REST, RxNav, SIDER `drug_names.tsv`, live 2026-10-29]

| NHANES `IN` | Mentions | Decagon CID | Verdict | Evidence |
|---|---|---|---|---|
| rosuvastatin | 104 | 4475485 | **Not the same compound** | Formula C22H30FN3O6S vs rosuvastatin C22H28FN3O6S (CID 446157); title "Rosuvastatin EP Impurity L". SIDER's rosuvastatin is CID 129806, which is **not** among the 645. Stays "not in the model". |
| nystatin | 10 | 3468412 | Likely the same | Formula C47H75NO17 equals nystatin's; synonym "Nystatin (JP18/USP/INN)". The connectivity differs from PubChem's nystatin records (drawn differently). No SIDER entry. |
| tenofovir | 9 | 5481350 | Same active drug (prodrug ester) | The Decagon CID is tenofovir disoproxil. NHANES reaches `IN` 117466 `tenofovir`, whose only PIN is `tenofovir anhydrous`. NHANES combination names ("emtricitabine; tenofovir") don't say whether the disoproxil or the alafenamide form was taken. |
| sevelamer | 8 | 3085017 | Same | SIDER name "sevelamer". The PubChem record is the polymer's two monomers (C6H12ClNO). |
| lactulose | 3 | 87177 | **Unclear** | PubChem title "Maltulose", an isomer of lactulose. It has the synonym "LACTULOSE [INN]". Connectivity differs from lactulose (CID 11333). No SIDER entry. |
| acarbose | 3 | 41774 | Same | SIDER name "acarbose"; formula C25H43NO18 equals acarbose's. |
| epinephrine | 3 | 838 | Same (racemate vs active enantiomer) | SIDER name "epinephrine"; same connectivity as L-epinephrine (CID 5816). The Decagon side reached `IN` `racepinephrine`, NHANES reached `IN` 3992 `epinephrine`. |
| ceftazidime | 1 | 2650 | Same | SIDER name "ceftazidime"; same connectivity as ceftazidime (CID 5481173). |
| risedronate | 3 | 5245 | Same (acid vs anion) | SIDER name "risedronate". RxNorm keeps `IN` 73056 `risedronate` and `IN` 55685 `risedronic acid` as two separate active concepts; neither is remapped to the other. |

- **Effect if the 8 "same / likely same" rows were joined:** 40 mentions (0.26% of 15,560) would move into "in the 645". Rosuvastatin (104 mentions, 1.31% weighted prevalence) stays outside the model. Lactulose (3 mentions) is unresolved. [derived]
- These are facts. Whether to join them (and by which rule: same connectivity, prodrug to active moiety, acid to salt) is for the map's open item "Handling failed drug-name mappings". **Decided on #26:** join acid / salt, prodrug ester ↔ active drug and racemate ↔ enantiomer; never epimers or diastereomers sold as separate drugs.

## 8. Most prevalent unmapped drugs (top 30)

Weighted prevalence is the share of the weighted population (all 9,254 participants) that reported the unit. All 30 are `not_in_model`, and none has an exact name-check hit. The full list (216 units) is in `unmapped_units_J.csv`. [derived]

| # | Ingredient | Participants | Weighted prev. | Loose name hit |
|---|---|---|---|---|
| 1 | ethinyl estradiol | 138 | 2.32% | |
| 2 | montelukast | 149 | 1.80% | |
| 3 | escitalopram | 84 | 1.65% | |
| 4 | hydrocodone | 139 | 1.65% | |
| 5 | rosuvastatin | 103 | 1.31% | yes (section 7) |
| 6 | insulin glargine | 118 | 1.13% | |
| 7 | allopurinol | 108 | 1.12% | |
| 8 | esomeprazole | 58 | 1.05% | |
| 9 | sitagliptin | 67 | 0.82% | |
| 10 | insulin aspart, human | 52 | 0.65% | |
| 11 | dextroamphetamine | 40 | 0.60% | |
| 12 | thyroid (USP) | 18 | 0.53% | |
| 13 | norgestimate | 45 | 0.53% | |
| 14 | sulfamethoxazole | 35 | 0.49% | |
| 15 | insulin lispro | 48 | 0.48% | |
| 16 | insulin detemir | 41 | 0.46% | |
| 17 | apixaban | 50 | 0.40% | |
| 18 | rivaroxaban | 34 | 0.39% | |
| 19 | clavulanate | 31 | 0.39% | |
| 20 | lisdexamfetamine | 22 | 0.38% | |
| 21 | vilanterol | 25 | 0.38% | |
| 22 | butalbital | 20 | 0.37% | |
| 23 | codeine | 28 | 0.36% | |
| 24 | triamterene | 43 | 0.32% | |
| 25 | benzonatate | 15 | 0.30% | |
| 26 | colchicine | 14 | 0.26% | |
| 27 | insulin isophane | 36 | 0.25% | |
| 28 | omega-3 fatty acids | 20 | 0.25% | |
| 29 | pseudoephedrine | 14 | 0.25% | |
| 30 | drospirenone | 11 | 0.24% | |

- The most prevalent **mapping failure** is `FLUORIDE TOPICAL` (d04397), at rank 40 with 0.16%. [derived]
- The 134 uncoded rows come from 130 participants, 1.29% of the weighted population. [derived]

## 9. Carried-over checks from #11

**(a) `RXDDRGID` across cycles** (`cycle_checks.py`). Code → name tables were built from `RXQ_RX_I` (2015–16) and `RXQ_RX_J` (2017–18). Results [derived]:
- 682 codes in I, 666 in J; 547 in both, 135 only in I, 119 only in J.
- No code has more than one name within a cycle.
- **No shared code changes its name** between I and J (0 of 547).
- **No name changes its code** (0 of 547 names found in both).
- Prefix counts: I `d` 655, `c` 21, `h` 5, `a` 1; J `d` 649, `c` 13, `h` 3, `a` 1.
- This covers only these two adjacent cycles (Multum Dec 2016 vs Dec 2018). Older cycles were not checked. [unverified: other cycles]

**(b) `SDMVSTRA` across cycles.** DEMO_I strata are 119–133 and DEMO_J strata are 134–148: **no overlap**. Both cycles have 15 strata × 2 PSUs (`SDMVPSU` 1, 2) = 30 cells. `SEQN` does not overlap either (I 83,732–93,702; J 93,703–102,956). Because the strata do not overlap, stacking the two cycles keeps strata distinct. PSU numbers (1, 2) repeat in both cycles. [derived]

## 10. Blockers and what is left

| What | Status |
|---|---|
| Hand check of the 34 unmatched Decagon drugs (28 with a rejected `approximateTerm` candidate, listed in section 4) | **not done**. Rules set on #26 (Decagon-side override file, half-day time box); done when the mapping pipeline is built. |
| Hand check of the 17 Decagon drugs with more than one `IN` | **not done**. All `IN`s are kept in the numbers; effect bounded at 231 mentions (section 7). Same as above: #26 alias rule. |
| Hand check of the loose name-check hits | **done** (section 7.1): 7 same or likely the same, rosuvastatin not the same, lactulose unclear. |
| `tenofovir` vs `tenofovir disoproxil`, `risedronate` vs `risedronic acid` | **done** (section 7.1): the same drug under different RxNorm `IN`s. |
| Variance / standard errors for the weighted shares (Taylor linearization with `SDMVSTRA` / `SDMVPSU`) | **not computed** (not asked by #16). Only weighted point estimates are given. |
| `RXDDRGID` stability for cycles other than I and J | not checked; relevant to #15 (which cycles to pool) |
| API failures | none left: 3 transient PubChem HTTP 502s, all retried successfully (`api_failures.log`). PubChem 404s ("no parent", "no UNII heading") are results, not failures. |

## Files written (all git-ignored except this note)

- `local/data/raw/nhanes/`: `DEMO_J.xpt`, `DEMO_I.xpt`, `RXQ_RX_I.xpt`
- `local/data/cache/`: `pubchem_decagon.json`, `rxnav_decagon.json`, `rxnav_nhanes_2017_pilot.json` (new NHANES lookups; `rxnav_nhanes_2017.json` was read and not changed)
- `local/data/scratch/coverage-pilot/`:
  - scripts: `cache.py`, `decagon_map.py`, `nhanes_map.py`, `coverage.py`, `extra_checks.py`, `cycle_checks.py`
  - tables: `decagon_rxnorm.csv` (645 rows), `decagon_names.json`, `nhanes_codes_J.csv` (666 rows), `unmapped_units_J.csv`
  - outputs: `coverage.out`, `extra_checks.out`, `cycle_checks.out`, `stereo_only.out`, `decagon_map.out`, `api_failures.log`
