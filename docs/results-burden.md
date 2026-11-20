# Burden results

Results of the burden pipeline of `docs/spec-burden.md` (Track 2, map #51): exposure (RQ2), burden rankings and how much of the burden is documented (RQ3), patient groups (RQ4), validation, and the limitations of §12.

All population figures are a share of the whole pooled population (all ages, non-users included), per 100,000 people (§1). The pool's weights sum to 320.1 million people. The exception is §9's band prevalence and band burden, which are per 100,000 band members.

**Scores (§2).** RQ3, RQ4 and the validation use the final score table of #40: hand-off `additive_free_w-none_lr0.01` (additive free encoder, MLP decoder, 10 cross-fitting folds per regime, Platt calibration per side effect), **`meta.json` at commit `a5e3077`**, dated 2026-11-13. The calibration report is in `docs/results-model.md`, "Final scores and calibration". The burden code reads only the score table; it never refits or rescores (§2).

## RQ2: exposure

How many people in the US take each drug pair, and how much of their co-medication the 645 Decagon drugs cover (§6, §7; #56). No scores are needed.

- **Data:** NHANES 2015–2016 (`I`) + 2017–March 2020 pre-pandemic (`P_`), interview weights `× 2/5.2` and `× 3.2/5.2` (§3). 25,531 participants; 10,873 reported a prescription drug in the past 30 days, 7,875 at least two. The 23 who refused or didn't know (`RXDUSE` 7 / 9) stay in every denominator with no drugs.
- **Mapping:** §4, with the hand fixes of `data/decagon_overrides.csv` (#55). 642 of the 645 Decagon drugs reach an RxNorm ingredient.
- **Run:** `python -m sief_druco.burden.rq2`, 2026-11-11. The tables, with the prevalence of every pair, are in `local/data/burden/rq2/I_P/`.

### Coverage

Share of each row by status, in %. Weighted rows use the pooled weights (§6). Point estimates; coverage is descriptive.

| Row | Total | Covered | Not in the model | Mapping failed |
|---|---|---|---|---|
| Distinct ingredients | 725 | 58.3 | 34.5 | 7.2 |
| Distinct ingredients, weighted by their users | 576.2 M | 86.5 | 13.1 | 0.4 |
| Mentions | 41,023 | 86.3 | 12.5 | 1.2 |
| Mentions, weighted | 583.4 M | 85.9 | 13.0 | 1.0 |
| Participants with ≥ 2 drugs | 7,875 | 58.0 | 39.5 | 2.5 |
| Participants with ≥ 2 drugs, weighted | 115.0 M | 57.1 | 40.7 | 2.3 |
| Person-pairs | 107,182 | 75.4 | 23.9 | 0.7 |
| Person-pairs, weighted | 1,480.7 M | 75.1 | 24.3 | 0.7 |

- **Units:** an ingredient is an RxNorm `IN`, or a source code that reached none; a combination product counts once per ingredient (a mention is a drug record × ingredient). The 339 records with no drug code (unknown, refused, don't know) are mentions and drugs of their participant, but no identifiable ingredient, so the first two rows leave them out.
- **Participants:** *covered* if all their drugs are in the 645, *mapping failed* if any drug failed, else *not in the model*. Two aliases of one Decagon drug count as two drugs here, though they never form a pair (§5).
- Three quarters of the weighted co-prescribed person-pairs can be scored. What limits coverage is drugs outside the 645, not the mapping: failures are about 1% of mentions and person-pairs.

**Supplement: mapping failed by side and cause.**

| Side | Cause | Drugs | Mentions | % of mentions | % weighted |
|---|---|---|---|---|---|
| Source | no drug code (55555 / 77777 / 99999) | 339 records | 339 | 0.83 | 0.64 |
| Source | `a` / `c` / `h` code | 35 | 91 | 0.22 | 0.20 |
| Source | retired RxNorm concept | 17 | 55 | 0.13 | 0.18 |
| Decagon | Decagon drug reached no `IN` | 0 | 0 | 0 | 0 |

- **Decagon side:** 3 Decagon drugs reach no `IN` (#55). No source ingredient matches their names exactly, so no source mention fails on the Decagon side.
- **What the hand fixes and aliases moved:** against the same run without the overrides, 1,251 of 41,023 mentions (3.0%) moved into *covered* (1,226 from *not in the model*, led by montelukast 457, rosuvastatin 320 and allopurinol 228; 25 from Decagon-side *mapping failed*: sevelamer, acarbose, ceftazidime), and 28 moved out (betamethasone and isotretinoin, dropped as provisional joins). 153 covered mentions join through an alias `IN`.
- **Footnotes:** NHANES "tenofovir" joins Decagon's tenofovir disoproxil, because NHANES doesn't separate TDF from TAF (§4). One `IN` is claimed by two Decagon drugs: hydrocortisone joins `CID000003640` (396 labelled pairs), not the butyrate ester `CID000003642` (26).

### Pair prevalence

Prevalence of a pair = pooled weight of the persons taking both ÷ the pooled weight of everyone, with a Korn–Graubard 95% interval from the Taylor variance (§7). Pairs holding a record with no drug code have no prevalence.

**How many covered pairs pass the minimum count.** A pair is ranked only with ≥ 30 unweighted participants; every covered pair still counts in coverage and in later totals.

| Participants | Covered pairs | Share of covered pair weight (%) |
|---|---|---|
| any | 15,296 | 100 |
| ≥ 10 | 1,644 | 60.0 |
| ≥ 20 | 713 | 44.6 |
| ≥ 30 (ranked) | **429** | 36.0 |

- **Ranked pairs, 429 against the rough estimate of ~500 (§7):** the estimate scaled 2017–18 alone up by the pool's 2.75× participants. On that scale, ≥ 30 in the pool means ≥ 11 in 2017–18, where 513 covered pairs pass. Only 398 of them reach 30 in the pool, and 31 others come up from below 11. Pairs just above a count cut in one cycle are there partly by chance, and 2015–16 and the 2019–20 part of `P_` don't add as many of them as straight scaling assumes.
- The ranked 429 pairs hold 36% of the covered pair weight. The other 64% sits in 14,867 pairs that are each too rare to rank; they still enter the side-effect totals and burden shares of §8 and §10.

**The 20 most prevalent pairs among all co-prescribed pairs**, per 100,000 people, with each cycle's own share (a survey domain, as a share of that cycle's population; §3 supplement).

| # | Pair | Status | Participants | Prevalence (95% CI) | 2015–16 | `P_` |
|---|---|---|---|---|---|---|
| 1 | lisinopril + metformin | covered | 554 | 2,224 (1,948–2,527) | 2,195 | 2,242 |
| 2 | hydrochlorothiazide + lisinopril | covered | 437 | 1,814 (1,571–2,083) | 1,814 | 1,814 |
| 3 | atorvastatin + lisinopril | covered | 421 | 1,623 (1,354–1,928) | 1,194 | 1,886 |
| 4 | atorvastatin + metoprolol | covered | 305 | 1,510 (1,208–1,865) | 1,227 | 1,685 |
| 5 | atorvastatin + metformin | covered | 409 | 1,499 (1,287–1,735) | 974 | 1,820 |
| 6 | hydrochlorothiazide + losartan | covered | 334 | 1,419 (1,180–1,690) | 1,532 | 1,349 |
| 7 | lisinopril + metoprolol | covered | 269 | 1,257 (968–1,606) | 1,393 | 1,174 |
| 8 | hydrochlorothiazide + metformin | covered | 314 | 1,096 (925–1,289) | 1,116 | 1,083 |
| 9 | metformin + simvastatin | covered | 249 | 1,079 (867–1,327) | 1,273 | 960 |
| 10 | metformin + metoprolol | covered | 248 | 1,072 (901–1,266) | 850 | 1,208 |
| 11 | lisinopril + simvastatin | covered | 239 | 1,069 (882–1,283) | 1,563 | 766 |
| 12 | amlodipine + hydrochlorothiazide | covered | 309 | 1,069 (850–1,326) | 1,134 | 1,029 |
| 13 | acetaminophen + hydrocodone | not in the model | 255 | 1,028 (822–1,269) | 1,253 | 890 |
| 14 | amlodipine + lisinopril | covered | 290 | 1,021 (841–1,227) | 1,123 | 958 |
| 15 | amlodipine + atorvastatin | covered | 299 | 1,012 (822–1,232) | 723 | 1,189 |
| 16 | levothyroxine + lisinopril | covered | 186 | 964 (752–1,216) | 1,281 | 769 |
| 17 | amlodipine + metformin | covered | 270 | 937 (756–1,148) | 705 | 1,079 |
| 18 | amlodipine + metoprolol | covered | 235 | 918 (745–1,118) | 963 | 890 |
| 19 | atorvastatin + hydrochlorothiazide | covered | 274 | 909 (673–1,200) | 745 | 1,009 |
| 20 | albuterol + fluticasone | covered | 229 | 893 (715–1,102) | 924 | 875 |

- Nineteen of the 20 are covered and ranked; the 20th ranked pair is losartan + metformin (887 per 100,000). The one uncovered pair is acetaminophen + hydrocodone, since hydrocodone is not in the model (#55).
- The most common pair, lisinopril + metformin, is taken by about 2.2% of the population: ≈ 7.1 million US residents. The list is the cardiometabolic core of US prescribing: ACE inhibitors and ARBs, thiazides, statins, metformin and beta blockers.
- **Cycle to cycle:** the pooled estimate is the result (§3). The two cycles agree for most pairs. The biggest shifts are consistent with the move from simvastatin to atorvastatin: atorvastatin pairs are higher in `P_` (atorvastatin + metformin 974 → 1,820), and simvastatin pairs are lower (lisinopril + simvastatin 1,563 → 766).

### Development check on 2017–2018 alone

Per §3's order of work, the pipeline was built on NHANES 2017–2018 (`J`) and checked against the coverage pilot (#16, `research/nhanes-coverage-pilot.md`) before the switch to the main pool. `tests/test_burden_rq2.py` holds these numbers.

| 2017–18 | Covered pairs | ≥ 10 participants | ≥ 30 participants |
|---|---|---|---|
| Pilot (#16; spec §7) | 8,884 | 542 (41% of covered pair weight) | 114 (21%) |
| This pipeline, without the hand fixes | 8,856 | 542 (41.0%) | 114 (20.5%) |
| This pipeline, with the hand fixes | 9,609 | 585 (40.5%) | 118 (19.5%) |

- **Without the hand fixes** the pipeline reproduces the pilot. The only difference: the pilot counted an alias `IN` and its main `IN` (methylphenidate and dexmethylphenidate) as two drugs, so a third drug paired with each gave two pairs. De-duplicating covered drugs by STITCH ID (§5) merges 28 such pairs. The covered column of the coverage table matches the pilot's too: 15,560 ingredient mentions (84.0% covered), 2,941 participants with ≥ 2 drugs (50.7% covered, 47.6% weighted) and 69.2% of weighted person-pairs. The split of the rest between *not in the model* and *mapping failed* differs slightly, since the Decagon-side name check differs: the pilot's flagged risedronate, a matched Decagon drug, while this one adds SIDER names and flags only Decagon drugs that reached no `IN` (sevelamer, acarbose and ceftazidime here).
- **With the hand fixes,** 753 more covered pairs, 43 more pass 10 participants and 4 more pass 30. They are the pairs of the drugs #55 joined, mostly montelukast, rosuvastatin and allopurinol. Participants with all drugs covered rise from 50.7% to 57.9% and weighted person-pairs covered from 69.2% to 74.5%, so the main pool's coverage table can't be compared with the pilot's directly.

## RQ3: burden

Burden weights each covered pair's prevalence by the model's calibrated `p_r(a, b)` (§8): a pair's burden is prevalence × Σ_r `p_r`, the expected number of side-effect signals among the people who take it, per 100,000 people. It measures exposure to a signal (reported together in TWOSIDES), not cases of harm. All 963 side effects count the same, and labelled pairs keep their cross-fitted prediction (S regime) for every side effect.

- **Data:** the main NHANES pool and mapping of RQ2, with the scores of #40.
- **Run:** `python -m sief_druco.burden.rq3 local/crossfit/additive_free_w-none_lr0.01`, 2026-11-20. Tables in `local/data/burden/rq3/additive_free_w-none_lr0.01/I_P/`.
- **Total:** the 15,296 covered pairs carry a burden of 46.96 million per 100,000 people: about 470 expected signals per US resident, summed over every covered pair a person takes. The 429 ranked pairs (≥ 30 participants) hold 43.5% of it.
- **Uncertainty:** intervals cover prevalence only, with `p` held fixed (§8). Model uncertainty shows only in the side-effect columns (held-out AUPRC, predicted vs observed positives).

### Pairs

The 20 pairs with the highest burden, of the 429 ranked. DDInter: the level DDInter lists for the pair (§10). *Unknown* is a pair DDInter lists without a grade, which counts as not documented; *not listed* is a pair DDInter doesn't list.

| # | Pair | Participants | Prevalence (95% CI) | Σ_r p_r | Burden (95% CI) | Known to Decagon | DDInter |
|---|---|---|---|---|---|---|---|
| 1 | lisinopril + metformin | 554 | 2,224 (1,948–2,527) | 174.8 | 388,745 (340,579–441,690) | yes | Moderate |
| 2 | hydrochlorothiazide + lisinopril | 437 | 1,814 (1,571–2,083) | 179.0 | 324,732 (281,255–372,927) | yes | Moderate |
| 3 | atorvastatin + lisinopril | 421 | 1,623 (1,354–1,928) | 168.1 | 272,726 (227,584–324,079) | yes | Unknown |
| 4 | atorvastatin + metoprolol | 305 | 1,510 (1,208–1,865) | 167.5 | 253,044 (202,316–312,467) | yes | Unknown |
| 5 | atorvastatin + metformin | 409 | 1,499 (1,287–1,735) | 164.2 | 246,087 (211,325–284,865) | yes | Unknown |
| 6 | hydrochlorothiazide + losartan | 334 | 1,419 (1,180–1,690) | 171.3 | 243,003 (202,207–289,517) | yes | Unknown |
| 7 | lisinopril + metoprolol | 269 | 1,257 (968–1,606) | 180.3 | 226,702 (174,493–289,466) | yes | Unknown |
| 8 | lisinopril + simvastatin | 239 | 1,069 (882–1,283) | 180.0 | 192,487 (158,869–231,053) | yes | Unknown |
| 9 | metformin + metoprolol | 248 | 1,072 (901–1,266) | 177.0 | 189,693 (159,381–224,041) | yes | Unknown |
| 10 | hydrochlorothiazide + metformin | 314 | 1,096 (925–1,289) | 173.1 | 189,638 (160,066–223,040) | yes | Moderate |
| 11 | metformin + simvastatin | 249 | 1,079 (867–1,327) | 175.1 | 188,853 (151,728–232,231) | yes | Unknown |
| 12 | levothyroxine + lisinopril | 186 | 964 (752–1,216) | 182.5 | 175,867 (137,192–221,973) | yes | Unknown |
| 13 | amlodipine + hydrochlorothiazide | 309 | 1,069 (850–1,326) | 164.0 | 175,233 (139,326–217,486) | yes | Minor |
| 14 | amlodipine + lisinopril | 290 | 1,021 (841–1,227) | 166.0 | 169,406 (139,592–203,648) | yes | Minor |
| 15 | albuterol + fluticasone | 229 | 893 (715–1,102) | 179.0 | 159,901 (128,012–197,264) | yes | Minor |
| 16 | amlodipine + atorvastatin | 299 | 1,012 (822–1,232) | 155.7 | 157,571 (128,019–191,842) | yes | Moderate |
| 17 | levothyroxine + metformin | 170 | 870 (672–1,107) | 178.8 | 155,479 (120,126–197,905) | yes | Moderate |
| 18 | amlodipine + metoprolol | 235 | 918 (745–1,118) | 166.6 | 152,892 (124,133–186,265) | yes | Moderate |
| 19 | atorvastatin + hydrochlorothiazide | 274 | 909 (673–1,200) | 166.2 | 151,025 (111,877–199,332) | yes | Unknown |
| 20 | amlodipine + metformin | 270 | 937 (756–1,148) | 160.7 | 150,520 (121,407–184,455) | yes | Unknown |

- **The top pair,** lisinopril + metformin, is taken by ≈ 7.1 million US residents, each carrying about 175 expected signals across the 963 side effects; hydrochlorothiazide + lisinopril follows with ≈ 5.8 million.
- **The burden ranking is mostly the prevalence ranking.** Among the 429 ranked pairs, prevalence spans 36× (62 to 2,224 per 100,000), but Σ_r `p_r` stays within a narrow band (median 167, interquartile range 149–177). The burden and prevalence ranks correlate at Spearman 0.95, and 48 of the top 50 by burden are also in the top 50 by prevalence. The model separates labelled from unlabelled pairs (mean Σ_r `p_r` 108 against 12 over all covered pairs), but 428 of the 429 ranked pairs are labelled: pairs common enough to rank were reported together in FAERS. Among labelled pairs the calibrated scores sit close together, so prevalence decides the order.
- **The one unlabelled ranked pair** is hydralazine + losartan (rank 429 of 429, Σ_r `p_r` 8.8); see the novelty candidates below.

### Triples

Triples (pair × side effect) of the 429 ranked pairs: 413,127, ranked by burden_r = prevalence × `p_r`. "Triple known to Decagon" means that exact triple is a training label. DDInter is pair-level only (§10).

| # | Pair | Side effect | p_r | Burden | Triple known to Decagon | Pair in DDInter |
|---|---|---|---|---|---|---|
| 1 | lisinopril + metformin | Fatigue | 0.622 | 1,384 | no | Moderate |
| 2 | lisinopril + metformin | Back Ache | 0.618 | 1,374 | no | Moderate |
| 3 | lisinopril + metformin | asthenia | 0.615 | 1,367 | no | Moderate |
| 4 | lisinopril + metformin | pleural pain | 0.612 | 1,361 | no | Moderate |
| 5 | atorvastatin + lisinopril | arterial pressure NOS decreased | 0.838 | 1,360 | yes | Unknown |
| 6 | lisinopril + metformin | anaemia | 0.609 | 1,355 | no | Moderate |
| 7 | lisinopril + metformin | Aching joints | 0.583 | 1,297 | no | Moderate |
| 8 | lisinopril + metformin | Hypoventilation | 0.572 | 1,273 | no | Moderate |
| 9 | hydrochlorothiazide + lisinopril | Difficulty breathing | 0.693 | 1,256 | no | Moderate |
| 10 | hydrochlorothiazide + lisinopril | neumonia | 0.687 | 1,247 | no | Moderate |
| 11 | lisinopril + metformin | chest pain | 0.560 | 1,245 | no | Moderate |
| 12 | lisinopril + metformin | apoplexy | 0.549 | 1,221 | yes | Moderate |
| 13 | atorvastatin + lisinopril | Pain | 0.750 | 1,218 | yes | Unknown |
| 14 | lisinopril + metformin | diarrhea | 0.542 | 1,206 | no | Moderate |
| 15 | lisinopril + metformin | Difficulty breathing | 0.541 | 1,203 | yes | Moderate |
| 16 | lisinopril + metformin | Infection Urinary Tract | 0.540 | 1,202 | yes | Moderate |
| 17 | atorvastatin + metformin | Pain | 0.791 | 1,186 | yes | Unknown |
| 18 | lisinopril + metformin | Extremity pain | 0.529 | 1,177 | no | Moderate |
| 19 | lisinopril + metformin | Mod | 0.518 | 1,152 | no | Moderate |
| 20 | hydrochlorothiazide + lisinopril | Fatigue | 0.634 | 1,150 | yes | Moderate |

- The most prevalent pair holds 14 of the top 20 triples: with `p_r` between 0.5 and 0.85 for its common side effects, prevalence again decides.
- Side-effect names are Decagon's own (`bio-decagon-combo.csv`), truncations included ("neumonia", "Mod").
- Most top triples are predictions for a labelled pair's *unlabelled* side effects: the pair is in TWOSIDES, but not with that side effect.

### Side effects

Each side effect's burden summed over **all** 15,296 covered pairs, with a Taylor interval (`p` fixed), next to the model's held-out quality for it (`side_effects.csv`, #40).

| # | Side effect | Burden (95% CI) | Labelled pairs | AUPRC U / S | Predicted / observed, U | Predicted / observed, S |
|---|---|---|---|---|---|---|
| 1 | arterial pressure NOS decreased | 188,733 (174,843–202,623) | 28,568 | 0.60 / 0.63 | 28,579 / 28,568 | 28,418 / 28,568 |
| 2 | Pain | 188,244 (174,115–202,374) | 23,894 | 0.62 / 0.64 | 23,905 / 23,894 | 23,932 / 23,894 |
| 3 | neumonia | 183,886 (169,986–197,786) | 24,430 | 0.60 / 0.61 | 24,439 / 24,430 | 24,432 / 24,430 |
| 4 | anaemia | 177,247 (164,109–190,386) | 27,006 | 0.58 / 0.61 | 27,006 / 27,006 | 27,055 / 27,006 |
| 5 | Fatigue | 176,652 (163,679–189,625) | 24,260 | 0.57 / 0.59 | 24,266 / 24,260 | 24,130 / 24,260 |
| 6 | Difficulty breathing | 174,671 (161,799–187,543) | 26,037 | 0.57 / 0.59 | 26,044 / 26,037 | 25,980 / 26,037 |
| 7 | asthenia | 173,918 (161,142–186,695) | 23,515 | 0.56 / 0.58 | 23,509 / 23,515 | 23,406 / 23,515 |
| 8 | chest pain | 170,160 (157,561–182,759) | 21,013 | 0.58 / 0.60 | 21,013 / 21,013 | 21,078 / 21,013 |
| 9 | Hypoventilation | 169,970 (157,301–182,640) | 21,322 | 0.60 / 0.62 | 21,321 / 21,322 | 21,291 / 21,322 |
| 10 | nausea | 166,101 (153,730–178,472) | 25,190 | 0.56 / 0.58 | 25,193 / 25,190 | 25,158 / 25,190 |
| 11 | edema extremities | 165,457 (153,070–177,844) | 21,981 | 0.58 / 0.59 | 21,973 / 21,981 | 21,943 / 21,981 |
| 12 | Back Ache | 165,111 (152,849–177,373) | 19,930 | 0.60 / 0.61 | 19,939 / 19,930 | 19,988 / 19,930 |
| 13 | diarrhea | 164,567 (152,440–176,695) | 23,848 | 0.53 / 0.56 | 23,847 / 23,848 | 24,067 / 23,848 |
| 14 | pleural pain | 162,125 (150,160–174,089) | 21,781 | 0.57 / 0.58 | 21,781 / 21,781 | 21,764 / 21,781 |
| 15 | emesis | 157,813 (146,270–169,355) | 23,043 | 0.50 / 0.53 | 23,033 / 23,043 | 23,071 / 23,043 |
| 16 | abdominal pain | 156,410 (144,791–168,029) | 21,410 | 0.55 / 0.56 | 21,410 / 21,410 | 21,298 / 21,410 |
| 17 | Anxiety | 152,089 (140,688–163,489) | 18,980 | 0.55 / 0.57 | 18,976 / 18,980 | 18,986 / 18,980 |
| 18 | Aching joints | 150,577 (139,342–161,812) | 18,061 | 0.56 / 0.57 | 18,057 / 18,061 | 18,032 / 18,061 |
| 19 | weight gain | 149,953 (138,868–161,038) | 17,112 | 0.55 / 0.57 | 17,106 / 17,112 | 17,119 / 17,112 |
| 20 | dizziness | 148,198 (137,350–159,046) | 20,204 | 0.50 / 0.52 | 20,204 / 20,204 | 20,057 / 20,204 |

- **The order follows the labels.** The side-effect ranking correlates with each side effect's number of labelled pairs at Spearman 0.99: calibrated scores reproduce how often a side effect is reported, and the most reported ones are broad FAERS terms (pain, fatigue, low blood pressure). Predicted and observed positives agree within 1% for these side effects in both regimes.
- **Burden is spread thin:** the top 10 side effects hold 3.8% of the summed side-effect burden, the top 100 hold 27.6%. The last-ranked, impetigo, still has 2,619 per 100,000.

### How much is undocumented (§10)

Every covered pair and triple has two labels: **known to Decagon** (pair: any training label; triple: that triple is a label) and **clinically documented** (DDInter 1.0 lists the pair as Major, Moderate or Minor; triples inherit the pair's status). Shares of total burden over all covered pairs and triples, with no count floor, in %:

| Pairs | DDInter documented | Not documented | Not checkable | Total |
|---|---|---|---|---|
| Known to Decagon | 33.16 | 65.97 | 0.51 | 99.64 |
| Not known to Decagon | 0.08 | **0.27** | 0.00 | 0.35 |

| Triples | DDInter documented (pair) | Not documented | Not checkable | Total |
|---|---|---|---|---|
| Known to Decagon | 10.20 | 19.64 | 0.15 | 29.99 |
| Not known to Decagon | 23.03 | **46.61** | 0.37 | 70.01 |

- **Pairs:** 99.6% of the pair burden is in pairs Decagon already knows, so the burden of pairs that are new to both sources (the bold cell) is 0.27%. The model scores unlabelled pairs roughly ten times lower than labelled ones (mean Σ_r `p_r` 12 against 108), and the unlabelled pairs are rarer.
- **Triples:** 70% of the triple burden is in side effects a pair has no label for, and 46.6% is new to both sources. The DDInter half of that is pair-level only: a documented pair says nothing about which side effect it causes.
- **Footnote, DDInter "Unknown":** 64.1% of the pair burden (61.4% of covered pair weight) is in pairs DDInter lists without a grade. They count as not documented, so almost all of the 66% "not documented" pair burden is DDInter-listed but ungraded. Checkable pairs DDInter doesn't list hold 2.2% of the burden (2.7% with the not-checkable pairs).
- **Checkable:** DDInter can check 99.5% of the pair burden (99.3% of covered pair weight). 1,896 of DDInter's 1,971 drug names reach an `IN`, and 630 of the 645 Decagon drugs appear in DDInter.

**Top 50 against the base rate.** The 2×2 for the top 50 ranked pairs and top 50 ranked triples, as counts, next to the share of all covered pairs and triples in each cell (a count share, not a burden share).

| Level | Known to Decagon / clinically documented | Top 50 | Top 50 (%) | All covered (%) |
|---|---|---|---|---|
| Pairs | known / documented | 20 | 40 | 26.3 |
| Pairs | known / not documented | 30 | 60 | 59.6 |
| Pairs | known / not checkable | 0 | 0 | 1.6 |
| Pairs | not known / documented | 0 | 0 | 2.4 |
| Pairs | not known / not documented | 0 | 0 | 9.8 |
| Pairs | not known / not checkable | 0 | 0 | 0.3 |
| Triples | known / documented | 11 | 22 | 3.7 |
| Triples | known / not documented | 6 | 12 | 7.3 |
| Triples | known / not checkable | 0 | 0 | 0.2 |
| Triples | not known / documented | 32 | 64 | 25.0 |
| Triples | not known / not documented | 1 | 2 | 62.1 |
| Triples | not known / not checkable | 0 | 0 | 1.7 |

- All top 50 pairs are known to Decagon, against 87.5% of covered pairs. DDInter grades 20 of them (14 Moderate, 5 Minor, 1 Major: amlodipine + simvastatin); the other 30 are DDInter "Unknown".
- Among the top 50 triples, the most common pattern (32) is a documented pair with a side effect Decagon has no label for.

### Novelty candidates

Pairs neither known to Decagon nor clinically documented, with ≥ 10 participants, by pair burden. Under 30 participants is flagged as imprecise.

| # | Pair | Participants | Prevalence (95% CI) | Σ_r p_r | Burden | Imprecise (< 30) |
|---|---|---|---|---|---|---|
| 1 | hydralazine + losartan | 37 | 87.7 (53.8–135.1) | 8.8 | 773 | no |
| 2 | hydralazine + metformin | 28 | 68.9 (37.3–116.4) | 10.7 | 740 | yes |
| 3 | hydralazine + hydrochlorothiazide | 21 | 46.0 (25.3–76.9) | 11.0 | 507 | yes |
| 4 | gabapentin + hydralazine | 17 | 30.7 (15.0–55.7) | 13.8 | 423 | yes |
| 5 | hydralazine + potassium chloride | 14 | 26.2 (12.0–49.5) | 12.1 | 317 | yes |
| 6 | hydralazine + allopurinol | 11 | 24.5 (11.7–45.1) | 9.3 | 228 | yes |
| 7 | glipizide + hydralazine | 13 | 21.0 (12.1–34.0) | 5.5 | 116 | yes |
| 8 | hydralazine + isosorbide dinitrate | 10 | 20.4 (9.5–38.3) | 2.1 | 44 | yes |

- **8 candidates, all with hydralazine.** Hydralazine has only 38 labelled pairs in Decagon (rank 552 of 645 drugs by labelled pairs; median 167), so its pairs with common drugs are unlabelled. Even its labelled pairs score low (Σ_r `p_r` ≈ 35, against ≈ 170 for the top pairs). None of the 8 is listed in DDInter.
- The top candidate, hydralazine + losartan, is taken by ≈ 0.28 million US residents. Its burden is 773 per 100,000, against 388,745 for the top pair: the candidates are a short list of hypotheses, not a large hidden burden.
- Hydralazine + isosorbide dinitrate is a fixed-dose combination (BiDil) as well as a co-prescription. Combination pills count as co-prescribed (§5).

**DrugBank hand check (§10).** All 8 candidates (fewer than 20), looked up by hand on 2026-11-20. DrugBank's free interaction checker is gone, so each pair was looked up in the drug-interaction list of hydralazine's DrugBank page (signed in). Pairs not found there were checked again from the partner drug's page. Listed yes/no only; no text copied. Source: DrugBank 6.0 (Knox et al., *Nucleic Acids Research*, 2024), read as within DrugBank's terms for non-clinical research.

| # | Pair | Listed in DrugBank | Severity |
|---|---|---|---|
| 1 | hydralazine + losartan | yes | not shown |
| 2 | hydralazine + metformin | yes | not shown |
| 3 | hydralazine + hydrochlorothiazide | no | — |
| 4 | gabapentin + hydralazine | yes | not shown |
| 5 | hydralazine + potassium chloride | no | — |
| 6 | hydralazine + allopurinol | yes | not shown |
| 7 | glipizide + hydralazine | yes | not shown |
| 8 | hydralazine + isosorbide dinitrate | no | — |

- **5 of 8 are listed:** losartan, metformin, gabapentin, allopurinol and glipizide with hydralazine.
- **Severity: partial.** The spec records severity, but the signed-in drug-page view shows none, so the column can't be filled. The check is otherwise complete.
- **Footnote:** the partner pages of hydrochlorothiazide and isosorbide dinitrate list dihydralazine, a different drug, but not hydralazine.

## RQ4: patient groups

Burden by age band and by drug-count band on the main pool (§9). A person's burden is the sum over their covered pairs of Σ_r `p_r`; a drug counts again in each pair it is part of, so this is not a patient risk. Bands are survey domains. The drug count is a person's distinct ingredients, including those outside the 645, plus one per drug that reached no ingredient.

- **Run:** `python -m sief_druco.burden.rq4 local/crossfit/additive_free_w-none_lr0.01`, 2026-11-20. Tables in `local/data/burden/rq4/additive_free_w-none_lr0.01/`.

| Band | Participants | Population (%) | Share of burden (%) | Burden per 100,000 band members (95% CI) | Covered pairs per person | Rankable pairs |
|---|---|---|---|---|---|---|
| Age 0–17 | 9,846 | 23.4 | 1.3 | 2.6 M (2.0–3.3 M) | 0.23 | 5 |
| Age 18–44 | 6,693 | 35.1 | 9.9 | 13.2 M (11.0–15.3 M) | 1.05 | 2 |
| Age 45–64 | 5,203 | 26.0 | 41.4 | 74.8 M (67.2–82.4 M) | 5.46 | 101 |
| Age 65+ | 3,789 | 15.5 | 47.4 | 143.8 M (132.2–155.5 M) | 10.55 | 209 |
| Drugs 0–1 | 17,656 | 64.1 | 0 | — | — | — |
| Drugs 2 | 2,075 | 9.6 | 1.9 | 9.4 M (8.9–9.9 M) | 0.67 | 2 |
| Drugs 3–4 | 2,489 | 11.8 | 11.5 | 45.5 M (43.2–47.7 M) | 3.27 | 32 |
| Drugs 5–9 | 2,652 | 11.8 | 49.3 | 195.9 M (188.0–203.7 M) | 14.32 | 166 |
| Drugs 10+ | 659 | 2.7 | 37.3 | 660.0 M (620.1–699.9 M) | 50.08 | 111 |

- **Age:** people aged 65+, 15.5% of the population (≈ 49.6 million US residents), hold 47.4% of the burden; with the 45–64 band, people aged 45+ hold 88.8%. A person aged 65+ carries on average ≈ 1,440 expected signals, one aged 0–17 ≈ 26.
- **Drug count:** people on 5 or more drugs, 14.5% of the population (≈ 46.4 million), hold 86.6% of the burden. People on 10 or more, 2.7% (≈ 8.5 million), hold 37.3%.
- **The rise is the pair count.** Burden per person divided by covered pairs per person stays between 114 and 140 in every band, so almost all of the rise with age and drug count comes from the number of pairs (n drugs give n(n−1)/2 pairs). The burden per pair barely changes.
- "Rankable pairs" counts pairs with ≥ 30 participants within the band. The under-45 bands and the 2-drug band have fewer than 10 (§9 expected a handful).

**Top pairs per band** by band burden = band prevalence × Σ_r `p_r`, per 100,000 band members, only pairs with ≥ 30 participants in the band. Overall rank is the rank in the whole pool (§8, "—" if it has fewer than 30 participants there).

**Age 0–17** (fewer than 10 rankable pairs)

| # | Pair | Participants | Band prevalence (95% CI) | Band burden | Overall rank | Known to Decagon | DDInter |
|---|---|---|---|---|---|---|---|
| 1 | albuterol + fluticasone | 80 | 950 (621–1,390) | 170,065 | 15 | yes | Minor |
| 2 | albuterol + montelukast | 68 | 815 (518–1,218) | 135,361 | 61 | yes | Unknown |
| 3 | fluticasone + montelukast | 51 | 578 (370–860) | 93,586 | 52 | yes | Unknown |
| 4 | cetirizine + fluticasone | 36 | 360 (214–569) | 61,724 | 208 | yes | Unknown |
| 5 | albuterol + cetirizine | 31 | 346 (189–581) | 60,723 | 308 | yes | Unknown |

**Age 18–44** (fewer than 10 rankable pairs)

| # | Pair | Participants | Band prevalence (95% CI) | Band burden | Overall rank | Known to Decagon | DDInter |
|---|---|---|---|---|---|---|---|
| 1 | hydrochlorothiazide + lisinopril | 36 | 460 (308–660) | 82,290 | 2 | yes | Moderate |
| 2 | lisinopril + metformin | 33 | 398 (222–658) | 69,547 | 1 | yes | Moderate |

**Age 45–64**

| # | Pair | Participants | Band prevalence (95% CI) | Band burden | Overall rank | Known to Decagon | DDInter |
|---|---|---|---|---|---|---|---|
| 1 | lisinopril + metformin | 279 | 4,514 (3,697–5,450) | 788,983 | 1 | yes | Moderate |
| 2 | hydrochlorothiazide + lisinopril | 220 | 3,717 (2,982–4,573) | 665,464 | 2 | yes | Moderate |
| 3 | atorvastatin + lisinopril | 184 | 3,035 (2,347–3,857) | 510,107 | 3 | yes | Unknown |
| 4 | atorvastatin + metformin | 187 | 2,742 (2,048–3,589) | 450,225 | 5 | yes | Unknown |
| 5 | lisinopril + simvastatin | 106 | 2,127 (1,613–2,751) | 382,942 | 8 | yes | Unknown |
| 6 | hydrochlorothiazide + losartan | 147 | 2,159 (1,648–2,775) | 369,780 | 6 | yes | Unknown |
| 7 | levothyroxine + lisinopril | 86 | 1,991 (1,344–2,838) | 363,395 | 12 | yes | Unknown |
| 8 | metformin + simvastatin | 100 | 2,071 (1,596–2,642) | 362,648 | 11 | yes | Unknown |
| 9 | hydrochlorothiazide + metformin | 137 | 1,949 (1,477–2,520) | 337,288 | 10 | yes | Moderate |
| 10 | lisinopril + metoprolol | 106 | 1,870 (1,364–2,500) | 337,208 | 7 | yes | Unknown |

**Age 65+**

| # | Pair | Participants | Band prevalence (95% CI) | Band burden | Overall rank | Known to Decagon | DDInter |
|---|---|---|---|---|---|---|---|
| 1 | lisinopril + metformin | 242 | 5,884 (4,624–7,364) | 1,028,560 | 1 | yes | Moderate |
| 2 | atorvastatin + metoprolol | 199 | 5,791 (4,831–6,876) | 970,084 | 4 | yes | Unknown |
| 3 | hydrochlorothiazide + losartan | 167 | 4,805 (3,730–6,080) | 823,073 | 6 | yes | Unknown |
| 4 | hydrochlorothiazide + lisinopril | 180 | 4,427 (3,508–5,505) | 792,590 | 2 | yes | Moderate |
| 5 | atorvastatin + lisinopril | 223 | 4,638 (3,753–5,661) | 779,473 | 3 | yes | Unknown |
| 6 | lisinopril + metoprolol | 150 | 4,196 (3,188–5,410) | 756,437 | 7 | yes | Unknown |
| 7 | atorvastatin + metformin | 200 | 4,569 (3,667–5,617) | 750,276 | 5 | yes | Unknown |
| 8 | metformin + metoprolol | 150 | 4,051 (2,968–5,385) | 716,856 | 9 | yes | Unknown |
| 9 | levothyroxine + atorvastatin | 119 | 3,942 (3,039–5,021) | 662,809 | 24 | yes | Unknown |
| 10 | losartan + metoprolol | 132 | 3,765 (2,692–5,107) | 655,554 | 31 | yes | Unknown |

**Drugs 2** (fewer than 10 rankable pairs)

| # | Pair | Participants | Band prevalence (95% CI) | Band burden | Overall rank | Known to Decagon | DDInter |
|---|---|---|---|---|---|---|---|
| 1 | hydrochlorothiazide + lisinopril | 61 | 3,606 (2,201–5,542) | 645,500 | 2 | yes | Moderate |
| 2 | albuterol + fluticasone | 35 | 870 (557–1,294) | 155,699 | 15 | yes | Minor |

**Drugs 3–4**

| # | Pair | Participants | Band prevalence (95% CI) | Band burden | Overall rank | Known to Decagon | DDInter |
|---|---|---|---|---|---|---|---|
| 1 | lisinopril + metformin | 155 | 4,969 (3,779–6,400) | 868,591 | 1 | yes | Moderate |
| 2 | hydrochlorothiazide + lisinopril | 117 | 4,085 (3,188–5,147) | 731,356 | 2 | yes | Moderate |
| 3 | hydrochlorothiazide + losartan | 91 | 3,574 (2,575–4,820) | 612,159 | 6 | yes | Unknown |
| 4 | atorvastatin + lisinopril | 108 | 3,373 (2,292–4,773) | 566,924 | 3 | yes | Unknown |
| 5 | lisinopril + simvastatin | 67 | 2,295 (1,752–2,949) | 413,158 | 8 | yes | Unknown |
| 6 | levothyroxine + lisinopril | 48 | 2,259 (1,543–3,186) | 412,217 | 12 | yes | Unknown |
| 7 | atorvastatin + metoprolol | 44 | 2,453 (1,552–3,676) | 411,018 | 4 | yes | Unknown |
| 8 | metformin + simvastatin | 62 | 2,228 (1,504–3,173) | 390,086 | 11 | yes | Unknown |
| 9 | atorvastatin + metformin | 89 | 2,238 (1,545–3,128) | 367,473 | 5 | yes | Unknown |
| 10 | amlodipine + hydrochlorothiazide | 73 | 2,178 (1,517–3,023) | 357,071 | 13 | yes | Minor |

**Drugs 5–9**

| # | Pair | Participants | Band prevalence (95% CI) | Band burden | Overall rank | Known to Decagon | DDInter |
|---|---|---|---|---|---|---|---|
| 1 | lisinopril + metformin | 289 | 10,443 (8,692–12,409) | 1,825,316 | 1 | yes | Moderate |
| 2 | atorvastatin + lisinopril | 226 | 7,989 (6,156–10,159) | 1,342,683 | 3 | yes | Unknown |
| 3 | atorvastatin + metoprolol | 182 | 7,998 (6,158–10,176) | 1,339,853 | 4 | yes | Unknown |
| 4 | lisinopril + metoprolol | 164 | 7,355 (5,458–9,656) | 1,326,115 | 7 | yes | Unknown |
| 5 | atorvastatin + metformin | 223 | 7,705 (6,186–9,458) | 1,265,197 | 5 | yes | Unknown |
| 6 | hydrochlorothiazide + lisinopril | 211 | 6,943 (5,821–8,206) | 1,242,972 | 2 | yes | Moderate |
| 7 | hydrochlorothiazide + metformin | 206 | 6,493 (5,073–8,164) | 1,123,857 | 10 | yes | Moderate |
| 8 | hydrochlorothiazide + losartan | 176 | 6,377 (5,058–7,915) | 1,092,350 | 6 | yes | Unknown |
| 9 | metformin + metoprolol | 144 | 5,581 (4,504–6,825) | 987,608 | 9 | yes | Unknown |
| 10 | metformin + simvastatin | 142 | 4,984 (3,808–6,393) | 872,563 | 11 | yes | Unknown |

**Drugs 10+**

| # | Pair | Participants | Band prevalence (95% CI) | Band burden | Overall rank | Known to Decagon | DDInter |
|---|---|---|---|---|---|---|---|
| 1 | lisinopril + metformin | 90 | 13,195 (8,976–18,462) | 2,306,320 | 1 | yes | Moderate |
| 2 | metformin + omeprazole | 71 | 12,337 (8,171–17,631) | 2,210,842 | 33 | yes | Unknown |
| 3 | albuterol + fluticasone | 74 | 11,767 (7,578–17,182) | 2,106,173 | 15 | yes | Minor |
| 4 | gabapentin + metformin | 78 | 11,403 (7,914–15,742) | 2,002,066 | 39 | yes | Unknown |
| 5 | lisinopril + omeprazole | 61 | 10,243 (6,394–15,335) | 1,856,156 | 22 | yes | Unknown |
| 6 | amlodipine + metoprolol | 69 | 10,640 (7,147–15,074) | 1,772,772 | 18 | yes | Moderate |
| 7 | gabapentin + omeprazole | 60 | 9,708 (6,419–13,934) | 1,762,459 | 51 | yes | Unknown |
| 8 | metformin + metoprolol | 70 | 9,892 (6,971–13,512) | 1,750,494 | 9 | yes | Unknown |
| 9 | metoprolol + omeprazole | 56 | 9,479 (6,215–13,698) | 1,745,754 | 26 | yes | Unknown |
| 10 | atorvastatin + metformin | 84 | 10,580 (7,004–15,158) | 1,737,463 | 5 | yes | Unknown |

- Every adult band's list is the cardiometabolic core of the overall ranking. Children's pairs are asthma and allergy drugs (albuterol, fluticasone, montelukast, cetirizine). In the 10+ band, omeprazole and gabapentin pairs come up from overall ranks 22–51: about one in eight people on 10+ drugs takes metformin + omeprazole.
- All 59 listed band pairs are known to Decagon, matching the overall top 50.

## Validation

### Positive controls (§11.1)

The AGS 2023 Beers Criteria Table 5 list of `data/positive_controls.csv` (#60), fixed before any score was seen. **Model view:** a control is *recovered* if its `p_r` is in the top 10% of all 207,690 pairs for at least one mapped side effect (ties count against it). **Burden view:** its overall rank and whether it reaches the top 50.

- **Run:** `python -m sief_druco.burden.controls local/crossfit/additive_free_w-none_lr0.01`, 2026-11-20. Tables in `local/data/burden/controls/additive_free_w-none_lr0.01/`. The full per-control report of §11.1 (each of the 787 controls with its labelled mark, best side effect, percentile, participants and rank) is `controls.csv` there; below are the per-entry counts and the controls common enough to matter for burden.

| Entry | Beers interaction → outcome | Controls | Recovered | Known to Decagon | Ranked (≥ 30) | Top 50 |
|---|---|---|---|---|---|---|
| 1 | RAS inhibitor or K-sparing diuretic + another → hyperkalaemia | 153 | 76 | 105 | 0 | 0 |
| 2 | opioid + benzodiazepine → overdose | 182 | 46 | 95 | 0 | 0 |
| 3 | opioid + gabapentinoid → sedation, respiratory depression | 28 | 12 | 19 | 2 | 0 |
| 4 | anticholinergic + anticholinergic → cognitive decline, delirium, falls, fractures | 378 | 108 | 189 | 0 | 0 |
| 6 | lithium + ACEI / ARB / ARNI → lithium toxicity | 16 | 11 | 11 | 0 | 0 |
| 7 | lithium + loop diuretic → lithium toxicity | 4 | 2 | 2 | 0 | 0 |
| 8 | non-selective peripheral alpha-1 blocker + loop diuretic → urinary incontinence | 12 | 1 | 8 | 0 | 0 |
| 9 | phenytoin + trimethoprim-sulfamethoxazole → phenytoin toxicity | 1 | 0 | 1 | 0 | 0 |
| 10 | theophylline + cimetidine → theophylline toxicity | 1 | 0 | 0 | 0 | 0 |
| 11 | theophylline + ciprofloxacin → theophylline toxicity | 1 | 0 | 1 | 0 | 0 |
| 12 | warfarin + amiodarone / ciprofloxacin / macrolide / TMP-SMX / SSRI → bleeding | 11 | 9 | 11 | 0 | 0 |
| **All** | | **787** | **265** | **442** | **2** | **0** |

- **265 of 787 controls are recovered.** Labelled controls: 253 of 442. Unlabelled: 12 of 345. The split reflects the two regimes, since labelled pairs fill most of each side effect's top 10%. For the 13 mapped side effects, 32% of labelled pairs reach the top 10%, but only 0.4% of unlabelled pairs do. Twelve of 345 unlabelled controls is above that rate, but S-regime scores are partly expected to be high (§11.1). No overall pass/fail.
- **Every control is checkable:** each entry has at least one mapped side effect. Entry 5 ("≥ 3 CNS-active drugs") is a drug set and out of scope. Beers drugs not in the 645 give no control (e.g. hydrocodone, codeine, diphenhydramine, escitalopram, sulfamethoxazole, triamterene; the full list is `unchecked.csv`). Phenytoin + TMP-SMX reduces to phenytoin + trimethoprim, since sulfamethoxazole is not in the 645.
- **Burden view:** no control reaches the top 50. Only two are common enough to rank: gabapentin + oxycodone (39 participants, rank 143) and gabapentin + tramadol (47, rank 331). Both are recovered (drowsiness, 95.4th percentile; hypoventilation, 96.6th). The Beers interactions are ones prescribers avoid, so they are rarely co-taken.
- **Footnote, conditions ignored:** RAS inhibitor + potassium-sparing diuretic applies in Beers only with chronic kidney disease stage 3a or worse, and alpha-1 blocker + loop diuretic only in older women. The controls test `p` for the combination, not patient eligibility (§11.1).

The controls with ≥ 10 participants in the pool:

| Entry | Pair | Known to Decagon | Best side effect | Percentile | Recovered | Participants | Rank |
|---|---|---|---|---|---|---|---|
| 3 | gabapentin + tramadol | yes | Hypoventilation | 96.6 | yes | 47 | 331 |
| 3 | gabapentin + oxycodone | yes | drowsiness | 95.4 | yes | 39 | 143 |
| 2 | alprazolam + oxycodone | yes | respiratory failure | 97.5 | yes | 24 | — |
| 2 | alprazolam + tramadol | yes | respiratory failure | 98.7 | yes | 15 | — |
| 3 | tramadol + pregabalin | yes | Hypoventilation | 94.4 | yes | 13 | — |
| 8 | doxazosin + furosemide | yes | Stress incontinence | 88.9 | no | 13 | — |
| 3 | gabapentin + morphine | yes | Hypoventilation | 99.5 | yes | 12 | — |
| 4 | cyclobenzaprine + hydroxyzine | yes | confusion | 99.7 | yes | 11 | — |
| 12 | sertraline + warfarin | yes | cerebral haemorrhage | 99.2 | yes | 11 | — |

### Re-runs: top-50 stability (§11.2)

Each re-run forms its own top 50 by the main rule (≥ 30 participants, MEPS: person-years; the 50 highest pair burdens) and is compared with the main pool's top 50. Every re-run uses the same scores, so a difference comes from exposure alone.

- **Run:** `python -m sief_druco.burden.stability local/crossfit/additive_free_w-none_lr0.01`, 2026-11-20. Tables in `local/data/burden/stability/additive_free_w-none_lr0.01/`.

| Re-run | Tests | Ranked pairs | Overlap of the top 50 | Not observed | Below 30 participants | Ranked lower |
|---|---|---|---|---|---|---|
| 2015–16 share | cycle to cycle | 103 | **44** of 50 | 0 | 0 | 6 |
| `P_` share | cycle to cycle | 244 | **47** of 50 | 0 | 0 | 3 |
| FAERS-era pool (2005–08) | time gap to the labels | 223 | **28** of 50 | 0 | 4 | 18 |
| MEPS 2015–17 | purchases instead of 30-day use | 3,589 | **36** of 50 | 0 | 0 | 14 |

Each main-pool top-50 pair missing from a re-run's top 50 gets one reason: *not observed* (no participant takes it in the re-run), *below 30 participants* (too rare there to rank), or *ranked lower* (ranked, but below 50). The drop-outs, with main rank → re-run rank or reason:

- **2015–16:** glipizide + metformin (35 → 54), gabapentin + metformin (39 → 63), atorvastatin + clopidogrel (43 → 68), gabapentin + lisinopril (44 → 55), losartan + omeprazole (46 → 56), clopidogrel + metoprolol (50 → 69).
- **`P_`:** omeprazole + simvastatin (45 → 52), levothyroxine + amlodipine (48 → 54), amlodipine + simvastatin (49 → 55).
- **FAERS-era:** metformin + metoprolol (9 → 84), levothyroxine + metformin (17 → 53), amlodipine + metformin (20 → 101), losartan + metformin (21 → 186), lisinopril + omeprazole (22 → 63), metoprolol + omeprazole (26 → 100), atorvastatin + losartan (28 → 119), amlodipine + losartan (29 → 155), losartan + metoprolol (31 → 173), hydrochlorothiazide + omeprazole (32 → 60), metformin + omeprazole (33 → 107), atorvastatin + omeprazole (34 → 124), levothyroxine + losartan (36, below 30), levothyroxine + omeprazole (38 → 72), gabapentin + metformin (39, below 30), losartan + simvastatin (42 → 150), atorvastatin + clopidogrel (43 → 73), gabapentin + lisinopril (44, below 30), omeprazole + simvastatin (45 → 51), losartan + omeprazole (46, below 30), amlodipine + omeprazole (47 → 156), amlodipine + simvastatin (49 → 67).
- **MEPS:** losartan + metoprolol (31 → 55), metformin + omeprazole (33 → 53), glipizide + metformin (35 → 52), levothyroxine + losartan (36 → 81), gabapentin + metformin (39 → 69), furosemide + metoprolol (41 → 57), losartan + simvastatin (42 → 62), atorvastatin + clopidogrel (43 → 82), omeprazole + simvastatin (45 → 60), losartan + omeprazole (46 → 83), amlodipine + omeprazole (47 → 68), levothyroxine + amlodipine (48 → 59), amlodipine + simvastatin (49 → 64), clopidogrel + metoprolol (50 → 85).

- **Cycle to cycle and MEPS:** the top 30 holds in every one of these re-runs. Only pairs near the cut at 50 drop out (main ranks 31–50), and they land just below it. Measuring exposure by 4-month purchase rounds instead of 30-day use gives the same core.
- **FAERS-era:** 28 of 50. The drop-outs reach up to rank 9 and are mostly pairs with omeprazole (9 of 22) and losartan (7), which rank far lower in 2005–08 (losartan + metformin 21 → 186). Since `p` is the same in both, the difference comes from exposure, not the scores: the burden ranking moves with what people took over the decade that separates the labels from the main pool.

## Re-run data: MEPS 2015–2017 (light)

The MEPS re-run tests exposure measured by purchases rather than 30-day use (§11.2; #62). This section holds its data and coverage; its top-50 overlap with the main result is under "Re-runs: top-50 stability". MEPS prevalence is used only for that overlap and is not reported as a finding.

- **Data:** Prescribed Medicines files HC-178A (2015), HC-188A (2016) and HC-197A (2017), one record per purchase. 2018–2019 are left out (the spring 2018 questionnaire redesign). 960,625 records of 56,724 person-years (42,007 distinct persons); a person present in two years counts as two person-years. Weights `PERWTyyF / 3`. Point estimates only: no HC-036 linkage, no BRR.
- **Co-prescribed** = both drugs bought in the same purchase round (`PURCHRD`, about 4 months), counted once per person-year. 43,071 person-years have ≥ 2 distinct ingredients, 41,678 of them within one round.
- **Denominator:** the people on the PMED files (195.6 million, weighted), since people with no purchase aren't listed. This rescales every pair equally, so the top-50 overlap doesn't change.
- **Mapping:** `RXDRGNAM` → `IN` by name under §4's rules, combinations split at `-`, with the hand fixes of `data/decagon_overrides.csv`. **92.9% of records reach an `IN`** (research note: 93.1%; its pilot also split names at `/` and ` AND `).
- **Run:** `python -m sief_druco.burden.meps`, 2026-11-16. Tables in `local/data/burden/meps/2015_2016_2017/`.

### Coverage

Share of each row by status, in %, as in §6. A mention here is a purchase record × ingredient, not a 30-day drug mention, so the rows can't be compared with NHANES one to one.

| Row | Total | Covered | Not in the model | Mapping failed |
|---|---|---|---|---|
| Distinct ingredients | 690 | 44.6 | 31.0 | 24.3 |
| Distinct ingredients, weighted by their users | 932.1 M | 77.2 | 16.7 | 6.1 |
| Mentions | 1,053,431 | 79.3 | 14.3 | 6.5 |
| Mentions, weighted | 3,619.1 M | 79.0 | 14.8 | 6.2 |
| Person-years with ≥ 2 drugs | 44,215 | 34.7 | 37.5 | 27.8 |
| Person-years with ≥ 2 drugs, weighted | 154.2 M | 33.9 | 39.0 | 27.0 |
| Person-pairs | 844,726 | 60.3 | 27.2 | 12.5 |
| Person-pairs, weighted | 2,928.0 M | 60.4 | 27.4 | 12.2 |

- *Mapping failed* is far larger than in NHANES (6.5% against 1.2% of mentions), and it comes from the files, not the mapping: MEPS replaces drugs used by fewer than 200,000 people with their therapeutic class, and diabetic supplies carry no name. Over a year of purchases, one masked or nameless record is enough to put a person-year in *mapping failed*, so only a third of person-years have every drug covered.
- 60% of weighted person-pairs can be scored (NHANES main pool: 75%).

**Supplement: mapping failed by side and cause.**

| Side | Cause | Drugs | Mentions | % of mentions | % weighted |
|---|---|---|---|---|---|
| Source | no name (`RXDRGNAM` = -9; almost all diabetic supplies) | 1 | 22,890 | 2.17 | 2.06 |
| Source | name not found (class labels, multivitamins, topicals) | 167 | 45,125 | 4.28 | 4.16 |
| Decagon | Decagon drug reached no `IN` | 0 | 0 | 0 | 0 |

- **What the hand fixes and aliases moved:** against the same run without the overrides, 26,196 of 1,053,431 mentions (2.5%) moved into *covered*. 25,862 came from *not in the model*, led by montelukast 8,109, rosuvastatin 5,346 and allopurinol 4,266. 334 came from Decagon-side *mapping failed* (sevelamer, vitamin K1). 981 moved out: betamethasone and isotretinoin, dropped as provisional joins, plus heparin and sorbitol. 2,294 covered mentions join through an alias `IN`.
- The most common names that don't map are therapeutic-class labels (psychotherapeutic agents 2,586 records, antiviral agents 2,471, CNS agents 1,749, narcotic analgesics 1,729) and multivitamins (1,658). No source-side hand fixes (§4).
- **Footnote:** the file documentation counts 5,496 / 4,473 / 4,226 fewer records with a positive person weight in 2015 / 2016 / 2017 than the released files, where every record has one (324,957 against 330,453 in 2015). The released weights are used.

## Re-run data: FAERS-era pool (NHANES 2005–2008)

The FAERS-era re-run tests the time gap between the labels (TWOSIDES covers AERS Q1 2004 – Q1 2009) and the exposure (the main pool centres on ~2017) (§11.2; #22, #63). This section holds its data and coverage; its top-50 overlap with the main result is under "Re-runs: top-50 stability".

- **Data:** NHANES 2005–2006 (`D`) + 2007–2008 (`E`), interview weights `WTINT2YR × 1/2`, Taylor variance on the stacked strata (44–74; the cycles don't share one). 20,497 participants; 8,166 reported a prescription drug in the past 30 days, 5,830 at least two ingredients (#22 counted 5,258 with ≥ 2 drug codes; a combination product is one code but several ingredients). 15 refused or didn't know (`RXDUSE` 7 / 9). The weights sum to 294.4 million people.
- **Mapping:** §4 unchanged, with the hand fixes of `data/decagon_overrides.csv`. **249 drug codes are not in the main pool's code table: 815 mentions, 3.1% of coded mentions**, 27 of them `a` / `c` / `h` codes; this run reproduces #22's counts exactly. Estriol is the class code `a70248` in 2005–06 (`d05027` in the main pool), so it stays *mapping failed* (no source-side hand fixes).
- **Run:** `python -m sief_druco.burden.stability <hand-off>`, 2026-11-17 on the no-encoder hand-off, repeated on the final scores on 2026-11-20 with the same coverage. The coverage below doesn't depend on the scores.

### Coverage

Share of each row by status, in %, as in §6.

| Row | Total | Covered | Not in the model | Mapping failed |
|---|---|---|---|---|
| Distinct ingredients | 690 | 63.3 | 28.0 | 8.7 |
| Distinct ingredients, weighted by their users | 494.3 M | 84.6 | 14.7 | 0.7 |
| Mentions | 29,471 | 85.5 | 13.3 | 1.2 |
| Mentions, weighted | 500.5 M | 84.2 | 14.7 | 1.1 |
| Participants with ≥ 2 drugs | 5,830 | 54.9 | 42.2 | 2.9 |
| Participants with ≥ 2 drugs, weighted | 102.2 M | 52.8 | 44.5 | 2.8 |
| Person-pairs | 73,209 | 75.4 | 23.7 | 0.9 |
| Person-pairs, weighted | 1,213.3 M | 73.3 | 25.8 | 0.9 |

- Coverage is close to the main pool's (weighted person-pairs covered: 73.3% against 75.1%): fifteen years earlier, the 645 Decagon drugs cover about the same share of co-medication.

**Supplement: mapping failed by side and cause.**

| Side | Cause | Drugs | Mentions | % of mentions | % weighted |
|---|---|---|---|---|---|
| Source | `a` / `c` / `h` code | 43 | 134 | 0.45 | 0.40 |
| Source | no drug code (55555 / 77777 / 99999) | 122 records | 122 | 0.41 | 0.36 |
| Source | retired RxNorm concept | 17 | 103 | 0.35 | 0.34 |
| Decagon | Decagon drug reached no `IN` | 0 | 0 | 0 | 0 |

## Limitations

The limitations of §12, with what this run measured.

- **Coverage:** the model is transductive, so drugs outside Decagon's 645 can't be scored. The 645 cover 75.1% of the weighted co-prescribed person-pairs; 24.3% involve a drug outside them. The gap is drugs outside the 645 (ethinyl estradiol, hydrocodone, escitalopram, …), not mapping, which fails for about 1% (RQ2).
- **Mapping:** 11 Decagon drugs are joined by their Decagon single-drug side effects, not their PubChem structure (structure slips, §4; lithium among them). Three Decagon drugs reach no ingredient.
- **Labels:** `p_r(a, b)` estimates "reported together with r in TWOSIDES (FAERS)", not causation. Unlabelled means "not reported", not "safe". The top side effects are broad reporting terms (pain, fatigue, low blood pressure).
- **Exposure:** co-listing in NHANES means co-medication in the past 30 days, not that a side effect occurred. Burden is predicted exposure to a signal, not cases of harm.
- **Model uncertainty:** burden intervals cover prevalence only (`p` fixed); there is no per-score uncertainty. The S-regime scores (labelled pairs, which carry 99.6% of the burden) are calibrated above p ≈ 0.1 and in per-side-effect totals, but off by up to 1.6× between 0.001 and 0.1 (#40, `docs/results-model.md`).
- **Ranking:** among pairs common enough to rank, nearly all are labelled and their Σ_r `p_r` varies little, so the burden ranking is close to the prevalence ranking (Spearman 0.95). The model's main effect on these results is the gap between labelled and unlabelled pairs.
- **Time gap:** TWOSIDES covers AERS Q1 2004 – Q1 2009, and the pool centres on ~2017. On 2005–2008 exposure, 28 of the top 50 pairs stay; the drop-outs follow changes in co-prescription, mostly omeprazole and losartan pairs.
- **Sums count drugs repeatedly:** a person's burden adds up pairs that share a drug, so it is not a patient risk, and its rise with drug count is mostly the number of pairs.
- **Documentation:** DDInter lists 64% of the pair burden without a grade ("Unknown"), which counts as not documented. "Not clinically documented" therefore mostly means "listed by DDInter without a severity", and only 2.2% of the pair burden is in checkable pairs DDInter doesn't list.
- **No causal claims:** the rankings and the novelty candidates are prioritised hypotheses.
- **Novelty of the method:** no precedent was found for combining a learned, calibrated risk with survey-weighted prevalence (#10; confidence medium–high). This file makes no priority claim. The Google Scholar citing-articles search of `research/prior-coprescription-burden.md` and the check of the Parker et al. 2017 wording behind the 30-participant rule (§7) are left for the paper. Related work: known-risk × prevalence (Sutherland 2015 on NHANES, Holm 2014, Guthrie 2015) and the closest ML work (Kee 2024, Hecker 2023, Dagenais 2024, Patrick 2021).
