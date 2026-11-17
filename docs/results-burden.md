# Burden results

Results of the burden pipeline of `docs/spec-burden.md` (Track 2). This file is filled ticket by ticket on map #51; the final run (#64) adds the rankings, patient groups, novelty tables, positive controls and re-runs, and states the limitations of §12.

All population figures are a share of the whole pooled population (all ages, non-users included), per 100,000 people (§1). The pool's weights sum to 320.1 million people.

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

## Re-run data: MEPS 2015–2017 (light)

The MEPS re-run tests exposure measured by purchases rather than 30-day use (§11.2; #62). This section holds its data and coverage; the top-50 overlap with the main result comes with the final run (#64). MEPS prevalence is used only for that overlap and is not reported as a finding.

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

The FAERS-era re-run tests the time gap between the labels (TWOSIDES covers AERS Q1 2004 – Q1 2009) and the exposure (the main pool centres on ~2017) (§11.2; #22, #63). This section holds its data and coverage; the top-50 overlap with the main result comes with the final run (#64).

- **Data:** NHANES 2005–2006 (`D`) + 2007–2008 (`E`), interview weights `WTINT2YR × 1/2`, Taylor variance on the stacked strata (44–74; the cycles don't share one). 20,497 participants; 8,166 reported a prescription drug in the past 30 days, 5,830 at least two ingredients (#22 counted 5,258 with ≥ 2 drug codes; a combination product is one code but several ingredients). 15 refused or didn't know (`RXDUSE` 7 / 9). The weights sum to 294.4 million people.
- **Mapping:** §4 unchanged, with the hand fixes of `data/decagon_overrides.csv`. **249 drug codes are not in the main pool's code table: 815 mentions, 3.1% of coded mentions**, 27 of them `a` / `c` / `h` codes; this run reproduces #22's counts exactly. Estriol is the class code `a70248` in 2005–06 (`d05027` in the main pool), so it stays *mapping failed* (no source-side hand fixes).
- **Run:** `python -m sief_druco.burden.stability <hand-off>`, 2026-11-17, dry run on the no-encoder hand-off. The coverage below doesn't depend on the scores.

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
