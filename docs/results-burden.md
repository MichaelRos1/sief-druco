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

- **Without the hand fixes** the pipeline reproduces the pilot. The only difference: the pilot counted an alias `IN` and its main `IN` (methylphenidate and dexmethylphenidate) as two drugs, so a third drug paired with each gave two pairs. De-duplicating covered drugs by STITCH ID (§5) merges 28 such pairs. The coverage table matches the pilot's too (15,560 ingredient mentions, 2,941 participants with ≥ 2 drugs, 69.2% of weighted person-pairs covered).
- **With the hand fixes,** 753 more covered pairs, 43 more pass 10 participants and 4 more pass 30. They are the pairs of the drugs #55 joined, mostly montelukast, rosuvastatin and allopurinol. Weighted person-pairs covered rise from 69.2% to 74.5%.
