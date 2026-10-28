# NHANES Prescription Medications files: cycles, drug ID fields, pooling rules

Research for issue #11 (map #1). It collects facts for choosing which cycles to pool (`10.1` in `PROJECT_CONTEXT.md`). It makes no decisions.

Sources are CDC/NCHS pages on `wwwn.cdc.gov`, fetched 2026-10-27: the per-cycle data documentation (codebooks), the NHANES weighting and variance tutorials, and the cycle overview briefs. Each claim links to the page it comes from. Numbers marked **derived** are arithmetic on published codebook counts, not figures NCHS publishes.

## 1. Short answer

- **Twelve releases have a Prescription Medications file.** Eleven of them have drug-level records with names and Multum IDs: 1999–2000 through 2015–2016, 2017–2018, and the 2017–March 2020 pre-pandemic file `P_RXQ_RX`.
- **August 2021–August 2023 (`RXQ_RX_L`) has no drug names.** It has only "any Rx in past 30 days" and a banded count (1, 2, 3, 4, 5+). NCHS says names, duration and reason "were not collected in this survey cycle". That cycle cannot give pair-level co-prescription data.
- **The drug fields we need are stable from 2003–2004 onward:** `RXDDRUG` (generic name; multi-ingredient products are joined with `; ` in alphabetical order), `RXDDRGID` (Multum generic drug code, `d#####`), `RXDDAYS` and `RXDCOUNT`. ICD-10-CM reason codes (`RXDRSC1–3`) exist only from 2013–2014. 1999–2002 use older variable names, and `RXDDRGID` was added to them later.
- **Weights for Rx questionnaire data are the interview weights:** `WTINT2YR`, `WTINTPRP` for the pre-pandemic file, and `WTINT4YR` when 1999–2002 are included. Pooling rule: divide 2-year weights by the number of 2-year cycles. Three exceptions: 1999–2002 must start from the 4-year weights; `P_` files count as 3.2 years; and NCHS says combining 2021–2023 with other cycles is "generally not recommended".
- **The 2017–2018 file and the 2017–March 2020 file overlap.** `P_` already contains the 2017–2018 participants, so you use one or the other, never both.
- **Design variables are `SDMVSTRA` and `SDMVPSU`** in every cycle's demographics file. NCHS recommends Taylor-series linearization.
- **Roughly 3,600–4,500 participants per 2-year cycle report at least one prescription drug** (derived). The pre-pandemic file has about 6,800 (derived).

## 2. Which cycles have the file

The file lists come from the Questionnaire data page of each cycle ([example: 2021–2023](https://wwwn.cdc.gov/Nchs/Nhanes/Search/DataPage.aspx?Component=Questionnaire&Cycle=2021-2023)). Every file is SAS transport (`.xpt`) with an HTML codebook at `https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/<begin-year>/DataFiles/<FILE>.htm`. The codebook header gives the date first published.

| Cycle | File | First published | Drug-level? | Multum DB used for release |
|---|---|---|---|---|
| 1999–2000 | [`RXQ_RX`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/RXQ_RX.htm) | Sep 2003 (rev. Jun 2009) | yes | added later (see §3) |
| 2001–2002 | [`RXQ_RX_B`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/RXQ_RX_B.htm) | Apr 2005 (rev. Jun 2009) | yes | added later |
| 2003–2004 | [`RXQ_RX_C`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2003/DataFiles/RXQ_RX_C.htm) | Aug 2007 (rev. Jun 2009) | yes | Lexicon Plus |
| 2005–2006 | [`RXQ_RX_D`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2005/DataFiles/RXQ_RX_D.htm) | Sep 2008 (rev. Jun 2009) | yes | Lexicon Plus |
| 2007–2008 | [`RXQ_RX_E`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2007/DataFiles/RXQ_RX_E.htm) | Apr 2010 | yes | Lexicon Plus |
| 2009–2010 | [`RXQ_RX_F`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2009/DataFiles/RXQ_RX_F.htm) | May 2012 | yes | Dec 2010 |
| 2011–2012 | [`RXQ_RX_G`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/RXQ_RX_G.htm) | Jul 2014 | yes | Dec 2012 |
| 2013–2014 | [`RXQ_RX_H`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/RXQ_RX_H.htm) | Dec 2016 | yes | Dec 2014 |
| 2015–2016 | [`RXQ_RX_I`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2015/DataFiles/RXQ_RX_I.htm) | Jan 2019 | yes | Dec 2016 |
| 2017–2018 | [`RXQ_RX_J`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/RXQ_RX_J.htm) | Mar 2020 | yes | Dec 2018 |
| 2017–Mar 2020 | [`P_RXQ_RX`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_RXQ_RX.htm) | Sep 2021 | yes | Dec 2019 |
| Aug 2021–Aug 2023 | [`RXQ_RX_L`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/RXQ_RX_L.htm) | Sep 2024 | **no** | — |

The drug information file [`RXQ_DRUG`](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1988/DataFiles/RXQ_DRUG.htm) is listed on every cycle's page through 2017–March 2020. It holds one row per `RXDDRGID` with up to 4 drug-level and 6 ingredient-level Multum therapeutic categories (3 nested levels each) and a single/multi-ingredient flag (`RXDINGFL`). It was last revised Sep 2021 using the Dec 2019 Multum classification. NCHS says it links to the Rx files of NHANES III and every cycle from 1999–2000 to 2017–March 2020 through `RXDDRGID`, and advises: "Analyst should use the most recently revised version."

## 3. Drug fields per cycle

Variable lists are taken from each codebook's table of contents.

| Field | 1999–2002 (`RXQ_RX`, `_B`) | 2003–2012 (`_C`…`_G`) | 2013–Mar 2020 (`_H`, `_I`, `_J`, `P_`) | 2021–2023 (`_L`) |
|---|---|---|---|---|
| Any Rx past 30 days / month | `RXD030` | `RXDUSE` | `RXDUSE` | `RXQ033` |
| Generic name | `RXD240B` ("standard generic ingredient name") | `RXDDRUG` | `RXDDRUG` | — |
| Multum generic drug code | `RXDDRGID` (added later) | `RXDDRGID` | `RXDDRGID` | — |
| Old NCHS ingredient code | `NHCODE` | — | — | — |
| Drug class on record | `FDACODE1–6` (FDA/NDC class) | — (use `RXQ_DRUG`) | — (use `RXQ_DRUG`) | — |
| Container seen | `RXQ250` | `RXQSEEN` | `RXQSEEN` | — |
| Days taken | `RXD260` | `RXDDAYS` | `RXDDAYS` | — |
| Reason for use (ICD-10-CM) | — | — | `RXDRSC1–3` codes, `RXDRSD1–3` text | — |
| Number of Rx drugs | `RXD295` | `RXDCOUNT` | `RXDCOUNT` | `RXQ050` (banded 1/2/3/4/5+) |
| Other | `RXQ280` beta-agonist canisters (1999–2000 only) | | | |

### Points that matter for drug mapping (`4.2`)

- **Name format.** "For multi-ingredient products, the ingredients were listed in alphabetical order (i.e., Acetaminophen; Codeine)." ([P_RXQ_RX](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_RXQ_RX.htm)). Products are named at the generic level, with no brand, dose or form, so duplicates can appear per person: the same generic name and code from different brands, forms or doses.
- **`RXDDRGID` prefixes.** Multum codes "always begin with a 'd'". Names not found in Multum get codes starting "a" or "h". Unspecified products with a known class start with "c" followed by the Multum category code. `RXDDRGID` is missing for the special name codes 55555 (unknown drug), 77777 (refused) and 99999 (don't know) ([P_RXQ_RX](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_RXQ_RX.htm); same text in `_C` and `_G`).
- **Unspecified products** are released as "<Multum category> - unspecified", for example "anti-infective - unspecified". In 1999–2002 the form was "<reported name>-unspecified".
- **1999–2002 retrofit.** These cycles were collected with the Master Drug Database (MDDB, Facts and Comparisons). `RXDDRGID` was added later by matching `RXD240B` names to Multum. That mapping is many-to-many: one name can map to several codes (topical vs. systemic betaxolol), and several names can map to one code (four penicillin names → `d00116`). Also "RXDDRGID is not equal to the variable NHCODE" ([RXQ_RX 1999–2000](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/RXQ_RX.htm)).
- **Collection matching quality differs by cycle:** 92% of drugs matched during the interview in 1999–2000, 70% automatically matched in 2003–2004, 92% in 2011–2012, and about 96% in 2017–March 2020 (83% exact, 13% similar). Unmatched entries were edited by NCHS afterwards (QA sections of each codebook).
- **What the Rx file leaves out.** Prescription vitamins, minerals and most prescription supplements are in the Dietary Supplement files (DSQ), not in RXQ. OTC-only products are removed, except OTC insulin and pseudoephedrine products. `RXDCOUNT` includes unknown, refused and don't-know entries ([P_RXQ_RX](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_RXQ_RX.htm)).
- **Days taken** (`RXDDAYS`) is total duration of use, converted to days (7 per week, 30.4 per month, 365 per year). 77777/99999 mean refused/don't know. It is not "days in the past 30".
- **ICD-10-CM reasons** come from what participants reported, "not … diagnoses from physicians". There are up to 3 per drug, and a trailing "P" marks prevention (e.g. `D75.9P`).
- **Licence.** NCHS outputs must cite the Multum Lexicon as the basis of the drug coding. The licence is reproduced in the [RXQ_DRUG](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1988/DataFiles/RXQ_DRUG.htm) appendix.

### File layout

There is one record per reported drug. Participants with no drug have one record with `RXDUSE` = 2. "SEQN … is not a unique identifier"; NCHS tells analysts to convert drug-level data to person-level before merging with DEMO by `SEQN` ([P_RXQ_RX](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_RXQ_RX.htm), Analytic Notes). Interviewers could record up to 20 drugs in 1999–2000. Later codebooks show `RXDCOUNT` maxima of 20–23.

## 4. 2017–March 2020 pre-pandemic vs. 2021–2023

**2017–March 2020 (`P_` files).**
- Field work stopped in March 2020 after 18 of 30 locations in 2019–2020. Those data are not nationally representative on their own, so NCHS combined them with all of 2017–2018 into a pre-pandemic file with special weights ([P_DEMO](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DEMO.htm)).
- The weights "are not appropriate for independent analyses of the 2019-2020 data and will not yield nationally representative results for either the 2017-2018 data alone or the 2019-March 2020 data alone" ([P_DEMO](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/P_DEMO.htm)).
- `SDDSRVYR` = 66. Weights are `WTINTPRP` and `WTMECPRP`: "should be used for all NHANES 2017-March 2020 pre-pandemic data analyses".
- The release has 24 masked strata and 49 masked PSUs. PSU and strata are masked so that the two parts cannot be told apart, and trend comparisons between 2017–2018 and 2019–March 2020 "are not possible and should not be conducted" ([overview brief 2017–2020](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/overviewbrief.aspx?Cycle=2017-2020)).
- The Rx variables are the same as in 2017–2018, coded with Multum Dec 2019.
- `P_RXQ_RX` includes the 2017–2018 participants. Pooling `RXQ_RX_J` together with `P_RXQ_RX` would count them twice. (This follows from the P_DEMO description; NCHS gives no explicit sentence about it.)

**August 2021–August 2023 (`_L` files).**
- `RXQ_RX_L` has only `SEQN`, `RXQ033` and `RXQ050`. "The detailed prescription medication name, duration of use, and reason for use were not collected in this survey cycle." In year 2, participants were asked about COVID-19 medications, but those data are only available through the NCHS Research Data Center ([RXQ_RX_L](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/RXQ_RX_L.htm)).
- Interviews took place in the home or by phone, and a proxy answered for participants under 16 ([RXQ_RX_L](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/RXQ_RX_L.htm)).
- Response rates were 34.6% for the interview and 25.7% for the exam (unweighted). Household screening was multimode ([overview brief 2021–2023](https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/overviewbrief.aspx?Cycle=2021-2023)).
- The interview-only share is much larger than in earlier cycles: 3,073 of 11,933 participants, versus 300–1,300 before ([DEMO_L](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/DEMO_L.htm)).
- The design has 15 strata and 30 PSUs. Weights are `WTINT2YR` and `WTMEC2YR`.

## 5. Weights, pooling and design variables

**Which weight.**
- NCHS: "use the weight that is appropriate for the variable of interest that was collected on the smallest number of respondents" ([weighting tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx)).
- Rx data come from the home interview. The `_L` codebook states the rule outright: interview weights if RXQ is analysed alone, MEC weights if it is merged with exam or lab data, subsample weights for subsample data ([RXQ_RX_L](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2021/DataFiles/RXQ_RX_L.htm)).
- So the interview weights are `WTINT2YR` (2-year cycles), `WTINT4YR` (1999–2002 combined) and `WTINTPRP` (2017–March 2020). All are on the cycle's DEMO file.

**Pooling rules** ([weighting tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/weighting.aspx), "Combining survey cycles"). The tutorial shows the formulas for MEC weights and says "the same structure would apply" to interview weights.
- **2001–2002 onward, without 1999–2000:** multiply each 2-year weight by 1 / (number of 2-year cycles). Example: `if sddsrvyr in (2,3,4) then MEC6YR = 1/3 * WTMEC2YR`.
- **Any pool that includes 1999–2000:** you "must start by using the 4-year weights provided by NCHS for 1999-2002, then include the 2-year weights for each additional 2-year cycle". Example 1999–2004: `2/3 * WTMEC4YR` for cycles 1–2 and `1/3 * WTMEC2YR` for cycle 3. Using the 2-year weights for 1999–2000 and 2001–2002 "would produce incorrect estimates".
- **With 2017–March 2020:** weight each part by its share of the total years, counting `P_` as 3.2 years. Example 2015–March 2020 (5.2 years): `2/5.2 * WTMEC2YR` for `sddsrvyr = 9` and `3.2/5.2 * WTMECPRP` for `sddsrvyr = 66`. The same logic extends to longer pools, e.g. 2013–Mar 2020 = 7.2 years (that extension is our arithmetic; the tutorial's worked example covers only 5.2 years).
- **2021–2023:** "It is generally not recommended to combine the August 2021-August 2023 cycle with other cycles given the 1.5-year gap" (April 2020–July 2021 unobserved, plus pandemic disruptions to health care).
- **Population totals.** Weight scaling affects estimated population totals, and "analysts are advised NOT to use the sum of weights to determine population estimates for a given health condition" because of exclusions and missing data. The tutorial points to its Sample Code module for the recommended method. This matters for `7` step 4, where prevalence is defined as a sum of weights.

**Design variables** ([variance tutorial](https://wwwn.cdc.gov/nchs/nhanes/tutorials/VarianceEstimation.aspx)).
- Every DEMO file carries `SDMVSTRA` (masked variance pseudo-stratum) and `SDMVPSU` (masked variance pseudo-PSU). These masked variance units stand in for the true design units, which are not released.
- NCHS "recommends the use of the Taylor Series Linearization methods for variance estimation in all NHANES surveys" (e.g. SAS `STRATA sdmvstra; CLUSTER sdmvpsu;`, SUDAAN `NEST sdmvstra sdmvpsu`).
- For subgroups (e.g. an age band), do not subset the data before analysis. Use the software's domain or subpopulation option instead, so variance estimation keeps the full design.
- Masked strata and PSUs per release (DEMO Analytic Notes):

| Release | Strata | PSUs |
|---|---|---|
| 2009–2010 | 15 | 31 |
| 2011–2012 | 14 | 31 |
| 2013–2014, 2015–2016, 2017–2018 | 15 | 30 |
| 2017–Mar 2020 | 24 | 49 |
| 2021–2023 | 15 | 30 |

## 6. Demographics fields needed (age bands, drug-count bands)

All are in the cycle's DEMO file (`DEMO`, `DEMO_B` … `DEMO_J`, `P_DEMO`, `DEMO_L`) and are merged on `SEQN`.

| Need | Variable | Notes |
|---|---|---|
| Cycle ID | `SDDSRVYR` | 1 = 1999–2000 … 10 = 2017–2018; 66 = 2017–Mar 2020 |
| Age | `RIDAGEYR` (age at screening) | Top-coded at **85** in 1999–2006 and at **80** from 2007–2008 onward (DEMO_E … DEMO_L, P_DEMO). NCHS says age at screening "should be used for most analyses". |
| Sex | `RIAGENDR` | |
| Interview/exam status | `RIDSTATR` | 1 = interviewed only, 2 = interviewed + MEC examined |
| Interview weight | `WTINT2YR` / `WTINT4YR` (1999–2002 only) / `WTINTPRP` (P only) | |
| MEC weight | `WTMEC2YR` / `WTMEC4YR` / `WTMECPRP` | only if exam or lab data are merged in |
| Design | `SDMVSTRA`, `SDMVPSU` | |
| Drug count | `RXDCOUNT` (Rx file, repeated on every drug record) or a count of rows per `SEQN` | `RXDCOUNT` includes unknown, refused and don't-know drugs, so it can exceed the number of usable `RXDDRGID`s |

## 7. Sample sizes per cycle

- **Participants:** `RIDSTATR` frequencies from each DEMO codebook.
- **Drug records:** total rows in the Rx codebook frequency table.
- **Rx users (derived):** participants minus non-user records (No + Refused + Don't know + Missing on the lead-in item, which equals the "Missing" count of `RXDCOUNT`/`RXD295`). This assumes every DEMO participant has at least one Rx record, as the "eligible sample = all participants" statement implies.
- **Drug records per user (derived)** = drug records ÷ users. Drug records include unknown, refused and don't-know entries.

| Release | Participants (interviewed) | of whom MEC examined | Rx drug records (`RXDUSE`=1) | Rx users (derived) | Records per user (derived) |
|---|---|---|---|---|---|
| 1999–2000 | 9,965 | 9,282 | 9,257 | 3,579 | 2.6 |
| 2001–2002 | 11,039 | 10,477 | 11,558 | 4,175–4,179 | 2.8 |
| 2003–2004 | 10,122 | 9,643 | 12,308 | 4,029 | 3.1 |
| 2005–2006 | 10,348 | 9,950 | 11,431 | 3,774 | 3.0 |
| 2007–2008 | 10,149 | 9,762 | 14,759 | 4,392 | 3.4 |
| 2009–2010 | 10,537 | 10,253 | 14,889 | 4,521 | 3.3 |
| 2011–2012 | 9,756 | 9,338 | 12,957 | 4,009 | 3.2 |
| 2013–2014 | 10,175 | 9,813 | 14,100 | 4,081 | 3.5 |
| 2015–2016 | 9,971 | 9,544 | 13,705 | 4,029 | 3.4 |
| 2017–2018 | 9,254 | 8,704 | 14,434 | 4,045 | 3.6 |
| 2017–Mar 2020 (incl. 2017–2018) | 15,560 | 14,300 | 24,246 | 6,844 | 3.5 |
| Aug 2021–Aug 2023 | 11,933 | 8,860 | — (no drug records) | 6,255 (`RXQ033`=1, direct) | — |

- For 2001–2002 the codebook shows 11,558 "Yes" records but only 11,554 records with an `RXD295` value. Four "Yes" records have no count, so the derived user count is a range.
- 2021–2023 `RXQ050` (persons): 1 drug 1,621; 2 drugs 1,183; 3 drugs 835; 4 drugs 690; 5+ drugs 1,917.
- Example pools (derived user counts, no overlap): 2015–2016 + 2017–Mar 2020 ≈ 10,900 users (5.2 years). The eight 2-year cycles 2003–2018 ≈ 32,900 users (16 years). Swapping 2017–2018 for `P_` in that pool gives ≈ 35,700 users (17.2 years).

## 8. Open points not settled by these sources

- **Exact Multum versions for 2003–2008 releases.** Those codebooks name Lexicon Plus but not a month/year. All codes are re-classified in the current `RXQ_DRUG` (Dec 2019).
- **Whether `RXDDRGID` is stable across cycles.** The same generic drug is expected to keep the same Multum `d` code, but no NCHS sentence guarantees it across database versions. Checking it means comparing `RXDDRUG`↔`RXDDRGID` pairs across the downloaded files.
- **Strata/PSU numbering across cycles** (whether `SDMVSTRA` values are disjoint between cycles, which matters when pooling). Check this in the data. NCHS's pooling examples simply stack cycles and use `SDMVSTRA`/`SDMVPSU` as given.
- **Reliability thresholds for small cells** (NCHS data presentation standards for proportions) are not covered here.
