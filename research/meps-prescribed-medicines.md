# MEPS Prescribed Medicines files 2015–2019: drug fields, concurrent use, weights and pooling

Research for issue #23 (map #1). It collects facts for the MEPS robustness check (`7.1`, `2.5` in `PROJECT_CONTEXT.md`), the counterpart of `research/nhanes-prescription-files.md` (#11). It makes no decisions.

**Status: complete for the ticket's four questions.** Questions 1–3 are answered from AHRQ primary sources. The data-dependent parts (persons with ≥ 2 medicines, how many records reach an RxNorm ingredient, name stability) were counted from the five PMED data files on 2026-10-30 (§2.4, §5; marked **[data 2026-10-30]**). Three side questions stay open (*Left open*); none blocks the ticket's questions.

Sources are AHRQ pages on `meps.ahrq.gov`, read on 2026-10-30 (file documentation as HTML, codebooks as PDF, appendix and pooling documentation, one methodology report). Markers: **[verified 2026-10-30]** = read in the cited source on that date; **[unverified]** = not read in a source, reason given; **derived** = arithmetic on published counts. The 2015–2017 codebooks were read with a stdlib text extractor and the 2018–2019 codebooks with a glyph-offset decoder (no PDF tool is installed); counts quoted from them are the codebook's own unweighted frequencies.

## 1. Short answer

- **Files.** The five PMED files are HC-178A (2015), HC-188A (2016), HC-197A (2017), HC-206A (2018) and HC-213A (2019). The file numbers in the ticket are right. One record = one household-reported purchase (fill or refill) of a prescribed medicine in that calendar year. Persons with no purchase are not on the file. 1,573,416 records in total over the five years (derived).
- **Drug fields.** The PUF has **no Multum drug code**. It has `RXDRGNAM` (Multum generic name, a string; since the 2013 file), `RXNAME` (pharmacy-reported name), `RXNDC` (11-digit NDC, partly imputed) and Multum therapeutic-class codes (`TC1…TC3`, with sub-classes; `PREGCAT`). The drug fields were the same in all five years. What changed is the missing-value code (-9 in 2015–2017, -15 in 2018–2019), the removal of `RXCCC1X–3X` after 2015, and one new flag (`DiabEquip`) in 2018. Masked drugs (orphan drugs and drugs used by fewer than 200,000 people) have `RXNDC` and `RXNAME` set to missing and `RXDRGNAM` replaced by a therapeutic class. About 2–3% of records have `RXDRGNAM` missing and 5–6% `RXNDC` missing.
- **Route to RxNorm.** There is no Multum-code route (#12's `GN` + code) because the code is not released. Two routes remain: `RXDRGNAM` → RxCUI by name, and `RXNDC` → RxCUI. Counted on the data: **93.1% of records reach an ingredient by name** (841 distinct names; combinations are written `A-B`), and the NDC route adds 0.3%, for **93.4%**. What does not map is mostly masking (a therapeutic-class label instead of a name, ~3.8%), diabetic supplies with no name (2.3%), and multivitamins, topicals and vaccines (~0.8%) (§2.4).
- **Concurrent use.** MEPS gives only the **round** of purchase (`PURCHRD`, 1–5; three rounds per person-year), not a purchase date. The PMED file has no fill date, no refill count and no "currently taking" flag. `RXDAYSUP` (days supplied) is "do not know" (-8) on 27–39% of records and was not imputed. The date the person first used the medicine (`RXBEGMM`/`RXBEGYRX`) is a first-use date, not a purchase date. Condition links exist through `LINKIDX` and the link files (`CLNK`/`RXLK`) to the Conditions file.
- **Weights and pooling.** The weight is `PERWT15F … PERWT19F` (full-year person weight), design variables `VARSTR` and `VARPSU`. AHRQ's rule: divide the weight by the number of years pooled. **Pooling 2015–2019 needs the HC-036 linkage file**, because HC-036 states that the 2019 variance structure is incompatible with 2002–2018 (the HC-213A text says otherwise; see §4).
- **Sample size.** No published count exists; counted from the data. Persons with a positive weight and **≥ 2 distinct RxNorm ingredients: 15,033 / 14,256 / 13,782 / 13,655 / 13,120** (2015–2019); with ≥ 2 in the same round: 14,631 / 13,735 / 13,312 / 13,218 / 12,670. Weighted, 142–145M people per year have ≥ 2 distinct drug names. Persons recur in adjacent years (shared panels, about 7,000–7,700 per pair of years), so the five years hold 70,185 distinct persons in 92,035 person-years (all persons on the PMED files) (§5).
- **Instrument change in 2018.** AHRQ states that the MEPS instrument was redesigned in spring 2018, including the section on prescribed medicines "obtained earlier … but no longer taking". The 2015–2017 and 2018–2019 files are therefore not collected the same way.

## 2. Q1: Files and drug fields

### 2.1 Which file per year

File numbers and titles are from each file's documentation page.

| Year | File | Documentation | Records | Records, weight > 0 | Variables | Panels (weight factor) |
|---|---|---|---|---|---|---|
| 2015 | HC-178A | [doc](https://meps.ahrq.gov/data_stats/download_data/pufs/h178a/h178adoc.shtml), [codebook](https://meps.ahrq.gov/data_stats/download_data/pufs/h178a/h178acb.pdf) | 330,453 | 324,957 | 69 | 19 (.46), 20 (.54) |
| 2016 | HC-188A | [doc](https://meps.ahrq.gov/data_stats/download_data/pufs/h188a/h188adoc.shtml), [codebook](https://meps.ahrq.gov/data_stats/download_data/pufs/h188a/h188acb.pdf) | 319,685 | 315,212 | 67 | 20 (.51), 21 (.49) |
| 2017 | HC-197A | [doc](https://meps.ahrq.gov/data_stats/download_data/pufs/h197a/h197adoc.shtml), [codebook](https://meps.ahrq.gov/data_stats/download_data/pufs/h197a/h197acb.pdf) | 310,487 | 306,261 | 63 | 21 (.50), 22 (.50) |
| 2018 | HC-206A | [doc](https://meps.ahrq.gov/data_stats/download_data/pufs/h206a/h206adoc.shtml), [codebook](https://meps.ahrq.gov/data_stats/download_data/pufs/h206a/h206acb.pdf) | 319,666 | 316,110 | 66 | 22 (.49), 23 (.51) |
| 2019 | HC-213A | [doc](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml), [codebook](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213acb.pdf) | 293,125 | 289,675 | 62 | 23 (.50), 24 (.50) |

[verified 2026-10-30] for all cells (documentation sections 1.0–2.0 and 3.2, "Weights"). Release months: HC-178A July 2017, HC-213A July 2021 (document headers); the others were not recorded. The 2015 download page lists ASCII (8.7 MB), SAS transport (9.8 MB) and programming statements ([HC-178A detail page](https://meps.ahrq.gov/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-178A)). The data files were downloaded later the same day for the counts in §2.4 and §5 (see Sources).

**Discrepancy in the "weight > 0" column [data 2026-10-30].**
- 2018 and 2019 match the documentation exactly: 3,556 and 3,450 records with `PERWT` = 0.
- 2015–2017 do not: **every record has `PERWT` > 0** (330,453 / 319,685 / 310,487), against the documented 324,957 / 315,212 / 306,261. The 5,496 / 4,473 / 4,226 records the documentation sets apart carry a positive weight in the released file.
- The documentation read does not explain this. The person counts for 2015–2017 in §5 therefore include those records' persons.

Each file covers calendar year *t* using Rounds 3–5 of the panel that started in *t*−1 and Rounds 1–3 of the panel that started in *t* ([HC-213A §2.6.1.4](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). Consecutive years therefore share a panel, so the same persons appear in two adjacent years: Panel 20 in 2015 and 2016, Panel 21 in 2016 and 2017, and so on.

Person-level companion files (full-year consolidated): HC-181 (2015), HC-192 (2016), HC-201 (2017), HC-209 (2018), HC-216 (2019) [verified 2026-10-30 from the document titles, e.g. [HC-216](https://meps.ahrq.gov/data_stats/download_data/pufs/h216/h216doc.shtml)]. Linking appendices with the condition link files: HC-178I, HC-188I, HC-197I, HC-206I, HC-213I; Conditions files HC-180, HC-190, HC-199, HC-207, HC-214 [verified 2026-10-30, appendix documents, e.g. [HC-213I](https://meps.ahrq.gov/data_stats/download_data/pufs/h213i/h213idoc.shtml)].

### 2.2 Drug identifier fields

Variable lists are taken from each file's Variable-Source Crosswalk ([e.g. HC-213A](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)) and codebook. [verified 2026-10-30]

| Field | Content | 2015 | 2016 | 2017 | 2018 | 2019 |
|---|---|---|---|---|---|---|
| `RXDRGNAM` | Multum Lexicon generic name "most commonly used by prescribing physicians", 60 chars | yes | yes | yes | yes | yes |
| `RXNAME` | medicine name as reported by the pharmacy (imputed), 50 chars | yes | yes | yes | yes | yes |
| `RXNDC` | NDC, 11 characters (imputed) | yes | yes | yes | yes | yes |
| `PREGCAT` | Multum FDA pregnancy category | yes | yes | yes | yes | yes |
| `TC1, TC1S1, TC1S1_1, TC1S1_2, TC1S2, TC1S2_1, TC1S3, TC1S3_1, TC2, TC2S1, TC2S1_1, TC2S1_2, TC2S2, TC3, TC3S1, TC3S1_1` | Multum therapeutic class, up to 3 classes, with sub- and sub-sub-classes | yes | yes | yes | yes | yes |
| `RXCCC1X–RXCCC3X` | CCS condition category linked to the purchase | yes | no | no | no | no |
| `DiabEquip` | record is diabetic equipment/supplies first reported at question PM40 | no | no | no | yes | yes |
| `DRUGIDX` | drug-level ID (person + counter), links the same drug across rounds | yes | yes | yes | yes | yes |
| `LINKIDX` | drug-round-level ID for the condition and event link files | 12 chars | 12 | 12 | 16 | 16 |
| `RXRECIDX` | unique purchase-record ID | 15 chars | 15 | 15 | 19 | 19 |

- **No Multum drug ID.** The crosswalk lists only the name and class variables above as "Cerner Multum, Inc." sources; there is no variable like NHANES `RXDDRGID`. The documentation says `RXDRGNAM` was added with the 2013 file and that earlier names are in separate "Multum Lexicon Addendum Files … 1996–2013" ([HC-213A §2.6.2.2](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). The class codes (`TC1` etc.) are therapeutic categories, not drug codes.
- **Multum version.** The documentation gives no Lexicon version or date. "The Multum variables on each year … reflect the most recent classification available in the year the data were released", and the classification has changed over the years (examples given: antidiabetics and antihyperlipidemics moved between classes between the 1996–2004 and later files), so users "should carefully review the data when … pooling years" ([HC-213A §2.6.3](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). This concerns the class variables; whether `RXDRGNAM` strings are stable across 2015–2019 is not stated in the documentation.
  - **In the data they are stable [data 2026-10-30].** 632–654 distinct names per year, 531 present in both 2015 and 2019. No name is written two ways after normalising case and whitespace.
  - Of the 27 ingredient sets reached by more than one name, nearly all are **route variants** of one drug: `MOMETASONE` / `MOMETASONE NASAL` / `MOMETASONE TOPICAL`; `CIPROFLOXACIN` / `CIPROFLOXACIN OPHTHALMIC`.
  - Only three look like renames between years:
    - `MAGNESIUM/POTASSIUM/SODIUM SULFATES` (2015–17) → `MAGNESIUM SULFATE/POTASSIUM SULFATE/SODIUM SULFATE` (2018–19);
    - `GENTAMICIN OPHTHALMIC` (2015–18) → `GENTAMICIN TOPICAL` (2019);
    - `GLYCERIN` → `GLYCERIN TOPICAL` (2019).
  - Route suffixes collapse to the same ingredient, so topical and systemic forms of a drug are not told apart at the `IN` level. NHANES names (`… TOPICAL`, `… OPHTHALMIC`) have the same issue.
- **Known Multum linkage errors** (diabetic supplies in `TC1S1` = 101 for 1996–2002 and = 37 for 2003–2010; NDC 00169750111 since 2006) are all older than 2015 and affect class codes, not names ([HC-213A §2.6.3](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)).
- **Where the names come from.** Households report a name in each round; the pharmacy follow-back (Pharmacy Component) supplies NDC, name, strength, quantity and days supplied. Household drug mentions are matched to pharmacy records "at a drug (active ingredient) level", after all drugs were coded "using a proprietary database", and unmatched household mentions are statistically matched to donor pharmacy records. So `RXNDC`, `RXNAME` and `RXDRGNAM` can come from a pharmacy record of another person (`PCIMPFLG` = 2, "with-replacement match") ([HC-213A §2.5, §2.6.2.4, §4.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). The share of such records is in the table below.
- **NDC imputation.** `RXFLG` = 1 no imputation, 2 from another pharmacy record, 3 from a secondary proprietary database because the pharmacy NDC was invalid. For 1996–2004 the imputed NDC was not released; for the years here it is part of `RXNDC`. [verified 2026-10-30]

### 2.3 Masking and missing values

Confidentiality masking: "orphan drugs and drugs AHRQ estimated were used by fewer than 200,000 people are masked … unless use of the drug does not reveal specific information about the condition treated (for example, cold remedies)". For these, details are "generally recoded as missing" and `RXNAME` is recoded to the remaining therapeutic class information. When a name is a confidentiality risk, `RXNDC` is replaced by the missing code, `RXDRGNAM` by the Multum therapeutic class, and `RXNAME` is set to the missing code ([HC-213A §2.2, §2.6.2.2](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml); same text in the 2015–2018 documents). [verified 2026-10-30]

Reserved codes: -1 inapplicable, -7 refused, -8 don't know, -14 not yet taken (valid only for `RXBEGYRX`), plus **-9 not ascertained in 2015–2017** and **-15 cannot be computed in 2018–2019** (the -9 code was replaced by -15; both are used as the masking code in the respective years). [verified 2026-10-30, §2.2 of each document]

Unweighted counts of missing/masked values, from the codebooks (frequency tables; the missing value is -9 in 2015–2017 and -15 in 2018–2019):

| Year | Records | `RXDRGNAM` missing | `RXNDC` missing | `RXNAME` missing | `RXFLG` = 3 (NDC from secondary source) |
|---|---|---|---|---|---|
| 2015 | 330,453 | 8,405 (2.5%) | 18,996 (5.7%) | 657 (0.2%) | not read |
| 2016 | 319,685 | 7,504 (2.3%) | 19,363 (6.1%) | 358 (0.1%) | 18,010 (5.6%) |
| 2017 | 310,487 | 6,981 (2.2%) | 19,392 (6.2%) | 337 (0.1%) | 18,473 (5.9%) |
| 2018 | 319,666 | 7,145 (2.2%) | 18,145 (5.7%) | 217 (0.1%) | 19,368 (6.1%) |
| 2019 | 293,125 | 6,379 (2.2%) | 15,152 (5.2%) | 273 (0.1%) | not read |

[verified 2026-10-30; percentages derived.] Masked and truly unascertained cases cannot be told apart in the codebooks. The documentation also says imputed records "may still have missing data" for strength, quantity and similar items, because those are embedded in the NDC.

### 2.4 Can a field reach an RxNorm ingredient (route of #12)?

The route of #12 is "Multum code → RxNorm" (`MMSL_CODE = GN + RXDDRGID`) with name fallback (`research/nhanes-stitch-mapping.md`). For MEPS:

- **Multum code route: not available**, because no Multum drug code is on the PUF (see §2.2).
- **Name route.** `RXDRGNAM` is a Multum generic name, so it is the closest equivalent to NHANES `RXDDRUG`. The #12 name fallback (`rxcui.json?name=…&search=2`, `approximateTerm`) applies. Live probe 2026-10-30 against RxNav (data version `05-Oct-2026`, API `3.1.356`): `name=amoxicillin&search=2` → RxCUI 723; `name=acetaminophen; codeine&search=2` → 817579 (`acetaminophen / codeine`). These two names were chosen by hand; they show only that the call works, not how well it matches MEPS strings. [verified 2026-10-30, live RxNav]
  - Whether `RXDRGNAM` uses a fixed convention for combination products (NHANES: alphabetical, joined by `; `) is not described in the MEPS documentation [unverified: needs a look at the strings in the data].
  - Masked records carry a therapeutic class in place of a drug name, so they cannot be mapped by name (§2.3).
- **NDC route.** `findRxcuiById` with `idtype=NDC` accepts "the CMS 11-digit NDC derivative (example: 16571043111)" and 10-digit forms; it "retrieves only the corresponding RxCUI". The page does not say whether historical NDCs are covered ([RxNav findRxcuiById](https://lhncbc.nlm.nih.gov/RxNav/APIs/api-RxNorm.findRxcuiById.html), [verified 2026-10-30]). Live probe with the two example NDCs printed in the HC-213A ID example (hypothetical records, used only as test input): `00093310905` → RxCUI 308191 (`amoxicillin 500 MG Oral Capsule`, status ACTIVE); `00003010955` → empty result from `rxcui.json?idtype=NDC`, while `ndcstatus.json` reports status OBSOLETE with RxCUI 94406 (`amoxicillin 500 MG Oral Capsule [Trimox]`) and an end date of 2011-04. So obsolete NDCs need the `ndcstatus` / history call. Then the RxCUI goes through `related.json?tty=IN` as in #12. [verified 2026-10-30, live RxNav]
  - The NDC identifies the product (manufacturer, strength, package), including combinations, so it maps to a clinical drug rather than directly to an ingredient.
  - `RXNDC` is imputed for part of the records (above) and missing for 5–6%.
- **Coverage of both routes on the data [data 2026-10-30].** Records of persons with `PERWT` > 0, all five years: 1,566,410. Scripts and output: `local/data/scratch/meps-pmed/` (`meps_map.py`, `meps_ndc.py`, `meps_map.out`, `meps_ndc.out`, `unmapped_check.out`; git-ignored). Files: SAS transport for 2015–2017 and Stata for 2018–2019, because the 2018–2019 `.ssp` files are SAS CPORT, which pandas cannot read.
  - **Name route:**
    - 841 distinct `RXDRGNAM` values; 664 reach ≥ 1 `IN`. 654 match as a whole name (`rxcui.json?name=…&search=2`, then `related.json?tty=IN`, as in #12), and 10 more after splitting at `-`, `/` or ` AND `.
    - Combinations are written **`A-B`** with a hyphen, not NHANES's `A; B`. Examples: `ACETAMINOPHEN-HYDROCODONE`, `HYDROCHLOROTHIAZIDE-LISINOPRIL`. RxNorm's normalised name search matches most of them whole.
    - **93.1%** of records reach an ingredient this way.
  - **What the name route misses:**
    - 2.3% of records have no name. These are almost all **diabetic supplies**: `RXNAME` ONETOUCH, FREESTYLE, BD PEN NEEDL, ACCU-CHEK, …; their NDCs are unknown to RxNorm (`ndcstatus` = ALIEN).
    - 4.6% have a name that does not map:
      - 82.9% of those are **therapeutic-class labels**, the masking of §2.3 (`PSYCHOTHERAPEUTIC AGENTS`, `ANTIVIRAL AGENTS`, `NARCOTIC ANALGESICS`, …). The masking is therefore about 3.8% of records, larger than the "missing" counts in §2.3 suggest.
      - The rest is 54 names: multivitamins, topicals, vaccines, probiotics, and a few drug names with a route suffix that RxNorm does not match (`LOTEPREDNOL OPHTHALMIC`, `CONJUGATED ESTROGENS TOPICAL`).
  - **NDC route as a fallback:**
    - 42,510 of the unreached records carry an NDC (1,322 distinct). `ndcstatus.json` status by record: ALIEN 36,268, ACTIVE 3,669, UNKNOWN 1,359, OBSOLETE 1,214.
    - It rescues 4,649 records (0.3% of all), almost all with a name that did not match.
    - **Total reaching an ingredient: 93.4%.**
  - **Three multi-ingredient names** are not written as combinations: `PANCRELIPASE` → amylase, lipase, protease; and two herbal/allergenic double hits (`CRANBERRY`, `GINGER`). These need the same rule as NHANES multi-`IN` codes.
  - **Join to Decagon's 645** (`decagon_rxnorm.csv` from #16; a first look, not a coverage analysis): 86.7–87.2% of records per year reach an ingredient that is one of Decagon's drugs. The NHANES 2017–18 mention-level figure is 82.6% (#16); the units differ (fills vs. mentions).
- **Licence.** MEPS asks users of the Multum variables to cite the Multum Lexicon as the data source ([HC-213A §2.6.3](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)).

## 3. Q2: Concurrent use

### 3.1 What time information the PMED file has

| Field | What it contains | Caveats |
|---|---|---|
| `PURCHRD` | round (1–5) in which the medicine was obtained | Rounds are panel-specific. A person-year has three rounds: 3, 4, 5 (earlier panel) or 1, 2, 3 (later panel). `PANEL` says which. Round numbers are not calendar intervals; calendar position needs `PANEL` and the reference-period dates (below). |
| `PANEL` | panel number (19/20 in 2015 … 23/24 in 2019) | |
| `RXBEGMM`, `RXBEGYRX` | month and year the person first used the medicine | "first taken" questions "are only asked the first time a prescription is mentioned"; "Values are carried forward from prior rounds". It is **not a purchase date**. |
| `DRUGIDX` | links the same drug across rounds | introduced for the 2009 file (earlier years: `RXNDC`) |
| `LINKIDX` | drug-round-level ID | Same within a round for several fills of one drug; varies across rounds. |
| `RXQUANTY` | quantity dispensed (imputed) | |
| `RXDAYSUP` | days supplied | see below |
| `RXFORM`, `RXFRMUNT`, `RXSTRENG`, `RXSTRUNT` | dosage form, quantity unit, strength | |

[verified 2026-10-30, [HC-213A §2.6.1–2.6.2](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml) and the same sections in the other four documents.]

- **No fill date, no refill field.** The variable list has no date of purchase. The Pharmacy Component collects "the dates the prescription was filled or refilled" ([Methodology Report 37, Introduction](https://meps.ahrq.gov/data_files/publications/mr37/mr37.pdf), 2021 data), but the PUF does not release it. AHRQ: "the data collection design for the HC does not allow separate records to be created for multiple acquisitions" and the file "does not distinguish originals from refills": "all records in the HC were 'unfolded' to create separate records for each original purchase and each refill" ([HC-213A §2.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). The refill information is therefore the **number of records** for one `DRUGIDX` in one `PURCHRD`. [verified 2026-10-30]
- **Number of fills is partly imputed.** The household reports the number of times each drug was obtained in the round. Outliers (checked against round length) were replaced, and where the respondent did not know the number, it was imputed; from the 2011 data a hot-deck conditioned on insurance, age, geography and drug was used. In 2017 (Panel 22 Round 3, Panel 21 Round 5) "more implausibly high numbers of fills were reported than in prior years, and so there was more extensive imputation of number of fills". [verified 2026-10-30, [HC-213A §2.5.1, §4.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml); 2017 sentence from [HC-197A §4.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h197a/h197adoc.shtml)]
- **Calendar window of a round.** The reference period of a round runs from the previous round's interview date to the current interview date (Round 1: from 1 January). The dates are on the consolidated person file: `BEGRFM31`, `BEGRFY31`, `BEGRFM42`, `BEGRFY42`, `BEGRFM53`, `BEGRFY53` and the end-date counterparts, "recoded for confidentiality" in a small number of cases. The same wording is in HC-181, HC-192, HC-201, HC-209 and HC-216 (variable names read in HC-216; the sentence "reference period for Round 1 … began on January 1" read in all five). [verified 2026-10-30, [HC-216](https://meps.ahrq.gov/data_stats/download_data/pufs/h216/h216doc.shtml)] The length of a round in days or months is not stated in the PMED or consolidated documentation read for this note [unverified: would be computed from the begin and end dates in the data].
- **Crossover rounds.** The third round of each panel spans two calendar years. The respondent reports how many times a drug was obtained since the last interview and how many of these were in the previous year; if missing, acquisitions are allocated to the years from the start date, pharmacy records, or in proportion to the round's days in each year ([MR 37, "Editing year in crossover rounds"](https://meps.ahrq.gov/data_files/publications/mr37/mr37.pdf); [HC-213A §2.5.1](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). [verified 2026-10-30]

Records per round (unweighted; from the codebooks; [verified 2026-10-30]):

| Year | R1 | R2 | R3 | R4 | R5 |
|---|---|---|---|---|---|
| 2015 | 37,893 | 77,554 | 83,710 | 72,796 | 58,500 |
| 2016 | 33,074 | 68,397 | 79,151 | 74,902 | 64,161 |
| 2017 | 35,925 | 67,792 | 79,406 | 64,939 | 62,425 |
| 2018 | 36,507 | 68,592 | 85,319 | 68,939 | 60,309 |
| 2019 | 34,117 | 65,512 | 73,812 | 69,523 | 50,161 |

### 3.2 `RXDAYSUP` (days supplied)

- Definition: "Days supplied was first collected and released to the public on the 2010 Prescribed Medicines file. Many pharmacies did not provide this information, and imputation was not attempted in these cases. A value of 999 indicates the medication is to be taken as needed. No edits were implemented to impose consistency between the quantity and days supplied, and no edits were implemented for very high values." ([HC-213A §2.6.2.2](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml); same text in 2015–2018) [verified 2026-10-30]
- Codebook values: -8 (don't know, i.e. not supplied), -7, 0, 1–990, 999.

| Year | Records | -8 (not known) | share | 1–990 | 999 (as needed) | 0 |
|---|---|---|---|---|---|---|
| 2015 | 330,453 | 128,401 | 38.9% | 201,868 | 57 | 10 |
| 2016 | 319,685 | 92,250 | 28.9% | 227,179 | 62 | 26 |
| 2017 | 310,487 | 90,679 | 29.2% | 219,696 | 16 | 21 |
| 2018 | 319,666 | 86,596 | 27.1% | 233,000 | 25 | 7 |
| 2019 | 293,125 | 85,026 | 29.0% | 207,981 | 1 | none listed |

[verified 2026-10-30, codebook frequency tables; shares derived.] `-7` (refused) is 38–168 records per year (117 in 2015 and 2019). The missingness is not a random subset by design: it is "not provided by the pharmacy", so the records missing it are mainly those without a matched pharmacy record (this inference is ours; the documentation does not break `RXDAYSUP` down by `PCIMPFLG`).

### 3.3 Condition links

- `LINKIDX` is the HC-generated drug-round-level ID. The condition-event link file `CLNK` links Conditions-file records (`CONDIDX`) to event records (`EVNTIDX`; for the PMED file the matching variable is `LINKIDX`, since `EVNTIDX` is not on the PMED file). The prescribed-medicines link file `RXLK` links PMED records to other event files (office visits, ER, outpatient, inpatient) ([HC-213I §2.3–2.4](https://meps.ahrq.gov/data_stats/download_data/pufs/h213i/h213idoc.shtml)). [verified 2026-10-30]
- Caveats from AHRQ: conditions are self-reported; a purchase can link to several conditions; one condition can link to several purchases; "not all prescribed medicine purchases link to the condition file" ([HC-213A §6.2](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). One RXLK event record can link to several PMED records and the reverse; allocation is left to the analyst ([HC-213I §3.2](https://meps.ahrq.gov/data_stats/download_data/pufs/h213i/h213idoc.shtml)). [verified 2026-10-30]
- **Condition coding changes inside the 2015–2019 window.** 2015: household-reported conditions coded to ICD-9-CM, aggregated to CCS and released on the PMED file as `RXCCC1X–3X` (up to three per purchase, "in the order in which the conditions were reported", not by importance; 99.72% of events had 0–3 condition records linked) ([HC-178A §2.6.2.6](https://meps.ahrq.gov/data_stats/download_data/pufs/h178a/h178adoc.shtml)). 2016 onward: ICD-10-CM; `RXCCC` is not provided on the 2016 PMED file because "the ICD-10 mapping to CCS codes is still under review" ([HC-188A §2.6.2.6](https://meps.ahrq.gov/data_stats/download_data/pufs/h188a/h188adoc.shtml)). From 2018 CCSR categories (up to three per ICD-10-CM code, v2019.1; v2020.3 in 2019) are on the Conditions file ([HC-206I](https://meps.ahrq.gov/data_stats/download_data/pufs/h206i/h206idoc.shtml), [HC-213I §2.5](https://meps.ahrq.gov/data_stats/download_data/pufs/h213i/h213idoc.shtml)). Condition-based grouping therefore needs the link files for all five years. [verified 2026-10-30]
- Condition names are interviewer-recorded text coded by professional coders; AHRQ: "analysts should not presume this level of precision" ([HC-178A §2.6.2.6](https://meps.ahrq.gov/data_stats/download_data/pufs/h178a/h178adoc.shtml)). [verified 2026-10-30]

### 3.4 Comparability inside 2015–2019

- **Spring 2018 instrument redesign.** "Major modifications to the MEPS instrument design and data collection process, particularly in the events sections of the instrument … introduced in the Spring of 2018 and thus affected data beginning with Round 1 of Panel 23, Round 3 of Panel 22, and Round 5 of Panel 21." To limit the effect on the 2017 file, Panel 22 Round 3 and Panel 21 Round 5 data were "transformed to make them as consistent as possible with data collected under the previous design". AHRQ expects "more complete" event data afterwards ([HC-213I §3.4](https://meps.ahrq.gov/data_stats/download_data/pufs/h213i/h213idoc.shtml)). Methodology Report 37 states that the 2018 redesign was aimed at, among other things, "prescribed medicines that sample members obtained earlier in the reference period but were no longer taking" and diabetes and asthma medications, and says "2017 data to 2018 data (instrument redesign might affect trends)" ([MR 37](https://meps.ahrq.gov/data_files/publications/mr37/mr37.pdf)). [verified 2026-10-30]
- **Other changes:** pharmacy types (`PHARTPn`) are those used by the person from 2018, by anyone in the reporting unit before; higher imputed prices from 2017 (an expenditure-only change); charge/payment no longer collected from households for Rx from Panel 21 Round 5 / Panel 22 Round 3; the diabetic-supplies text also refers to Panel 21 Round 5 / Panel 22 Round 3 ([HC-206A §2.5.1, §2.6.2.3](https://meps.ahrq.gov/data_stats/download_data/pufs/h206a/h206adoc.shtml), [HC-213A §4.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). [verified 2026-10-30]
- **`RXBEGMM`/`RXBEGYRX` coding changed.** `RXBEGMM` = -15 ("cannot be computed") on 179,618 of 293,125 records in 2019, versus 114 of 319,666 in 2018 and, in the earlier files, -1 (inapplicable) on 194,780 (2016), 190,610 (2017) and 200,496 (2015) records. The documentation does not explain this change [unverified: no explanation found in HC-213A]. In the data, 2019 has -15 on 61.3% of records and -1 on only 1.9%, while 2015–2018 have -1 on 60.7–64.7%. So -15 takes the place of "inapplicable" in 2019, a recode of the same share of records; the reason is undocumented [data 2026-10-30]. Bottom-coding of `RXBEGYRX`: 1940 (2015–2018 documents) and 1934 (2019 document); the 2019 codebook range starts at 1938. [verified 2026-10-30, codebooks and documents]

## 4. Q3: Weights, pooling, design variables

### 4.1 Which weight

- `PERWT15F`, `PERWT16F`, `PERWT17F`, `PERWT18F`, `PERWT19F`: the "final person-level weight", one value per person, repeated on every PMED record of that person. Each year's weight is a composite of the two panels (factors in the table in §2.1). The weights were built from nonresponse adjustment over the rounds and raking to CPS control totals for the civilian noninstitutionalized population on 31 December ([HC-213A §3.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). [verified 2026-10-30]
- **Event-file versus person-level use.** The PMED file contains only persons with ≥ 1 purchase. AHRQ: "if the estimate relates to the entire population, this file cannot be used to calculate the denominator … the full year consolidated file must be used for person-level analyses that include both persons with and without prescribed medicine events" ([HC-213A §5.2](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). Summing `PERWT` over records estimates the number of purchases; person-level estimates need one record per person. [verified 2026-10-30]
- The same weight variable on the consolidated file is `PERWT15F…PERWT19F`. The PMED file's own weight is the full-year weight, so no special Rx-subsample weight exists in these documents.
- Persons on the PMED file: key in-scope persons with a positive weight, or members of a family whose key members all have a positive weight (via `FAMWT19F`); 289,675 of 293,125 records (2019) belong to persons with `PERWT19F` > 0 ([HC-213A §2.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). [verified 2026-10-30]

### 4.2 Pooling rule (weights)

HC-036 §4.0: "It is generally recommended that analysts adjust the analytic weight variable by dividing it by the number of years being pooled. The sum of these adjusted weights represents the average annual population size for the pooled period … this adjustment will have no effect on estimated means, proportions or regression coefficients … if the objective is to produce an estimated total for the entire pooled period … then the analytic weight variable should not be divided by the number of years" ([HC-036 documentation 1996–2019](https://meps.ahrq.gov/data_stats/download_data/pufs/h036/h36u19doc.shtml)). For 2015–2019 that is `PERWT*F / 5`. [verified 2026-10-30]

Other pooling statements in the same document: MEPS samples are "not completely independent because households are drawn from the same sample geographic areas and many persons are in the sample for two consecutive years"; "it is valid to pool multiple years … and keep all observations in the analysis because each year … is designed to be nationally representative". Weight variable names differ by year (`PERWT15F` … `PERWT19F`) and must be renamed to a common name before stacking (HC-036 §3.0 step 2). [verified 2026-10-30]

### 4.3 Design variables and the HC-036 pooled linkage file

- `VARSTR` and `VARPSU` are on every PMED file. AHRQ recommends Taylor-series linearization with a with-replacement design; "one can generally expect to have at least 100 degrees of freedom" for national estimates ([HC-213A §5.4.1](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). [verified 2026-10-30]
- **Strata by year.** 2015 and 2016: the 165 strata of the 2006 NHIS design, `VARSTR` 1001–1165 (documents HC-178A, HC-188A). 2017: Panel 21 still under the old design, Panel 22 under the new one (117 new strata); the 2017 file therefore has 282 strata, new-design strata numbered 2001–2117 ([HC-197A §5.4.1](https://meps.ahrq.gov/data_stats/download_data/pufs/h197a/h197adoc.shtml)). 2018: both panels in the new design, 117 strata, 2001–2117 ([HC-206A §5.4.1](https://meps.ahrq.gov/data_stats/download_data/pufs/h206a/h206adoc.shtml)). 2019: the new design was modified in 2018 and the number of strata fell to 105 (numbered within 2001–2117, with gaps) ([HC-213A §5.4.1](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)). [verified 2026-10-30]
- **Do the annual variance structures pool?** The HC-213A text says "when pooling any year from 2002 or later, one can use the variance strata numbering as is" and that HC-036 is needed for 1996–2001. **The HC-036 documentation says the opposite for 2019.** It "contains the proper variance structure to use … where one or more years are from 1996-2001 or 2019", because "starting in 2002, the annual MEPS public use files were released with a common variance structure (except in 2019)", which "is not compatible with the structure on the annual PUFs released prior to 2002 and also in 2019". Its scenario list ends: "In the last scenario, data from year 2019 are pooled with other years. The variance structure from HC-036 must be used in this circumstance as well. In no circumstance should the variance structure on the annual PUFs be combined with the variance structure on the HC-036 dataset" ([HC-036 documentation 1996–2019, §1.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h036/h36u19doc.shtml)). Both documents are dated July 2021. The later HC-036 releases keep this rule for 2019 and later: "the common variance structure for 2002 to 2018 and the common variance structure for 2019 and later years are not compatible with each other" ([HC-036 documentation 1996–2024, §1.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h036/h36u24doc.shtml)). So for 2015–2019, which includes 2019 together with 2015–2018, the dedicated pooling document requires HC-036; the HC-213A sentence conflicts with it and is not the pooling document. [verified 2026-10-30]
- **HC-036 contents.** One record per person on the full-year person files, keyed on `DUPERSID` + `PANEL`, with the pooled `STRA9619` and `PSU9619` (1996–2019 release; its source files include HC-181, HC-192, HC-201, HC-209 and HC-216). 433,143 records. Panel 22 persons appear twice because the `DUPERSID` format changed in 2018 (panel prefix), once with the 2017 ID and once with the 2018 ID. Procedure: stack the yearly files, rename variables, merge `STRA9619`/`PSU9619` by `DUPERSID` and `PANEL`, then use them in place of `VARSTR`/`VARPSU` ([HC-036 §2.0, §3.0](https://meps.ahrq.gov/data_stats/download_data/pufs/h036/h36u19doc.shtml)). [verified 2026-10-30]
- **Which release.** The current HC-036 page is the 1996–2024 release (August 2026, 483,346 records, variables `STRA9624`/`PSU9624`) ([HC-036 detail page](https://meps.ahrq.gov/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-036)); the variable names change with each release. Documents for the 1996–2019, 1996–2020, 1996–2022 and 1996–2024 releases exist on the site; the 1996–2019, 1996–2020 and 1996–2024 ones were read. Whether the pooled stratum/PSU values for 2015–2019 are identical across releases is not stated [unverified: would need to compare the files].
- **Subpopulations.** Do not subset to the analysis group before variance estimation; create a flag and use the software's domain/subpopulation option on the full file (HC-036 §5.0, [verified 2026-10-30]). For this project, "persons with ≥ 2 prescribed medicines" is such a subpopulation, so the person file with all respondents (the consolidated files, including persons without any Rx) is needed for the variance, not only the PMED file.
- **BRR alternative.** The HC-036BRR file provides 128 half-sample flags (`BRR1–BRR128`) for replicate weights across years ([HC-213A §5.4.2](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml), [HC-036BRR 1996–2019 page](https://meps.ahrq.gov/data_stats/download_data/pufs/h036brr/h36brr19doc.shtml) was listed by search but not read). [verified for HC-213A text; HC-036BRR page unverified: only seen in a search result]

### 4.4 Comparison with the NHANES side

NHANES pools 2-year cycles with weights multiplied by 1/(number of cycles) and uses `SDMVSTRA`/`SDMVPSU` as given (`research/nhanes-prescription-files.md`, §5). MEPS pools single years with weights divided by the number of years and, for 2019, replaces the design variables with those of HC-036. The person unit differs too: an NHANES participant is one interview; a MEPS person can appear in two adjacent years (shared panel), and `DUPERSID` changes format between 2017 and 2018 for the same Panel 22 persons.

## 5. Q4: Sample size

**Counted from the data [data 2026-10-30].** No published source gives the figure (see below), so it was counted on the five PMED files. Script: `local/data/scratch/meps-pmed/meps_counts.py`, `meps_map.py`. Counts are unweighted persons with `PERWT` > 0, except the weighted column. "Names" = distinct `RXDRGNAM`, excluding missing codes. "INs" = distinct RxNorm ingredients after the name route of §2.4, so a combination product counts as two. "Same round" = at least two in one `PURCHRD`, the closest MEPS gets to concurrent use (§3.1).

| Year | Persons on PMED | ≥ 2 `DRUGIDX` | ≥ 2 names | ≥ 2 INs | ≥ 2 INs, same round | ≥ 2 Decagon INs | Weighted persons (M) | Weighted, ≥ 2 names (M) |
|---|---|---|---|---|---|---|---|---|
| 2015 | 19,604 | 14,271 | 14,094 | 15,033 | 14,631 | 13,032 | 197.4 | 145.2 |
| 2016 | 19,013 | 13,889 | 13,692 | 14,256 | 13,735 | 12,515 | 194.9 | 142.0 |
| 2017 | 18,107 | 13,422 | 13,226 | 13,782 | 13,312 | 12,123 | 194.5 | 142.0 |
| 2018 | 17,753 | 13,450 | 13,254 | 13,655 | 13,218 | 12,398 | 196.8 | 144.2 |
| 2019 | 16,856 | 12,986 | 12,779 | 13,120 | 12,670 | 11,964 | 195.3 | 144.8 |

- **Persons repeat across years.** Shared persons in adjacent years: 7,660 (2015/16), 7,057 (2016/17), 7,214 (2017/18, after removing the 2-digit panel prefix that `DUPERSID` gained in 2018), 7,133 (2018/19). All five files together hold 70,185 distinct persons in 92,035 person-years. A raw `DUPERSID` join finds 0 persons shared by 2017 and 2018, so any cross-year person link needs `PANEL` plus the prefix rule (§4.3).
- **Rounds:** `PURCHRD` 1–5 counts over all five years are 177,516 / 347,847 / 401,398 / 351,099 / 295,556 records.
- **Compared with NHANES:** the main NHANES pool has ~10,900 Rx users over 5.2 years (#15). MEPS has ~13,000–15,000 persons with ≥ 2 ingredients *per year*: a much larger sample, with a weaker exposure definition.

Published figures, kept for reference:

- **Not found.** No AHRQ page read for this note gives the number (weighted or unweighted) of persons with ≥ 2 prescribed medicines in a year. The MEPS summary-tables tool (`https://meps.ahrq.gov/mepstrends/hc_pmed/`) is described by a search result as offering "number of persons with purchases", but a direct request returned **HTTP 401 Unauthorized** (both through the fetch tool and `curl`), so it could not be read. The statistical briefs found (#522, #532, #541, #542, #535) were not found to contain the ≥ 2 split. A PMED record is one fill, so record counts are not person counts.
- **What counting would need.** Persons with ≥ 2 distinct prescribed medicines need the PMED data (distinct `DRUGIDX`, or distinct ingredients after mapping), not the published `RXTOT` totals: `RXTOT` on the consolidated file "is a count of all fills and refills" ([HC-216, prescribed medicines variables](https://meps.ahrq.gov/data_stats/download_data/pufs/h216/h216doc.shtml)), so ≥ 2 fills can be one drug filled twice. The download sizes are modest (HC-178A ASCII 8.7 MB, SAS transport 9.8 MB) but were not downloaded here, per the task rules.
- **Closest published figures** [verified 2026-10-30]:

| Year | PMED records | of which `PERWT` > 0 | Persons on the consolidated file | of which `PERWT` > 0 |
|---|---|---|---|---|
| 2015 | 330,453 | 324,957 | 35,427 | 33,893 |
| 2016 | 319,685 | 315,212 | 34,655 | 33,259 |
| 2017 | 310,487 | 306,261 | 31,880 | 30,716 |
| 2018 | 319,666 | 316,110 | 30,461 | 29,415 |
| 2019 | 293,125 | 289,675 | 28,512 | 27,648 |
| Sum | 1,573,416 (derived) | 1,552,215 (derived) | 160,935 (derived) | 154,931 (derived) |

  Sources: PMED documents, section 2.0; consolidated documents [HC-181](https://meps.ahrq.gov/data_stats/download_data/pufs/h181/h181doc.shtml), [HC-192](https://meps.ahrq.gov/data_stats/download_data/pufs/h192/h192doc.shtml), [HC-201](https://meps.ahrq.gov/data_stats/download_data/pufs/h201/h201doc.shtml), [HC-209](https://meps.ahrq.gov/data_stats/download_data/pufs/h209/h209doc.shtml), [HC-216](https://meps.ahrq.gov/data_stats/download_data/pufs/h216/h216doc.shtml), section 2.0 ("These N persons were part of one of the two MEPS panels …").
- **Share with any Rx.** AHRQ Statistical Brief #522: "the percentage of the population that obtained prescribed medicines was fairly constant (ranging from 60.3 to 62.8 percent)" over 2009–2016 ([SB #522](https://meps.ahrq.gov/data_files/publications/st522/stat522.shtml)). This covers 2009–2016, "not shown" by year, so it only brackets 2015–2016 and says nothing about 2017–2019 or about ≥ 2 medicines. [verified 2026-10-30]
- **No per-person average is given:** the number of persons on the PMED file itself is not published, and dividing records by all persons on the consolidated file would mix users and non-users.

## Left open

Resolved on 2026-10-30 from the data: persons with ≥ 2 medicines (§5), the share of records reaching an RxNorm ingredient and the combination format (§2.4), name stability (§2.2), and the `RXBEGMM` -15 recode, whose reason is still undocumented (§3.4). The MEPS summary-tables tool (HTTP 401) is no longer needed for the count.

Still open. None of these is one of #23's four questions; they matter when the MEPS check is run:

1. **Length of a round in days.** Needs the begin/end dates on the consolidated person files (`BEGRFM31` … `ENDRFY53`, HC-181 … HC-216), which were not downloaded. It bounds how "same round" (§5) approximates concurrent use.
2. **HC-036 releases.** Are the pooled stratum/PSU values for 2015–2019 identical across the HC-036 releases (§4.3)? Needs the HC-036 file itself, which the variance step will download anyway.
3. **Weight > 0 discrepancy, 2015–2017** (§2.1): the released files give every record a positive `PERWT`, but the documentation counts 4,226–5,496 fewer such records per year. It could change the 2015–2017 person counts in §5 by up to ~1–2%.
4. **HC-036BRR documentation** was seen only in a search result, not read. Needed only if replicate weights are chosen over Taylor linearization.

## Sources

AHRQ MEPS (all read 2026-10-30 unless noted)

- PMED file documentation: [HC-178A](https://meps.ahrq.gov/data_stats/download_data/pufs/h178a/h178adoc.shtml), [HC-188A](https://meps.ahrq.gov/data_stats/download_data/pufs/h188a/h188adoc.shtml), [HC-197A](https://meps.ahrq.gov/data_stats/download_data/pufs/h197a/h197adoc.shtml), [HC-206A](https://meps.ahrq.gov/data_stats/download_data/pufs/h206a/h206adoc.shtml), [HC-213A](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213adoc.shtml)
- PMED codebooks (PDF, text decoded locally): [h178acb.pdf](https://meps.ahrq.gov/data_stats/download_data/pufs/h178a/h178acb.pdf), [h188acb.pdf](https://meps.ahrq.gov/data_stats/download_data/pufs/h188a/h188acb.pdf), [h197acb.pdf](https://meps.ahrq.gov/data_stats/download_data/pufs/h197a/h197acb.pdf), [h206acb.pdf](https://meps.ahrq.gov/data_stats/download_data/pufs/h206a/h206acb.pdf), [h213acb.pdf](https://meps.ahrq.gov/data_stats/download_data/pufs/h213a/h213acb.pdf)
- Download page, HC-178A: https://meps.ahrq.gov/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-178A
- Full-year consolidated documentation: [HC-181](https://meps.ahrq.gov/data_stats/download_data/pufs/h181/h181doc.shtml), [HC-192](https://meps.ahrq.gov/data_stats/download_data/pufs/h192/h192doc.shtml), [HC-201](https://meps.ahrq.gov/data_stats/download_data/pufs/h201/h201doc.shtml), [HC-209](https://meps.ahrq.gov/data_stats/download_data/pufs/h209/h209doc.shtml), [HC-216](https://meps.ahrq.gov/data_stats/download_data/pufs/h216/h216doc.shtml)
- Linking appendices (CLNK, RXLK): [HC-178I](https://meps.ahrq.gov/data_stats/download_data/pufs/h178i/h178idoc.shtml), [HC-188I](https://meps.ahrq.gov/data_stats/download_data/pufs/h188i/h188idoc.shtml), [HC-197I](https://meps.ahrq.gov/data_stats/download_data/pufs/h197i/h197idoc.shtml), [HC-206I](https://meps.ahrq.gov/data_stats/download_data/pufs/h206i/h206idoc.shtml), [HC-213I](https://meps.ahrq.gov/data_stats/download_data/pufs/h213i/h213idoc.shtml)
- Pooled linkage variance file: [HC-036 detail page (1996–2024)](https://meps.ahrq.gov/data_stats/download_data_files_detail.jsp?cboPufNumber=HC-036), [HC-036 documentation 1996–2019](https://meps.ahrq.gov/data_stats/download_data/pufs/h036/h36u19doc.shtml), [HC-036 documentation 1996–2020](https://meps.ahrq.gov/data_stats/download_data/pufs/h036/h36u20doc.shtml), [HC-036 documentation 1996–2024](https://meps.ahrq.gov/data_stats/download_data/pufs/h036/h36u24doc.shtml)
- Methodology Report 37, "Outpatient Prescription Drugs: Data Collection and Editing in the 2021 MEPS" (January 2024): https://meps.ahrq.gov/data_files/publications/mr37/mr37.pdf (PDF text decoded locally; used for the pharmacy fill-date collection, crossover-round allocation and the 2018 redesign statements)
- Statistical Brief #522: https://meps.ahrq.gov/data_files/publications/st522/stat522.shtml
- Not readable: MEPS summary tables https://meps.ahrq.gov/mepstrends/hc_pmed/ (HTTP 401). Listed by a search only, not read: [HC-036BRR 1996–2019](https://meps.ahrq.gov/data_stats/download_data/pufs/h036brr/h36brr19doc.shtml).

AHRQ MEPS data files (downloaded 2026-10-30 to `local/data/raw/meps/`, git-ignored; SHA-256 in `SHA256SUMS` there)

- SAS transport: `https://meps.ahrq.gov/data_files/pufs/h178assp.zip` (2015), `…/pufs/h188assp.zip` (2016), `…/pufs/h197a/h197assp.zip` (2017)
- Stata: `…/pufs/h206a/h206adta.zip` (2018), `…/pufs/h213a/h213adta.zip` (2019). The 2018–2019 `.ssp` files are SAS CPORT, which pandas cannot read.

NLM RxNav (for the mapping route; fetched/called 2026-10-30)

- findRxcuiById documentation: https://lhncbc.nlm.nih.gov/RxNav/APIs/api-RxNorm.findRxcuiById.html
- Live calls to `https://rxnav.nlm.nih.gov/REST/` (`version.json`, `rxcui.json?idtype=NDC`, `ndcstatus.json`, `rxcui.json?name=…&search=2`, `approximateTerm.json`), data version `05-Oct-2026`.

Project documents: `research/nhanes-prescription-files.md`, `research/nhanes-stitch-mapping.md`, `PROJECT_CONTEXT.md` (`4`, `4.1`, `4.4`, `7.1`).
