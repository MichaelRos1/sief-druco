# Severity weights for Decagon side effects

Research for issue #14 (point `10.2` in `PROJECT_CONTEXT.md`): can each of Decagon's polypharmacy side effects get a severity weight for `burden(a, b) = Σ_r prevalence(a, b) × p_r(a, b) × severity_r`? Which sources exist, in what form, with what IDs, and how many of the side effects do they cover?

All coverage numbers below were computed on 2026-10-27 from the downloaded files (method in [How the numbers were computed](#how-the-numbers-were-computed)).

**Second pass, 2026-10-28.** This pass did the three "Still to do" items: the hand check of the synonym matches ([section 7](#7-hand-check-of-the-synonym-matches)), the EMA IME list ([section 4](#4-meddra-based-severity-and-the-ema-ime-list)) and 963 vs 964 ([section 6](#6-data-note-963-not-964)). `bio-decagon-combo.tar.gz` is 35.7 MB, more than the 5 MB download limit, so it was not downloaded. Without it, the exact 329 matches and the exact 963 set cannot be rebuilt (see [Blockers](#blockers)). In this pass each new claim is marked [verified: source], [derived] (computed or reasoned here) or [unverified]. Corrections to the first pass are marked **Correction**.

## Short answer

- **The Gottlieb et al. 2015 ranking can be downloaded.** It is an Excel file in the paper's open-access supplement (CC BY 2.0). It gives 2,929 adverse drug reactions (ADRs), each with a continuous severity score from 0.01 to 1. The ADRs have **names only (upper-case MedDRA terms), no IDs**.
- **Coverage of the 963 side effects that pass Decagon's ≥500-pair filter:** 522 (54%) by exact name match. That rises to **851 (88%)** when each Decagon UMLS CUI is expanded to its MedDRA synonyms through SIDER's CUI→MedDRA table. The 851 covered side effects hold **92%** of the filtered positive (pair, side effect) labels.
- **CTCAE (v5.0 / v6.0)** uses MedDRA LLT codes and names. It covers far fewer side effects: 287 of 963 (30%) with synonyms. It also gives **no per-term severity**: every term has its own grade 1–5 scale. At best it supplies a weak proxy, such as "can this term reach grade 4/5?"
- **MedDRA has no severity attribute.** It is a terminology. MedDRA also needs a licence. The **EMA Important Medical Event (IME) list**, version 29.1 (September 2026), is an Excel list of 7,846 MedDRA **PT codes** with names. It can be joined to Decagon CUIs **by code** through SIDER's `meddra.tsv.gz`. On the Decagon side effects that could be checked without the combo file, it flags 46% (260 of 561 in `effectcategories`). It agrees well with Gottlieb: an IME ADR outranks a non-IME ADR with probability 0.85. Its count on the 963 is blocked (see [section 4](#4-meddra-based-severity-and-the-ema-ime-list)).
- **Hand check of the synonym matches:** the same-CUI synonym path made **no wrong matches** in 167 checked combo side effects. All wrong or neighbouring matches came from the second path, SIDER's LLT→PT parent mapping. That path is rare among combo side effects (6 of 173) but makes up 49 of 121 synonym matches in a control set, of which 5 were wrong and 5 were neighbours ([section 7](#7-hand-check-of-the-synonym-matches)). The exact 329 could not be rebuilt.
- **963 vs 964:** still open (blocked). The Decagon authors' own loader counts **distinct ordered (drug 1, drug 2) keys** per side effect, not unordered pairs or raw rows. Its printed output gives 1,318 side effect types, not the 1,317 found in the first pass ([section 6](#6-data-note-963-not-964)).
- **`bio-decagon-effectcategories.csv`** gives a disease class, not a severity. It covers only 399 of the 963 filtered side effects.

## 1. Gottlieb et al. 2015 crowdsourced ranking

**Source:** Gottlieb A, Hoehndorf R, Dumontier M, Altman RB. "Ranking Adverse Drug Reactions With Crowdsourcing." *J Med Internet Res* 2015;17(3):e80. doi:10.2196/jmir.3962, PMC4387295. The paper is open access under CC BY 2.0, which covers its multimedia appendices.

**How it was made (from the paper):**
- ADR list: the 2,929 "common" ADRs from **SIDER2** (October 2012 version) that also appear in FDA AERS. The paper says ADRs from SIDER2, AERS and OFFSIDES "are coded using the medical dictionary for regulatory activities (MedDRA) terminology". The paper never mentions UMLS CUIs.
- 2,589 Amazon Mechanical Turk workers (US only) made 126,512 pairwise "which is worse?" judgements. After quality control, 124,513 comparisons remained. A linear-programming model (CPLEX) turned them into one global ranking.
- Validation: Spearman ρ = .53 with the relative death rate per ADR in AERS (P < .001); ρ = .35 with "life-threatening" and .34 with "hospitalization" outcomes. Reproducibility across duplicate batches: ρ = .71.
- Most severe: cardiac arrest. Least severe: elevated mood.

**Downloadable files (JMIR article page, Multimedia Appendices):**

| File | Content |
|---|---|
| `jmir_v17i3e80_app4.xlsx` | **Table S2, the ranked ADR list.** Sheet "Ranked ADRs", 2,929 rows, columns `Name`, `Rank score`, `Rank Stdev (% out 2929)`. |
| `jmir_v17i3e80_app2.xlsx` | Table S1, the raw pairwise worker comparisons. |
| `jmir_v17i3e80_app3.pdf` | Supplementary methods, including the LP formulation. |
| app1 / app5 (png), app6 / app7 (xlsx) | Example task screen; similarity plot; novel severe OFFSIDES ADRs; genes vs. severe ADRs. |

Direct download that worked: `https://www.jmir.org/api/download?alt_name=jmir_v17i3e80_app4.xlsx`. The PMC and Europe PMC `bin/` links returned HTML instead of the file.

**Score format, as observed in the file:**
- `Rank score` is a float in [0.01, 1.0]. It is 1.0 for CARDIAC ARREST, about 0.98 for BONE CANCER METASTATIC, and 0.01 for ELEVATED MOOD. The median over all 2,929 ADRs is 0.494, so the scores are spread roughly uniformly, like a rank percentile.
- `Rank Stdev` is the uncertainty of the rank as a percentage of 2,929, for example 0.9% for cardiac arrest and 4.8% for elevated mood.
- Names are upper case (e.g. `LEFT VENTRICULAR FAILURE`, `HIV INFECTION`). There is no MedDRA code and no UMLS CUI, so the file can only be joined **by name**.
- The paper's caption calls Table S2 a list "with their reported frequency", but the file has no frequency column, and the text says frequency was not used.

**Coverage of Decagon side effects:**

| Decagon set | Exact name match | Name + SIDER MedDRA synonyms of the CUI |
|---|---|---|
| 963 side effects with ≥500 pairs (the modelled set) | 522 (54.2%) | **851 (88.4%)** |
| All 1,317 side effect types in `bio-decagon-combo.csv` | 710 (53.9%) | 1,100 (83.5%) |

- The 851 covered side effects account for **91.9%** of all (pair, side effect) positives among the 963. Mostly rare side effects are missing.
- Covered scores span 0.04–1.0 with median 0.494, the same as the whole Gottlieb list, so the covered subset is not biased towards mild or severe ADRs.
- Highest covered: asystole (1.0, via its CUI's MedDRA synonym "cardiac arrest"), HIV disease 0.97, intracranial haemorrhage 0.95, coma 0.94, breast cancer 0.94, pulmonary embolism 0.93. Lowest: night sweats 0.07, nasal congestion 0.07, agitation 0.07, aptyalism 0.04.
- **What is missing:** mostly terms that are not drug reactions and that SIDER drug labels do not list. Examples are procedures (appendectomy, colostomy, caesarean section, coronary artery bypass graft), injuries (femur fracture, head injury, animal bite, bruise), social or behavioural items (alcohol consumption, abuse) and some chronic diseases (coeliac disease, primary biliary cirrhosis, diabetic nephropathy). These come from FAERS reports via TWOSIDES. Some are probably indications or confounders rather than drug effects (`3.4.1`).
- Caveat: the ranking dates from 2015, MedDRA versions differ, and the crowd was lay people, not clinicians. The paper states this limitation itself.

## 2. ID systems and how to join them

| Source | Side-effect identifier |
|---|---|
| Decagon (`bio-decagon-combo.csv`) | UMLS CUI (e.g. `C0151714`) plus a free-text name (mixed case, often a MedDRA-like term). |
| Gottlieb 2015 Table S2 | **Name only**, upper-case MedDRA term from SIDER2. |
| CTCAE v5.0 | MedDRA v20.1 **LLT code** plus term (column "CTCAE Term/MedDRA LLT v20.1"). |
| CTCAE v6.0 | MedDRA 28.0 LLT code plus term. |
| SIDER 4.1 `meddra.tsv.gz` | UMLS CUI → MedDRA concept type (`LT`/`PT`), MedDRA ID, name. Free download, **CC BY-SA 4.0** (see the correction below). |
| SIDER 4.1 `meddra_all_se.tsv.gz` | Per drug: label CUI, `LLT`/`PT`, CUI of the MedDRA term, name. "All side effects found on the labels are given as LLT. Additionally, the PT is shown." The LLT→PT relations were taken from UMLS [verified: https://sideeffects.embl.de/media/download/README]. |
| EMA IME list 29.1 | MedDRA **PT code** plus PT name and SOC name (see [section 4](#4-meddra-based-severity-and-the-ema-ime-list)). |

**Correction (licence):** the first pass listed SIDER as CC BY-NC-SA. The SIDER download page says the MedDRA files "are available under a Creative Commons Attribution-Share Alike 4.0 License" and links `creativecommons.org/licenses/by-sa/4.0/`. There is no non-commercial clause [verified: https://sideeffects.embl.de/download/]. `https://sideeffects.embl.de` worked over HTTPS (HTTP/2 200) on 2026-10-28. The TLS error seen on stitch.embl.de did not occur here.

- **Bridge used here:** Decagon CUI → all MedDRA names for that CUI in SIDER 4.1 (`meddra.tsv.gz` and `meddra_all_se.tsv.gz`) → case-insensitive name match. 948 of the 963 Decagon CUIs appear in SIDER's table; 896 have a MedDRA PT there.
- Decagon's own name matches a Gottlieb name exactly for 522 side effects. The SIDER synonyms add 329 more (e.g. "asystole" → "cardiac arrest").
- **Risk of many-to-one hits:** a synonym can pull a Decagon CUI onto a broader or neighbouring Gottlieb term. The hand check in [section 7](#7-hand-check-of-the-synonym-matches) found that this happens only on the LLT→PT parent path (`meddra_all_se.tsv.gz`), not on the same-CUI path (`meddra.tsv.gz`) [derived].
- **Decagon's names are lossy; the CUI is the key.** Several Decagon names do not say what their CUI means. Examples: `Clotting` is CUI C0005779, whose MedDRA PT is Coagulopathy. `bone marrow failure` is C0030312, which holds the PTs Pancytopenia and Bone marrow failure. In `bio-decagon-mono.csv`, `alanine aminotransferase` is C0855618, which is "Alanine aminotransferase decreased" [verified: SIDER `meddra.tsv.gz` rows for these CUIs]. A match through the CUI's MedDRA names can be right even when it looks wrong next to the Decagon name.
- **One CUI can hold several MedDRA PTs.** Examples: C0018790 holds Cardiac arrest and Ventricular asystole, and C0234376 holds Action tremor and Intention tremor [verified: SIDER `meddra.tsv.gz`]. When several of these PTs are in Gottlieb, the CUI gets more than one score. This happened for 3 of 173 combo side effects, with score spreads of 0.02 to 0.12 [derived].
- UMLS itself (the official CUI ↔ MedDRA crosswalk) needs a UMLS licence, and MedDRA needs an MSSO licence (free for academic or non-commercial users). Neither was needed for the numbers above.

## 3. CTCAE grades

**Source:** NCI DCTD / CTEP "Trial Development" page, which offers `ctcae-v6.0.xlsx` (2025), `ctcae-v5.0.xlsx` (2017), v4.03 and v3.0 as Excel. CTCAE v6 applies to newly opened CTEP/DCP studies from 1 January 2026.

**Format, as observed in the files:**
- v5.0 sheet "CTCAE v5.0 Clean Copy": 837 terms. Columns: MedDRA v20.1 code, SOC, term (= MedDRA LLT), Grade 1 to Grade 5 text, definition, navigational note. v6.0 has 850 terms with MedDRA 28.0 LLT codes, plus v5→v6 and v4→v5 mapping sheets.
- **A grade belongs to an event in a patient, not to a term.** Each term defines what grade 1 (mild) to grade 5 (death) means for that event. A "–" marks grades that do not apply to the term. So CTCAE gives no single severity per side effect. Possible derived proxies:
  - highest grade the term allows: among the 287 covered side effects, 140 allow grade 5, 22 stop at grade 4, 96 at 3, 25 at 2 and 4 at 1;
  - lowest grade the term allows: 230 start at grade 1, 34 at 2, 17 at 3 and 6 at 4.
  Both are coarse: about half the covered terms allow grade 5.

**Coverage of the 963 filtered Decagon side effects:** v5.0: 149 by exact name, **287 (29.8%)** with SIDER synonyms; v6.0: 151 / 293. CTCAE is an oncology-trial vocabulary and has only ~840 terms, so many FAERS-style terms are absent.

Gottlieb and CTCAE together cover 866 of the 963 (90%), only 15 more than Gottlieb alone.

## 4. MedDRA-based severity and the EMA IME list

- MedDRA is a hierarchical terminology (SOC > HLGT > HLT > PT > LLT). It has **no severity or seriousness field**. Its SOC can group side effects, much like Decagon's disease classes.
- The closest MedDRA-native option is the **EMA Important Medical Event (IME) terms list**. It gives a binary weight: important medical event or not. It was checked on 2026-10-28, with the results below.

**What it is.**
- The EudraVigilance Expert Working Group developed the list to "facilitate the classification of suspected adverse reactions as well as aggregated data analysis and case assessment" in EU pharmacovigilance. EMA says the list "is for guidance purposes only" [verified: https://www.ema.europa.eu/en/human-regulatory-overview/research-development/pharmacovigilance-research-development/eudravigilance/eudravigilance-system-overview].
- Inclusion is based on the ICH definition of a serious adverse reaction (death, life-threatening, hospitalisation, disability, congenital anomaly) and of an "important medical event". The list explicitly leaves out trauma, exposures and poisonings by external agents, and (now) procedure terms such as PT Liver transplant. Criteria by category include embolism, infarction, haemorrhage, cytopenias, organ failure, seizures and coma, and serious infections [verified: EMA/126913/2021 Rev., 22 September 2025, https://www.ema.europa.eu/en/documents/other/inclusion-exclusion-criteria-important-medical-events-ime-list_en.pdf].
- **It is a seriousness flag, not a severity scale.** It does not grade within a term. A non-IME term can still be serious in a given case: the criteria document gives the example of anaemia that needs a transfusion [verified: same PDF].

**Latest version and maintenance.**
- The latest version is **29.1**: "MedDRA important medical event terms list (IME) - version 29.1", reference EMA/211097/2026, XLSX, 432.24 KB, last updated 18/09/2026. EMA updates it twice a year with MedDRA, X.0 by 31 March and X.1 by 30 September [verified: EudraVigilance system overview page above].
- Download that worked: `https://www.ema.europa.eu/en/documents/other/meddra-important-medical-event-terms-list-ime-version-291_en.xlsx` (HTTP 200, 442,616 bytes).

**Format, as observed in the file** [verified: the 29.1 xlsx]:
- Sheet "29.1 IME List": a header block (ref EMA/211097/2026, "Important medical event terms list (MedDRA version 29.1)"), then the columns `MedDRA Code`, `PT Name`, `SOC Name`, `Comment`, `Added in 29.1`, `Primary SOC Change`. There are **7,846 rows, all 8-digit MedDRA PT codes**. There are no LLTs, no CUIs and no grades.
- Sheet "Deleted PTs": 8 PTs dropped in 29.1 because MedDRA demoted them to LLTs (e.g. Iridocyclitis → under PT Anterior uveitis).

**Licence / terms.**
- The criteria PDF carries "© European Medicines Agency, 2025. Reproduction is authorised provided the source is acknowledged." [verified: criteria PDF above]. EMA's legal notice allows reproduction with acknowledgement but excludes "content supplied by third parties" [verified: https://www.ema.europa.eu/en/about-us/legal-notice].
- The xlsx has no licence text of its own. Its PT codes and names are MedDRA content, and MedDRA is licensed separately by the MSSO, so using them in bulk is subject to MedDRA's terms [unverified: the meddra.org pages could not be read; see [Blockers](#blockers)].

**Coverage of Decagon side effects** (join: Decagon CUI → SIDER `meddra.tsv.gz` MedDRA codes → IME PT code; fallbacks are an exact name match and the SIDER LLT→PT parent name):

| Set | Size | IME-flagged | of which: exact name / PT code of same CUI / LT code of same CUI / PT parent |
|---|---|---|---|
| `effectcategories` CUIs (all are combo side effects; 399 of them are in the 963) | 561 | **260 (46.3%)** | 160 / 85 / 4 / 11 |
| `bio-decagon-mono.csv` side effects (control set, not combo) | 10,184 | 3,074 (30.2%) | 2,952 / 74 / 14 / 34 |

[derived; scripts in `local/data/scratch/severity-weights/scripts/ime_match.py`]
- **Not measured on the 963:** that needs `bio-decagon-combo.csv` (blocked). The 46% above is for a subset of combo side effects that is not random. It is not an estimate for the 963.
- **Version gap:** the SIDER 4.1 codes date from 2015 (an older MedDRA version). A PT that was demoted since then will not match a 29.1 PT code. The 8 "Deleted PTs" in 29.1 show that this happens. Fixing it needs the MedDRA LLT→PT tables, which need a licence [derived].
- **Agreement with Gottlieb:** 927 of the 2,929 Gottlieb names are IME PT names. Their median `Rank score` is 0.657, against 0.424 for the rest. The probability that a random IME ADR outranks a random non-IME ADR (AUC) is **0.845** [derived]. So the binary IME flag and the crowd ranking agree to a large extent.
- FAERS seriousness outcomes (death, hospitalization, life-threatening) could give a data-driven severity per side effect, such as a death rate. Gottlieb used this as validation (ρ = .53 with mortality). It would mean extra FAERS processing.

## 5. Decagon effect categories (`bio-decagon-effectcategories.csv`)

- Columns: `Side Effect` (UMLS CUI), `Side Effect Name`, `Disease Class`. 561 rows.
- Covers **399 of the 963** filtered side effects (41%); its 561 rows include CUIs outside the filter.
- 30 disease classes among the 399. Largest: nervous system disease 60, cardiovascular 52, gastrointestinal 43, musculoskeletal 30, respiratory 27, integumentary 27, cancer 25. The labels have inconsistencies (e.g. both "hematopoietic system disease" and "hematopoietic system diseases"; "psoriatic arthritis" as its own class).
- **It is a body-system grouping, not a severity scale.** It could only feed a hand-made per-class weight, and it covers fewer than half of the side effects.

## 6. Data note: 963, not 964

Counting unique unordered drug pairs per side effect in `bio-decagon-combo.csv` and keeping those with ≥500 pairs gives **963** side effects, not 964. Positions 960–968 of the sorted counts are 508, 507, 507, 504, 502, 498, 492, … The next one, avascular necrosis (C0221247), has 498 pairs. So the filtered set has 963 side effects, and the "964" in the paper and on the SNAP page does not come out of this rule. It matters for the data layer (`5.1`), not for severity. Not checked: whether counting raw rows instead (the file may list a pair more than once) or a different file version gives 964.

**Second pass (2026-10-28): still open, blocked.**
- `http://snap.stanford.edu/decagon/bio-decagon-combo.tar.gz` is **35,688,667 bytes** (HEAD request, HTTP 200). That is over the 5 MB limit, so it was not downloaded. `bio-decagon-combo.csv` cannot be fetched on its own (HTTP 404). So raw rows were **not counted**.
- **What the paper says:** "we focus on predicting the 964 commonly occurring types of polypharmacy side effects that each occurred in at least 500 drug combinations". TWOSIDES "details 1318 side effects types across 63 473 drug combinations". "Overall, it contains 4 651 131 drug combination-side effect associations." [verified: Zitnik et al. 2018, https://pmc.ncbi.nlm.nih.gov/articles/PMC6022705/]. The paper does not say whether a "drug combination" is ordered or unordered, or how duplicates are handled.
- **What the SNAP page says:** only that the graph has "964 different polypharmacy side effects". It says nothing about the filter [verified: http://snap.stanford.edu/decagon/].
- **What the authors' code does:** `polypharmacy/utility.py::load_combo_se` keys each row by the **ordered** string `stitch_id1 + '_' + stitch_id2` and adds side effects to a **set** for that key. The notebook `polypharmacy/Exploratory Analysis.ipynb` then counts, per side effect, the number of such keys (`get_se_counter`). Its saved output reads "Drug combinations: 63473 Side effects: 1318" and "Drug-drug interactions: 4651131" [verified: https://github.com/mims-harvard/decagon/blob/master/polypharmacy/utility.py and https://github.com/mims-harvard/decagon/blob/master/polypharmacy/Exploratory%20Analysis.ipynb]. The notebook does not apply the ≥500 filter. The filter code that gives 964 is not in the public repo [verified: repo tree, `gh api repos/mims-harvard/decagon/git/trees/master?recursive=1`].
- **Consequences** [derived]:
  - The authors' unit is the **distinct ordered (drug 1, drug 2) key** per side effect. Duplicate rows collapse, and A–B and B–A count separately if both appear. "Raw rows" is therefore not the authors' unit. The test to run once the file is available is "distinct ordered keys per side effect ≥ 500", plus raw rows for comparison. The ordered count is never below the unordered count, and it is higher whenever a pair appears in both orders. So avascular necrosis (498 unordered pairs) could reach 500 that way. This is a hypothesis, not a result.
  - **Conflict:** the authors' loader prints 1,318 side effect types. The first pass found 1,317 in the same file. The difference is unexplained and was not checked [unverified].

## 7. Hand check of the synonym matches

**Goal:** check the 329 Gottlieb matches among the 963 that come only from SIDER synonyms.

**Limitation:** the 329 could not be rebuilt. That needs the 963 CUIs and their Decagon names from `bio-decagon-combo.csv`, which is blocked (35.7 MB). Instead, the same matching rule (same normalisation, then SIDER synonyms) was rerun on the two Decagon side-effect sets that are small enough to download:
- **Main sample: `bio-decagon-effectcategories.csv`.** Its 561 CUIs are combo side effects, and 399 of them are among the 963 (first pass). It gives 331 exact matches, **173 synonym-only matches** and 57 unmatched. The synonym share is 31%, close to the 34% (329/963) of the first pass. **All 173 were checked by eye**, which covers more than half of the 329 if the overlap is proportional [derived].
- **Control sample: `bio-decagon-mono.csv`.** These 10,184 CUIs are single-drug side effects, not combo side effects. They give 121 synonym-only matches, all checked by eye. They show how the method fails on term types that are rarer in `effectcategories` [derived].

**Classification.** Each match was put in one of two classes by path:
- **same-CUI:** the Gottlieb name is one of the MedDRA names that `meddra.tsv.gz` lists for the same CUI.
- **PT parent:** the Gottlieb name is the MedDRA PT that `meddra_all_se.tsv.gz` gives for the CUI as an LLT, which is a different CUI.

Each match then got one verdict, judged against the CUI's meaning as shown by its MedDRA names:
- *equivalent*;
- *name lossy*: correct for the CUI, but the Decagon name says something else;
- *broader*: the Gottlieb term is a parent;
- *neighbour*: a related but different concept;
- *wrong*.

The PT/LLT status of the matched Gottlieb term was also recorded.

| Sample | Path | n | equivalent | name lossy | broader | neighbour | wrong |
|---|---|---|---|---|---|---|---|
| effectcategories (combo) | same-CUI | 167 | 160 | 7 | 0 | 0 | **0** |
| effectcategories (combo) | PT parent | 6 | 3 | 0 | 3 | 0 | **0** |
| mono (control) | same-CUI | 72 | 56 | 16 | 0 | 0 | **0** |
| mono (control) | PT parent | 49 | 16 | 0 | 23 | 5 | **5** |

[derived: hand verdicts in `local/data/scratch/severity-weights/scripts/verdicts.tsv`, tallies from `classify.py`]

**Findings** [derived]:
- **Error rate among combo side effects: 0 of 173 wrong; 3 of 173 (1.7%) broader.** The broader ones are acute bronchitis → BRONCHITIS, causalgia → COMPLEX REGIONAL PAIN SYNDROME, and allergic vasculitis → LEUKOCYTOCLASTIC VASCULITIS. All three come from the PT parent path.
- **The same-CUI path made no wrong matches** in either sample (0 of 239). In the combo sample, 160 of its 167 matched Gottlieb terms are a MedDRA **PT** of that same CUI, and 7 are LLTs only. What looks wrong on this path is the Decagon name, not the match. Examples: `Clotting` → COAGULOPATHY, `Mod` → TYPE 2 DIABETES MELLITUS, `peliosis` → PURPURA, and in mono `pregnancy test` → PREGNANCY TEST NEGATIVE, where the CUI really is "pregnancy test negative".
- **All wrong and neighbouring matches come from the PT parent path.** Wrong: thinking abnormal → MENTAL DISABILITY, musculoskeletal discomfort → BACK PAIN, emprosthotonus → OPISTHOTONUS (the opposite posture), Budd-Chiari syndrome → VENOOCCLUSIVE LIVER DISEASE, acid base balance → ELECTROLYTE IMBALANCE. Neighbours are mostly injection-site and mucosal terms. Broader matches drop a qualifier (localised/generalised, facial, unilateral, exacerbated, steroid-induced).
- **Multiple scores per CUI:** 3 combo CUIs hit more than one Gottlieb term: acute bronchitis 0.31/0.43, colpitis 0.33–0.43, retinal disorder 0.44/0.46.
- **What this means for the 329:** if the 329 look like the combo sample, nearly all are usable. The PT parent path is the one to review by hand or drop. Dropping it costs 6 of 173 (3.5%) synonym matches in the combo sample. This could not be confirmed on the 329 themselves (blocked).
- Caveat: the verdicts are one non-clinician reviewer's reading of MedDRA terms. Their boundary cases (e.g. neighbour vs broader) are a judgement call.

## Options for the user (facts only)

| Option | Coverage of 963 | Granularity | IDs / join | Cost |
|---|---|---|---|---|
| A. Gottlieb 2015 `Rank score` | 851 (88%), 92% of labels | Continuous 0.01–1 | Name; needs SIDER CUI→MedDRA bridge. Hand check: same-CUI path 0 wrong of 239 checked; LLT→PT parent path is where the errors are | Low; open licence |
| B. CTCAE max/min allowed grade | 287 (30%) | 1–5 ordinal, coarse | MedDRA LLT; same bridge | Low; weak signal |
| C. EMA IME list 29.1 (binary) | not measured on the 963 (blocked); 260 of 561 `effectcategories` CUIs flagged | 0/1 (seriousness, not severity) | MedDRA PT **code**; join by code via SIDER `meddra.tsv.gz` | Low; EMA reuse with acknowledgement, MedDRA terms unverified |
| D. FAERS outcome rates (death, hospitalization) | potentially all | Continuous | CUI/MedDRA via FAERS | High (FAERS processing) |
| E. Decagon disease class + manual weights | 399 (41%) | Per class | UMLS CUI (native) | Manual; not a severity |
| F. No weights (unweighted sum) | 963 | n/a | n/a | None |

Uncovered side effects need a rule whichever source is used: drop them, give them the median weight, or report weighted and unweighted burden side by side.

## How the numbers were computed

- Downloaded on 2026-10-27:
  - `bio-decagon-combo.tar.gz` and `bio-decagon-effectcategories.tar.gz` from `http://snap.stanford.edu/decagon/`
  - `jmir_v17i3e80_app4.xlsx` from JMIR
  - `ctcae-v5.0.xlsx` and `ctcae-v6.0.xlsx` from `https://dctd.cancer.gov/research/ctep-trials/trial-development/`
  - `meddra.tsv.gz` and `meddra_all_se.tsv.gz` (SIDER 4.1) from `http://sideeffects.embl.de/media/download/`
- Filter: unique unordered (STITCH 1, STITCH 2) pairs per `Polypharmacy Side Effect`, keep ≥500, giving 963.
- Name normalisation: lower case, with runs of spaces, hyphens, commas, apostrophes, slashes and parentheses collapsed to one space. Exact match after normalisation.
- Synonym step: if the Decagon name does not match, try every MedDRA name SIDER lists for that CUI.
- The first pass did not hand-check the matches: 522 come from exact names and 329 from SIDER synonyms. The second pass checked a proxy sample instead (see [section 7](#7-hand-check-of-the-synonym-matches)).
- Scripts were throwaway and are not committed. Re-running takes the downloads above plus about 150 lines of stdlib Python.

**Second pass (2026-10-28).** Downloads and scripts are in `local/data/scratch/severity-weights/`, which is not committed (it holds a `.gitignore` of `*`).
- Downloaded into `local/data/scratch/severity-weights/data/`:
  - `jmir_v17i3e80_app4.xlsx` (151,885 B)
  - SIDER `meddra.tsv.gz` (1,084,933 B) and `meddra_all_se.tsv.gz` (2,381,171 B) from `https://sideeffects.embl.de/media/download/`
  - `bio-decagon-mono.tar.gz` (1,928,066 B) and `bio-decagon-effectcategories.tar.gz` (7,742 B) from `http://snap.stanford.edu/decagon/`
  - IME 29.1 xlsx (442,616 B) and the IME criteria PDF (214,165 B) from EMA
  - `polypharmacy/utility.py` and `Exploratory Analysis.ipynb` from GitHub `mims-harvard/decagon`
- Scripts in `local/data/scratch/severity-weights/scripts/`, stdlib only and run with `python3 -I`:
  - `match.py`: the matching rule; output in `syn_matches.tsv`
  - `verdicts.tsv`: hand verdicts; the default is "equivalent"
  - `classify.py`: tallies
  - `ime_match.py`: IME join
  - `pdftext.py`: crude PDF text extraction, because no PDF tool was installed
- Same name normalisation as the first pass. Synonym order: same-CUI names from `meddra.tsv.gz` first, then the LLT→PT parent names from `meddra_all_se.tsv.gz`. The first pass used both files, but its order of preference was not recorded [unverified].

## Not checked / open

- ~~The EMA IME list: its current download location, its MedDRA version, and its coverage of Decagon.~~ Done in the second pass, except its coverage of the 963 (blocked).
- Hand validation of the exact 329 synonym-based Gottlieb matches. Not possible without `bio-decagon-combo.csv`. A proxy of 173 combo side effects plus 121 control matches was checked instead ([section 7](#7-hand-check-of-the-synonym-matches)).
- 963 vs 964: counting distinct ordered keys and raw rows per side effect (blocked; [section 6](#6-data-note-963-not-964)). Also the 1,317 vs 1,318 side effect types.
- Whether the Gottlieb scores still hold up against newer severity sources. No newer crowdsourced ranking was looked for here.
- MedDRA licence terms for reusing the IME list's PT codes and names.

## Blockers

| What | Exact URL / command | Exact error or size | Blocks |
|---|---|---|---|
| Decagon combo file over the 5 MB limit | `curl -sSIL http://snap.stanford.edu/decagon/bio-decagon-combo.tar.gz` | `HTTP/1.1 200 OK`, `Content-Length: 35688667` (35.7 MB); not downloaded | rebuilding the exact 329; 963 vs 964 (ordered keys, raw rows); IME count on the 963; 1,317 vs 1,318 |
| Uncompressed combo CSV not offered | `curl -sSIL http://snap.stanford.edu/decagon/bio-decagon-combo.csv` | `HTTP/1.1 404 Not Found` | same as above |
| MedDRA licence pages unreadable | `curl -sSL https://www.meddra.org/subscription/subscription-rate` (also `/subscription/process`, `/about-meddra/organisation/ownership`); WebFetch of the first | HTTP 200 but no body text, only the page title "MedDRA" (page is rendered by script) | verifying MedDRA terms for IME reuse |
| Decagon article at the publisher | `curl -sSL https://academic.oup.com/bioinformatics/article/34/13/i457/5045770` | `HTTP 403` | none: PMC copy PMC6022705 used instead |
| IME criteria PDF text | WebFetch of the PDF URL; Read tool on the local PDF | WebFetch: "The text is stored as compressed PDF streams"; Read: "pdftoppm is not installed" | none: read with the stdlib extractor `pdftext.py` |

## Follow-up 2026-10-28: 963 vs 964 on the downloaded Decagon file

`bio-decagon-combo.tar.gz` (35,688,667 bytes) was downloaded to `local/data/raw/decagon/` in the main checkout (git-ignored) and counted with a stdlib script.

- The file has **1,317** side-effect types (not 1,318). [verified: `bio-decagon-combo.csv`]
- With the ≥ 500 filter, **963** types remain, whether counted as raw rows, unordered pairs or ordered (drug 1, drug 2) keys: each pair appears once per side effect, so all three counts are equal. The next type below the cut is avascular necrosis with 498. [verified: same file] **964 cannot be reproduced from the public file** with a ≥ 500 rule; it would need a threshold of 498 or a different data version. The paper's 964 stays unexplained.
- The exact 329 synonym-only matches can now be rebuilt from this file; that rerun is still to do.
