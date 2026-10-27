# Severity weights for Decagon side effects

Research for issue #14 (point `10.2` in `PROJECT_CONTEXT.md`): can each of Decagon's polypharmacy side effects get a severity weight for `burden(a, b) = Σ_r prevalence(a, b) × p_r(a, b) × severity_r`? Which sources exist, in what form, with what IDs, and how many of the side effects do they cover?

All coverage numbers below were computed on 2026-10-27 from the downloaded files (method in [How the numbers were computed](#how-the-numbers-were-computed)).

## Short answer

- **The Gottlieb et al. 2015 ranking can be downloaded.** It is an Excel file in the paper's open-access supplement (CC BY 2.0). It gives 2,929 adverse drug reactions (ADRs), each with a continuous severity score from 0.01 to 1. The ADRs have **names only (upper-case MedDRA terms), no IDs**.
- **Coverage of the 963 side effects that pass Decagon's ≥500-pair filter:** 522 (54%) by exact name match. That rises to **851 (88%)** when each Decagon UMLS CUI is expanded to its MedDRA synonyms through SIDER's CUI→MedDRA table. The 851 covered side effects hold **92%** of the filtered positive (pair, side effect) labels.
- **CTCAE (v5.0 / v6.0)** uses MedDRA LLT codes and names. It covers far fewer side effects: 287 of 963 (30%) with synonyms. It also gives **no per-term severity**: every term has its own grade 1–5 scale. At best it supplies a weak proxy, such as "can this term reach grade 4/5?"
- **MedDRA has no severity attribute.** It is a terminology. MedDRA also needs a licence. The Important Medical Event (IME) list is the EMA's MedDRA-based list of medically serious terms. It is a candidate binary weight but was not checked here (see [Not checked](#not-checked--open)).
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
| SIDER 4.1 `meddra.tsv.gz` | UMLS CUI → MedDRA concept type (`LT`/`PT`), MedDRA ID, name. Free download, CC BY-NC-SA. |

- **Bridge used here:** Decagon CUI → all MedDRA names for that CUI in SIDER 4.1 (`meddra.tsv.gz` and `meddra_all_se.tsv.gz`) → case-insensitive name match. 948 of the 963 Decagon CUIs appear in SIDER's table; 896 have a MedDRA PT there.
- Decagon's own name matches a Gottlieb name exactly for 522 side effects. The SIDER synonyms add 329 more (e.g. "asystole" → "cardiac arrest").
- **Risk of many-to-one hits:** a synonym can pull a Decagon CUI onto a broader or neighbouring Gottlieb term. Hand-check the synonym-based matches before using them (see [How the numbers were computed](#how-the-numbers-were-computed)).
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

## 4. MedDRA-based severity

- MedDRA is a hierarchical terminology (SOC > HLGT > HLT > PT > LLT). It has **no severity or seriousness field**. Its SOC can group side effects, much like Decagon's disease classes.
- The closest MedDRA-native option is the **EMA Important Medical Event (IME) terms list**: MedDRA PTs that the EMA considers medically important. It would give a binary weight (serious or not serious). It was **not downloaded or checked here**; see [Not checked](#not-checked--open).
- FAERS seriousness outcomes (death, hospitalization, life-threatening) could give a data-driven severity per side effect, such as a death rate. Gottlieb used this as validation (ρ = .53 with mortality). It would mean extra FAERS processing.

## 5. Decagon effect categories (`bio-decagon-effectcategories.csv`)

- Columns: `Side Effect` (UMLS CUI), `Side Effect Name`, `Disease Class`. 561 rows.
- Covers **399 of the 963** filtered side effects (41%); its 561 rows include CUIs outside the filter.
- 30 disease classes among the 399. Largest: nervous system disease 60, cardiovascular 52, gastrointestinal 43, musculoskeletal 30, respiratory 27, integumentary 27, cancer 25. The labels have inconsistencies (e.g. both "hematopoietic system disease" and "hematopoietic system diseases"; "psoriatic arthritis" as its own class).
- **It is a body-system grouping, not a severity scale.** It could only feed a hand-made per-class weight, and it covers fewer than half of the side effects.

## 6. Data note: 963, not 964

Counting unique unordered drug pairs per side effect in `bio-decagon-combo.csv` and keeping those with ≥500 pairs gives **963** side effects, not 964. Positions 960–968 of the sorted counts are 508, 507, 507, 504, 502, 498, 492, … The next one, avascular necrosis (C0221247), has 498 pairs. So the filtered set has 963 side effects, and the "964" in the paper and on the SNAP page does not come out of this rule. It matters for the data layer (`5.1`), not for severity. Not checked: whether counting raw rows instead (the file may list a pair more than once) or a different file version gives 964.

## Options for the user (facts only)

| Option | Coverage of 963 | Granularity | IDs / join | Cost |
|---|---|---|---|---|
| A. Gottlieb 2015 `Rank score` | 851 (88%), 92% of labels | Continuous 0.01–1 | Name; needs SIDER CUI→MedDRA bridge plus a hand check | Low; open licence |
| B. CTCAE max/min allowed grade | 287 (30%) | 1–5 ordinal, coarse | MedDRA LLT; same bridge | Low; weak signal |
| C. EMA IME list (binary) | not measured | 0/1 | MedDRA PT; same bridge | Low; unchecked |
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
- The matches were not hand-checked. 522 come from exact names and 329 from SIDER synonyms; the synonym matches may include broader terms or many-to-one hits.
- Scripts were throwaway and are not committed. Re-running takes the downloads above plus about 150 lines of stdlib Python.

## Not checked / open

- The EMA IME list: its current download location, its MedDRA version, and its coverage of Decagon.
- Hand validation of all 329 synonym-based Gottlieb matches.
- Whether the Gottlieb scores still hold up against newer severity sources. No newer crowdsourced ranking was looked for here.
