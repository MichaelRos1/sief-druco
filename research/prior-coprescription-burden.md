# Prior work: co-prescription prevalence × drug-pair risk

Ticket: #10 (map #1). Answers `10.4` in `PROJECT_CONTEXT.md`: has anyone weighted **predicted** drug-pair (polypharmacy side-effect) risk by real-world co-prescription prevalence to rank pairs by population burden? Also covers work on how common known risky pairs are in NHANES / MEPS / EHR / registries, and on how much real co-medication DDI benchmarks cover.

Date of search: 2026-10-27.

## Short answer

- **Known risk × prevalence is well established.** Many studies take a national survey, EHR or prescription registry, look up each co-medicated pair in a curated interaction source (Drugs.com, Micromedex, Lexicomp, DrugBank, Liverpool checker), and report the share of people exposed to a potential DDI, often by severity, age and sex. Some also rank the most common interacting pairs. NHANES has been used this way at least since Sutherland et al. 2015.
- **Predicted risk × prevalence: no direct precedent found.** No paper turned up that scores co-prescribed pairs with a *learned* side-effect model (Decagon-style or other), multiplies by survey-weighted prevalence, and ranks pairs and side effects by burden. The closest are (a) Sutherland et al. 2015, who show that pairs co-prescribed by ≥0.1% of NHANES participants agree poorly with the pairs the literature discusses, and (b) Correia et al. 2019 / the 2024 three-population follow-up, which compute population interaction risk from co-administration but only for **known** DrugBank interactions.
- **Coverage of real co-medication by DDI benchmarks:** no paper found that measures what share of NHANES / MEPS co-medicated pairs fall inside Decagon / TWOSIDES. Sutherland et al. 2015 is the nearest result (a mismatch between real co-prescription and what is studied). This project's coverage figure (`4.3`, `8.1.2`) looks like a new number.
- **Novelty of `1.1`:** the combination (learned polypharmacy side-effect model + calibrated probabilities + survey-weighted prevalence + "undocumented" share) appears new. The ingredients are not new, and the "share of population exposed to known DDIs" part has a large literature to compare against.

Confidence: medium. The search below was narrower than planned because the search tools failed repeatedly during the session (see *Limits of this search*). A focused follow-up search (Google Scholar, Semantic Scholar) for "co-prescription" + "drug-drug interaction prediction" is cheap and worth doing before the paper's related-work section is written.

## 1. Closest prior work

### 1.1 Sutherland et al. 2015, PLoS One (NHANES, known DDIs) — the closest single paper

Sutherland JJ, Daly TM, Liu X, Goldstein K, Johnston JA, Ryan TP. "Co-prescription trends in a large cohort of subjects predict substantial drug-drug interactions." *PLoS One* 10(3):e0118991, 2015. doi:10.1371/journal.pone.0118991. PMID 25739022. (Authors at Eli Lilly and Sano Informed Prescribing.)

- **Data:** NHANES self-reported prescription data 1999–2010; detailed analysis of the 2009–2010 cycle (690 drugs, 10,537 participants).
- **Risk source:** Drugs.com interaction alerts (known interactions) plus CYP450 substrate / inhibitor annotations.
- **How combined:** per participant, count alerts among co-listed drugs; separately, take every drug pair co-prescribed in ≥0.1% of the population and compare its co-prescription rate with how often the pair is co-discussed in the literature.
- **Findings (from the abstract):** alerts rise with drug count, from 3.3 alerts at 5 medications to 11.7 at 10; 22% of elderly participants take a substrate and an inhibitor of the same CYP enzyme, 4% take several inhibitors of one enzyme; ≥98 unique regimens per 100 elderly people on 3+ drugs; the most-prescribed drugs are not the most co-prescribed in any of 16 classes; **low agreement between co-prescription rate and literature co-discussion**.
- **Relevance:** same data source (NHANES) and the same motivation as `1.6` (what people actually take ≠ what is studied). It does not use predicted risk, does not weight by survey design in the pair ranking as far as the abstract shows, and does not rank side effects. Cite as the main NHANES precedent and as support for the "undocumented burden" argument.

### 1.2 Correia et al. 2019, npj Digital Medicine (city-wide EHR, known DDIs)

Correia RB, de Araújo Kohler L, Mattos MM, Rocha LM. "City-wide electronic health records reveal gender and age biases in administration of known drug–drug interactions." *npj Digital Medicine* 2:74, 2019. doi:10.1038/s41746-019-0141-x. arXiv:1803.03571.

- **Data:** 18 months of public-health-system dispensing records, Blumenau, Brazil (~340k inhabitants).
- **Risk source:** known DDIs from DrugBank, with a subset linked to major adverse reactions.
- **How combined:** co-administration network of drug pairs dispensed to the same patient in overlapping time windows; "interaction risk" = probability that a patient is co-administered a known interacting pair, stratified by sex and age, with a null model controlling for number of drugs.
- **Findings:** about 12% of patients were co-administered a known interacting pair; women had ~60% higher risk than men (~90% for major interactions); patients aged 70–79 dispensed 2+ drugs had 34% risk; polypharmacy alone did not explain the sex/age gaps.
- **Relevance:** closest methodological template for RQ4 (`1.5`) and for "prevalence of risky pairs" weighting. Known interactions only.

### 1.3 Three-population follow-up, BMC Medicine 2024 (EHR, known DDIs)

"Prevalence and differences in the co-administration of drugs known to interact: an analysis of three distinct and large populations." *BMC Medicine*, 2024. doi:10.1186/s12916-024-03384-1. (Same Rocha group line of work; author list not re-checked.)

- **Data:** Blumenau (Brazil), Indianapolis (US), Catalonia (Spain) EHR / dispensing data.
- **How combined:** defines prevalence of co-administration and prevalence of interaction over known DDIs.
- **Findings:** co-administration of known interacting drugs affects 12.51% (Blumenau), 12.12% (Indianapolis) and 10.06% (Catalonia) of individuals.
- **Relevance:** cross-population comparison figures for the "share exposed" headline number.

### 1.4 Machine-learning medication recommendation on MIMIC (EHR, known DDI graph)

GAMENet (Shang et al., AAAI 2019) and SafeDrug (Yang et al., IJCAI 2021) recommend medication sets from MIMIC-III and report a "DDI rate" of real and recommended prescriptions, using a DDI graph built from TWOSIDES, the same FAERS-mined source behind Decagon (`3.4.1`). (From prior knowledge; not re-fetched in this session.)

- **Relevance:** the only line found that puts TWOSIDES-style polypharmacy labels next to real co-prescriptions. It uses known TWOSIDES pairs, not predictions, and the goal is recommendation, not population burden. It shows that real ICU co-prescriptions already contain many TWOSIDES-flagged pairs. Worth citing if MIMIC-IV is used (`2.5`).

## 2. Prevalence of known risky pairs in NHANES / MEPS / registries

All of these combine prevalence with a **known** interaction source. None uses predicted risk.

| Study | Data | Risk source | Headline finding |
|---|---|---|---|
| Sutherland et al. 2015, *PLoS One*, PMID 25739022 | NHANES 1999–2010 | Drugs.com, CYP450 | See 1.1 |
| Igho-Osagie et al. 2023, *Clin Ther* 45(5):390–399, PMID 37032225 | NHANES 2015–Mar 2020, 15,685 adults | Liverpool COVID-19 checker, Lexicomp, FDA fact sheets | 29.3% (survey-weighted) of US adults would have a major-to-contraindicated pDDI if given ritonavir; 60.2% at age ≥60 |
| Meng et al. 2025, *J Am Acad Child Adolesc Psychiatry*, PMID 41423041 | NHANES 2001–Mar 2020, ages 6–24 | Micromedex | 26.0% of psychotropic users exposed to a contraindicated or major DDI |
| Qato et al. 2008, *JAMA* 300(24):2867–78; Qato et al. 2016, *JAMA Intern Med* 176(4):473–82 | NSHAP (US adults 57–85, in-home medication inspection), 2005–06 and 2010–11 | Micromedex / Lexi-Interact | Share at risk of a major DDI rose from about 1 in 25 (2005–06) to about 15% (2010–11). (Prior knowledge; figures not re-fetched.) |
| Correia et al. 2019 and BMC Med 2024 | Blumenau, Indianapolis, Catalonia EHR | DrugBank | 10–12.5% of patients co-administered a known interacting pair |
| PLOS One 2026 national outpatient cross-sectional study, doi:10.1371/journal.pone.0345076 | National outpatient prescriptions (country not re-checked) | Interaction checker with severity grades | Nearly half of prescriptions contained ≥1 interaction; co-prescription network with edges weighted by co-prescription frequency and coloured by severity; most common major pairs involve aspirin with metformin, clopidogrel or hydrochlorothiazide |

MEPS: a PubMed search for MEPS + drug interaction returned only 3 hits, none of which estimates pair-level DDI prevalence (Murphy et al. 2018 *Cancer*, PMID 29645083, describes polypharmacy patterns in cancer survivors; Vyas et al. 2012 *Int J Family Med*, PMID 22900173, links chronic-condition clusters to polypharmacy). MEPS appears rarely used for pair-level DDI exposure, which fits `4.1`'s concern that MEPS annual purchase records do not show concurrent use.

Scandinavian and UK prescription-registry studies (for example the Swedish Prescribed Drug Register and Scottish primary-care data) also report potential-DDI prevalence rising with polypharmacy. They are known from background knowledge only and were not re-checked here, so no figures are given.

## 3. How risk and prevalence are combined in prior work

| Approach | Used by | Comparable to this project's |
|---|---|---|
| Share of people exposed to ≥1 known DDI (binary per person, weighted) | Igho-Osagie 2023, Meng 2025, Qato 2008/2016, Correia 2019, BMC Med 2024 | Patient-level risk (`7.6`), but with binary known labels instead of calibrated probabilities |
| Count of alerts per person vs. number of drugs | Sutherland 2015 | RQ4 by number of concurrent drugs (`1.5`) |
| Rank pairs by co-prescription frequency, then label by severity | Sutherland 2015 (pairs ≥0.1%), PLOS One 2026 network | Burden ranking (`7.5`) with severity (`10.2`), but no probability term |
| Network of co-administered pairs, interaction risk by sex and age with a null model | Correia 2019 | RQ4 stratification |
| prevalence × **predicted, calibrated** per-side-effect probability, summed over side effects | none found | `7.5` itself |

## 4. What this means for `1.1`

- **New:** using a learned, calibrated polypharmacy side-effect model on population co-medication; burden broken down by side effect; the share of top-burden pairs that no known source lists; the coverage of NHANES co-medication by Decagon's 645 drugs.
- **Not new:** "X% of the population takes a known interacting pair", and its rise with age and drug count. The project should compare its known-label numbers with these papers rather than present them as findings.
- **Comparison targets:**
  1. Sutherland 2015: same NHANES source; reproduce "pairs ≥0.1% prevalence" and check how many are Decagon-labelled.
  2. Igho-Osagie 2023 and Meng 2025: survey-weighted NHANES DDI prevalence methods (cycle pooling, weights), useful for `10.1`.
  3. Correia 2019 / BMC Med 2024: sex and age stratification method and the 10–12.5% known-DDI exposure figures.
- **Sanity checks (`7.1`):** the PLOS One 2026 study's most common major pairs (aspirin + clopidogrel, aspirin + metformin, aspirin + hydrochlorothiazide) are candidates for "should rank high" if those drugs are in Decagon.

## Limits of this search

- The search tools (web search, web fetch, shell) failed repeatedly mid-session with "no verdict" errors, so searching stopped early.
- Done: PubMed searches for NHANES + DDI (10 hits, all reviewed by title), NHANES + co-prescription terms (30 hits), MEPS + DDI (3 hits); abstracts read for Sutherland 2015, Igho-Osagie 2023, Meng 2025; web search for predicted DDI × co-prescription prevalence and for Correia 2019.
- Not done: a full search of ML venues (NeurIPS, KDD, AMIA, Bioinformatics, JAMIA) for "Decagon" or "polypharmacy side effect" applied to real medication lists; Scandinavian and UK registry studies; DrugBank / TWOSIDES coverage studies. Items marked "prior knowledge" (GAMENet, SafeDrug, Qato) are from background knowledge, not re-read sources.
- The "no direct precedent" finding means "none found in this search", not "none exists".
