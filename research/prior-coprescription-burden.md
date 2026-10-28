# Prior work: co-prescription prevalence × drug-pair risk

Ticket: #10 (map #1). Answers `10.4` in `PROJECT_CONTEXT.md`: has anyone weighted **predicted** drug-pair (polypharmacy side-effect) risk by real-world co-prescription prevalence to rank pairs by population burden? Also covers work on how common known risky pairs are in NHANES / MEPS / EHR / registries, and on how much real co-medication DDI benchmarks cover.

Date of search: 2026-10-27 (first pass); 2026-10-28 (follow-up: ML / informatics venues, Nordic and UK registries, check of items that came from prior knowledge). Exact queries are in *Search log*. Failed requests are in *Blockers*.

Markers: **[verified: URL]** means the claim was checked against that primary source (abstract or full text, as stated). **[verified 2026-10-27]** means it was checked in the first pass and not re-read on 2026-10-28. **[unverified]** means it was not checked against a source.

## Short answer

- **Known risk × prevalence is well established.** Many studies take a national survey, EHR or prescription registry, look up each co-medicated pair in a curated interaction source (Drugs.com, Micromedex, Lexicomp, DrugBank, SFINX, BNF, Liverpool checker), and report the share of people exposed to a potential DDI, often by severity, age and sex. Some also rank the most common interacting pairs. NHANES has been used this way at least since Sutherland et al. 2015. Swedish register studies (Johnell & Klarin 2007, Holm et al. 2014) and Scottish dispensing data (Guthrie et al. 2015) do the same at whole-population scale (section 2).
- **Predicted risk × prevalence: still no direct precedent found.** The follow-up search found no paper that scores co-prescribed pairs with a *learned* side-effect model, multiplies the score by survey-weighted (or any population) prevalence, and ranks pairs and side effects by burden. Learned DDI models *have* been applied to real medication data, but only in three narrower ways:
  1. **Kee et al. 2024** (Heliyon): ranked DeSIDE-DDI predictions (a model trained on TWOSIDES-derived BioSNAP data, the same data family as Decagon), kept top predictions that had enough co-users in one hospital's EHR, and tested one pair (cefpodoxime + chlorpheniramine → lung edema) in a retrospective cohort. Co-use counts were a filter for feasibility, not a weight. [verified: https://doi.org/10.1016/j.heliyon.2024.e24620, full text PMC10831713]
  2. **Hecker et al. 2023** (Pharmaceutics): ran a retrained DeepDDI model (DrugBank DDI types, not side effects) over the medication plans of 627 MS patients and reported the share of patients with ≥1 predicted DDI (81.2%). This is the nearest "predicted DDI prevalence" figure. It is one clinical cohort, with no survey weights, no calibrated probabilities and no burden ranking. [verified: https://doi.org/10.3390/pharmaceutics16010003, abstract]
  3. **Dagenais et al. 2024** (Clin Transl Sci, Pfizer): a "predicted DDI impact score" built from claims-based (Optum) prevalence of exposure to CYP victim substrates for hypothetical new perpetrator drugs. The word "predicted" is used, but the risk side is a known PK classification (UW Drug Interaction Database) scored in points, not a learned model. [verified: https://doi.org/10.1111/cts.13741, full text PMC10915735]
- **Coverage of real co-medication by DDI benchmarks:** still no paper found that measures what share of NHANES / MEPS co-medicated pairs fall inside Decagon / TWOSIDES. The nearest numbers are the GAMENet / SafeDrug "DDI rate" on MIMIC-III: 7.77% and 8.08% of within-visit drug pairs (ATC level 3) are TWOSIDES top-40-severity pairs (section 1.4). That measures the share of real pairs that are *flagged*, not the share of real pairs the benchmark *covers*. This project's coverage figure (`4.3`, `8.1.2`) still looks like a new number.
- **Novelty of `1.1`:** the combination (learned polypharmacy side-effect model + calibrated probabilities + survey-weighted prevalence + "undocumented" share) still appears new. The ingredients are not new. The related-work section must cite Kee 2024 and Hecker 2023 as the closest ML precedents and say how this project differs (population-weighted, nationally representative, per-side-effect burden, calibrated).

Confidence: medium. It is higher than after the first pass because PubMed, Europe PMC, arXiv and OpenAlex (citers of Zitnik et al. 2018) were searched. It is not high: Semantic Scholar and Google Scholar could not be queried (see *Blockers*), and OpenAlex indexes conference papers (NeurIPS, KDD, AMIA) only partly. "No direct precedent" still means "none found in this search", not "none exists".

## 1. Closest prior work

### 1.1 Sutherland et al. 2015, PLoS One (NHANES, known DDIs) — the closest single paper on the data side

Sutherland JJ, Daly TM, Liu X, Goldstein K, Johnston JA, Ryan TP. "Co-prescription trends in a large cohort of subjects predict substantial drug-drug interactions." *PLoS One* 10(3):e0118991, 2015. doi:10.1371/journal.pone.0118991. PMID 25739022. (Authors at Eli Lilly and Sano Informed Prescribing.) [verified 2026-10-27: https://pubmed.ncbi.nlm.nih.gov/25739022/]

- **Data:** NHANES self-reported prescription data 1999–2010; detailed analysis of the 2009–2010 cycle (690 drugs, 10,537 participants).
- **Risk source:** Drugs.com interaction alerts (known interactions) plus CYP450 substrate / inhibitor annotations.
- **How combined:** per participant, count alerts among co-listed drugs; separately, take every drug pair co-prescribed in ≥0.1% of the population and compare its co-prescription rate with how often the pair is co-discussed in the literature.
- **Findings (from the abstract):** alerts rise with drug count, from 3.3 alerts at 5 medications to 11.7 at 10; 22% of elderly participants take a substrate and an inhibitor of the same CYP enzyme, 4% take several inhibitors of one enzyme; ≥98 unique regimens per 100 elderly people on 3+ drugs; the most-prescribed drugs are not the most co-prescribed in any of 16 classes; **low agreement between co-prescription rate and literature co-discussion**.
- **Relevance:** same data source (NHANES) and the same motivation as `1.6` (what people actually take ≠ what is studied). It does not use predicted risk, does not weight by survey design in the pair ranking as far as the abstract shows, and does not rank side effects. Cite as the main NHANES precedent and as support for the "undocumented burden" argument.

### 1.2 Correia et al. 2019, npj Digital Medicine (city-wide EHR, known DDIs)

Correia RB, de Araújo Kohler L, Mattos MM, Rocha LM. "City-wide electronic health records reveal gender and age biases in administration of known drug–drug interactions." *npj Digital Medicine* 2:74, 2019. doi:10.1038/s41746-019-0141-x. arXiv:1803.03571. [verified 2026-10-27 (web search); arXiv ID confirmed 2026-10-28 in arXiv API results: https://arxiv.org/abs/1803.03571]

- **Data:** 18 months of public-health-system dispensing records, Blumenau, Brazil (~340k inhabitants).
- **Risk source:** known DDIs from DrugBank, with a subset linked to major adverse reactions.
- **How combined:** co-administration network of drug pairs dispensed to the same patient in overlapping time windows; "interaction risk" = probability that a patient is co-administered a known interacting pair, stratified by sex and age, with a null model controlling for number of drugs.
- **Findings:** about 12% of patients were co-administered a known interacting pair; women had ~60% higher risk than men (~90% for major interactions); patients aged 70–79 dispensed 2+ drugs had 34% risk; polypharmacy alone did not explain the sex/age gaps.
- **Relevance:** closest methodological template for RQ4 (`1.5`) and for "prevalence of risky pairs" weighting. Known interactions only.

### 1.3 Three-population follow-up, BMC Medicine 2024 (EHR, known DDIs)

Sánchez-Valle J, Correia RB, Camacho-Artacho M, Lepore R, Mattos MM, Rocha LM, Valencia A. "Prevalence and differences in the co-administration of drugs known to interact: an analysis of three distinct and large populations." *BMC Medicine* 22:166, 2024. doi:10.1186/s12916-024-03384-1. PMID 38637816. Sánchez-Valle and Correia contributed equally. Update of medRxiv 10.1101/2023.02.06.23285566. [verified: https://pubmed.ncbi.nlm.nih.gov/38637816/ and Crossref https://api.crossref.org/works/10.1186/s12916-024-03384-1]

*Correction:* the first pass gave no author list and called this "the Rocha group line of work". The first author is Sánchez-Valle (Barcelona Supercomputing Center); it is a joint BSC (Valencia) and Rocha / Correia paper.

- **Data:** population-wide EHR drug administrations in Blumenau, Brazil (133K individuals), Catalonia, Spain (5.5M) and Indianapolis, USA (264K). [verified: PubMed abstract]
- **Risk source:** DrugBank (Blumenau drugs matched to DrugBank IDs; Catalonia ATC codes mapped to DrugBank). [verified: full text PMC11027217]
- **How combined:** DDI prevalence stratified by severity, sex and age at the time of administration; null models for the expected effect of polypharmacy; simulated swaps to alternative drugs. [verified: PubMed abstract]
- **Findings:** co-administration of known interacting drugs affects 12.51% (Blumenau), 12.12% (Indianapolis) and 10.06% (Catalonia) of individuals; swapping omeprazole for other proton pump inhibitors would cut the share exposed by up to 21% (Blumenau, Catalonia) and 2% (Indianapolis). [verified: PubMed abstract]
- **Relevance:** cross-population comparison figures for the "share exposed" headline number.

### 1.4 Machine-learning medication recommendation on MIMIC-III (EHR, known TWOSIDES DDI graph)

GAMENet: Shang J, Xiao C, Ma T, Li H, Sun J. *AAAI 2019*. arXiv:1809.01852. SafeDrug: Yang C, Xiao C, Ma F, et al. *IJCAI 2021*. arXiv:2105.02711.

*The first pass had these from prior knowledge only. They are now checked against the full texts and the claim holds, with details added.*

- **GAMENet:** "we used DDI knowledge from TWOSIDES dataset", keeping the top-40 severity DDI types and mapping NDC codes to ATC level 3 to match MIMIC-III. 123 of 145 medications are in the DDI graph. "DDI rate" = share of drug pairs within a visit's medication set that are edges of the TWOSIDES graph. Base DDI rate in the EHR test data: **0.0777**. The DDI graph is a known-edge adjacency matrix, not a learned prediction. [verified: https://arxiv.org/abs/1809.01852, full text v3, pp. 5–7]
- **SafeDrug:** "We extract the DDI information by Top-40 severity types from TWOSIDES", reported at ATC level 3. Ground-truth DDI rate on MIMIC-III: **0.0808** ("around 8%"). [verified: https://arxiv.org/abs/2105.02711, full text v2, Table 2 and Appendix B]
- **Relevance:** the only line found that puts TWOSIDES polypharmacy labels next to real co-prescriptions at scale. About 8% of real ICU within-visit drug-class pairs are TWOSIDES-flagged. Known TWOSIDES pairs, not predictions; ATC-3 classes, not drugs; ICU, not a population; goal is recommendation, not burden. Worth citing if MIMIC-IV is used (`2.5`).

### 1.5 Kee et al. 2024, Heliyon (learned polypharmacy side-effect model, single-hospital EHR) — the closest ML paper

Kee G, Kang HJ, Ahn I, et al. (14 authors; last author Kim YH). "Are polypharmacy side effects predicted by public data still valid in real-world data?" *Heliyon* e24620, 2024. doi:10.1016/j.heliyon.2024.e24620. PMID 38304832. [verified: full text https://pmc.ncbi.nlm.nih.gov/articles/PMC10831713/]

- **Model:** DeSIDE-DDI, trained on 4,576,287 drug–drug–side-effect triples (63,472 drug pairs, 963 side effects) from BioSNAP, "refined from TWOSIDES". The paper mentions Decagon as related work but does not use it.
- **Data:** Asan Medical Center (Seoul) EHR, 2000–2020, adults.
- **How combined:** sort predicted drug–drug–side-effect triples that are not already known by score; for the top ones, count patients at Asan who used both drugs; with clinicians, pick a pair that had both a high score and enough co-users; run an IPTW cohort study for that one pair.
- **Findings:** cefpodoxime + chlorpheniramine users (n = 1,628) had a higher 1-year risk of lung edema than cefpodoxime-only users (HR 2.10, 95% CI 1.26–3.52).
- **Difference from this project:** co-use count is used to decide whether a prediction is *testable*, not to weight burden. Single hospital, no survey weights, no calibration, one pair validated, no ranking of pairs by population burden. It is, however, direct precedent for "take a TWOSIDES-trained model, look up how many real patients take the pair" and for the validation idea in `7.1`.

### 1.6 Hecker, Frahm & Zettl 2023, Pharmaceutics (learned DDI model on real medication plans)

Hecker M, Frahm N, Zettl UK. "Update and Application of a Deep Learning Model for the Prediction of Interactions between Drugs Used by Patients with Multiple Sclerosis." *Pharmaceutics* 16(1):3, 2023. doi:10.3390/pharmaceutics16010003. PMID 38276481. [verified: abstract via Europe PMC]

- **Model:** DeepDDI (multi-label classifier of DDI types) retrained on >1 million DrugBank DDIs; DrugBank interaction types, not TWOSIDES side effects.
- **Data:** medication plans of 627 patients with MS, 312 small-molecule drugs.
- **Findings:** 3,748 DDIs found in DrugBank vs 13,365 predicted by DeepDDI in the medication plans; ≥1 predicted DDI in 509 patients (81.2%).
- **Difference from this project:** the share of a cohort with ≥1 predicted DDI is reported, which is the closest existing "prevalence of predicted DDIs" number. One disease cohort, not population-representative, no survey weights, binary predicted labels, no per-side-effect burden ranking.

### 1.7 Dagenais et al. 2024, Clin Transl Sci (claims prevalence × known PK risk class, "predicted DDI impact score")

Dagenais S, Lee C, Cronenberger C, Wang E, Sahasrabudhe V. "Proposing a framework to quantify the potential impact of pharmacokinetic drug-drug interactions caused by a new drug candidate by using real world data about the target patient population." *Clin Transl Sci* 17(3):e13741, 2024. doi:10.1111/cts.13741. PMID 38445532. [verified: full text https://pmc.ncbi.nlm.nih.gov/articles/PMC10915735/]

- **Data:** Optum Market Clarity (US real-world data), 14,163,271 adults; target populations T2DM, obesity, migraine.
- **How combined:** points for the prevalence of use of CYP3A / CYP2D6 / CYP2C9 victim substrates (classes from the UW Drug Interaction Database), plus points for risk factors among those exposed (age, polypharmacy, duration, number of prescribers), summed into a "predicted DDI impact score".
- **Difference from this project:** it is a prevalence-weighted population-impact score, the closest conceptual analogue of `7.5`, but the risk side is a known PK category scored in points, not a learned or calibrated probability, and it scores perpetrator mechanisms, not drug pairs or side effects.

### 1.8 Other learned-model papers that touch real co-prescription (minor)

- **DCSE, Jimenez & Paccanaro 2026** (*PLoS Comput Biol*, doi:10.1371/journal.pcbi.1013619): new polypharmacy side-effect model with prospective (time-split) TWOSIDES evaluation. One case study (fluoxetine + naratriptan) cites a registry study's triptan–antidepressant co-prescription rate (21–29% per year) to argue the pair matters. Anecdotal, one pair, not a weighting. [verified: full text PMC13657403]
- **Schwarz et al. 2022** (arXiv:2207.08665): UK Biobank prescription clusters and the number / types of *expected* DDIs per cluster, from DrugBank (via Therapeutics Data Commons). Known DDIs, not predictions. [verified: https://arxiv.org/abs/2207.08665, full text v2]

## 2. Prevalence of known risky pairs in NHANES / MEPS / registries

All of these combine prevalence with a **known** interaction source. None uses predicted risk.

| Study | Data | Risk source | Headline finding |
|---|---|---|---|
| Sutherland et al. 2015, *PLoS One*, PMID 25739022 | NHANES 1999–2010 | Drugs.com, CYP450 | See 1.1 [verified 2026-10-27] |
| Igho-Osagie et al. 2023, *Clin Ther* 45(5):390–399, PMID 37032225 | NHANES 2015–Mar 2020, 15,685 adults | Liverpool COVID-19 checker, Lexicomp, FDA fact sheets | 29.3% (survey-weighted) of US adults would have a major-to-contraindicated pDDI if given ritonavir; 60.2% at age ≥60 [verified 2026-10-27] |
| Meng et al. 2025, *J Am Acad Child Adolesc Psychiatry*, PMID 41423041 | NHANES 2001–Mar 2020, ages 6–24 | Micromedex | 26.0% of psychotropic users exposed to a contraindicated or major DDI [verified 2026-10-27] |
| Qato et al. 2008, *JAMA* 300(24):2867–78, PMID 19109115 | NSHAP 2005–06, 3,005 community-dwelling adults aged 57–85, in-home medication logs | Thomson Micromedex, applied to the 20 most common prescription / OTC drugs and 20 most common supplements | 4% ("nearly 1 in 25") potentially at risk of a major DDI; half involved nonprescription products; 10% of men aged 75–85; no contraindicated combinations [verified: https://pubmed.ncbi.nlm.nih.gov/19109115/ and full text PMC2702513] |
| Qato et al. 2016, *JAMA Intern Med* 176(4):473–82, PMID 26998708 | NSHAP longitudinal cohort aged 62–85, 2,351 (2005–06) and 2,206 (2010–11), direct medication inspection | Micromedex | Share at risk of a potential major DDI rose from **8.4%** (2005–06) to **15.1%** (2010–11) [verified: https://pubmed.ncbi.nlm.nih.gov/26998708/] |
| Correia et al. 2019 and Sánchez-Valle et al. 2024 | Blumenau, Indianapolis, Catalonia EHR | DrugBank | 10–12.5% of patients co-administered a known interacting pair [verified: see 1.2, 1.3] |
| Aarabi et al. 2026, *PLoS One* 21(4):e0345076, PMID 41950170 | 2,365,811 outpatient prescriptions (982,102 patients), Tehran, Iran | Micromedex, top 100 most-prescribed drugs | Potential DDIs in 46.1% of prescriptions and 57.8% of patients; network with edges weighted by co-prescription frequency and coloured by severity; most common major pairs: aspirin + metformin (2.36% of interactions), + clopidogrel (1.46%), + hydrochlorothiazide (1.36%) [verified: abstract and full text PMC13061183] |

*Corrections to the first pass:*
- **Qato row was wrong.** It said the share at risk "rose from about 1 in 25 (2005–06) to about 15% (2010–11)". The 1-in-25 (4%) figure is from the 2008 paper (ages 57–85, Micromedex on the top-20 lists only). The 2016 paper re-estimated 2005–06 at **8.4%** in its own cohort (ages 62–85) and compared that with 15.1%. The two 2005–06 figures are not comparable, so 4% → 15% overstates the rise.
- **Qato risk source was wrong.** It said "Micromedex / Lexi-Interact". Both papers used Micromedex; no Lexi-Interact was found in either.
- **PLOS One 2026 row:** the country (not re-checked before) is Iran (Tehran); the interaction checker is Micromedex; the aspirin pairs are confirmed in the full text.

MEPS: a PubMed search for MEPS + drug interaction returned only 3 hits, none of which estimates pair-level DDI prevalence (Murphy et al. 2018 *Cancer*, PMID 29645083, describes polypharmacy patterns in cancer survivors; Vyas et al. 2012 *Int J Family Med*, PMID 22900173, links chronic-condition clusters to polypharmacy). [verified 2026-10-27] Europe PMC full-text search for "polypharmacy side effects" AND (NHANES OR MEPS) also found nothing relevant (2026-10-28). MEPS appears rarely used for pair-level DDI exposure, which fits `4.1`'s concern that MEPS annual purchase records do not show concurrent use.

### 2.1 Nordic and UK prescription-registry studies

*The first pass listed these from background knowledge only and gave no figures. All rows below were read on 2026-10-28.*

| Study | Data | Risk source | Headline finding |
|---|---|---|---|
| Johnell & Klarin 2007, *Drug Saf* 30(10):911–8, PMID 17867728 | Swedish Prescribed Drug Register, Oct–Dec 2005, 630,743 people aged ≥75 with ≥2 drugs | Type C / D interaction classes (database not named in the abstract) | Type C potential DDIs 26%, type D 5%; strong association with number of drugs; type D less likely with older age and in women [verified: https://pubmed.ncbi.nlm.nih.gov/17867728/] |
| Holm et al. 2014, *Eur J Clin Pharmacol* 70(11):1375–83, PMID 25190295 | Swedish national register of dispensed drugs, Jan–Apr 2010, whole population | SFINX | 38% dispensed C-interactions, 3.8% D-interactions; **a top-15 list of D-interactions covers 80% of all interacting combinations; 10 drugs are involved in 94% of D-interactions** [verified: https://pubmed.ncbi.nlm.nih.gov/25190295/] |
| Holm et al. 2019, *PLoS One* 14(8):e0220685, PMID 31381591 | Same register, children 0–17 | SFINX | Among children with ≥2 drugs: D-interactions 0.14%, C-interactions 1.3% [verified: https://pubmed.ncbi.nlm.nih.gov/31381591/] |
| Åstrand et al. 2007, *Eur J Clin Pharmacol* 63(9):851–9, PMID 17598096 | Jämtland county (Sweden) cohort, ~8,000 people, 1983–84, 1993–94, 2003–04 | Computerised interaction system (C / D classes) | Relative risk of C-interactions rose and of D-interactions fell over 20 years; few combinations account for much of the risk [verified: https://pubmed.ncbi.nlm.nih.gov/17598096/] |
| Guthrie et al. 2015, *BMC Med* 13:74, PMID 25889849 | All ~310,000 adults in Tayside, Scotland, community dispensing, 1995 and 2010 | BNF (March 2010) "potentially serious" interactions | Share with a potentially serious DDI more than doubled to 13% of adults in 2010; 10.9% at 2–4 drugs vs 80.8% at ≥15 drugs [verified: https://pubmed.ncbi.nlm.nih.gov/25889849/, DDI definition from full text PMC4417329] |

- **Not found in this search:** a whole-population DDI-prevalence study from the Danish, Norwegian or Finnish prescription registries (only drug-specific studies, e.g. warfarin or metoprolol, turned up), and a CPRD / THIN study of population-wide DDI prevalence across all drug pairs (only pair-specific cohort studies, e.g. DOAC + statin). This means "not found in this search", not "does not exist".
- **Relevance:** Holm et al. 2014 is the strongest support for a burden-ranking claim: a few pairs make up most of the population exposure. This project should test whether the same concentration holds for *predicted* burden. Guthrie 2015 and Johnell & Klarin 2007 give the "rises steeply with drug count" curve to compare against for RQ4 (`1.5`).

## 3. How risk and prevalence are combined in prior work

| Approach | Used by | Comparable to this project's |
|---|---|---|
| Share of people exposed to ≥1 known DDI (binary per person, weighted) | Igho-Osagie 2023, Meng 2025, Qato 2008/2016, Correia 2019, Sánchez-Valle 2024, Johnell & Klarin 2007, Holm 2014, Guthrie 2015 | Patient-level risk (`7.6`), but with binary known labels instead of calibrated probabilities |
| Count of alerts per person vs. number of drugs | Sutherland 2015, Guthrie 2015 | RQ4 by number of concurrent drugs (`1.5`) |
| Rank pairs by co-prescription frequency, then label by severity | Sutherland 2015 (pairs ≥0.1%), Holm 2014 (top-15 D-interactions), Aarabi 2026 network | Burden ranking (`7.5`) with severity (`10.2`), but no probability term |
| Network of co-administered pairs, interaction risk by sex and age with a null model | Correia 2019, Sánchez-Valle 2024 | RQ4 stratification |
| Share of within-visit drug pairs that are TWOSIDES-flagged (DDI rate) | GAMENet 2019, SafeDrug 2021 (MIMIC-III) | Coverage / overlap of real co-medication with TWOSIDES (`4.3`), but flagged share, not coverage |
| Share of a cohort with ≥1 **predicted** DDI (binary, learned model) | Hecker 2023 (DeepDDI, 627 MS patients) | Patient-level predicted risk (`7.6`), but one clinical cohort, binary labels, DDI types not side effects |
| Rank **predicted** triples by score, keep those with enough local co-users, validate one | Kee 2024 (DeSIDE-DDI, one hospital EHR) | Validation of top predictions (`7.1`); co-use as a filter, not a weight |
| Prevalence of exposure × known PK risk class, summed into an impact score | Dagenais 2024 (Optum claims) | Population-impact score (`7.5`), but no learned or calibrated probability |
| prevalence × **predicted, calibrated** per-side-effect probability, summed over side effects | none found | `7.5` itself |

## 4. What this means for `1.1`

- **New:** using a learned, calibrated polypharmacy side-effect model on nationally representative, survey-weighted co-medication; burden broken down by side effect; the share of top-burden pairs that no known source lists; the coverage of NHANES co-medication by Decagon's 645 drugs.
- **Not new:** "X% of the population takes a known interacting pair", its rise with age and drug count, and its concentration in a few pairs (Holm 2014). Also not new: applying a learned DDI model to real patients' medication lists (Hecker 2023) or checking real-world co-use of top predictions (Kee 2024). The project should compare its known-label numbers with these papers rather than present them as findings, and should state its difference from Kee and Hecker explicitly.
- **Comparison targets:**
  1. Sutherland 2015: same NHANES source; reproduce "pairs ≥0.1% prevalence" and check how many are Decagon-labelled.
  2. Igho-Osagie 2023 and Meng 2025: survey-weighted NHANES DDI prevalence methods (cycle pooling, weights), useful for `10.1`.
  3. Correia 2019 / Sánchez-Valle 2024: sex and age stratification method and the 10–12.5% known-DDI exposure figures.
  4. Holm 2014: check whether predicted burden is as concentrated as known-DDI exposure (top-15 pairs = 80%).
  5. Hecker 2023: their 81.2% "≥1 predicted DDI" vs this project's patient-level predicted risk; explains why a binary predicted label inflates prevalence and why calibration matters.
  6. GAMENet / SafeDrug: their ~8% TWOSIDES-flagged share of MIMIC-III pairs vs this project's NHANES figure.
- **Sanity checks (`7.1`):** the Aarabi 2026 most common major pairs (aspirin + metformin, aspirin + clopidogrel, aspirin + hydrochlorothiazide) are candidates for "should rank high" if those drugs are in Decagon. Kee 2024's cefpodoxime + chlorpheniramine → lung edema is a prediction with real-world support and could serve as a positive control if both drugs are in Decagon.

## Search log

### 2026-10-27 (first pass)

PubMed searches for NHANES + DDI (10 hits, all reviewed by title), NHANES + co-prescription terms (30 hits), MEPS + DDI (3 hits); abstracts read for Sutherland 2015, Igho-Osagie 2023, Meng 2025; web search for predicted DDI × co-prescription prevalence and for Correia 2019. Exact query strings were not recorded in that session.

### 2026-10-28 (follow-up)

**OpenAlex, works citing Zitnik et al. 2018 (Decagon, OpenAlex W2786016794, 1,499 citers)**, `https://api.openalex.org/works?filter=cites:W2786016794&search=<q>`:

| `search=` | Hits | Relevant results read |
|---|---|---|
| `NHANES` | 2 | none relevant (a 2025 review, a mental-health network paper) |
| `MEPS OR "Medical Expenditure Panel Survey"` | 1 | none relevant |
| `co-prescription OR coprescription OR "co-prescribed" OR "concomitant"` | 109 | Hecker 2023, DCSE 2026, Saha et al. 2024 (abstract not reachable, see *Blockers*) |
| `claims` | 216 | Neural bandits (Larouche et al., simulated claims data, arXiv:2212.05190), none relevant |
| `prescription registry OR "prescription register"` | 116 | none beyond the above |
| `prevalence population burden` | 37 | Zhou et al. 2025 *PLoS One* (MDFLDRR model; predicts cardiovascular × antidepressant DDIs, no patient data) |
| `"electronic health records" prevalence drug pairs` | 26 | **Kee et al. 2024** |
| `"real-world" prescriptions polypharmacy side effects` | 61 | none beyond the above |

**OpenAlex, Decagon citers filtered by venue** (`primary_location.source.id`), then keyword filter on title + abstract (EHR, claims, NHANES, MEPS, prescription, co-prescri, real-world, registry, cohort, medication list, MIMIC, prevalence): JAMIA 1 citer / 0 hits; JAMIA Open 1 / 0; Bioinformatics 37 / 7; BMC Bioinformatics 25 / 4; J Biomed Inform 7 / 1; IEEE JBHI 11 / 1; KDD 9 / 4; NeurIPS 0. All hits were "real-world" used loosely in benchmark papers, or EHR used for representation learning; none applies a polypharmacy model to population co-medication. These citer counts are far below the true numbers for those venues, so OpenAlex venue coverage is incomplete (especially conferences: NeurIPS 0, AMIA proceedings not indexed as a source).

**PubMed E-utilities** (`esearch.fcgi?db=pubmed&term=…`):

| Query | Hits | Relevant |
|---|---|---|
| `(decagon OR "polypharmacy side effect*" OR "polypharmacy side-effect*") AND ("electronic health record*" OR EHR OR claims OR NHANES OR MEPS OR "real-world" OR prescription* OR registry)` | 7 | Kee 2024, DCSE 2026 |
| `("drug-drug interaction*"[tiab] AND (predict*[ti] OR "deep learning"[tiab] OR "machine learning"[tiab] OR "graph neural"[tiab])) AND (prevalence[tiab] OR "co-prescri*"[tiab] OR coprescri*[tiab] OR "medication list*"[tiab] OR "medication plan*"[tiab]) AND (NHANES[tiab] OR "national survey"[tiab] OR claims[tiab] OR registry[tiab] OR "electronic health record*"[tiab] OR "real-world"[tiab])` | 5 | Sutherland 2015, Aarabi 2026, Lalitkanjanakul 2026 (RA, Micromedex/WebMD vs co-prescription, LLMs; known DDIs) |
| `(DeepDDI OR "DeSIDE-DDI" OR decagon[tiab] OR TWOSIDES[tiab]) AND (patient*[tiab] OR cohort[tiab] OR population[tiab])` | 28 | Hecker 2023 |
| `("predicted drug-drug interaction*"[tiab] OR "predicted interactions"[tiab] OR "predicted DDI*"[tiab]) AND (prevalence[tiab] OR "co-prescri*"[tiab] OR exposure[tiab] OR "population"[tiab])` | 83 | Dagenais 2024 (rest mostly PBPK and network-pharmacology papers) |
| `Hecker M[au] AND (DeepDDI OR "drug-drug interaction*" OR "deep learning")` | 8 | only Hecker 2023 uses a learned model; the others use DDI databases |
| `("deep learning"[tiab] OR "machine learning"[tiab] OR "neural network"[tiab]) AND ("drug-drug interaction*"[tiab] OR "polypharmacy"[tiab]) AND ("medication plans"[tiab] OR "medication lists"[tiab] OR "prescription data"[tiab] OR "dispensing data"[tiab] OR "claims data"[tiab] OR NHANES[tiab]) AND (predict*[tiab])` | 6 | Hecker 2023 only |
| `("drug interaction*"[tiab] OR "drug-drug interaction*"[tiab] OR DDI[tiab] OR DDIs[tiab]) AND ("Prescribed Drug Register"[tiab] OR "Swedish Prescribed Drug"[tiab] OR "Norwegian Prescription Database"[tiab] OR "Danish National Prescription Registry"[tiab] OR "Danish National Prescription Register"[tiab] OR "Finnish Prescription Register"[tiab] OR "prescription register*"[tiab] OR "prescription registr*"[tiab]) AND (prevalence[tiab] OR exposure[tiab] OR frequency[tiab] OR "potential"[tiab])` | 17 | Johnell & Klarin 2007, Åstrand 2007 |
| `("drug-drug interaction*"[ti] OR "drug interaction*"[ti] OR "interacting drug*"[ti]) AND (Sweden[tiab] OR Swedish[tiab] OR Norway[tiab] OR Norwegian[tiab] OR Denmark[tiab] OR Danish[tiab] OR Finland[tiab] OR Finnish[tiab]) AND (register*[tiab] OR registry[tiab] OR nationwide[tiab]) AND (prevalence[tiab] OR "potential"[tiab])` | 17 | Holm 2014, Holm 2019 |
| `("drug interaction*"[tiab] OR "drug-drug interaction*"[tiab] OR DDIs[tiab]) AND (prevalence[tiab] OR "potential"[tiab]) AND (CPRD[tiab] OR "Clinical Practice Research Datalink"[tiab] OR "The Health Improvement Network"[tiab] OR QResearch[tiab] OR Scotland[tiab] OR Scottish[tiab] OR Tayside[tiab] OR "UK Biobank"[tiab] OR "Wales"[tiab] OR "England"[tiab]) AND (prescri*[tiab] OR dispens*[tiab] OR "primary care"[tiab])` | 34 | pair-specific UK cohort studies only |
| `(polypharmacy[ti] OR "drug interaction*"[ti] OR "drug-drug interaction*"[ti]) AND (Scotland[tiab] OR Tayside[tiab] OR UK[tiab] OR "United Kingdom"[tiab]) AND (population[tiab] OR "primary care"[tiab]) AND (prevalence[tiab] OR trend*[tiab])` | 17 | polypharmacy-only prevalence (e.g. Payne 2014, PMID 24487416); no DDI figures |
| `"rising tide of polypharmacy"` | 3 | Guthrie 2015 (searched by title because the two queries above missed it) |
| `Qato DM[au] AND 2016[dp] AND JAMA Intern Med[ta]` | 2 | Qato 2016 |

**Europe PMC full-text search** (`webservices/rest/search?query=…`):
- `"polypharmacy side effects" AND (NHANES OR "Medical Expenditure Panel Survey")`: 2 hits, none relevant.
- `(decagon AND polypharmacy) AND ("co-prescription" OR "co-prescribed" OR coprescription)`: 5 hits. Relevant: DCSE 2026 (case study, see 1.8); Fast polypharmacy side-effect prediction using tensor factorisation (*Bioinformatics* 2024, doi:10.1093/bioinformatics/btae706), which mentions co-prescription only as a limitation (models handle pairs, patients take 3+ drugs).
- `("predicted" AND "polypharmacy side effect") AND ("survey weights" OR "survey-weighted")`: 0 hits.

**arXiv API** (`export.arxiv.org/api/query?search_query=…`):
- `abs:"polypharmacy" AND (abs:"electronic health records" OR abs:claims OR abs:NHANES OR abs:"prescription" OR abs:"co-prescription")`: 9 hits. Relevant: Schwarz et al. 2022 (UK Biobank, see 1.8), Correia et al. 2018 preprint.
- `abs:decagon AND abs:polypharmacy`: 2 hits (Decagon; DrIVeNN 2023, whose abstract describes a cardiovascular "treatment cohort" case study with drug-pair AUROC; whether patient medication data are used is not clear from the abstract).
- `abs:"drug-drug interaction" AND (abs:prevalence OR abs:NHANES OR abs:"co-prescribed" OR abs:"co-prescription")`: 7 hits, none relevant.

**Web search** (two "extended" searches and one "standard" search on predicted DDI × co-prescription / survey prevalence / ML venues): only Decagon, Sutherland 2015, Hecker 2023 and this repo's own issue #10 came back. The issue is not a source.

**Prior-knowledge checks:** GAMENet and SafeDrug full texts (arXiv PDFs, 0.46 MB and 1.2 MB); Qato 2008 (PubMed + PMC full text), Qato 2016 (PubMed); BMC Med 2024 (PubMed, Crossref, PMC full text).

## Blockers

1. **Semantic Scholar API: rate limited.** `curl "https://api.semanticscholar.org/graph/v1/paper/search?query=polypharmacy%20side%20effect%20prediction%20co-prescription%20prevalence&limit=30&fields=title,year,venue,externalIds"`: HTTP 429 on 4 attempts (back-off 8, 16, 24 s). A second query (`predicted%20polypharmacy%20side%20effects%20real-world%20prescriptions%20prevalence`) about 20 minutes later: HTTP 429 on 3 attempts 30 s apart. Body: `{"message": "Too Many Requests. Please wait and try again or apply for a key for higher rate limits. https://www.semanticscholar.org/product/api#api-key-form", "code": "429"}`. Partly replaced by OpenAlex citers, PubMed and Europe PMC.
2. **Google Scholar: no results returned.** WebFetch `https://scholar.google.com/scholar?q=%22polypharmacy+side+effect%22+decagon+co-prescription+prevalence` returned a page saying no matching articles and "Das System kann den Vorgang jetzt nicht ausführen. Versuchen Sie es später erneut" ("The system cannot perform the operation right now. Please try again later"). Not retried. Google Scholar was not searched.
3. **Saha et al. 2024, "Adverse Event Classification from Co-prescribed Drugs by Integrating Chemical, Phenotypic and Graph Embedding Features"** (LNCS, PReMI, doi:10.1007/978-3-031-12700-7_35): no abstract in OpenAlex, Crossref or Europe PMC. WebFetch of `https://link.springer.com/chapter/10.1007/978-3-031-12700-7_35` returned `303 See Other` to `https://idp.springer.com/authorize?response_type=cookie&client_id=springerlink&…`. A curl retry with a cookie jar was refused by the session's permission classifier ("Bot Check Circumvention"). Not read. From the title it may use co-prescription data, so it is an open lead. [unverified]
4. **Europe PMC full-text XML for Qato 2008:** `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC2702513/fullTextXML` returned HTTP 500 (`"Internal Server Error"`). Worked around via NCBI `efetch.fcgi?db=pmc&id=2702513`, so it has no effect on the findings.
5. **Venue coverage gap (not an error):** OpenAlex lists 0 NeurIPS and 9 KDD citers of Decagon, and AMIA proceedings are not indexed as a source. NeurIPS, KDD and AMIA proceedings were therefore not searched completely. AMIA Joint Summits papers in PubMed were covered by the PubMed queries (e.g. PMID 32477666, "Predicting Polypharmacy Side-effects Using Knowledge Graph Embeddings", came up; title only, abstract not read).

## Limits of this search

- Not done: Semantic Scholar and Google Scholar (see *Blockers*); full-venue browsing of NeurIPS / KDD / AMIA proceedings; DrugBank / TWOSIDES coverage studies beyond the GAMENet / SafeDrug DDI rate; Danish, Norwegian and Finnish population-wide studies (none found by the queries above).
- Abstract-only reads (not full text): Hecker 2023, Qato 2016, Johnell & Klarin 2007, Holm 2014, Holm 2019, Åstrand 2007, Lalitkanjanakul 2026, Zhou 2025, DrIVeNN 2023. The first-pass rows marked [verified 2026-10-27] were not re-read.
- The "no direct precedent" finding means "none found in this search", not "none exists".
