# Spec: co-prescription burden (Track 2, RQ2–RQ4)

Build-ready spec for the pipeline that weights the model's calibrated `p_r(a, b)` by how many people in the US take each drug pair, ranks pairs, side effects and patient groups by burden, and checks how much of that burden is documented. Every decision below comes from a closed ticket of the [wayfinder map #1](https://github.com/MichaelRos1/sief-druco/issues/1); the ticket is named in each section and holds the reasoning. Background: `PROJECT_CONTEXT.md`. Vocabulary: `GLOSSARY.md`. Model side: `docs/spec-model.md`.

The pipeline treats the model as a black box: it reads only the score table (§2). Pair-level only; drug sets (3+ drugs) are out of scope, but every table keeps each person's full drug list so they can be added later (#21).

---

## 1. What the pipeline delivers

| Research question | Output | Section |
|---|---|---|
| RQ2 (exposure) | Pair prevalence for every covered pair; coverage of real co-medication by the 645 drugs | §5, §6, §7 |
| RQ3 (burden) | Rankings of pairs, triples (pair × side effect) and side effects by burden; the share of burden that is undocumented; a novelty candidate list | §8, §10 |
| RQ4 (patient groups) | Burden by age band and by drug-count band, with each band's top pairs | §9 |
| Validation | Positive controls, top-50 stability across re-runs (MEPS, FAERS-era, per cycle) | §11 |

All population figures are a **share of the whole pooled population** (all ages, non-users included), written **per 100,000 people** (#18). Top items are translated into "≈ N million US residents" only in the text. One exception: band prevalence and band burden in §9 use the band's own population.

Decisions marked "spec audit" were made after the tickets closed, when both specs were checked against the tickets; they are recorded on map #1.

## 2. Input from Track 1 (#17)

Read from the git-ignored `local/` folder, never recomputed here:

- `scores.npy`: float32 `[207690, 963]`, cross-fitted, Platt-calibrated `p_r(a, b)`.
- `pairs.csv`: `pair_id` (row of `scores.npy`), `stitch_a`, `stitch_b` (unordered, a < b by index), `labelled`, `regime` (U / S).
- `labels.npz`: sparse 0/1 `[207690, 963]`, same order; the Decagon training labels, used only for "known to Decagon" per triple (§10).
- `side_effects.csv`: `cui` (column of `scores.npy`), `name`, number of labelled pairs, held-out AUPRC per regime, predicted vs observed positives.
- `meta.json`: git commit, model setup, folds, seed, date.

`p_r(a, b)` estimates P(pair reported with r in TWOSIDES), not P(pair causes r) (#7, #9). There is no per-score uncertainty and no live scoring function.

## 3. Co-prescription data: NHANES (#11, #15)

- **Files:** `RXQ_RX` (one row per drug; `SEQN`, `RXDUSE`, `RXDDRUG`, `RXDDRGID`, ...) and `DEMO` (`SEQN`, `RIDAGEYR`, `RIAGENDR`, weights, `SDMVSTRA`, `SDMVPSU`) per cycle. Non-users have one row with `RXDUSE` = 2.
- **Main pool:** NHANES **2015–2016 (`_I`) + 2017–March 2020 pre-pandemic (`P_`)**: 5.2 years, ~10,900 Rx users, middle year ~2017.
  - **Never** 2017–2018 (`_J`) next to `P_`: `P_` already holds those participants under new `SEQN`s.
  - Out: 1999–2002 (retrofitted Multum codes, 4-year weights) and 2021–2023 (no drug names).
- **Weights:** interview weights only: `WTINT2YR × 2/5.2` for 2015–16, `WTINTPRP × 3.2/5.2` for `P_`. Raw weight sums are never used as population estimates.
- **Variance:** Taylor linearization on stacked `SDMVSTRA` / `SDMVPSU` (strata 119–133 and 149–172 don't overlap; no single-PSU stratum). Subgroups (age bands, drug-count bands) are survey **domains**, never subsets of the data.
- **Proportion intervals:** Korn–Graubard 95% (#18).
- **Software (spec audit):** a small module of our own in numpy / scipy: Taylor variance for ratio and mean estimators from stratum and PSU totals, domains, and Korn–Graubard as a Clopper–Pearson interval at the effective sample size. Checked once against R's `survey` package (`svyciprop(method = "beta")`, `svymean`) on a few 2017–18 estimates; R is not part of the pipeline.
- **Order of work:** the pipeline takes a list of cycles. Develop on 2017–2018 alone (already mapped in the pilot #16), then switch to the main pool.
- **Supplement:** the top pairs also get separate 2015–16 and `P_` shares; the pooled estimate is the result (#15).

## 4. Drug mapping (#12, #16, #26)

All sources meet the model at the **RxNorm ingredient (`IN`)**, on the reverse route.

1. **Decagon table (once):** each of the 645 STITCH IDs → PubChem CID (strip `CID` and leading zeros) → UNII → RxNorm `IN`. Most IDs are flat CIDs without a UNII, so the parent / same-connectivity CID and name-match fallbacks are the normal path. A bare UNII listed by PubChem is accepted only if PUG-View confirms it; CIDs without PubChem relations get a SMILES-based same-connectivity step (#16). The automatic table has `stitch_id`, `cid`, `rxcui`, `name`, `match_method`; each of the 645 is matched or listed unmatched with a reason (#12). With the overrides applied, the result is STITCH ID → **one main `IN` + hand-checked aliases**, each alias tagged with its relation.
2. **Source → `IN`:** NHANES `RXDDRGID` (Multum) → RxNorm → `IN` (PIN → `IN`; multi-ingredient products split into their `IN`s). MEPS `RXDRGNAM` and DDInter names go → `IN` by name (MEPS combinations split at `-`).
3. **Join** source `IN`s against the Decagon table.

**Rules (#26):**
- **Alias rule:** a Decagon drug joins an `IN` other than its own only when they are the same active drug: acid / salt / anion forms, a prodrug ester and its active drug, or a racemate and its enantiomer. Epimers and diastereomers sold as separate drugs (quinine / quinidine, dexamethasone / betamethasone) never join. NHANES / MEPS "tenofovir" joins Decagon's tenofovir disoproxil, footnoted (sources don't separate TDF from TAF).
- **One `IN`, several STITCH IDs:** the `IN` points to the STITCH ID with the **most labelled pairs** (hydrocortisone → `CID000003640`, 396 pairs, not the butyrate ester `CID000003642`, 26). Collisions are footnoted.
- **Hand fixes, Decagon side only, about half a day:** the 34 unmatched Decagon drugs, the 17 that reach more than one `IN`, and the hits of the exact / loose name checks on the main pool and MEPS. The **name check** compares the most prevalent unmatched source ingredients by name (exact, then loose) against the 645 Decagon names, to catch missed mappings (pilot: 2 exact, 8 loose hits; `research/nhanes-coverage-pilot.md` §7). Unclear cases stay unmatched with a reason. Fixes live in a committed override file `data/decagon_overrides.csv` (`stitch_id`, `main_in`, `aliases`, `verdict`, `reason`, `evidence`), holding IDs, names and verdicts only, never survey rows. `aliases` are each tagged acid/salt, prodrug or racemate-enantiomer; `verdict` is matched / dropped `IN` / unmatched; `evidence` is SIDER name, PubChem formula / InChIKey match, RxNav status. A script builds the evidence sheet, Claude proposes verdicts, the user reviews the *unclear* and *dropped* rows.
- **No source-side hand fixes:** source codes and names that fail stay `mapping_failed`, logged with their cause (no code `55555` / `77777` / `99999`, `a` / `c` / `h` code, retired code, name not found).
- **Two failure reasons, never merged:** **not in the model** (reached an `IN` outside the 645; reported only) and **mapping failed** (never reached an `IN`, or the Decagon drug never reached RxNorm). For DDInter, a failed name makes the pair *not checkable* instead (§10).
- Cache every RxNav / PubChem lookup to disk.

Pilot (NHANES 2017–18, #16): 99.4% of mentions reach an `IN`; 611 of 645 Decagon drugs reach RxNorm; 69.2% of weighted co-prescribed person-pairs have both drugs in the 645, limited by drugs outside Decagon (ethinyl estradiol, montelukast, escitalopram, rosuvastatin, ...), not by mapping (≈ 1% failures).

## 5. Data layer (#17)

Two long tables, the same for NHANES and MEPS:

- **persons:** `source` (NHANES / MEPS), `cycle`, `seqn`, pooled weight, stratum, PSU, age, sex. Key (`source`, `cycle`, `seqn`).
- **person_drugs:** one row per person × RxNorm `IN`, with `stitch_id` or the reason `not_in_model` / `mapping_failed`.

Pairs (and later drug sets) are **derived on the fly** from `person_drugs`, never stored as counts:
- all unordered pairs of distinct `IN`s within a person; covered drugs are de-duplicated **by STITCH ID** first, so two aliases of one Decagon drug (e.g. methylphenidate and dexmethylphenidate) never form a pair (spec audit);
- fixed-dose combinations are split into ingredients and the ingredient pair **counts as co-prescribed**;
- the same ingredient twice in one person counts once; a drug paired with itself is dropped;
- a pair is **covered** (scorable) if both drugs have a `stitch_id`; otherwise it keeps its reason(s). A pair with one failed and one not-in-model drug counts as *mapping failed*.

Sex stays in `persons` but is not used in this spec (#25).

## 6. Coverage reporting (`4.3`, #26)

- **Main table:** three columns adding to 100%: *covered* / *not in the model* / *mapping failed*. Rows: distinct ingredients, mentions, participants with ≥ 2 drugs, person-pairs; each unweighted and weighted. Point estimates (coverage is descriptive).
- **Supplement:** *mapping failed* split by side (source / Decagon) and cause; one line on how many mentions the hand fixes and aliases moved. All numbers are after the fixes.

## 7. Prevalence (#15, #18)

- **Prevalence of a pair** = Σ pooled weights of persons taking both ÷ Σ all pooled weights, per 100,000. Korn–Graubard 95% interval from the Taylor variance.
- **Headline population:** all ages (#18). Children rarely take two prescription drugs, so the ranking barely changes; adults appear through the RQ4 age bands.
- **Minimum count for ranking:** a pair, and its triples, is ranked only if **≥ 30 unweighted participants** in the pool take it (in the spirit of NCHS's presentation standard for proportions, Parker et al. 2017; check its wording before citing it in the paper). Every covered pair still counts in coverage, base rates, side-effect totals and burden shares.
- **Most common pairs overall (RQ2, spec audit):** one table of the 20 most prevalent pairs among **all** co-prescribed pairs, covered or not, each with its status (covered / not in the model / mapping failed). Everything else in this spec uses covered pairs only.
- The pipeline reports how many covered pairs pass 10, 20 and 30 participants. NHANES 2017–18 alone: 8,884 covered pairs; ≥ 10: 542 (41% of covered pair weight); ≥ 30: 114 (21%). The pool has ~2.75× the participants; rough scaling gives ~500 ranked pairs.

## 8. Burden (#18)

```
burden_r(a, b) = prev(a, b) × p_r(a, b)            triple (pair × side effect r)
burden(a, b)   = prev(a, b) × Σ_r p_r(a, b)        pair, all 963 side effects
burden_r       = Σ_pairs prev(a, b) × p_r(a, b)    side effect, over all covered pairs
```

- **Meaning:** a pair's burden is the expected number of side-effect signals per person in the population who takes the pair. It measures exposure to a signal, not cases of harm.
- **Equal weights:** all 963 side effects count the same. No probability threshold, no subset, no disease classes (`effectcategories.csv` covers only 399 of 963). Severity weights are a later extension (#20).
- **Labelled pairs:** `p` is the cross-fitted prediction for **every** triple, including a labelled pair's labelled side effects (S regime). Labels never replace predictions, so labelled and unlabelled pairs compete on the same footing.
- **Three rankings:** pairs and triples (only pairs with ≥ 30 participants, §7), side effects (all covered pairs).
- **Uncertainty:** the prevalence interval is carried into burden with `p` held fixed. Model uncertainty is shown only through per-side-effect columns from `side_effects.csv` (held-out AUPRC, predicted vs observed positives) and stated as a limitation. No intervals on ranks.

## 9. Patient groups (RQ4, #25)

- **Pair level only, main NHANES pool only** (the re-runs in §11 don't repeat RQ4). A person's burden = Σ over their covered pairs of Σ_r `p_r` (expected signals across their pairs; a drug counts again in each pair it is part of). No per-patient risk: aggregating a person's pairs (noisy-OR / max) is drug-set work (#21).
- **Age bands:** 0–17, 18–44, 45–64, 65+ (`RIDAGEYR`, top-coded at 80).
- **Drug-count bands:** 2, 3–4, 5–9, 10+. The count is a person's distinct `IN`s, including those outside the 645, **plus one per drug that never reached an `IN`** (uncoded row, class code), as in the pilot. Not `RXDCOUNT`, which counts products, so a combination pill is one drug there but two ingredients and one pair here (spec audit). People with 0–1 drugs appear only as a population share.
- No sex split; no age × drug-count cross-table.
- **Per band, both band types:**
  - (a) the band's **share of total burden** = the band's weighted sum of person burdens ÷ the same sum over everyone;
  - (b) **burden per person** in the band (a domain mean), per 100,000 band members, with a 95% Taylor interval (`p` fixed), next to the mean number of covered pairs per person (n drugs give n(n−1)/2 pairs, so burden rises steeply with drug count partly by construction);
  - (c) the band's **top 10 pairs** by band burden = band prevalence × Σ_r `p_r`, only pairs with ≥ 30 participants **within the band**. Fewer than 10: shown and marked "fewer than 10 rankable pairs"; none: "too few participants to rank".
- **Band denominators (spec audit):** band prevalence and band burden are shares of the **band's own population** (e.g. per 100,000 people aged 65+); the ranking within a band is the same either way, and whole-population shares would make small bands look tiny. (a) stays relative to everyone.
- **Band list columns:** pair; band prevalence per 100,000 band members with Korn–Graubard interval; band burden; overall rank (§8); known to Decagon; clinically documented (§10). No per-band novelty tables.
- NHANES 2017–18: people aged 45+ hold 92% of covered person-pairs, people on 5+ drugs 89%. The under-45 bands and the 2-drug band will have only a handful of rankable pairs.

## 10. Novelty: what is undocumented (#19, #18, #13)

Every covered pair and triple gets **two labels, never merged**:

- **Known to Decagon:** pair = has any label in the 963-side-effect training table; triple = that triple is a label.
- **Clinically documented:** pair = **DDInter 1.0** (all 14 ATC CSVs, de-duplicated, highest level kept, names → `IN` by §4's rules) lists it as Major, Moderate or Minor. "Unknown" (description incomplete) counts as not documented; its share goes in one footnote. Triples inherit the pair's status, marked "pair-level only".
- Each label is documented / not documented / **not checkable** (a drug is not in the source, or its name doesn't reach an `IN`).
- **Novelty candidate** = not known to Decagon and not clinically documented.

**Reporting:**
- **Headline ("how much is undocumented"):** the 2×2 as **shares of total burden** over all covered pairs and over all covered triples, with no count floor.
- **Top-k tables:** the 2×2 for the **top 50 pairs** and the **top 50 triples** of the ranked set (≥ 30 participants), read against the base rate over all covered pairs and triples.
- **Novelty candidate list:** candidates with **≥ 10 participants** in the pool, ordered by pair burden; those under 30 flagged as imprecise. NHANES 2017–18: only 3 of 682 candidate pairs reach 10 participants, so the ≥ 30 rule would leave almost none.
- **DrugBank hand check:** the top 20 of the candidate list (all of it if shorter), looked up by hand in DrugBank's interaction checker once, at the end. Record listed yes/no + severity, copy no text, cite DrugBank 6.0 (read as within DrugBank's terms for non-clinical research; not legal advice).

**Not used:** TWOSIDES / nSIDES (FAERS-derived, circular with the labels); CRESCENDDI (4.2% of covered pair weight, not a complete list); Drugs.com (terms forbid scraping and need written consent for research); the prospective FAERS check (cut on #19). DDInter is CC BY-NC-SA 4.0: its files stay in `local/`, only derived labels are reported.

NHANES 2017–18, covered pair weight: DDInter can check 99.5%; graded-documented 31.9% (in Decagon) and 1.0% (not in Decagon); "Unknown" 60.3%; novelty pool 3.4%.

## 11. Validation (#27, #22, #15)

Positive controls test the model; re-runs on other exposure data test the exposure measure.

### 11.1 Positive controls

- **Source:** AGS 2023 Beers Criteria, Table 5 (drug–drug interactions to avoid in older adults; [PMC12478568](https://pmc.ncbi.nlm.nih.gov/articles/PMC12478568/)), fixed before any results are seen. No hand-picked extras. 12 entries; "≥ 3 CNS-active drugs" is a drug set (out of scope); "anticholinergic + anticholinergic" is expanded from Beers Table 7.
- **Conditions ignored (spec audit):** two entries apply only to subgroups (RAS inhibitor + potassium-sparing diuretic: chronic kidney disease stage 3a or worse; alpha-1 blocker + loop diuretic: older women). The controls test the model's `p` for the combination, not patient eligibility, so the conditions are ignored and footnoted.
- **Before any score is looked up**, write down: each class expanded into drugs in the 645, and each outcome mapped to Decagon side effects by name (e.g. hyperkalaemia → "Excess potassium"; bleeding → "Bleeding", "gastrointestinal bleed", "cerebral haemorrhage"; respiratory depression → "Hypoventilation", "respiratory failure").
- **Model view:** for each drug pair and mapped side effect r, the percentile of `p_r(a, b)` among **all 207,690 pairs'** scores for r (a model question, so it doesn't depend on who takes which pair; spec audit). **Recovered** = top 10% for at least one mapped side effect.
- **Burden view:** the pair's overall rank and whether it reaches the top 50 (≥ 30 participants only).
- **Report:** "k of n recovered", each control's percentile and rank. Controls labelled in Decagon are marked (S-regime scores are partly expected to be high); controls with no matching side effect are listed as not checkable. No overall pass/fail.
- NHANES 2017–18: only opioid + gabapentinoid is common enough for the burden view (gabapentin + tramadol 17, gabapentin + oxycodone 14 participants); lithium and hydrocodone didn't match the 645 by name and are worth a look during the hand fixes (§4). Phansalkar et al. 2012 was rejected: its interactions are contraindicated, so almost never co-taken.

### 11.2 Re-runs and the stability measure

Every re-run is compared with the main result by **overlap of the top 50 pairs**, with the pairs that drop out named. Pairs only; no pass/fail threshold; no rank correlations.

- **Forming each top 50 (spec audit):** the same rule as the main result, applied within each dataset: ≥ 30 unweighted participants (MEPS: person-years), then the 50 highest pair burdens. Overlap = shared pairs out of 50. Each main-pool top-50 pair that drops out gets one reason: not observed in the re-run, below 30 participants, or ranked lower.

| Re-run | Tests | Setup |
|---|---|---|
| 2015–16 vs `P_` shares | cycle-to-cycle | §3 supplement |
| **FAERS-era pool** (#22) | time gap between labels and exposure | NHANES 2005–2006 + 2007–2008 (`_D`, `_E`), `WTINT2YR × 1/2`, Taylor on stacked strata (44–74, no overlap). 249 new codes (3.1% of mentions) go through §4 unchanged (estriol is the class code `a70248` there, so it stays `mapping_failed`). Also report its coverage. |
| **MEPS, light** (#27, #23) | exposure measured by purchases, not 30-day use | below |

**MEPS (light):**
- Files: PMED HC-178A (2015), HC-188A (2016), HC-197A (2017); **2015–2017 only**, before the spring 2018 questionnaire redesign. ≈ 43,000 person-years with ≥ 2 ingredients.
- **Co-prescribed = both drugs bought in the same round** (`PURCHRD`; three rounds per calendar year, about 4 months each), counted once per person-year.
- Weights `PERWT15F / 3`, `PERWT16F / 3`, `PERWT17F / 3`. A person present in two years counts as two person-years. **Point estimates only:** MEPS prevalence is used only for the top-50 overlap and is never reported as a finding; no HC-036 linkage, no BRR.
- Mapping by `RXDRGNAM` (93% of records reach an `IN`) under §4's rules; masked drugs and unmapped records get `mapping_failed`; coverage reported as in §6. Same `persons` / `person_drugs` tables (`source` = MEPS).
- Footnote: the documentation counts 4,226–5,496 fewer positive-weight records per year in 2015–2017 than the released files; the released weights are used.

## 12. Limitations to state

- **Coverage:** drugs outside the 645 can't be scored (transductive model); coverage is reported, not hidden.
- **Labels:** `p` is "reported together in FAERS", not causation; unlabelled means "not reported", not "safe".
- **Exposure:** co-listing means co-medication, not that a side effect occurred. Burden is predicted exposure to a signal, not cases.
- **Model uncertainty:** burden intervals cover prevalence only.
- **Time gap:** TWOSIDES covers AERS Q1 2004 – Q1 2009; the pool centres on ~2017. Measured by the FAERS-era re-run.
- **Sums count drugs repeatedly:** a person's burden adds up pairs that share a drug; it is not a patient risk.
- **No causal claims:** results are prioritised hypotheses.
- **Novelty claim:** no precedent was found for learned, calibrated risk × survey-weighted prevalence (#10; confidence medium–high). Before writing "first to", run the Google Scholar citing-articles search listed in `research/prior-coprescription-burden.md`. Cite the known-risk × prevalence work (Sutherland 2015 on NHANES, Holm 2014, Guthrie 2015) and the closest ML work: Kee 2024 (co-use as a filter), Hecker 2023 (share of a cohort with ≥ 1 predicted DDI), Dagenais 2024 (prevalence × a known CYP class), Patrick 2021 (claims data as a model input).

## 13. Not in this spec

- **The model and the score table:** `docs/spec-model.md`.
- **Drug sets of 3+ drugs** (noisy-OR / max, per-patient risk): #21. **Severity weights:** #20. **HGT, TIP:** #24.
- **True higher-order (3+ drug) models:** Decagon has no 3+ drug data.
- **Application layer** (`PROJECT_CONTEXT.md` `9`): belongs with drug sets.
- **MIMIC-IV:** different population, needs PhysioNet credentialing.
- **Inductive embedding of drugs outside the 645** (`10.5`): out of scope; coverage is reported instead.
