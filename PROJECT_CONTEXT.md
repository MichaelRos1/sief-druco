# Project Context: Polypharmacy Side-Effect Model × Real-World Co-Prescription Burden

> **For the next model reading this:** This file holds everything decided so far. Use it to map out how to build the project. Every point has a stable ID (e.g. `3.2`, `5.4`); the user refers to points by these IDs.

---

## 0. Project summary

| | |
|---|---|
| **Idea** | Build a model that predicts adverse side effects of drug combinations, then combine its predictions with real-world co-prescription data to rank which **commonly co-prescribed** drug pairs carry the highest predicted side-effect burden, and how much of that burden is not yet documented. |
| **Research question** | See `1`. |
| **Model data** | Decagon dataset (Stanford SNAP / BioSNAP), see `3`. |
| **Co-prescription data** | NHANES Prescription Medications (primary), see `4`. |
| **Core tech** | PyTorch + PyTorch Geometric. |
| **Status** | Planning stage. No code written yet. |

---

## 1. Research question

- **1.1 Main question:** Which drug combinations that are commonly taken together in the population carry the highest predicted polypharmacy side-effect burden, and how much of that burden is not yet documented?
- **1.2 RQ1 (model):** How well does a graph neural network trained on Decagon predict polypharmacy side effects under pair-level and drug-level splits, compared with simple baselines (additive per-drug model, tensor factorisation)?
- **1.3 RQ2 (exposure):** Which drug pairs and sets are taken together in the population, how often, and what share of real co-medication is covered by Decagon's 645 drugs?
- **1.4 RQ3 (burden):** When predicted risk is weighted by co-prescription prevalence, which pairs and which side effects account for most of the burden? How many top-burden pairs are neither known to Decagon (its labels) nor clinically documented (DDInter graded levels, DrugBank hand check)? Defined in `7` step 7 (#19).
- **1.5 RQ4 (patient groups):** How is predicted burden distributed across patient groups (age, number of concurrent drugs)?
- **1.6 Motivation:** Existing work (Decagon, `10.1`) asks "is this pair risky?" on a fixed benchmark and ignores which pairs people actually take. A risky pair that is almost never co-prescribed matters less for public health than a moderately risky pair taken by millions. Risk × prevalence gives a prioritised list for pharmacovigilance and prescribers.

---

## 2. Decisions made so far

| ID | Decision | Choice |
|---|---|---|
| **2.1** | Model dataset | Decagon (`3`) |
| **2.2** | Graph view | Nodes = drugs, proteins. Edges = protein–protein, drug–protein targets, drug–drug (one edge type per side effect). Single-drug side effects = drug node features. See `3.2`. |
| **2.3** | Pairs vs. sets | Build and evaluate the model on **pairs** first. Extend to N drugs by **pairwise aggregation** (max / noisy-OR, `6.7`). True higher-order (3+ drug) models are out of scope: Decagon has no 3+ drug data. |
| **2.4** | Research direction | Real-world co-prescription burden (`1`), not a pure reproduction of Decagon. |
| **2.5** | Co-prescription source | NHANES (primary), MEPS (robustness check), MIMIC-IV (optional extension). See `4`. |
| **2.6** | Architecture style | Modular: data → input layer → encoder → decoder → training → evaluation → inference. Modules talk only through fixed interfaces so the encoder can be swapped and the N-drug aggregator added without touching the rest (`5`). |
| **2.7** | Encoder strategy | RGCN is the **reference** encoder (it is what Decagon used), not assumed to be best. It is compared against no encoder (mono input only), collapsed GraphSAGE and collapsed GAT, plus tensor factorisation and per-drug additive baselines (`5.3`, `6.5`). Comparison procedure: #8. Answers RQ1. |

---

## 3. Model data: the Decagon dataset

### 3.1 Files and verified statistics

Source: Stanford SNAP (`snap.stanford.edu/decagon`) and BioSNAP. Column names are from memory and must be checked on download.

| File | One row means | Columns (verify) | Size (from SNAP pages) |
|---|---|---|---|
| `bio-decagon-combo.csv` | Drugs a + b together were reported with side effect r | `STITCH 1, STITCH 2, Polypharmacy Side Effect, Side Effect Name` | 4,649,441 triples; 645 drugs; 63,473 drug pairs |
| `bio-decagon-mono.csv` | Drug a alone causes side effect s | `STITCH, Individual Side Effect, Side Effect Name` | 174,977 rows; 640 drugs; 10,185 side effects |
| `bio-decagon-targets.csv` | Drug a binds protein p | `STITCH, Gene` | 18,690 edges; **only 284 drugs**; 3,648 proteins |
| `bio-decagon-targets-all.csv` | Same, merged from several curated databases | `STITCH, Gene` | 131,034 edges |
| `bio-decagon-ppi.csv` | Protein p interacts with protein q | `Gene 1, Gene 2` | 715,612 edges; 19,081 proteins |
| `bio-decagon-effectcategories.csv` | Side effect → disease class | `Side Effect, Side Effect Name, Disease Class` | small |

**IDs:** drugs = STITCH CIDs (e.g. `CID000002173`), proteins = NCBI Entrez Gene IDs, side effects = UMLS concept IDs (e.g. `C0151714`) plus a name.

### 3.2 The graph

| Element | Count | Meaning |
|---|---|---|
| `drug` nodes | 645 | drugs in the combo file |
| `protein` nodes | 19,081 | human proteins |
| protein – protein edges | 715,612 | protein interactions (PPI) |
| drug → protein edges (+ reverse) | 18,690 | drug binds target; 284 drugs have any |
| drug – drug edges | 964 edge types (one per side effect) | **the edges to predict** |
| single-drug side effects | 10,185-dim multi-hot per drug | drug node features (as in Decagon) |

- **Side effects are edge types, not nodes.**
- PyG layout:
  ```python
  data = HeteroData()
  data['drug'].x       = mono_multi_hot          # [645, 10185]
  data['protein'].x    = None                    # learned embedding [19081, d]
  data['protein', 'ppi', 'protein'].edge_index   # both directions
  data['drug', 'targets', 'protein'].edge_index  # + ('protein','rev_targets','drug')
  data['drug', 'combo', 'drug'].edge_index       # drug pairs, both directions
  data['drug', 'combo', 'drug'].edge_label       # [num_pairs, 964] multi-hot
  ```
- For RGCN with one relation per side effect: `data.to_homogeneous()` + `RGCNConv(..., num_relations=967, num_bases=B)` (964 side effects + PPI + targets + reverse targets).

### 3.3 Properties that shape the modelling

- **3.3.1 Side-effect filter:** the Decagon paper kept side effects occurring in at least 500 drug pairs, giving 964 types. The raw combo file has more types, so apply this filter to match the paper.
- **3.3.2 Dense drug–drug graph:** 63,473 of 207,690 possible pairs (about 31%) are labelled; clustering coefficient 0.81; diameter 3. Two layers of drug–drug message passing already reach almost all drugs, so **oversmoothing** is a real risk.
- **3.3.3 Labels per pair:** 73.25 side effects on average (median 53, standard deviation 66). Labels are about 7.6% ones within labelled pairs; rare side effects have few positives.
- **3.3.4 Target coverage:** 361 of 645 drugs have no target edges in `targets.csv`; they get no biological signal. Check how many the fuller `targets-all.csv` covers before relying on it.
- **3.3.5 Drug graph is small, protein graph is large:** the full graph fits on one GPU; the cost is in decoding millions of drug–drug triples plus negatives.

### 3.4 Data caveats

- **3.4.1 Noisy labels.** Decagon comes from TWOSIDES (Tatonetti et al. 2012), mined from FDA spontaneous reports (FAERS). A positive means "statistically over-reported together", not proven causation. Known biases: indication bias, notoriety bias, the Weber effect.
- **3.4.2 Missing ≠ safe.** Unlabelled pairs include real interactions nobody reported. This is **positive-unlabeled (PU)** data.
- **3.4.3 Single-drug leakage.** Many pair labels are really side effects of one of the two drugs. A model can score well by learning per-drug tendencies, so always compare against a per-drug additive baseline (`6.5`).

---

## 4. Co-prescription data

| Source | Role | Access |
|---|---|---|
| **NHANES Prescription Medications** (primary) | Per participant: all prescription drugs taken in the past 30 days, with demographics and survey weights | Public, no registration |
| MEPS Prescribed Medicines | Robustness check: prescribed medicines purchased during a year, with survey weights | Public |
| MIMIC-IV prescriptions | Optional extension: hospital inpatient co-medication | Credentialed via PhysioNet |
| DrugBank interactions | Reference for "already documented" pairs (RQ3) | Free for academic use |

- **4.1 Why NHANES first:** it records drugs taken *currently* (past 30 days), so co-listed drugs are very likely taken at the same time. MEPS records purchases over a year, where co-listing does not guarantee concurrent use. Survey weights turn sample counts into US population estimates; several survey cycles can be pooled.
- **4.2 Drug-name mapping (main data-engineering task):** reverse route, meeting at the RxNorm ingredient (details and open checks: issue #12). The model is transductive, so only Decagon's 645 drugs need the hard hops.
  1. **Decagon side (once):** STITCH ID → PubChem CID (strip the `CID` prefix and leading zeros: `CID000002173` → PubChem `2173`) → RxNorm ingredient. Save as a table; every one of the 645 is matched or listed with a reason.
  2. **NHANES side:** Multum code (`RXDDRGID`) → RxNorm → ingredient level (salt forms stripped, e.g. metoprolol succinate → metoprolol; combination products such as "amlodipine; benazepril" split into their ingredients).
  3. **Join** each NHANES ingredient against the Decagon table.
  4. **Two failure reasons, kept apart:** *not in the model* (reached an ingredient outside the 645; reported only) and *mapping failed* (never reached an ingredient, or the Decagon drug never reached RxNorm; fixed where possible).
  Cache every lookup to disk; log failures with their reason.
- **4.3 Coverage must be reported:** the share of participants, drugs and co-prescribed pairs that map to Decagon's 645 drugs. Drugs outside Decagon cannot be scored by a transductive model.
- **4.4 Pooled cycles (decided on #15):** main pool = NHANES 2015–2016 + 2017–March 2020 pre-pandemic file (`P_`), 5.2 years, ~10,900 Rx users, middle year ~2017.
  - **Never** 2017–2018 next to `P_` (`P_` contains those participants under new `SEQN`s, so a double count can't be caught by ID). 1999–2002 and 2021–2023 (no drug names) are out.
  - **Weights:** interview weights only, `WTINT2YR × 2/5.2` and `WTINTPRP × 3.2/5.2`. Prevalence is a weighted **share** of the population, not a sum of weights.
  - **Variance:** Taylor linearization on stacked `SDMVSTRA` / `SDMVPSU` (strata 119–133 and 149–172, no single-PSU stratum). The minimum-count rule for ranking a pair is set in #18.
  - **Population:** all ages stay in the pool; age bands are survey domains. The headline population is set in #18.
  - **Order of work:** the pipeline takes a list of cycles; develop on 2017–2018 (already mapped, #16), then switch to the main pool.
  - **Reporting:** pooled estimate is the result; top pairs also get separate 2015–16 and `P_` shares (supplement).
  - **Robustness:** FAERS-era pool 2005–2008 (#22, see `8.2`); MEPS 2015–2019 without 2020 (#23).

---

## 5. Model architecture

```
 1. DATA         CSVs → filter → index → HeteroData → split by drug pair
       ▼
 2. INPUT        drug: mono profile (10,185) → Linear → d
                 protein: Embedding(19,081, d)
       ▼
 3. ENCODER      2 × relational message passing over PPI, targets, drug–drug
                 → z_drug [645, d]
       ▼
 4. DECODER      (z_a, z_b, r) → logit,  e.g. z_aᵀ · D_r · R · D_r · z_b
       ▼
 5. TRAINING     BCE on positives + degree-matched negatives
 6. EVALUATION   AUROC / AUPRC / AP@50 per side effect, pair- and drug-level splits
       ▼
 7. INFERENCE    precompute z → score pairs → calibrate → aggregate → burden (`7`)
```

### 5.1 Data layer
1. Filter to 964 side effects (`3.3.1`).
2. Index mappings: STITCH → 0..644, Entrez → 0..19,080, UMLS → 0..963 (keep reverse mappings for names).
3. Drug features: [645, 10,185] multi-hot from `mono.csv`.
4. Build `HeteroData` (`3.2`).
5. **Split by drug pair** (all side effects of a pair in the same set), plus a **drug-level (cold-drug)** split.
6. Split training pairs into **message edges** (used in the encoder) and **supervision edges** (used only in the loss). Validation/test pairs never appear in the message graph, including reverse direction. PyG: `RandomLinkSplit(..., disjoint_train_ratio=...)`.

### 5.2 Input layer
- Drugs: `Linear(10185 → d)` on the mono profile (reduces 10k inputs; 10k features for 645 drugs overfits otherwise).
- Proteins: `Embedding(19081, d)`.

### 5.3 Encoder (compared, answers RQ1)

| Encoder | How it works | Why compare it |
|---|---|---|
| **None (no encoder)** | Mono input → `Linear` → decoder, no message passing | Same input as the GNNs, so it tests whether message passing helps at all. Tensor factorisation (free per-drug vectors, Malone et al. 2018) is a separate baseline (`6.5.2`). |
| **RGCN** (reference) | One weight matrix per edge type, built from ~30 shared bases (basis decomposition) | What Decagon used; comparable to published numbers |
| **Collapsed drug–drug edges + GraphSAGE / GAT** (`to_hetero`) | One drug–drug edge type for message passing; GAT adds attention over neighbours | Much cheaper; attention may help on the dense drug graph |
| ~~**HGT**~~ (future work, #24) | Attention with type-specific projections | Learns which relations and neighbours matter |
| ~~**TIP-style staged propagation**~~ (future work, #24) | Protein → protein, then protein → drug, then drug → drug | Reported to beat Decagon in accuracy and speed (Xu et al. 2019) |

- Known RGCN weaknesses here: equal weighting of neighbours, oversmoothing on the dense drug graph (`3.3.2`), high cost (millions of typed messages per layer).
- RGCN uses the dense drug–drug rewrite (30 shared bases), not `RGCNConv`: 11× faster per epoch, half the memory (#8 timing).
- Setup: 2 layers (fixed), ReLU, dropout and d tuned (`6.4`), degree normalisation.
- Output: z_drug [645, d]. Protein vectors are intermediate.

### 5.4 Decoder
- Default: Decagon's decoder `z_aᵀ · D_r · R · D_r · z_b + b_r`, D_r diagonal per side effect ([963, d]), R shared [d, d], constrained symmetric so f(a, b) = f(b, a). Computed as one matmul, `((z_a ⊗ z_b) ⊙ R).flatten() @ (D ⊗ D)ᵀ` (8× faster on MPS, #8).
- Alternatives: DistMult (always symmetric), multi-label MLP on `[z_a + z_b, z_a ⊙ z_b]` → 963 logits. Chosen in Stage A of the comparison (#8).
- **Per-side-effect bias `b_r`** on every decoder, started at side effect r's training base rate. Deviation from Decagon (no bias): with full rows ~98.6% of labels are 0, and without it the untrained model predicts 0.5 everywhere.
- Start values: D_r ~ N(0, 1), R ~ N(0, 1/√d). Smaller start values give near-zero gradients through the four-factor product (#8 pilot).
- Must stay a separate module taking only (z_a, z_b, r), so the N-drug aggregator reuses it unchanged.

### 5.5 Training
1. Encoder forward on the full graph → z for all drugs.
2. Supervision pairs, each with its **full row** of 963 labels (decided in #7). Every unlabelled triple is a plain 0, with no positive-unlabeled correction, so the model estimates P(pair reported with r).
3. Negatives: none sampled by default; the full rows hold every zero. Switchable fallback: zeros drawn uniformly per side effect at one rate β, independent of the positives. Decagon's degree-proportional corruption is not used: it can never draw 29% of the zeros and distorts each pair differently.
4. Decoder → logits; plain `BCEWithLogitsLoss` (no `pos_weight`, no focal), per-side-effect weights {none, `1/√n_pos(r)`} chosen once in Stage A and fixed for all models (#8).
5. Backpropagate; full-batch Adam; early stopping on validation macro AUPRC (patience 30, cap 500 epochs), test score at the best epoch.

---

## 6. Evaluation

- **6.1 Splits:** pair-level (main for the paper-style comparison) and drug-level / cold-drug (does the model generalise to drugs never seen in a pair). The Decagon paper used a random split by triple, which is optimistic because a pair's other side effects are in training.
- **6.2 Metrics:** AUROC, AUPRC, AP@50 per side effect, macro-averaged over 963 (Decagon reported AUROC 0.872, AUPRC 0.832, AP@50 0.803 on its random split, over 964). Also micro averages. AUPRC must be read against each side effect's prevalence. **Model selection uses macro AUPRC over the full validation rows of the pair-level split** (#8, ADR 0001); the other metrics are reported only.
- **6.3 Statistics:** 3 folds × 3 seeds over pairs/drugs, the same for every model (ADR 0001, #8); report mean ± spread. Pairs sharing a drug are not independent, so resample by drug.
- **6.4 Hyperparameters:** tune each encoder once, on pair-level fold 0 with one seed: learning rate {3e-3, 1e-2} × d {32, 64, 128} × dropout {0, 0.3}, 2 layers fixed; reuse the winning settings on the random and cold-drug splits (#8). Limitation: cold-drug may prefer more regularisation than pair-level picks.
- **6.4.1 Run order (#8):** Stage A picks the decoder and loss weights on pair-level (no encoder + RGCN, 1 fold × 3 seeds, default settings); then tuning; then Stage B compares all encoders and baselines on all three splits; then the no-protein ablation. ≈ 50 h of MacBook compute at an estimated ~300 epochs per run; Stage A measures the real count.
- **6.5 Baselines and ablations:**
  1. Per-drug additive model: σ(u_r(a) + u_r(b)). Two versions: a free per-drug table (random and pair-level) and from the mono input (all splits). Shows how much is per-drug tendency.
  2. Tensor factorisation: free per-drug vectors + the chosen decoder (random and pair-level; not applicable on cold-drug). Shows whether drug features beat drug IDs. Whether message passing helps is shown by encoder "none" (`5.3`).
  3. The best GNN from Stage B on drug–drug edges only (no proteins). Shows whether biology helps.
  4. Full graph; report separately for the 284 drugs with targets vs. the 361 without.
- **6.6 Calibration (required for burden):** measure with reliability diagrams / ECE; fix on validation data with temperature scaling, Platt scaling or isotonic regression. Burden multiplies probabilities by prevalence, so probabilities must be meaningful, not just correct rankings. Full-row training (#7) needs no subsampling correction; runs that use the uniform-sampling fallback at rate β are corrected with p = β·p_s / (β·p_s − p_s + 1) (Dal Pozzolo et al. 2015). PU bias is not fixed by either: `p` estimates the reported probability.
- **6.7 N-drug aggregation:** per side effect over all pairs of a drug set, max or noisy-OR: P = 1 − ∏(1 − p_ij). Noisy-OR assumes independent, calibrated pairs; pairs sharing a drug are correlated and count that drug's own tendency repeatedly. Report both, and remember which pair drives each risk.

---

## 7. Burden analysis pipeline (answers RQ2–RQ4)

1. **Train and evaluate the model** (`5`, `6`). → RQ1
2. **Calibrate** probabilities (`6.6`).
3. **Map NHANES drugs** to Decagon IDs (`4.2`); report coverage (`4.3`).
4. **Extract co-prescriptions:** for each participant, all pairs of mapped drugs. Weighted prevalence of a pair = sum of survey weights of participants taking both. → RQ2
5. **Compute burden:**
   ```
   burden_r(a, b) = prevalence(a, b) × p_r(a, b)                        per side effect r
   burden(a, b)   = Σ_r  prevalence(a, b) × p_r(a, b) [× severity_r]     per pair
   ```
   Rank pairs and side effects by burden. → RQ3
6. **Patient-level risk:** per participant, aggregate over all their pairs (`6.7`); summarise by age group and number of drugs using survey weights. → RQ4
7. **Novelty check (decided on #19):** every covered pair and pair × side effect gets two labels, never merged:
   - **Known to Decagon:** pair = has any label in the 963-side-effect table the burden model is trained on; triple = that triple is a label.
   - **Clinically documented:** pair = DDInter 1.0 (all 14 CSVs, de-duplicated, highest level kept, joined by name → RxNorm IN) lists it as Major, Moderate or Minor. "Unknown" (description incomplete) counts as not documented; its share goes in a footnote. Triples inherit the pair's status, marked "pair-level only"; no source documents a specific side effect.
   - Each source gives documented / not documented / **not checkable** (a drug is not in the source).
   - Top-burden items are reported in the 2×2 table; **novelty candidates** = not known to Decagon and not clinically documented. Shares in the top-k (cut-offs from #18) are read against the base rate over all covered pairs.
   - **DrugBank hand check:** the top ~20 novelty candidates are looked up by hand in DrugBank's interaction checker (listed yes/no + severity, no text copied, cited as DrugBank 6.0), once, at the end.
   - Not used: TWOSIDES / nSIDES (FAERS-derived, circular); CRESCENDDI (cut: documents 4.2% of covered pair weight and is not a complete list).
   - NHANES 2017–18, covered pair weight: DDInter can check 99.5%; graded-documented 31.9% (in Decagon) and 1.0% (not in Decagon); "Unknown" 60.3%; novelty pool (in neither) 3.4%.

### 7.1 Validation of the burden ranking
- **Sanity check:** well-known risky, common combinations (e.g. anticoagulant + antiplatelet → bleeding) should rank high.
- **External evidence** for undocumented top pairs: the DrugBank hand check (`7` step 7). The prospective FAERS / nSIDES check was cut on #19 (still FAERS, licence unclear, little value).
- **Robustness:** repeat with MEPS instead of NHANES, and with max instead of noisy-OR; report how stable the top of the ranking is.

---

## 8. Expected contributions and limitations

### 8.1 Contributions
1. A population-weighted ranking of polypharmacy risk: from "which pairs are risky" to "where the real-world risk is".
2. A measured coverage figure: how much of real co-medication a Decagon-trained model can score.
3. A list of high-burden combinations not documented in standard interaction sources, as hypotheses for pharmacovigilance.
4. An honest model comparison on pair-level and drug-level splits, with baselines.
5. An application: screen a real medication list (`9`).

### 8.2 Limitations (state up front)
- **Coverage:** drugs outside Decagon's 645 cannot be scored by a transductive model; coverage is reported, not hidden.
- **Labels:** the model learns "over-reported together in FAERS", not causation (`3.4.1`).
- **Exposure:** co-listing in a survey means co-medication, not that a side effect occurred; burden is predicted, not observed.
- **Time gap:** Decagon's labels come from FAERS reports filed Q1 2004 – Q1 2009 (TWOSIDES); the main NHANES pool centres on ~2017 (`4.4`). The gap mostly shows as lower coverage (newer drugs are outside the 645). Measured by re-running on a FAERS-era pool, NHANES 2005–2008, and reporting coverage and top-50 burden overlap with the main pool (#22).
- **No causal claims:** results are prioritised hypotheses.

---

## 9. Application layer

- **9.1 Inference:** precompute z_drug once; per request run only the decoder on all N·(N−1)/2 pairs (milliseconds), calibrate, aggregate (`6.7`).
- **9.2 Backend:** e.g. FastAPI `POST /predict {drugs: [...]}` → top-k side effects with probability and the responsible pair.
- **9.3 Frontend:** Streamlit or Gradio for a prototype.
- **9.4 Explanations:** which pair drives each risk; optionally which proteins (PyG `torch_geometric.explain`, GNNExplainer).
- **9.5 Disclaimer:** research / educational tool, not medical advice.

---

## 10. Open items

- **10.1** ~~Which NHANES cycles to pool.~~ Decided on #15, see `4.4`.
- **10.2** Whether to add severity weights to the burden score (e.g. Gottlieb et al. 2015 crowdsourced severity rankings).
- **10.3** Final encoder: decided by the comparison procedure locked in #8 (`5.3`, `6.4`, `6.5`).
- **10.4** Literature check: prior work combining co-prescription or prescription-survey data with predicted DDI risk.
- **10.5** Optional extension for coverage: embed drugs outside Decagon inductively (from single-drug side-effect profiles, known targets, or molecular structure via SMILES → molecular GNN, fused with the graph embedding).

---

## 11. Background knowledge the project relies on

Explained in depth in an earlier session; re-explain on request, grounded in Decagon.

- **11.1 Representation:** embeddings and embedding tables; node features vs. learned embeddings; heterogeneous graphs; reverse edges.
- **11.2 Message passing on Decagon:** one layer step by step; layers = hops; aggregation (sum / mean / max) and degree normalisation; oversmoothing.
- **11.3 RGCN:** per-relation weights; basis decomposition (at d = 64: 964 full matrices ≈ 3.9M parameters vs. ≈ 152k with 30 bases); block-diagonal decomposition.
- **11.4 Attention:** GAT, HGT.
- **11.5 Decoders as tensor factorisation:** dot product, DistMult, RESCAL, Decagon (DEDICOM-style), MLP head; symmetry; low rank.
- **11.6 Training:** BCE; why negatives are needed; negative sampling and degree bias; PU learning (Elkan & Noto 2008; reporting bias violates "selected completely at random", see Bekker & Davis 2019); class imbalance (`pos_weight`, focal loss); message vs. supervision edges.
- **11.7 Transductive vs. inductive.**
- **11.8 Evaluation:** train/val/test; split types and what each measures; AUROC, AUPRC, AP@k; macro vs. micro; variance and resampling by drug; baselines and ablations.
- **11.9 Calibration and aggregation:** reliability, ECE, temperature / Platt / isotonic; noisy-OR derivation and its independence assumption.
- **11.10 Higher-order (out of scope, for context):** DeepSets, Set Transformer, hypergraph convolution; need 3+ drug data, which Decagon lacks.

---

## 12. Reading list

- **12.1** Zitnik, Agrawal & Leskovec (2018), "Modeling polypharmacy side effects with graph convolutional networks" (Decagon). Dataset, RGCN encoder, decoder, evaluation protocol.
- **12.2** Tatonetti et al. (2012), "Data-driven prediction of drug effects and interactions". Origin of TWOSIDES / OFFSIDES; how the labels were made.
- **12.3** Schlichtkrull et al. (2018), "Modeling Relational Data with Graph Convolutional Networks". RGCN and basis decomposition.
- **12.4** Malone, García-Durán & Niepert (2018), "Knowledge Graph Completion to Predict Polypharmacy Side Effects". Knowledge-graph embedding beating Decagon.
- **12.5** Xu et al. (2019), TIP: "Tri-graph Information Propagation for Polypharmacy Side Effect Prediction". Staged propagation, faster and more accurate than Decagon.
- **12.6** Elkan & Noto (2008), "Learning classifiers from only positive and unlabeled data". PU learning.
- **12.7** Bekker & Davis (2019), "Beyond the Selected Completely At Random Assumption for Learning from Positive and Unlabeled Data". PU with biased reporting.
- **12.8** Dal Pozzolo et al. (2015), "Calibrating Probability with Undersampling for Unbalanced Classification".
- **12.9** Gottlieb et al. (2015), "Ranking adverse drug reactions with crowdsourcing". Severity weights.
- **12.10** Cheng, Kovács & Barabási (2019), "Network-based prediction of drug combinations" (Nat Commun). Why protein-network proximity matters.
- **12.11** NHANES Prescription Medications documentation (CDC); MEPS Prescribed Medicines documentation (AHRQ).
