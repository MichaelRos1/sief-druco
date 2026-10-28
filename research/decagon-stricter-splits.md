# Decagon-style results under pair-level and cold-drug splits

Ticket: #2 (map #1). Question: what numbers exist for Decagon-style models on the Decagon dataset under a **pair-level split** and a **drug-level / cold-drug split**, instead of Decagon's random per-triple split? How large is the drop, and does any work report a per-drug additive baseline (`6.5.1`)?

Researched 2026-10-27 from the papers themselves (journal/PMC/arXiv full text). Items marked *not verified* were not checked against the source in this session.

Follow-up 2026-10-28 (PARTIAL items): Malone et al. 2018 and TIP (Xu et al. 2019) checked against the arXiv PDFs, the TIP workshop PDF and the TIP code on GitHub. Claims from that pass are tagged **[verified: source]** or **[unverified]**. Remaining gaps are listed under **Blockers**.

## Answer in brief

- Only **one** paper found reports a true **pair-level** split on Decagon's data: **Dewulf, Stock & De Baets 2021**. It uses a kernel ridge model, not a GNN, and scores the **full label tensor** (every unlabelled triple counts as a negative, about 2% positives). So its AUPRC values cannot be compared with Decagon's balanced-negative AUPRC.
- **Cold-drug** numbers exist in three forms, and each defines "new drug" differently:
  - Dewulf 2021: one or both drugs unknown, by drug-level k-fold.
  - SimVec 2022: the 98 lowest-degree drugs, with all their triples held out.
  - EmerGNN 2023: a random drug split, run on a *different* TWOSIDES subset (604 drugs, 200 side effects).
- Under pair-level and cold-drug splits, AUROC drops modestly and AUPRC drops a lot (Dewulf). Knowledge-graph embedding models such as Decagon and TriVec fall to **AUROC 0.30–0.44** on SimVec's weak-node split, which is below chance.
- **No paper found reports a per-drug additive baseline σ(u_r(a)+u_r(b)).** The closest are:
  - Decagon's "concatenated drug features" baseline (random split only).
  - Malone et al. 2018's concatenated mono-side-effect + target indicator baseline, which *beats* Decagon on the random split (AUROC 0.896 vs 0.872).
  - Dewulf's trivial "any labelled pair gets every effect" baseline (AUROC_dd 0.85).
  - EmerGNN's MLP baseline, which beats Decagon in the one-new-drug setting.
- No paper re-runs **Decagon itself** under a pair-level split on its 645-drug / 964-effect data.
- **Malone et al. 2018 and TIP (Xu et al. 2019) report random per-triple splits only** (now verified, see below). Neither has a pair-level or cold-drug result. Both copy Decagon's 0.872 / 0.832 / 0.803 from the Decagon paper instead of re-running it on their own split.

## Reference point: the random split everyone compares against

### Decagon: Zitnik, Agrawal & Leskovec 2018, *Bioinformatics* 34(13):i457

Source: <https://academic.oup.com/bioinformatics/article/34/13/i457/5045770>

- **Split:** "For each side effect type, we use 80% of drug pairs to train a model, and 10% of drug pairs to select model parameters"; the remaining 10% are test. The split is made separately for each side effect, so one pair can be in training for side effect r1 and in test for side effect r2. This is the per-triple split that `6.1` calls optimistic.
- **Negatives:** in training, each positive (v_i, r, v_j) is paired with a corrupted (v_i, r, v_n), with v_n drawn from a sampling distribution P_r. The page does not describe how negatives are built for evaluation.
- **Metrics:** AUROC, AUPRC and AP@50, computed per side effect and averaged over 964 side effects (Table 2).

Table 2 (average over 964 side effects):

| Approach | AUROC | AUPRC | AP@50 |
|---|---|---|---|
| Decagon | 0.872 | 0.832 | 0.803 |
| RESCAL | 0.693 | 0.613 | 0.476 |
| DEDICOM | 0.705 | 0.637 | 0.567 |
| DeepWalk | 0.761 | 0.737 | 0.658 |
| Concatenated drug features | 0.793 | 0.764 | 0.712 |

- The page contains no analysis by drug degree or held-out drug.
- "Concatenated drug features" is a per-drug feature baseline, but it is not an additive per-side-effect model; its exact construction was *not verified* in this session.

### Tensor factorisation on the same random split

- **TriVec** (Nováček & Mohamed): AUROC 0.975, AUPRC 0.966, as quoted in SimVec's Table 1. The earlier draft called TriVec "the successor of Malone et al. 2018's work". Neither the Malone paper nor the TIP paper supports that link, and the TriVec paper was not read, so it is **[unverified]** and has been dropped.
- **Malone et al. 2018** and **TIP** now have their own subsections below.
- **Lloyd, Liu & Gaunt 2024**, *Bioinformatics* btae706 (arXiv:2404.11374, <https://arxiv.org/html/2404.11374v2>):
  - Split: Decagon's protocol, "randomly removing 10% of the edges belonging to each PSE" (458,061 held-out edges). No pair-level or drug-level split.
  - Evaluation negatives: not described. Training uses 1vsAll / KvsAll.
  - Metrics: AUROC, AUPRC and AP@50 per side effect, reported as the **median** over 963 side effects. Decagon reported a mean.
  - SimplE (self-loops variant): **0.978 / 0.971 / 1.000**. ComplEx: 0.970 / 0.963 / 0.978. DistMult: 0.946 / 0.930 / 0.941.
  - The authors say tensor factorisation methods "simply cannot make inferences about unseen nodes".
  - The paper's literature table lists random-split AUPRC on Decagon data of up to 0.998 (Carletti et al.).
### Malone, García-Durán & Niepert 2018, DILS 2018 (LNCS), pp. 144–149; arXiv:1810.09227

Source: <https://arxiv.org/pdf/1810.09227> (arXiv v1, 22 Oct 2018, 6-page "Poster" version). The published chapter (DOI 10.1007/978-3-030-06016-9_14) is behind a login, so it was not compared with the arXiv version (see Blockers).

- **Data:** the Decagon preprocessed dataset: 645 drugs, 19,089 proteins, 4,649,441 drug-drug interactions, **963** distinct polypharmacy side effects (Table 1, p. 2). [verified: arXiv PDF, p. 2, Table 1]
- **Model:** DistMult, and KBlrn (a product of experts over DistMult embeddings and rule-based "relational features" from shared or interacting drug targets). Training uses negative sampling. [verified: arXiv PDF, §3, pp. 2–3]
- **Split:** the paper says it follows Decagon's "common experimental design". It builds one fixed record set of positives plus negatives, then uses "stratified sampling to split the records in training, validation and testing sets". The split proportions are **not stated**. Table 2 is on "the pre-defined test set". This is a per-triple random split: there is **no pair-level or drug-level holdout**. [verified: arXiv PDF, §4 "Dataset construction", p. 3; Table 2 caption, p. 5]
- **Negatives (evaluation):** made "by randomly selecting a pair of drugs and a polypharmacy side effect which does not exist in the knowledge graph", with equal numbers of positives and negatives for each side effect. So the evaluation is balanced 1:1, random-pair negatives, and the no-skill AUPRC is about 0.5. [verified: arXiv PDF, §4, p. 3]
- **Metrics:** AuROC, AuPR, and AP@50 ("average precision for the top 50 predictions for each polypharmacy side effect"). The paper does not say how per-side-effect values are combined (mean or median). [verified: Table 2 caption, p. 5; aggregation: unverified]
- **Decagon row:** "values reported in [14]". It is copied from the Decagon paper, not re-run on Malone's split. [verified: Table 2, p. 5]

Table 2 (p. 5), in the paper's column order:

| Method | AuROC | AuPR | AP@50 |
|---|---|---|---|
| Baseline (concatenated drug indicator vectors) | 0.896 | 0.859 | 0.812 |
| Decagon (copied from [14]) | 0.872 | 0.832 | 0.803 |
| DistMult | 0.923 | 0.898 | 0.899 |
| KBlrn | 0.899 | 0.878 | 0.857 |
| DistMult (drug-drug interactions only) | 0.931 | 0.909 | 0.919 |
| KBlrn (drug-drug interactions only) | 0.894 | 0.886 | 0.892 |
| DistMult (drugs with protein targets only) | 0.534 | 0.545 | 0.394 |
| KBlrn (drugs with protein targets only) | 0.829 | 0.797 | 0.774 |

[verified: arXiv PDF, Table 2, p. 5]

- **Beyond the random split:** none. The two extra settings are **graph ablations, not stricter splits**:
  - "Drug-drug interactions only" removes the protein and target edges.
  - "Drugs with protein targets only" removes every drug without an annotated target, together with its triples.
  - Both keep the same kind of record-level split. [verified: arXiv PDF, §4, p. 4]
- **Baseline construction:** each drug becomes "a binary vector of indicators for each mono side effect and gene target". The pair's two vectors are concatenated, and the model predicts each side effect from the concatenation. The classifier is not named. [verified: arXiv PDF, §4 "Baselines", p. 3] This per-drug-feature baseline **beats the Decagon numbers it quotes** on all three metrics. It differs from Decagon's own "concatenated drug features" baseline (0.793 / 0.764 / 0.712), and the two were run on different splits and negatives.
- **Code:** the paper links no repository. A GitHub search (`gh search repos kblrn`, plus "KBlrn polypharmacy") returned no results. [verified: arXiv PDF; gh search, 2026-10-28]

### TIP: Xu, Sang & Lu 2019, NeurIPS 2019 Workshop on Graph Representation Learning; arXiv:2001.10516

Sources: <https://arxiv.org/pdf/2001.10516> (arXiv v1, 28 Jan 2020). The workshop PDF <https://grlearning.github.io/papers/94.pdf> has the same text and tables. Code: <https://github.com/NYXFLOWER/TIP>, read at HEAD (`cfe3df4`, 2022-12-02) and at the paper-era commit `82e63fa` (2020-01-07).

**What the paper says**

- **Data:** the BioSNAP-Decagon files, with 645 drugs. The D-D graph has 63,473 edges and 1,317 labels before filtering. Side effects with fewer than 500 occurrences are removed, described as "the same pre-processing as in Zitnik et al.". [verified: arXiv PDF, Table 1, p. 2; Supp. §8, p. 7]
- **Split:** "For each side effect, we use 80% of the total edges in D-D graph for model training and the remaining 20% for testing." This is a per-triple random split. **No validation set** is described. [verified: arXiv PDF, Supp. §8, p. 7]
- **Negatives:** "the negative samples E′n are sampled randomly from R until E′n∩E′p = ∅". The paper does not separate training negatives from evaluation negatives. [verified: arXiv PDF, Supp. §8, p. 7]
- **Metrics:** AUPRC, AUROC, and "AP@k: average precision for the top k predictions for each side effect". [verified: arXiv PDF, Supp. §8, p. 7]
- **Decagon row:** "we quote the accuracy score in [16]". It is copied, not re-run. [verified: Table 2 caption, p. 4]

Table 2 (p. 4), SNAP-Decagon dataset:

| Model | ARCT | AUPRC | AUROC | AP@50 |
|---|---|---|---|---|
| Decagon (copied from [16]) | n/a | 0.832 | 0.872 | 0.803 |
| DistMult | DF | 0.835 | 0.859 | 0.834 |
| R-GCN | DDM-DF | 0.882 | 0.908 | 0.883 |
| dTIP_D | DDM-NN | 0.791 | 0.847 | 0.792 |
| dTIP_P | PPM-GGM-NN | 0.746 | 0.743 | 0.733 |
| TIP-cat | PPM-GGM-DDM-DF | 0.889 | 0.913 | 0.890 |
| TIP-sum | PPM-GGM-DDM-DF | **0.890** | **0.914** | **0.890** |

[verified: arXiv PDF and workshop PDF, Table 2, p. 4]

- **Beyond the random split:** none. The only extra analyses are model ablations (dTIP_D, dTIP_P) and the best and worst side effects by AUPRC (Supp. Figs. 4–5). [verified: arXiv PDF, pp. 4, 8]
- **"7.2+%" improvement over Decagon in AUPRC:** by my arithmetic, (0.890 − 0.832) / 0.832 = 7.0%, and the absolute gain is 0.058. The paper's 7.2% could not be reproduced from Table 2.
- **The DistMult row is not Malone's number.** TIP labels DistMult "(reported by [9])", where [9] is Malone et al. 2018. TIP's row (0.835 / 0.859 / 0.834) differs from Malone's DistMult (AuPR 0.898 / AuROC 0.923 / AP@50 0.899). The repository has a commit "add distmult model evaluation" (`5d41e5b`, 2020-01-06), which suggests the authors re-ran it, but this is **[unverified]**. Malone's reported DistMult AUPRC (0.898) is higher than TIP-sum's (0.890), but the two come from different splits and negatives.

**What the code actually does** (paper-era commit `82e63fa`; the split and negative-sampling functions are unchanged at HEAD)

- **Split:** `tip-cat.py` calls `process_edges(data.dd_edge_index)` with the default `p=0.9`. For each side effect, `process_edges` draws `np.random.binomial(1, p, n_edges)` and sends each (pair, side effect) edge to train or test. That gives about **90% train / 10% test, not the paper's 80 / 20**, and no validation set. One drug pair can be in training for one side effect and in test for another, so this is a per-triple split, as in Decagon. [verified: `src/utils.py` `process_edges` at `cc7f616` / `82e63fa` / HEAD; `tip-cat.py` line 35 at `82e63fa`; HEAD `src/layers.py` `Setting(sp_rate=0.9)`]
- **Negatives:** `typed_negative_sampling` works one side effect at a time. It draws as many ordered drug pairs, uniformly from all N_d², as there are positives, and rejects only pairs found among that side effect's positives **in the same slice**. [verified: `src/neg_sampling.py`, unchanged since `cc7f616`]
  - Test negatives are drawn once, before training. They exclude only that side effect's *test* positives, so a test "negative" can be a training positive for the same side effect.
  - Self-pairs (i, i) are not excluded.
  - Positives are stored in both directions (`to_bidirection`), so each test pair is scored twice.
  - Evaluation is therefore balanced 1:1 with random-pair negatives (no-skill AUPRC about 0.5).
- **Metrics:** for each side effect, `auprc_auroc_ap` returns AUPRC (trapezoidal area under the precision-recall curve), AUROC, and `sklearn.metrics.average_precision_score` over **all** test positives and negatives. The script prints that last value as "ap@50". It is **not a top-50 measure**, so the code does not match the paper's "AP@k … top k" definition. Results are averaged as `record.sum(axis=1) / n_dd_et`, a plain mean over side effects. [verified: `src/utils.py` `auprc_auroc_ap` and `tip-cat.py` lines 233–275 at `82e63fa`] Whether Table 2 was produced by exactly this code is **[unverified]**. The README says the default settings reproduce the paper's results. [verified: README at HEAD]
- **Side-effect set:** the committed list `data/decagon_et.pkl` (identical to `data/dd_edge_type_list_500.pkl`) holds **1,097** side-effect IDs out of 1,317. That is not Decagon's 964 (or 963). The paper's dTIP_D decoder has 964 outputs (Supp. §7), and the repo's `data.ipynb` says Decagon's subset has 963 types. [verified: opcode count of the pickle via `pickletools`, without unpickling; `data.ipynb` cell 2 at `82e63fa`] Which side-effect set produced Table 2, and why the list has 1,097 entries, is **[unverified]**. One possible cause is that `get_edge_list` thresholds `adj.nnz` on the *symmetric* matrix, which counts every pair twice. That was not checked (see Blockers).

**Implication for Malone / TIP:** Malone's and TIP's numbers are random per-triple results with balanced 1:1 random-pair negatives. Each paper puts its own re-runs in the same table as Decagon's numbers copied from the Decagon paper, which used a different split. TIP's "AP@50" column, at least in the released code, is full average precision. It is not comparable with Decagon's or Malone's AP@50.

**Implication:** published random-split results are saturated (AUROC 0.97–0.998), so the leaderboard can no longer tell methods apart.

## Pair-level split

### Dewulf, Stock & De Baets 2021, *Pharmaceuticals* 14(5):429 (PMC8147651)

Source: <https://pmc.ncbi.nlm.nih.gov/articles/PMC8147651/>

- **Data:** the Zitnik et al. (Decagon) dataset, 645 drugs and 963 effects.
- **Model:** "three-step kernel ridge regression". The drug kernel is a Gaussian over single-drug (mono) side-effect vectors. The effect kernel is a Tanimoto similarity over training labels.
- **Four settings** (10-fold CV; regularisation "tuned separately for each prediction task"):
  - **dde^**: random triples. This is Decagon-like tensor completion.
  - **dd^e**: pair-level. Each drug pair goes to a fold *with all its effects*.
  - **d^de**: one unknown drug. Drugs are assigned to folds together with all their combinations.
  - **d^d^e**: two unknown drugs. Only pairs whose two drugs are both test drugs are evaluated.
- **Negatives:** none are sampled. "Any other triplet is assumed to have a zero-label", which gives about 2% positives, so the AUPRC no-skill level is about 0.02, not 0.5.
- **Metrics:** AUROC and AUPRC (no AP@50), averaged two ways after pooling all 10 folds:
  - **E_dd:** per effect, ranking all test pairs, then averaged over effects. This is the macro average comparable to Decagon's.
  - **E_e:** per pair, ranking effects, then averaged over pairs.

Table 1:

| Setting | AUROC E_dd | AUROC E_e | AUPRC E_dd | AUPRC E_e |
|---|---|---|---|---|
| dde^ (random triple) | 0.957 | 0.888 | 0.557 | 0.257 |
| dd^e (**pair-level**) | 0.919 | 0.865 | 0.286 | 0.179 |
| d^de (**one new drug**) | 0.910 | 0.859 | 0.221 | 0.176 |
| d^d^e (**two new drugs**) | 0.843 | 0.834 | 0.112 | 0.144 |

- **Drop from random to pair-level (E_dd):** AUROC falls by 0.038 and AUPRC roughly halves (0.557 → 0.286). From random to two new drugs, AUROC falls by 0.114 and AUPRC falls about 5× (0.557 → 0.112).
- **Rescaling check:** under a balanced (no-skill 0.5) convention, the paper recalculates its dde^ AUPRC as **0.957**. It sets this against the state of the art at the time, AUROC 0.965 and AUPRC 0.938. The authors note that balanced sampling in prior work is "possibly indicating optimistic results".
- **Trivial baseline:** "predicts all the effects as occurring for each drug–drug pair that has at least one effect in the database" scores **AUROC_dd = 0.85** and AUROC_e = 0.5. Under the full-tensor negative scheme, per-effect AUROC is therefore largely earned by separating labelled pairs from unlabelled ones.
- **Per-drug additive baseline:** not reported. Mono side effects are used only as drug features.

## Cold-drug / new-drug splits

### SimVec: Lukashina et al. 2022, *J Cheminform* 14:49 (PMC9327181)

Source: <https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9327181/fullTextXML> (Europe PMC full text)

- **Split ("weak nodes"):**
  - The **98 lowest-degree drugs** are the "new drugs", about 1/6 of the 645.
  - 20,000 "weak triples" (triples with at least one new drug) are split half into test and half into validation.
  - The remaining triples are used for training.
  - This is not a random drug split: it selects drugs with few edges.
- **Comparison split:** a "uniform" split that keeps each node's degree across train and test. The original-model values in this column are copied from the TriVec paper, so they are Decagon's random split.
- **Negatives:** in training, one endpoint is replaced by a uniformly sampled vertex. Evaluation negatives are not described.
- **Metrics:** ROC AUC and AUC-PR (no AP@50), averaged over all relations, with mean ± sd over 10 runs.

Table 1, selected rows:

| Model | Uniform AUROC | Uniform AUPRC | Weak-node AUROC | Weak-node AUPRC |
|---|---|---|---|---|
| RESCAL (original) | 0.693 | 0.613 | 0.44 ± 0.015 | 0.535 ± 0.08 |
| Decagon (original) | 0.872 | 0.832 | **0.30 ± 0.02** | 0.45 ± 0.008 |
| TriVec (original) | 0.975 | 0.966 | 0.44 ± 0.016 | 0.535 ± 0.07 |
| SimVec_full | 0.975 ± 0.002 | 0.968 ± 0.003 | 0.773 ± 0.006 | 0.745 ± 0.003 |

- **Drops (my arithmetic, not the paper's):**
  - Decagon AUROC: 0.872 → 0.30.
  - TriVec AUROC: 0.975 → 0.44.
  - SimVec AUROC: 0.975 → 0.773. SimVec gets its signal from mono side effects and chemical-structure similarity edges.
- **Below-chance results:** AUROC below 0.5 means the trained-graph models rank weak-drug positives *below* the sampled negatives. Embeddings of drugs with few or no training edges are barely trained.
- **Negative sampling changes the cold result substantially** (Table 3, SimVec_SE_weighted, weak split, AUROC):
  - Uniform 6:1: 0.753.
  - NSCaching 6:1: 0.780.
  - "Stay Positive": 0.849.
- **Per-drug additive baseline:** none.

### EmerGNN: Zhang et al. 2023 (Nature Computational Science; arXiv:2311.09261)

Source: <https://arxiv.org/html/2311.09261v1>

- **Different dataset:** this is SumGNN's TWOSIDES subset, with **604 drugs and 200 side effects**, not Decagon's 645 / 964. The numbers are close in spirit but not directly comparable.
- **S1 (one new drug):** drugs are split at random into disjoint sets (seed 1: 514 train / 30 valid / 60 test). Test pairs are (train ∪ valid drug) × test drug.
- **S2 (two new drugs):** both drugs of a test pair are test drugs (2,466 test edges for seed 1).
- **Negatives:** 1:1 for each positive. The new drug is kept and its partner is replaced by a *training* drug v′ such that the triple is not in the network.
- **Metrics:** PR-AUC (primary), ROC-AUC and accuracy, with each side effect treated as a binary task. Mean ± sd over 5 runs. How results are averaged across side effects is not clearly defined.

TWOSIDES results (%):

| Method | S1 PR-AUC | S1 ROC-AUC | S2 PR-AUC | S2 ROC-AUC |
|---|---|---|---|---|
| Decagon | 79.0 ± 2.0 | 78.5 ± 2.3 | n/a | n/a |
| MLP (drug-feature baseline) | 81.5 ± 1.5 | 81.2 ± 1.9 | n/a | n/a |
| SumGNN | 80.3 ± 1.1 | 81.4 ± 1.0 | n/a | n/a |
| CompGCN | 80.3 ± 3.2 | 79.4 ± 4.0 | n/a | n/a |
| HIN-DDI | 81.9 ± 0.6 | 83.8 ± 0.9 | 64.8 ± 2.3 | 58.5 ± 1.6 |
| KG-DDI | 79.1 ± 0.9 | 77.7 ± 1.0 | 53.9 ± 3.9 | 47.0 ± 5.5 |
| EmerGNN | **90.6 ± 0.7** | **91.5 ± 1.0** | **81.4 ± 7.4** | **79.6 ± 7.9** |

- In S1, a plain MLP on drug features matches or beats Decagon and the other graph neural networks, except EmerGNN and HIN-DDI.
- S2 numbers have a large spread (sd about 7–8 points).
- Lloyd 2024 lists EmerGNN's AUPRC on the Decagon data as 0.976. That figure is a random / existing-drug result.
- **Per-drug additive baseline:** none.

## Per-drug additive baseline (`6.5.1`)

- No paper found reports σ(u_r(a)+u_r(b)) or an equivalent "sum of two single-drug scores per side effect" baseline on Decagon data, under any split.
- The nearest evidence that single-drug signal dominates comes from three results:
  - Decagon's concatenated drug-feature baseline: AUROC 0.793 on the random split.
  - Malone et al. 2018's concatenated mono-side-effect + gene-target indicator baseline: AUROC / AuPR / AP@50 of 0.896 / 0.859 / 0.812 on its own random split with balanced negatives. That is above the Decagon numbers it quotes (0.872 / 0.832 / 0.803). [verified: arXiv:1810.09227, Table 2, p. 5] The baseline concatenates features rather than adding per-drug scores, and it was not tested under a stricter split.
  - TIP's dTIP_P (protein and target information only, NN decoder; no drug-drug message passing): AUROC 0.743. [verified: arXiv:2001.10516, Table 2, p. 4]
  - Dewulf's trivial baseline: AUROC_dd 0.85, from the pair-has-any-label signal alone.
  - EmerGNN's MLP baseline: beats Decagon in S1.
- SimVec's cold-drug gains come from mono side-effect features. SimVec_SE reaches AUROC 0.742 against TriVec_SE's 0.397.

## Implications for the split-protocol decision

- **No published number is a direct comparison point for "Decagon-style GNN, Decagon 645/964 data, pair-level split, Decagon's negatives".** Any pair-level or cold-drug result in this project needs its own baselines run under the same protocol.
- **The negative scheme moves AUPRC more than the split does.** Balanced 1:1 sampling gives a no-skill level of 0.5, and full-tensor scoring gives about 0.02. A pair-level AUPRC is comparable to Dewulf only under full-tensor negatives and per-effect (E_dd) averaging.
- **Comparable external numbers by protocol:**
  - Pair-level, full-tensor zeros, per-effect mean, 10-fold over pairs: Dewulf, AUROC 0.919 / AUPRC 0.286.
  - One new drug: Dewulf, 0.910 / 0.221.
  - Two new drugs: Dewulf, 0.843 / 0.112.
  - Lowest-degree-drug split with sampled negatives: SimVec, Decagon 0.30 and TriVec 0.44 AUROC.
  - Random drug split, 604 / 200 subset, 1:1 negatives with a training-drug partner: EmerGNN, Decagon 78.5 and MLP 81.2 ROC-AUC (S1).
- **Aggregation differs between papers:**
  - Decagon: mean over 964 side effects.
  - Lloyd: median over 963.
  - Dewulf: mean over effects after pooling folds.
  - SimVec: mean over relations, then over 10 runs.
- **AP@50 is reported only on the random split.** No stricter-split AP@50 was found. Even on the random split, "AP@50" is not one metric everywhere: TIP's released code computes full average precision and labels it "ap@50" (see the TIP subsection).
- **Copied baselines:** Malone 2018 and TIP 2019 both quote Decagon's paper numbers rather than re-running Decagon on their own split and negatives. Their "beats Decagon" margins therefore compare different protocols.
- **How "cold drug" is defined changes the result:**
  - The lowest-degree drugs (SimVec) give much lower numbers than a random drug split (EmerGNN, Dewulf).
  - Scoring only the pairs whose two drugs are both new (Dewulf d^d^e, EmerGNN S2) is the hardest setting.

## Sources

- Zitnik, Agrawal, Leskovec (2018). Modeling polypharmacy side effects with graph convolutional networks. *Bioinformatics* 34(13):i457. <https://academic.oup.com/bioinformatics/article/34/13/i457/5045770>
- Dewulf, Stock, De Baets (2021). Cold-start problems in data-driven prediction of drug–drug interaction effects. *Pharmaceuticals* 14(5):429. <https://pmc.ncbi.nlm.nih.gov/articles/PMC8147651/>
- Lukashina et al. (2022). SimVec: predicting polypharmacy side effects for new drugs. *J Cheminform* 14:49. <https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9327181/fullTextXML>
- Zhang et al. (2023). Emerging drug interaction prediction enabled by a flow-based graph neural network with biomedical network. arXiv:2311.09261. <https://arxiv.org/html/2311.09261v1>
- Lloyd, Liu, Gaunt (2024). Fast polypharmacy side effect prediction using tensor factorisation. *Bioinformatics* btae706; arXiv:2404.11374. <https://arxiv.org/html/2404.11374v2>
- Malone, García-Durán, Niepert (2018). Knowledge graph completion to predict polypharmacy side effects. DILS 2018, LNCS, pp. 144–149, DOI 10.1007/978-3-030-06016-9_14. Read via arXiv:1810.09227v1: <https://arxiv.org/pdf/1810.09227>
- Xu, Sang, Lu (2019). Tri-graph Information Propagation for polypharmacy side effect prediction. NeurIPS 2019 Workshop on Graph Representation Learning. arXiv:2001.10516v1: <https://arxiv.org/pdf/2001.10516>; workshop PDF: <https://grlearning.github.io/papers/94.pdf>
- TIP code: <https://github.com/NYXFLOWER/TIP>. Files read: `src/utils.py`, `src/neg_sampling.py`, `src/layers.py`, `data/utils.py`, `data/preprocess_data.py`, `prepare.py`, `tip.py`, `model/ddm-df_rgcn.py`, `README.md` at HEAD `cfe3df4`; `tip-cat.py`, `src/utils.py`, `src/neg_sampling.py`, `data.ipynb`, `data/decagon_et.pkl`, `data/dd_edge_type_list_500.pkl` at `82e63fa`; `src/utils.py` at `cc7f616`.

## Blockers

- **Published Malone chapter not read.** `curl -I https://link.springer.com/content/pdf/10.1007/978-3-030-06016-9_14.pdf` returned `HTTP/2 303` with `location: https://idp.springer.com/authorize?...` (a login redirect). The chapter page <https://link.springer.com/chapter/10.1007/978-3-030-06016-9_14> shows "Access this chapter" / "Log in". Any differences between the DILS version and the arXiv poster version are unchecked.
- **Malone code not found.** The paper links no repository, and `gh search repos kblrn` and `gh search repos "KBlrn polypharmacy"` returned no results. Malone's split proportions and the classifier behind its baseline therefore cannot be checked against code.
- **TIP side-effect count (1,097 vs 964) not reconciled.** Counting edges per side effect would mean downloading the 1,317 files `data/sym_adj/drug-sparse-adj/type_*.npz` (about 12.5 MB in total) and reading them with numpy/scipy. `python3 -I -c "import numpy, scipy"` failed with `ModuleNotFoundError: No module named 'numpy'`, and installing packages is out of scope.
- **Tooling note (not a content gap):** `pdftotext` is not installed (`pdftotext not found`), and a page-range PDF read failed with "pdftoppm is not installed". Both PDFs were read whole instead, so nothing is missing. Semantic Scholar was not used.
