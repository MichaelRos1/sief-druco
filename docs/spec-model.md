# Spec: drug-pair side-effect model (Track 1, RQ1)

Build-ready spec for the model that outputs a calibrated `p_r(a, b)` for any pair of Decagon's 645 drugs and 963 side effects, and for the evaluation protocol that answers RQ1. Every decision below comes from a closed ticket of the [wayfinder map #1](https://github.com/MichaelRos1/sief-druco/issues/1); the ticket is named in each section and holds the reasoning. Background: `PROJECT_CONTEXT.md`. Vocabulary: `GLOSSARY.md`.

Choices that depend on experiment results (the final encoder, the decoder, hyperparameters) are fixed as **procedures**, not outcomes.

---

## 1. What the model delivers

- **Output:** for a drug pair `(a, b)`, 963 probabilities `p_r(a, b)`, one per side effect r. Symmetric: `p_r(a, b) = p_r(b, a)`.
- **Meaning:** `p_r(a, b)` estimates **P(pair reported with side effect r in TWOSIDES)**, not P(pair causes r). Unlabelled triples are treated as plain zeros, with no positive-unlabeled correction; reporting bias is a stated limitation (#7).
- **Hand-off to the burden track:** cross-fitted, Platt-calibrated scores for all 207,690 pairs × 963 side effects (`9`). The exact interface (format, IDs, metadata) is settled in #17.
- **Interfaces kept ready for drug sets (#21):** the decoder takes only `(z_a, z_b, r)`; the model exposes `score(pairs) → [n, 963]` calibrated probabilities.

## 2. Data (#4, #5)

**Files** (from `https://snap.stanford.edu/decagon/`, SHA-256 checksums in `local/data/raw/decagon/SHA256SUMS`; facts in `research/decagon-files.md`): `bio-decagon-combo.csv`, `bio-decagon-mono.csv`, `bio-decagon-ppi.csv`, `bio-decagon-targets.csv`. `targets-all.csv` is not used: for the 645 drugs it is identical to `targets.csv`.

| Item | Rule | Result |
|---|---|---|
| Side effects | keep side effects with ≥ 500 distinct pairs in `combo` | **963** (not the paper's 964; state this when comparing) |
| Drugs | every drug in `combo` | 645 |
| Labelled pairs | unordered pairs with ≥ 1 kept side effect | 63,472 of 207,690 possible pairs |
| Triples | (pair, side effect) after the filter | 4,576,287 |
| Proteins | all proteins in `ppi.csv` plus the 8 target proteins missing from it | **19,089** |
| Targets | `targets.csv`, 645 drugs only | 284 drugs, 18,690 edges; 361 drugs have none |
| Mono profiles | `mono.csv` | 639 drugs; 6 have none; 5 have neither profile nor targets |

- **Indexing:** STITCH ID → 0..644, Entrez ID → 0..19,088, side-effect CUI → 0..962. Keep the reverse mappings (with names).
- **Direction:** pairs and PPI edges are stored once, in one order. Split on **unordered** pairs; add reverse edges in code **after** splitting.
- Mono and combo side-effect CUIs are disjoint: never join them.

**Graph** (`HeteroData`): node types `drug` (645), `protein` (19,089). Edge types: protein–protein (PPI), drug → protein (targets), protein → drug (reverse targets), and drug–drug edges, one type per side effect (963) for RGCN or a single collapsed type for the other encoders. RGCN: `num_relations = 966`.

## 3. Input layer (#5)

- **Drugs:** mono profile only, with **no learned per-drug embedding**, for every encoder and every split.
  - Columns: mono side effects present in **≥ 2 training drugs**, recomputed per split (8,152 when all 645 drugs are training drugs). Fixed, not tuned.
  - Vector: the **mean** of the drug's mono side-effect indicators over those columns, plus **`log(1 + n)`** as one extra input (n = number of mono side effects).
  - The 6 drugs without a profile get all zeros (n = 0); they stay in the graph and in every split.
  - → `Linear(n_cols + 1 → d)`.
- **Proteins:** no input features; `Embedding(19089, d)`.

## 4. Models

### 4.1 Encoders (#8 Q1)

All have 2 layers, ReLU, dropout (tuned), output `z_drug [645, d]`.

| Encoder | Definition |
|---|---|
| **None** | `z = Linear(drug input)`; no message passing. Tests whether the graph helps on the same input. |
| **RGCN** (reference) | Relational GCN over all 966 relations, basis decomposition with **30 shared bases**, mean aggregation, root weight. Implemented as the **dense drug–drug rewrite** (drug–drug relations as a dense `[963, 645, 645]` mean adjacency combined through the bases; protein relations by `index_add_`), not `RGCNConv` (11× slower, 11.8 GB). |
| **GraphSAGE** | Collapsed graph: one drug–drug edge type ("labelled together"), plus PPI and target edges; `SAGEConv` via `to_hetero`. |
| **GAT** | As GraphSAGE, with `GATConv` (4 heads × d/4). |

HGT and TIP-style staged propagation are future work (#24).

### 4.2 Decoders (#8 Q3, Q9, Q10)

All take only `(z_a, z_b)` and return 963 logits, and all add a **per-side-effect bias `b_r`** started at the logit of side effect r's training base rate. The bias is a deviation from Decagon: with full rows ~97.7% of labels are 0, and without it the untrained model predicts 0.5 everywhere.

| Decoder | Logit for side effect r |
|---|---|
| **Decagon** (default) | `z_aᵀ D_r R D_r z_b + b_r`, D_r diagonal ([963, d]), R shared [d, d] and made symmetric as `(R₀ + R₀ᵀ) / 2`. Compute as one matmul: `((z_a ⊗ z_b) ⊙ R).flatten(1) @ (D ⊗ D)ᵀ` (8× faster on MPS than broadcasting to `[n, 963, d]`). |
| **DistMult** | `Σ_k z_a[k] · w_r[k] · z_b[k] + b_r` |
| **MLP head** | `MLP([z_a + z_b, z_a ⊙ z_b]) → 963`, + `b_r` |

**Start values** (Decagon decoder): `D_r ~ N(0, 1)`, `R₀ ~ N(0, 1/√d)`. Fixed, not tuned. With 0.1-scale start values the four-factor product has near-zero gradients and training stalls for 30+ epochs.

### 4.3 Baselines and ablations (#8 Q2, Q4, Q14; #5)

| Model | Definition | Splits |
|---|---|---|
| **Tensor factorisation** | free per-drug vectors `Embedding(645, d)` + the chosen decoder, no mono input, no graph | random, pair-level (cold-drug: "not applicable") |
| **Per-drug additive, free table** | `σ(u_r(a) + u_r(b))`, `u` a learned `[645, 963]` table | random, pair-level |
| **Per-drug additive, mono input** | as above, `u(drug) = Linear(drug input) → 963` | all three |
| **No-protein ablation** | the best GNN from Stage B, retrained on drug–drug edges only (no PPI, no targets) | all three |

Every metric is also reported separately for pairs among the **284 drugs with targets** vs the **361 without** (#5). If the no-protein ablation shows that proteins clearly help on the 284, look for a free target source (STITCH, ChEMBL, DGIdb, Guide to Pharmacology); DrugBank is not an option.

## 5. Training (#7, #8)

- **Full rows:** every supervision pair enters the loss with all 963 labels. No negatives are sampled.
- **Loss:** `BCEWithLogitsLoss`, no `pos_weight`, no focal loss. Per-side-effect weights `{none, 1/√n_pos(r)}` (normalised to mean 1), chosen once in Stage A (`8`) and then fixed for every model.
- **Optimiser:** full-batch Adam, one encoder pass and one step per epoch. The decoder runs in chunks of pairs; pass its gradient back exactly through a detached `z` (`z.detach().requires_grad_()`, then `z.backward(z_leaf.grad)`).
- **Early stopping:** on validation macro AUPRC (`7`), patience 30, cap 500 epochs; keep the best epoch's weights.
- **Sampling switch** (fallback only): zeros drawn uniformly per side effect at one rate β, redrawn every epoch, corrected with `p = β·p_s / (β·p_s − p_s + 1)` before calibration. Use it only if a model exceeds **12 GB MPS or 10 s per epoch**; none of the compared models does. Within a split, all compared models use the same variant; each run records it. Decagon's degree-based sampler is not used (it can never draw 29% of the zeros).

**Message vs supervision edges:**
- **Pair-level and cold-drug:** each epoch, training pairs (labelled and unlabelled) are redrawn into 70% message / 30% supervision. Only labelled message pairs create edges. Only supervision pairs' full rows enter the loss. Validation and test pairs never enter the message graph (in either direction). At validation and test time, the message graph is all training labelled pairs.
- **Random split:** as Decagon, every training triple is both a message edge and in the loss.

**Measured cost on the M4** (d = 64, 2 layers, pair-level / random epoch, without validation): none 0.8 / 3.4 s; RGCN dense 1.3 / 4.1 s; SAGE 2.1 / 5.4 s; GAT 2.3 / 5.6 s; ≤ 5.6 GB MPS. About 1 s of the collapsed encoders' time is building the collapsed edge index, which can be precomputed. Validation over full rows adds ~0.3 s.

**MPS notes (#3):** a minimal `torch` + `torch_geometric` install (no extensions) runs everything. Avoid a sorted `EdgeIndex` / `SparseTensor` for GAT (needs `aten::segment_reduce`, no MPS kernel). Avoid `FastRGCNConv` (~51 GB).

## 6. Splits (#6, #7, ADR 0001)

Three splits; every model runs on all three (except where `4.3` says otherwise). **Models are selected on pair-level.**

| Split | Role | Construction |
|---|---|---|
| **Random per-triple** | co-headline; "seen pair, missing side effect" | Decagon's recipe: per side effect, 5% of its labelled pairs to validation and 5% to test (≥ 50 each). Additionally 5% / 5% of each side effect's **zero** triples (from labelled and unlabelled pairs) go to validation / test and are never trained on. |
| **Pair-level** | co-headline and **model selection**; "unseen pair" | All 207,690 unordered pairs (labelled and unlabelled together), each with all 963 labels: 80% train / 10% validation / 10% test. |
| **Cold-drug** | secondary: memorisation check, forecast for new drugs | Drugs at random: 10% test, a separate 10% validation; the 5 drugs with neither a mono profile nor targets always stay in training. Training = pairs of two training drugs. Scored on **one-new-drug pairs** (macro over 963). Two-new-drug pairs are reported only as a pooled micro average. |

- **Repeats:** 3 folds × 3 seeds per model and split, the **same folds and seeds for every model**, so comparisons are paired run by run. Folds are never re-drawn; each fold records its smallest number of test positives per side effect.

## 7. Evaluation (#6, #7, #8)

- **Metrics per side effect:** AUPRC (average precision), AUROC, AP@50; macro-averaged over 963, plus micro averages. AUPRC is read against each side effect's base rate.
- **Evaluation rows:**
  - Pair-level and cold-drug: the **full rows** of the held-out pairs, at real base rates (~20M triples per fold).
  - Random split, two numbers: Decagon's 1:1 uniform negatives per side effect (comparable to the published 0.872 / 0.832 / 0.803), and the full held-out number on the held-out positives + held-out zeros.
- **Selection metric:** **macro AUPRC over the full validation rows of the pair-level split.** Every other metric is reported only.
- **Uncertainty:** mean ± spread over the 9 runs; paired differences between models; per run a 95% interval from a **cluster bootstrap over test drugs** (1,000 resamples).
- **Subgroups:** 284 drugs with targets vs 361 without (`4.3`).

## 8. Comparison procedure (#8)

Run in this order:

1. **Stage A: decoder and loss weights** (~3 h). Decagon vs DistMult vs MLP head, each with encoder None and RGCN, on pair-level, 1 fold × 3 seeds, at default settings (learning rate 1e-2, d = 64, dropout 0). The loss-weight switch `{none, 1/√n_pos}` is tested on None. Pick both by pair-level validation macro AUPRC. **Also record epochs to convergence** and revisit the budget below.
2. **Tuning** (~10 h). Each encoder with the chosen decoder, once, on pair-level fold 0, one seed: learning rate {3e-3, 1e-2} × d {32, 64, 128} × dropout {0, 0.3}. Layers fixed at 2. The winning settings are reused unchanged on the random and cold-drug splits (limitation: cold-drug may prefer more regularisation).
3. **Stage B: main comparison** (~33 h). The 4 encoders + tensor factorisation + both additive baselines, 3 splits × 3 folds × 3 seeds.
4. **No-protein ablation** (~7 h). The best GNN from Stage B on drug–drug edges only.

The **final model** is the Stage B winner on pair-level validation macro AUPRC.

**Budget:** ≈ 50 h of MacBook compute (5–6 nights), assuming ~300 epochs per run (from a pilot stopped at epoch 180, still improving). It scales with that number: 150 epochs → ~25 h, 500 → ~85 h. Cheapest cuts if needed: drop GAT (~11 h), fix d = 64 (~5 h).

## 9. Calibration and burden scores (#9)

1. **Cross-fitting:** the final model, with its tuned settings, is retrained so that no pair is scored by a model that trained on it. One seed per fold model; ~5–7 h.
   - **U regime** (unseen pair) → scores for **unlabelled pairs**: pair-level protocol, all 207,690 pairs in 10 folds. Each run holds out one fold (scored), takes a validation slice from the rest, trains on the remainder.
   - **S regime** (seen pair, missing side effect) → scores for **every triple of labelled pairs**: random-split protocol, triples in 10 folds. Each run holds out one fold of triples (positives and zeros) plus a validation slice; the pair's other side effects stay in the message graph.
2. **Platt scaling per side effect:** `p = σ(a_r · logit + c_r)`, fitted per fold model and regime on that run's early-stopping validation slice. Plain PyTorch: all 963 side effects at once, minimising BCE (seconds). Not temperature scaling, not isotonic, not a library.
3. **Reporting** (final model only, per regime, on the held-out folds): reliability diagrams pooled and for rare / median / common side effects; ECE with 15 equal-count bins; predicted vs observed number of positives per side effect.

Whether the burden uses a labelled side effect's label or its S-regime prediction is decided in #18.

## 10. Deviations from Decagon

| Decagon | This project | Why | Ticket |
|---|---|---|---|
| 964 side effects | 963 | the public file gives 963 under every counting rule | #4 |
| random per-triple split only | + pair-level (selection) and cold-drug | the random split rewards knowing a pair is labelled | #6, ADR 0001 |
| every training edge is message and supervision | 70/30 redraw on pair-level and cold-drug | a held-out pair is never in the graph | #6 |
| degree-based negative sampling, 1:1 | full rows, no sampling | the sampler never draws 29% of the zeros | #7 |
| hinge loss (code) / cross-entropy (paper) | plain BCE | probabilities are needed | #7 |
| decoder without bias | per-side-effect bias `b_r` | full rows are ~97.7% zeros | #8 |
| uncalibrated scores | cross-fitted, Platt per side effect | burden multiplies `p` by prevalence | #9 |

## 11. Not in this spec

- **Model ↔ burden interface:** #17.
- **Burden formula, label vs prediction for labelled pairs:** #18.
- **HGT, TIP:** #24. **Drug sets of 3+ drugs:** #21. **Severity weights:** #20.
- **Drugs outside the 645** (inductive embedding, `PROJECT_CONTEXT.md` `10.5`): out of scope; coverage is reported instead.
