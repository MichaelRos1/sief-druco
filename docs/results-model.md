# Model results: tuning, Stage B and the no-protein ablation

Results of the comparison procedure of `docs/spec-model.md` §8, steps 2–4 (#39). Stage A (step 1) is in #37; the cross-fitted scores and the calibration report follow in #40. Metrics as in §7: per side effect, macro-averaged over the 963 side effects, on the full held-out rows unless named. All runs use the Stage A choices (MLP-head decoder, no loss weights), 2 layers, at most 500 epochs, early stopping with patience 30 on validation macro AUPRC.

**Final model: the per-drug additive model with a free table** (`σ(u_r(a) + u_r(b) + b_r)`, §4.3), the Stage B winner on pair-level validation macro AUPRC (0.326 ± 0.003, ahead of every other model in all 9 paired runs). The best GNN is GAT (0.264).

Commits: the runs record the commits they ran at, `af6cfa2` (Stage B) and `6eefbbe` (ablation), before the branch was rebased onto `main`; on `main` they are `2e03433` and `ecd4789`, with the same model code (the rebase only brought in burden-track files).

## Tuning

Each model once on pair-level fold 0, seed 0, over learning rate {3e-3, 1e-2} × d {32, 64, 128} × dropout {0, 0.3}, over the settings that apply to it (the additive models have no d; the free table has no dropout). The winner by validation macro AUPRC is used unchanged on all three splits. Code at `fdd0cf8`; 66 runs, 6.1 h.

| model | lr | d | dropout | val macro AUPRC | runner-up | settings tried |
|---|---|---|---|---|---|---|
| None | 0.01 | 128 | 0.3 | 0.221 | 0.212 | 12 |
| RGCN | 0.003 | 128 | 0 | 0.250 | 0.246 | 12 |
| GraphSAGE | 0.003 | 128 | 0 | 0.259 | 0.251 | 12 |
| GAT | 0.003 | 128 | 0 | 0.269 | 0.256 | 12 |
| Tensor factorisation | 0.01 | 128 | 0 | 0.314 | 0.313 | 12 |
| Additive, free table | 0.01 | – | – | 0.329 | 0.323 | 2 |
| Additive, mono input | 0.01 | – | 0.3 | 0.188 | 0.188 | 4 |

d = 128, the top of the grid, wins for every model that has a d, and most winners were still improving at the 500-epoch cap (47 of 66 runs had their best epoch at ≥ 490). A larger d or more epochs may help every model; the grid and cap are fixed by §8 and were not widened.

## Stage B

Every model at its tuned settings, 3 splits × 3 folds × 3 seeds, the same folds and seeds for every model so runs pair up. Mean ± SD over the 9 runs; test rows unless named. "284 with targets" and "361 without" are the macro AUPRC on the test pairs of drugs with and without protein targets (§4.3, §7). Code at `af6cfa2` (the 7 pair-level fold 0 seed 0 runs are the tuning winners, at `fdd0cf8` / `c1144a3`; the commits between differ only in burden-track files). 171 runs, 22.6 h.

Per-run 95% cluster-bootstrap intervals over test drugs (1,000 resamples) are in the local tables (`local/results/stage_b/results.md`); they are about ±0.02–0.03 macro AUPRC on pair-level and random and ±0.02–0.05 on cold-drug.

### Pair-level

| model | val macro AUPRC | macro AUPRC | macro AUROC | macro AP@50 | 284 with targets | 361 without |
|---|---|---|---|---|---|---|
| None | 0.202 ± 0.017 | 0.206 ± 0.015 | 0.923 ± 0.009 | 0.196 ± 0.034 | 0.242 ± 0.012 | 0.193 ± 0.014 |
| RGCN | 0.247 ± 0.003 | 0.251 ± 0.003 | 0.942 ± 0.001 | 0.267 ± 0.007 | 0.288 ± 0.007 | 0.236 ± 0.006 |
| GraphSAGE | 0.249 ± 0.017 | 0.253 ± 0.017 | 0.941 ± 0.006 | 0.276 ± 0.028 | 0.296 ± 0.018 | 0.232 ± 0.017 |
| GAT | 0.264 ± 0.003 | 0.267 ± 0.002 | 0.945 ± 0.001 | 0.298 ± 0.005 | 0.310 ± 0.005 | 0.244 ± 0.006 |
| Tensor factorisation | 0.310 ± 0.002 | 0.316 ± 0.003 | 0.946 ± 0.001 | 0.405 ± 0.009 | 0.345 ± 0.008 | 0.304 ± 0.008 |
| Additive, free table | 0.326 ± 0.003 | 0.331 ± 0.002 | 0.950 ± 0.000 | 0.394 ± 0.006 | 0.365 ± 0.004 | 0.315 ± 0.006 |
| Additive, mono input | 0.183 ± 0.005 | 0.185 ± 0.005 | 0.891 ± 0.003 | 0.187 ± 0.005 | 0.221 ± 0.006 | 0.174 ± 0.006 |

Reference: Additive, free table. Model − reference over the 9 runs with the same fold and seed; (runs where the model is ahead / 9).

| model | Δ val macro AUPRC | Δ test macro AUPRC |
|---|---|---|
| None | -0.124 ± 0.017 (0/9) | -0.126 ± 0.016 (0/9) |
| RGCN | -0.078 ± 0.002 (0/9) | -0.081 ± 0.004 (0/9) |
| GraphSAGE | -0.076 ± 0.016 (0/9) | -0.078 ± 0.017 (0/9) |
| GAT | -0.062 ± 0.001 (0/9) | -0.064 ± 0.002 (0/9) |
| Tensor factorisation | -0.016 ± 0.002 (0/9) | -0.016 ± 0.003 (0/9) |
| Additive, mono input | -0.143 ± 0.004 (0/9) | -0.147 ± 0.006 (0/9) |

### Random

| model | val macro AUPRC | macro AUPRC | macro AUROC | macro AP@50 | 284 with targets | 361 without |
|---|---|---|---|---|---|---|
| None | 0.200 ± 0.032 | 0.200 ± 0.032 | 0.911 ± 0.027 | 0.160 ± 0.041 | 0.239 ± 0.031 | 0.189 ± 0.035 |
| RGCN | 0.265 ± 0.008 | 0.265 ± 0.008 | 0.943 ± 0.002 | 0.244 ± 0.010 | 0.310 ± 0.011 | 0.250 ± 0.006 |
| GraphSAGE | 0.248 ± 0.027 | 0.248 ± 0.027 | 0.936 ± 0.010 | 0.221 ± 0.036 | 0.298 ± 0.029 | 0.227 ± 0.026 |
| GAT | 0.281 ± 0.011 | 0.281 ± 0.011 | 0.947 ± 0.002 | 0.266 ± 0.016 | 0.334 ± 0.014 | 0.258 ± 0.009 |
| Tensor factorisation | 0.373 ± 0.001 | 0.372 ± 0.001 | 0.962 ± 0.000 | 0.401 ± 0.004 | 0.398 ± 0.003 | 0.371 ± 0.002 |
| Additive, free table | 0.334 ± 0.000 | 0.334 ± 0.001 | 0.944 ± 0.000 | 0.334 ± 0.002 | 0.371 ± 0.002 | 0.322 ± 0.002 |
| Additive, mono input | 0.193 ± 0.003 | 0.193 ± 0.003 | 0.891 ± 0.002 | 0.169 ± 0.005 | 0.235 ± 0.003 | 0.185 ± 0.003 |

Reference: Additive, free table. Model − reference over the 9 runs with the same fold and seed; (runs where the model is ahead / 9).

| model | Δ val macro AUPRC | Δ test macro AUPRC |
|---|---|---|
| None | -0.133 ± 0.032 (0/9) | -0.133 ± 0.032 (0/9) |
| RGCN | -0.068 ± 0.008 (0/9) | -0.068 ± 0.008 (0/9) |
| GraphSAGE | -0.086 ± 0.028 (0/9) | -0.086 ± 0.027 (0/9) |
| GAT | -0.052 ± 0.011 (0/9) | -0.052 ± 0.010 (0/9) |
| Tensor factorisation | +0.040 ± 0.001 (9/9) | +0.039 ± 0.002 (9/9) |
| Additive, mono input | -0.141 ± 0.003 (0/9) | -0.140 ± 0.003 (0/9) |

### Cold-drug

| model | val macro AUPRC | macro AUPRC | macro AUROC | macro AP@50 | 284 with targets | 361 without |
|---|---|---|---|---|---|---|
| None | 0.142 ± 0.012 | 0.129 ± 0.014 | 0.861 ± 0.015 | 0.121 ± 0.040 | 0.151 ± 0.016 | 0.116 ± 0.015 |
| RGCN | 0.106 ± 0.011 | 0.105 ± 0.013 | 0.838 ± 0.007 | 0.090 ± 0.050 | 0.125 ± 0.020 | 0.097 ± 0.005 |
| GraphSAGE | 0.135 ± 0.011 | 0.120 ± 0.014 | 0.842 ± 0.009 | 0.138 ± 0.043 | 0.142 ± 0.030 | 0.110 ± 0.003 |
| GAT | 0.125 ± 0.008 | 0.117 ± 0.010 | 0.832 ± 0.021 | 0.123 ± 0.039 | 0.133 ± 0.024 | 0.117 ± 0.005 |
| Additive, mono input | 0.124 ± 0.008 | 0.119 ± 0.012 | 0.842 ± 0.011 | 0.131 ± 0.029 | 0.143 ± 0.014 | 0.109 ± 0.011 |

Reference: None. Model − reference over the 9 runs with the same fold and seed; (runs where the model is ahead / 9).

| model | Δ val macro AUPRC | Δ test macro AUPRC |
|---|---|---|
| RGCN | -0.036 ± 0.016 (0/9) | -0.024 ± 0.018 (1/9) |
| GraphSAGE | -0.008 ± 0.016 (1/9) | -0.009 ± 0.011 (3/9) |
| GAT | -0.018 ± 0.016 (1/9) | -0.012 ± 0.013 (1/9) |
| Additive, mono input | -0.018 ± 0.013 (1/9) | -0.009 ± 0.015 (2/9) |

### Random split, 1:1 negatives

| model | AUROC | AUPRC | AP@50 |
|---|---|---|---|
| Decagon, published (964 side effects) | 0.872 | 0.832 | 0.803 |
| None | 0.911 ± 0.026 | 0.897 ± 0.024 | 0.873 ± 0.028 |
| RGCN | 0.943 ± 0.002 | 0.930 ± 0.003 | 0.917 ± 0.004 |
| GraphSAGE | 0.936 ± 0.010 | 0.922 ± 0.011 | 0.906 ± 0.015 |
| GAT | 0.947 ± 0.002 | 0.934 ± 0.003 | 0.924 ± 0.004 |
| Tensor factorisation | 0.962 ± 0.000 | 0.952 ± 0.000 | 0.950 ± 0.000 |
| Additive, free table | 0.944 ± 0.000 | 0.936 ± 0.000 | 0.945 ± 0.001 |
| Additive, mono input | 0.890 ± 0.002 | 0.883 ± 0.003 | 0.875 ± 0.003 |

The published Decagon numbers are over 964 side effects; ours are over 963 (§2: side effects with ≥ 500 distinct pairs in `combo`), so the comparison is close but not exact. On 1:1 negatives every model except the mono-input additive one is above the published Decagon numbers, including the no-encoder model.

### What Stage B shows

- **Free per-drug parameters win wherever a test drug was seen in training.** On pair-level, the additive free table is ahead of tensor factorisation by 0.016 and of GAT by 0.062 val macro AUPRC, in all 9 paired runs each. On random, tensor factorisation is first (0.373) and the additive table second (0.334). A model of each drug's own tendency to cause each side effect explains most of what the GNNs learn on these splits.
- **Among the GNNs, GAT is best on pair-level and random** (ahead of RGCN in 9 / 8 of the 9 paired runs on pair-level / random). GraphSAGE and RGCN are close on pair-level (0.249 / 0.247 val); GraphSAGE is the least stable (SD 0.017–0.027: some runs stop early at a lower level).
- **On cold-drug every model is weak (0.10–0.14) and no GNN beats the no-encoder model** (0.142 val, ahead of RGCN in all 9 runs). The free-table models are not applicable there (they learn nothing for an unseen drug). The GNNs stop early (best epoch 50–160), so they overfit the seen drugs; the tuned settings come from pair-level and are reused unchanged (§8 limitation: cold-drug may prefer more regularisation).
- **The 284 drugs with targets score higher than the 361 without in every model**, also in the models that use no proteins (additive table: 0.365 vs 0.315 on pair-level). The gap therefore reflects the drugs (targets are known for better-studied drugs with more labels), not what the proteins add; the ablation below isolates that.
- **The additive free table stops at epoch 1–2 on the random split.** It starts at u = 0, and after one full-batch step each drug's row is about its own excess frequency per side effect, which already ranks well (val macro AUPRC 0.3332); further training lowers it by 0.0002 over 30 epochs, so patience stops it. Its raw scores then sit close to the base rate, hence its low micro AUPRC (0.078). Macro AUPRC, the selection metric, ranks within each side effect and is unaffected; on pair-level it trains normally to ~500 epochs. Relevant for #40: the S regime uses the random-split protocol, so its scores rely on Platt scaling for calibration.
- **The no-encoder model and GraphSAGE vary most across runs** (val macro AUPRC ranges up to 0.08 on random); the additive table and tensor factorisation vary by ≤ 0.006 across folds and seeds.

## No-protein ablation

The best GNN of Stage B (GAT) at the same settings, on drug–drug edges only: no protein nodes, no drug–protein or protein–protein edges (§4.3), 3 splits × 3 folds × 3 seeds, paired with the Stage B GAT runs. Code at `6eefbbe`; 27 runs, 3.2 h.

### Pair-level

| model | val macro AUPRC | macro AUPRC | macro AUROC | macro AP@50 | 284 with targets | 361 without |
|---|---|---|---|---|---|---|
| GAT (no proteins) | 0.201 ± 0.017 | 0.206 ± 0.018 | 0.925 ± 0.008 | 0.190 ± 0.030 | 0.240 ± 0.016 | 0.195 ± 0.021 |
| GAT | 0.264 ± 0.003 | 0.267 ± 0.002 | 0.945 ± 0.001 | 0.298 ± 0.005 | 0.310 ± 0.005 | 0.244 ± 0.006 |

Reference: GAT. Model − reference over the 9 runs with the same fold and seed; (runs where the model is ahead / 9).

| model | Δ test macro AUPRC | Δ 284 with targets | Δ 361 without |
|---|---|---|---|
| GAT (no proteins) | -0.061 ± 0.018 (0/9) | -0.070 ± 0.017 (0/9) | -0.049 ± 0.020 (0/9) |

### Random

| model | val macro AUPRC | macro AUPRC | macro AUROC | macro AP@50 | 284 with targets | 361 without |
|---|---|---|---|---|---|---|
| GAT (no proteins) | 0.231 ± 0.009 | 0.232 ± 0.010 | 0.933 ± 0.004 | 0.196 ± 0.012 | 0.265 ± 0.009 | 0.227 ± 0.012 |
| GAT | 0.281 ± 0.011 | 0.281 ± 0.011 | 0.947 ± 0.002 | 0.266 ± 0.016 | 0.334 ± 0.014 | 0.258 ± 0.009 |

Reference: GAT. Model − reference over the 9 runs with the same fold and seed; (runs where the model is ahead / 9).

| model | Δ test macro AUPRC | Δ 284 with targets | Δ 361 without |
|---|---|---|---|
| GAT (no proteins) | -0.050 ± 0.016 (0/9) | -0.069 ± 0.017 (0/9) | -0.031 ± 0.017 (0/9) |

### Cold-drug

| model | val macro AUPRC | macro AUPRC | macro AUROC | macro AP@50 | 284 with targets | 361 without |
|---|---|---|---|---|---|---|
| GAT (no proteins) | 0.146 ± 0.013 | 0.133 ± 0.003 | 0.868 ± 0.007 | 0.137 ± 0.024 | 0.155 ± 0.013 | 0.120 ± 0.008 |
| GAT | 0.125 ± 0.008 | 0.117 ± 0.010 | 0.832 ± 0.021 | 0.123 ± 0.039 | 0.133 ± 0.024 | 0.117 ± 0.005 |

Reference: GAT. Model − reference over the 9 runs with the same fold and seed; (runs where the model is ahead / 9).

| model | Δ test macro AUPRC | Δ 284 with targets | Δ 361 without |
|---|---|---|---|
| GAT (no proteins) | +0.016 ± 0.009 (8/9) | +0.023 ± 0.020 (8/9) | +0.003 ± 0.006 (7/9) |

### What the ablation shows

- **Proteins help GAT on pair-level and random, most on the 284 drugs with targets.** Without them GAT loses 0.061 / 0.050 test macro AUPRC on pair-level / random (behind in 9 of 9 runs on both); on the 284 drugs with targets it loses 0.070 / 0.069, on the 361 without 0.049 / 0.031. Without proteins GAT falls to the level of the no-encoder model (pair-level val 0.201 vs 0.202).
- **On cold-drug, proteins hurt:** without them GAT gains 0.016 (ahead in 8 of 9 runs), mostly on the drugs with targets (+0.023). For a drug never seen in training, its targets do not carry over what the model learnt for the seen drugs; on this split GAT stops early (best epoch ~160) and the protein graph adds parameters to overfit.
- **Proteins do not reach the final model.** The additive free table uses no proteins and is ahead of GAT with proteins on the 284 drugs with targets too (pair-level 0.365 vs 0.310). Better target data could narrow the GNNs' gap, but only on the drugs that already have targets and only for a model that is not the final one.

