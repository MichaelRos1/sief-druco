# Decagon-style results under pair-level and cold-drug splits

DRAFT - in progress.

## Dewulf, Stock & De Baets 2021 (Pharmaceuticals 14(5):429, PMC8147651)

- Data: Zitnik et al. (Decagon) dataset, 645 drugs, 963 effects.
- Settings: dde^ (random triple), dd^e (pair-level), d^de (one unknown drug), d^d^e (both drugs unknown). 10-fold CV, hyperparameters tuned per task.
- Negatives: none sampled; every unlabelled triple = 0 (about 2% positives).
- Table 1 (pooled; AUROC E_dd / AUROC E_e / AUPRC E_dd / AUPRC E_e):
  - dde^ 0.957 / 0.888 / 0.557 / 0.257
  - dd^e 0.919 / 0.865 / 0.286 / 0.179
  - d^de 0.910 / 0.859 / 0.221 / 0.176
  - d^d^e 0.843 / 0.834 / 0.112 / 0.144
- Trivial baseline (predict all effects for any pair with >=1 known effect): AUROC_dd 0.85, AUROC_e 0.5.
