---
status: accepted
---

# Three evaluation splits; select models on the pair-level split

Every model is evaluated on three splits of the Decagon data:
- **random per-triple**, Decagon's recipe;
- **pair-level**, where a drug pair leaves training with all its side effects;
- **cold-drug**, where whole drugs leave training.

Random per-triple and pair-level are co-headline results, and **models are selected on pair-level**. Each split serves a different part of the burden pipeline. On NHANES 2017–18, 95.6% of the covered co-prescription weight sits on drug pairs Decagon already labels (#16). Scoring those pairs means recovering missing side effects of a known pair, which is what the random split tests. The remaining pairs are unlabelled: they are the candidates for undocumented burden, and only the pair-level split tests them. Selection uses pair-level because the random split rewards knowing that a pair is labelled. In Dewulf 2021, a rule that gives every labelled pair every side effect already scores AUROC 0.85. Cold-drug is a secondary result: a memorisation check and a forecast for drugs outside the 645. Full protocol: #6.

## Considered options

- **Random per-triple only** (as Decagon, Malone and TIP): comparable to published numbers, but saturated (AUROC up to ~0.98). It never tests a pair the model hasn't seen.
- **Pair-level only:** tests novelty, but ignores the task behind most of the burden weight.

## Consequences

- **Deliberate deviation from Decagon:** on the pair-level and cold-drug splits, training pairs are split into 70% message and 30% supervision edges, redrawn every epoch. Decagon used every training edge for both. That works on its own split, where a test pair is still in the graph through its other side effects. On ours, a test pair is never in the graph, so sharing edges would teach the model a cue that is always missing at test time. The random split keeps Decagon's shared edges, so its numbers stay comparable.
- **Zeros are held out too** (#7). On pair-level and cold-drug, unlabelled pairs are split with the labelled ones (pair-level 80/10/10 on the same folds; cold-drug by drug membership), and only this epoch's supervision pairs (30% of all training pairs, labelled and unlabelled) enter the loss. On the random split, 5% / 5% of each side effect's zero triples go to validation / test and are never trained on. Without this, test zeros would be training zeros. The random split's Decagon-style 1:1 number is unchanged; the held-out zeros add a full-row number for the seen-pair regime.
- **Only the random split is comparable to published numbers,** and even there with 963 side effects, not 964 (#4).
- **Cost:** 3 splits × each compared model × at least 3 folds × 3 seeds. All models share the same folds and seeds, so comparisons are paired, run by run.
