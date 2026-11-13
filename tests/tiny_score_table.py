"""A tiny score table in the §2 format (`docs/spec-burden.md`), shared by the burden tests.

4 drugs, 6 pairs, 3 side effects. `tests/` is on pytest's `pythonpath`, so test
files import it as `tiny_score_table`.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

A, B, C, D = "CID000000001", "CID000000002", "CID000000003", "CID000000004"
DRUG_PAIRS = [(A, B), (A, C), (A, D), (B, C), (B, D), (C, D)]
SCORES = np.array(
    [
        [0.1, 0.2, 0.3],  # AB
        [0.5, 0.0, 0.5],  # AC
        [0.9, 0.9, 0.9],  # AD: nobody takes it
        [0.2, 0.2, 0.2],  # BC
        [0.4, 0.4, 0.4],  # BD: nobody takes it
        [0.0, 0.1, 0.0],  # CD
    ],
    dtype=np.float32,
)
LABELS = np.zeros(SCORES.shape, np.uint8)
LABELS[0, 0] = 1  # AB is labelled with side effect 0, at a predicted 0.1


def write_table(folder: Path, scores=SCORES, labels=LABELS, drug_pairs=DRUG_PAIRS) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    np.save(folder / "scores.npy", scores)
    labelled = labels.sum(axis=1) > 0
    pd.DataFrame({
        "pair_id": np.arange(len(drug_pairs)),
        "stitch_a": [a for a, _ in drug_pairs],
        "stitch_b": [b for _, b in drug_pairs],
        "labelled": labelled.astype(np.uint8),
        "regime": np.where(labelled, "S", "U"),
    }).to_csv(folder / "pairs.csv", index=False)
    sp.save_npz(folder / "labels.npz", sp.csr_matrix(labels))
    k = labels.shape[1]
    pd.DataFrame({
        "cui": [f"C{j:07d}" for j in range(k)], "name": [f"effect {j}" for j in range(k)],
        "labelled_pairs": labels.sum(axis=0),
        "auprc_U": 0.5, "predicted_U": 1.0, "observed_U": 1, "auprc_S": 0.6, "predicted_S": 1.0, "observed_S": 1,
    }).to_csv(folder / "side_effects.csv", index=False)
    (folder / "meta.json").write_text(json.dumps({"name": "tiny"}))
    return folder
