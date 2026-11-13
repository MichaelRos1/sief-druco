"""Reader for the score table handed over by Track 1.

Spec: `docs/spec-burden.md` §2 (and `docs/spec-model.md` §1, which writes it).

`read(folder)` opens a hand-off folder: `scores.npy` is memory-mapped (float32
[207690, 963], ~800 MB), so only the rows a caller indexes are read. The other files
are checked against it: `pairs.csv` must list one row per score row in order
(`pair_id` 0, 1, …), `labels.npz` must have the same shape, and `side_effects.csv`
one row per column. The scores are the cross-fitted, Platt-calibrated `p_r(a, b)`;
the labels serve only "known to Decagon" (§10), never the burden.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

FILES = ("scores.npy", "pairs.csv", "labels.npz", "side_effects.csv", "meta.json")
PAIR_COLUMNS = ["pair_id", "stitch_a", "stitch_b", "labelled", "regime"]
SIDE_EFFECT_COLUMNS = [
    "cui", "name", "labelled_pairs",
    "auprc_U", "predicted_U", "observed_U", "auprc_S", "predicted_S", "observed_S",
]


def _unordered(a, b) -> pd.MultiIndex:
    """Pairs of STITCH IDs as (lower, higher), whatever order they come in."""
    a, b = np.asarray(a, dtype=object), np.asarray(b, dtype=object)
    return pd.MultiIndex.from_arrays([np.where(a < b, a, b), np.where(a < b, b, a)])


@dataclass(frozen=True)
class ScoreTable:
    scores: np.ndarray  # memory-mapped float32 [pairs, side effects]
    pairs: pd.DataFrame  # `PAIR_COLUMNS`, row i describes scores[i]
    labels: sp.csr_array  # 0/1 [pairs, side effects]
    side_effects: pd.DataFrame  # `SIDE_EFFECT_COLUMNS`, row j describes scores[:, j]
    meta: dict
    _row_of: pd.Series  # score row by the unordered pair (lower STITCH ID first)

    def rows(self, a, b) -> np.ndarray:
        """Score rows of the unordered pairs (`a[i]`, `b[i]`) of STITCH IDs."""
        rows = self._row_of.reindex(_unordered(a, b))
        if rows.isna().any():
            missing = rows.index[rows.isna()][:5].tolist()
            raise KeyError(f"pairs not in the score table, e.g. {missing}")
        return rows.to_numpy(dtype=np.int64)


def read(folder) -> ScoreTable:
    """Open a hand-off folder and check that its files line up."""
    folder = Path(folder)
    missing = [f for f in FILES if not (folder / f).exists()]
    if missing:
        raise FileNotFoundError(f"{folder} lacks {missing}")

    scores = np.load(folder / "scores.npy", mmap_mode="r")
    pairs = pd.read_csv(folder / "pairs.csv", dtype={"stitch_a": str, "stitch_b": str, "regime": str})
    labels = sp.csr_array(sp.load_npz(folder / "labels.npz"))
    side_effects = pd.read_csv(folder / "side_effects.csv", dtype={"cui": str, "name": str})
    meta = json.loads((folder / "meta.json").read_text())

    if scores.ndim != 2 or scores.dtype != np.float32:
        raise ValueError(f"scores.npy is {scores.dtype} {scores.shape}, not a float32 matrix")
    for name, df, cols in [("pairs.csv", pairs, PAIR_COLUMNS), ("side_effects.csv", side_effects, SIDE_EFFECT_COLUMNS)]:
        if list(df.columns[: len(cols)]) != cols:
            raise ValueError(f"{name} has columns {list(df.columns)}, expected {cols}")
    if scores.shape != (len(pairs), len(side_effects)):
        raise ValueError(
            f"scores.npy is {scores.shape}, but pairs.csv has {len(pairs)} rows "
            f"and side_effects.csv {len(side_effects)}"
        )
    if labels.shape != scores.shape:
        raise ValueError(f"labels.npz is {labels.shape}, scores.npy {scores.shape}")
    if not np.array_equal(pairs.pair_id.to_numpy(), np.arange(len(pairs))):
        raise ValueError("pairs.csv is out of order: pair_id must be 0, 1, … in row order")
    if (pairs.stitch_a == pairs.stitch_b).any():
        raise ValueError("pairs.csv pairs a drug with itself")
    if side_effects.cui.duplicated().any():
        raise ValueError("side_effects.csv lists a side effect twice")
    labelled = np.diff(labels.indptr) > 0
    if not np.array_equal(pairs.labelled.to_numpy().astype(bool), labelled):
        raise ValueError("pairs.csv `labelled` disagrees with labels.npz")
    if not np.array_equal(pairs.regime.to_numpy() == "S", labelled):
        raise ValueError("pairs.csv `regime` must be S exactly for labelled pairs")

    key = _unordered(pairs.stitch_a, pairs.stitch_b)
    if key.duplicated().any():
        raise ValueError("pairs.csv lists a pair twice")
    row_of = pd.Series(pairs.pair_id.to_numpy(), index=key)
    return ScoreTable(scores, pairs, labels, side_effects, meta, row_of)
