"""The three evaluation splits of ADR 0001, drawn once per fold and saved in `local/splits/`.

Spec: `docs/spec-model.md` §6; held-out zeros: #7.

Fold k is drawn from seed k and saved; `load` reads the saved draw and never redraws
it. Roles are uint8 codes: the random split gives one per triple `[207690, 963]`, the
pair-level and cold-drug splits one per pair `[207690]` (cold-drug also keeps its
per-drug roles). Shares are rounded down, as in Decagon's code.

Draw all folds once with `python -m sief_druco.splits`.
"""

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from sief_druco import decagon

SPLITS_DIR = decagon.LOCAL_DIR / "splits"
SPLITS = ("random", "pair", "cold_drug")
N_FOLDS = 3

TRAIN, VAL, TEST, UNUSED, TEST_TWO_NEW = 0, 1, 2, 3, 4
ROLE_NAMES = {TRAIN: "train", VAL: "val", TEST: "test", UNUSED: "unused", TEST_TWO_NEW: "test_two_new"}

RANDOM_PERCENT = 5  # per side effect, to validation and again to test
RANDOM_MIN_POSITIVES = 50
PAIR_PERCENT = 10
COLD_DRUG_PERCENT = 10


@dataclass
class Fold:
    split: str
    fold: int
    roles: np.ndarray  # random: [207690, 963] per triple; pair, cold_drug: [207690] per pair
    drug_roles: np.ndarray | None  # cold_drug only: [645] per drug
    meta: dict


def _share(n: int, percent: int) -> int:
    return n * percent // 100


# --- Drawing ---


def random_roles(labels: sp.spmatrix, rng: np.random.Generator) -> np.ndarray:
    """Per side effect: 5% of its positives (≥ 50) and 5% of its zeros each to val and test.

    The held-out zeros come from labelled and unlabelled pairs alike and are never trained on.
    """
    n_pairs, n_side_effects = labels.shape
    cols = sp.csc_matrix(labels)
    # Side-effect-major while drawing, so each side effect writes one contiguous row.
    roles = np.full((n_side_effects, n_pairs), TRAIN, dtype=np.uint8)
    for r in range(n_side_effects):
        is_pos = np.zeros(n_pairs, dtype=bool)
        is_pos[cols.indices[cols.indptr[r] : cols.indptr[r + 1]]] = True
        order = rng.permutation(n_pairs)
        n_pos = int(is_pos.sum())
        for group, k in (
            (order[is_pos[order]], max(RANDOM_MIN_POSITIVES, _share(n_pos, RANDOM_PERCENT))),
            (order[~is_pos[order]], _share(n_pairs - n_pos, RANDOM_PERCENT)),
        ):
            roles[r, group[:k]] = VAL
            roles[r, group[k : 2 * k]] = TEST
    return np.ascontiguousarray(roles.T)


def pair_roles(rng: np.random.Generator, n_pairs: int = decagon.N_PAIRS) -> np.ndarray:
    """All pairs, labelled and unlabelled together: 80% train, 10% val, 10% test."""
    order = rng.permutation(n_pairs)
    k = _share(n_pairs, PAIR_PERCENT)
    roles = np.full(n_pairs, TRAIN, dtype=np.uint8)
    roles[order[:k]] = VAL
    roles[order[k : 2 * k]] = TEST
    return roles


def cold_drug_drug_roles(always_train: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """10% of drugs to test, a separate 10% to val; `always_train` drugs are never drawn."""
    eligible = np.setdiff1d(np.arange(decagon.N_DRUGS), always_train)
    order = rng.permutation(eligible)
    k = _share(decagon.N_DRUGS, COLD_DRUG_PERCENT)
    roles = np.full(decagon.N_DRUGS, TRAIN, dtype=np.uint8)
    roles[order[:k]] = VAL
    roles[order[k : 2 * k]] = TEST
    return roles


def cold_drug_pair_roles(drug_roles: np.ndarray) -> np.ndarray:
    """A pair's role from its two drugs' roles.

    train + train → train; val + train → val; test + train → test (one new drug);
    test + test → test_two_new; val + val and val + test → unused.
    """
    ra, rb = drug_roles[decagon.all_pairs()].T
    lo, hi = np.minimum(ra, rb), np.maximum(ra, rb)
    roles = np.full(len(ra), UNUSED, dtype=np.uint8)
    roles[(lo == TRAIN) & (hi == TRAIN)] = TRAIN
    roles[(lo == TRAIN) & (hi == VAL)] = VAL
    roles[(lo == TRAIN) & (hi == TEST)] = TEST
    roles[(lo == TEST) & (hi == TEST)] = TEST_TWO_NEW
    return roles


def no_input_drugs(data) -> np.ndarray:
    """Drugs with neither a mono profile nor targets (5 in Decagon)."""
    has_mono = data.mono.getnnz(axis=1) > 0
    has_targets = np.isin(np.arange(decagon.N_DRUGS), data.targets[:, 0])
    return np.flatnonzero(~has_mono & ~has_targets)


def draw(split: str, fold: int, data) -> Fold:
    """Draw fold `fold` of `split` from seed `fold`. Use `load` to get the saved draw."""
    rng = np.random.default_rng(fold)
    drug_roles = None
    if split == "random":
        roles = random_roles(data.labels, rng)
    elif split == "pair":
        roles = pair_roles(rng, data.labels.shape[0])
    elif split == "cold_drug":
        drug_roles = cold_drug_drug_roles(no_input_drugs(data), rng)
        roles = cold_drug_pair_roles(drug_roles)
    else:
        raise ValueError(f"unknown split {split!r}; expected one of {SPLITS}")
    return Fold(split, fold, roles, drug_roles, _meta(split, fold, roles, drug_roles, data.labels))


def _meta(split, fold, roles, drug_roles, labels) -> dict:
    labels = sp.csr_matrix(labels)
    if roles.ndim == 2:  # per triple
        n = np.bincount(roles.ravel(), minlength=len(ROLE_NAMES))
        rows, cols = labels.nonzero()
        pos_roles = roles[rows, cols]
        positives = np.bincount(pos_roles, minlength=len(ROLE_NAMES))
        test_pos = np.bincount(cols[pos_roles == TEST], minlength=labels.shape[1])
    else:  # per pair
        n = np.bincount(roles, minlength=len(ROLE_NAMES))
        per_pair = np.asarray(labels.sum(axis=1)).ravel()
        positives = np.bincount(roles, weights=per_pair, minlength=len(ROLE_NAMES))
        test_pos = np.asarray(labels[roles == TEST].sum(axis=0)).ravel()
    counts = {
        name: {"n": int(n[code]), "positives": int(positives[code])}
        for code, name in ROLE_NAMES.items()
        if n[code]
    }
    meta = {
        "split": split,
        "fold": fold,
        "seed": fold,
        "unit": "triple" if roles.ndim == 2 else "pair",
        "counts": counts,
        "min_test_positives_per_side_effect": int(test_pos.min()),
    }
    if drug_roles is not None:
        meta["drugs"] = {
            name: int((drug_roles == code).sum())
            for code, name in ROLE_NAMES.items()
            if (drug_roles == code).any()
        }
    return meta


# --- Saving and loading ---


def _paths(split: str, fold: int, splits_dir: Path) -> tuple[Path, Path]:
    base = splits_dir / split / f"fold{fold}"
    return base.with_suffix(".npz"), base.with_suffix(".json")


def save(f: Fold, splits_dir: Path = SPLITS_DIR) -> None:
    arrays_path, meta_path = _paths(f.split, f.fold, splits_dir)
    arrays_path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {"roles": f.roles}
    if f.drug_roles is not None:
        arrays["drug_roles"] = f.drug_roles
    # Arrays last: `load` treats the fold as saved once they exist.
    meta_path.write_text(json.dumps(f.meta, indent=2) + "\n")
    np.savez_compressed(arrays_path, **arrays)


def load(split: str, fold: int, data=None, splits_dir: Path = SPLITS_DIR) -> Fold:
    """The saved fold; drawn and saved first if it doesn't exist yet. Never redrawn."""
    arrays_path, meta_path = _paths(split, fold, splits_dir)
    if not arrays_path.exists():
        save(draw(split, fold, data if data is not None else decagon.load()), splits_dir)
    with np.load(arrays_path) as z:
        roles = z["roles"]
        drug_roles = z["drug_roles"] if "drug_roles" in z else None
    return Fold(split, fold, roles, drug_roles, json.loads(meta_path.read_text()))


if __name__ == "__main__":
    data = decagon.load()
    for split in SPLITS:
        for k in range(N_FOLDS):
            f = load(split, k, data)
            print(split, k, f.meta["counts"], "min test positives:",
                  f.meta["min_test_positives_per_side_effect"])
