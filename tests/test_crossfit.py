import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import torch

from sief_druco import crossfit, decagon
from sief_druco.splits import TEST, TRAIN, VAL

needs_raw = pytest.mark.skipif(
    not decagon.RAW_DIR.is_dir(), reason=f"raw Decagon files not in {decagon.RAW_DIR}"
)

# --- A tiny synthetic dataset (no real data needed) ---

N_DRUGS, N_SIDE_EFFECTS = 40, 7


def tiny_data(seed=0):
    """Decagon-shaped: 40 drugs, 7 side effects, about a third of the pairs labelled."""
    rng = np.random.default_rng(seed)
    pairs = decagon.all_pairs(N_DRUGS)
    labelled = rng.random(len(pairs)) < 0.3
    dense = (rng.random((len(pairs), N_SIDE_EFFECTS)) < 0.3) & labelled[:, None]
    dense[labelled & ~dense.any(axis=1), 0] = True  # a labelled pair has ≥ 1 side effect
    return SimpleNamespace(
        drug_ids=np.array([f"CID{i:09d}" for i in range(N_DRUGS)]),
        side_effect_ids=np.array([f"C{r:07d}" for r in range(N_SIDE_EFFECTS)]),
        side_effect_names=np.array([f"effect {r}" for r in range(N_SIDE_EFFECTS)]),
        pairs=pairs,
        labels=sp.csr_matrix(dense.astype(np.uint8)),
    )


def folds_of(regime, data, seed=0):
    rng = np.random.default_rng(seed)
    if regime == "U":
        return crossfit.u_folds(len(data.pairs), rng)
    return crossfit.s_folds(data.labels, rng)


# --- Folds: no pair or triple is scored by a model that trained on it ---


@pytest.mark.parametrize("regime", crossfit.REGIMES)
def test_run_k_scores_fold_k_validates_on_k_plus_1_and_trains_on_the_rest(regime):
    folds = folds_of(regime, tiny_data())
    for k in range(crossfit.N_FOLDS):
        roles = crossfit.roles(folds, k)
        assert roles.shape == folds.shape
        np.testing.assert_array_equal(roles == TEST, folds == k)
        np.testing.assert_array_equal(roles == VAL, folds == (k + 1) % crossfit.N_FOLDS)
        np.testing.assert_array_equal(roles == TRAIN, (folds != k) & (folds != (k + 1) % 10))


@pytest.mark.parametrize("regime", crossfit.REGIMES)
def test_every_scored_entry_comes_from_the_one_model_that_held_it_out(regime):
    data = tiny_data()
    folds = folds_of(regime, data)
    shape = data.labels.shape
    # Fold model k "scores" every entry it is asked for with the value k.
    per_fold = (np.full(crossfit.held_out(folds, k, shape[1]).sum(), k, dtype=np.float32)
                for k in range(crossfit.N_FOLDS))
    oof = crossfit.out_of_fold(folds, per_fold, shape[1])
    scored = ~np.isnan(oof)
    rows, cols = np.nonzero(scored)
    for k in range(crossfit.N_FOLDS):
        is_k = oof[rows, cols] == k
        roles = np.broadcast_to(crossfit.roles(folds, k).reshape(len(folds), -1), shape)
        assert (roles[rows[is_k], cols[is_k]] == TEST).all()  # never TRAIN or VAL there
    labelled = data.labels.getnnz(axis=1) > 0
    if regime == "U":
        assert scored.all()  # all pairs, full rows
    else:
        np.testing.assert_array_equal(scored, np.broadcast_to(labelled[:, None], shape))


def test_u_folds_split_all_pairs_evenly():
    folds = crossfit.u_folds(1003, np.random.default_rng(0))
    counts = np.bincount(folds, minlength=crossfit.N_FOLDS)
    assert counts.sum() == 1003 and counts.max() - counts.min() <= 1


def test_s_folds_split_each_side_effects_positives_and_zeros_of_labelled_pairs_evenly():
    data = tiny_data()
    folds = crossfit.s_folds(data.labels, np.random.default_rng(0))
    labelled = data.labels.getnnz(axis=1) > 0
    assert (folds[~labelled] == crossfit.ALWAYS_TRAIN).all()  # unlabelled pairs: trained on in every run
    dense = data.labels.toarray().astype(bool)
    for r in range(N_SIDE_EFFECTS):
        for is_pos in (True, False):
            counts = np.bincount(folds[labelled & (dense[:, r] == is_pos), r], minlength=10)
            assert counts.max() - counts.min() <= 1


@needs_raw
def test_folds_on_the_real_data():
    data = decagon.load()
    u = crossfit.u_folds(decagon.N_PAIRS, np.random.default_rng(0))
    assert np.bincount(u).min() == decagon.N_PAIRS // 10
    s = crossfit.s_folds(data.labels, np.random.default_rng(0))
    labelled = data.labels.getnnz(axis=1) > 0
    assert labelled.sum() == 63_472
    assert (s[labelled] >= 0).all() and (s[~labelled] == crossfit.ALWAYS_TRAIN).all()


# --- Platt scaling and ECE ---


def test_platt_recovers_known_parameters_and_lowers_ece():
    rng = np.random.default_rng(0)
    n, k = 40_000, 5
    a_true = np.array([0.5, 0.8, 1.0, 1.5, 2.0])
    c_true = np.array([-1.0, 1.5, 1.0, 1.0, 0.5])
    cols = np.tile(np.arange(k), n)
    x = rng.normal(-2.5, 1.5, n * k)  # miscalibrated logits
    y = rng.random(n * k) < 1 / (1 + np.exp(-(a_true[cols] * x + c_true[cols])))

    a, c = crossfit.fit_platt(x, y, cols, k)
    np.testing.assert_allclose(a, a_true, atol=0.06)
    np.testing.assert_allclose(c, c_true, atol=0.12)

    raw = crossfit.platt(x, np.ones(k), np.zeros(k), cols)
    calibrated = crossfit.platt(x, a, c, cols)
    ece_raw, _ = crossfit.ece(raw, y)
    ece_cal, bins = crossfit.ece(calibrated, y)
    assert ece_cal < ece_raw / 3
    assert ece_cal < 0.01
    for r in range(k):  # and per side effect
        col = cols == r
        assert crossfit.ece(calibrated[col], y[col])[0] < crossfit.ece(raw[col], y[col])[0]
    assert len(bins) == crossfit.N_BINS
    assert sum(b["n"] for b in bins) == n * k


def test_platt_gives_float32_strictly_inside_0_1_even_for_extreme_logits():
    x = np.array([-1e4, -200.0, -50.0, 0.0, 50.0, 200.0, 1e4])
    p = crossfit.platt(x, np.ones(1), np.zeros(1), np.zeros(len(x), dtype=np.int64))
    assert p.dtype == np.float32
    assert ((p > 0) & (p < 1)).all()
    assert (np.diff(p) >= 0).all()


def test_ece_uses_equal_count_bins():
    p = np.linspace(0.01, 0.99, 150)
    y = np.zeros(150)
    value, bins = crossfit.ece(p, y, n_bins=15)
    assert [b["n"] for b in bins] == [10] * 15
    assert value == pytest.approx(p.mean())  # every bin over-predicts by its mean p


def test_calibration_report_groups_side_effects_by_training_positives():
    rng = np.random.default_rng(0)
    p = rng.random((300, 9)).astype(np.float32)
    y = (rng.random((300, 9)) < p).astype(np.uint8)
    n_pos = np.array([5, 900, 50, 7, 300, 60, 2, 1000, 55])
    report, predicted, observed = crossfit.calibration_report(p, y, n_pos, raw=p)
    assert report["groups"]["rare"]["side_effects"] == [6, 0, 3]
    assert report["groups"]["median"]["side_effects"] == [2, 8, 5]
    assert report["groups"]["common"]["side_effects"] == [4, 1, 7]
    np.testing.assert_allclose(predicted, p.sum(axis=0, dtype=np.float64), rtol=1e-6)
    np.testing.assert_array_equal(observed, y.sum(axis=0))
    assert report["ece"] == pytest.approx(report["ece_raw"])
    assert 0 <= report["ece"] < 0.1


# --- The hand-off files ---


def write_tiny(tmp_path):
    data = tiny_data()
    n_pairs = len(data.pairs)
    rng = np.random.default_rng(1)
    u = rng.random((n_pairs, N_SIDE_EFFECTS)).astype(np.float32) * 0.5 + 0.01  # in (0, 0.51)
    s = u + 0.4  # distinguishable: S scores lie above 0.41
    labelled = data.labels.getnnz(axis=1) > 0
    scores = crossfit.combine(u, s, labelled)
    side_effects = pd.DataFrame({"cui": data.side_effect_ids})
    crossfit.write_handoff(tmp_path, data, scores, side_effects, {"seed": 0})
    return data, u, s, labelled


def test_handoff_rows_follow_the_canonical_pair_order_and_regimes(tmp_path):
    data, u, s, labelled = write_tiny(tmp_path)
    assert {p.name for p in tmp_path.iterdir()} >= {
        "scores.npy", "pairs.csv", "labels.npz", "side_effects.csv", "meta.json"}
    pairs = pd.read_csv(tmp_path / "pairs.csv")
    np.testing.assert_array_equal(pairs["pair_id"], np.arange(len(data.pairs)))
    a, b = decagon.all_pairs(N_DRUGS).T
    np.testing.assert_array_equal(pairs["stitch_a"], data.drug_ids[a])
    np.testing.assert_array_equal(pairs["stitch_b"], data.drug_ids[b])
    np.testing.assert_array_equal(pairs["labelled"], labelled)
    np.testing.assert_array_equal(pairs["regime"], np.where(labelled, "S", "U"))

    labels = sp.load_npz(tmp_path / "labels.npz")
    assert (labels != data.labels).nnz == 0 and labels.shape == data.labels.shape

    scores = np.load(tmp_path / "scores.npy")
    np.testing.assert_array_equal(scores[labelled], s[labelled])  # predictions, never labels
    np.testing.assert_array_equal(scores[~labelled], u[~labelled])
    assert json.loads((tmp_path / "meta.json").read_text())["seed"] == 0


def test_scores_have_one_row_per_unordered_pair_and_lie_in_0_1(tmp_path):
    data, *_ = write_tiny(tmp_path)
    scores = np.load(tmp_path / "scores.npy")
    assert scores.dtype == np.float32
    assert scores.shape == (N_DRUGS * (N_DRUGS - 1) // 2, N_SIDE_EFFECTS)
    assert ((scores > 0) & (scores < 1)).all()
    pairs = pd.read_csv(tmp_path / "pairs.csv")
    a = np.searchsorted(data.drug_ids, pairs["stitch_a"])
    b = np.searchsorted(data.drug_ids, pairs["stitch_b"])
    assert (a < b).all()
    # Each unordered pair once, so (a, b) and (b, a) share a row: symmetric by construction.
    np.testing.assert_array_equal(np.stack([a, b], axis=1), decagon.all_pairs(N_DRUGS))
