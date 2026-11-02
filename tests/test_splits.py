from types import SimpleNamespace

import numpy as np
import pytest
import scipy.sparse as sp

from sief_druco import decagon, splits
from sief_druco.splits import TEST, TEST_TWO_NEW, TRAIN, UNUSED, VAL

needs_raw = pytest.mark.skipif(
    not decagon.RAW_DIR.is_dir(), reason=f"raw Decagon files not in {decagon.RAW_DIR}"
)

N = decagon.N_PAIRS


def rng(seed=0):
    return np.random.default_rng(seed)


@pytest.fixture(scope="module")
def toy_labels():
    # Three side effects: one below the 50-positive floor (5% of 600 = 30), two above it.
    g = rng(99)
    cols = [g.choice(N, size=n, replace=False) for n in (600, 2_000, 30_000)]
    rows = np.concatenate(cols)
    col_ids = np.repeat(np.arange(3), [len(c) for c in cols])
    return sp.csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, col_ids)), shape=(N, 3))


@pytest.fixture(scope="module")
def toy_data(toy_labels):
    # Drugs 0..4 have no mono profile and no targets; every other drug has a profile.
    mono = sp.csr_matrix(np.r_[np.zeros(5), np.ones(decagon.N_DRUGS - 5)][:, None])
    return SimpleNamespace(labels=toy_labels, mono=mono, targets=np.empty((0, 2), dtype=int))


# --- Random per-triple ---


def test_random_holds_out_5_percent_of_positives_and_zeros(toy_labels):
    roles = splits.random_roles(toy_labels, rng())
    assert roles.shape == (N, 3)
    dense = toy_labels.toarray().astype(bool)
    for r, n_pos in enumerate((600, 2_000, 30_000)):
        k_pos = max(50, n_pos * 5 // 100)
        k_zero = (N - n_pos) * 5 // 100
        pos, zero = roles[dense[:, r], r], roles[~dense[:, r], r]
        assert (pos == VAL).sum() == (pos == TEST).sum() == k_pos
        assert (zero == VAL).sum() == (zero == TEST).sum() == k_zero
        # Everything not held out is trained on; nothing is unused.
        assert (roles[:, r] == TRAIN).sum() == N - 2 * k_pos - 2 * k_zero


def test_random_floor_of_50_positives(toy_labels):
    roles = splits.random_roles(toy_labels, rng())
    pos = roles[toy_labels.toarray()[:, 0].astype(bool), 0]
    assert (pos == TEST).sum() == 50  # not 5% of 600 = 30


def test_random_rejects_side_effects_too_small_for_the_floor():
    labels = sp.csr_matrix((np.ones(80, dtype=np.uint8), (np.arange(80), np.zeros(80))), shape=(N, 1))
    with pytest.raises(ValueError):
        splits.random_roles(labels, rng())


# --- Pair-level ---


def test_pair_level_80_10_10():
    roles = splits.pair_roles(rng())
    assert roles.shape == (N,)
    assert set(np.unique(roles)) == {TRAIN, VAL, TEST}
    assert (roles == VAL).sum() == (roles == TEST).sum() == N // 10
    assert (roles == TRAIN).sum() == N - 2 * (N // 10)


# --- Cold-drug ---


def test_cold_drug_drug_fractions_and_always_train():
    always = np.array([3, 100, 200, 300, 644])
    drug_roles = splits.cold_drug_drug_roles(always, rng())
    assert (drug_roles == VAL).sum() == (drug_roles == TEST).sum() == 64
    assert (drug_roles == TRAIN).sum() == 645 - 128
    assert np.all(drug_roles[always] == TRAIN)


def test_cold_drug_no_held_out_drug_in_a_training_pair():
    drug_roles = splits.cold_drug_drug_roles(np.array([0]), rng())
    roles = splits.cold_drug_pair_roles(drug_roles)
    a, b = decagon.all_pairs().T
    ra, rb = drug_roles[a], drug_roles[b]
    assert np.all((ra[roles == TRAIN] == TRAIN) & (rb[roles == TRAIN] == TRAIN))


def test_cold_drug_pair_roles_follow_drug_roles():
    drug_roles = splits.cold_drug_drug_roles(np.array([0]), rng())
    roles = splits.cold_drug_pair_roles(drug_roles)
    a, b = decagon.all_pairs().T
    pair = np.sort(np.stack([drug_roles[a], drug_roles[b]], axis=1), axis=1)
    expected = {
        (TRAIN, TRAIN): TRAIN,
        (TRAIN, VAL): VAL,
        (TRAIN, TEST): TEST,
        (TEST, TEST): TEST_TWO_NEW,
        (VAL, VAL): UNUSED,
        (VAL, TEST): UNUSED,
    }
    for (lo, hi), role in expected.items():
        assert np.all(roles[(pair[:, 0] == lo) & (pair[:, 1] == hi)] == role)
    n_train, n_held = 645 - 128, 64
    assert (roles == TRAIN).sum() == n_train * (n_train - 1) // 2
    assert (roles == VAL).sum() == (roles == TEST).sum() == n_held * n_train
    assert (roles == TEST_TWO_NEW).sum() == n_held * (n_held - 1) // 2


# --- Seeds, saving and loading ---


def test_same_seed_same_draw_other_seed_other_draw():
    assert np.array_equal(splits.pair_roles(rng(0)), splits.pair_roles(rng(0)))
    assert not np.array_equal(splits.pair_roles(rng(0)), splits.pair_roles(rng(1)))


@pytest.mark.parametrize("split", splits.SPLITS)
def test_load_returns_the_saved_draw_and_never_redraws(split, toy_data, tmp_path, monkeypatch):
    first = splits.load(split, 1, toy_data, tmp_path)

    def no_redraw(*args):
        raise AssertionError("a saved fold was drawn again")

    monkeypatch.setattr(splits, "draw", no_redraw)
    again = splits.load(split, 1, toy_data, tmp_path)
    assert np.array_equal(first.roles, again.roles)
    assert first.meta == again.meta
    if split == "cold_drug":
        assert np.array_equal(first.drug_roles, again.drug_roles)


@pytest.mark.parametrize("split", splits.SPLITS)
@pytest.mark.parametrize("fold", range(splits.N_FOLDS))
def test_fold_k_is_drawn_from_its_split_seed(split, fold, toy_data):
    g = rng(splits.SEEDS[split][fold])
    if split == "random":
        expected = splits.random_roles(toy_data.labels, g)
    elif split == "pair":
        expected = splits.pair_roles(g)
    else:
        expected = splits.cold_drug_pair_roles(
            splits.cold_drug_drug_roles(splits.no_input_drugs(toy_data), g)
        )
    f = splits.draw(split, fold, toy_data)
    assert np.array_equal(f.roles, expected)
    assert f.meta["seed"] == splits.SEEDS[split][fold]


def test_no_two_splits_share_a_seed():
    seeds = [s for split in splits.SPLITS for s in splits.SEEDS[split]]
    assert len(set(seeds)) == len(seeds) == len(splits.SPLITS) * splits.N_FOLDS


@pytest.mark.parametrize("split", splits.SPLITS)
def test_meta_matches_an_independent_count(split, toy_data):
    f = splits.draw(split, 0, toy_data)
    dense = toy_data.labels.toarray().astype(bool)
    in_role = lambda code: f.roles == code
    if split == "random":  # roles per triple
        positives = lambda code: (in_role(code) & dense).sum()
        test_pos = (in_role(TEST) & dense).sum(axis=0)
    else:  # roles per pair; for cold-drug, "test" is the one-new-drug pairs only
        positives = lambda code: dense[in_role(code)].sum()
        test_pos = dense[in_role(TEST)].sum(axis=0)
    # Names spelled out here, not taken from the module, so a swapped name fails.
    names = {TRAIN: "train", VAL: "val", TEST: "test", UNUSED: "unused", TEST_TWO_NEW: "test_two_new"}
    expected = {
        name: {"n": int(in_role(code).sum()), "positives": int(positives(code))}
        for code, name in names.items()
        if in_role(code).any()
    }
    assert f.meta["counts"] == expected
    assert f.meta["min_test_positives_per_side_effect"] == test_pos.min()
    if split == "cold_drug":
        assert f.meta["drugs"] == {"train": 645 - 128, "val": 64, "test": 64}


def test_unknown_split_raises(toy_data):
    with pytest.raises(ValueError):
        splits.draw("scaffold", 0, toy_data)


# --- The real dataset ---


@pytest.fixture(scope="module")
def data():
    return decagon.load()


@pytest.fixture(scope="module")
def real_folds(data, tmp_path_factory):
    out = tmp_path_factory.mktemp("splits")
    return {s: splits.load(s, 0, data, out) for s in splits.SPLITS}


@needs_raw
def test_real_random_has_50_val_and_test_positives_per_side_effect(data, real_folds):
    roles = real_folds["random"].roles
    rows, cols = data.labels.nonzero()
    for role in (VAL, TEST):
        per_side_effect = np.bincount(cols[roles[rows, cols] == role], minlength=963)
        assert per_side_effect.min() >= 50
    assert real_folds["random"].meta["min_test_positives_per_side_effect"] >= 50


@needs_raw
def test_real_cold_drug_keeps_no_input_drugs_in_training(data, real_folds):
    no_input = splits.no_input_drugs(data)
    assert len(no_input) == 5
    assert np.all(real_folds["cold_drug"].drug_roles[no_input] == TRAIN)


@needs_raw
def test_real_meta_covers_every_split(real_folds):
    for split, f in real_folds.items():
        assert f.meta["seed"] == splits.SEEDS[split][0]
        assert f.meta["counts"]["test"]["positives"] > 0
        assert f.meta["min_test_positives_per_side_effect"] > 0
