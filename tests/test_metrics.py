import numpy as np
import pytest
import scipy.sparse as sp
import torch
from sklearn.metrics import average_precision_score, roc_auc_score

from sief_druco import metrics


def synthetic(n=300, k=6, seed=0, decimals=None):
    """Scores that carry some signal, labels at a few percent; rounding makes ties."""
    rng = np.random.default_rng(seed)
    labels = (rng.random((n, k)) < 0.08).astype(np.uint8)
    labels[:2] = 1  # every side effect has positives
    scores = rng.normal(size=(n, k)) + 1.5 * labels
    if decimals is not None:
        scores = scores.round(decimals)
    return scores.astype(np.float32), labels


# --- Per side effect and micro vs sklearn ---


@pytest.mark.parametrize("decimals", [None, 1])
def test_per_side_effect_auprc_auroc_match_sklearn(decimals):
    scores, labels = synthetic(decimals=decimals)
    out = metrics.evaluate(scores, labels)
    for r in range(labels.shape[1]):
        assert out["auprc"][r] == pytest.approx(average_precision_score(labels[:, r], scores[:, r]))
        assert out["auroc"][r] == pytest.approx(roc_auc_score(labels[:, r], scores[:, r]))
    assert out["macro_auprc"] == pytest.approx(out["auprc"].mean())
    assert np.array_equal(out["n_pos"], labels.sum(0))
    assert np.allclose(out["base_rate"], labels.mean(0))


@pytest.mark.parametrize("decimals", [None, 1])
def test_micro_matches_sklearn_on_pooled_triples(decimals):
    scores, labels = synthetic(decimals=decimals)
    out = metrics.evaluate(scores, labels)
    assert out["micro_auprc"] == pytest.approx(average_precision_score(labels.ravel(), scores.ravel()))
    assert out["micro_auroc"] == pytest.approx(roc_auc_score(labels.ravel(), scores.ravel()))
    assert out["micro_base_rate"] == pytest.approx(labels.mean())


def test_mask_scores_only_the_selected_triples():
    scores, labels = synthetic(decimals=1)
    mask = np.random.default_rng(1).random(labels.shape) < 0.5
    mask[:2] = True
    out = metrics.evaluate(scores, labels, mask)
    for r in range(labels.shape[1]):
        y, s = labels[mask[:, r], r], scores[mask[:, r], r]
        assert out["auprc"][r] == pytest.approx(average_precision_score(y, s))
        assert out["auroc"][r] == pytest.approx(roc_auc_score(y, s))
        assert out["n_pos"][r] == y.sum()
    assert out["micro_auprc"] == pytest.approx(average_precision_score(labels[mask], scores[mask]))


@pytest.mark.parametrize("decimals", [None, 1])
def test_weighted_metrics_match_sklearn_sample_weight(decimals):
    # The bootstrap weighs each pair by its drugs' multiplicities.
    scores, labels = synthetic(decimals=decimals)
    w = np.random.default_rng(2).integers(0, 4, size=len(labels))
    auprc, auroc = metrics._Weighted(metrics._by_side_effect(scores, labels, None))(w)
    for r in range(labels.shape[1]):
        assert auprc[r] == pytest.approx(average_precision_score(labels[:, r], scores[:, r], sample_weight=w))
        assert auroc[r] == pytest.approx(roc_auc_score(labels[:, r], scores[:, r], sample_weight=w))


# --- Edge cases ---


def test_perfect_scores_give_one():
    _, labels = synthetic()
    out = metrics.evaluate(labels.astype(np.float32), labels)
    for name in ("auprc", "auroc", "ap50"):
        assert np.allclose(out[name], 1.0)
    assert out["micro_auprc"] == pytest.approx(1.0)


def test_constant_scores_give_the_base_rate():
    _, labels = synthetic()
    out = metrics.evaluate(np.zeros(labels.shape, np.float32), labels)
    assert np.allclose(out["auprc"], labels.mean(0))
    assert np.allclose(out["auroc"], 0.5)
    assert out["micro_auprc"] == pytest.approx(labels.mean())


def test_side_effect_without_positives_is_nan_and_left_out_of_macro():
    scores, labels = synthetic()
    labels[:, 0] = 0
    out = metrics.evaluate(scores, labels)
    assert np.isnan(out["auprc"][0]) and np.isnan(out["auroc"][0]) and np.isnan(out["ap50"][0])
    assert out["macro_auprc"] == pytest.approx(out["auprc"][1:].mean())


def test_sparse_labels_are_accepted():
    scores, labels = synthetic()
    dense = metrics.evaluate(scores, labels)
    sparse = metrics.evaluate(scores, sp.csr_matrix(labels))
    assert dense["macro_auprc"] == sparse["macro_auprc"]


@pytest.mark.parametrize("masked", [False, True])
def test_torch_macro_auprc_matches_evaluate_without_ties(masked):
    scores, labels = synthetic(k=40)
    mask = np.random.default_rng(1).random(scores.shape) < 0.5 if masked else None
    mask_t = None if mask is None else torch.from_numpy(mask)
    fast = metrics.macro_auprc_torch(torch.from_numpy(scores), torch.from_numpy(labels), mask_t, chunk=16)
    assert fast == pytest.approx(metrics.evaluate(scores, labels, mask)["macro_auprc"], rel=1e-5)


# --- AP@50 (Decagon's apk, hand-checked) ---


def ranked_column(hits_in_rank_order):
    """One side effect whose scores rank the given labels top to bottom."""
    y = np.array(hits_in_rank_order, dtype=np.uint8)[:, None]
    s = -np.arange(len(y), dtype=np.float32)[:, None]
    return s, y


def test_ap50_hand_checked():
    # Hits at ranks 1 and 3: (1/1 + 2/3) / min(2, 50).
    s, y = ranked_column([1, 0, 1, 0, 0])
    assert metrics.evaluate(s, y)["ap50"][0] == pytest.approx((1 + 2 / 3) / 2)


def test_ap50_ignores_hits_below_rank_50_and_caps_the_denominator():
    # 60 positives, the top 50 ranks all hits: precision 1 at each, / min(60, 50).
    s, y = ranked_column([1] * 60 + [0] * 40)
    assert metrics.evaluate(s, y)["ap50"][0] == pytest.approx(1.0)
    # Two positives, one at rank 2, one at rank 51: only the first counts.
    s, y = ranked_column([0, 1] + [0] * 48 + [1])
    assert metrics.evaluate(s, y)["ap50"][0] == pytest.approx((1 / 2) / 2)


def test_ap50_counts_only_masked_triples():
    s, y = ranked_column([1, 0, 1, 0])
    mask = np.array([[True], [False], [True], [True]])
    assert metrics.evaluate(s, y, mask)["ap50"][0] == pytest.approx(1.0)


# --- Decagon's 1:1 negatives on the random split ---


def held_out(seed=0):
    rng = np.random.default_rng(seed)
    labels = (rng.random((400, 5)) < 0.1).astype(np.uint8)
    mask = rng.random(labels.shape) < 0.3
    return labels, mask


def test_uniform_negatives_are_one_to_one_per_side_effect_from_its_held_out_zeros():
    labels, mask = held_out()
    sel = metrics.uniform_negatives(labels, mask, seed=0)
    pos, zeros = mask & (labels == 1), mask & (labels == 0)
    assert np.array_equal(sel & (labels == 1), pos)  # every held-out positive kept
    assert not np.any(sel & ~mask)  # nothing outside the held-out triples
    assert np.array_equal((sel & zeros).sum(0), pos.sum(0))  # 1:1 within each side effect


def test_uniform_negatives_are_reproducible_from_the_seed():
    labels, mask = held_out()
    a = metrics.uniform_negatives(labels, mask, seed=0)
    assert np.array_equal(a, metrics.uniform_negatives(labels, mask, seed=0))
    assert not np.array_equal(a, metrics.uniform_negatives(labels, mask, seed=1))


def test_uniform_negatives_need_enough_zeros():
    labels = np.ones((10, 1), np.uint8)
    labels[0] = 0
    with pytest.raises(ValueError):
        metrics.uniform_negatives(labels, np.ones_like(labels, bool), seed=0)


# --- Subgroups ---


def test_subgroups_split_pairs_by_how_many_drugs_have_targets():
    pairs = np.array([[0, 1], [0, 2], [2, 3], [1, 3]])
    has_targets = np.array([True, True, False, False])
    g = metrics.subgroups(pairs, has_targets)
    assert g["targets"].tolist() == [True, False, False, False]
    assert g["no_targets"].tolist() == [False, False, True, False]
    assert g["mixed"].tolist() == [False, True, False, True]


def test_new_drug_groups_split_cold_drug_pairs_by_number_of_new_drugs():
    pairs = np.array([[0, 1], [0, 5], [1, 6], [5, 6]])
    g = metrics.new_drug_groups(pairs, new_drugs=np.array([5, 6]))
    assert g["one_new"].tolist() == [False, True, True, False]
    assert g["two_new"].tolist() == [False, False, False, True]


# --- Cluster bootstrap over drugs ---


def test_resample_weights_are_products_of_drug_multiplicities():
    # Drug 4 is not resampled, so a pair (d, 4) reads off drug d's multiplicity.
    drugs = np.array([0, 1, 2, 3])
    pairs = np.array([[0, 4], [1, 4], [2, 4], [3, 4], [0, 1], [0, 2], [1, 3], [2, 3]])
    rng = np.random.default_rng(0)
    for _ in range(20):
        w = metrics.resample_weights(pairs, drugs, rng)
        m = np.append(w[:4], 1.0)
        assert m[:4].sum() == len(drugs)  # drugs drawn with replacement, as many as there are
        assert np.array_equal(w[4:], m[pairs[4:, 0]] * m[pairs[4:, 1]])


def test_bootstrap_resamples_drugs_not_triples():
    # Cold-drug shape: one test drug paired with many training drugs. Resampling that one
    # drug always gives the same rows back, so the interval has zero width; a bootstrap
    # over triples or pairs would not.
    scores, labels = synthetic(n=200)
    pairs = np.stack([np.zeros(200, int), np.arange(1, 201)], axis=1)
    ci = metrics.cluster_bootstrap(scores, labels, pairs, drugs=np.array([0]), n_resamples=50)
    point = metrics.evaluate(scores, labels)["macro_auprc"]
    assert ci["macro_auprc"] == pytest.approx((point, point))


def test_bootstrap_interval_brackets_the_point_estimate_and_is_reproducible():
    scores, labels = synthetic(n=400)
    rng = np.random.default_rng(3)
    pairs = np.sort(np.stack([rng.choice(40, 400), rng.choice(40, 400)], axis=1), axis=1)
    pairs[pairs[:, 0] == pairs[:, 1], 1] += 1
    ci = metrics.cluster_bootstrap(scores, labels, pairs, n_resamples=200, seed=0)
    point = metrics.evaluate(scores, labels)
    for name in ("macro_auprc", "macro_auroc", "micro_auprc", "micro_auroc"):
        lo, hi = ci[name]
        assert lo < point[name] < hi
    # Same seed, same drugs drawn; MPS sums in float32 in no fixed order, so only ~1e-7 apart.
    again = metrics.cluster_bootstrap(scores, labels, pairs, n_resamples=200, seed=0)
    on_cpu = metrics.cluster_bootstrap(scores, labels, pairs, n_resamples=200, seed=0, device="cpu")
    for name in ci:
        assert ci[name] == pytest.approx(again[name], abs=1e-6)
        assert ci[name] == pytest.approx(on_cpu[name], abs=1e-6)
