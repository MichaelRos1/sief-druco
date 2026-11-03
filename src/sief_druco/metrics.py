"""Evaluation metrics: every number spec §7 asks for, from scores and labels.

Spec: `docs/spec-model.md` §7 (metrics, evaluation rows, bootstrap), §4.3 (subgroups),
§6 (cold-drug scoring).

Independent of any model or split. Inputs are plain arrays over the evaluated pairs:
`scores` `[n_pairs, n_side_effects]` (logits or probabilities; only the ranking within
a side effect matters), `labels` 0/1 of the same shape (dense or sparse), an optional
boolean `mask` of the triples that count (the random split's held-out triples; `None` =
full rows), and `pairs` `[n_pairs, 2]` with each row's drugs, for subgroups and the
bootstrap. How the splits use it:

- Pair-level: `evaluate` on the full rows of the held-out pairs.
- Cold-drug: `new_drug_groups(pairs, test_drugs)` splits the test pairs. One-new-drug
  pairs: `evaluate` as above. Two-new-drug pairs: report only the `micro_*` numbers.
- Random split: `evaluate` with `mask` = held-out triples (full held-out number), and with
  `mask = uniform_negatives(labels, mask, seed)` (Decagon's 1:1 number).
- Subgroups: `subgroups(pairs, has_targets)` gives row masks; evaluate each subset.
- Early stopping: `macro_auprc_torch` gives `evaluate`'s macro AUPRC on the GPU, every epoch.

AUPRC is average precision with ties grouped, as in `sklearn.metrics.average_precision_score`;
AUROC as `roc_auc_score`; AP@50 is Decagon's `apk(k=50)` on the ranking (ties keep row
order). AP@50 has no micro average: the top 50 of millions of pooled triples says little.
A side effect without positives (or, for AUROC, without zeros) gets NaN and is left out
of the macro average.
"""

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import torch

AP_AT = 50


@dataclass
class _Ranked:
    """Triples ranked by score, one row per group (a side effect, or all pooled).

    Rows are padded to equal length; padding has `pair = n_pairs` (weight 0) and comes last.
    """

    pair: np.ndarray  # [groups, m], the evaluated pair at each rank
    y: np.ndarray  # [groups, m] float 0/1, label at each rank
    first: np.ndarray | None  # [groups, m], first rank of each rank's tie group; None = no ties
    last: np.ndarray | None  # [groups, m], last rank of the tie group


def _dense(labels) -> np.ndarray:
    return labels.toarray() if sp.issparse(labels) else np.asarray(labels)


def _inputs(scores, labels, mask):
    scores = np.asarray(scores)
    if not np.issubdtype(scores.dtype, np.floating):
        scores = scores.astype(np.float64)
    labels = _dense(labels)
    if labels.shape != scores.shape:
        raise ValueError(f"labels {labels.shape} and scores {scores.shape} differ in shape")
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != scores.shape:
            raise ValueError(f"mask {mask.shape} and scores {scores.shape} differ in shape")
    return scores, labels, mask


def _rank(s: np.ndarray, y: np.ndarray, pair: np.ndarray) -> _Ranked:
    """Sort each row by descending score; padding (score -inf) sorts last."""
    order = np.argsort(-s, axis=1, kind="stable")
    s, y, pair = (np.take_along_axis(x, order, axis=1) for x in (s, y, pair))
    m = s.shape[1]
    same = s[:, 1:] == s[:, :-1]
    if not np.any(same & np.isfinite(s[:, 1:])):
        return _Ranked(pair, y, None, None)
    idx = np.broadcast_to(np.arange(m), s.shape)
    starts = np.ones(s.shape, dtype=bool)
    starts[:, 1:] = ~same
    ends = np.ones(s.shape, dtype=bool)
    ends[:, :-1] = ~same
    first = np.maximum.accumulate(np.where(starts, idx, 0), axis=1)
    last = np.minimum.accumulate(np.where(ends, idx, m)[:, ::-1], axis=1)[:, ::-1].copy()
    return _Ranked(pair, y, first, last)


def _by_side_effect(scores, labels, mask) -> _Ranked:
    """One ranked row per side effect, over the triples that count."""
    n, k = scores.shape
    mask_t = np.ones((k, n), dtype=bool) if mask is None else np.ascontiguousarray(mask.T)
    m = int(mask_t.sum(axis=1).max())
    pair = np.full((k, m), n)
    for r in range(k):
        idx = np.flatnonzero(mask_t[r])
        pair[r, : len(idx)] = idx
    cols = np.arange(k)[:, None]
    s = np.vstack([scores, np.full((1, k), -np.inf)])[pair, cols]
    y = np.vstack([labels, np.zeros((1, k), labels.dtype)])[pair, cols].astype(np.float64)
    return _rank(s, y, pair)


def _pooled(scores, labels, mask) -> _Ranked:
    """All counted triples in one ranked row, for the micro averages."""
    mask = np.ones(scores.shape, dtype=bool) if mask is None else mask
    rows = np.nonzero(mask)[0]
    return _rank(scores[mask][None], labels[mask][None].astype(np.float64), rows[None])


class _Weighted:
    """AUPRC and AUROC per ranked row, each triple weighted by its pair's integer weight.

    Torch, so the bootstrap can run on MPS. Weights are summed as integers, so the counts
    stay exact in float32. After the two running sums only the positives matter: each one
    adds its precision (AUPRC) and the zeros it beats (AUROC).
    """

    def __init__(self, ranked: _Ranked, device: str = "cpu"):
        self.device = torch.device(device)
        cpu = self.device.type == "cpu"
        self.float = torch.float64 if cpu else torch.float32
        self.int = np.int64 if cpu else np.int32  # MPS is faster on 32 bits; counts fit
        groups, m = ranked.y.shape
        rows, cols = np.nonzero(ranked.y)
        first = cols if ranked.first is None else ranked.first[rows, cols]
        last = cols if ranked.last is None else ranked.last[rows, cols]
        t = lambda x: torch.as_tensor(x, device=self.device)
        self.pair, self.y = t(ranked.pair), t(ranked.y.astype(self.int))
        self.n_groups, self.group = groups, t(rows)
        # Flat positions of each positive, its tie group's last rank, and the rank before the group.
        self.pos, self.end = t(rows * m + cols), t(rows * m + last)
        self.before, self.has_before = t(rows * m + np.maximum(first - 1, 0)), t(first > 0)

    def __call__(self, w_pair: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        w = torch.as_tensor(np.append(w_pair, 0).astype(self.int), device=self.device)[self.pair]
        cw, ctp = w.cumsum(1), (w * self.y).cumsum(1)
        f = lambda x: x.to(self.float)
        p = f(ctp[:, -1])
        n = f(cw[:, -1]) - p
        cw, ctp = cw.flatten(), ctp.flatten()
        w_pos = f(w.flatten()[self.pos])
        # Within a tie group every triple shares the group's end point.
        tp_end, w_end = f(ctp[self.end]), f(cw[self.end])
        fp_before = torch.where(self.has_before, f(cw[self.before] - ctp[self.before]), 0)
        auprc = self._sum(w_pos * tp_end / w_end.clamp(min=1)) / p
        # A positive beats the zeros ranked below it and half of those tied with it.
        auroc = self._sum(w_pos * (n[self.group] - (w_end - tp_end + fp_before) / 2)) / (p * n)
        return auprc.cpu().double().numpy(), auroc.cpu().double().numpy()

    def _sum(self, x: torch.Tensor) -> torch.Tensor:
        out = torch.zeros(self.n_groups, dtype=self.float, device=self.device)
        return out.index_add_(0, self.group, x)


def _ap_at(ranked: _Ranked, k: int = AP_AT) -> np.ndarray:
    """Decagon's `apk`: precision at each hit in the top k, over min(positives, k)."""
    hits = ranked.y[:, :k]
    precision = hits.cumsum(axis=1) / np.arange(1, hits.shape[1] + 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (hits * precision).sum(axis=1) / np.minimum(ranked.y.sum(axis=1), k)


def evaluate(scores, labels, mask=None) -> dict:
    """Per side effect AUPRC, AUROC, AP@50, base rate and positives; macro and micro averages."""
    scores, labels, mask = _inputs(scores, labels, mask)
    cols = _by_side_effect(scores, labels, mask)
    ones = np.ones(len(scores), dtype=np.int64)
    auprc, auroc = _Weighted(cols)(ones)
    ap50 = _ap_at(cols)
    n_pos = cols.y.sum(axis=1).astype(np.int64)
    n_triples = (cols.pair < len(scores)).sum(axis=1)
    micro_auprc, micro_auroc = _Weighted(_pooled(scores, labels, mask))(ones)
    return {
        "auprc": auprc,
        "auroc": auroc,
        "ap50": ap50,
        "base_rate": n_pos / n_triples,
        "n_pos": n_pos,
        "macro_auprc": float(np.nanmean(auprc)),
        "macro_auroc": float(np.nanmean(auroc)),
        "macro_ap50": float(np.nanmean(ap50)),
        "micro_auprc": float(micro_auprc[0]),
        "micro_auroc": float(micro_auroc[0]),
        "micro_base_rate": float(n_pos.sum() / n_triples.sum()),
    }


def macro_auprc_torch(scores: torch.Tensor, labels: torch.Tensor, mask=None, chunk: int = 128) -> float:
    """Macro AUPRC on the scores' device, for early stopping every epoch.

    Same as `evaluate`'s `macro_auprc` except that ties are broken in arbitrary order, not
    grouped; with float logits ties are rare. `labels` 0/1 and `mask` (bool, or `None` =
    full rows) are tensors of the scores' shape. Side effects in chunks of `chunk` columns.
    """
    n = scores.shape[0]
    precision_at = 1 / torch.arange(1, n + 1, device=scores.device, dtype=torch.float32)[:, None]
    aps = []
    for s in range(0, scores.shape[1], chunk):
        sc, y = scores[:, s : s + chunk], labels[:, s : s + chunk].float()
        if mask is not None:
            # Uncounted triples sort last with label 0, so they add to no positive's precision.
            m = mask[:, s : s + chunk]
            sc, y = sc.masked_fill(~m, -torch.inf), y * m
        y = y.gather(0, sc.argsort(dim=0, descending=True))
        aps.append((y * y.cumsum(0) * precision_at).sum(0) / y.sum(0))
    return float(torch.cat(aps).nanmean())


def uniform_negatives(labels, mask, seed: int) -> np.ndarray:
    """Decagon's 1:1 evaluation set: per side effect, every held-out positive plus as many
    held-out zeros of that side effect, drawn uniformly without replacement.

    Returns a new mask for `evaluate`. Zeros come only from `mask` (never trained on).
    """
    labels = _dense(labels)
    mask = np.asarray(mask, dtype=bool)
    rng = np.random.default_rng(seed)
    out = mask & (labels == 1)
    for r in range(labels.shape[1]):
        zeros = np.flatnonzero(mask[:, r] & (labels[:, r] == 0))
        n_pos = int(out[:, r].sum())
        if len(zeros) < n_pos:
            raise ValueError(f"side effect {r}: {n_pos} held-out positives, {len(zeros)} zeros")
        out[rng.choice(zeros, n_pos, replace=False), r] = True
    return out


def subgroups(pairs: np.ndarray, has_targets: np.ndarray) -> dict[str, np.ndarray]:
    """Row masks: pairs among drugs with targets, among drugs without, and mixed pairs.

    `has_targets` is a boolean per drug, e.g. `np.isin(np.arange(645), data.targets[:, 0])`.
    """
    n_with = has_targets[pairs].sum(axis=1)
    return {"targets": n_with == 2, "no_targets": n_with == 0, "mixed": n_with == 1}


def new_drug_groups(pairs: np.ndarray, new_drugs: np.ndarray) -> dict[str, np.ndarray]:
    """Cold-drug row masks: pairs with one new drug (scored, macro) and with two (micro only)."""
    n_new = np.isin(pairs, new_drugs).sum(axis=1)
    return {"one_new": n_new == 1, "two_new": n_new == 2}


def resample_weights(pairs: np.ndarray, drugs: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """One cluster-bootstrap draw: weight per pair from drugs resampled with replacement.

    `drugs` are the clusters (the test drugs). A drug drawn m times counts m times, so a
    pair weighs the product of its two drugs' counts; drugs outside `drugs` count once.
    """
    n = max(int(pairs.max()), int(drugs.max())) + 1
    counts = np.ones(n)
    counts[drugs] = 0
    counts += np.bincount(rng.choice(drugs, len(drugs)), minlength=n)
    return counts[pairs[:, 0]] * counts[pairs[:, 1]]


def cluster_bootstrap(
    scores,
    labels,
    pairs,
    mask=None,
    drugs=None,
    n_resamples: int = 1000,
    seed: int = 0,
    device: str | None = None,
) -> dict[str, tuple[float, float]]:
    """95% intervals of macro and micro AUPRC / AUROC from resampling drugs.

    `drugs` defaults to every drug in `pairs` (pair-level, random split); for cold-drug
    pass the test drugs, so each one-new-drug pair follows its new drug. `device`
    defaults to MPS when available.

    Normal interval `θ ± 1.96 · sd` over the resamples, not percentiles: a resample
    repeats drugs, which shrinks the sample and lifts AUPRC by more than its spread, so
    percentile intervals, plain or reflected, can miss the point estimate. The spread
    itself is what the drug clusters set.
    """
    scores, labels, mask = _inputs(scores, labels, mask)
    pairs = np.asarray(pairs)
    drugs = np.unique(pairs) if drugs is None else np.asarray(drugs)
    device = device or ("mps" if torch.backends.mps.is_available() else "cpu")
    cols = _Weighted(_by_side_effect(scores, labels, mask), device)
    pooled = _Weighted(_pooled(scores, labels, mask), device)

    def stats(w):
        auprc, auroc = cols(w)
        micro_auprc, micro_auroc = pooled(w)
        return np.nanmean(auprc), np.nanmean(auroc), micro_auprc[0], micro_auroc[0]

    point = np.array(stats(np.ones(len(pairs), dtype=np.int64)))
    rng = np.random.default_rng(seed)
    draws = np.array(
        [stats(resample_weights(pairs, drugs, rng).astype(np.int64)) for _ in range(n_resamples)]
    )
    half = 1.96 * np.nanstd(draws, axis=0)
    lo, hi = point - half, point + half
    names = ("macro_auprc", "macro_auroc", "micro_auprc", "micro_auroc")
    return {name: (float(lo[j]), float(hi[j])) for j, name in enumerate(names)}
