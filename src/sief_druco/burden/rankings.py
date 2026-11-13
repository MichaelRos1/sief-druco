"""Burden of pairs, triples and side effects, and their rankings.

Spec: `docs/spec-burden.md` §8 (and §7 for the minimum count).

    burden_r(a, b) = prev(a, b) × p_r(a, b)            triple (pair × side effect r)
    burden(a, b)   = prev(a, b) × Σ_r p_r(a, b)        pair, all side effects
    burden_r       = Σ_pairs prev(a, b) × p_r(a, b)    side effect, over all covered pairs

`prev` is per 100,000 people, so a burden is expected side-effect signals per 100,000
people. All side effects weigh the same, and `p` is the cross-fitted prediction for
every triple, a labelled one included: the labels are never read here.

Pairs and triples are ranked only for pairs with ≥ `RANK_MIN` participants; side
effects sum over every covered pair. The prevalence interval is carried into burden
with `p` held fixed: a pair's or triple's Korn–Graubard bounds times its `p`, and for
a side effect a Taylor interval of the person-level outcome Σ over the person's covered
pairs of `p_r` (its weighted mean is exactly the sum over pairs above).
"""

import numpy as np
import pandas as pd
from scipy import stats

from sief_druco.burden import survey
from sief_druco.burden.person_drugs import COVERED
from sief_druco.burden.prevalence import PER, RANK_MIN, pair_matrix
from sief_druco.burden.scores import ScoreTable

PAIR_COLUMNS = [
    "unit_a", "unit_b", "participants", "prevalence", "ci_low", "ci_high",
    "labelled", "regime", "sum_p", "burden", "burden_low", "burden_high", "rank",
]
TRIPLE_COLUMNS = [
    "unit_a", "unit_b", "cui", "side_effect", "participants", "pair_rank",
    "p", "burden", "burden_low", "burden_high", "rank",
]
SIDE_EFFECT_COLUMNS = ["cui", "name", "burden", "ci_low", "ci_high", "rank"]


def _rank(df: pd.DataFrame, eligible=None) -> pd.Series:
    """Rank by burden from 1 among `eligible` rows (all if None); NA for the others."""
    eligible = np.ones(len(df), bool) if eligible is None else np.asarray(eligible, bool)
    rank = pd.Series(pd.NA, index=df.index, dtype="Int64")
    order = df.burden[eligible].sort_values(ascending=False, kind="stable").index
    rank[order] = np.arange(1, len(order) + 1)
    return rank


def pair_burden(prevalence: pd.DataFrame, table: ScoreTable, min_participants: int = RANK_MIN) -> pd.DataFrame:
    """Every covered pair of a `prevalence.pair_prevalence` table with its burden
    (`PAIR_COLUMNS`), highest first. `rank` only for pairs with ≥ `min_participants`."""
    c = prevalence[prevalence.status == COVERED].reset_index(drop=True)
    rows = table.rows(c.unit_a, c.unit_b)
    sum_p = np.asarray(table.scores[rows], dtype=np.float64).sum(axis=1)
    out = c.assign(
        labelled=table.pairs.labelled.to_numpy()[rows],
        regime=table.pairs.regime.to_numpy()[rows],
        sum_p=sum_p,
        burden=c.prevalence * sum_p,
        burden_low=c.ci_low * sum_p,
        burden_high=c.ci_high * sum_p,
    )
    out["rank"] = _rank(out, out.participants >= min_participants)
    return out[PAIR_COLUMNS].sort_values("burden", ascending=False, kind="stable").reset_index(drop=True)


def triple_burden(pairs: pd.DataFrame, table: ScoreTable) -> pd.DataFrame:
    """Every side effect of the ranked pairs of `pair_burden` (`TRIPLE_COLUMNS`),
    highest burden first, ranked from 1."""
    ranked = pairs[pairs["rank"].notna()].reset_index(drop=True)
    p = np.asarray(table.scores[table.rows(ranked.unit_a, ranked.unit_b)], dtype=np.float64)
    k = p.shape[1]
    rep = lambda col: np.repeat(ranked[col].to_numpy(), k)  # noqa: E731
    out = pd.DataFrame({
        "unit_a": rep("unit_a"),
        "unit_b": rep("unit_b"),
        "cui": np.tile(table.side_effects.cui.to_numpy(), len(ranked)),
        "side_effect": np.tile(table.side_effects.name.to_numpy(), len(ranked)),
        "participants": rep("participants"),
        "pair_rank": rep("rank"),
        "p": p.ravel(),
        "burden": (ranked.prevalence.to_numpy()[:, None] * p).ravel(),
        "burden_low": (ranked.ci_low.to_numpy()[:, None] * p).ravel(),
        "burden_high": (ranked.ci_high.to_numpy()[:, None] * p).ravel(),
    })
    out = out.sort_values("burden", ascending=False, kind="stable").reset_index(drop=True)
    out["rank"] = np.arange(1, len(out) + 1)
    return out[TRIPLE_COLUMNS]


def side_effect_burden(
    persons: pd.DataFrame, pairs: pd.DataFrame, table: ScoreTable, design=None, level: float = 0.95
) -> pd.DataFrame:
    """Each side effect's burden over all covered pairs, with its Taylor interval,
    highest first, next to the model columns of `side_effects.csv` (held-out AUPRC,
    predicted vs observed positives per regime).

    `persons`, `pairs`: the `persons` table and the person-pairs of
    `person_drugs.pairs` (spec §5).
    """
    design = design or survey.Design.from_persons(persons)
    y, distinct = pair_matrix(persons, pairs)
    covered = (distinct.status == COVERED).to_numpy()
    p = np.asarray(table.scores[table.rows(distinct.unit_a[covered], distinct.unit_b[covered])])
    est = survey.mean(design, y[:, np.flatnonzero(covered)] @ p)  # [persons, side effects]
    half = stats.t.ppf(0.5 + level / 2, est.df) * est.se
    out = table.side_effects.assign(
        burden=PER * est.value, ci_low=PER * (est.value - half), ci_high=PER * (est.value + half)
    )
    out["rank"] = _rank(out)
    model = [c for c in table.side_effects.columns if c not in ("cui", "name")]
    return out[SIDE_EFFECT_COLUMNS + model].sort_values("rank").reset_index(drop=True)
