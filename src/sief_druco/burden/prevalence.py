"""Pair prevalence per 100,000 people, and the minimum count for ranking.

Spec: `docs/spec-burden.md` §7 (and §1: shares of the whole pooled population).

Prevalence of a pair = Σ pooled weights of persons taking both ÷ Σ all pooled
weights, so the denominator is everyone in `persons`, non-users included. Its 95%
interval is Korn–Graubard from the Taylor variance (`survey`). A subgroup (one
cycle, an age band) is a survey domain: the same estimator with `domain`, so the
share is of the domain's own population.

A pair is **ranked** only with ≥ `RANK_MIN` unweighted participants; every covered
pair still counts everywhere else (coverage, totals, burden shares).

Pairs holding an NHANES unknown / refused / don't know row are no identifiable
pair and get no prevalence; they still count as person-pairs in coverage (§6).
"""

import numpy as np
import pandas as pd
import scipy.sparse as sp

from sief_druco.burden import survey
from sief_druco.burden.person_drugs import COVERED, PERSON, UNCODED_PREFIX

PER = 100_000
RANK_MIN = 30
THRESHOLDS = (10, 20, 30)

COLUMNS = ["unit_a", "unit_b", "status", "participants", "prevalence", "ci_low", "ci_high"]


def pair_matrix(persons: pd.DataFrame, pairs: pd.DataFrame) -> tuple[sp.csr_array, pd.DataFrame]:
    """Sparse 0/1 [persons, distinct pairs] and the distinct pairs (`unit_a`, `unit_b`, `status`)."""
    keep = ~(pairs.unit_a.str.startswith(UNCODED_PREFIX) | pairs.unit_b.str.startswith(UNCODED_PREFIX))
    p = pairs[keep]
    row_of = pd.Series(np.arange(len(persons)), index=pd.MultiIndex.from_frame(persons[PERSON]))
    rows = row_of.reindex(pd.MultiIndex.from_frame(p[PERSON])).to_numpy()
    if np.isnan(rows.astype(float)).any():
        raise ValueError("pairs of persons missing from the persons table")
    g = p.groupby(["unit_a", "unit_b"], sort=True)
    codes = g.ngroup().to_numpy()
    out = g.status.agg(["first", "nunique"]).reset_index()
    if (out["nunique"] > 1).any():
        raise ValueError("a pair has more than one status")
    out = out.drop(columns="nunique").rename(columns={"first": "status"})
    y = sp.csr_array((np.ones(len(p)), (rows.astype(int), codes)), shape=(len(persons), len(out)))
    if y.nnz and y.max() > 1:
        raise ValueError("a person holds the same pair twice")
    return y, out


def pair_prevalence(persons: pd.DataFrame, pairs: pd.DataFrame, domain=None, design=None) -> pd.DataFrame:
    """One row per distinct pair (`COLUMNS`), most prevalent first.

    `participants`: unweighted persons taking the pair (within `domain`).
    `prevalence`, `ci_low`, `ci_high`: per 100,000 people (of the domain).
    """
    design = design or survey.Design.from_persons(persons)
    y, out = pair_matrix(persons, pairs)
    dom = np.ones(len(persons), bool) if domain is None else np.asarray(domain, bool)
    est = survey.ratio(design, y, domain=dom)
    lo, hi = est.korn_graubard()
    out["participants"] = np.asarray(y.T @ dom.astype(float)).astype(int)
    out["prevalence"] = PER * np.asarray(est.value)
    out["ci_low"], out["ci_high"] = PER * lo, PER * hi
    return out[COLUMNS].sort_values(["prevalence", "unit_a", "unit_b"], ascending=[False, True, True]).reset_index(drop=True)


def ranked(prevalence: pd.DataFrame, min_participants: int = RANK_MIN) -> pd.DataFrame:
    """Covered pairs with ≥ `min_participants`, by prevalence, with a `rank` from 1."""
    r = prevalence[(prevalence.status == COVERED) & (prevalence.participants >= min_participants)]
    r = r.sort_values("prevalence", ascending=False, kind="stable").reset_index(drop=True)
    return r.assign(rank=np.arange(1, len(r) + 1))


def thresholds(prevalence: pd.DataFrame, levels=THRESHOLDS) -> pd.DataFrame:
    """Covered pairs passing each participant count, and their share of covered pair weight (%)."""
    c = prevalence[prevalence.status == COVERED]
    total = c.prevalence.sum()
    rows = [[0, len(c), 100.0]] + [
        [k, int((c.participants >= k).sum()), 100 * c.prevalence[c.participants >= k].sum() / total]
        for k in levels
    ]
    return pd.DataFrame(rows, columns=["min_participants", "pairs", "weight_pct"])
