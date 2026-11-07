"""Survey estimators: Taylor-linearized variance, domains, Korn–Graubard intervals.

Spec: `docs/spec-burden.md` §3 (variance, proportion intervals, software). Checked
once against R's `survey` package; the reference numbers are in
`tests/test_burden_survey.py`.

The design is stratified, with PSUs nested in strata and treated as sampled with
replacement, the setup NCHS recommends for NHANES's masked variance units
(`svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, nest = TRUE)` in R). A subgroup is
a **domain**: its rows get weight 0 in the sums but stay in the design, so every
stratum keeps its full PSU count in the variance. Subsetting the data first would
understate the variance.
"""

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
from scipy import stats


@dataclass(frozen=True)
class Design:
    """Weights, strata and PSUs of one survey (or one pool of stacked cycles)."""

    weight: np.ndarray
    stratum: np.ndarray
    psu: np.ndarray
    # PSU index of each row, and the sparse 0/1 [n_psu, n_rows] that sums rows into PSUs.
    _psu_of_row: np.ndarray
    _rows_to_psu: sp.csr_array
    # Stratum index of each PSU, and each stratum's PSU count in the full design.
    _psu_stratum: np.ndarray
    _n_psu_in_stratum: np.ndarray

    @classmethod
    def from_arrays(cls, weight, stratum, psu) -> "Design":
        weight = np.asarray(weight, dtype=float)
        stratum, psu = np.asarray(stratum), np.asarray(psu)
        if not len(weight) == len(stratum) == len(psu):
            raise ValueError("weight, stratum and psu need one entry per row")
        if np.any(weight < 0) or not np.all(np.isfinite(weight)):
            raise ValueError("weights must be finite and non-negative")

        _, stratum_idx = np.unique(stratum, return_inverse=True)
        # PSUs are nested: PSU 1 of two strata are two PSUs.
        cells, psu_idx = np.unique(
            np.stack([stratum_idx, np.unique(psu, return_inverse=True)[1]], axis=1),
            axis=0,
            return_inverse=True,
        )
        psu_idx = psu_idx.ravel()
        psu_stratum = cells[:, 0]
        n_psu_in_stratum = np.bincount(psu_stratum)
        if np.any(n_psu_in_stratum < 2):
            bad = np.unique(stratum)[n_psu_in_stratum < 2].tolist()
            raise ValueError(f"single PSU in strata {bad}: no variance can be estimated")

        n = len(weight)
        rows_to_psu = sp.csr_array(
            (np.ones(n), (psu_idx, np.arange(n))), shape=(len(cells), n)
        )
        return cls(weight, stratum, psu, psu_idx, rows_to_psu, psu_stratum, n_psu_in_stratum)

    @classmethod
    def from_persons(cls, persons) -> "Design":
        """Design of a `persons` table (spec §5): pooled weight, stratum, PSU."""
        return cls.from_arrays(persons["weight"], persons["stratum"], persons["psu"])

    @property
    def n_psu(self) -> int:
        return len(self._psu_stratum)

    @property
    def df(self) -> int:
        """Design degrees of freedom: PSUs minus strata."""
        return self.n_psu - len(self._n_psu_in_stratum)

    def _domain_df(self, domain: np.ndarray) -> int:
        # As R's degf() on a subset: PSUs minus strata among rows in the domain.
        psus = np.unique(self._psu_of_row[domain & (self.weight > 0)])
        return len(psus) - len(np.unique(self._psu_stratum[psus]))

    def _variance(self, psu_totals: np.ndarray) -> np.ndarray:
        """With-replacement variance from linearized PSU totals [n_psu, k]."""
        n_h = self._n_psu_in_stratum
        stratum_sum = np.zeros((len(n_h), psu_totals.shape[1]))
        np.add.at(stratum_sum, self._psu_stratum, psu_totals)
        dev = psu_totals - (stratum_sum / n_h[:, None])[self._psu_stratum]
        scale = (n_h / (n_h - 1))[self._psu_stratum]
        return (scale[:, None] * dev**2).sum(axis=0)


@dataclass(frozen=True)
class Estimate:
    """A ratio or mean: scalar fields for one outcome, arrays for a matrix of them."""

    value: np.ndarray | float
    variance: np.ndarray | float
    n: int  # rows in the domain
    df: int  # design degrees of freedom of the domain

    @property
    def se(self):
        return np.sqrt(self.variance)

    def korn_graubard(self, level: float = 0.95):
        """Korn–Graubard interval for a proportion, as R's `svyciprop(method = "beta")`.

        Clopper–Pearson at the effective sample size p(1-p)/var, scaled by
        (t_{n-1} / t_df)^2. With zero variance (p = 0 or 1) the effective sample
        size is the domain's row count.
        """
        p, var = np.asarray(self.value, float), np.asarray(self.variance, float)
        alpha = 1 - level
        with np.errstate(divide="ignore", invalid="ignore"):
            n_eff = np.where(var > 0, p * (1 - p) / var, self.n)
        n_eff = n_eff * (stats.t.ppf(alpha / 2, self.n - 1) / stats.t.ppf(alpha / 2, self.df)) ** 2
        with np.errstate(invalid="ignore"):
            lo = np.where(p > 0, stats.beta.ppf(alpha / 2, n_eff * p, n_eff * (1 - p) + 1), 0.0)
            hi = np.where(p < 1, stats.beta.ppf(1 - alpha / 2, n_eff * p + 1, n_eff * (1 - p)), 1.0)
        if np.ndim(self.value) == 0:
            return float(lo), float(hi)
        return lo, hi


def ratio(design: Design, y, x=None, domain=None) -> Estimate:
    """Σ w·y / Σ w·x over the domain, with its Taylor variance.

    `y`: one outcome [n] or several [n, k], dense or scipy sparse.
    `x`: the denominator [n]; None means 1 for every row (a mean).
    `domain`: bool [n]; None means everyone.
    """
    n = len(design.weight)
    domain = np.ones(n, bool) if domain is None else np.asarray(domain, bool)
    x = np.ones(n) if x is None else np.asarray(x, float)
    if len(domain) != n or len(x) != n or y.shape[0] != n:
        raise ValueError("y, x and domain need one entry per row of the design")

    one_column = not sp.issparse(y) and np.ndim(y) == 1
    y = sp.csr_array(y) if sp.issparse(y) else np.asarray(y, float).reshape(n, -1)

    w = design.weight * domain
    x_total = w @ x
    if x_total == 0:
        raise ValueError("the denominator sums to zero in this domain")
    wy = (y.T @ w) if sp.issparse(y) else (w @ y)
    r = np.asarray(wy).ravel() / x_total

    # Linearized values z = w (y - r x) / X, summed per PSU without forming [n, k].
    agg = design._rows_to_psu
    wy_psu = agg @ (sp.diags_array(w) @ y) if sp.issparse(y) else agg @ (w[:, None] * y)
    wy_psu = wy_psu.toarray() if sp.issparse(wy_psu) else wy_psu
    psu_totals = (wy_psu - np.outer(agg @ (w * x), r)) / x_total

    var = design._variance(psu_totals)
    est_df = design.df if domain.all() else design._domain_df(domain)
    if one_column:
        return Estimate(float(r[0]), float(var[0]), int(domain.sum()), est_df)
    return Estimate(r, var, int(domain.sum()), est_df)


def mean(design: Design, y, domain=None) -> Estimate:
    """Weighted mean of `y` over the domain; a proportion when `y` is 0/1."""
    return ratio(design, y, None, domain)
