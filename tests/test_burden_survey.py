import numpy as np
import pytest
import scipy.sparse as sp

from sief_druco.burden import nhanes, survey

needs_raw = pytest.mark.skipif(
    not nhanes.RAW_DIR.is_dir(), reason=f"raw NHANES files not in {nhanes.RAW_DIR}"
)


def tiny_design():
    """Two strata of two PSUs, one person each, all weights 1."""
    return survey.Design.from_arrays(
        weight=[1.0, 1.0, 1.0, 1.0], stratum=[1, 1, 2, 2], psu=[1, 2, 1, 2]
    )


def synthetic_survey(seed=0):
    """Six strata of three PSUs, 20 people per PSU, uneven weights."""
    rng = np.random.default_rng(seed)
    stratum = np.repeat(np.arange(6), 60)
    psu = np.tile(np.repeat(np.arange(3), 20), 6)
    weight = rng.uniform(0.5, 3.0, size=360)
    age = rng.integers(0, 90, size=360)
    # The outcome depends on the PSU, so clustering shows in the variance.
    psu_effect = rng.normal(0, 1, size=18)[stratum * 3 + psu]
    y = (rng.normal(0, 1, size=360) + psu_effect > 0.8).astype(float)
    return survey.Design.from_arrays(weight, stratum, psu), y, age


# --- Design ---


def test_design_rejects_a_stratum_with_one_psu():
    with pytest.raises(ValueError, match="single PSU"):
        survey.Design.from_arrays(weight=[1, 1, 1], stratum=[1, 1, 2], psu=[1, 2, 1])


def test_psus_are_nested_in_strata():
    # PSU 1 of stratum 1 and PSU 1 of stratum 2 are different PSUs.
    d = tiny_design()
    assert d.n_psu == 4
    assert d.df == 2


# --- Estimators ---


def test_mean_and_taylor_variance_match_hand_computation():
    # Mean 3/4. Linearized values (y - 3/4) / 4 per PSU: stratum 1 has 1/16 and
    # -3/16 (deviations ±1/8), stratum 2 has 1/16 twice (no deviation).
    # Variance = 2/(2-1) * (2 * (1/8)^2) = 1/16.
    est = survey.mean(tiny_design(), np.array([1.0, 0.0, 1.0, 1.0]))
    assert est.value == pytest.approx(0.75)
    assert est.variance == pytest.approx(1 / 16)
    assert est.n == 4 and est.df == 2


def test_ratio_is_weighted_sum_over_weighted_sum():
    d = survey.Design.from_arrays(weight=[1, 2, 3, 4], stratum=[1, 1, 2, 2], psu=[1, 2, 1, 2])
    y, x = np.array([2.0, 0.0, 1.0, 3.0]), np.array([1.0, 1.0, 2.0, 2.0])
    est = survey.ratio(d, y, x)
    assert est.value == pytest.approx((2 + 0 + 3 + 12) / (1 + 2 + 6 + 8))


def test_columns_of_a_matrix_are_estimated_independently():
    d, y, age = synthetic_survey()
    ys = np.stack([y, 1 - y], axis=1)
    both = survey.mean(d, ys)
    for k, col in enumerate([y, 1 - y]):
        one = survey.mean(d, col)
        assert both.value[k] == pytest.approx(one.value)
        assert both.variance[k] == pytest.approx(one.variance)


def test_sparse_and_dense_outcomes_give_the_same_estimate():
    d, y, age = synthetic_survey()
    ys = np.stack([y, 1 - y], axis=1)
    dense, sparse = survey.mean(d, ys), survey.mean(d, sp.csr_array(ys))
    np.testing.assert_allclose(sparse.value, dense.value)
    np.testing.assert_allclose(sparse.variance, dense.variance)


def test_domain_keeps_all_rows_and_differs_from_the_naive_subset():
    d, y, age = synthetic_survey()
    # Older people, none of whom live in PSU 0 of the first two strata: a PSU
    # with no domain members, as happens for small age bands.
    old = (age >= 65) & ~((d.psu == 0) & (d.stratum < 2))
    domain = survey.mean(d, y, domain=old)

    naive_design = survey.Design.from_arrays(d.weight[old], d.stratum[old], d.psu[old])
    naive = survey.mean(naive_design, y[old])

    # Same point estimate. The domain keeps the empty PSUs (each a zero total),
    # the subset drops them and so changes those strata's PSU counts and means.
    assert domain.value == pytest.approx(naive.value)
    assert domain.variance != pytest.approx(naive.variance, rel=1e-3)
    assert domain.n == old.sum()


def test_domain_counts_strata_without_domain_members():
    d, y, age = synthetic_survey()
    # A domain that sits in stratum 0 only.
    only_first = d.stratum == d.stratum[0]
    est = survey.mean(d, y, domain=only_first)
    assert est.df == 3 - 1  # PSUs - strata among rows in the domain


def test_korn_graubard_brackets_the_estimate_and_is_wider_than_wald():
    d, y, age = synthetic_survey()
    est = survey.mean(d, y)
    lo, hi = est.korn_graubard()
    assert lo < est.value < hi
    # Wider than the Wald interval under a design effect > 1 at these sizes.
    assert hi - lo > 2 * 1.96 * est.se * 0.9


def test_korn_graubard_handles_a_zero_proportion():
    d, y, age = synthetic_survey()
    est = survey.mean(d, np.zeros_like(y))
    lo, hi = est.korn_graubard()
    assert lo == 0.0
    assert 0.0 < hi < 0.05


# --- Reference values from R's `survey` package ---
#
# Computed once (2026-11-07) with R 4.5.3 and survey 4.5 on NHANES 2017-2018,
# run in local/data/raw/nhanes/. R is not part of the pipeline; rerun this only if
# the estimators change. The domain is a subset() of the design, which R treats
# as a domain (it keeps each stratum's PSU count).
#
#   library(survey); library(foreign)
#   demo <- read.xport("DEMO_J.xpt"); rx <- read.xport("RXQ_RX_J.xpt")
#   rx <- rx[!duplicated(rx$SEQN), c("SEQN", "RXDUSE")]
#   d <- merge(demo, rx, by = "SEQN")
#   d$any_rx <- as.numeric(d$RXDUSE == 1)
#   des <- svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTINT2YR,
#                    nest = TRUE, data = d)
#   old <- subset(des, RIDAGEYR >= 65)
#   svymean(~any_rx, des); svyciprop(~any_rx, des, method = "beta")
#   svymean(~any_rx, old); svyciprop(~any_rx, old, method = "beta")
#   svymean(~RIDAGEYR, des); svymean(~RIDAGEYR, old)

R_REFERENCE = {
    # name: (value, se, ci_low, ci_high); ci from svyciprop(method = "beta")
    "any_rx": (0.489569, 0.012519, 0.462564, 0.516618),
    "any_rx_65plus": (0.894260, 0.012312, 0.864850, 0.919161),
    "age": (38.423973, 0.524406, None, None),
    "age_65plus": (73.042543, 0.337875, None, None),
}


@needs_raw
def test_estimates_agree_with_r_survey():
    persons = nhanes.load_persons(["J"])
    d = survey.Design.from_persons(persons)
    any_rx = (persons.rx_use == 1).to_numpy(float)
    age = persons.age.to_numpy(float)
    old = age >= 65

    for name, est, is_prop in [
        ("any_rx", survey.mean(d, any_rx), True),
        ("any_rx_65plus", survey.mean(d, any_rx, domain=old), True),
        ("age", survey.mean(d, age), False),
        ("age_65plus", survey.mean(d, age, domain=old), False),
    ]:
        value, se, lo, hi = R_REFERENCE[name]
        assert est.value == pytest.approx(value, abs=5e-4), name
        assert est.se == pytest.approx(se, abs=5e-4), name
        if is_prop:
            ours = est.korn_graubard()
            assert ours[0] == pytest.approx(lo, abs=5e-4), name
            assert ours[1] == pytest.approx(hi, abs=5e-4), name
