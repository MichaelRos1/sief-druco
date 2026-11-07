import pytest

from sief_druco.burden import nhanes, survey

needs_raw = pytest.mark.skipif(
    not nhanes.RAW_DIR.is_dir(), reason=f"raw NHANES files not in {nhanes.RAW_DIR}"
)


# --- Weight factors (no data needed) ---


def test_main_pool_weights_by_years():
    assert nhanes.weight_factors(["I", "P_"]) == pytest.approx({"I": 2 / 5.2, "P_": 3.2 / 5.2})


def test_faers_era_pool_halves_the_two_year_weights():
    assert nhanes.weight_factors(["D", "E"]) == pytest.approx({"D": 0.5, "E": 0.5})


@pytest.mark.parametrize("cycle", ["I", "J", "P_"])
def test_a_single_cycle_keeps_its_weights(cycle):
    assert nhanes.weight_factors([cycle]) == {cycle: 1.0}


def test_2017_18_next_to_the_pre_pandemic_file_raises():
    with pytest.raises(ValueError, match="P_ already holds"):
        nhanes.weight_factors(["J", "P_"])
    with pytest.raises(ValueError, match="P_ already holds"):
        nhanes.load_persons(["I", "J", "P_"])


@pytest.mark.parametrize("cycles", [["A"], ["B", "C"], ["I", "I"], ["L"], []])
def test_unknown_repeated_or_missing_cycles_raise(cycles):
    with pytest.raises(ValueError):
        nhanes.weight_factors(cycles)


# --- Loaded pools (raw files needed) ---


@pytest.fixture(scope="module")
def main_pool():
    return nhanes.load_persons(nhanes.MAIN_POOL)


@needs_raw
def test_main_pool_holds_every_participant(main_pool):
    assert len(main_pool) == 25_531
    assert (main_pool.rx_use == 1).sum() == 10_873
    # Refused / don't know stay in the denominators.
    assert main_pool.rx_use.isin([7, 9]).sum() == 23
    assert set(main_pool.columns) >= {
        "source", "cycle", "seqn", "weight", "stratum", "psu", "age", "sex", "rx_use"
    }
    assert (main_pool.source == "NHANES").all()
    assert not main_pool.duplicated(["source", "cycle", "seqn"]).any()


@needs_raw
def test_main_pool_strata_do_not_overlap_and_have_two_psus(main_pool):
    strata = {c: set(g.stratum) for c, g in main_pool.groupby("cycle")}
    assert strata["I"] == set(range(119, 134))
    assert strata["P_"] == set(range(149, 173))
    psus = main_pool.groupby("stratum").psu.nunique()
    assert psus.min() >= 2
    survey.Design.from_persons(main_pool)  # would raise on a single-PSU stratum


@needs_raw
def test_main_pool_weights_sum_to_one_population(main_pool):
    # Each cycle's pooled weights carry its share of the 5.2 years.
    by_cycle = main_pool.groupby("cycle").weight.sum()
    total = by_cycle.sum()
    assert 300e6 < total < 340e6
    assert by_cycle["I"] / total == pytest.approx(2 / 5.2, abs=0.02)


@needs_raw
def test_faers_era_pool_strata_do_not_overlap():
    persons = nhanes.load_persons(nhanes.FAERS_ERA_POOL)
    strata = {c: set(g.stratum) for c, g in persons.groupby("cycle")}
    assert not strata["D"] & strata["E"]
    assert strata["D"] | strata["E"] == set(range(44, 75))
    survey.Design.from_persons(persons)
