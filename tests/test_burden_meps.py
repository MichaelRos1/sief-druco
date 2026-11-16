import numpy as np
import pandas as pd
import pytest

from sief_druco import decagon
from sief_druco.burden import coverage as cov
from sief_druco.burden import decagon_table as dt
from sief_druco.burden import handfix, lookup, meps, survey
from sief_druco.burden import ingredients as ing
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden import prevalence as prv
from sief_druco.burden.ingredients import Resolution

# --- Synthetic records: drug names -> resolutions, and IN -> STITCH ID ---

METF, LISIN, HCTZ, ATOR = "6809", "29046", "5487", "83367"
STITCH_OF = {METF: "CID000004091", LISIN: "CID000003961", HCTZ: "CID000003639", ATOR: "CID000060823"}


def _found(*ins):
    return Resolution(tuple((i, f"in{i}") for i in ins), None, ing.NAME)


RESOLUTIONS = {
    "METFORMIN": _found(METF),
    "LISINOPRIL": _found(LISIN),
    "HYDROCHLOROTHIAZIDE-LISINOPRIL": _found(HCTZ, LISIN),
    "ATORVASTATIN": _found(ATOR),
    "ANTIVIRAL AGENTS": Resolution((), ing.NAME_NOT_FOUND, ""),  # a masked drug
    "-9": Resolution((), meps.NO_NAME, ""),
}


def records_of(rows):
    """`rows`: (year, dupersid, round, name, weight)."""
    return pd.DataFrame(rows, columns=meps.RECORD_COLUMNS)


def drugs_of(records):
    return pdr.build(meps.drug_rows(records), lambda code, name: RESOLUTIONS[name], STITCH_OF)


def covered_pairs(records):
    p = pdr.pairs(drugs_of(records))
    return p[p.status == pdr.COVERED]


# --- Years ---


@pytest.mark.parametrize("years", [[2018], [2019], [2017, 2018], [2015, 2016, 2017, 2018, 2019]])
def test_the_2018_and_2019_files_are_refused(years):
    with pytest.raises(ValueError, match="2018 questionnaire redesign"):
        meps.check_years(years)
    with pytest.raises(ValueError, match="2018 questionnaire redesign"):
        meps.load_records(years)  # before any file is read


@pytest.mark.parametrize("years", [[], [2014], [2016, 2016]])
def test_unknown_repeated_or_missing_years_raise(years):
    with pytest.raises(ValueError):
        meps.check_years(years)


def test_the_pool_is_2015_to_2017():
    assert meps.check_years(meps.YEARS) == (2015, 2016, 2017)


# --- Rounds and person-years ---


def test_drugs_bought_in_different_rounds_dont_pair():
    r = records_of([
        (2015, 1, 1, "METFORMIN", 900.0),
        (2015, 1, 2, "LISINOPRIL", 900.0),
    ])
    assert covered_pairs(r).empty


def test_drugs_bought_in_the_same_round_pair_once_per_person_year():
    r = records_of([
        (2015, 1, 1, "METFORMIN", 900.0),
        (2015, 1, 1, "METFORMIN", 900.0),  # a refill
        (2015, 1, 1, "LISINOPRIL", 900.0),
        (2015, 1, 3, "METFORMIN", 900.0),  # the same pair again in round 3
        (2015, 1, 3, "LISINOPRIL", 900.0),
    ])
    p = covered_pairs(r)
    assert len(p) == 1
    assert {p.unit_a.iloc[0], p.unit_b.iloc[0]} == {STITCH_OF[METF], STITCH_OF[LISIN]}


def test_a_combination_bought_alone_pairs_its_ingredients():
    p = covered_pairs(records_of([(2016, 7, 2, "HYDROCHLOROTHIAZIDE-LISINOPRIL", 500.0)]))
    assert {p.unit_a.iloc[0], p.unit_b.iloc[0]} == {STITCH_OF[HCTZ], STITCH_OF[LISIN]}


def test_a_person_in_two_years_counts_as_two_person_years():
    r = records_of([
        (2015, 60001101, 3, "METFORMIN", 3000.0),
        (2015, 60001101, 3, "LISINOPRIL", 3000.0),
        (2016, 60001101, 1, "METFORMIN", 6000.0),
        (2016, 60001101, 1, "LISINOPRIL", 6000.0),
        (2016, 60002201, 2, "ATORVASTATIN", 3000.0),
    ])
    persons = meps.persons(r)
    assert len(persons) == 3
    one = persons[persons.seqn == 60001101].set_index("cycle")
    # Weight PERWTyyF ÷ the two years pooled, each year its own.
    assert one.weight.to_dict() == {"2015": 1500.0, "2016": 3000.0}

    p = covered_pairs(r)
    assert sorted(p.cycle) == ["2015", "2016"]
    prev = prv.point_prevalence(persons, pdr.pairs(drugs_of(r)))
    row = prev[prev.status == pdr.COVERED].iloc[0]
    assert row.participants == 2
    assert row.prevalence == pytest.approx(prv.PER * 4500 / 6000)


def test_the_pool_divides_each_weight_by_three():
    r = records_of([(y, 1, 1, "METFORMIN", 300.0) for y in meps.YEARS])
    assert meps.persons(r).weight.tolist() == [100.0, 100.0, 100.0]


def test_one_person_year_has_one_weight():
    r = records_of([(2015, 1, 1, "METFORMIN", 300.0), (2015, 1, 2, "LISINOPRIL", 301.0)])
    with pytest.raises(ValueError, match="weight differs"):
        meps.persons(r)


# --- Masked and unmapped records ---


def test_masked_and_nameless_records_are_mapping_failed(tmp_path):
    def offline(url, form):
        raise OSError(f"no request expected: {url}")

    resolve = meps.resolver(lookup.Cache(tmp_path / "rxnav.json", fetch=offline))
    assert resolve("", "-9") == Resolution((), meps.NO_NAME, "")
    assert resolve("", "") == Resolution((), meps.NO_NAME, "")

    r = records_of([
        (2017, 5, 1, "METFORMIN", 100.0),
        (2017, 5, 1, "-9", 100.0),
        (2017, 5, 1, "-9", 100.0),
        (2017, 5, 1, "ANTIVIRAL AGENTS", 100.0),
    ])
    d = drugs_of(r)
    failed = d[d.status == pdr.MAPPING_FAILED].set_index("cause")
    assert failed.mentions.to_dict() == {meps.NO_NAME: 2, ing.NAME_NOT_FOUND: 1}
    assert set(pdr.pairs(d).status) == {pdr.MAPPING_FAILED}


def test_counts_records_reaching_an_ingredient():
    r = records_of([
        (2017, 5, 1, "METFORMIN", 100.0),
        (2017, 5, 1, "HYDROCHLOROTHIAZIDE-LISINOPRIL", 100.0),
        (2017, 5, 2, "-9", 100.0),
        (2017, 6, 1, "ANTIVIRAL AGENTS", 100.0),
    ])
    n = meps.counts(r, meps.persons(r), drugs_of(r))
    assert (n["records"], n["records_reaching_in"], n["no_name_records"]) == (4, 2, 1)
    assert (n["person_years"], n["person_years_2_ins"], n["person_years_2_ins_same_round"]) == (2, 1, 1)


# --- No variance or interval ---


@pytest.fixture
def meps_tables():
    r = records_of([
        (2015, i, 1, name, 100.0 + i)
        for i in range(6)
        for name in (["METFORMIN", "LISINOPRIL"] if i % 2 else ["ATORVASTATIN", "METFORMIN"])
    ])
    return meps.persons(r), pdr.pairs(drugs_of(r))


def test_meps_gets_no_variance_or_interval(meps_tables):
    persons, pairs = meps_tables
    assert persons[["stratum", "psu"]].isna().all().all()
    with pytest.raises(ValueError, match="point estimates only"):
        survey.Design.from_persons(persons)
    with pytest.raises(ValueError, match="point estimates only"):
        prv.pair_prevalence(persons, pairs)
    # Not even with a design built some other way.
    design = survey.Design.from_arrays(persons.weight, np.arange(6) % 2, (np.arange(6) // 2) % 2)
    with pytest.raises(ValueError, match="point estimates only"):
        prv.pair_prevalence(persons, pairs, design=design)

    prev = prv.point_prevalence(persons, pairs)
    assert prev[["ci_low", "ci_high"]].isna().all().all()
    assert prev.prevalence.notna().all()


def test_point_prevalence_is_the_survey_point_estimate(meps_tables):
    persons, pairs = meps_tables
    as_nhanes = persons.assign(source="NHANES", stratum=np.arange(6) % 2, psu=(np.arange(6) // 2) % 2)
    a = prv.point_prevalence(persons, pairs)
    b = prv.pair_prevalence(as_nhanes, pairs.assign(source="NHANES"))
    pd.testing.assert_frame_equal(a.drop(columns=["ci_low", "ci_high"]), b.drop(columns=["ci_low", "ci_high"]))


# --- MEPS 2015–2017 (raw files and lookup caches needed) ---

LOCAL_INPUTS = [
    *(meps.RAW_DIR / f / f"{f}.ssp" for f in meps.FILES.values()),
    decagon.RAW_DIR,
    handfix.SIDER_NAMES,
    lookup.CACHE_DIR / meps.RXNAV_CACHE,
    lookup.CACHE_DIR / dt.PUBCHEM_CACHE,
    lookup.CACHE_DIR / dt.RXNAV_CACHE,
]

needs_local = pytest.mark.skipif(
    not all(p.exists() for p in LOCAL_INPUTS),
    reason="raw MEPS / Decagon / SIDER files or the lookup caches not in local/",
)


def _offline(url, form):
    raise OSError(f"tests make no network calls: {url}")


@needs_local
def test_meps_2015_17_matches_the_research_counts(tmp_path):
    out = meps.run(out_dir=tmp_path, fetch=_offline)
    n = out["counts"].set_index("what").n

    # research/meps-prescribed-medicines.md §2.4, §5: 93% of records reach an IN, and
    # 15,033 + 14,256 + 13,782 person-years have ≥ 2 distinct ingredients.
    assert n["records"] == 330_453 + 319_685 + 310_487
    assert n["records_reaching_in_pct"] == pytest.approx(93, abs=0.5)
    assert n["person_years_2_ins"] == 43_071
    assert n["person_years_2_ins_same_round"] == 14_631 + 13_735 + 13_312

    # Every released record has a positive weight, against the documented counts.
    rc = out["record_counts"].set_index("year")
    assert rc.difference.to_dict() == {2015: 5_496, 2016: 4_473, 2017: 4_226}

    c = out["coverage"].set_index(["level", "weighted"])
    assert c[list(pdr.STATUSES)].sum(axis=1).to_numpy() == pytest.approx(100)
    assert (out["failures"].side == cov.SOURCE).all()
    assert out["prevalence"][["ci_low", "ci_high"]].isna().all().all()
