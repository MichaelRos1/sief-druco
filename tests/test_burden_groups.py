import numpy as np
import pandas as pd
import pytest
from scipy import stats

from sief_druco.burden import groups as gr
from sief_druco import decagon
from sief_druco.burden import decagon_table as dt
from sief_druco.burden import handfix, lookup, nhanes, rq2
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden import scores as sc
from sief_druco.burden import survey
from sief_druco.burden.person_drugs import COVERED, MAPPING_FAILED, NOT_IN_MODEL
from tiny_score_table import A, B, C, D, write_table

# --- 16 persons in two strata of two PSUs, with drug lists in the `person_drugs` format ---
#
# (age, stratum, PSU, drugs). A drug is a STITCH ID (covered), "in:<n>" (not in the
# model) or a key that never reached an `IN`: "code:<class code>" or an uncoded row
# "row:<n>". A tuple is one combination pill: one source record, two ingredients.
# Every band has persons in both PSUs of some stratum, so it has a design df ≥ 1.
OTHERS = [f"in:{900 + i}" for i in range(6)]
PEOPLE = [
    (70, 0, 0, [A, B]),
    (70, 1, 0, [A, B, C]),
    (50, 0, 0, [A, B, "in:999", "code:a123", "row:7"]),  # 5 drugs, only AB covered
    (50, 1, 0, [C, D]),
    (10, 0, 0, []),                                    # non-user
    (80, 1, 1, [A]),                                   # one drug
    (50, 0, 1, [(A, C)]),                              # a combination pill: 2 ingredients, pair AC
    (20, 1, 1, []),
    (65, 0, 1, [A, B, C, D] + OTHERS),                 # 10 drugs
    (40, 0, 1, []),
    (12, 0, 1, [A, D]),
    (30, 1, 0, [B, C, D]),
    (68, 0, 0, [A, B, C, D] + OTHERS),
    (55, 1, 1, [A, B, C, D, "in:999"]),
    (75, 1, 1, [A, C, D]),
    (60, 0, 1, [A, B, C, "in:998", "in:999"]),
]
WEIGHTS = np.arange(1, len(PEOPLE) + 1, dtype=float)
AGE = np.array([p[0] for p in PEOPLE])
INGREDIENT_OF = {A: "1", B: "2", C: "3", D: "4"}


def person_drugs_rows():
    rows = []
    for seqn, (*_, drugs) in enumerate(PEOPLE):
        for record, drug in enumerate(drugs):
            for d in drug if isinstance(drug, tuple) else (drug,):
                if d.startswith("CID"):
                    row = dict(drug=f"in:{INGREDIENT_OF[d]}", stitch_id=d, status=COVERED, cause=None)
                elif d.startswith("in:"):
                    row = dict(drug=d, stitch_id=None, status=NOT_IN_MODEL, cause=None)
                else:
                    row = dict(drug=d, stitch_id=None, status=MAPPING_FAILED, cause="class code")
                rows.append(dict(source="NHANES", cycle="J", seqn=seqn, round=pd.NA, code=f"c{record}", **row))
    t = pd.DataFrame(rows)
    t["round"] = t["round"].astype("Int64")
    return t


def tiny():
    persons = pd.DataFrame({
        "source": "NHANES", "cycle": "J", "seqn": np.arange(len(PEOPLE)), "weight": WEIGHTS,
        "stratum": [p[1] for p in PEOPLE], "psu": [p[2] for p in PEOPLE], "age": AGE,
        "rx_use": [1 if p[3] else 2 for p in PEOPLE],
    })
    drugs = person_drugs_rows()
    return persons, drugs, pdr.pairs(drugs)


@pytest.fixture
def table(tmp_path):
    return sc.read(write_table(tmp_path / "handoff"))


# Σ_r p of each pair in the tiny score table, and each person's covered pairs by hand.
SUM_P = {(A, B): 0.6, (A, C): 1.0, (A, D): 2.7, (B, C): 0.6, (B, D): 1.2, (C, D): 0.1}


def covered_pairs(drugs):
    flat = sorted({d for x in drugs for d in (x if isinstance(x, tuple) else (x,)) if d.startswith("CID")})
    return [(a, b) for i, a in enumerate(flat) for b in flat[i + 1:]]


PERSON_PAIRS = [covered_pairs(p[3]) for p in PEOPLE]
PERSON_BURDEN = np.array([sum(SUM_P[p] for p in pairs) for pairs in PERSON_PAIRS])


# --- Bands ---


def test_drug_count_counts_ingredients_not_products():
    persons, drugs, _ = tiny()
    # Drugs outside the 645 count, one per drug that never reached an `IN`, and the
    # combination pill of person 6 counts as two ingredients.
    assert gr.drug_count(persons, drugs).tolist() == [2, 3, 5, 2, 0, 1, 2, 0, 10, 0, 2, 3, 10, 5, 3, 5]


def test_people_who_refused_or_did_not_know_are_in_the_0_1_group(table):
    persons, drugs, pairs = tiny()
    # RXDUSE 7 / 9: no drug rows, but they stay in every denominator (spec §3).
    refused = pd.DataFrame({
        "source": "NHANES", "cycle": "J", "seqn": [16, 17], "weight": [3.0, 4.0],
        "stratum": [0, 1], "psu": [0, 1], "age": [40, 70], "rx_use": [7, 9],
    })
    persons = pd.concat([persons, refused], ignore_index=True)
    assert gr.drug_count(persons, drugs)[-2:].tolist() == [0, 0]
    b = gr.bands(persons, drugs, pairs, table).set_index(["band_type", "band"])
    total = WEIGHTS.sum() + 7
    assert b.loc[("drugs", "0–1"), "population_pct"] == pytest.approx(100 * (5 + 6 + 8 + 10 + 7) / total)
    assert b.loc[("drugs", "0–1"), "persons"] == 6
    old = AGE >= 65
    assert b.loc[("age", "65+"), "burden_per_person"] == pytest.approx(
        1e5 * (WEIGHTS * PERSON_BURDEN)[old].sum() / (WEIGHTS[old].sum() + 4)
    )


def test_band_edges():
    assert gr.age_band(np.array([0, 17, 18, 44, 45, 64, 65, 80])).tolist() == [
        "0–17", "0–17", "18–44", "18–44", "45–64", "45–64", "65+", "65+",
    ]
    assert gr.drug_band(np.array([0, 1, 2, 3, 4, 5, 9, 10, 25])).tolist() == [
        "0–1", "0–1", "2", "3–4", "3–4", "5–9", "5–9", "10+", "10+",
    ]


def test_person_burden_sums_their_covered_pairs(table):
    persons, _, pairs = tiny()
    burden, n_pairs = gr.person_burden(persons, pairs, table)
    assert burden == pytest.approx(PERSON_BURDEN)
    assert n_pairs.tolist() == [len(p) for p in PERSON_PAIRS]


# --- (a) shares and (b) burden per person ---


def bands(table, **kw):
    persons, drugs, pairs = tiny()
    return gr.bands(persons, drugs, pairs, table, **kw).set_index(["band_type", "band"])


def test_shares_of_total_burden_add_to_100(table):
    b = bands(table)
    for band_type, labels in [("age", gr.AGE_BANDS), ("drugs", gr.DRUG_BANDS)]:
        assert b.loc[band_type].index.tolist() == list(labels)
        assert b.loc[band_type].burden_pct.sum() == pytest.approx(100)
        assert b.loc[band_type].population_pct.sum() == pytest.approx(100)
    # The 0–1 drug group holds no burden and appears only as a population share.
    zero = b.loc[("drugs", "0–1")]
    assert zero.burden_pct == 0
    assert zero[["burden_per_person", "ci_low", "ci_high", "covered_pairs_per_person"]].isna().all()
    assert zero.population_pct == pytest.approx(100 * (5 + 6 + 8 + 10) / WEIGHTS.sum())


def test_share_and_burden_per_person_by_hand(table):
    b = bands(table).loc[("age", "65+")]
    w, old = WEIGHTS, AGE >= 65
    assert b.burden_pct == pytest.approx(100 * (w * PERSON_BURDEN)[old].sum() / (w * PERSON_BURDEN).sum())
    # Burden per person uses the band's own population, per 100,000 band members.
    assert b.burden_per_person == pytest.approx(1e5 * (w * PERSON_BURDEN)[old].sum() / w[old].sum())
    n_pairs = np.array([len(p) for p in PERSON_PAIRS])
    assert b.covered_pairs_per_person == pytest.approx((w * n_pairs)[old].sum() / w[old].sum())
    assert b.ci_low < b.burden_per_person < b.ci_high


def test_bands_are_survey_domains_not_subsets(table):
    persons, _, pairs = tiny()
    design = survey.Design.from_persons(persons)
    burden, _ = gr.person_burden(persons, pairs, table)
    domain = (AGE >= 18) & (AGE <= 44)
    est = survey.mean(design, burden, domain=domain)
    b = bands(table).loc[("age", "18–44")]
    assert b.burden_per_person == pytest.approx(1e5 * est.value)
    assert b.ci_high - b.ci_low == pytest.approx(2e5 * est.se * stats.t.ppf(0.975, est.df))
    assert est.se > 0

    # Subsetting the data first loses the band's empty PSU in stratum 0.
    with pytest.raises(ValueError, match="single PSU"):
        survey.Design.from_persons(persons[domain])


# --- (c) top pairs within a band ---


def tops(table, **kw):
    persons, drugs, pairs = tiny()
    return gr.band_tops(persons, drugs, pairs, table, **kw)


def test_band_prevalence_and_burden_use_the_bands_own_population(table):
    top = tops(table, min_participants=1)
    old = top[(top.band_type == "age") & (top.band == "65+")].set_index(["unit_a", "unit_b"])
    takers = np.array([(A, B) in p for p in PERSON_PAIRS]) & (AGE >= 65)
    prev = 1e5 * WEIGHTS[takers].sum() / WEIGHTS[AGE >= 65].sum()
    assert old.loc[(A, B), "prevalence"] == pytest.approx(prev)
    assert old.loc[(A, B), "participants"] == takers.sum() == 4
    assert old.loc[(A, B), "burden"] == pytest.approx(prev * SUM_P[A, B])
    assert (old.ci_low <= old.prevalence).all() and (old.prevalence <= old.ci_high).all()
    # Ranked within the band, with the overall rank and "known to Decagon" next to it.
    assert old.band_rank.tolist() == list(range(1, len(old) + 1))
    assert old.burden.is_monotonic_decreasing
    assert old.loc[(A, B), "labelled"] == 1 and old.loc[(A, C), "labelled"] == 0
    assert old.overall_rank.notna().all()


def test_band_burdens_add_up_to_the_bands_burden_per_person(table):
    top = tops(table, min_participants=1, top=100)
    b = bands(table)
    for (band_type, band), g in top.groupby(["band_type", "band"]):
        assert g.burden.sum() == pytest.approx(b.loc[(band_type, band), "burden_per_person"])


def test_the_rank_rule_applies_within_the_band(table):
    # Overall, AB has 7 participants; within 65+ only 4, within 45–64 only 3.
    top = tops(table, min_participants=4)
    assert set(zip(top.band_type, top.band, top.unit_a, top.unit_b)) == {("age", "65+", A, B), ("age", "65+", A, C)}
    notes = bands(table, min_participants=4).note
    assert notes[("age", "65+")] == gr.FEWER == "fewer than 10 rankable pairs"
    assert notes[("age", "45–64")] == gr.TOO_FEW == "too few participants to rank"
    assert notes[("drugs", "10+")] == gr.TOO_FEW
    assert notes[("drugs", "0–1")] == ""


def test_a_band_with_enough_rankable_pairs_has_no_note(table):
    assert bands(table, min_participants=1, top=2).note[("age", "65+")] == ""
    assert len(tops(table, min_participants=1, top=2).query("band == '65+'")) == 2


# --- NHANES 2017–2018 (raw files and lookup caches needed) ---

LOCAL_INPUTS = [
    nhanes.RAW_DIR,
    decagon.RAW_DIR,
    handfix.SIDER_NAMES,
    lookup.CACHE_DIR / pdr.NHANES_RXNAV_CACHE,
    lookup.CACHE_DIR / dt.PUBCHEM_CACHE,
    lookup.CACHE_DIR / dt.RXNAV_CACHE,
]

needs_local = pytest.mark.skipif(
    not all(p.exists() for p in LOCAL_INPUTS),
    reason="raw NHANES / Decagon / SIDER files or the lookup caches not in local/",
)


def _offline(url, form):
    raise OSError(f"tests make no network calls: {url}")



@needs_local
def test_2017_18_older_and_polypharmacy_bands_hold_most_covered_pairs():
    e = rq2.exposure(["J"], overrides_path=None, fetch=_offline)
    share = gr.person_pair_shares(e.persons, e.person_drugs, e.pairs)
    # #25's pilot (no hand fixes) counted 29,977 covered person-pairs by band.
    pilot = {
        ("age", "0–17"): 569, ("age", "18–44"): 2_014, ("age", "45–64"): 10_699, ("age", "65+"): 16_695,
        ("drugs", "0–1"): 0, ("drugs", "2"): 480, ("drugs", "3–4"): 2_982, ("drugs", "5–9"): 13_455,
        ("drugs", "10+"): 13_060,
    }
    assert share.to_dict() == pytest.approx({k: 100 * v / 29_977 for k, v in pilot.items()})
    # So aged 45+ hold about 92% of them, people on 5+ drugs about 89% (spec §9).
    assert share[("age", "45–64")] + share[("age", "65+")] == pytest.approx(92, abs=1)
    assert share[("drugs", "5–9")] + share[("drugs", "10+")] == pytest.approx(89, abs=1)
