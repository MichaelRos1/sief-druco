import numpy as np
import pandas as pd
import pytest
from scipy import stats

from sief_druco import decagon
from sief_druco.burden import coverage as cov
from sief_druco.burden import decagon_table as dt
from sief_druco.burden import handfix, lookup, nhanes, rq2
from sief_druco.burden import ingredients as ing
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden import prevalence as prv
from sief_druco.burden.ingredients import Resolution

# --- Synthetic drugs (as in test_burden_person_drugs) ---

MPH, DEXMPH, METF, LISIN, ESCIT = "6901", "1091", "6809", "29046", "321988"
STITCH_OF = {MPH: "CID000004158", DEXMPH: "CID000004158", METF: "CID000004091", LISIN: "CID000003961"}


def _found(*ins):
    return Resolution(tuple((i, f"in{i}") for i in ins), None, ing.MMSL)


RESOLUTIONS = {
    "d_mph": _found(MPH),
    "d_dexmph": _found(DEXMPH),
    "d_metf": _found(METF),
    "d_lisin": _found(LISIN),
    "d_combo": _found(METF, LISIN),
    "d_escit": _found(ESCIT),
    "c00049": Resolution((), ing.ACH_CODE, ""),
    "": Resolution((), ing.NO_CODE, ""),
}


def drugs_of(people: dict[int, list[str]], cycle="J"):
    rows = [
        dict(source="NHANES", cycle=cycle, seqn=s, round=pd.NA, code=c, name="55555" if c == "" else c.upper())
        for s, codes in people.items()
        for c in codes
    ]
    return pdr.build(pd.DataFrame(rows), lambda code, name: RESOLUTIONS[code], STITCH_OF)


def persons_of(weights, stratum, psu, cycle="J"):
    n = len(weights)
    return pd.DataFrame({
        "source": "NHANES", "cycle": cycle, "seqn": np.arange(n), "weight": np.asarray(weights, float),
        "stratum": stratum, "psu": psu, "age": 50, "sex": "female", "rx_use": 2,
    })


def survey_of(people: dict[int, list[str]], n: int, weights=None):
    """`n` persons in two strata of two PSUs; `people` gives the drugs of some."""
    weights = np.ones(n) if weights is None else weights
    persons = persons_of(weights, stratum=np.arange(n) % 2, psu=(np.arange(n) // 2) % 2)
    persons.loc[list(people), "rx_use"] = 1
    drugs = drugs_of(people)
    return persons, drugs, pdr.pairs(drugs)


# --- Coverage ---


def test_coverage_columns_add_to_100_in_every_row():
    persons, drugs, pairs = survey_of(
        {
            0: ["d_metf", "d_lisin"],  # covered
            1: ["d_metf", "d_escit"],  # not in the model
            2: ["d_metf", "c00049", ""],  # mapping failed, one uncoded row
            3: ["d_combo", "d_mph"],  # a combination: two covered ingredients
            4: ["d_metf"],
        },
        n=8,
        weights=[1, 2, 3, 4, 5, 6, 7, 8],
    )
    t = cov.table(persons, drugs, pairs)
    assert t.level.tolist() == [lv for lv in cov.LEVELS for _ in (0, 1)]
    assert t[list(pdr.STATUSES)].sum(axis=1).to_numpy() == pytest.approx(100)

    unweighted = t[~t.weighted].set_index("level")
    # metf, lisin, escit, mph, c00049; the uncoded row is no identifiable drug.
    assert unweighted.total[cov.INGREDIENTS] == 5
    assert unweighted.loc[cov.INGREDIENTS, pdr.COVERED] == pytest.approx(60)
    assert unweighted.total[cov.PARTICIPANTS] == 4
    # Pairs: 0 one, 1 one, 2 three, 3 three.
    assert unweighted.total[cov.PERSON_PAIRS] == 8
    assert unweighted.loc[cov.PERSON_PAIRS, pdr.COVERED] == pytest.approx(100 * 4 / 8)

    # Mentions per person 2, 2, 3, 3, 1 with weights 1, 2, 3, 4, 5: 2+4+9+12+5 = 32.
    # Person 2's failed code and uncoded row: 2 mentions × weight 3.
    weighted = t[t.weighted].set_index("level")
    assert weighted.total[cov.MENTIONS] == 32
    assert weighted.loc[cov.MENTIONS, pdr.MAPPING_FAILED] == pytest.approx(100 * 6 / 32)


def test_participants_with_a_failed_drug_are_mapping_failed():
    persons, drugs, _ = survey_of({0: ["d_metf", "d_escit", "c00049"], 1: ["d_mph", "d_dexmph"]}, n=4)
    s = cov.participant_status(drugs).set_index("seqn")
    assert s.status.to_dict() == {0: pdr.MAPPING_FAILED, 1: pdr.COVERED}
    # Two aliases of one Decagon drug are two drugs, though they never pair.
    assert s.n_drugs[1] == 2


def test_decagon_side_failures_are_split_from_source_failures():
    persons, drugs, _ = survey_of({0: ["d_metf", "d_escit", "c00049"], 1: ["d_escit"]}, n=4)
    drugs = cov.mark_decagon_failures(drugs, {ESCIT})
    escit = drugs[drugs.rxcui == ESCIT]
    assert (escit.status == pdr.MAPPING_FAILED).all() and (escit.cause == cov.DECAGON_UNMATCHED).all()

    f = cov.failures(persons, drugs).set_index(["side", "cause"])
    assert f.mentions.to_dict() == {
        (cov.DECAGON, cov.DECAGON_UNMATCHED): 2, (cov.SOURCE, ing.ACH_CODE): 1,
    }
    assert f.mentions_pct.sum() == pytest.approx(100 * 3 / 4)


def test_moved_counts_status_changes_and_alias_mentions():
    rows = pd.DataFrame([
        dict(source="NHANES", cycle="J", seqn=0, round=pd.NA, code=c, name=c.upper())
        for c in ["d_metf", "d_dexmph", "d_escit"]
    ])
    resolve = lambda code, name: RESOLUTIONS[code]  # noqa: E731
    before = pdr.build(rows, resolve, STITCH_OF)
    # A hand fix joins escitalopram to a Decagon drug.
    after = pdr.build(rows, resolve, STITCH_OF | {ESCIT: "CID000009999"})
    m = cov.moved(before, after, alias_rxcuis=[DEXMPH])
    assert m["changed"] == {(pdr.NOT_IN_MODEL, pdr.COVERED): 1}
    assert m["via_alias"] == 1 and m["mentions"] == 3


# --- Prevalence ---


def test_prevalence_and_interval_match_hand_computation_with_non_users():
    # Strata A (PSUs 1, 2) and B (PSUs 1, 2), two persons per PSU, weights
    # 1, 1, 2, 2, 1, 1, 2, 2. Persons 0, 2 and 4 take the pair; 6 and 7 take nothing.
    persons = persons_of([1, 1, 2, 2, 1, 1, 2, 2], stratum=[1, 1, 1, 1, 2, 2, 2, 2], psu=[1, 1, 2, 2, 1, 1, 2, 2])
    drugs = drugs_of({0: ["d_metf", "d_lisin"], 1: ["d_metf"], 2: ["d_combo"], 3: ["d_lisin"], 4: ["d_metf", "d_lisin"]})
    prev = prv.pair_prevalence(persons, pdr.pairs(drugs)).iloc[0]

    # p = (1 + 2 + 1) / 12 = 1/3: the denominator is all 12 of weight, non-users included.
    assert prev.prevalence == pytest.approx(1e5 / 3)
    assert prev.participants == 3
    # Linearized PSU totals w (y - p) / 12: A1 1/36, A2 2/36, B1 1/36, B2 -4/36.
    # Stratum A: 2 * ((1/72)^2 * 2) = 1/1296; stratum B: 2 * ((5/72)^2 * 2) = 25/1296.
    var = 26 / 1296
    p = 1 / 3
    # Korn–Graubard: effective n = p(1-p)/var, scaled by (t_7 / t_2)^2 (8 rows; 4 PSUs - 2 strata).
    n_eff = p * (1 - p) / var * (stats.t.ppf(0.025, 7) / stats.t.ppf(0.025, 2)) ** 2
    lo = stats.beta.ppf(0.025, n_eff * p, n_eff * (1 - p) + 1)
    hi = stats.beta.ppf(0.975, n_eff * p + 1, n_eff * (1 - p))
    assert (prev.ci_low, prev.ci_high) == pytest.approx((1e5 * lo, 1e5 * hi))


def test_a_pair_under_30_participants_is_not_ranked_but_counts_in_coverage():
    people = {i: ["d_metf", "d_lisin"] for i in range(30)} | {30: ["d_mph", "d_lisin"], 31: ["d_mph", "d_lisin"]}
    persons, drugs, pairs = survey_of(people, n=40)
    prev = prv.pair_prevalence(persons, pairs)
    ranked = prv.ranked(prev)
    assert ranked[["unit_a", "unit_b", "rank"]].values.tolist() == [[STITCH_OF[LISIN], STITCH_OF[METF], 1]]

    small = prev[prev.participants == 2]
    assert small.status.tolist() == [pdr.COVERED]
    t = cov.table(persons, drugs, pairs).set_index(["level", "weighted"])
    assert t.total[(cov.PERSON_PAIRS, False)] == 32
    assert t.loc[(cov.PERSON_PAIRS, False), pdr.COVERED] == pytest.approx(100)

    th = prv.thresholds(prev).set_index("min_participants")
    assert th.pairs.to_dict() == {0: 2, 10: 1, 20: 1, 30: 1}
    assert th.weight_pct[30] == pytest.approx(100 * 30 / 32)


def test_pairs_with_an_uncoded_row_get_no_prevalence():
    persons, drugs, pairs = survey_of({0: ["d_metf", ""], 1: ["d_metf", "d_lisin"]}, n=4)
    prev = prv.pair_prevalence(persons, pairs)
    assert len(prev) == 1 and prev.status.tolist() == [pdr.COVERED]


def test_a_cycle_share_is_a_domain_of_its_own_population():
    persons = pd.concat([
        persons_of([1, 1, 1, 1], [1, 1, 2, 2], [1, 2, 1, 2], cycle="I"),
        persons_of([3, 3, 3, 3], [3, 3, 4, 4], [1, 2, 1, 2], cycle="P_"),
    ], ignore_index=True)
    persons["seqn"] = np.arange(8)
    drugs = pd.concat([
        drugs_of({0: ["d_metf", "d_lisin"]}, cycle="I"),
        drugs_of({4: ["d_metf", "d_lisin"], 5: ["d_metf", "d_lisin"]}, cycle="P_"),
    ], ignore_index=True)
    pairs = pdr.pairs(drugs)
    pooled = prv.pair_prevalence(persons, pairs).prevalence[0]
    assert pooled == pytest.approx(1e5 * 7 / 16)
    by = {c: prv.pair_prevalence(persons, pairs, domain=(persons.cycle == c).to_numpy()) for c in ("I", "P_")}
    assert by["I"].prevalence[0] == pytest.approx(1e5 / 4)
    assert by["P_"].prevalence[0] == pytest.approx(1e5 / 2)
    assert by["I"].participants[0] == 1


# --- NHANES 2017–2018 and the main pool (raw files and lookup caches needed) ---

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
def test_2017_18_covered_pairs_match_the_pilot(tmp_path):
    out = rq2.run(["J"], out_dir=tmp_path, fetch=_offline)

    # Without the hand fixes: the pilot's 8,884 distinct covered pairs, less the 28
    # that de-duplicating by STITCH ID merges; ≥ 10: 542, ≥ 30: 114 (spec §7).
    before = out["thresholds_no_overrides"].set_index("min_participants").pairs
    assert before.to_dict() == {0: 8_856, 10: 542, 20: 194, 30: 114}

    # With them, montelukast, rosuvastatin, allopurinol and the others join (#55).
    after = out["thresholds"].set_index("min_participants").pairs
    assert after.to_dict() == {0: 9_609, 10: 585, 20: 207, 30: 118}

    c = out["coverage"].set_index(["level", "weighted"])
    assert c[list(pdr.STATUSES)].sum(axis=1).to_numpy() == pytest.approx(100)
    assert c.total[(cov.MENTIONS, False)] == 15_560
    assert c.total[(cov.PARTICIPANTS, False)] == 2_941
    # No source drug is one of the 3 Decagon drugs left unmatched.
    assert (out["failures"].side == cov.SOURCE).all()
