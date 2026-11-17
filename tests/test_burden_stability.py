from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from sief_druco import decagon
from sief_druco.burden import decagon_table as dt
from sief_druco.burden import handfix, lookup, meps, nhanes, rq2
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden import scores as sc
from sief_druco.burden import stability as st
from sief_druco.burden.ingredients import ACH_CODE, Resolver
from tiny_score_table import A, B, C, D, write_table

# --- Overlap and the three drop-out reasons, on hand-made rankings ---


def ranking(rows):
    """A `rankings.pair_burden`-like table from (pair, participants, rank) rows."""
    return pd.DataFrame(
        [(a, b, n, r) for (a, b), n, r in rows], columns=["unit_a", "unit_b", "participants", "rank"]
    ).astype({"rank": "Int64"})


P = {i: (f"a{i}", f"b{i}") for i in range(1, 11)}

# Main top 5: pairs 1–5; pair 6 is ranked 6th, outside it.
MAIN = ranking([(P[i], 100, i) for i in range(1, 7)])
RERUN = ranking([
    (P[1], 80, 3),  # stays in the top 5
    # P[2]: nobody takes it
    (P[3], 0, pd.NA),  # in the table (a domain lists every pair), but nobody takes it
    (P[4], 12, pd.NA),  # taken, but below 30 participants
    (P[5], 40, 7),  # ranked, but lower than 5th
    (P[6], 90, 1),  # enters the top 5
    (P[7], 70, 2), (P[8], 60, 4), (P[9], 50, 5), (P[10], 45, 6),
])


def test_overlap_and_each_drop_out_reason_match_the_synthetic_case():
    c = st.compare(MAIN, RERUN, k=5)
    assert list(c.main_rank) == [1, 2, 3, 4, 5]
    assert list(c.reason) == ["", st.NOT_OBSERVED, st.NOT_OBSERVED, st.BELOW_MIN, st.RANKED_LOWER]
    assert list(c.in_top) == [True, False, False, False, False]
    assert list(c.rerun_participants) == [80, 0, 0, 12, 40]

    o = st.overlap(RERUN, c, k=5)
    assert o == {
        "ranked": 7, "top": 5, "overlap": 1, st.NOT_OBSERVED: 2, st.BELOW_MIN: 1, st.RANKED_LOWER: 1,
    }


def test_a_rerun_with_fewer_ranked_pairs_than_the_top_still_counts_out_of_k():
    rerun = ranking([(P[1], 80, 1), (P[2], 40, 2)])
    overlaps, dropouts = st.stability(MAIN, {"small": rerun}, k=5)
    row = overlaps.iloc[0]
    assert (row.ranked, row.top, row.overlap) == (2, 2, 2)
    assert dropouts.reason.tolist() == ["", "", st.NOT_OBSERVED, st.NOT_OBSERVED, st.NOT_OBSERVED]


def test_the_below_30_reason_uses_the_ranking_minimum():
    # 29 participants: below 30; a pair with 30 is ranked and can only be "ranked lower".
    main = ranking([(P[1], 100, 1), (P[2], 100, 2)])
    rerun = ranking([(P[3], 99, 1), (P[1], 29, pd.NA), (P[2], 30, 2)])
    assert st.compare(main, rerun, k=1).reason.tolist() == [st.BELOW_MIN]
    assert st.compare(main.assign(rank=[2, 1]).astype({"rank": "Int64"}), rerun, k=1).reason.tolist() == [
        st.RANKED_LOWER
    ]


# --- A cycle re-run through the same code as the main result ---


def _cycle_survey():
    """Cycle I: 40 take AB, 35 take AC. Cycle P_: 40 take AB, 31 take CD. 10 non-users each."""
    takes = {"I": [(A, B)] * 40 + [(A, C)] * 35 + [None] * 10, "P_": [(A, B)] * 40 + [(C, D)] * 31 + [None] * 10}
    persons, pairs = [], []
    for c, (cycle, rows) in enumerate(takes.items()):
        for s, pair in enumerate(rows):
            persons.append(dict(source="NHANES", cycle=cycle, seqn=s, weight=1.0 + s % 3,
                                stratum=10 * c + s % 2, psu=(s // 2) % 2))
            if pair:
                pairs.append(dict(source="NHANES", cycle=cycle, seqn=s, unit_a=pair[0], unit_b=pair[1],
                                  status=pdr.COVERED))
    return SimpleNamespace(persons=pd.DataFrame(persons), pairs=pd.DataFrame(pairs))


def test_each_cycle_ranks_its_own_pairs_as_a_domain(tmp_path):
    table = sc.read(write_table(tmp_path / "handoff"))
    e = _cycle_survey()
    main = st.burden(e, table)
    reruns = {
        name: st.burden(e, table, domain=(e.persons.cycle == cycle).to_numpy())
        for name, cycle in st.CYCLE_OF.items()
    }
    # Pooled, AB (80), AC (35) and CD (31) are ranked; AC scores higher than CD.
    assert main.dropna(subset=["rank"])[["unit_a", "unit_b"]].apply(tuple, axis=1).tolist() == [
        (A, B), (A, C), (C, D)
    ]
    overlaps, dropouts = st.stability(main, reruns, k=2)
    d = dropouts.set_index(["rerun", "unit_a", "unit_b"])
    assert d.reason[(st.CYCLE_I, A, C)] == ""
    # Nobody in P_ takes AC: listed in the domain's table, but with no participant.
    assert d.reason[(st.CYCLE_P, A, C)] == st.NOT_OBSERVED
    assert d.rerun_participants[(st.CYCLE_P, A, B)] == 40
    assert overlaps.set_index("rerun").overlap.to_dict() == {st.CYCLE_I: 2, st.CYCLE_P: 1}
    # A domain's prevalence is a share of the cycle's own population: AB is 40 of 85 there.
    w = e.persons.weight[e.persons.cycle == "P_"]
    ab = reruns[st.CYCLE_P].set_index(["unit_a", "unit_b"]).loc[(A, B)]
    assert ab.prevalence == pytest.approx(1e5 * w[:40].sum() / w.sum())


# --- FAERS-era codes ---


def _offline(url, form=None):
    raise OSError(f"tests make no network calls: {url}")


def test_estriol_class_code_stays_mapping_failed(tmp_path):
    # Estriol is `a70248` in 2005–06 (a class code), `d05027` in the main pool. Even
    # with estriol's `IN` in the model, the class code reaches no `IN`: no source-side
    # hand fixes (§4).
    estriol_in = "4099"
    rows = pd.DataFrame({
        "source": "NHANES", "cycle": "D", "seqn": [1, 1], "round": pd.array([pd.NA] * 2, dtype="Int64"),
        "code": [st.ESTRIOL_CODE, ""], "name": ["ESTRIOL", "99999"],
    })
    resolve = Resolver(lookup.Cache(tmp_path / "rxnav.json", fetch=_offline)).multum
    drugs = pdr.build(rows, resolve, {estriol_in: "CID000005756"})
    got = st.code_status(drugs, st.ESTRIOL_CODE)
    assert got.status.tolist() == [pdr.MAPPING_FAILED]
    assert got.cause.tolist() == [ACH_CODE]


def test_new_codes_counts_codes_missing_from_the_main_pool():
    rows = pd.DataFrame({"code": ["d1", "d1", "d2", "a9", "", "d3"]})
    main_rows = pd.DataFrame({"code": ["d1", "d3", ""]})
    assert st.new_codes(rows, main_rows) == {
        "new_codes": 2, "new_codes_ach": 1, "new_code_mentions": 2,
        "coded_mentions": 5, "new_code_mentions_pct": 40.0,
    }


# --- The FAERS-era pool (raw files and lookup caches needed) ---

FAERS_INPUTS = [
    *(nhanes.RAW_DIR / f"{kind}_{c}.xpt" for kind in ("DEMO", "RXQ_RX") for c in nhanes.FAERS_ERA_POOL),
    decagon.RAW_DIR,
    handfix.SIDER_NAMES,
    lookup.CACHE_DIR / pdr.NHANES_RXNAV_CACHE,
    lookup.CACHE_DIR / dt.PUBCHEM_CACHE,
    lookup.CACHE_DIR / dt.RXNAV_CACHE,
]
needs_faers = pytest.mark.skipif(
    not all(p.exists() for p in FAERS_INPUTS),
    reason="raw NHANES 2005-2008 / Decagon / SIDER files or the lookup caches not in local/",
)


@needs_faers
def test_faers_era_new_codes_match_the_research_counts():
    # #22: 249 codes, 815 mentions = 3.1% of coded mentions, not in the main pool.
    codes = st.new_codes(nhanes.load_drug_rows(nhanes.FAERS_ERA_POOL), nhanes.load_drug_rows(nhanes.MAIN_POOL))
    assert codes["new_codes"] == 249
    assert codes["new_code_mentions"] == 815
    assert codes["new_code_mentions_pct"] == pytest.approx(3.1, abs=0.05)


@needs_faers
def test_the_faers_era_run_keeps_estriol_as_mapping_failed():
    f = rq2.exposure(nhanes.FAERS_ERA_POOL, fetch=_offline)
    got = st.code_status(f.person_drugs, st.ESTRIOL_CODE)
    assert len(got) == 1
    assert got.status.iloc[0] == pdr.MAPPING_FAILED
    assert got.cause.iloc[0] == ACH_CODE
    # Pooled with × 1/2 each, on strata that don't overlap.
    assert set(f.persons.cycle) == set(nhanes.FAERS_ERA_POOL)
    assert f.persons.groupby("stratum").cycle.nunique().max() == 1


# The dry-run on the no-encoder hand-off: needs every input of the main pool, the
# FAERS-era pool and MEPS. Its numbers aren't results.
HANDOFF = lookup.LOCAL_DIR / "crossfit" / "none_mlp_w-none_d64_lr0.01_do0_ep10"


@pytest.mark.skipif(
    not ((HANDOFF / "scores.npy").exists() and meps.RAW_DIR.exists() and all(p.exists() for p in FAERS_INPUTS)),
    reason="no-encoder hand-off, MEPS or NHANES files not in local/",
)
def test_dry_run_of_all_reruns_on_the_no_encoder_hand_off(tmp_path):
    out = st.run(HANDOFF, out_dir=tmp_path, fetch=_offline)
    o = out["overlap"].set_index("rerun")
    assert list(o.index) == list(st.RERUNS)
    assert (o.top == st.TOP).all()
    # Every main top-50 pair is either kept or has exactly one reason.
    assert (o.overlap + o[list(st.REASONS)].sum(axis=1) == st.TOP).all()
    d = out["dropouts"]
    assert (d.groupby("rerun").size() == st.TOP).all()
    assert ((d.reason == "") == d.in_top).all()
    c = out["faers_counts"].set_index("what").n
    assert c["new_codes"] == 249
    assert c["estriol_a70248_status"] == pdr.MAPPING_FAILED
    for name in out:
        assert (tmp_path / f"{name}.csv").exists()
