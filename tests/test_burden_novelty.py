import numpy as np
import pandas as pd
import pytest

from sief_druco import decagon
from sief_druco.burden import decagon_table as dt
from sief_druco.burden import handfix, lookup, nhanes, rq2
from sief_druco.burden import novelty as nv
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden import prevalence as prv
from sief_druco.burden import rankings as rk
from sief_druco.burden import scores as sc
from sief_druco.burden.ingredients import NAME, NAME_NOT_FOUND, Resolution
from sief_druco.burden.novelty import DOCUMENTED, NOT_CHECKABLE, NOT_DOCUMENTED
from tiny_score_table import A, B, C, D, write_table

# --- DDInter on the four drugs of `tiny_score_table` ---

# IN rxcui -> STITCH ID; "15" is an alias (a salt) of A, "99" is outside the 645.
STITCH_OF = {"1": A, "15": A, "2": B, "3": C, "4": D}
NAMES = {
    "Drug A": ["1"], "Drug A salt": ["15"], "Drug B": ["2"], "Drug C": ["3"], "Other": ["99"],
    "Drug D (topical)": [],  # D's only DDInter name reaches no IN
}


def resolve(name):
    rxcuis = NAMES[name]
    if not rxcuis:
        return Resolution((), NAME_NOT_FOUND, "")
    return Resolution(tuple((r, name) for r in rxcuis), None, NAME)


ROWS = pd.DataFrame([
    ("Drug A", "Drug B", "Minor"),  # AB at two levels, in either order: Major is kept
    ("Drug B", "Drug A", "Major"),
    ("Drug A salt", "Drug C", "Unknown"),  # AC via A's alias, only "Unknown"
    ("Drug C", "Other", "Moderate"),
    ("Drug D (topical)", "Drug A", "Major"),  # D never reaches an IN
], columns=["drug_a", "drug_b", "level"])


@pytest.fixture
def clinical():
    return nv.map_ddinter(ROWS, resolve, STITCH_OF)


def label(clinical, a, b):
    return clinical.label([a], [b]).iloc[0].to_dict()


def test_a_pair_listed_at_two_levels_keeps_the_highest(clinical):
    assert label(clinical, A, B) == {"clinically_documented": DOCUMENTED, "ddinter_level": "Major"}
    assert label(clinical, B, A) == label(clinical, A, B)


def test_unknown_is_not_documented(clinical):
    # Reached through A's alias: the alias joins its Decagon drug (§4).
    assert label(clinical, A, C) == {"clinically_documented": NOT_DOCUMENTED, "ddinter_level": "Unknown"}


def test_a_pair_of_two_listed_drugs_that_is_not_listed_is_not_documented(clinical):
    assert label(clinical, B, C) == {"clinically_documented": NOT_DOCUMENTED, "ddinter_level": ""}


def test_a_name_that_reaches_no_in_makes_the_pair_not_checkable(clinical):
    # DDInter lists D + A as Major, but "Drug D (topical)" reaches no IN.
    assert label(clinical, A, D)["clinically_documented"] == NOT_CHECKABLE
    assert label(clinical, C, D)["clinically_documented"] == NOT_CHECKABLE
    assert clinical.drugs == {A, B, C}
    assert clinical.names.set_index("name").cause["Drug D (topical)"] == NAME_NOT_FOUND


def test_load_ddinter_reads_all_14_files(tmp_path):
    for i, f in enumerate(nv.ddinter_files(tmp_path)):
        pd.DataFrame({
            "DDInterID_A": ["DDInter1"], "Drug_A": [f"Drug {i}"], "DDInterID_B": ["DDInter2"],
            "Drug_B": ["Drug X"], "Level": [nv.LEVELS[i % 4]],
        }).to_csv(f, index=False)
    rows = nv.load_ddinter(tmp_path)
    assert len(rows) == 14 and list(rows.columns) == ["drug_a", "drug_b", "level"]
    nv.ddinter_files(tmp_path)[3].unlink()
    with pytest.raises(FileNotFoundError, match="code_D"):
        nv.load_ddinter(tmp_path)


# --- Labels and tables on the tiny score table ---

# Covered pairs with hand-set prevalence (per 100,000) and participants.
PAIRS = pd.DataFrame({
    "unit_a": [A, A, B, C], "unit_b": [B, C, C, D],
    "participants": [40, 35, 12, 5], "prevalence": [600.0, 300.0, 300.0, 400.0],
    "ci_low": [500.0, 250.0, 200.0, 300.0], "ci_high": [700.0, 350.0, 400.0, 500.0],
})


@pytest.fixture
def table(tmp_path):
    return sc.read(write_table(tmp_path / "handoff"))


@pytest.fixture
def labelled(table, clinical):
    pb = PAIRS.assign(status=pdr.COVERED)
    pb = rk.pair_burden(pb, table, min_participants=30)
    return nv.label_pairs(pb, table, clinical)


def test_known_to_decagon_is_any_label_for_pairs_and_the_label_for_triples(table, clinical, labelled):
    # `tiny_score_table`: AB is labelled, with side effect 0 only.
    pairs = labelled.set_index(["unit_a", "unit_b"])
    assert pairs.known_to_decagon.to_dict() == {
        (A, B): DOCUMENTED, (A, C): NOT_DOCUMENTED, (B, C): NOT_DOCUMENTED, (C, D): NOT_DOCUMENTED,
    }
    t = nv.label_triples(rk.triple_burden(labelled, table), table, clinical).set_index(["unit_a", "unit_b", "cui"])
    assert t.known_to_decagon[(A, B, "C0000000")] == DOCUMENTED
    assert t.known_to_decagon[(A, B, "C0000001")] == NOT_DOCUMENTED
    assert t.known_to_decagon[(A, C, "C0000000")] == NOT_DOCUMENTED
    # Triples inherit the pair's clinical status, marked as such.
    assert (t.clinically_documented.xs((A, B), level=[0, 1]) == DOCUMENTED).all()
    assert (t.clinical_scope == nv.PAIR_PART).all()


def test_candidates_are_neither_known_nor_documented(labelled):
    got = labelled.set_index(["unit_a", "unit_b"]).candidate.to_dict()
    # AC (Unknown) and BC (not listed) are candidates; CD is not checkable.
    assert got == {(A, B): False, (A, C): True, (B, C): True, (C, D): False}


def test_shares_of_burden_add_to_100(table, labelled):
    s = nv.burden_shares(labelled, table)
    assert s.groupby("level").burden_pct.sum().to_dict() == pytest.approx({"pairs": 100, "triples": 100})
    assert len(s) == 2 * 2 * 3  # every cell, empty ones included
    # Pairs and triples hold the same total burden.
    assert s.groupby("level").burden.sum().to_dict() == pytest.approx({"pairs": labelled.burden.sum()} | {"triples": labelled.burden.sum()})
    cell = s.set_index(["level", "known_to_decagon", "clinically_documented"]).burden
    # Burden: AB 600 × 0.6, AC 300 × 1.0, BC 300 × 0.6, CD 400 × 0.1.
    assert cell[("pairs", DOCUMENTED, DOCUMENTED)] == pytest.approx(360)
    assert cell[("pairs", NOT_DOCUMENTED, NOT_DOCUMENTED)] == pytest.approx(300 + 180)
    assert cell[("pairs", NOT_DOCUMENTED, NOT_CHECKABLE)] == pytest.approx(40)
    # Triples: only AB × effect 0 (p = 0.1) is known to Decagon.
    assert cell[("triples", DOCUMENTED, DOCUMENTED)] == pytest.approx(60)
    assert cell[("triples", NOT_DOCUMENTED, DOCUMENTED)] == pytest.approx(300)


def test_top_k_counts_against_the_base_rate(table, clinical, labelled):
    triples = nv.label_triples(rk.triple_burden(labelled, table), table, clinical)
    t = nv.top_k(labelled, triples, table, k=4).set_index(["level", "known_to_decagon", "clinically_documented"])
    assert t.groupby("level").top_pct.sum().to_dict() == pytest.approx({"pairs": 100, "triples": 100})
    assert t.groupby("level").base_pct.sum().to_dict() == pytest.approx({"pairs": 100, "triples": 100})
    # Ranked (≥ 30 participants): AB, AC. Base: all 4 covered pairs, and their 12 triples.
    assert t.top_n.xs("pairs").sum() == 2 and t.base_n.xs("pairs").sum() == 4
    assert t.top_n[("pairs", DOCUMENTED, DOCUMENTED)] == 1
    assert t.top_n[("pairs", NOT_DOCUMENTED, NOT_DOCUMENTED)] == 1
    assert t.top_n.xs("triples").sum() == 4 and t.base_n.xs("triples").sum() == 12
    assert t.base_n[("triples", DOCUMENTED, DOCUMENTED)] == 1
    assert t.base_n[("triples", NOT_DOCUMENTED, DOCUMENTED)] == 2
    assert t.base_n[("triples", NOT_DOCUMENTED, NOT_CHECKABLE)] == 3


def test_candidate_list_has_its_own_floor_and_flags_the_imprecise(labelled):
    c = nv.candidates(labelled)
    # AC (35 participants) and BC (12); CD (5) would be under the floor anyway.
    assert list(zip(c.unit_a, c.unit_b)) == [(A, C), (B, C)]
    assert c.imprecise.tolist() == [False, True]
    assert c.burden.is_monotonic_decreasing


def test_notes_give_the_unknown_share(labelled, clinical):
    n = nv.notes(labelled, clinical).set_index("what").value
    total = labelled.burden.sum()
    assert n["ddinter_unknown_burden_pct"] == pytest.approx(100 * 300 / total)
    assert n["checkable_weight_pct"] == pytest.approx(100 * 1200 / 1600)
    assert n["ddinter_names_reaching_an_in"] == 5 and n["ddinter_names_reaching_decagon"] == 4


# --- NHANES 2017–2018 and the real DDInter (raw files, lookup caches and the final hand-off needed) ---

HANDOFF = lookup.LOCAL_DIR / "crossfit" / "additive_free_w-none_lr0.01"
LOCAL_INPUTS = [
    nhanes.RAW_DIR,
    decagon.RAW_DIR,
    nv.DDINTER_DIR,
    handfix.SIDER_NAMES,
    HANDOFF / "scores.npy",
    lookup.CACHE_DIR / pdr.NHANES_RXNAV_CACHE,
    lookup.CACHE_DIR / dt.PUBCHEM_CACHE,
    lookup.CACHE_DIR / dt.RXNAV_CACHE,
    lookup.CACHE_DIR / nv.DDINTER_RXNAV_CACHE,
]

needs_local = pytest.mark.skipif(
    not all(p.exists() for p in LOCAL_INPUTS),
    reason="raw NHANES / Decagon / DDInter / SIDER files, the lookup caches or the final hand-off not in local/",
)


def _offline(url, form):
    raise OSError(f"tests make no network calls: {url}")


def _2017_18_shares(overrides_path):
    e = rq2.exposure(["J"], overrides_path=overrides_path, fetch=_offline)
    clinical = nv.ddinter(e.table.stitch_of, fetch=_offline)
    prev = prv.pair_prevalence(e.persons, e.pairs)
    pairs = nv.label_pairs(prev[prev.status == pdr.COVERED], table=sc.read(HANDOFF), clinical=clinical)
    share = pairs.groupby(["known_to_decagon", "ddinter_level"]).prevalence.sum() / pairs.prevalence.sum() * 100
    by_state = nv.weight_shares(pairs).set_index(["known_to_decagon", "clinically_documented"]).pct
    return {
        "checkable": 100 - by_state.xs(NOT_CHECKABLE, level=1).sum(),
        "graded_in_decagon": share[DOCUMENTED].reindex(list(nv.GRADED)).sum(),
        "graded_not_in_decagon": share[NOT_DOCUMENTED].reindex(list(nv.GRADED)).sum(),
        "unknown": share.xs(nv.UNKNOWN, level=1).sum(),
        "novelty_pool": by_state[(NOT_DOCUMENTED, NOT_DOCUMENTED)],
    }


@needs_local
def test_2017_18_covered_pair_weight_matches_the_pilot():
    # Spec §10, from #19's pilot, which had no hand fixes: the same run without them.
    got = _2017_18_shares(overrides_path=None)
    pilot = {"checkable": 99.5, "graded_in_decagon": 31.9, "graded_not_in_decagon": 1.0, "unknown": 60.3, "novelty_pool": 3.4}
    assert got == pytest.approx(pilot, abs=0.2)


@needs_local
def test_2017_18_shares_stay_close_after_the_hand_fixes():
    # The hand fixes add covered pairs (8,856 → 9,609); the shares move by up to 2 points.
    got = _2017_18_shares(overrides_path=dt.OVERRIDES_PATH)
    pilot = {"checkable": 99.5, "graded_in_decagon": 31.9, "graded_not_in_decagon": 1.0, "unknown": 60.3, "novelty_pool": 3.4}
    assert got == pytest.approx(pilot, abs=2.0)
