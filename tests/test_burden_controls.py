from itertools import combinations

import numpy as np
import pandas as pd
import pytest

from sief_druco.burden import controls as ct
from sief_druco.burden import lookup
from sief_druco.burden import prevalence as prv
from sief_druco.burden import rankings as rk
from sief_druco.burden import scores as sc
from sief_druco.burden.person_drugs import COVERED
from tiny_score_table import write_table

# --- A synthetic score table: 5 drugs, 10 pairs, 2 side effects ---

D1, D2, D3, D4, D5 = (f"CID00000000{i}" for i in range(1, 6))
DRUG_PAIRS = list(combinations([D1, D2, D3, D4, D5], 2))  # row i = DRUG_PAIRS[i]
ROW = {p: i for i, p in enumerate(DRUG_PAIRS)}
SCORES = np.full((10, 2), 0.5, np.float32)
SCORES[:, 0] = np.arange(10) / 10  # effect 0: (D4, D5) highest, (D3, D5) second
SCORES[[ROW[D1, D2], ROW[D2, D3]], 1] = 0.95  # effect 1: two pairs tie at the top
LABELS = np.zeros(SCORES.shape, np.uint8)
LABELS[ROW[D4, D5], 1] = 1  # (D4, D5) is labelled, but not with its entry's effect 0
E0, E1 = "C0000000", "C0000001"

ROWS = [
    # entry, kind, name, id
    (1, "a", "drug 4", D4), (1, "b", "drug 5", D5), (1, "side_effect", "effect 0", E0),
    (2, "a", "drug 3", D3), (2, "b", "drug 5", D5), (2, "side_effect", "effect 0", E0),
    (3, "member", "drug 1", D1), (3, "member", "drug 2", D2), (3, "member", "drug 3", D3),
    (3, "side_effect", "effect 1", E1),
    (4, "a", "drug 1", D1), (4, "b", "drug 2", D2), (4, "b", "outside", ""),
    (4, "side_effect", "unmodelled effect", "C9999999"),
    (5, "out_of_scope", "", ""),
]


@pytest.fixture
def table(tmp_path):
    return sc.read(write_table(tmp_path / "handoff", SCORES, LABELS, DRUG_PAIRS))


@pytest.fixture
def controls(tmp_path):
    df = pd.DataFrame(ROWS, columns=["entry", "kind", "name", "id"])
    df = df.assign(label=[f"entry {e}" for e in df.entry], **{"class": "", "note": ""})
    path = tmp_path / "controls.csv"
    df[["entry", "label", "kind", "class", "name", "id", "note"]].to_csv(path, index=False)
    return ct.read(path)


def by_pair(df):
    return {(r.entry, r.unit_a, r.unit_b): r for r in df.itertuples()}


# --- Model view ---


def test_control_pairs_are_a_times_b_or_two_distinct_members(controls, table):
    got = set(map(tuple, ct.pairs(controls, table)[["entry", "unit_a", "unit_b"]].to_numpy()))
    assert got == {
        (1, D4, D5), (2, D3, D5),
        (3, D1, D2), (3, D1, D3), (3, D2, D3),
        (4, D1, D2),  # "outside" has no id, so it forms no pair
    }


def test_percentile_and_recovered_match_hand_computed_values(controls, table):
    _, report = ct.model_view(controls, table)
    r = by_pair(report)
    # Highest of 10 on effect 0: 9 pairs below → 90; 1 of 10 scores as high = exactly 10%.
    assert r[1, D4, D5].best_percentile == pytest.approx(90.0) and r[1, D4, D5].recovered
    # Second highest: 80; 2 of 10 score as high → 20%, not in the top 10%.
    assert r[2, D3, D5].best_percentile == pytest.approx(80.0) and not r[2, D3, D5].recovered
    # A tie at the top counts against the control: 2 of 10 score ≥ 0.95.
    for pair in [(D1, D2), (D2, D3)]:
        assert r[(3, *pair)].best_percentile == pytest.approx(80.0) and not r[(3, *pair)].recovered
    # 0.5 on effect 1: the two 0.95s are above, 8 pairs are at 0.5, none below.
    assert r[3, D1, D3].best_percentile == pytest.approx(0.0)


def test_labelled_controls_are_marked(controls, table):
    triples, report = ct.model_view(controls, table)
    r = by_pair(report)
    assert r[1, D4, D5].labelled and not r[1, D4, D5].labelled_mapped  # labelled, but not with effect 0
    assert not r[2, D3, D5].labelled
    assert not triples.labelled.any()


def test_a_control_with_no_matching_side_effect_is_not_checkable(controls, table):
    _, report = ct.model_view(controls, table)
    c = by_pair(report)[4, D1, D2]
    assert not c.checkable and pd.isna(c.recovered) and pd.isna(c.best_percentile)

    summary = ct.entries(controls, ct.burden_view(report, pb_frame(table, {})))
    row = summary.set_index("entry").loc[4]
    assert (row.controls, row.checkable, row.recovered) == (1, 0, 0)
    overall = summary.set_index("entry").loc["all"]
    assert overall.recovered_text == "1 of 5 recovered"  # entry 4's control is not in n

    un = ct.unchecked(controls, table, report)
    assert set(map(tuple, un[["entry", "reason"]].to_numpy())) == {
        (4, ct.NOT_IN_645), (4, ct.NO_SIDE_EFFECT), (5, ct.OUT_OF_SCOPE),
    }
    assert un[un.reason == ct.NOT_IN_645].name.tolist() == ["outside"]


# --- Burden view ---


def pb_frame(table, participants: dict) -> pd.DataFrame:
    """`rankings.pair_burden` with the default minimum of 30, from a prevalence table of
    the given participants per pair (prevalence = participants, so burden follows it)."""
    prev = pd.DataFrame(
        [(a, b, COVERED, n, float(n), float(n), float(n)) for (a, b), n in participants.items()],
        columns=prv.COLUMNS,
    )
    return rk.pair_burden(prev, table)


def test_burden_view_ranks_only_pairs_with_30_participants(controls, table):
    _, report = ct.model_view(controls, table)
    # Burden = participants × Σp: (D4, D5) 29 × 1.4 = 40.6 is the highest, but has only
    # 29 participants; (D3, D5) 30 × 1.3 = 39; (D1, D2) 31 × 0.95 = 29.45.
    pb = pb_frame(table, {(D4, D5): 29, (D3, D5): 30, (D1, D2): 31})
    r = by_pair(ct.burden_view(report, pb))
    assert r[1, D4, D5].participants == 29 and pd.isna(r[1, D4, D5].rank) and not r[1, D4, D5].top50
    assert r[2, D3, D5].rank == 1 and r[2, D3, D5].top50
    assert r[3, D1, D2].rank == 2 and r[3, D1, D2].top50
    assert r[3, D1, D3].participants == 0 and pd.isna(r[3, D1, D3].rank)  # nobody takes it


def test_burden_view_top_k_boundary(controls, table):
    _, report = ct.model_view(controls, table)
    pb = pb_frame(table, {(D3, D5): 30, (D1, D2): 31})  # ranks 1 and 2, as above
    r = by_pair(ct.burden_view(report, pb, top=1))
    assert r[2, D3, D5].top50 and not r[3, D1, D2].top50


def test_model_view_does_not_depend_on_prevalence(controls, table):
    _, report = ct.model_view(controls, table)
    model_cols = ["best_side_effect", "best_p", "best_percentile", "recovered", "checkable"]
    views = [
        ct.burden_view(report, pb_frame(table, p))
        for p in [{}, {(D4, D5): 1000, (D3, D5): 30}, {(a, b): 30 + i for i, (a, b) in enumerate(DRUG_PAIRS)}]
    ]
    for v in views[1:]:
        pd.testing.assert_frame_equal(v[model_cols], views[0][model_cols])
    # Percentiles are among all 10 pairs, including those nobody takes.
    assert by_pair(views[0])[1, D4, D5].best_percentile == pytest.approx(90.0)


# --- The committed list on the no-encoder hand-off (dry run: `python -m sief_druco.burden.controls`) ---

HANDOFF = lookup.LOCAL_DIR / "crossfit" / "none_mlp_w-none_d64_lr0.01_do0_ep10"


def test_committed_list_reads():
    c = ct.read()
    assert set(c.kind) == {"a", "b", "member", "side_effect", "note", "out_of_scope"}
    assert c.entry.nunique() == 12


@pytest.mark.skipif(not (HANDOFF / "scores.npy").exists(), reason="no-encoder hand-off not in local/")
def test_committed_list_on_the_real_hand_off():
    t = sc.read(HANDOFF)
    controls = ct.read()
    triples, report = ct.model_view(controls, t)
    assert report.checkable.all()  # every entry maps to a side effect of the 963
    assert report.best_percentile.between(0, 100).all()
    assert set(report.entry) == set(range(1, 13)) - {5}
    un = ct.unchecked(controls, t, report)
    assert set(un.reason) == {ct.OUT_OF_SCOPE, ct.NOT_IN_645}
