import subprocess
import sys

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from sief_druco.burden import lookup
from sief_druco.burden import prevalence as prv
from sief_druco.burden import rankings as rk
from sief_druco.burden import scores as sc
from sief_druco.burden.person_drugs import COVERED, NOT_IN_MODEL
from tiny_score_table import A, B, C, D, LABELS, SCORES, write_table

# --- The tiny score table of `tiny_score_table` ---

# 8 persons, weights 1..8 (Σ 36), two strata of two PSUs.
PERSON_UNITS = {0: [A, B], 1: [A, B], 2: [A, B, C], 3: [C, D], 4: [A, "in:999"]}
WEIGHTS = np.arange(1, 9, dtype=float)


def survey():
    n = len(WEIGHTS)
    persons = pd.DataFrame({
        "source": "NHANES", "cycle": "J", "seqn": np.arange(n), "weight": WEIGHTS,
        "stratum": np.arange(n) % 2, "psu": (np.arange(n) // 2) % 2,
    })
    rows = []
    for s, units in PERSON_UNITS.items():
        for i, a in enumerate(units):
            for b in units[i + 1:]:
                status = COVERED if a.startswith("CID") and b.startswith("CID") else NOT_IN_MODEL
                rows.append(dict(source="NHANES", cycle="J", seqn=s, unit_a=min(a, b), unit_b=max(a, b), status=status))
    return persons, pd.DataFrame(rows)


@pytest.fixture
def table(tmp_path):
    return sc.read(write_table(tmp_path / "handoff"))


def burdens(table, min_participants=2):
    persons, pairs = survey()
    pb = rk.pair_burden(prv.pair_prevalence(persons, pairs), table, min_participants)
    return persons, pairs, pb


# Hand-computed prevalence per 100,000: AB 6/36, AC 3/36, BC 3/36, CD 4/36.
PREV = {(A, B): 6 / 36, (A, C): 3 / 36, (B, C): 3 / 36, (C, D): 4 / 36}
PREV = {k: 1e5 * v for k, v in PREV.items()}


# --- The three formulas ---


def test_pair_burden_is_prevalence_times_summed_scores(table):
    _, _, pb = burdens(table)
    got = dict(zip(zip(pb.unit_a, pb.unit_b), pb.burden))
    assert got == pytest.approx({
        (A, B): PREV[A, B] * 0.6, (A, C): PREV[A, C] * 1.0, (B, C): PREV[B, C] * 0.6, (C, D): PREV[C, D] * 0.1,
    })
    # The interval is the prevalence interval with p fixed.
    assert pb.burden_low.to_numpy() == pytest.approx(pb.ci_low * pb.sum_p)
    assert pb.burden_high.to_numpy() == pytest.approx(pb.ci_high * pb.sum_p)


def test_only_pairs_with_enough_participants_are_ranked(table):
    _, _, pb = burdens(table, min_participants=2)
    # AB (3 participants) is the only pair with ≥ 2; the not-in-model pair is left out.
    assert len(pb) == 4
    ranked = pb[pb["rank"].notna()]
    assert list(zip(ranked.unit_a, ranked.unit_b, ranked["rank"])) == [(A, B, 1)]


def test_triple_burden_is_prevalence_times_score(table):
    _, _, pb = burdens(table, min_participants=1)
    t = rk.triple_burden(pb, table).set_index(["unit_a", "unit_b", "cui"])
    assert t.burden[(A, B, "C0000001")] == pytest.approx(PREV[A, B] * 0.2)
    assert t.burden[(A, C, "C0000002")] == pytest.approx(PREV[A, C] * 0.5)
    assert t.burden[(C, D, "C0000000")] == 0
    assert t["rank"].tolist() == list(range(1, len(t) + 1))
    assert t.burden.is_monotonic_decreasing


def test_side_effect_burden_sums_over_all_covered_pairs(table):
    persons, pairs, _ = burdens(table)
    se = rk.side_effect_burden(persons, pairs, table).set_index("cui")
    hand = {
        "C0000000": PREV[A, B] * 0.1 + PREV[A, C] * 0.5 + PREV[B, C] * 0.2 + PREV[C, D] * 0.0,
        "C0000001": PREV[A, B] * 0.2 + PREV[A, C] * 0.0 + PREV[B, C] * 0.2 + PREV[C, D] * 0.1,
        "C0000002": PREV[A, B] * 0.3 + PREV[A, C] * 0.5 + PREV[B, C] * 0.2 + PREV[C, D] * 0.0,
    }
    assert se.burden.to_dict() == pytest.approx(hand)
    assert (se.ci_low < se.burden).all() and (se.burden < se.ci_high).all()
    assert se.loc["C0000002", "rank"] == 1
    assert {"auprc_U", "auprc_S", "predicted_U", "observed_U"} <= set(se.columns)


# --- Sums ---


def test_triples_sum_to_their_pair(table):
    _, _, pb = burdens(table, min_participants=1)
    t = rk.triple_burden(pb, table)
    by_pair = t.groupby(["unit_a", "unit_b"]).burden.sum()
    assert by_pair.to_dict() == pytest.approx(pb.set_index(["unit_a", "unit_b"]).burden.to_dict())


def test_side_effects_sum_to_the_total_over_covered_pairs(table):
    persons, pairs, pb = burdens(table, min_participants=10)  # nothing ranked: totals still count every pair
    se = rk.side_effect_burden(persons, pairs, table)
    assert pb["rank"].isna().all()
    assert se.burden.sum() == pytest.approx(pb.burden.sum())


def test_a_labelled_pair_uses_its_prediction_not_its_label(table, tmp_path):
    assert table.labels[[0], [0]][0] == 1  # AB × effect 0 is a label, predicted at 0.1
    unlabelled = sc.read(write_table(tmp_path / "unlabelled", labels=np.zeros_like(LABELS)))
    persons, pairs = survey()

    def all_burdens(t):
        pb = rk.pair_burden(prv.pair_prevalence(persons, pairs), t, min_participants=1)
        return [
            pb.burden.to_numpy(),
            rk.triple_burden(pb, t).burden.to_numpy(),
            rk.side_effect_burden(persons, pairs, t).burden.to_numpy(),
        ]

    # Every burden is the same with or without the label: only the prediction counts.
    for with_label, without in zip(all_burdens(table), all_burdens(unlabelled)):
        assert with_label == pytest.approx(without)
    pb = rk.pair_burden(prv.pair_prevalence(persons, pairs), table, min_participants=1).set_index(["unit_a", "unit_b"])
    assert pb.loc[(A, B), "regime"] == "S" and pb.loc[(A, B), "sum_p"] == pytest.approx(0.6)


# --- The reader ---


def test_reader_finds_rows_in_either_order(table):
    assert table.rows([B, D], [A, C]).tolist() == [0, 5]
    with pytest.raises(KeyError):
        table.rows([A], ["CID000000099"])


def test_reader_rejects_a_table_of_the_wrong_shape(tmp_path):
    write_table(tmp_path, scores=SCORES[:5])
    with pytest.raises(ValueError, match="scores.npy"):
        sc.read(tmp_path)


def test_reader_rejects_labels_of_the_wrong_shape(tmp_path):
    write_table(tmp_path)
    sp.save_npz(tmp_path / "labels.npz", sp.csr_matrix(np.zeros((6, 4), np.uint8)))
    with pytest.raises(ValueError, match="labels.npz"):
        sc.read(tmp_path)


def test_reader_rejects_side_effects_that_do_not_match_the_columns(tmp_path):
    write_table(tmp_path)
    se = pd.read_csv(tmp_path / "side_effects.csv")
    se.iloc[:2].to_csv(tmp_path / "side_effects.csv", index=False)
    with pytest.raises(ValueError, match="side_effects.csv"):
        sc.read(tmp_path)


def test_reader_rejects_pairs_out_of_order(tmp_path):
    write_table(tmp_path)
    pairs = pd.read_csv(tmp_path / "pairs.csv")
    pairs.iloc[[1, 0, 2, 3, 4, 5]].to_csv(tmp_path / "pairs.csv", index=False)
    with pytest.raises(ValueError, match="out of order"):
        sc.read(tmp_path)


def test_reader_needs_every_file(tmp_path):
    write_table(tmp_path)
    (tmp_path / "labels.npz").unlink()
    with pytest.raises(FileNotFoundError):
        sc.read(tmp_path)


# --- No torch ---


def test_importing_the_burden_package_does_not_load_torch():
    code = (
        "import pkgutil, importlib, sys, sief_druco.burden as b\n"
        "for m in pkgutil.iter_modules(b.__path__): importlib.import_module(f'sief_druco.burden.{m.name}')\n"
        "assert 'torch' not in sys.modules, 'torch was imported'\n"
    )
    subprocess.run([sys.executable, "-c", code], check=True)


# --- The real format: the no-encoder hand-off (the dry run is `python -m sief_druco.burden.rq3`) ---

HANDOFF = lookup.LOCAL_DIR / "crossfit" / "none_mlp_w-none_d64_lr0.01_do0_ep10"


@pytest.mark.skipif(not (HANDOFF / "scores.npy").exists(), reason="no-encoder hand-off not in local/")
def test_reader_opens_the_real_hand_off():
    t = sc.read(HANDOFF)
    assert t.scores.shape == (207690, 963) and isinstance(t.scores, np.memmap)
    assert t.pairs.regime.isin(["U", "S"]).all()
