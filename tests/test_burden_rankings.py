import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from sief_druco.burden import lookup
from sief_druco.burden import prevalence as prv
from sief_druco.burden import rankings as rk
from sief_druco.burden import scores as sc
from sief_druco.burden.person_drugs import COVERED, NOT_IN_MODEL

# --- A tiny score table in the §2 format: 4 drugs, 6 pairs, 3 side effects ---

A, B, C, D = "CID000000001", "CID000000002", "CID000000003", "CID000000004"
DRUG_PAIRS = [(A, B), (A, C), (A, D), (B, C), (B, D), (C, D)]
SCORES = np.array(
    [
        [0.1, 0.2, 0.3],  # AB
        [0.5, 0.0, 0.5],  # AC
        [0.9, 0.9, 0.9],  # AD: nobody takes it
        [0.2, 0.2, 0.2],  # BC
        [0.4, 0.4, 0.4],  # BD: nobody takes it
        [0.0, 0.1, 0.0],  # CD
    ],
    dtype=np.float32,
)
LABELS = np.zeros(SCORES.shape, np.uint8)
LABELS[0, 0] = 1  # AB is labelled with side effect 0, at a predicted 0.1


def write_table(folder: Path, scores=SCORES, labels=LABELS, drug_pairs=DRUG_PAIRS) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    np.save(folder / "scores.npy", scores)
    labelled = labels.sum(axis=1) > 0
    pd.DataFrame({
        "pair_id": np.arange(len(drug_pairs)),
        "stitch_a": [a for a, _ in drug_pairs],
        "stitch_b": [b for _, b in drug_pairs],
        "labelled": labelled.astype(np.uint8),
        "regime": np.where(labelled, "S", "U"),
    }).to_csv(folder / "pairs.csv", index=False)
    sp.save_npz(folder / "labels.npz", sp.csr_matrix(labels))
    k = labels.shape[1]
    pd.DataFrame({
        "cui": [f"C{j:07d}" for j in range(k)], "name": [f"effect {j}" for j in range(k)],
        "labelled_pairs": labels.sum(axis=0),
        "auprc_U": 0.5, "predicted_U": 1.0, "observed_U": 1, "auprc_S": 0.6, "predicted_S": 1.0, "observed_S": 1,
    }).to_csv(folder / "side_effects.csv", index=False)
    (folder / "meta.json").write_text(json.dumps({"name": "tiny"}))
    return folder


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


def test_a_labelled_pair_uses_its_prediction_not_its_label(table):
    persons, pairs, pb = burdens(table, min_participants=1)
    t = rk.triple_burden(pb, table).set_index(["unit_a", "unit_b", "cui"])
    assert table.labels[[0], [0]][0] == 1  # AB × effect 0 is a label ...
    assert t.p[(A, B, "C0000000")] == pytest.approx(0.1)  # ... but its p is the prediction
    assert pb.set_index(["unit_a", "unit_b"]).loc[(A, B), "regime"] == "S"


# --- The reader ---


def test_reader_finds_rows_in_either_order(table):
    assert table.rows([B, D], [A, C]).tolist() == [0, 5]
    with pytest.raises(KeyError):
        table.rows([A], ["CID000000099"])


def test_reader_rejects_a_table_of_the_wrong_shape(tmp_path):
    write_table(tmp_path, scores=SCORES[:5])
    with pytest.raises(ValueError, match="scores.npy"):
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
