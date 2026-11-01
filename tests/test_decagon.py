import hashlib

import numpy as np
import pytest

from sief_druco import decagon

needs_raw = pytest.mark.skipif(
    not decagon.RAW_DIR.is_dir(), reason=f"raw Decagon files not in {decagon.RAW_DIR}"
)


# --- Pair order (no data needed) ---


def test_pair_order_is_upper_triangle_row_by_row():
    pairs = decagon.all_pairs(4)
    assert pairs.tolist() == [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]


def test_pair_id_round_trips_for_all_pairs():
    pairs = decagon.all_pairs(decagon.N_DRUGS)
    assert len(pairs) == 207_690
    ids = decagon.pair_id(pairs[:, 0], pairs[:, 1])
    assert np.array_equal(ids, np.arange(len(pairs)))
    assert np.array_equal(decagon.pair_of(ids), pairs)


def test_pair_id_ignores_drug_order():
    a, b = np.array([3, 10]), np.array([10, 3])
    assert np.array_equal(decagon.pair_id(a, b), decagon.pair_id(b, a))


def test_pair_id_rejects_self_pairs():
    with pytest.raises(ValueError):
        decagon.pair_id(np.array([5]), np.array([5]))


@pytest.mark.parametrize("a, b", [(0, 645), (644, 645), (-1, 3)])
def test_pair_id_rejects_drugs_out_of_range(a, b):
    with pytest.raises(ValueError):
        decagon.pair_id(np.array([a]), np.array([b]))


@pytest.mark.parametrize("pid", [-1, 207_690])
def test_pair_of_rejects_ids_out_of_range(pid):
    with pytest.raises(ValueError):
        decagon.pair_of(np.array([pid]))


# --- Checksums (no data needed) ---


def test_checksum_mismatch_raises(tmp_path):
    (tmp_path / "a.csv").write_text("x\n")
    (tmp_path / "SHA256SUMS").write_text(f"{'0' * 64}  a.csv\n")
    with pytest.raises(ValueError, match="a.csv"):
        decagon.verify_checksums(tmp_path, ["a.csv"])


def test_checksum_match_passes(tmp_path):
    (tmp_path / "a.csv").write_text("x\n")
    digest = hashlib.sha256(b"x\n").hexdigest()
    (tmp_path / "SHA256SUMS").write_text(f"{digest}  a.csv\n")
    decagon.verify_checksums(tmp_path, ["a.csv"])


# --- The real dataset ---


@pytest.fixture(scope="module")
def cache_path(tmp_path_factory):
    return tmp_path_factory.mktemp("cache") / "decagon.npz"


@pytest.fixture(scope="module")
def data(cache_path):
    # Parse from the raw CSVs into a fresh cache, so the parser itself is tested.
    return decagon.load(cache_path=cache_path)


@needs_raw
def test_counts(data):
    assert len(data.drug_ids) == 645
    assert len(data.side_effect_ids) == 963
    assert len(data.protein_ids) == 19_089
    assert data.labels.shape == (207_690, 963)
    assert data.labels.nnz == 4_576_287
    assert (data.labels.getnnz(axis=1) > 0).sum() == 63_472


@needs_raw
def test_targets(data):
    assert len(data.targets) == 18_690
    with_targets = np.unique(data.targets[:, 0])
    assert len(with_targets) == 284
    assert decagon.N_DRUGS - len(with_targets) == 361


@needs_raw
def test_mono_profiles(data):
    has_mono = data.mono.getnnz(axis=1) > 0
    assert data.mono.shape[0] == 645
    assert has_mono.sum() == 639
    has_targets = np.isin(np.arange(645), data.targets[:, 0])
    assert (~has_mono & ~has_targets).sum() == 5


@needs_raw
def test_mono_and_combo_cuis_are_disjoint(data):
    assert not set(data.mono_ids) & set(data.side_effect_ids)


@needs_raw
def test_ids_are_sorted_and_unique(data):
    for ids in (data.drug_ids, data.protein_ids, data.side_effect_ids, data.mono_ids):
        assert np.all(ids[:-1] < ids[1:])


@needs_raw
def test_labelled_pairs_use_canonical_order(data):
    assert np.array_equal(data.pairs, decagon.all_pairs(645))
    assert np.all(data.pairs[:, 0] < data.pairs[:, 1])
    labelled = np.flatnonzero(data.labels.getnnz(axis=1))
    a, b = decagon.pair_of(labelled).T
    assert np.array_equal(decagon.pair_id(a, b), labelled)


@needs_raw
def test_ppi_edges_stored_once(data):
    assert len(data.ppi) == 715_612
    assert np.all(data.ppi[:, 0] != data.ppi[:, 1])
    undirected = np.sort(data.ppi, axis=1)
    assert len(np.unique(undirected, axis=0)) == len(data.ppi)


@needs_raw
def test_known_pair_and_side_effect(data):
    # A row of bio-decagon-combo.csv: CID000003929, CID000150610, C0002871, anaemia.
    a, b = data.drug_index(["CID000003929", "CID000150610"])
    r = data.side_effect_index(["C0002871"])[0]
    assert data.labels[decagon.pair_id(a, b), r] == 1
    assert data.side_effect_names[r] == "anaemia"
    # Hypermagnesemia (387 pairs) is below the 500-pair cut.
    with pytest.raises(KeyError):
        data.side_effect_index(["C0151714"])


@needs_raw
def test_lookup_rejects_ids_longer_than_a_real_one(data):
    # Each is a real ID with one character appended.
    with pytest.raises(KeyError):
        data.side_effect_index(["C00028711"])
    with pytest.raises(KeyError):
        data.drug_index(["CID0000039299"])


@needs_raw
def test_cache_round_trips(data, cache_path):
    cached = decagon.load(cache_path=cache_path)
    for field in decagon.Decagon.__dataclass_fields__:
        x, y = getattr(data, field), getattr(cached, field)
        if hasattr(x, "toarray"):
            assert (x != y).nnz == 0
        else:
            assert np.array_equal(x, y)
