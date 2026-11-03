import numpy as np
import pytest
import scipy.sparse as sp

from sief_druco import decagon, graph, splits
from sief_druco.graph import DRUG, LABELLED_WITH, PPI, PROTEIN, TARGETED_BY, TARGETS
from sief_druco.splits import TRAIN

needs_raw = pytest.mark.skipif(
    not decagon.RAW_DIR.is_dir(), reason=f"raw Decagon files not in {decagon.RAW_DIR}"
)


# --- Drug input (no data needed) ---


def test_drug_input_keeps_columns_of_two_training_drugs_and_takes_the_mean():
    mono = sp.csr_matrix(
        np.array(
            [
                # col 0: 3 training drugs; col 1: 2; col 2: 1; col 3: only drug 3, not training
                [1, 1, 1, 0],
                [1, 1, 0, 0],
                [1, 0, 0, 0],
                [1, 0, 0, 1],
                [0, 0, 0, 0],  # no profile
            ],
            dtype=np.uint8,
        )
    )
    x = graph.drug_input(mono, np.array([0, 1, 2, 4]))
    assert x.dtype == np.float32
    expected_mean = np.array([[0.5, 0.5], [0.5, 0.5], [1, 0], [1, 0], [0, 0]])
    expected_n = np.array([3, 2, 1, 2, 0])  # all of a drug's side effects, kept or not
    np.testing.assert_allclose(x[:, :2], expected_mean)
    np.testing.assert_allclose(x[:, 2], np.log1p(expected_n), rtol=1e-6)


# --- The real dataset ---


@pytest.fixture(scope="module")
def data():
    return decagon.load()


@pytest.fixture(scope="module")
def folds(data):
    return {s: splits.draw(s, 0, data) for s in splits.SPLITS}


def edge_pair_ids(edge_index) -> np.ndarray:
    a, b = edge_index.numpy()
    return decagon.pair_id(a, b)


def drug_drug_types(g):
    return [t for t in g.edge_types if t[0] == t[2] == DRUG]


@needs_raw
def test_8152_mono_columns_when_all_drugs_are_training_drugs(data, folds):
    x = graph.drug_input(data.mono, np.arange(645))
    assert x.shape == (645, 8152 + 1)
    for split in ("random", "pair"):
        assert graph.MessageGraphs(data, folds[split]).x.shape == (645, 8152 + 1)


@needs_raw
def test_cold_drug_input_uses_its_training_drugs_only(data, folds):
    f = folds["cold_drug"]
    x = graph.MessageGraphs(data, f).x.numpy()
    np.testing.assert_array_equal(x, graph.drug_input(data.mono, f.drug_roles == TRAIN))
    assert x.shape[1] < 8152 + 1


@needs_raw
def test_drugs_without_a_profile_have_zero_input_and_stay_in_the_graph(data, folds):
    no_profile = np.flatnonzero(data.mono.getnnz(axis=1) == 0)
    assert len(no_profile) == 6
    g = graph.MessageGraphs(data, folds["pair"]).eval_graph
    assert g[DRUG].num_nodes == 645
    assert np.all(g[DRUG].x.numpy()[no_profile] == 0)
    assert np.all(np.isin(no_profile, g[LABELLED_WITH].edge_index.numpy()[0]))


@needs_raw
@pytest.mark.parametrize("split", ["pair", "cold_drug"])
@pytest.mark.parametrize("drug_edges", graph.DRUG_EDGES)
def test_drug_edges_are_exactly_the_labelled_message_pairs(data, folds, split, drug_edges):
    # So no validation or test pair is an edge, in either direction, in any message graph.
    f = folds[split]
    labelled = data.labels.getnnz(axis=1) > 0
    m = graph.MessageGraphs(data, f, drug_edges)
    rng = np.random.default_rng(0)
    epochs = [m.epoch(rng) for _ in range(3 if drug_edges == "collapsed" else 1)]
    cases = [(m.eval_graph, np.flatnonzero(f.roles == TRAIN))] + [
        (e.graph, e.message) for e in epochs
    ]
    for g, message in cases:
        ids = np.concatenate([edge_pair_ids(g[t].edge_index) for t in drug_drug_types(g)])
        assert np.all(f.roles[ids] == TRAIN)
        np.testing.assert_array_equal(np.unique(ids), message[labelled[message]])


@needs_raw
def test_per_side_effect_edges_follow_the_labels(data, folds):
    m = graph.MessageGraphs(data, folds["pair"], "per_side_effect")
    e = m.epoch(np.random.default_rng(0))
    labels = sp.csc_matrix(data.labels)
    for r in (0, 481, 962):
        positives = labels.indices[labels.indptr[r] : labels.indptr[r + 1]]
        ids = edge_pair_ids(e.graph[graph.side_effect_edge(r)].edge_index)
        np.testing.assert_array_equal(np.unique(ids), np.intersect1d(positives, e.message))


@needs_raw
@pytest.mark.parametrize("split", ["pair", "cold_drug"])
def test_message_and_supervision_split_the_training_pairs_70_30(data, folds, split):
    m = graph.MessageGraphs(data, folds[split])
    rng = np.random.default_rng(0)
    first, second = m.epoch(rng), m.epoch(rng)
    for e in (first, second):
        assert len(np.intersect1d(e.message, e.supervision)) == 0
        np.testing.assert_array_equal(np.union1d(e.message, e.supervision), m.train_pairs)
        assert len(e.message) / len(m.train_pairs) == pytest.approx(0.7, abs=1e-4)
    assert not np.array_equal(first.message, second.message)  # redrawn every epoch


@needs_raw
def test_random_split_edges_are_the_training_triples_every_epoch(data, folds):
    f = folds["random"]
    m = graph.MessageGraphs(data, f, "per_side_effect")
    e = m.epoch(np.random.default_rng(0))
    assert e.graph is m.eval_graph and e.message is None
    np.testing.assert_array_equal(e.supervision, np.arange(decagon.N_PAIRS))
    labels = sp.csc_matrix(data.labels)
    for r in (0, 481, 962):
        positives = labels.indices[labels.indptr[r] : labels.indptr[r + 1]]
        trained = positives[f.roles[positives, r] == TRAIN]
        ids = edge_pair_ids(e.graph[graph.side_effect_edge(r)].edge_index)
        np.testing.assert_array_equal(np.unique(ids), trained)


@needs_raw
def test_edge_types_and_reverses(data, folds):
    rgcn = graph.MessageGraphs(data, folds["pair"], "per_side_effect").eval_graph
    assert len(rgcn.edge_types) == 966

    g = graph.MessageGraphs(data, folds["pair"], "collapsed").eval_graph
    assert set(g.edge_types) == {PPI, TARGETS, TARGETED_BY, LABELLED_WITH}
    assert g[PROTEIN].num_nodes == 19_089
    assert g[PPI].num_edges == 2 * len(data.ppi)
    assert g[TARGETS].num_edges == 18_690
    assert (g[TARGETED_BY].edge_index == g[TARGETS].edge_index.flip(0)).all()
    for t in (PPI, LABELLED_WITH):  # same-type edges hold both directions
        src, dst = g[t].edge_index.numpy()
        n = g[t[0]].num_nodes
        np.testing.assert_array_equal(np.sort(src * n + dst), np.sort(dst * n + src))


@needs_raw
@pytest.mark.parametrize("drug_edges, n_types", [("per_side_effect", 963), ("collapsed", 1)])
def test_no_protein_graph_has_drug_drug_edges_only(data, folds, drug_edges, n_types):
    m = graph.MessageGraphs(data, folds["pair"], drug_edges, proteins=False)
    for g in (m.eval_graph, m.epoch(np.random.default_rng(0)).graph):
        assert g.node_types == [DRUG]
        assert len(g.edge_types) == len(drug_drug_types(g)) == n_types


def test_unknown_drug_edges_raises():
    with pytest.raises(ValueError):
        graph.MessageGraphs(None, None, "dense")
