from types import SimpleNamespace

import pytest
import torch
from torch_geometric.data import HeteroData
from torch_geometric.nn import RGCNConv

from sief_druco.encoders import ENCODERS, N_BASES, RGCN, dense_adjacency
from sief_druco.graph import DRUG, LABELLED_WITH, PPI, PROTEIN, TARGETED_BY, TARGETS, side_effect_edge

N_DRUGS, N_PROTEINS, N_SIDE_EFFECTS, N_IN, D = 14, 20, 6, 9, 8
GNNS = ("rgcn", "sage", "gat")


def random_pairs(n: int, k: int, gen) -> torch.Tensor:
    """`k` distinct unordered pairs of `n` nodes, as an edge index with both directions."""
    a, b = torch.triu_indices(n, n, 1)[:, torch.randperm(n * (n - 1) // 2, generator=gen)[:k]]
    return torch.stack([torch.cat([a, b]), torch.cat([b, a])])


def small_graphs(drug_edges: str, proteins: bool = True, seed: int = 0):
    """A stand-in for `MessageGraphs`: a small random message graph of either kind."""
    gen = torch.Generator().manual_seed(seed)
    g = HeteroData()
    g[DRUG].x = torch.rand(N_DRUGS, N_IN, generator=gen)
    if proteins:
        g[PROTEIN].num_nodes = N_PROTEINS
        g[PPI].edge_index = random_pairs(N_PROTEINS, 40, gen)
        targets = torch.stack([torch.randint(0, N_DRUGS - 3, (25,), generator=gen),  # 3 drugs: none
                               torch.randint(0, N_PROTEINS, (25,), generator=gen)]).unique(dim=1)
        g[TARGETS].edge_index = targets
        g[TARGETED_BY].edge_index = targets.flip(0)
    if drug_edges == "collapsed":
        g[LABELLED_WITH].edge_index = random_pairs(N_DRUGS, 20, gen)
    else:
        for r in range(N_SIDE_EFFECTS):
            g[side_effect_edge(r)].edge_index = random_pairs(N_DRUGS, 3 + 2 * r, gen)
    return SimpleNamespace(x=g[DRUG].x, proteins=proteins, n_proteins=N_PROTEINS,
                           n_side_effects=N_SIDE_EFFECTS, eval_graph=g)


def test_plain_dense_adjacency_is_the_mean_over_each_relations_neighbours():
    g = small_graphs("per_side_effect").eval_graph
    A = dense_adjacency(g, N_SIDE_EFFECTS, over_relations=False)
    assert A.shape == (N_SIDE_EFFECTS, N_DRUGS, N_DRUGS)
    for r in range(N_SIDE_EFFECTS):
        source, target = g[side_effect_edge(r)].edge_index
        expected = torch.zeros(N_DRUGS, N_DRUGS)
        expected[target, source] = 1
        expected /= expected.sum(1, keepdim=True).clamp(min=1)
        assert torch.equal(A[r], expected)


def test_the_models_adjacency_also_takes_the_mean_over_relations():
    g = small_graphs("per_side_effect").eval_graph
    plain = dense_adjacency(g, N_SIDE_EFFECTS, over_relations=False)
    A = dense_adjacency(g, N_SIDE_EFFECTS)
    k = (plain.sum(2) > 0).sum(0)  # relations in which each drug has a neighbour
    assert (k > 1).any()
    assert torch.allclose(A, plain / k.clamp(min=1)[None, :, None])
    assert torch.allclose(A.sum((0, 2)), (k > 0).float())  # each drug's weights sum to 1


@pytest.mark.parametrize("proteins", [True, False])
def test_dense_rgcn_matches_rgcnconv(proteins):
    """Same bases, mean aggregation and root weight: the rewrite is RGCNConv, two layers."""
    graphs = small_graphs("per_side_effect", proteins)
    g = graphs.eval_graph
    torch.manual_seed(0)
    encoder = RGCN(graphs, D, dropout=0.0).eval()

    # The homogeneous graph: drugs 0..N_DRUGS-1, then proteins; relations in the layer's order.
    edges = [g[side_effect_edge(r)].edge_index for r in range(N_SIDE_EFFECTS)]
    if proteins:
        shift = {DRUG: 0, PROTEIN: N_DRUGS}
        for t in (PPI, TARGETS, TARGETED_BY):
            edges.append(g[t].edge_index + torch.tensor([[shift[t[0]]], [shift[t[2]]]]))
    edge_type = torch.cat([torch.full((e.shape[1],), r) for r, e in enumerate(edges)])
    edge_index = torch.cat(edges, 1)
    x = encoder.inputs(g)
    x = torch.cat([x[DRUG], x[PROTEIN]]) if proteins else x[DRUG]

    h = x
    for k, layer in enumerate(encoder.layers):
        conv = RGCNConv(D, D, len(edges), num_bases=N_BASES, aggr="mean", root_weight=True)
        with torch.no_grad():
            conv.weight.copy_(layer.basis)
            conv.comp.copy_(layer.comp)
            conv.root.copy_(layer.root)
            conv.bias.copy_(layer.bias)
        h = conv(h, edge_index, edge_type)
        if k == 0:
            h = h.relu()
    with torch.no_grad():
        plain = dense_adjacency(g, N_SIDE_EFFECTS, over_relations=False)
        assert torch.allclose(encoder(g, plain), h[:N_DRUGS], atol=1e-5)


@pytest.mark.parametrize("name", GNNS)
@pytest.mark.parametrize("proteins", [True, False])
def test_every_gnn_gives_one_embedding_per_drug_with_or_without_proteins(name, proteins):
    Encoder = ENCODERS[name]
    graphs = small_graphs(Encoder.drug_edges, proteins)
    encoder = Encoder(graphs, D, dropout=0.0)
    z = encoder(graphs.eval_graph)
    assert z.shape == (N_DRUGS, D)
    z.sum().backward()
    assert encoder.inputs.drug.weight.grad.abs().sum() > 0
    if proteins:
        assert encoder.inputs.protein.weight.grad.abs().sum() > 0
    else:
        assert encoder.inputs.protein is None


@pytest.mark.parametrize("name", GNNS)
def test_dropout_acts_in_training_only(name):
    Encoder = ENCODERS[name]
    graphs = small_graphs(Encoder.drug_edges)
    encoder = Encoder(graphs, D, dropout=0.5)
    g = graphs.eval_graph
    encoder.train()
    assert not torch.equal(encoder(g), encoder(g))
    encoder.eval()
    assert torch.equal(encoder(g), encoder(g))


@pytest.mark.parametrize("name", GNNS)
def test_a_drug_with_no_edges_still_sees_its_own_input(name):
    """A new drug (cold-drug test) without targets has no messages; its z must still follow x."""
    Encoder = ENCODERS[name]
    graphs = small_graphs(Encoder.drug_edges, proteins=False)
    g = graphs.eval_graph
    lonely = N_DRUGS  # one more drug, in no edge
    g[DRUG].x = torch.cat([g[DRUG].x, torch.zeros(1, N_IN)])
    graphs.x = g[DRUG].x
    encoder = Encoder(graphs, D, dropout=0.0).eval()
    z0 = encoder(g)[lonely]
    g[DRUG].x[lonely] = 1.0
    assert not torch.allclose(encoder(g)[lonely], z0)


def test_gat_refuses_d_not_divisible_by_its_heads():
    with pytest.raises(ValueError, match="divisible"):
        ENCODERS["gat"](small_graphs("collapsed"), 10, dropout=0.0)
