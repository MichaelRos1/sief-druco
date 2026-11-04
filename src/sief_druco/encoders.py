"""Encoders: a message graph → drug embeddings `z_drug [645, d]`.

Spec: `docs/spec-model.md` §3 (input layer), §4.1 (encoders), §5 (MPS notes).

Every encoder is built as `Encoder(graphs, d, dropout)` from the fold's `MessageGraphs`
and called as `encoder(graph) → z_drug`, so the training loop never needs to know
which one it runs. `drug_edges` and `uses_proteins` say which message graph it needs; with
`graphs.proteins` False (the no-protein ablation) the GNNs use drug–drug edges only.

The GNNs share the input layer of §3 (`Linear(n_cols + 1 → d)` on the drug input,
`Embedding(19089, d)` for proteins), then two message-passing layers with ReLU and dropout
between them.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv, SAGEConv, to_hetero

from sief_druco.graph import DRUG, PPI, PROTEIN, TARGETED_BY, TARGETS, side_effect_edge

N_BASES = 30  # RGCN: bases shared by all relations
GAT_HEADS = 4  # GAT: 4 heads × d/4


class NoEncoder(nn.Module):
    """No message passing: two Linear layers on the drug input, ReLU and dropout between.

    The first is the input layer `Linear(n_cols + 1 → d)` of §3. Tests whether the graph
    helps on the same input.
    """

    drug_edges = "collapsed"
    uses_proteins = False

    def __init__(self, graphs, d: int, dropout: float):
        super().__init__()
        self.input = nn.Linear(graphs.x.shape[1], d)
        self.out = nn.Linear(d, d)
        self.dropout = dropout

    def forward(self, graph) -> torch.Tensor:
        h = F.relu(self.input(graph[DRUG].x))
        return self.out(F.dropout(h, self.dropout, self.training))


class _Inputs(nn.Module):
    """The input layer of §3: drug input → d, and a free `d`-vector per protein."""

    def __init__(self, graphs, d: int):
        super().__init__()
        self.drug = nn.Linear(graphs.x.shape[1], d)
        self.protein = nn.Embedding(graphs.n_proteins, d) if graphs.proteins else None

    def forward(self, graph) -> dict[str, torch.Tensor]:
        x = {DRUG: self.drug(graph[DRUG].x)}
        if self.protein is not None:
            x[PROTEIN] = self.protein.weight
        return x


# --- RGCN: the dense drug–drug rewrite ---


def dense_adjacency(graph, n_side_effects: int, over_relations: bool = True) -> torch.Tensor:
    """`[n_side_effects, n_drugs, n_drugs]` weights of drug i's drug–drug messages from j.

    `over_relations=False`: `A[r, i, j] = 1 / deg_r(i)` if `j → i` is an edge of side effect
    r, else 0, PyG's mean aggregation at the target as a matrix. With `over_relations` (the
    model's), each row is further divided by `k_i`, the number of relations in which drug i
    has a neighbour, so a drug's drug–drug messages are a mean over its relations, not a sum
    of hundreds of terms (see `RGCN`). 1.6 GB float32 for Decagon, built in place.
    """
    n_drugs = graph[DRUG].num_nodes
    edges = [graph[side_effect_edge(r)].edge_index for r in range(n_side_effects)]
    counts = torch.tensor([e.shape[1] for e in edges], device=edges[0].device)
    relation = torch.repeat_interleave(torch.arange(n_side_effects, device=counts.device), counts)
    source, target = torch.cat(edges, 1)
    adjacency = torch.zeros(n_side_effects, n_drugs, n_drugs, device=source.device)
    adjacency[relation, target, source] = 1.0
    degree = adjacency.sum(2, keepdim=True)  # [n_side_effects, n_drugs, 1]
    if over_relations:
        degree *= (degree > 0).sum(0)  # × k_i
    return adjacency.div_(degree.clamp_(min=1))


def _mean_messages(x_source, edge_index, n_targets: int, weight) -> torch.Tensor:
    """Mean over each target's incoming `x_source @ weight`; zero for a target with none."""
    source, target = edge_index
    messages = (x_source @ weight)[source]
    total = torch.zeros(n_targets, weight.shape[1], device=x_source.device).index_add_(0, target, messages)
    degree = torch.zeros(n_targets, device=x_source.device).index_add_(
        0, target, torch.ones_like(target, dtype=x_source.dtype))
    return total / degree.clamp(min=1)[:, None]


class DenseRGCNLayer(nn.Module):
    """`RGCNConv` with basis decomposition, mean aggregation and root weight, rewritten.

    Relations: the drug–drug ones first (one per side effect), then PPI, targets and
    targeted-by. The drug–drug sum `Σ_r A_r x W_r = Σ_b (Σ_r comp[r, b] A_r) x basis_b`
    runs on the dense adjacency as two matmuls; the protein relations go by `index_add_`.
    Parameters have the shapes and start values of `RGCNConv`, so its weights load here.
    """

    def __init__(self, d: int, n_side_effects: int, proteins: bool, n_bases: int = N_BASES):
        super().__init__()
        self.n_side_effects = n_side_effects
        self.basis = nn.Parameter(torch.empty(n_bases, d, d))  # RGCNConv.weight
        self.comp = nn.Parameter(torch.empty(n_side_effects + 3 * proteins, n_bases))
        self.root = nn.Parameter(torch.empty(d, d))
        self.bias = nn.Parameter(torch.zeros(d))
        for p in (self.basis, self.comp, self.root):
            nn.init.xavier_uniform_(p)  # as RGCNConv (PyG's glorot)

    def forward(self, x, adjacency, protein_edges=None) -> dict[str, torch.Tensor]:
        x_drug = x[DRUG]
        n_bases, d, _ = self.basis.shape
        n_drugs = x_drug.shape[0]
        mixed = (self.comp[: self.n_side_effects].T @ adjacency.flatten(1)).view(n_bases, n_drugs, n_drugs)
        neighbours = (mixed @ x_drug).permute(1, 0, 2).reshape(n_drugs, n_bases * d)  # [n_drugs, B·d]
        out = {DRUG: neighbours @ self.basis.reshape(n_bases * d, d)}
        if protein_edges is not None:
            x_protein = x[PROTEIN]
            n_proteins = x_protein.shape[0]
            w_ppi, w_targets, w_targeted_by = torch.einsum(
                "rb,bio->rio", self.comp[self.n_side_effects :], self.basis)
            out[PROTEIN] = (_mean_messages(x_protein, protein_edges[PPI], n_proteins, w_ppi)
                            + _mean_messages(x_drug, protein_edges[TARGETS], n_proteins, w_targets))
            out[DRUG] = out[DRUG] + _mean_messages(
                x_protein, protein_edges[TARGETED_BY], n_drugs, w_targeted_by)
        return {k: v + x[k] @ self.root + self.bias for k, v in out.items()}


class RGCN(nn.Module):
    """Relational GCN over every relation (§4.1, the reference), on the dense rewrite.

    Two `DenseRGCNLayer`s (each exactly `RGCNConv`) on the dense adjacency, with one
    deviation: a drug's drug–drug messages are averaged over its relations (`k_i`), where
    RGCN sums them. A drug has neighbours in hundreds of the 963 relations, and Adam moves
    every basis coefficient by about the learning rate per step, so the plain sum shifts
    `z` several-fold in one step and diverges at lr 1e-2 (loss 0.098 → 4,125 in two
    steps); the mean trains stably.

    Not `RGCNConv` (11× slower, 11.8 GB) and never `FastRGCNConv` (~51 GB).
    """

    drug_edges = "per_side_effect"
    uses_proteins = True

    def __init__(self, graphs, d: int, dropout: float):
        super().__init__()
        self.n_side_effects = graphs.n_side_effects
        self.inputs = _Inputs(graphs, d)
        self.layers = nn.ModuleList(
            DenseRGCNLayer(d, graphs.n_side_effects, graphs.proteins) for _ in range(2))
        self.dropout = dropout

    def forward(self, graph, adjacency: torch.Tensor | None = None) -> torch.Tensor:
        """`adjacency`: from `dense_adjacency` (default); tests pass the plain-RGCN one."""
        if adjacency is None:
            adjacency = dense_adjacency(graph, self.n_side_effects)
        protein_edges = None
        if self.inputs.protein is not None:
            protein_edges = {t: graph[t].edge_index for t in (PPI, TARGETS, TARGETED_BY)}
        h = self.layers[0](self.inputs(graph), adjacency, protein_edges)
        h = {k: F.dropout(F.relu(v), self.dropout, self.training) for k, v in h.items()}
        return self.layers[1](h, adjacency, protein_edges)[DRUG]


# --- GraphSAGE and GAT: the collapsed graph through to_hetero ---


class _TwoLayers(nn.Module):
    """conv → ReLU → dropout → conv, on one edge type; `to_hetero` copies it per edge type."""

    def __init__(self, conv1, conv2, dropout: float):
        super().__init__()
        self.conv1, self.conv2 = conv1, conv2
        self.dropout = nn.Dropout(dropout)  # a module, so to_hetero's trace keeps train/eval

    def forward(self, x, edge_index):
        return self.conv2(self.dropout(F.relu(self.conv1(x, edge_index))), edge_index)


class _Collapsed(nn.Module):
    """A homogeneous two-layer GNN made heterogeneous: messages summed over edge types."""

    drug_edges = "collapsed"
    uses_proteins = True

    def __init__(self, graphs, d: int, dropout: float):
        super().__init__()
        self.inputs = _Inputs(graphs, d)
        layers = _TwoLayers(self.conv(d), self.conv(d), dropout)
        self.gnn = to_hetero(layers, graphs.eval_graph.metadata(), aggr="sum")

    def forward(self, graph) -> torch.Tensor:
        return self.gnn(self.inputs(graph), graph.edge_index_dict)[DRUG]


class GraphSAGE(_Collapsed):
    """`SAGEConv` (mean aggregation, root weight) on the collapsed graph."""

    @staticmethod
    def conv(d: int):
        return SAGEConv(d, d)


class GAT(_Collapsed):
    """`GATConv`, 4 heads × d/4 concatenated, on the collapsed graph.

    No self-loops (they are undefined on drug → protein edges), so `residual=True` gives
    each node its own term, as GraphSAGE's root weight does; without it a drug with no
    message edges and no targets (every cold-drug test drug without targets) gets the same
    `z` whatever its input. Plain edge-index tensors only: a sorted `EdgeIndex` or
    `SparseTensor` needs `segment_reduce`, which has no MPS kernel.
    """

    @staticmethod
    def conv(d: int):
        return GATConv(d, d // GAT_HEADS, heads=GAT_HEADS, add_self_loops=False, residual=True)


ENCODERS = {"none": NoEncoder, "rgcn": RGCN, "sage": GraphSAGE, "gat": GAT}
