"""Encoders: a message graph → drug embeddings `z_drug [645, d]`.

Spec: `docs/spec-model.md` §3 (input layer), §4.1 (encoders).

Every encoder is built as `Encoder(graphs, d, dropout)` from the fold's `MessageGraphs`
and called as `encoder(graph) → z_drug`, so the training loop never needs to know
which one it runs. `drug_edges` and `uses_proteins` say which message graph it needs.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

from sief_druco.graph import DRUG


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


ENCODERS = {"none": NoEncoder}
