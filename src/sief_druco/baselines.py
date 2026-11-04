"""The baselines: tensor factorisation and the two per-drug additive models.

Spec: `docs/spec-model.md` §4.3; tuned like every other model: §8 step 2.

They plug into the training loop as encoders, `Model(graphs, d, dropout)` called as
`model(graph) → z_drug`, and ignore the graph:

- **Tensor factorisation**: free per-drug vectors `Embedding(645, d)`, then the chosen
  decoder. No mono input, no graph. Start values N(0, 0.1²), about the scale of the
  no-encoder model's start `z` (~0.07), so it too starts at the training base rate.
- **Per-drug additive**: `σ(u_r(a) + u_r(b))`, with `u` a free `[645, 963]` table (free
  table) or `u(drug) = Linear(drug input → 963)` (mono input). The "encoder" gives `u`;
  its own decoder `Additive` adds the two drugs' rows and the per-side-effect bias `b_r`
  of §4.2, which adds nothing a free table or a Linear bias couldn't learn. Both start
  `u` at 0 (no hidden units, so a zero start is no symmetry problem): with `b_r` at the
  base-rate logit, the untrained model predicts exactly the base rate.

The additive models have no `d` (ignored) and use no other decoder. The free table has no
dropout (it has no layer to drop out of); the mono-input model drops out its input.
Free-table models learn nothing for a drug they never saw in training, so on cold-drug
they are "not applicable" and not run.
"""

import torch
import torch.nn as nn
import torch.nn.functional as F

from sief_druco.decoders import Decoder
from sief_druco.graph import DRUG

FREE_TABLE = ("tensor_factorisation", "additive_free")  # cold-drug: "not applicable", not run
ADDITIVE = ("additive_free", "additive_mono")  # own decoder, no d
NO_DROPOUT = ("additive_free",)
TF_INIT_STD = 0.1


class TensorFactorisation(nn.Module):
    drug_edges = "collapsed"
    uses_proteins = False

    def __init__(self, graphs, d: int, dropout: float):
        super().__init__()
        self.z = nn.Embedding(graphs.x.shape[0], d)
        nn.init.normal_(self.z.weight, std=TF_INIT_STD)
        self.dropout = dropout

    def forward(self, graph) -> torch.Tensor:
        return F.dropout(self.z.weight, self.dropout, self.training)


class Additive(Decoder):
    """`u_r(a) + u_r(b) + b_r`, where the "embeddings" are the drugs' rows of `u`."""

    def __init__(self, d: int, base_rate: torch.Tensor):  # d: unused, the decoders' signature
        super().__init__(base_rate)

    def forward(self, u_a, u_b):
        return u_a + u_b + self.b

    def triples(self, u_a, u_b, r):
        rows = torch.arange(len(r), device=r.device)
        return u_a[rows, r] + u_b[rows, r] + self.b[r]


class AdditiveFree(nn.Module):
    drug_edges = "collapsed"
    uses_proteins = False
    decoder = Additive

    def __init__(self, graphs, d: int, dropout: float):
        super().__init__()
        self.u = nn.Parameter(torch.zeros(graphs.x.shape[0], graphs.n_side_effects))

    def forward(self, graph) -> torch.Tensor:
        return self.u


class AdditiveMono(nn.Module):
    drug_edges = "collapsed"
    uses_proteins = False
    decoder = Additive

    def __init__(self, graphs, d: int, dropout: float):
        super().__init__()
        self.u = nn.Linear(graphs.x.shape[1], graphs.n_side_effects, bias=False)  # b_r is the bias
        nn.init.zeros_(self.u.weight)
        self.dropout = dropout

    def forward(self, graph) -> torch.Tensor:
        return self.u(F.dropout(graph[DRUG].x, self.dropout, self.training))


BASELINES = {"tensor_factorisation": TensorFactorisation, "additive_free": AdditiveFree,
             "additive_mono": AdditiveMono}
