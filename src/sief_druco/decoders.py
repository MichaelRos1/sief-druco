"""The three decoders: drug embeddings `(z_a, z_b)` → one logit per side effect.

Spec: `docs/spec-model.md` §4.2; the bias and start values: #8.

Every decoder is symmetric in `a, b` and adds a per-side-effect bias `b_r`, started at the
logit of side effect r's training base rate. `forward(z_a, z_b)` gives full rows
`[n, n_side_effects]`; `triples(z_a, z_b, r)` gives one logit per `(pair, side effect)`, for
the sampling switch of §5, which scores only the drawn triples.
"""

import torch
import torch.nn as nn


class Decoder(nn.Module):
    def __init__(self, base_rate: torch.Tensor):
        super().__init__()
        self.b = nn.Parameter(torch.logit(base_rate.double()).float())

    def triples(self, z_a, z_b, r) -> torch.Tensor:
        """Default: the full row, then one entry, so the sampling switch saves no work (MLP head)."""
        return self(z_a, z_b).gather(1, r[:, None]).squeeze(1)


class Decagon(Decoder):
    """`z_aᵀ D_r R D_r z_b + b_r`, D_r diagonal, R = (R₀ + R₀ᵀ) / 2 shared.

    One matmul, `((z_a ⊗ z_b) ⊙ R).flatten(1) @ (D ⊗ D)ᵀ`: 8× faster on MPS than
    broadcasting to `[n, n_side_effects, d]`. Start values D_r ~ N(0, 1), R₀ ~ N(0, 1/√d)
    (standard deviation); smaller ones stall training.
    """

    def __init__(self, d: int, base_rate: torch.Tensor):
        super().__init__(base_rate)
        self.D = nn.Parameter(torch.randn(len(base_rate), d))
        self.R0 = nn.Parameter(torch.randn(d, d) / d**0.5)

    @property
    def R(self) -> torch.Tensor:
        return (self.R0 + self.R0.T) / 2

    def forward(self, z_a, z_b):
        DD = (self.D[:, :, None] * self.D[:, None, :]).flatten(1)
        return ((z_a[:, :, None] * z_b[:, None, :]) * self.R).flatten(1) @ DD.T + self.b

    def triples(self, z_a, z_b, r):
        D = self.D[r]
        return ((z_a * D) @ self.R * (z_b * D)).sum(1) + self.b[r]


class DistMult(Decoder):
    """`Σ_k z_a[k] · w_r[k] · z_b[k] + b_r`. Start values w_r ~ N(0, 1), as the Decagon
    decoder's D_r (the spec fixes start values for the Decagon decoder only)."""

    def __init__(self, d: int, base_rate: torch.Tensor):
        super().__init__(base_rate)
        self.w = nn.Parameter(torch.randn(len(base_rate), d))

    def forward(self, z_a, z_b):
        return (z_a * z_b) @ self.w.T + self.b

    def triples(self, z_a, z_b, r):
        return (z_a * self.w[r] * z_b).sum(1) + self.b[r]


class MLPHead(Decoder):
    """`MLP([z_a + z_b, z_a ⊙ z_b]) + b_r`.

    The spec leaves the MLP's shape open: here one hidden layer of width d with ReLU, and
    no bias in the last layer, since `b_r` is the bias.
    """

    def __init__(self, d: int, base_rate: torch.Tensor):
        super().__init__(base_rate)
        self.mlp = nn.Sequential(
            nn.Linear(2 * d, d), nn.ReLU(), nn.Linear(d, len(base_rate), bias=False)
        )

    def forward(self, z_a, z_b):
        return self.mlp(torch.cat([z_a + z_b, z_a * z_b], dim=1)) + self.b


DECODERS = {"decagon": Decagon, "distmult": DistMult, "mlp": MLPHead}
