"""Model input for one fold: the drug input and the message graphs.

Spec: `docs/spec-model.md` §2 (Graph), §3 (Input layer), §4.3 (no-protein ablation) and
§5 (message vs supervision edges).

The drug input is recomputed from each split's training drugs. Drug–drug edges come only
from labelled training pairs (random split: training triples), so validation and test pairs
never enter a message graph. Reverse edges are added here, after splitting: every edge type
between nodes of the same type holds both directions, and targets get a reverse type.
"""

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import torch
from torch_geometric.data import HeteroData

from sief_druco import decagon
from sief_druco.splits import TRAIN

MIN_MONO_DRUGS = 2  # keep a mono column only if ≥ 2 training drugs have it
MESSAGE_PERCENT = 70  # pair-level and cold-drug: share of training pairs drawn as message pairs

DRUG, PROTEIN = "drug", "protein"
PPI = (PROTEIN, "ppi", PROTEIN)
TARGETS = (DRUG, "targets", PROTEIN)
TARGETED_BY = (PROTEIN, "targeted_by", DRUG)
LABELLED_WITH = (DRUG, "labelled_with", DRUG)  # the collapsed drug–drug type
DRUG_EDGES = ("per_side_effect", "collapsed")  # RGCN; GraphSAGE and GAT


def side_effect_edge(r: int) -> tuple[str, str, str]:
    """The drug–drug edge type of side effect `r` in a per-side-effect graph."""
    return (DRUG, f"se{r}", DRUG)


def drug_input(mono: sp.spmatrix, train_drugs: np.ndarray) -> np.ndarray:
    """Float32 `[n_drugs, n_cols + 1]`: the mean of a drug's mono indicators, then `log(1 + n)`.

    The columns are the mono side effects of ≥ 2 training drugs (`train_drugs`: indices or
    a mask); the mean is over the drug's side effects among them. n counts all of the drug's
    mono side effects, kept columns or not. A drug with no profile gets all zeros.
    """
    mono = sp.csr_matrix(mono, dtype=np.float32)
    per_column = np.asarray(mono[train_drugs].sum(axis=0)).ravel()
    x = mono[:, per_column >= MIN_MONO_DRUGS].toarray()
    n_kept = x.sum(axis=1, keepdims=True)
    mean = np.divide(x, n_kept, out=np.zeros_like(x), where=n_kept > 0)
    n = np.asarray(mono.sum(axis=1), dtype=np.float32)
    return np.hstack([mean, np.log1p(n)])


@dataclass
class Epoch:
    graph: HeteroData
    supervision: np.ndarray  # pair IDs whose full rows enter the loss
    message: np.ndarray | None  # pair IDs drawn as message pairs; None on the random split


class MessageGraphs:
    """The message graphs of one fold: a fresh one per training epoch, one for evaluation.

    Pair-level and cold-drug: each `epoch` redraws the training pairs 70 / 30 into message
    and supervision pairs; only labelled message pairs become edges. `eval_graph` holds
    every labelled training pair. Random split: every training triple is an edge and every
    pair a supervision pair (the loss masks the held-out triples), so each epoch reuses
    `eval_graph`.

    `drug_edges` is "per_side_effect" (963 drug–drug types, for RGCN) or "collapsed" (one
    "labelled together" type). `proteins=False` drops proteins, PPI and targets (the
    no-protein ablation). Tensors are built on the CPU.
    """

    def __init__(self, data, fold, drug_edges: str = "collapsed", proteins: bool = True):
        if drug_edges not in DRUG_EDGES:
            raise ValueError(f"unknown drug_edges {drug_edges!r}; expected one of {DRUG_EDGES}")
        self.split = fold.split
        self.drug_edges = drug_edges
        self.proteins = proteins
        self.n_side_effects = data.labels.shape[1]
        self._pairs = data.pairs

        if fold.split == "cold_drug":
            train_drugs = fold.drug_roles == TRAIN
        else:
            train_drugs = np.arange(decagon.N_DRUGS)
        self.x = torch.from_numpy(drug_input(data.mono, train_drugs))

        if proteins:
            self.n_proteins = len(data.protein_ids)
            self._ppi = _both_ways(*data.ppi.T)
            self._targets = torch.from_numpy(data.targets.T.astype(np.int64))

        # Labelled triples sorted by side effect, so each side effect is one slice.
        labels = sp.csc_matrix(data.labels)
        rows = labels.indices.astype(np.int64)
        cols = np.repeat(np.arange(self.n_side_effects), np.diff(labels.indptr))
        if fold.split == "random":
            self.train_pairs = None
            keep = fold.roles[rows, cols] == TRAIN
        else:
            self.train_pairs = np.flatnonzero(fold.roles == TRAIN)
            keep = fold.roles[rows] == TRAIN
        self._rows, self._cols = rows[keep], cols[keep]
        self.eval_graph = self._graph(self._rows, self._cols)

    def epoch(self, rng: np.random.Generator) -> Epoch:
        if self.train_pairs is None:
            return Epoch(self.eval_graph, np.arange(decagon.N_PAIRS), None)
        order = rng.permutation(self.train_pairs)
        k = len(order) * MESSAGE_PERCENT // 100
        message, supervision = np.sort(order[:k]), np.sort(order[k:])
        is_message = np.zeros(decagon.N_PAIRS, dtype=bool)
        is_message[message] = True
        keep = is_message[self._rows]
        return Epoch(self._graph(self._rows[keep], self._cols[keep]), supervision, message)

    def _graph(self, rows: np.ndarray, cols: np.ndarray) -> HeteroData:
        """The graph whose drug–drug edges are the triples `(rows, cols)`, sorted by `cols`."""
        g = HeteroData()
        g[DRUG].x = self.x
        if self.proteins:
            g[PROTEIN].num_nodes = self.n_proteins
            g[PPI].edge_index = self._ppi
            g[TARGETS].edge_index = self._targets
            g[TARGETED_BY].edge_index = self._targets.flip(0)
        if self.drug_edges == "collapsed":
            g[LABELLED_WITH].edge_index = _both_ways(*self._pairs[np.unique(rows)].T)
        else:
            a, b = self._pairs[rows].T
            ends = np.cumsum(np.bincount(cols, minlength=self.n_side_effects))
            for r, (s, e) in enumerate(zip(np.r_[0, ends[:-1]], ends)):
                g[side_effect_edge(r)].edge_index = _both_ways(a[s:e], b[s:e])
        return g


def _both_ways(a: np.ndarray, b: np.ndarray) -> torch.Tensor:
    """Edge index `[2, 2n]` with each edge `a → b` and its reverse `b → a`."""
    return torch.from_numpy(np.stack([np.r_[a, b], np.r_[b, a]]).astype(np.int64))
