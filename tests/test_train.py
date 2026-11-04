import csv
import json
from types import SimpleNamespace

import numpy as np
import pytest
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import HeteroData

from sief_druco import decagon, splits, train
from sief_druco.decoders import DECODERS
from sief_druco.encoders import NoEncoder
from sief_druco.graph import DRUG, MessageGraphs
from sief_druco.splits import TRAIN

needs_raw = pytest.mark.skipif(
    not decagon.RAW_DIR.is_dir(), reason=f"raw Decagon files not in {decagon.RAW_DIR}"
)


# --- Chunked gradients (no data needed) ---

N_DRUGS, N_IN, N_SIDE_EFFECTS, D = 30, 12, 20, 8


def tiny(decoder: str, seed=0):
    """A no-encoder model on random drug input, with random labels over all pairs."""
    torch.manual_seed(seed)
    graph = HeteroData()
    graph[DRUG].x = torch.rand(N_DRUGS, N_IN)
    model = nn.ModuleDict({
        "encoder": NoEncoder(SimpleNamespace(x=graph[DRUG].x), D, dropout=0.0),
        "decoder": DECODERS[decoder](D, torch.full((N_SIDE_EFFECTS,), 0.1)),
    })
    pairs = torch.from_numpy(decagon.all_pairs(N_DRUGS))
    labels = (torch.rand(len(pairs), N_SIDE_EFFECTS) < 0.1).to(torch.uint8)
    return model, graph, pairs, labels


def grads(model):
    return {k: p.grad.clone() for k, p in model.named_parameters()}


@pytest.mark.parametrize("decoder", DECODERS)
@pytest.mark.parametrize("masked", [False, True])
def test_chunked_detached_gradients_equal_unchunked(decoder, masked):
    model, graph, pairs, labels = tiny(decoder)
    ids = torch.arange(0, len(pairs), 2)  # a supervision subset
    weights = torch.rand(N_SIDE_EFFECTS) + 0.5
    mask = torch.rand(labels.shape) < 0.7 if masked else None

    # Reference: one plain forward and backward of the mean weighted loss.
    model.zero_grad()
    z = model["encoder"](graph)
    a, b = pairs[ids].T
    loss = F.binary_cross_entropy_with_logits(
        model["decoder"](z[a], z[b]), labels[ids].float(), reduction="none") * weights
    if mask is not None:
        loss = loss[mask[ids]]
    loss.mean().backward()
    expected = grads(model)

    model.zero_grad()
    value = train.train_step(model, graph, ids, pairs, labels, weights, mask, chunk=17)
    assert value == pytest.approx(loss.mean().item(), rel=1e-5)
    for k, g in grads(model).items():
        torch.testing.assert_close(g, expected[k], rtol=1e-4, atol=1e-7, msg=k)


def test_sampling_switch_draws_every_positive_and_beta_of_the_zeros():
    torch.manual_seed(0)
    y = (torch.rand(400, 250) < 0.05).float()
    drawn = train.draw_triples(y, beta=0.2)
    assert drawn[y == 1].all()
    assert drawn[y == 0].float().mean().item() == pytest.approx(0.2, abs=0.01)
    keep = torch.rand(y.shape) < 0.5
    assert not train.draw_triples(y, 0.2, keep)[~keep].any()


def test_sampling_shift_is_the_spec_correction():
    beta = 0.1
    p_s = torch.linspace(0.01, 0.99, 50, dtype=torch.float64)
    corrected = torch.sigmoid(torch.logit(p_s) + train.sampling_shift(beta))
    torch.testing.assert_close(corrected, beta * p_s / (beta * p_s - p_s + 1))
    assert train.sampling_shift(None) == 0


@pytest.mark.parametrize("decoder", DECODERS)
def test_sampled_step_gives_gradients_for_every_parameter(decoder):
    model, graph, pairs, labels = tiny(decoder)
    model.zero_grad()
    loss = train.train_step(model, graph, torch.arange(len(pairs)), pairs, labels, beta=0.2, chunk=50)
    assert loss > 0
    assert all(p.grad is not None for p in model.parameters())


# --- The real dataset ---


@pytest.fixture(scope="module")
def data():
    return decagon.load()


@needs_raw
@pytest.mark.parametrize("decoder", DECODERS)
def test_mean_probability_at_start_is_the_training_base_rate(data, decoder):
    fold = splits.load("pair", 0, data)
    graphs = MessageGraphs(data, fold, proteins=False)
    train_ids = np.flatnonzero(fold.roles == TRAIN)
    labels = data.labels[train_ids]
    base_rate = torch.from_numpy(np.asarray(labels.mean(axis=0)).ravel()).float()
    torch.manual_seed(0)
    model = nn.ModuleDict({
        "encoder": NoEncoder(graphs, 64, dropout=0.0),
        "decoder": DECODERS[decoder](64, base_rate),
    })
    logits = train.score(model, graphs.eval_graph, torch.from_numpy(data.pairs[train_ids]))
    error = (torch.sigmoid(logits).mean(0) / base_rate - 1).abs()
    # z has a nonzero mean, so the decoder adds a small fixed offset per side effect; without
    # the bias every prediction would be 0.5, ~25× the typical base rate.
    assert error.median() < 0.05
    assert error.max() < 0.25


@needs_raw
@pytest.mark.parametrize("split", splits.SPLITS)
def test_a_short_run_writes_its_files(data, split, tmp_path):
    config = train.Config(split, 0, seed=0, max_epochs=2, n_resamples=2)
    summary = train.run(config, data, runs_dir=tmp_path)
    out = config.run_dir(tmp_path)
    assert {p.name for p in out.iterdir()} == {
        "config.json", "log.csv", "weights.pt", "metrics.npz", "metrics.json"}
    assert json.loads((out / "metrics.json").read_text())["epochs"] == 2
    assert 0 < summary["test"]["test"]["macro_auprc"] < 1
    if split == "random":
        assert "test_uniform_negatives" in summary["test"]
    if split == "cold_drug":
        assert "micro_auprc" in summary["test"]["test_two_new"]


@needs_raw
@pytest.mark.parametrize("encoder", ("rgcn", "sage", "gat"))
def test_a_short_run_of_each_gnn_improves_validation(data, encoder, tmp_path):
    config = train.Config("pair", 0, seed=0, encoder=encoder, max_epochs=8, n_resamples=0)
    summary = train.run(config, data, runs_dir=tmp_path)
    with open(config.run_dir(tmp_path) / "log.csv") as f:
        first = float(next(csv.DictReader(f))["val_macro_auprc"])
    assert summary["best_epoch"] > 1
    assert summary["val_macro_auprc"] > first
