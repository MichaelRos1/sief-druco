import json
from types import SimpleNamespace

import numpy as np
import pytest
import torch
import torch.nn.functional as F
from torch_geometric.data import HeteroData

from sief_druco import decagon, results, splits, sweep, train
from sief_druco.baselines import ADDITIVE, BASELINES, FREE_TABLE
from sief_druco.decoders import DECODERS
from sief_druco.graph import DRUG, MessageGraphs
from sief_druco.splits import TRAIN
from sief_druco.train import Config

needs_raw = pytest.mark.skipif(
    not decagon.RAW_DIR.is_dir(), reason=f"raw Decagon files not in {decagon.RAW_DIR}"
)

# Every baseline with each decoder it can take: tensor factorisation uses the chosen one.
MODELS = [("tensor_factorisation", dec) for dec in DECODERS] + [(m, "decagon") for m in ADDITIVE]


# --- Shape, symmetry, gradients (no data needed) ---

N_DRUGS, N_IN, N_SIDE_EFFECTS, D = 30, 12, 963, 8


def tiny(model: str, decoder: str, dropout=0.0, seed=0):
    """A baseline on random drug input, with random labels over all pairs."""
    torch.manual_seed(seed)
    graph = HeteroData()
    graph[DRUG].x = torch.rand(N_DRUGS, N_IN)
    graphs = SimpleNamespace(x=graph[DRUG].x, n_side_effects=N_SIDE_EFFECTS)
    config = Config("pair", 0, 0, encoder=model, decoder=decoder, d=D, dropout=dropout)
    net = train.build(config, graphs, torch.full((N_SIDE_EFFECTS,), 0.1))
    with torch.no_grad():  # move the free table off its all-zero start, so the test can see asymmetry
        for p in net["encoder"].parameters():
            p.add_(torch.randn_like(p))
    pairs = torch.from_numpy(decagon.all_pairs(N_DRUGS))
    labels = (torch.rand(len(pairs), N_SIDE_EFFECTS) < 0.1).to(torch.uint8)
    return net, graph, pairs, labels


@pytest.mark.parametrize("model,decoder", MODELS)
def test_baselines_are_symmetric_and_give_full_rows(model, decoder):
    net, graph, pairs, _ = tiny(model, decoder)
    ab = train.score(net, graph, pairs)
    ba = train.score(net, graph, pairs.flip(1))
    assert ab.shape == (len(pairs), N_SIDE_EFFECTS)
    torch.testing.assert_close(ab, ba)


@pytest.mark.parametrize("model", ADDITIVE)
def test_additive_models_add_the_two_drugs_rows(model):
    net, graph, pairs, _ = tiny(model, "decagon")
    u = net["encoder"].eval()(graph)
    a, b = pairs.T
    torch.testing.assert_close(train.score(net, graph, pairs), u[a] + u[b] + net["decoder"].b)
    r = torch.randint(N_SIDE_EFFECTS, (len(pairs),))
    torch.testing.assert_close(net["decoder"].triples(u[a], u[b], r), (u[a] + u[b] + net["decoder"].b)[torch.arange(len(r)), r])


@pytest.mark.parametrize("model", ADDITIVE)
def test_additive_models_start_at_the_base_rate(model):
    graph = HeteroData()
    graph[DRUG].x = torch.rand(N_DRUGS, N_IN)
    graphs = SimpleNamespace(x=graph[DRUG].x, n_side_effects=5)
    base_rate = torch.tensor([0.01, 0.02, 0.05, 0.1, 0.3])
    net = train.build(Config("pair", 0, 0, encoder=model), graphs, base_rate)
    logits = train.score(net, graph, torch.from_numpy(decagon.all_pairs(N_DRUGS)))
    torch.testing.assert_close(torch.sigmoid(logits), base_rate.expand_as(logits))


@pytest.mark.parametrize("model,decoder", MODELS)
def test_chunked_detached_gradients_equal_unchunked(model, decoder):
    net, graph, pairs, labels = tiny(model, decoder)
    ids = torch.arange(0, len(pairs), 2)

    net.zero_grad()
    net.train()
    z = net["encoder"](graph)
    a, b = pairs[ids].T
    F.binary_cross_entropy_with_logits(net["decoder"](z[a], z[b]), labels[ids].float()).backward()
    expected = {k: p.grad.clone() for k, p in net.named_parameters()}

    net.zero_grad()
    train.train_step(net, graph, ids, pairs, labels, chunk=17)
    for k, p in net.named_parameters():
        torch.testing.assert_close(p.grad, expected[k], rtol=1e-4, atol=1e-7, msg=k)


@pytest.mark.parametrize("model,decoder", MODELS)
def test_sampled_step_gives_gradients_for_every_parameter(model, decoder):
    net, graph, pairs, labels = tiny(model, decoder)
    net.zero_grad()
    assert train.train_step(net, graph, torch.arange(len(pairs)), pairs, labels, beta=0.2, chunk=50) > 0
    assert all(p.grad is not None for p in net.parameters())


def test_dropout_reaches_the_models_that_have_it():
    for model in ("tensor_factorisation", "additive_mono"):
        net, graph, _, _ = tiny(model, "decagon", dropout=0.5)
        net.train()
        assert not torch.equal(net["encoder"](graph), net["encoder"](graph)), model
        net.eval()
        assert torch.equal(net["encoder"](graph), net["encoder"](graph)), model


# --- Configs, sweeps and tables (no data needed) ---


def test_run_names_hold_only_the_settings_a_model_has():
    assert Config("pair", 0, 0, encoder="tensor_factorisation", d=32).name == \
        "tensor_factorisation_decagon_w-none_d32_lr0.01_do0"
    assert Config("pair", 0, 0, encoder="additive_mono", decoder="mlp", dropout=0.3).name == \
        "additive_mono_w-none_lr0.01_do0.3"
    assert Config("pair", 0, 0, encoder="additive_free", lr=3e-3).name == "additive_free_w-none_lr0.003"


@pytest.mark.parametrize("model", FREE_TABLE)
def test_free_table_models_refuse_cold_drug(model, tmp_path):
    with pytest.raises(ValueError, match="not applicable on cold-drug"):
        train.check(Config("cold_drug", 0, 0, encoder=model))
    calls = []
    with pytest.raises(ValueError, match="not applicable"):  # before the first run of the sweep
        sweep.run_sweep([Config("pair", 0, 0), Config("cold_drug", 0, 0, encoder=model)], tmp_path,
                        data="unused", run=lambda c, *_, **__: calls.append(c), log=lambda _: None)
    assert calls == []
    train.check(Config("pair", 0, 0, encoder=model))
    train.check(Config("random", 0, 0, encoder=model))


def test_additive_mono_runs_on_cold_drug_and_the_free_table_has_no_dropout():
    train.check(Config("cold_drug", 0, 0, encoder="additive_mono"))
    with pytest.raises(ValueError, match="no dropout"):
        train.check(Config("pair", 0, 0, encoder="additive_free", dropout=0.3))


def write_fake_run(config, runs_dir, val):
    out = config.run_dir(runs_dir)
    out.mkdir(parents=True)
    (out / "config.json").write_text(json.dumps(
        {"name": config.name, "git_commit": "abc1234", "git_dirty": False, "sampling": "full rows"}))
    (out / "metrics.json").write_text(json.dumps(
        {"best_epoch": 10, "seconds": 60.0, "val_macro_auprc": val, "test": {"test": {"macro_auprc": val}}}))


def test_the_tables_report_free_table_models_as_not_applicable_on_cold_drug(tmp_path):
    configs = [Config(split, 0, 0, encoder=model) for split in splits.SPLITS
               for model in ("none", *BASELINES) if not (split == "cold_drug" and model in FREE_TABLE)]
    for i, c in enumerate(configs):
        write_fake_run(c, tmp_path / "runs", 0.1 + i / 100)
    md = (results.write(configs, tmp_path / "results", tmp_path / "runs") / "results.md").read_text()
    sections = dict(part.split("\n", 1) for part in md.split("## ")[1:])
    names = [Config("pair", 0, 0, encoder=m).name for m in FREE_TABLE]
    assert f"Not applicable on cold-drug, not run (free per-drug parameters learn nothing for a new drug): " \
           f"`{names[0]}`, `{names[1]}`." in sections["cold_drug"]
    assert "Not applicable" not in sections["pair"] + sections["random"]


# --- The real dataset ---


@pytest.fixture(scope="module")
def data():
    return decagon.load()


@needs_raw
@pytest.mark.parametrize("decoder", DECODERS)
def test_tensor_factorisation_starts_near_the_training_base_rate(data, decoder):
    fold = splits.load("pair", 0, data)
    graphs = MessageGraphs(data, fold, proteins=False)
    train_ids = np.flatnonzero(fold.roles == TRAIN)
    base_rate = torch.from_numpy(np.asarray(data.labels[train_ids].mean(axis=0)).ravel()).float()
    torch.manual_seed(0)
    net = train.build(Config("pair", 0, 0, encoder="tensor_factorisation", decoder=decoder), graphs, base_rate)
    logits = train.score(net, graphs.eval_graph, torch.from_numpy(data.pairs[train_ids]))
    error = (torch.sigmoid(logits).mean(0) / base_rate - 1).abs()
    assert error.median() < 0.05  # the bounds of the no-encoder test in test_train.py
    assert error.max() < 0.25


@needs_raw
@pytest.mark.parametrize("model", BASELINES)
def test_a_short_run_of_each_baseline_writes_its_files(data, model, tmp_path):
    config = Config("pair", 0, seed=0, encoder=model, max_epochs=2, n_resamples=2)
    summary = train.run(config, data, runs_dir=tmp_path)
    out = config.run_dir(tmp_path)
    assert {p.name for p in out.iterdir()} == {
        "config.json", "log.csv", "weights.pt", "metrics.npz", "metrics.json"}
    assert json.loads((out / "metrics.json").read_text())["epochs"] == 2
    assert 0 < summary["test"]["test"]["macro_auprc"] < 1
