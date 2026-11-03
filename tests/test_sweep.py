import json

import numpy as np
import pandas as pd
import pytest

from sief_druco import decagon, results, splits, sweep, train
from sief_druco.train import Config

needs_raw = pytest.mark.skipif(
    not decagon.RAW_DIR.is_dir(), reason=f"raw Decagon files not in {decagon.RAW_DIR}"
)


# --- Running and restarting (no data needed) ---


class Killed(BaseException):
    """Stands in for a kill: like KeyboardInterrupt, it isn't an Exception, so it isn't caught."""


def fake_run(calls, kill_at=None, fail=()):
    """A `train.run` stand-in that writes the run files `train` writes first and last."""

    def run(config, data, runs_dir):
        calls.append(config)
        out = config.run_dir(runs_dir)
        out.mkdir(parents=True, exist_ok=True)
        (out / "config.json").write_text("{}")
        if len(calls) == kill_at:
            raise Killed  # mid-run: config.json written, metrics.json not
        if config in fail:
            raise RuntimeError("MPS out of memory")
        (out / "metrics.json").write_text(json.dumps({"seconds": 60.0}))

    return run


CONFIGS = [Config("pair", 0, seed, decoder=decoder) for decoder in ("decagon", "mlp") for seed in (0, 1, 2)]


def test_a_killed_sweep_restarts_with_the_remaining_runs(tmp_path):
    calls = []
    with pytest.raises(Killed):
        sweep.run_sweep(CONFIGS, tmp_path, data="unused", run=fake_run(calls, kill_at=3), log=lambda _: None)
    assert calls == CONFIGS[:3]

    calls, lines = [], []
    failed = sweep.run_sweep(CONFIGS, tmp_path, data="unused", run=fake_run(calls), log=lines.append)
    assert calls == CONFIGS[2:]  # the killed run again, then the ones never started
    assert failed == []
    assert lines[0] == "6 runs: 2 done, 4 to run"
    assert lines[1].endswith("ETA 0h04m")  # 4 runs left at the 60 s of the done ones
    assert all(sweep.is_done(c, tmp_path) for c in CONFIGS)

    calls = []
    sweep.run_sweep(CONFIGS, tmp_path, data="unused", run=fake_run(calls), log=lambda _: None)
    assert calls == []


def test_a_failing_run_is_logged_and_the_sweep_goes_on(tmp_path):
    calls, lines = [], []
    failed = sweep.run_sweep(CONFIGS, tmp_path, data="unused", run=fake_run(calls, fail=CONFIGS[1:2]),
                             log=lines.append)
    assert calls == CONFIGS
    assert failed == CONFIGS[1:2]
    assert any("MPS out of memory" in line for line in lines)


def test_a_typo_stops_the_sweep_before_the_first_run(tmp_path):
    calls = []
    with pytest.raises(ValueError, match="rgnc"):
        sweep.run_sweep([*CONFIGS, Config("pair", 0, 0, encoder="rgnc")], tmp_path, data="unused",
                        run=fake_run(calls), log=lambda _: None)
    assert calls == []


# --- Sweep definitions ---


def test_stage_b_runs_every_model_on_the_same_folds_and_seeds(monkeypatch):
    monkeypatch.setattr(sweep, "DECODER", "decagon")
    monkeypatch.setattr(sweep, "LOSS_WEIGHTS", "none")
    monkeypatch.setattr(sweep, "TUNED", {m: {"lr": 3e-3} for m in (*sweep.ENCODERS, *sweep.BASELINES)})
    runs = sweep.stage_b()
    cells = {}
    for c in runs:
        cells.setdefault((c.encoder, c.split), set()).add((c.fold, c.seed))
    full = {(f, s) for f in sweep.FOLDS for s in sweep.SEEDS}
    assert all(cell == full for cell in cells.values())
    for model in (*sweep.ENCODERS, *sweep.BASELINES):
        expected = set(splits.SPLITS) - ({"cold_drug"} if model in sweep.FREE_TABLE else set())
        assert {split for m, split in cells if m == model} == expected
    assert all(c.lr == 3e-3 and c.decoder == "decagon" for c in runs)


def test_tuning_covers_the_grid_without_d_for_the_additive_models(monkeypatch):
    monkeypatch.setattr(sweep, "DECODER", "decagon")
    monkeypatch.setattr(sweep, "LOSS_WEIGHTS", "none")
    per_model = pd.Series([c.encoder for c in sweep.tuning()]).value_counts()
    assert all(per_model[m] == 12 for m in (*sweep.ENCODERS, "tensor_factorisation"))
    assert all(per_model[m] == 4 for m in sweep.NO_D)


def test_later_stages_refuse_to_run_before_their_choices_are_filled_in():
    for stage in (sweep.tuning, sweep.stage_b, sweep.ablation):
        with pytest.raises(ValueError, match="#37|#39"):
            stage()


def test_the_ablation_switch_names_its_runs():
    assert Config("pair", 0, 0, proteins=False).name.endswith("_noprot")
    assert "noprot" not in Config("pair", 0, 0).name


# --- Results tables ---


def runs_frame(rows):
    """`(model, fold, seed, value)` → a runs table with that value as every metric."""
    return pd.DataFrame([{"split": "pair", "model": m, "fold": f, "seed": s,
                          "val_macro_auprc": v, "test_macro_auprc": v} for m, f, s, v in rows])


def test_paired_differences_use_matching_fold_and_seed_only():
    runs = runs_frame([
        ("a", 0, 0, 0.30), ("a", 0, 1, 0.32), ("a", 1, 0, 0.40),  # (1, 0): b has no such run
        ("b", 0, 0, 0.20), ("b", 0, 1, 0.34), ("b", 0, 2, 0.90),  # (0, 2): a has no such run
    ])
    pairs = results.paired(runs)
    row = pairs[(pairs.model_a == "a") & (pairs.model_b == "b") & (pairs.metric == "test_macro_auprc")].iloc[0]
    assert row.n == 2
    assert row.mean_diff == pytest.approx((0.10 - 0.02) / 2)
    assert row.sd_diff == pytest.approx(np.std([0.10, -0.02], ddof=1))
    assert row.a_ahead == 1
    back = pairs[(pairs.model_a == "b") & (pairs.model_b == "a") & (pairs.metric == "test_macro_auprc")].iloc[0]
    assert back.mean_diff == pytest.approx(-row.mean_diff)
    # Unpaired means would give a different answer: (0.30 + 0.32 + 0.40) / 3 − (0.20 + 0.34 + 0.90) / 3.
    summ = results.summary(runs).set_index("model")
    assert summ.at["a", "test_macro_auprc_mean"] - summ.at["b", "test_macro_auprc_mean"] != pytest.approx(row.mean_diff)


def test_paired_differences_never_cross_splits():
    runs = pd.concat([runs_frame([("a", 0, 0, 0.3)]), runs_frame([("b", 0, 0, 0.2)]).assign(split="random")])
    assert results.paired(runs).empty


# --- End to end on the real dataset ---


@needs_raw
def test_a_short_run_sweep_writes_the_results_table(tmp_path):
    configs = [Config("pair", 0, seed, decoder=decoder, max_epochs=2, n_resamples=2)
               for decoder in ("decagon", "distmult") for seed in (0, 1)]
    assert sweep.run_sweep(configs, tmp_path / "runs", log=lambda _: None) == []
    out = results.write(configs, tmp_path / "results", tmp_path / "runs", title="smoke")
    assert {p.name for p in out.iterdir()} == {"runs.csv", "summary.csv", "paired.csv", "results.md"}

    runs = pd.read_csv(out / "runs.csv")
    assert len(runs) == 4
    assert {"git_commit", "sampling", "seconds", "mean_train_s_per_epoch",
            "test_bootstrap_95_macro_auprc_lo", "test_subgroups_targets_macro_auprc",
            "test_subgroups_no_targets_macro_auprc"} <= set(runs.columns)
    summ = pd.read_csv(out / "summary.csv")
    assert list(summ["val_macro_auprc_n"]) == [2, 2]
    assert (summ["test_macro_auprc_sd"] > 0).all()
    pairs = pd.read_csv(out / "paired.csv")
    assert set(pairs["n"]) == {2}
    md = (out / "results.md").read_text()
    assert "4 of 4 runs finished" in md and "Paired difference to" in md
