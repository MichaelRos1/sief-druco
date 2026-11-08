"""Sweeps: run many configs one after another, unattended, then write the results tables.

Spec: `docs/spec-model.md` §8 (procedure, repeats), §6 (same folds and seeds for every model).

A run counts as done once its `metrics.json` exists (`train` writes it last), so a sweep
skips done runs: killing it, or a closed lid overnight, loses at most the run in progress,
and the same command finishes the rest. A run that raises is logged and skipped; the sweep
goes on and lists the failures at the end. After the runs, `results` writes the tables to
`local/results/<sweep>/`.

The sweeps of §8 are defined below as Python lists of `Config`s. Tuning, Stage B and the
ablation depend on earlier stages; fill in the choices under "Choices from earlier stages"
by hand once a stage is done (#37, #39). Run one with `python -m sief_druco.sweep stage_a`
(`--dry-run` lists the runs and which are done, `--tables-only` rewrites the tables).
"""

import argparse
import itertools
import json
import sys
import time
import traceback
from dataclasses import replace
from pathlib import Path

from sief_druco import decagon, results, splits, train
from sief_druco.baselines import ADDITIVE, BASELINES, FREE_TABLE, NO_DROPOUT
from sief_druco.decoders import DECODERS
from sief_druco.train import Config

FOLDS = tuple(range(splits.N_FOLDS))
SEEDS = (0, 1, 2)  # the same folds and seeds for every model, so runs pair up (§6)

GNNS = ("rgcn", "sage", "gat")
ENCODERS = ("none", *GNNS)
GRID = {"lr": (3e-3, 1e-2), "d": (32, 64, 128), "dropout": (0.0, 0.3)}  # tuning, §8 step 2

# --- Choices from earlier stages: fill in by hand once each stage is done ---

DECODER: str | None = "mlp"  # Stage A (#37): MLP head, best pair-level val macro AUPRC
LOSS_WEIGHTS: str | None = "none"  # Stage A (#37): 1/√n_pos did not help
TUNED: dict[str, dict] = {  # tuning (#39): best pair-level val macro AUPRC per model
    "none": {"lr": 1e-2, "d": 128, "dropout": 0.3},
    "rgcn": {"lr": 3e-3, "d": 128, "dropout": 0.0},
    "sage": {"lr": 3e-3, "d": 128, "dropout": 0.0},
    "gat": {"lr": 3e-3, "d": 128, "dropout": 0.0},
    "tensor_factorisation": {"lr": 1e-2, "d": 128, "dropout": 0.0},
    "additive_free": {"lr": 1e-2},
    "additive_mono": {"lr": 1e-2, "dropout": 0.3},
}
BEST_GNN: str | None = None  # Stage B winner among the GNNs (#39), for the ablation


def _stage_a_choices() -> dict:
    if DECODER is None or LOSS_WEIGHTS is None:
        raise ValueError("set DECODER and LOSS_WEIGHTS in sweep.py from Stage A first (#37)")
    return {"decoder": DECODER, "loss_weights": LOSS_WEIGHTS}


def _tuned(model: str) -> dict:
    if model not in TUNED:
        raise ValueError(f"no tuned settings for {model!r}: fill in TUNED in sweep.py (#39)")
    return TUNED[model]


# --- The sweeps of §8 ---


def smoke() -> list[Config]:
    """The ticket's small end-to-end check: no encoder × 3 decoders × fold 0 × 2 seeds, 5 epochs."""
    return [Config("pair", 0, seed, decoder=decoder, max_epochs=5, n_resamples=20)
            for decoder in DECODERS for seed in SEEDS[:2]]


def stage_a() -> list[Config]:
    """Decoders × encoders None and RGCN, then the loss weights on encoder None; pair-level fold 0."""
    runs = [Config("pair", 0, seed, encoder=encoder, decoder=decoder)
            for encoder in ("none", "rgcn") for decoder in DECODERS for seed in SEEDS]
    runs += [Config("pair", 0, seed, decoder=decoder, loss_weights="inv_sqrt")
             for decoder in DECODERS for seed in SEEDS]
    return runs


def tuning() -> list[Config]:
    """Every model once over the settings that apply to it, on pair-level fold 0, seed 0."""
    runs = []
    for model in (*ENCODERS, *BASELINES):
        grid = {k: v for k, v in GRID.items()
                if not (k == "d" and model in ADDITIVE or k == "dropout" and model in NO_DROPOUT)}
        for values in itertools.product(*grid.values()):
            runs.append(Config("pair", 0, SEEDS[0], encoder=model, **_stage_a_choices(),
                               **dict(zip(grid, values))))
    return runs


def stage_b() -> list[Config]:
    """Every model at its tuned settings, 3 splits × 3 folds × 3 seeds."""
    return [Config(split, fold, seed, encoder=model, **_stage_a_choices(), **_tuned(model))
            for model in (*ENCODERS, *BASELINES)
            for split in splits.SPLITS if not (split == "cold_drug" and model in FREE_TABLE)
            for fold in FOLDS for seed in SEEDS]


def ablation() -> list[Config]:
    """The best GNN of Stage B on drug–drug edges only, 3 splits × 3 folds × 3 seeds."""
    if BEST_GNN is None:
        raise ValueError("set BEST_GNN in sweep.py from Stage B first (#39)")
    return [Config(split, fold, seed, encoder=BEST_GNN, proteins=False, **_stage_a_choices(),
                   **_tuned(BEST_GNN))
            for split in splits.SPLITS for fold in FOLDS for seed in SEEDS]


def ablation_tables() -> list[Config]:
    """The ablation's runs next to the same GNN's Stage B runs, so the tables pair them (§4.3)."""
    runs = ablation()
    return runs + [replace(c, proteins=True) for c in runs]


SWEEPS = {"smoke": smoke, "stage_a": stage_a, "tuning": tuning, "stage_b": stage_b,
          "ablation": ablation}
TABLES = {"ablation": ablation_tables}  # runs a sweep's tables cover, if more than its own


# --- Running ---


def _duration(seconds: float) -> str:
    minutes = round(seconds / 60)
    return f"{minutes // 60}h{minutes % 60:02d}m"


def run_sweep(configs: list[Config], runs_dir: Path = train.RUNS_DIR, data=None,
              run=train.run, log=print) -> list[Config]:
    """Run every config that isn't done yet, in order; return the ones that raised.

    All configs are checked before the first run, so a typo doesn't stop a sweep overnight.
    The ETA is the mean duration of the done runs (earlier ones included) × the runs left.
    """
    for config in configs:
        train.check(config)
    todo = [c for c in configs if not c.is_done(runs_dir)]
    seconds = [json.loads((c.run_dir(runs_dir) / "metrics.json").read_text())["seconds"]
               for c in configs if c.is_done(runs_dir)]
    log(f"{len(configs)} runs: {len(configs) - len(todo)} done, {len(todo)} to run")
    if todo and data is None:
        data = decagon.load()
    failed = []
    for i, config in enumerate(todo):
        eta = f", ETA {_duration(sum(seconds) / len(seconds) * (len(todo) - i))}" if seconds else ""
        log(f"[{i + 1}/{len(todo)}] {config.label}{eta}")
        started = time.time()
        try:
            run(config, data, runs_dir=runs_dir)
        except Exception:
            log(traceback.format_exc())
            failed.append(config)
            continue
        seconds.append(time.time() - started)
    for config in failed:
        log(f"failed: {config.label}")
    return failed


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description="Run a sweep of §8 and write its results tables.")
    parser.add_argument("sweep", choices=SWEEPS)
    parser.add_argument("--dry-run", action="store_true", help="list the runs, run nothing")
    parser.add_argument("--tables-only", action="store_true", help="only rewrite the tables")
    args = parser.parse_args(argv)
    configs = SWEEPS[args.sweep]()
    if args.dry_run:
        for c in configs:
            print(f"{'done' if c.is_done() else '    '}  {c.label}")
        print(f"{len(configs)} runs, {sum(c.is_done() for c in configs)} done")
        return
    failed = [] if args.tables_only else run_sweep(configs)
    tabled = TABLES[args.sweep]() if args.sweep in TABLES else configs
    out = results.write(tabled, results.RESULTS_DIR / args.sweep, title=args.sweep)
    print(f"tables in {out}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
