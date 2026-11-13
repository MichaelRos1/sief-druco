"""Cross-fitting, Platt calibration and the score table handed to the burden track.

Spec: `docs/spec-model.md` §1 (hand-off files), §9 (cross-fitting, Platt, reporting).

Two regimes, 10 folds each; run k scores fold k, validates on fold k + 1 (wrapping) for
early stopping and Platt, and trains on the other 8, so nothing is scored by a model that
trained on it:

- **U** (unseen pair), pair-level protocol: all 207,690 pairs in 10 folds. Its scores
  serve the unlabelled pairs; its report covers the full rows of every held-out pair.
- **S** (seen pair, missing side effect), random-split protocol: the triples of labelled
  pairs, positives and zeros, in 10 folds per side effect. The pair's other side effects
  stay in the message graph. Unlabelled pairs' triples are training zeros in every S run:
  S scores only labelled pairs, so its validation slice, and with it Platt, holds only
  their triples.

Per fold model and regime, Platt per side effect `p = σ(a_r · logit + c_r)` is fitted on
the validation slice: all side effects at once, each by its own mean BCE (LBFGS).

Everything goes to `local/crossfit/`: the fold assignments in `folds/` (drawn once, as in
`splits`), and per model `<Config.name>/` with the fold runs (`runs/`, as in `train`), the
calibrated held-out scores of each fold (`U/fold<k>.npz`, `S/fold<k>.npz`), the hand-off
files and `report/`. `side_effects.csv` gives each side effect's held-out AUPRC per regime
as the mean over the 10 fold models (Platt keeps each run's ranking while `a_r > 0`).
Resumable: a done run or fold file is skipped, so the same command
finishes the rest. Run with `python -m sief_druco.crossfit` (the model's settings as in
`train`, e.g. `--max_epochs 10`).
"""

import argparse
import json
from dataclasses import asdict, replace
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import torch.nn.functional as F
from scipy.special import expit

from sief_druco import decagon, splits, train
from sief_druco.graph import MessageGraphs
from sief_druco.splits import TEST, TRAIN, VAL
from sief_druco.train import Config

CROSSFIT_DIR = decagon.LOCAL_DIR / "crossfit"
REGIMES = {"U": "pair", "S": "random"}  # regime → the split protocol it trains with
N_CROSSFIT_FOLDS = 10
FOLD_SEEDS = {"U": 30, "S": 31}  # the splits use 0–2, 10–12 and 20–22
SEED = 0  # one seed per fold model (#9)
ALWAYS_TRAIN = -1  # S: the triples of unlabelled pairs, in no fold
N_BINS = 15  # ECE: equal-count bins
# Terciles of side effects by labelled pairs: one grouping for every run, ranked as by any
# run's training positives (80% of them).
GROUPS = ("rare", "median", "common")

# Calibrated probabilities stay strictly inside (0, 1) in float32.
_LOW = np.finfo(np.float32).tiny
_HIGH = np.nextafter(np.float32(1), np.float32(0))


# --- Folds ---


def u_folds(n_pairs: int, rng: np.random.Generator) -> np.ndarray:
    """int8 fold per pair: all pairs dealt into 10 folds in random order."""
    folds = np.empty(n_pairs, dtype=np.int8)
    folds[rng.permutation(n_pairs)] = np.arange(n_pairs) % N_CROSSFIT_FOLDS
    return folds


def s_folds(labels: sp.spmatrix, rng: np.random.Generator) -> np.ndarray:
    """int8 fold per triple `[n_pairs, n_side_effects]`.

    Per side effect, the labelled pairs' positives and their zeros are each dealt into 10
    folds in random order. Unlabelled pairs get `ALWAYS_TRAIN`.
    """
    labels = sp.csc_matrix(labels)
    n_pairs, n_side_effects = labels.shape
    labelled = np.flatnonzero(sp.csr_matrix(labels).getnnz(axis=1) > 0)
    folds = np.full((n_side_effects, n_pairs), ALWAYS_TRAIN, dtype=np.int8)
    for r in range(n_side_effects):
        is_pos = np.zeros(n_pairs, dtype=bool)
        is_pos[labels.indices[labels.indptr[r] : labels.indptr[r + 1]]] = True
        for group in (labelled[is_pos[labelled]], labelled[~is_pos[labelled]]):
            folds[r, rng.permutation(group)] = np.arange(len(group)) % N_CROSSFIT_FOLDS
    return np.ascontiguousarray(folds.T)


def load_folds(regime: str, data, folds_dir: Path) -> np.ndarray:
    """The saved fold assignment of `regime`; drawn from its seed and saved if missing."""
    path = folds_dir / f"{regime}.npz"
    if not path.exists():
        rng = np.random.default_rng(FOLD_SEEDS[regime])
        folds = u_folds(len(data.pairs), rng) if regime == "U" else s_folds(data.labels, rng)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp.npz")
        np.savez_compressed(tmp, folds=folds)
        tmp.replace(path)
    with np.load(path) as z:
        return z["folds"]


def roles(folds: np.ndarray, k: int) -> np.ndarray:
    """Roles of run k: fold k test (scored), fold k + 1 validation, the rest training."""
    out = np.full(folds.shape, TRAIN, dtype=np.uint8)
    out[folds == (k + 1) % N_CROSSFIT_FOLDS] = VAL
    out[folds == k] = TEST
    return out


def split_fold(regime: str, folds: np.ndarray, k: int) -> splits.Fold:
    """Run k of `regime` as a fold of its split protocol, for `train.run`."""
    return splits.Fold(REGIMES[regime], k, roles(folds, k), None, {"crossfit": regime})


def held_out(folds: np.ndarray, k: int, n_columns: int) -> np.ndarray:
    """Boolean `[n_pairs, n_columns]` of the triples run k scores."""
    is_k = folds == k
    return np.broadcast_to(is_k[:, None], (len(folds), n_columns)) if is_k.ndim == 1 else is_k


def out_of_fold(folds: np.ndarray, values, n_columns: int) -> np.ndarray:
    """float32 `[n_pairs, n_columns]` from each run's held-out scores (NaN where none).

    `values` gives run k's scores in row-major order of `held_out(folds, k, …)`.
    """
    out = np.full((len(folds), n_columns), np.nan, dtype=np.float32)
    for k, v in enumerate(values):
        out[held_out(folds, k, n_columns)] = v
    return out


# --- Platt scaling and calibration metrics ---


def fit_platt(logits, labels, columns, n_columns: int, device="cpu", max_iter: int = 200):
    """`(a, c)` per column minimising each column's mean BCE of `σ(a_r · logit + c_r)`.

    Inputs are flat, one entry per triple, `columns` its side effect. All columns at once
    with LBFGS; the loss is a sum of independent per-column terms, so each column gets its
    own minimum. Starts at the identity, a = 1, c = 0.
    """
    device = torch.device(device)
    dtype = torch.float64 if device.type == "cpu" else torch.float32  # MPS has no float64
    x = torch.as_tensor(np.asarray(logits), dtype=dtype, device=device)
    y = torch.as_tensor(np.asarray(labels), dtype=dtype, device=device)
    r = torch.as_tensor(np.asarray(columns), dtype=torch.long, device=device)
    w = (1 / torch.bincount(r, minlength=n_columns).clamp(min=1).to(dtype))[r]
    a = torch.ones(n_columns, dtype=dtype, device=device, requires_grad=True)
    c = torch.zeros(n_columns, dtype=dtype, device=device, requires_grad=True)
    optimizer = torch.optim.LBFGS([a, c], max_iter=max_iter, line_search_fn="strong_wolfe",
                                  tolerance_grad=1e-10, tolerance_change=1e-14)

    def closure():
        optimizer.zero_grad()
        loss = F.binary_cross_entropy_with_logits(a[r] * x + c[r], y, w, reduction="sum")
        loss.backward()
        return loss

    optimizer.step(closure)
    return a.detach().cpu().double().numpy(), c.detach().cpu().double().numpy()


def platt(logits, a, c, columns) -> np.ndarray:
    """float32 `σ(a_r · logit + c_r)` per entry, strictly inside (0, 1)."""
    z = np.asarray(a)[columns] * np.asarray(logits, dtype=np.float64) + np.asarray(c)[columns]
    return np.clip(expit(z).astype(np.float32), _LOW, _HIGH)


def ece(p, y, n_bins: int = N_BINS) -> tuple[float, list[dict]]:
    """Expected calibration error over `n_bins` equal-count bins, and the bins themselves.

    Each bin: its size, mean predicted probability and observed positive rate; the
    reliability diagram plots the last two.
    """
    p, y = np.asarray(p).ravel(), np.asarray(y).ravel()
    bins = []
    for idx in np.array_split(np.argsort(p), n_bins):
        bins.append({"n": len(idx), "predicted": float(p[idx].mean(dtype=np.float64)),
                     "observed": float(y[idx].mean(dtype=np.float64))})
    value = sum(b["n"] * abs(b["predicted"] - b["observed"]) for b in bins) / len(p)
    return float(value), bins


def calibration_report(p, y, n_pos, raw=None) -> tuple[dict, np.ndarray, np.ndarray]:
    """ECE and reliability bins of `p` against `y` `[n, n_side_effects]`, pooled and per
    group of side effects by training positives `n_pos`; predicted and observed positives
    per side effect. `raw`: the uncalibrated probabilities, for their pooled ECE and bins."""
    value, bins = ece(p, y)
    report = {"ece": value, "bins": bins, "groups": {}}
    if raw is not None:
        report["ece_raw"], report["bins_raw"] = ece(raw, y)
    order = np.argsort(n_pos, kind="stable")
    for name, cols in zip(GROUPS, np.array_split(order, len(GROUPS))):
        value, bins = ece(p[:, cols], y[:, cols])
        report["groups"][name] = {"side_effects": cols.tolist(), "ece": value, "bins": bins}
    return report, p.sum(axis=0, dtype=np.float64), y.sum(axis=0, dtype=np.int64)


QUANTILES = (0.0, 0.01, 0.25, 0.5, 0.75, 0.99, 1.0)


def _quantiles(x) -> dict:
    return {str(q): float(v) for q, v in zip(QUANTILES, np.quantile(x, QUANTILES))}


def platt_summary(a: np.ndarray, c: np.ndarray) -> dict:
    """The fitted Platt parameters of one regime, `[n_folds, n_side_effects]` each:
    quantiles, and how many `a_r` are not positive (those would flip a ranking)."""
    return {"a": _quantiles(a), "c": _quantiles(c), "a_not_positive": int((a <= 0).sum()),
            "n": int(a.size)}


def score_summary(scores: np.ndarray, y: np.ndarray) -> dict:
    """The score table's rows of one regime against their labels: score quantiles over all
    triples and over the positive ones, and predicted positives per pair (mean row sum)."""
    out = {"pairs": len(scores), "all": _quantiles(scores),
           "mean_predicted_per_pair": float(scores.sum(axis=1, dtype=np.float64).mean()),
           "mean_observed_per_pair": float(y.sum(axis=1, dtype=np.float64).mean())}
    if y.any():
        out["positives"] = _quantiles(scores[y.astype(bool)])
    return out


def plot_reliability(reports: dict, path: Path) -> None:
    """Reliability diagrams, one panel per regime: pooled, the three groups and the pooled
    raw (pre-Platt) probabilities, log axes."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, len(reports), figsize=(5.5 * len(reports), 5), squeeze=False)
    for ax, (regime, report) in zip(axes[0], reports.items()):
        curves = [("pooled", report)] + [(g, report["groups"][g]) for g in GROUPS]
        for name, curve in curves:
            ax.plot([b["predicted"] for b in curve["bins"]], [b["observed"] for b in curve["bins"]],
                    marker="o", label=f"{name} (ECE {curve['ece']:.4f})")
        if "bins_raw" in report:
            ax.plot([b["predicted"] for b in report["bins_raw"]],
                    [b["observed"] for b in report["bins_raw"]], marker="x", linestyle=":",
                    color="grey", label=f"pooled, before Platt (ECE {report['ece_raw']:.4f})")
        ends = [1e-4, 1]
        ax.plot(ends, ends, "k--", linewidth=0.8)
        ax.set(xscale="log", yscale="log", xlim=ends, ylim=ends, xlabel="mean predicted p",
               ylabel="observed positive rate", title=f"{regime} regime, {N_BINS} equal-count bins")
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


# --- The hand-off files ---


def combine(u: np.ndarray, s: np.ndarray, labelled: np.ndarray) -> np.ndarray:
    """The score table: S scores for labelled pairs, U for the others; never labels."""
    return np.where(labelled[:, None], s, u).astype(np.float32, copy=False)


def write_handoff(out: Path, data, scores: np.ndarray, side_effects: pd.DataFrame,
                  meta: dict) -> None:
    """`scores.npy`, `pairs.csv`, `labels.npz`, `side_effects.csv`, `meta.json` (§1).

    Rows of the first three in the canonical pair order (`decagon.all_pairs`), columns in
    side-effect index order.
    """
    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "scores.npy", scores.astype(np.float32, copy=False))
    labelled = data.labels.getnnz(axis=1) > 0
    a, b = data.pairs.T
    pd.DataFrame({
        "pair_id": np.arange(len(data.pairs)),
        "stitch_a": data.drug_ids[a],
        "stitch_b": data.drug_ids[b],
        "labelled": labelled.astype(np.uint8),
        "regime": np.where(labelled, "S", "U"),
    }).to_csv(out / "pairs.csv", index=False)
    sp.save_npz(out / "labels.npz", sp.csr_matrix(data.labels))
    side_effects.to_csv(out / "side_effects.csv", index=False)
    (out / "meta.json").write_text(json.dumps(meta, indent=2) + "\n")


# --- Running ---


def _load_model(config: Config, data, f: splits.Fold, runs_dir: Path, device):
    """A finished run's best weights, and the evaluation graph of its fold."""
    Encoder = train.MODELS[config.encoder]
    graphs = MessageGraphs(data, f, Encoder.drug_edges, Encoder.uses_proteins and config.proteins)
    base_rate = torch.full((data.labels.shape[1],), 0.5)  # replaced by the saved weights
    model = train.build(config, graphs, base_rate)
    state = torch.load(config.run_dir(runs_dir) / "weights.pt", map_location="cpu")
    model.load_state_dict(state)
    return model.to(device), graphs.eval_graph.to(device)


def score_fold(config: Config, data, f: splits.Fold, runs_dir: Path, path: Path, device) -> None:
    """Fit Platt on the run's validation slice; write its calibrated held-out scores.

    `path` holds `logits` and `p` of the test triples in row-major order (see
    `out_of_fold`), and the fitted `a`, `c`.
    """
    model, graph = _load_model(config, data, f, runs_dir, device)
    pairs = torch.from_numpy(data.pairs).to(device)
    shift = train.sampling_shift(config.beta)

    def triples(role):
        ids, mask = train.eval_rows(f, role)
        x = train.score(model, graph, pairs[torch.from_numpy(ids).to(device)], config.chunk, shift)
        x, y = x.cpu().numpy(), data.labels[ids].toarray()
        if mask is None:  # pair protocol: full rows
            mask = np.ones(x.shape, dtype=bool)
        return x[mask], y[mask], np.nonzero(mask)[1]

    a, c = fit_platt(*triples(VAL), data.labels.shape[1], device)
    x, _, cols = triples(TEST)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, logits=x.astype(np.float32), p=platt(x, a, c, cols), a=a, c=c)
    tmp.replace(path)  # last and whole: the fold counts as done once this exists


def run(model: dict, data=None, root: Path = CROSSFIT_DIR, device=None, log=print) -> Path:
    """Cross-fit `model` (`Config` settings) in both regimes; write the hand-off and report."""
    data = data if data is not None else decagon.load()
    device = train._device(device)
    template = Config("pair", 0, SEED, **{**model, "n_resamples": 0})
    train.check(template)
    out = root / template.name
    runs_dir = out / "runs"
    configs = {}
    for regime, split in REGIMES.items():
        folds = load_folds(regime, data, root / "folds")
        configs[regime] = [replace(template, split=split, fold=k) for k in range(N_CROSSFIT_FOLDS)]
        for k, config in enumerate(configs[regime]):
            path = out / regime / f"fold{k}.npz"
            if path.exists():
                continue
            log(f"[{regime} {k + 1}/{N_CROSSFIT_FOLDS}] {config.label}")
            f = split_fold(regime, folds, k)
            if not config.is_done(runs_dir):
                train.run(config, data, runs_dir=runs_dir, device=device, fold=f)
            score_fold(config, data, f, runs_dir, path, device)
    log("writing the score table and the report")
    write_outputs(out, data, root, configs, template)
    return out


def _fold_values(out: Path, regime: str, key: str):
    for k in range(N_CROSSFIT_FOLDS):
        with np.load(out / regime / f"fold{k}.npz") as z:
            yield z[key]


def _out_of_fold(out: Path, regime: str, folds: np.ndarray, key: str, n_columns: int):
    return out_of_fold(folds, _fold_values(out, regime, key), n_columns)


def write_outputs(out: Path, data, root: Path, configs: dict, template: Config) -> None:
    n_columns = data.labels.shape[1]
    labels = data.labels.toarray()  # uint8, ~200 MB
    labelled = data.labels.getnnz(axis=1) > 0
    n_pos = data.labels.getnnz(axis=0)
    rows = {"U": slice(None), "S": labelled}  # the pairs each regime's report covers
    side_effects = pd.DataFrame({"cui": data.side_effect_ids, "name": data.side_effect_names,
                                 "labelled_pairs": n_pos})
    reports, calibrated = {}, {}
    for regime in REGIMES:
        folds = load_folds(regime, data, root / "folds")
        y = labels[rows[regime]]
        raw = expit(_out_of_fold(out, regime, folds, "logits", n_columns)[rows[regime]])
        p = _out_of_fold(out, regime, folds, "p", n_columns)
        reports[regime], predicted, observed = calibration_report(p[rows[regime]], y, n_pos, raw)
        del raw
        calibrated[regime] = p
        a, c = (np.stack(list(_fold_values(out, regime, key))) for key in ("a", "c"))
        reports[regime]["platt"] = platt_summary(a, c)
        auprc = [np.load(c.run_dir(out / "runs") / "metrics.npz")["test_auprc"] for c in configs[regime]]
        side_effects[f"auprc_{regime}"] = np.nanmean(auprc, axis=0)
        side_effects[f"predicted_{regime}"] = predicted
        side_effects[f"observed_{regime}"] = observed
    scores = combine(calibrated["U"], calibrated["S"], labelled)
    del calibrated
    for regime, rows in (("U", ~labelled), ("S", labelled)):  # the rows each regime serves
        reports[regime]["score_table"] = score_summary(scores[rows], labels[rows])

    runs = {regime: [{"fold": c.fold, **{k: v for k, v in _summary(c, out).items()
                                          if k in ("best_epoch", "epochs", "val_macro_auprc", "seconds")}}
                     for c in configs[regime]] for regime in REGIMES}
    commit, dirty = train._git_commit()
    meta = {
        "git_commit": commit,
        "git_dirty": dirty,
        "date": date.today().isoformat(),
        "model": {k: v for k, v in asdict(template).items() if k not in ("split", "fold")},
        "name": template.name,
        "sampling": train.sampling_variant(template.beta),
        "seed": SEED,
        "folds": N_CROSSFIT_FOLDS,
        "fold_rule": "run k scores fold k, validates on fold k + 1 (wrapping), trains on the rest",
        "fold_seeds": FOLD_SEEDS,
        "regimes": {
            "U": "pair-level protocol, all pairs in 10 folds; scores unlabelled pairs",
            "S": "random-split protocol, labelled pairs' triples in 10 folds per side effect; "
                 "scores every triple of labelled pairs",
        },
        "calibration": "Platt per side effect, fitted per fold model and regime on its validation slice",
        "runs": runs,
    }
    write_handoff(out, data, scores, side_effects, meta)
    report_dir = out / "report"
    report_dir.mkdir(exist_ok=True)
    (report_dir / "calibration.json").write_text(json.dumps(reports, indent=2) + "\n")
    plot_reliability(reports, report_dir / "reliability.png")


def _summary(config: Config, out: Path) -> dict:
    return json.loads((config.run_dir(out / "runs") / "metrics.json").read_text())


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description="Cross-fit a model and write the score table.")
    train.add_config_arguments(parser, skip=("split", "fold", "seed", "n_resamples"))
    out = run(vars(parser.parse_args(argv)))
    print(f"score table and report in {out}")


if __name__ == "__main__":
    main()
