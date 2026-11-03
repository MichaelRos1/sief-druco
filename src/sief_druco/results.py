"""Results tables: the run files of a sweep → CSV and Markdown under `local/results/`.

Spec: `docs/spec-model.md` §7 (uncertainty, subgroups), §8 (repeats), §6 (paired runs).

- `runs.csv`: one row per finished run: its config, git commit, sampling variant, timing
  and every test number of its `metrics.json` (flattened; the 95% cluster-bootstrap
  interval as `*_lo` / `*_hi`).
- `summary.csv`: per split and model, mean, SD (spread) and count over the runs.
- `paired.csv`: per split, each model against each other one, run by run: only runs with
  the same fold and seed are paired; mean and SD of the difference, and how often the first
  model is ahead.
- `results.md`: the above per split, with each model's paired difference to the split's
  reference model (the highest mean validation macro AUPRC).

A model is a config without its split, fold and seed: `Config.name`.
"""

import json
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd

from sief_druco import decagon
from sief_druco.train import RUNS_DIR

RESULTS_DIR = decagon.LOCAL_DIR / "results"
KEYS = ["split", "model", "fold", "seed"]
SELECTION = "val_macro_auprc"  # §7: pair-level validation macro AUPRC

# Columns of the summary and paired tables, with their Markdown labels. Only those a split
# has are shown (1:1 numbers exist on the random split, two-new-drug ones on cold-drug).
METRICS = {
    "val_macro_auprc": "val macro AUPRC",
    "test_macro_auprc": "macro AUPRC",
    "test_macro_auroc": "macro AUROC",
    "test_macro_ap50": "macro AP@50",
    "test_micro_auprc": "micro AUPRC",
    "test_subgroups_targets_macro_auprc": "macro AUPRC, 284 with targets",
    "test_subgroups_no_targets_macro_auprc": "macro AUPRC, 361 without",
    "test_uniform_negatives_macro_auprc": "1:1 macro AUPRC",
    "test_uniform_negatives_macro_auroc": "1:1 macro AUROC",
    "test_uniform_negatives_macro_ap50": "1:1 macro AP@50",
    "test_two_new_micro_auprc": "two-new micro AUPRC",
}
COSTS = {"best_epoch": "best epoch", "minutes": "minutes"}


def _flatten(d: dict, prefix: str = "") -> dict:
    out = {}
    for key, value in d.items():
        name = f"{prefix}{key}"
        if isinstance(value, dict):
            out.update(_flatten(value, f"{name}_"))
        elif isinstance(value, list):  # a bootstrap interval
            out[f"{name}_lo"], out[f"{name}_hi"] = value
        else:
            out[name] = value
    return out


def run_rows(configs, runs_dir: Path) -> pd.DataFrame:
    """One row per finished run among `configs`; runs not done yet are left out."""
    rows = []
    for config in configs:
        out = config.run_dir(runs_dir)
        if not (out / "metrics.json").exists():
            continue
        meta = json.loads((out / "config.json").read_text())
        summary = json.loads((out / "metrics.json").read_text())
        test = summary.pop("test")
        row = {"split": config.split, "model": config.name, "fold": config.fold,
               "seed": config.seed, **{k: v for k, v in meta.items() if k not in KEYS + ["name"]},
               **summary, "minutes": summary["seconds"] / 60}
        row.update(_flatten(test["test"], "test_"))
        for name, part in test.items():
            if name != "test":
                row.update(_flatten(part, f"{name}_"))
        rows.append(row)
    return pd.DataFrame(rows)


def summary(runs: pd.DataFrame) -> pd.DataFrame:
    """Per split and model: mean, SD and n of every metric and cost, models in run order."""
    cols = [c for c in (*METRICS, *COSTS) if c in runs]
    grouped = runs.groupby(["split", "model"], sort=False)[cols]
    out = pd.concat({"mean": grouped.mean(), "sd": grouped.std(), "n": grouped.count()}, axis=1)
    out.columns = [f"{metric}_{stat}" for stat, metric in out.columns]
    return out.reset_index()


def paired(runs: pd.DataFrame) -> pd.DataFrame:
    """Each ordered pair of models within a split, on runs with the same fold and seed.

    `diff` = model_a − model_b per run; `a_ahead` counts the paired runs where it is > 0.
    """
    cols = [c for c in METRICS if c in runs]
    rows = []
    for split, runs_s in runs.groupby("split", sort=False):
        by_model = {m: r.set_index(["fold", "seed"])[cols] for m, r in runs_s.groupby("model", sort=False)}
        for a, b in ((a, b) for a in by_model for b in by_model if a != b):
            diff = (by_model[a] - by_model[b]).dropna(how="all")  # index join: matching runs only
            for metric in cols:
                d = diff[metric].dropna()
                if len(d):
                    rows.append({"split": split, "model_a": a, "model_b": b, "metric": metric,
                                 "n": len(d), "mean_diff": d.mean(), "sd_diff": d.std(),
                                 "a_ahead": int((d > 0).sum())})
    return pd.DataFrame(rows, columns=["split", "model_a", "model_b", "metric", "n",
                                       "mean_diff", "sd_diff", "a_ahead"])


# --- Markdown ---


def _md(df: pd.DataFrame) -> str:
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def _pm(mean: float, sd: float, digits: int) -> str:
    if np.isnan(mean):
        return ""
    return f"{mean:.{digits}f}" if np.isnan(sd) else f"{mean:.{digits}f} ± {sd:.{digits}f}"


def _split_md(split: str, runs: pd.DataFrame, summ: pd.DataFrame, pairs: pd.DataFrame) -> str:
    runs, summ, pairs = (t[t["split"] == split] for t in (runs, summ, pairs))
    shown = {c: label for c, label in METRICS.items() if c in runs and runs[c].notna().any()}
    table = pd.DataFrame({"model": summ["model"], "runs": summ[f"{SELECTION}_n"]})
    for c, label in {**shown, **COSTS}.items():
        digits = 4 if c in METRICS else (0 if c == "best_epoch" else 1)
        table[label] = [_pm(m, s, digits) for m, s in zip(summ[f"{c}_mean"], summ[f"{c}_sd"])]
    parts = [f"## {split}", "Mean ± SD over runs (test rows unless named).", _md(table)]

    reference = summ.loc[summ[f"{SELECTION}_mean"].idxmax(), "model"]
    vs = pairs[pairs["model_b"] == reference]
    if len(vs):
        diff = pd.DataFrame({"model": vs["model_a"].unique()})
        for c, label in shown.items():
            d = vs[vs["metric"] == c].set_index("model_a")
            diff[f"Δ {label}"] = [
                f"{_pm(d.at[m, 'mean_diff'], d.at[m, 'sd_diff'], 4)} ({d.at[m, 'a_ahead']}/{d.at[m, 'n']})"
                if m in d.index else "" for m in diff["model"]]
        parts += [f"Paired difference to `{reference}` (highest mean {METRICS[SELECTION]}): "
                  "mean ± SD of model − reference over runs with the same fold and seed "
                  "(runs where the model is ahead / paired runs).", _md(diff)]

    per_run = runs[["model", "fold", "seed", "best_epoch"]].rename(columns={"best_epoch": "best epoch"})
    auprc = runs["test_macro_auprc"].map("{:.4f}".format)
    lo, hi = "test_bootstrap_95_macro_auprc_lo", "test_bootstrap_95_macro_auprc_hi"
    if lo in runs:  # absent when a sweep sets n_resamples = 0
        auprc += [f" [{l:.4f}, {h:.4f}]" for l, h in zip(runs[lo], runs[hi])]
    per_run["macro AUPRC [95% bootstrap]"] = auprc
    parts += ["Per run: test macro AUPRC with its 95% cluster-bootstrap interval over test drugs.",
              _md(per_run)]
    return "\n\n".join(parts)


def write(configs, out_dir: Path, runs_dir: Path = RUNS_DIR, title: str = "") -> Path:
    """Write `runs.csv`, `summary.csv`, `paired.csv` and `results.md` for the finished runs."""
    runs = run_rows(configs, runs_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    header = [f"# Results: {title}", f"{len(runs)} of {len(configs)} runs finished; written {date.today()}."]
    if runs.empty:
        (out_dir / "results.md").write_text("\n\n".join(header) + "\n")
        return out_dir
    summ, pairs = summary(runs), paired(runs)
    runs.to_csv(out_dir / "runs.csv", index=False)
    summ.to_csv(out_dir / "summary.csv", index=False)
    pairs.to_csv(out_dir / "paired.csv", index=False)
    commits = ", ".join(f"`{str(c)[:7]}`" + (" (dirty)" if dirty else "") for c, dirty in
                        runs[["git_commit", "git_dirty"]].drop_duplicates().itertuples(index=False))
    header.append(f"Git commits: {commits}. Sampling: {', '.join(runs['sampling'].unique())}.")
    sections = [_split_md(s, runs, summ, pairs) for s in runs["split"].unique()]
    (out_dir / "results.md").write_text("\n\n".join(header + sections) + "\n")
    return out_dir
