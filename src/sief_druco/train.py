"""Training: one run = a config → best weights, test metrics and a per-epoch log.

Spec: `docs/spec-model.md` §5 (training), §6 (splits), §7 (evaluation rows, selection).

Full rows: every supervision pair enters the loss with all its side effects (the random
split masks its held-out triples). Full-batch Adam: each epoch is one encoder pass and one
step; the decoder runs in chunks of pairs and its gradient goes back exactly through a
detached `z`. Early stopping on validation macro AUPRC; the best epoch's weights are kept
and scored on the test rows with `metrics`.

Run files go to `local/runs/<split>/fold<k>/seed<s>/<name>/`; `metrics.json` is written
last, so a run is done once it exists. Run one with
`python -m sief_druco.train --split pair --fold 0 --seed 0`.
"""

import argparse
import copy
import csv
import json
import math
import subprocess
import time
from dataclasses import MISSING, asdict, dataclass, fields
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from sief_druco import decagon, metrics, splits
from sief_druco.decoders import DECODERS
from sief_druco.encoders import ENCODERS
from sief_druco.graph import MessageGraphs
from sief_druco.splits import TEST, TEST_TWO_NEW, TRAIN, VAL

RUNS_DIR = decagon.LOCAL_DIR / "runs"
LOSS_WEIGHTS = ("none", "inv_sqrt")  # per side effect: 1, or 1/√n_pos(r) normalised to mean 1
CHUNK = 8192  # pairs per decoder chunk; doesn't change the result


@dataclass
class Config:
    split: str
    fold: int
    seed: int
    encoder: str = "none"
    decoder: str = "decagon"
    d: int = 64
    lr: float = 1e-2
    dropout: float = 0.0
    loss_weights: str = "none"
    proteins: bool = True  # False: the no-protein ablation (§4.3), drug–drug edges only
    # Sampling switch (§5, fallback only): keep each zero with probability beta, redrawn
    # every epoch, and add log(beta) to the logits at evaluation. None = full rows.
    beta: float | None = None
    patience: int = 30
    max_epochs: int = 500
    chunk: int = CHUNK
    n_resamples: int = 1000  # cluster bootstrap of the test metrics

    @property
    def name(self) -> str:
        parts = [self.encoder, self.decoder, f"w-{self.loss_weights}", f"d{self.d}",
                 f"lr{self.lr:g}", f"do{self.dropout:g}"]
        defaults = Config(self.split, self.fold, self.seed)
        for key, tag in (("beta", "beta"), ("patience", "pat"), ("max_epochs", "ep")):
            if getattr(self, key) != getattr(defaults, key):
                parts.append(f"{tag}{getattr(self, key):g}")
        if not self.proteins:
            parts.append("noprot")
        return "_".join(parts)

    def run_dir(self, runs_dir: Path = RUNS_DIR) -> Path:
        return runs_dir / self.split / f"fold{self.fold}" / f"seed{self.seed}" / self.name


# --- One step and scoring ---


def sampling_shift(beta: float | None) -> float:
    """What the sampling switch adds to a logit at evaluation: log(beta) (0 without it).

    Equals the spec's `p = β·p_s / (β·p_s − p_s + 1)`: it multiplies the odds by beta.
    """
    return 0.0 if beta is None else math.log(beta)


def draw_triples(y: torch.Tensor, beta: float, keep: torch.Tensor | None = None) -> torch.Tensor:
    """Sampling switch: mask of every positive and each zero with probability `beta`.

    Within `keep` if given. Redrawn on every call, so every epoch.
    """
    drawn = (y == 1) | (torch.rand(y.shape, device=y.device) < beta)
    return drawn if keep is None else keep & drawn


def train_step(model, graph, ids, pairs, labels, weights=None, mask=None, beta=None,
               chunk: int = CHUNK) -> float:
    """Gradients of the mean (weighted) BCE over the pairs `ids`; returns the loss.

    One encoder pass; the decoder runs on chunks of pairs against a detached copy of `z`
    whose gradient is passed back through the encoder at the end, so the result equals
    the unchunked gradient. `mask` [n_pairs, n_side_effects] marks the triples that count
    (random split); `beta` draws zeros for the sampling switch. Loss sums are rescaled to
    the mean once the number of counted triples is known.
    """
    z = model["encoder"](graph)
    z_leaf = z.detach().requires_grad_()
    decoder = model["decoder"]
    total, n_terms = 0.0, 0
    for s in range(0, len(ids), chunk):
        rows = ids[s : s + chunk]
        drug_a, drug_b = pairs[rows].T
        y = labels[rows].float()
        keep = None if mask is None else mask[rows]
        if beta is None:
            logits = decoder(z_leaf[drug_a], z_leaf[drug_b])
            loss = F.binary_cross_entropy_with_logits(logits, y, reduction="none")
            if weights is not None:
                loss = loss * weights
            if keep is not None:
                loss = loss * keep
            n_terms += y.numel() if keep is None else int(keep.sum())
        else:
            row, r = draw_triples(y, beta, keep).nonzero().T
            logits = decoder.triples(z_leaf[drug_a[row]], z_leaf[drug_b[row]], r)
            loss = F.binary_cross_entropy_with_logits(logits, y[row, r], reduction="none")
            if weights is not None:
                loss = loss * weights[r]
            n_terms += len(r)
        loss = loss.sum()
        loss.backward()
        total += loss.item()
    for param in [z_leaf, *decoder.parameters()]:
        if param.grad is not None:
            param.grad /= n_terms
    z.backward(z_leaf.grad)
    return total / n_terms


@torch.no_grad()
def score(model, graph, ab: torch.Tensor, chunk: int = CHUNK, shift: float = 0.0) -> torch.Tensor:
    """Logits `[n, n_side_effects]` for the drug pairs `ab` [n, 2], in eval mode."""
    model.eval()
    z = model["encoder"](graph)
    chunks = (ab[s : s + chunk].T for s in range(0, len(ab), chunk))
    return torch.cat([model["decoder"](z[drug_a], z[drug_b]) for drug_a, drug_b in chunks]) + shift


def eval_rows(fold, role: int) -> tuple[np.ndarray, np.ndarray | None]:
    """Pair IDs evaluated for `role`, and on the random split the mask of its triples."""
    if fold.split == "random":
        is_role = fold.roles == role
        ids = np.flatnonzero(is_role.any(axis=1))
        return ids, is_role[ids]
    return np.flatnonzero(fold.roles == role), None


# --- A run ---


def _device(device=None) -> torch.device:
    return torch.device(device or ("mps" if torch.backends.mps.is_available() else "cpu"))


def _on(ids: np.ndarray, device) -> torch.Tensor:
    return torch.from_numpy(ids).to(device)


def _sync(device: torch.device) -> None:
    if device.type == "mps":
        torch.mps.synchronize()


def _git_commit() -> tuple[str | None, bool]:
    try:
        head = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True, cwd=Path(__file__).parent).stdout.strip()
        dirty = subprocess.run(["git", "status", "--porcelain"], capture_output=True, text=True,
                               cwd=Path(__file__).parent).stdout.strip() != ""
        return head, dirty
    except (OSError, subprocess.CalledProcessError):
        return None, False


def check(config: Config) -> None:
    """Raise ValueError if the config names a split, encoder, decoder or loss weight that doesn't exist."""
    for value, allowed in ((config.split, splits.SPLITS), (config.encoder, ENCODERS),
                           (config.decoder, DECODERS), (config.loss_weights, LOSS_WEIGHTS)):
        if value not in allowed:
            raise ValueError(f"unknown setting {value!r}; expected one of {list(allowed)}")


def run(config: Config, data=None, runs_dir: Path = RUNS_DIR, device=None) -> dict:
    """Train one config, keep the best epoch, score its test rows and write the run files."""
    check(config)
    device = _device(device)
    data = data if data is not None else decagon.load()
    fold = splits.load(config.split, config.fold, data)
    torch.manual_seed(config.seed)
    rng = np.random.default_rng(config.seed)
    started = time.time()

    Encoder = ENCODERS[config.encoder]
    graphs = MessageGraphs(data, fold, Encoder.drug_edges, Encoder.uses_proteins and config.proteins)
    pairs = torch.from_numpy(data.pairs).to(device)
    labels = torch.from_numpy(data.labels.toarray()).to(device)  # uint8, ~200 MB
    is_train = torch.from_numpy(fold.roles).to(device) == TRAIN
    train_mask = is_train if fold.split == "random" else None  # per triple: the loss mask
    if train_mask is None:
        is_train = is_train[:, None]  # per pair

    # Base rate and positives of each side effect over the training triples.
    n_pos = (labels.bool() & is_train).sum(0).float()  # float32 counts are exact below 2^24
    base_rate = (n_pos / is_train.sum(0)).cpu()
    weights = None
    if config.loss_weights == "inv_sqrt":
        weights = n_pos.rsqrt()
        weights = weights / weights.mean()
    shift = sampling_shift(config.beta)
    if config.beta is not None:
        # Start at the sampled base rate, so the corrected model starts at the true one.
        base_rate = torch.sigmoid(torch.logit(base_rate.double()) - shift).float()

    model = nn.ModuleDict({
        "encoder": Encoder(graphs, config.d, config.dropout),
        "decoder": DECODERS[config.decoder](config.d, base_rate),
    }).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=config.lr)

    eval_graph = graphs.eval_graph.to(device)
    val_ids, val_mask = eval_rows(fold, VAL)
    val_rows = _on(val_ids, device)
    val_ab, val_labels = pairs[val_rows], labels[val_rows]
    val_mask = None if val_mask is None else _on(val_mask, device)

    log, best, best_epoch, best_state, stopped_early = [], -math.inf, 0, None, False
    for epoch in range(1, config.max_epochs + 1):
        t0 = time.time()
        this_epoch = graphs.epoch(rng)
        model.train()
        optimizer.zero_grad()
        loss = train_step(model, this_epoch.graph.to(device), _on(this_epoch.supervision, device),
                          pairs, labels, weights, train_mask, config.beta, config.chunk)
        optimizer.step()
        _sync(device)
        t1 = time.time()
        val_scores = score(model, eval_graph, val_ab, config.chunk, shift)
        val_auprc = metrics.macro_auprc_torch(val_scores, val_labels, val_mask)
        del val_scores
        t2 = time.time()
        # MPS memory held at the end of the epoch, not a peak (MPS has no peak counter).
        mps_mb = torch.mps.driver_allocated_memory() / 2**20 if device.type == "mps" else None
        log.append({"epoch": epoch, "train_loss": loss, "val_macro_auprc": val_auprc,
                    "train_s": t1 - t0, "val_s": t2 - t1, "mps_mb": mps_mb})
        print(f"[{config.split} {config.name}] epoch {epoch} loss {loss:.5f} "
              f"val macro AUPRC {val_auprc:.4f} ({t1 - t0:.2f} + {t2 - t1:.2f} s)", flush=True)
        if val_auprc > best:
            best, best_epoch = val_auprc, epoch
            best_state = {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
        if epoch - best_epoch >= config.patience:
            stopped_early = True
            break

    if best_state is None:
        raise RuntimeError("validation macro AUPRC was NaN in every epoch")
    model.load_state_dict(best_state)
    test, arrays = _test_metrics(config, model, eval_graph, pairs, data, fold, shift)
    train_s = [row["train_s"] for row in log]
    summary = {
        "best_epoch": best_epoch,
        "epochs": len(log),
        "stopped_early": stopped_early,
        "val_macro_auprc": best,
        "test": test,
        "seconds": time.time() - started,
        "mean_train_s_per_epoch": float(np.mean(train_s)),
        "max_mps_mb": max(row["mps_mb"] for row in log) if device.type == "mps" else None,
    }
    _write(config, runs_dir, model, log, arrays, summary, device)
    return summary


def _test_metrics(config, model, graph, pairs, data, fold, shift):
    """Summaries (JSON) and per-side-effect arrays (npz) of every test number §7 asks for."""
    has_targets = np.isin(np.arange(decagon.N_DRUGS), data.targets[:, 0])
    test, arrays = {}, {}

    def evaluate_rows(name, ids, mask=None, drugs=None, micro_only=False):
        scores = score(model, graph, pairs[_on(ids, pairs.device)], config.chunk, shift).cpu().numpy()
        labels = data.labels[ids].toarray()
        result = metrics.evaluate(scores, labels, mask)
        prefixes = ("micro_",) if micro_only else ("macro_", "micro_")
        averages = lambda r: {k: v for k, v in r.items() if k.startswith(prefixes)}
        test[name] = averages(result)
        if micro_only:
            return
        for k in ("auprc", "auroc", "ap50", "base_rate", "n_pos"):
            arrays[f"{name}_{k}"] = result[k]
        test[name]["subgroups"] = {}
        for group, rows in metrics.subgroups(data.pairs[ids], has_targets).items():
            sub = metrics.evaluate(scores[rows], labels[rows], None if mask is None else mask[rows])
            test[name]["subgroups"][group] = averages(sub)
        if config.n_resamples:
            test[name]["bootstrap_95"] = metrics.cluster_bootstrap(
                scores, labels, data.pairs[ids], mask, drugs, config.n_resamples, config.seed)
        if mask is not None:  # random split: also Decagon's 1:1 number
            uniform = metrics.uniform_negatives(labels, mask, config.seed)
            test[f"{name}_uniform_negatives"] = averages(metrics.evaluate(scores, labels, uniform))

    ids, mask = eval_rows(fold, TEST)
    if config.split == "cold_drug":
        evaluate_rows("test", ids, drugs=np.flatnonzero(fold.drug_roles == TEST))  # one new drug
        evaluate_rows("test_two_new", eval_rows(fold, TEST_TWO_NEW)[0], micro_only=True)
    else:
        evaluate_rows("test", ids, mask)
    return test, arrays


def _write(config, runs_dir, model, log, arrays, summary, device) -> None:
    out = config.run_dir(runs_dir)
    out.mkdir(parents=True, exist_ok=True)
    commit, dirty = _git_commit()
    meta = {**asdict(config), "name": config.name, "git_commit": commit, "git_dirty": dirty,
            "device": device.type, "torch": torch.__version__,
            "sampling": "full rows" if config.beta is None else f"sampled, beta={config.beta:g}"}
    (out / "config.json").write_text(json.dumps(meta, indent=2) + "\n")
    with open(out / "log.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(log[0]))
        writer.writeheader()
        writer.writerows(log)
    torch.save(model.state_dict(), out / "weights.pt")
    np.savez_compressed(out / "metrics.npz", **arrays)
    # Last: the run counts as done once this exists.
    (out / "metrics.json").write_text(json.dumps(summary, indent=2) + "\n")


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description="Train and evaluate one run.")
    for f in fields(Config):
        if f.type is bool:
            kind = lambda s: {"true": True, "false": False}[s.lower()]
        else:
            kind = f.type if f.type in (int, float, str) else float  # beta: float | None
        if f.default is MISSING:
            parser.add_argument(f"--{f.name}", type=kind, required=True)
        else:
            parser.add_argument(f"--{f.name}", type=kind, default=f.default)
    summary = run(Config(**vars(parser.parse_args(argv))))
    print(json.dumps({k: v for k, v in summary.items() if k != "test"}, indent=2))
    print(json.dumps(summary["test"]["test"], indent=2))


if __name__ == "__main__":
    main()
