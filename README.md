# Polypharmacy Side-Effect Burden

Predicting side effects of drugs taken together, and weighting those predictions by how often people actually take them together.

## Development

Python 3.13. From the repo root:

```bash
python3.13 -m venv .venv
.venv/bin/pip install -e '.[dev]'
.venv/bin/pytest
```

The raw Decagon files go in `local/data/raw/decagon/` (git-ignored, with `SHA256SUMS`). Tests that need them skip when they are missing. The parsed dataset is cached in `local/data/processed/decagon.npz` on first load.

One training run (config → best weights, test metrics and an epoch log in `local/runs/`):

```bash
.venv/bin/python -m sief_druco.train --split pair --fold 0 --seed 0
```
