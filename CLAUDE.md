# sief-druco

The repo is public: personal notes and learning files live in the git-ignored `local/`.

## Sources of truth

- `docs/spec-model.md` (Track 1, the model) and `docs/spec-burden.md` (Track 2, the burden): the build specs. Each section names the closed ticket holding the reasoning. Where a ticket and the spec disagree, the spec wins and the ticket gets fixed.
- `docs/adr/` for decisions, `GLOSSARY.md` for vocabulary.
- `docs/PROJECT_CONTEXT.md`: background; read the section a spec points to.
- `research/`: verified facts about external data (Decagon files, NHANES, MEPS, …). Check here before researching a source again.

## Tickets

- Issue #28 maps Track 1 and issue #51 maps Track 2: sub-issues in build order, linked by blocked-by. Each ticket names its spec sections and adds only acceptance criteria.
- One branch per ticket, `<number>-<slug>` from `main`, one PR; commit subjects end with `(#<number>)`. The user merges.

## Running

- `.venv/bin/pytest`. The full suite takes ~3 min, mostly the short training runs and the short sweep (`short_run` in their names); `-k "not short_run"` is the fast pass.
- Raw Decagon files sit in `local/data/raw/decagon/`; everything derived (splits, runs, scores) goes under `local/`.
- Training shares the one MPS GPU: run timed jobs alone, since two at once roughly double epoch times.
- Pilot scripts and logs from the design tickets live in `local/data/scratch/` (e.g. `convergence/converge.py`, the training pilot of #8).
