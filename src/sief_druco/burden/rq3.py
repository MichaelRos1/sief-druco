"""The RQ3 run: burden rankings of pairs, triples and side effects on NHANES.

Spec: `docs/spec-burden.md` §8 and §10 (with §2 for the score table, §7 for prevalence).

`python -m sief_druco.burden.rq3 <hand-off folder> [cycles ...]` (default cycles: the
main pool `I P_`) reads the score table, rebuilds the exposure tables as the RQ2 run
does, and writes to `local/data/burden/rq3/<hand-off name>/<cycles>/`:

- **pairs**: every covered pair with its prevalence, Σ_r `p_r` and burden, each with
  its interval; `rank` for the pairs with ≥ 30 participants.
- **triples**: every side effect of the ranked pairs, ranked by burden.
- **side_effects**: every side effect's burden over all covered pairs, with the
  model's held-out AUPRC and predicted vs observed positives per regime.
- **novelty_shares**, **novelty_top**, **novelty_candidates**, **novelty_notes**
  (§10): the 2×2 of known to Decagon × clinically documented (DDInter) as shares of
  total burden, as counts for the top 50 next to the base rate, the novelty candidate
  list, and the footnote numbers. Pairs and triples carry both labels.

Numbers from a hand-off other than the final one (#64) are no results.
"""

import sys
from pathlib import Path

import pandas as pd

from sief_druco.burden import lookup, nhanes, rq2
from sief_druco.burden import novelty as nv
from sief_druco.burden import prevalence as prv
from sief_druco.burden import rankings as rk
from sief_druco.burden import scores as sc
from sief_druco.burden import survey

OUT_DIR = rq2.OUT_DIR.parent / "rq3"
TOP = 20


def run(
    handoff: Path,
    cycles=nhanes.MAIN_POOL,
    out_dir: Path | None = None,
    ddinter_dir: Path = nv.DDINTER_DIR,
    **exposure_args,
) -> dict[str, pd.DataFrame]:
    handoff = Path(handoff)
    table = sc.read(handoff)
    e = rq2.exposure(cycles, **exposure_args)
    names = rq2.unit_names(e.table, e.person_drugs)
    design = survey.Design.from_persons(e.persons)

    clinical = nv.ddinter(
        e.table.stitch_of, ddinter_dir,
        exposure_args.get("cache_dir", lookup.CACHE_DIR), exposure_args.get("fetch", lookup.http_fetch),
    )
    pairs = rk.pair_burden(prv.pair_prevalence(e.persons, e.pairs, design=design), table)
    pairs = nv.label_pairs(pairs, table, clinical)
    triples = nv.label_triples(rk.triple_burden(pairs, table), table, clinical)
    out = {
        "pairs": rq2.named(pairs, names),
        "triples": rq2.named(triples, names),
        "side_effects": rk.side_effect_burden(e.persons, e.pairs, table, design),
        "novelty_shares": nv.burden_shares(pairs, table),
        "novelty_top": nv.top_k(pairs, triples, table),
        "novelty_candidates": rq2.named(nv.candidates(pairs), names),
        "novelty_notes": nv.notes(pairs, clinical),
    }
    out_dir = Path(out_dir or OUT_DIR / handoff.name / "_".join(c.rstrip("_") for c in e.cycles))
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in out.items():
        df.to_csv(out_dir / f"{name}.csv", index=False)
    return out


def main(argv) -> None:
    if not argv:
        sys.exit("usage: python -m sief_druco.burden.rq3 <hand-off folder> [cycles ...]")
    out = run(Path(argv[0]), tuple(argv[1:] or nhanes.MAIN_POOL))
    shown = {
        "pairs": ["name_a", "name_b", "participants", "prevalence", "sum_p", "burden", "burden_low", "burden_high", "rank"],
        "triples": ["name_a", "name_b", "side_effect", "participants", "p", "burden", "rank"],
        "side_effects": ["name", "burden", "ci_low", "ci_high", "auprc_U", "auprc_S", "rank"],
        "novelty_candidates": ["name_a", "name_b", "participants", "prevalence", "burden", "rank", "ddinter_level", "imprecise"],
    }
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.precision", 3):
        for name, cols in shown.items():
            df = out[name]
            if name == "pairs":
                print(f"\n{len(df)} covered pairs, {df['rank'].notna().sum()} ranked; "
                      f"total burden {df.burden.sum():.1f} per 100,000")
                df = df[df["rank"].notna()]
            print(f"\n== {name}, top {TOP}\n{df[cols].head(TOP).to_string(index=False)}")
        for name in ["novelty_shares", "novelty_top", "novelty_notes"]:
            print(f"\n== {name}\n{out[name].to_string(index=False)}")


if __name__ == "__main__":
    main(sys.argv[1:])
