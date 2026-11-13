"""The RQ4 run: burden by age band and drug-count band on the main NHANES pool.

Spec: `docs/spec-burden.md` §9.

`python -m sief_druco.burden.rq4 <hand-off folder>` reads the score table, rebuilds
the exposure tables of the main pool as the RQ2 run does, and writes to
`local/data/burden/rq4/<hand-off name>/`:

- **bands**: each band's share of the population and of total burden, burden per
  person per 100,000 band members with its interval, mean covered pairs per person,
  and the note its top list gets.
- **tops**: each band's top 10 pairs by band burden (≥ 30 participants within the band).

RQ4 runs on the main pool only; the re-runs of §11 don't repeat it. Numbers from a
hand-off other than the final one (#64) are no results.
"""

import sys
from pathlib import Path

import pandas as pd

from sief_druco.burden import groups as gr
from sief_druco.burden import nhanes, rq2
from sief_druco.burden import scores as sc
from sief_druco.burden import survey

OUT_DIR = rq2.OUT_DIR.parent / "rq4"


def run(handoff: Path, out_dir: Path | None = None, **exposure_args) -> dict[str, pd.DataFrame]:
    handoff = Path(handoff)
    table = sc.read(handoff)
    e = rq2.exposure(nhanes.MAIN_POOL, **exposure_args)
    design = survey.Design.from_persons(e.persons)
    args = (e.persons, e.person_drugs, e.pairs, table, design)
    out = {
        "bands": gr.bands(*args),
        "tops": rq2.named(gr.band_tops(*args), rq2.unit_names(e.table, e.person_drugs)),
    }
    out_dir = Path(out_dir or OUT_DIR / handoff.name)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in out.items():
        df.to_csv(out_dir / f"{name}.csv", index=False)
    return out


def main(argv) -> None:
    if len(argv) != 1:
        sys.exit("usage: python -m sief_druco.burden.rq4 <hand-off folder>")
    out = run(Path(argv[0]))
    cols = ["band", "name_a", "name_b", "participants", "prevalence", "burden", "band_rank", "overall_rank", "labelled"]
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.precision", 2):
        print(f"\n== bands\n{out['bands'].to_string(index=False)}")
        for (band_type, band), g in out["tops"].groupby(["band_type", "band"], sort=False):
            print(f"\n== {band_type} {band}\n{g[cols].to_string(index=False)}")


if __name__ == "__main__":
    main(sys.argv[1:])
