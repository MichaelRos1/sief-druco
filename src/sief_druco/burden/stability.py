"""Re-runs on other exposure data and the top-50 stability measure.

Spec: `docs/spec-burden.md` §11.2 (with §8 for burden, §7 for the minimum count).

Each re-run goes through the same code as the main result (`rankings.pair_burden`)
and forms its own top 50: pairs with ≥ 30 unweighted participants (MEPS: person-years),
then the 50 highest pair burdens. It is compared with the main top 50 by **overlap**,
the shared pairs out of 50, and each main top-50 pair that drops out gets one reason:

- **not observed in the re-run:** nobody in the re-run takes the pair;
- **below 30 participants:** taken, but too rarely to rank;
- **ranked lower:** ranked, but not in the re-run's top 50.

Pairs only; no pass/fail threshold and no rank correlations. RQ4 isn't repeated.

The re-runs (`RERUNS`):

- **2015–16** and **`P_`** separately: each cycle a survey domain of the main pool, so
  its prevalence is a share of that cycle's own population (§3 supplement).
- **FAERS-era pool:** NHANES 2005–2006 + 2007–2008 (`D`, `E`), through §4 unchanged.
  Its coverage (§6) and the codes it adds to the main pool's (`new_codes`) are
  reported too. Estriol is the class code `a70248` there and stays *mapping failed*.
- **MEPS 2015–2017 (light):** purchases in one round, point estimates (`meps`).

`python -m sief_druco.burden.stability <hand-off folder>` writes to
`local/data/burden/stability/<hand-off name>/`:

- **overlap**: per re-run, its ranked pairs, the overlap and the drop-outs by reason;
- **dropouts**: each main top-50 pair per re-run, with its rank and participants
  there and the reason it dropped out (empty if it stayed);
- **faers_coverage**, **faers_failures**, **faers_counts**: the FAERS-era coverage
  table and its supplement (§6), head counts and the new codes.

Numbers from a hand-off other than the final one (#64) are no results.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from sief_druco.burden import coverage as cov
from sief_druco.burden import lookup, meps, nhanes, rq2
from sief_druco.burden import prevalence as prv
from sief_druco.burden import rankings as rk
from sief_druco.burden import scores as sc
from sief_druco.burden import survey

OUT_DIR = rq2.OUT_DIR.parent / "stability"
TOP = 50

NOT_OBSERVED = "not observed in the re-run"
BELOW_MIN = "below 30 participants"
RANKED_LOWER = "ranked lower"
REASONS = (NOT_OBSERVED, BELOW_MIN, RANKED_LOWER)

CYCLE_I, CYCLE_P, FAERS_ERA, MEPS = "2015-16", "P_", "faers_era", "meps"
RERUNS = (CYCLE_I, CYCLE_P, FAERS_ERA, MEPS)
CYCLE_OF = {CYCLE_I: "I", CYCLE_P: "P_"}

ESTRIOL_CODE = "a70248"  # estriol in 2005–06: a class code, `d05027` in the main pool

KEY = ["unit_a", "unit_b"]
DROPOUT_COLUMNS = [
    "rerun", *KEY, "main_rank", "rerun_participants", "rerun_rank", "in_top", "reason",
]
OVERLAP_COLUMNS = ["rerun", "ranked", "top", "overlap", *REASONS]


def _in_top(rank: pd.Series, k: int) -> pd.Series:
    return rank.le(k).fillna(False).astype(bool)


def compare(main: pd.DataFrame, rerun: pd.DataFrame, k: int = TOP) -> pd.DataFrame:
    """Each pair of the main top `k` against a re-run (`DROPOUT_COLUMNS` less `rerun`).

    `main`, `rerun`: tables of `rankings.pair_burden`, each ranked within its own data.
    `reason` is empty for a pair in the re-run's top `k`.
    """
    m = main[_in_top(main["rank"], k)][[*KEY, "rank"]].rename(columns={"rank": "main_rank"})
    r = rerun[[*KEY, "participants", "rank"]].rename(
        columns={"participants": "rerun_participants", "rank": "rerun_rank"}
    )
    out = m.merge(r, on=KEY, how="left", validate="one_to_one")
    out["rerun_participants"] = out.rerun_participants.fillna(0).astype(int)
    out["in_top"] = _in_top(out.rerun_rank, k).to_numpy()
    out["reason"] = np.select(
        [out.in_top, out.rerun_participants == 0, out.rerun_participants < prv.RANK_MIN],
        ["", NOT_OBSERVED, BELOW_MIN],
        RANKED_LOWER,
    )
    return out.sort_values("main_rank").reset_index(drop=True)


def overlap(rerun: pd.DataFrame, comparison: pd.DataFrame, k: int = TOP) -> dict:
    """One `OVERLAP_COLUMNS` row (less `rerun`): the re-run's ranked pairs, the size of
    its top `k`, the shared pairs and the drop-outs by reason."""
    ranked = int(rerun["rank"].notna().sum())
    reasons = comparison.reason.value_counts()
    return {
        "ranked": ranked, "top": min(ranked, k), "overlap": int(comparison.in_top.sum()),
        **{r: int(reasons.get(r, 0)) for r in REASONS},
    }


def stability(
    main: pd.DataFrame, reruns: dict[str, pd.DataFrame], k: int = TOP
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The overlap table and the drop-out table of every re-run against `main`."""
    rows, dropouts = [], []
    for name, rerun in reruns.items():
        c = compare(main, rerun, k)
        rows.append({"rerun": name, **overlap(rerun, c, k)})
        dropouts.append(c.assign(rerun=name))
    return (
        pd.DataFrame(rows, columns=OVERLAP_COLUMNS),
        pd.concat(dropouts, ignore_index=True)[DROPOUT_COLUMNS],
    )


def new_codes(rows: pd.DataFrame, main_rows: pd.DataFrame) -> dict:
    """Drug codes of `rows` (`nhanes.load_drug_rows`) not in the main pool's
    `main_rows`: distinct codes, how many are `a` / `c` / `h` codes, their mentions
    and the share of all coded mentions in %. Records with no code are left out."""
    coded = rows[rows.code != ""]
    new = coded[~coded.code.isin(set(main_rows.code))]
    return {
        "new_codes": int(new.code.nunique()),
        "new_codes_ach": int(new.code[new.code.str[0].isin(list("ach"))].nunique()),
        "new_code_mentions": len(new),
        "coded_mentions": len(coded),
        "new_code_mentions_pct": round(100 * len(new) / len(coded), 2),
    }


def code_status(person_drugs: pd.DataFrame, code: str) -> pd.DataFrame:
    """Status and cause of the rows from one source drug code."""
    return person_drugs.loc[person_drugs.code == code, ["status", "cause", "name", "mentions"]]


def burden(
    exposure: rq2.Exposure, table: sc.ScoreTable, point: bool = False, domain=None, design=None
) -> pd.DataFrame:
    """`rankings.pair_burden` of one pool, ranked within it: survey prevalence (of a
    `domain`'s own population if given), or point prevalence for MEPS."""
    if point:
        prev = prv.point_prevalence(exposure.persons, exposure.pairs)
    else:
        prev = prv.pair_prevalence(exposure.persons, exposure.pairs, domain=domain, design=design)
    return rk.pair_burden(prev, table)


def run(
    handoff: Path,
    out_dir: Path | None = None,
    nhanes_raw_dir: Path = nhanes.RAW_DIR,
    meps_raw_dir: Path = meps.RAW_DIR,
    cache_dir: Path = lookup.CACHE_DIR,
    fetch: lookup.Fetch = lookup.http_fetch,
) -> dict[str, pd.DataFrame]:
    handoff = Path(handoff)
    table = sc.read(handoff)
    look = {"cache_dir": cache_dir, "fetch": fetch}

    pool = rq2.exposure(nhanes.MAIN_POOL, raw_dir=nhanes_raw_dir, **look)
    design = survey.Design.from_persons(pool.persons)
    main = burden(pool, table, design=design)
    reruns = {
        name: burden(pool, table, domain=(pool.persons.cycle == cycle).to_numpy(), design=design)
        for name, cycle in CYCLE_OF.items()
    }

    faers = rq2.exposure(nhanes.FAERS_ERA_POOL, raw_dir=nhanes_raw_dir, **look)
    reruns[FAERS_ERA] = burden(faers, table)
    purchases = meps.exposure(meps.load_records(meps.YEARS, meps_raw_dir), **look)
    reruns[MEPS] = burden(purchases, table, point=True)

    # A covered unit is a STITCH ID, named by its main `IN` in every pool.
    names = {}
    for x in (purchases, faers, pool):
        names |= rq2.unit_names(x.table, x.person_drugs)
    overlaps, dropouts = stability(main, reruns)

    codes = new_codes(
        nhanes.load_drug_rows(nhanes.FAERS_ERA_POOL, nhanes_raw_dir),
        nhanes.load_drug_rows(nhanes.MAIN_POOL, nhanes_raw_dir),
    )
    estriol = code_status(faers.person_drugs, ESTRIOL_CODE)
    if estriol.empty:
        raise ValueError(f"no FAERS-era row has the estriol code {ESTRIOL_CODE}")
    counts = cov.counts(faers.persons, faers.person_drugs) | codes | {
        "estriol_a70248_status": ";".join(sorted(set(estriol.status))),
        "population": round(faers.persons.weight.sum()),
    }
    out = {
        "overlap": overlaps,
        "dropouts": rq2.named(dropouts, names),
        "faers_coverage": cov.table(faers.persons, faers.person_drugs, faers.pairs),
        "faers_failures": cov.failures(faers.persons, faers.person_drugs),
        "faers_counts": pd.DataFrame(list(counts.items()), columns=["what", "n"]),
    }
    out_dir = Path(out_dir or OUT_DIR / handoff.name)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in out.items():
        df.to_csv(out_dir / f"{name}.csv", index=False)
    return out


def main(argv) -> None:
    if len(argv) != 1:
        sys.exit("usage: python -m sief_druco.burden.stability <hand-off folder>")
    out = run(Path(argv[0]))
    cols = ["main_rank", "name_a", "name_b", "rerun_participants", "rerun_rank", "reason"]
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.precision", 2):
        print(f"\n== overlap of the top {TOP}\n{out['overlap'].to_string(index=False)}")
        d = out["dropouts"]
        for name, g in d[~d.in_top].groupby("rerun", sort=False):
            print(f"\n== dropped out, {name}\n{g[cols].to_string(index=False)}")
        for name in ["faers_coverage", "faers_failures", "faers_counts"]:
            print(f"\n== {name}\n{out[name].to_string(index=False)}")


if __name__ == "__main__":
    main(sys.argv[1:])
