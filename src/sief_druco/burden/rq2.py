"""The RQ2 run: coverage and pair prevalence on NHANES. No scores needed.

Spec: `docs/spec-burden.md` §6, §7 and §3 (order of work, supplement).

`python -m sief_druco.burden.rq2 [cycles ...]` (default: the main pool `I P_`)
writes its tables to `local/data/burden/rq2/<cycles>/` and prints them; the
aggregate numbers go into `docs/results-burden.md` by hand.

- **coverage**, **failures**: the main table and the supplement of §6.
- **moved**, **coverage_no_overrides**, **thresholds_no_overrides**: what the hand
  fixes and aliases moved, against the same run without `data/decagon_overrides.csv`.
- **prevalence**: every distinct co-prescribed pair, per 100,000 people, with its
  Korn–Graubard interval (person-level data stays in `local/`).
- **thresholds**: covered pairs passing 10, 20 and 30 participants.
- **top_all**: the 20 most prevalent pairs among all co-prescribed pairs, with status
  and, when the pool has several cycles, each cycle's own share next to the pooled one.
"""

import sys
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from sief_druco.burden import coverage as cov
from sief_druco.burden import decagon_table as dt
from sief_druco.burden import handfix, lookup, nhanes
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden import prevalence as prv
from sief_druco.burden import survey
from sief_druco.burden.lookup import Cache

OUT_DIR = lookup.LOCAL_DIR / "data" / "burden" / "rq2"
TOP = 20


@dataclass
class Exposure:
    """One pool's exposure tables (spec §5), Decagon-side failures marked."""

    cycles: tuple[str, ...]
    table: dt.DrugTable
    persons: pd.DataFrame
    person_drugs: pd.DataFrame
    pairs: pd.DataFrame


def decagon_failures(
    table: dt.DrugTable, person_drugs: pd.DataFrame, cache_dir: Path, fetch: lookup.Fetch
) -> set[str]:
    """Source `IN`s whose name exactly matches a Decagon drug that reached no `IN`.

    The name check of §4 (`handfix.name_check`), restricted to the drugs still
    unmatched and to exact hits on `IN` names: the loose hits and the hits on
    matched drugs were all ruled on in #55 (`research/decagon-hand-fixes.md` §4).
    """
    if table.unmatched.empty:
        return set()
    pc = Cache(Path(cache_dir) / dt.PUBCHEM_CACHE, lookup.PUBCHEM_INTERVAL, fetch)
    rx = Cache(Path(cache_dir) / dt.RXNAV_CACHE, lookup.RXNAV_INTERVAL, fetch)
    try:
        names = handfix.decagon_names(table, handfix.sider_names(), pc, rx)
    finally:
        pc.save()
        rx.save()
    names = names[names.stitch_id.isin(set(table.unmatched.stitch_id))]
    hits = handfix.name_check(handfix.source_units(person_drugs), names)
    hits = hits[(hits.check == "exact") & (hits.kind == "in")]
    return set(hits.rxcui.dropna())


def exposure(
    cycles,
    overrides_path: Path | None = dt.OVERRIDES_PATH,
    raw_dir: Path = nhanes.RAW_DIR,
    cache_dir: Path = lookup.CACHE_DIR,
    fetch: lookup.Fetch = lookup.http_fetch,
) -> Exposure:
    cycles = tuple(cycles)
    table = dt.build(overrides_path=overrides_path, cache_dir=cache_dir, fetch=fetch)
    drugs = pdr.build_nhanes(cycles, table.stitch_of, raw_dir, cache_dir, fetch)
    drugs = cov.mark_decagon_failures(drugs, decagon_failures(table, drugs, cache_dir, fetch))
    persons = nhanes.load_persons(cycles, raw_dir)
    return Exposure(cycles, table, persons, drugs, pdr.pairs(drugs))


def unit_names(table: dt.DrugTable, person_drugs: pd.DataFrame) -> dict[str, str]:
    """Display name of each pair unit: a STITCH ID by its main `IN`, else the
    source's `IN` name or drug name."""
    main = table.ingredients[table.ingredients.role == "main"]
    names = dict(zip(main.stitch_id, main.name))
    other = person_drugs[person_drugs.status != pdr.COVERED]
    label = other.in_name.where(other.in_name.notna(), other.name)
    names.update(dict(zip(other.drug, label.fillna(""))))
    return names


def named(df: pd.DataFrame, names: dict[str, str]) -> pd.DataFrame:
    return df.assign(name_a=df.unit_a.map(names), name_b=df.unit_b.map(names))


def per_cycle(e: Exposure, top: pd.DataFrame) -> pd.DataFrame:
    """Each cycle's own share (a domain) of the `top` pairs, next to the pooled one."""
    out = top.copy()
    if len(e.cycles) < 2:
        return out
    design = survey.Design.from_persons(e.persons)
    keys = set(zip(top.unit_a, top.unit_b))
    sub = e.pairs[[k in keys for k in zip(e.pairs.unit_a, e.pairs.unit_b)]]
    for c in e.cycles:
        p = prv.pair_prevalence(e.persons, sub, domain=(e.persons.cycle == c).to_numpy(), design=design)
        p = p.set_index(["unit_a", "unit_b"])
        idx = pd.MultiIndex.from_frame(out[["unit_a", "unit_b"]])
        out[f"prevalence_{c}"] = p.prevalence.reindex(idx).to_numpy()
        out[f"participants_{c}"] = p.participants.reindex(idx).to_numpy()
    return out


def run(
    cycles=nhanes.MAIN_POOL,
    out_dir: Path | None = None,
    overrides_path: Path = dt.OVERRIDES_PATH,
    raw_dir: Path = nhanes.RAW_DIR,
    cache_dir: Path = lookup.CACHE_DIR,
    fetch: lookup.Fetch = lookup.http_fetch,
) -> dict[str, pd.DataFrame]:
    e = exposure(cycles, overrides_path, raw_dir, cache_dir, fetch)
    before = exposure(cycles, None, raw_dir, cache_dir, fetch)  # the same run without the overrides
    aliases = e.table.ingredients.rxcui[e.table.ingredients.role == "alias"]
    moved = cov.moved(before.person_drugs, e.person_drugs, aliases)
    names = unit_names(e.table, e.person_drugs)

    prev = prv.pair_prevalence(e.persons, e.pairs)
    ranked = prv.ranked(prev)
    counts = cov.counts(e.persons, e.person_drugs) | {
        "ranked_pairs": len(ranked), "population": round(e.persons.weight.sum()),
    }
    out = {
        "coverage": cov.table(e.persons, e.person_drugs, e.pairs),
        "coverage_no_overrides": cov.table(before.persons, before.person_drugs, before.pairs),
        "failures": cov.failures(e.persons, e.person_drugs),
        "moved": pd.DataFrame(
            [[b, a, m] for (b, a), m in moved["changed"].items()]
            + [["", "via_alias", moved["via_alias"]], ["", "all_mentions", moved["mentions"]]],
            columns=["status_before", "status_after", "mentions"],
        ),
        "counts": pd.DataFrame(list(counts.items()), columns=["what", "n"]),
        "prevalence": named(prev, names),
        "thresholds": prv.thresholds(prev),
        "thresholds_no_overrides": prv.thresholds(prv.pair_prevalence(before.persons, before.pairs)),
        "top_all": per_cycle(e, named(prev.head(TOP), names)),
        "collisions": e.table.collisions,
    }
    out_dir = Path(out_dir or OUT_DIR / "_".join(c.rstrip("_") for c in e.cycles))
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in out.items():
        df.to_csv(out_dir / f"{name}.csv", index=False)
    return out


def main(argv=None) -> None:
    cycles = tuple(argv if argv else nhanes.MAIN_POOL)
    out = run(cycles)
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.precision", 2):
        for name, df in out.items():
            if name != "prevalence":
                print(f"\n== {name}\n{df.to_string(index=False)}")


if __name__ == "__main__":
    main(sys.argv[1:])
