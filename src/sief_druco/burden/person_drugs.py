"""The `person_drugs` table and the drug pairs derived from it.

Spec: `docs/spec-burden.md` §4 (step 3, the two failure reasons) and §5.

`build` takes one row per source drug record and returns one row per person ×
drug × `round`:

- a drug that reached RxNorm ingredients gives one row per `IN` (`drug` = `in:<rxcui>`),
  with its `stitch_id` (*covered*) or *not in the model*;
- a drug that reached no `IN` keeps one row with its `cause` (*mapping failed*),
  because §9 counts it as a drug: one row per failed code (`code:<code>`), per
  failed name when the source has no code (`name:<name>`, MEPS and DDInter), and
  per record for NHANES's unknown / refused / don't know (`row:<n>`, as the pilot).

`mentions` counts the source records behind each row, so coverage (§6) can count
mentions after the de-duplication. Nothing is dropped: the table is each person's
full drug list (#21).

`pairs` derives the person-pairs on the fly; they are never stored as counts.
"""

from collections.abc import Callable
from pathlib import Path

import pandas as pd

from sief_druco.burden import lookup, nhanes
from sief_druco.burden.ingredients import NO_CODE, Resolution, Resolver

COVERED, NOT_IN_MODEL, MAPPING_FAILED = "covered", "not_in_model", "mapping_failed"
UNCODED_PREFIX = "row:"  # `drug` key of an NHANES unknown / refused / don't know row
STATUSES = (COVERED, NOT_IN_MODEL, MAPPING_FAILED)

NHANES_RXNAV_CACHE = "rxnav_nhanes.json"

PERSON = ["source", "cycle", "seqn"]
ROW_COLUMNS = [*PERSON, "round", "code", "name"]
COLUMNS = [
    *PERSON, "round", "drug", "rxcui", "in_name", "stitch_id", "status", "cause",
    "code", "name", "mentions",
]
PAIR_COLUMNS = [*PERSON, "unit_a", "unit_b", "status", "cause_a", "cause_b"]

# resolve(code, name) -> Resolution
Resolve = Callable[[str, str], Resolution]


def build(drug_rows: pd.DataFrame, resolve: Resolve, stitch_of: dict[str, str]) -> pd.DataFrame:
    """`person_drugs` from source drug records (`ROW_COLUMNS`).

    `resolve` is called once per distinct (`code`, `name`); `stitch_of` maps an `IN`
    RxCUI to the STITCH ID it joins (`decagon_table.DrugTable.stitch_of`).
    """
    missing = set(ROW_COLUMNS) - set(drug_rows.columns)
    if missing:
        raise ValueError(f"drug rows lack columns {sorted(missing)}")
    records = drug_rows[ROW_COLUMNS].reset_index(drop=True)
    records["record"] = records.index

    resolved = []
    for code, name in records[["code", "name"]].drop_duplicates().itertuples(index=False):
        res = resolve(code, name)
        if res.ingredients:
            resolved += [(code, name, rxcui, in_name, None) for rxcui, in_name in res.ingredients]
        else:
            resolved.append((code, name, None, None, res.cause))
    resolved = pd.DataFrame(resolved, columns=["code", "name", "rxcui", "in_name", "cause"])
    t = records.merge(resolved, on=["code", "name"], how="left")

    reached = t.rxcui.notna()
    by_name = ("name:" + t.name).where(t.cause != NO_CODE, UNCODED_PREFIX + t.record.astype(str))
    failed_key = ("code:" + t.code).where(t.code != "", by_name)
    t["drug"] = ("in:" + t.rxcui).where(reached, failed_key)
    t["stitch_id"] = t.rxcui.map(stitch_of)
    t["status"] = MAPPING_FAILED
    t.loc[reached, "status"] = NOT_IN_MODEL
    t.loc[t.stitch_id.notna(), "status"] = COVERED

    keys = [*PERSON, "round", "drug"]
    t["mentions"] = t.groupby(keys, dropna=False).record.transform("size")
    t = t.drop_duplicates(keys).reset_index(drop=True)
    t["round"] = t["round"].astype("Int64")
    return t[COLUMNS]


def pairs(person_drugs: pd.DataFrame) -> pd.DataFrame:
    """One row per person × unordered pair of distinct drugs (spec §5).

    A covered drug's unit is its STITCH ID, so aliases of one Decagon drug never
    pair; other units are the `drug` key. Pairs form within a `round` (MEPS), and a
    pair seen in several rounds counts once. `unit_a` < `unit_b`: for a covered pair
    both are STITCH IDs, in Decagon's index order; otherwise they are `drug` keys of
    `person_drugs`. `status`: *covered* if both drugs are covered, *mapping failed*
    if either failed, else *not in the model*. `cause_a` / `cause_b`: each side's
    failure cause, empty unless that side failed.
    """
    t = person_drugs.assign(
        unit=person_drugs.stitch_id.where(person_drugs.status == COVERED, person_drugs.drug),
        round_key=person_drugs["round"].astype("Int64").fillna(-1),
    )
    units = t.drop_duplicates([*PERSON, "round_key", "unit"])
    units = units[[*PERSON, "round_key", "unit", "status", "cause"]]
    p = units.merge(units, on=[*PERSON, "round_key"], suffixes=("_a", "_b"))
    p = p[p.unit_a < p.unit_b]

    status = pd.Series(NOT_IN_MODEL, index=p.index)
    status[(p.status_a == COVERED) & (p.status_b == COVERED)] = COVERED
    status[(p.status_a == MAPPING_FAILED) | (p.status_b == MAPPING_FAILED)] = MAPPING_FAILED
    p = p.assign(status=status).drop_duplicates([*PERSON, "unit_a", "unit_b"])
    return p[PAIR_COLUMNS].reset_index(drop=True)


def build_nhanes(
    cycles,
    stitch_of: dict[str, str],
    raw_dir: Path = nhanes.RAW_DIR,
    cache_dir: Path = lookup.CACHE_DIR,
    fetch: lookup.Fetch = lookup.http_fetch,
) -> pd.DataFrame:
    """`person_drugs` for NHANES cycles, codes resolved by `Resolver.multum`."""
    rx = lookup.Cache(Path(cache_dir) / NHANES_RXNAV_CACHE, lookup.RXNAV_INTERVAL, fetch)
    try:
        return build(nhanes.load_drug_rows(cycles, raw_dir), Resolver(rx).multum, stitch_of)
    finally:
        rx.save()
