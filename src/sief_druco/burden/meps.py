"""MEPS persons and drug records: the light re-run of `docs/spec-burden.md` §11.2.

Facts about the files: `research/meps-prescribed-medicines.md`. The PMED files
HC-178A (2015), HC-188A (2016) and HC-197A (2017) hold one record per purchase
(fill or refill); people with no purchase are not on them.

- **Years:** 2015–2017 only. 2018–2019 are refused: the spring 2018 questionnaire
  redesign changed how purchases are collected.
- **Persons:** one row per person-year (`cycle` = year, `seqn` = `DUPERSID`), so a
  person present in two years counts twice. Weight `PERWTyyF` ÷ the number of years
  pooled (÷ 3 for the pool). The denominator is the people on the PMED files, not
  the whole population: that rescales every pair equally, and MEPS prevalence is
  used only for the top-50 overlap.
- **Point estimates only:** no HC-036 linkage and no BRR, so `stratum` and `psu` are
  empty and `survey.Design.from_persons` refuses the table. Age and sex are not on
  the PMED files and stay empty too (MEPS doesn't repeat RQ4).
- **Drugs:** `RXDRGNAM` → `IN` by name, combinations split at `-` (§4). `round` is
  `PURCHRD`, so pairs form within a purchase round (`person_drugs.pairs`). A record
  with a missing-value code instead of a name (masked, or diabetic supplies) gets the
  cause `NO_NAME`; a name that reaches no `IN` (masking by therapeutic class,
  multivitamins, ...) `NAME_NOT_FOUND`. Both are *mapping failed*.

`python -m sief_druco.burden.meps [years ...]` (default 2015 2016 2017) writes the
coverage tables of §6 and the point prevalence of every pair to
`local/data/burden/meps/<years>/`.
"""

import re
import sys
from pathlib import Path

import pandas as pd

from sief_druco.burden import coverage as cov
from sief_druco.burden import decagon_table as dt
from sief_druco.burden import lookup
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden import prevalence as prv
from sief_druco.burden import rq2
from sief_druco.burden.ingredients import Resolution, Resolver

SOURCE = "MEPS"
RAW_DIR = lookup.LOCAL_DIR / "data" / "raw" / "meps"
OUT_DIR = lookup.LOCAL_DIR / "data" / "burden" / "meps"
RXNAV_CACHE = "rxnav_meps.json"

FILES = {2015: "h178a", 2016: "h188a", 2017: "h197a"}
YEARS = tuple(FILES)
REDESIGNED = {2018: "h206a", 2019: "h213a"}  # after the spring 2018 redesign

SEPARATORS = ("-",)  # §4: MEPS writes combinations as `A-B`
NO_NAME = "no_name"  # RXDRGNAM is a missing-value code (-9 in 2015–2017)
MISSING_NAME = re.compile(r"^-\d+$")

# Records with a positive person weight by the file documentation (section 2.0). In
# the released files every record has one (research note §2.1); the released
# weights are used.
DOCUMENTED_POSITIVE_WEIGHT = {2015: 324_957, 2016: 315_212, 2017: 306_261}

RECORD_COLUMNS = ["year", "dupersid", "round", "name", "weight"]


def check_years(years) -> tuple[int, ...]:
    years = tuple(int(y) for y in years)
    if not years:
        raise ValueError("no years given")
    late = [y for y in years if y in REDESIGNED]
    if late:
        raise ValueError(
            f"MEPS {late}: the spring 2018 questionnaire redesign changed how purchases "
            "are collected; use 2015-2017 only (spec §11.2)"
        )
    unknown = [y for y in years if y not in FILES]
    if unknown:
        raise ValueError(f"unknown MEPS years {unknown}; use {list(YEARS)}")
    if len(set(years)) != len(years):
        raise ValueError(f"a year is listed twice: {years}")
    return years


def load_records(years=YEARS, raw_dir: Path = RAW_DIR) -> pd.DataFrame:
    """PMED records with a positive person weight (`RECORD_COLUMNS`).

    `weight` is the year's raw `PERWTyyF`; `persons` divides it by the years pooled.
    """
    tables = []
    for year in check_years(years):
        f = FILES[year]
        d = pd.read_sas(Path(raw_dir) / f / f"{f}.ssp", format="xport", encoding="latin-1")
        weight = d[f"PERWT{year % 100}F"]
        d = d[weight > 0]
        ids = d.DUPERSID.astype(str).str.strip()
        if not ids.str.fullmatch(r"\d+").all():
            raise ValueError(f"MEPS {year}: DUPERSID is not numeric")
        tables.append(pd.DataFrame({
            "year": year,
            "dupersid": ids.astype("int64"),
            "round": d.PURCHRD.astype("int64"),
            "name": d.RXDRGNAM.fillna("").astype(str).str.strip(),
            "weight": weight[weight > 0],
        }))
    return pd.concat(tables, ignore_index=True)[RECORD_COLUMNS]


def record_counts(years=YEARS, raw_dir: Path = RAW_DIR) -> pd.DataFrame:
    """Per year: records in the released file, those with a positive weight, and the
    documented count of positive-weight records (the footnote of §11.2)."""
    rows = []
    for year in check_years(years):
        f = FILES[year]
        d = pd.read_sas(Path(raw_dir) / f / f"{f}.ssp", format="xport", encoding="latin-1")
        positive = int((d[f"PERWT{year % 100}F"] > 0).sum())
        documented = DOCUMENTED_POSITIVE_WEIGHT[year]
        rows.append([year, len(d), positive, documented, positive - documented])
    return pd.DataFrame(rows, columns=["year", "records", "positive_weight", "documented", "difference"])


def persons(records: pd.DataFrame) -> pd.DataFrame:
    """The `persons` table (spec §5): one row per person-year on the PMED files."""
    years = check_years(sorted(records.year.unique()))
    by = records.groupby(["year", "dupersid"]).weight
    if (by.nunique() > 1).any():
        raise ValueError("a person's weight differs between their records")
    w = by.first().reset_index()
    n = len(w)
    empty = pd.array([pd.NA] * n, dtype="Int64")
    return pd.DataFrame({
        "source": SOURCE,
        "cycle": w.year.astype(str),
        "seqn": w.dupersid.astype("int64"),
        "weight": w.weight / len(years),
        "stratum": empty,
        "psu": empty.copy(),
        "age": empty.copy(),
        "sex": pd.Series([None] * n, dtype="object"),
        "rx_use": empty.copy(),
    })


def drug_rows(records: pd.DataFrame) -> pd.DataFrame:
    """One row per purchase record, as `person_drugs.build` takes them: `code` is
    empty (the PUF has no Multum drug code), `name` is `RXDRGNAM`."""
    return pd.DataFrame({
        "source": SOURCE,
        "cycle": records.year.astype(str),
        "seqn": records.dupersid.astype("int64"),
        "round": records["round"].astype("Int64"),
        "code": "",
        "name": records.name,
    })


def resolver(rx: lookup.Cache) -> pdr.Resolve:
    """`resolve(code, name)` for MEPS: a missing-value code is `NO_NAME`, a name goes
    through `Resolver.name`, split at `-`."""
    names = Resolver(rx)

    def resolve(_code: str, name: str) -> Resolution:
        if not name or MISSING_NAME.match(name):
            return Resolution((), NO_NAME, "")
        return names.name(name, SEPARATORS)

    return resolve


def build_person_drugs(
    records: pd.DataFrame,
    stitch_of: dict[str, str],
    cache_dir: Path = lookup.CACHE_DIR,
    fetch: lookup.Fetch = lookup.http_fetch,
) -> pd.DataFrame:
    rx = lookup.Cache(Path(cache_dir) / RXNAV_CACHE, lookup.RXNAV_INTERVAL, fetch)
    try:
        return pdr.build(drug_rows(records), resolver(rx), stitch_of)
    finally:
        rx.save()


def counts(records: pd.DataFrame, persons_: pd.DataFrame, person_drugs: pd.DataFrame) -> dict:
    """Head counts for the text: records reaching an `IN`, person-years with ≥ 2
    ingredients (any round, and within one round), distinct persons."""
    failed = person_drugs[~person_drugs.drug.str.startswith("in:")]
    reached = len(records) - int(failed.mentions.sum())
    ins = person_drugs[person_drugs.drug.str.startswith("in:")]
    per_year = ins.groupby(pdr.PERSON).drug.nunique()
    per_round = ins.groupby([*pdr.PERSON, "round"]).drug.nunique().groupby(level=pdr.PERSON).max()
    return {
        "records": len(records),
        "records_reaching_in": reached,
        "records_reaching_in_pct": round(100 * reached / len(records), 2),
        "no_name_records": int(failed.mentions[failed.cause == NO_NAME].sum()),
        "person_years": len(persons_),
        "persons": int(persons_.seqn.nunique()),
        "person_years_2_ins": int((per_year >= 2).sum()),
        "person_years_2_ins_same_round": int((per_round >= 2).sum()),
        "population": round(persons_.weight.sum()),
    }


def exposure(
    years=YEARS,
    overrides_path: Path | None = dt.OVERRIDES_PATH,
    raw_dir: Path = RAW_DIR,
    cache_dir: Path = lookup.CACHE_DIR,
    fetch: lookup.Fetch = lookup.http_fetch,
    records: pd.DataFrame | None = None,
) -> rq2.Exposure:
    """The MEPS exposure tables (spec §5), Decagon-side failures marked as for NHANES."""
    years = check_years(years)
    records = load_records(years, raw_dir) if records is None else records
    table = dt.build(overrides_path=overrides_path, cache_dir=cache_dir, fetch=fetch)
    drugs = build_person_drugs(records, table.stitch_of, cache_dir, fetch)
    drugs = cov.mark_decagon_failures(drugs, rq2.decagon_failures(table, drugs, cache_dir, fetch))
    return rq2.Exposure(tuple(map(str, years)), table, persons(records), drugs, pdr.pairs(drugs))


def run(
    years=YEARS,
    out_dir: Path | None = None,
    overrides_path: Path = dt.OVERRIDES_PATH,
    raw_dir: Path = RAW_DIR,
    cache_dir: Path = lookup.CACHE_DIR,
    fetch: lookup.Fetch = lookup.http_fetch,
) -> dict[str, pd.DataFrame]:
    years = check_years(years)
    records = load_records(years, raw_dir)
    e = exposure(years, overrides_path, raw_dir, cache_dir, fetch, records=records)
    names = rq2.unit_names(e.table, e.person_drugs)
    prev = prv.point_prevalence(e.persons, e.pairs)
    n = counts(records, e.persons, e.person_drugs) | {"ranked_pairs": len(prv.ranked(prev))}
    out = {
        "coverage": cov.table(e.persons, e.person_drugs, e.pairs),
        "failures": cov.failures(e.persons, e.person_drugs),
        "counts": pd.DataFrame(list(n.items()), columns=["what", "n"]),
        "record_counts": record_counts(years, raw_dir),
        "prevalence": rq2.named(prev, names),
        "thresholds": prv.thresholds(prev),
    }
    out_dir = Path(out_dir or OUT_DIR / "_".join(map(str, years)))
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in out.items():
        df.to_csv(out_dir / f"{name}.csv", index=False)
    return out


def main(argv=None) -> None:
    out = run(tuple(argv) if argv else YEARS)
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.precision", 2):
        for name, df in out.items():
            if name != "prevalence":
                print(f"\n== {name}\n{df.to_string(index=False)}")


if __name__ == "__main__":
    main(sys.argv[1:])
