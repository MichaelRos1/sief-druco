"""NHANES persons: every interviewed participant, with pooled weights and the design.

Spec: `docs/spec-burden.md` §3 (cycles, weights, design, `RXDUSE` 7 / 9) and §5
(`persons`). Facts about the files: `research/nhanes-prescription-files.md`.

Cycles are named by their file suffix: "D" is 2005–2006 (`DEMO_D`), "I" 2015–2016,
"J" 2017–2018 and "P_" the 2017–March 2020 pre-pandemic file (`P_DEMO`).
"""

from pathlib import Path

import pandas as pd

LOCAL_DIR = Path(__file__).resolve().parents[3] / "local"
RAW_DIR = LOCAL_DIR / "data" / "raw" / "nhanes"

MAIN_POOL = ("I", "P_")
FAERS_ERA_POOL = ("D", "E")

# Years each cycle covers. 1999–2002 are left out: they need the 4-year weights
# and have retrofitted Multum codes; 2021–2023 has no drug names (spec §3).
YEARS = {c: 2.0 for c in "CDEFGHIJ"} | {"P_": 3.2}

# RXDUSE: any prescription drug in the past 30 days. One value per participant,
# repeated on each of their drug rows.
RX_USE_YES, RX_USE_NO, RX_USE_REFUSED, RX_USE_DONT_KNOW = 1, 2, 7, 9

SEX = {1: "male", 2: "female"}


def weight_factors(cycles) -> dict[str, float]:
    """Factor on each cycle's interview weight: its share of the pool's years.

    Main pool: 2/5.2 for "I", 3.2/5.2 for "P_". FAERS-era pool: 1/2 each. A single
    cycle keeps its weights. NCHS's rule for pooling cycles from 2001–2002 on.
    """
    cycles = list(cycles)
    if not cycles:
        raise ValueError("no cycles given")
    unknown = [c for c in cycles if c not in YEARS]
    if unknown:
        raise ValueError(f"unknown or unsupported cycles {unknown}; use {sorted(YEARS)}")
    if len(set(cycles)) != len(cycles):
        raise ValueError(f"a cycle is listed twice: {cycles}")
    if "J" in cycles and "P_" in cycles:
        raise ValueError(
            "P_ already holds the 2017-2018 participants under new SEQNs; "
            "pool J or P_, never both"
        )
    total = sum(YEARS[c] for c in cycles)
    return {c: YEARS[c] / total for c in cycles}


def _read(name: str, raw_dir: Path, columns: list[str]) -> pd.DataFrame:
    return pd.read_sas(raw_dir / f"{name}.xpt", format="xport", encoding="latin-1")[columns]


def _file(kind: str, cycle: str) -> str:
    return f"P_{kind}" if cycle == "P_" else f"{kind}_{cycle}"


def load_persons(cycles, raw_dir: Path = RAW_DIR) -> pd.DataFrame:
    """The `persons` table (spec §5) for the given cycles.

    One row per `DEMO` participant, non-users included: `source`, `cycle`, `seqn`,
    pooled `weight`, `stratum`, `psu`, `age` (`RIDAGEYR`, top-coded at 80 from
    2007–2008 and 85 before), `sex`, and `rx_use` (`RXDUSE`). People with
    `rx_use` 7 (refused) or 9 (don't know) stay in, with no drugs.
    """
    factors = weight_factors(cycles)
    tables = []
    for cycle, factor in factors.items():
        weight_col = "WTINTPRP" if cycle == "P_" else "WTINT2YR"
        demo = _read(
            _file("DEMO", cycle), raw_dir,
            ["SEQN", weight_col, "SDMVSTRA", "SDMVPSU", "RIDAGEYR", "RIAGENDR"],
        )
        rx = _read(_file("RXQ_RX", cycle), raw_dir, ["SEQN", "RXDUSE"])

        if rx.groupby("SEQN").RXDUSE.nunique().max() > 1:
            raise ValueError(f"cycle {cycle}: RXDUSE differs between one person's rows")
        rx = rx.drop_duplicates("SEQN")
        if set(rx.SEQN) != set(demo.SEQN):
            raise ValueError(f"cycle {cycle}: DEMO and RXQ_RX hold different participants")

        t = demo.merge(rx, on="SEQN", validate="one_to_one")
        tables.append(
            pd.DataFrame({
                "source": "NHANES",
                "cycle": cycle,
                "seqn": t.SEQN.astype("int64"),
                "weight": t[weight_col] * factor,
                "stratum": t.SDMVSTRA.astype("int64"),
                "psu": t.SDMVPSU.astype("int64"),
                "age": t.RIDAGEYR.astype("int64"),
                "sex": t.RIAGENDR.map(SEX),
                "rx_use": t.RXDUSE.astype("int64"),
            })
        )

    persons = pd.concat(tables, ignore_index=True)
    strata = persons.groupby("stratum").cycle.nunique()
    if (strata > 1).any():
        raise ValueError(f"strata shared between cycles: {strata[strata > 1].index.tolist()}")
    if persons[["weight", "sex"]].isna().any().any():
        raise ValueError("missing weight or sex")
    return persons
