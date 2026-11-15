"""Novelty: how much of the burden is undocumented.

Spec: `docs/spec-burden.md` §10 (with §8 for burden and §7 for the minimum count).

Every covered pair and triple gets two labels, never merged, each `DOCUMENTED`,
`NOT_DOCUMENTED` or `NOT_CHECKABLE`:

- **known to Decagon**: a pair has any label in the 963-side-effect training table
  (`pairs.csv` `labelled`); a triple is itself a label (`labels.npz`). Every covered
  pair is in the 645, so this label is always checkable.
- **clinically documented**: DDInter 1.0 lists the pair as Major, Moderate or Minor.
  All 14 ATC files are read and de-duplicated, keeping the highest level of a pair.
  DDInter names go to RxNorm `IN`s by name (`ingredients.Resolver.name`), and an `IN`
  joins the Decagon drug it joins in the drug table (main `IN` or alias, §4). A
  DDInter "Unknown" (description incomplete) counts as not documented and is kept in
  `ddinter_level` for the footnote. A pair is not checkable when one of its drugs is
  reached by no DDInter name: the drug is not in DDInter, or its name reaches no `IN`.
  Triples inherit their pair's status: pair-level only.

A **novelty candidate** is not known to Decagon and not clinically documented (not
checkable is neither). Reporting: the 2×2 (known to Decagon × clinically documented,
with a third column for not checkable) as shares of total burden over all covered
pairs and triples (`burden_shares`); the same 2×2 as counts for the top 50 ranked
pairs and triples next to the base rate over all covered ones (`top_k`); and the
candidate list, ≥ 10 participants, by pair burden, those under 30 flagged
(`candidates`).

DDInter is CC BY-NC-SA 4.0: its files stay in `local/`; only derived labels leave it.
"""

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from sief_druco.burden import lookup
from sief_druco.burden.ingredients import Resolution, Resolver
from sief_druco.burden.prevalence import RANK_MIN
from sief_druco.burden.scores import ScoreTable

DOCUMENTED, NOT_DOCUMENTED, NOT_CHECKABLE = "documented", "not documented", "not checkable"
STATES = (DOCUMENTED, NOT_DOCUMENTED, NOT_CHECKABLE)
DECAGON_STATES = (DOCUMENTED, NOT_DOCUMENTED)

LEVELS = ("Unknown", "Minor", "Moderate", "Major")  # lowest first
GRADED = frozenset(LEVELS[1:])
UNKNOWN = "Unknown"

ATC_CODES = "ABCDGHJLMNPRSV"  # the 14 DDInter files, 8 of them not on the download page
DDINTER_DIR = lookup.LOCAL_DIR / "data" / "raw" / "ddinter"
DDINTER_RXNAV_CACHE = "rxnav_ddinter.json"

TOP_K = 50
CANDIDATE_MIN = 10
PAIR_PART = "pair-level only"  # how a triple's clinical status is marked

LABEL_COLUMNS = ["known_to_decagon", "clinically_documented", "ddinter_level", "candidate"]
SHARE_COLUMNS = ["level", "known_to_decagon", "clinically_documented", "burden", "burden_pct"]
TOP_COLUMNS = ["level", "known_to_decagon", "clinically_documented", "top_n", "top_pct", "base_n", "base_pct"]
CANDIDATE_COLUMNS = [
    "unit_a", "unit_b", "participants", "prevalence", "ci_low", "ci_high",
    "sum_p", "burden", "rank", "ddinter_level", "imprecise",
]

# resolve(name) -> Resolution
ResolveName = Callable[[str], Resolution]


def ddinter_files(raw_dir: Path = DDINTER_DIR) -> list[Path]:
    return [Path(raw_dir) / f"ddinter_downloads_code_{c}.csv" for c in ATC_CODES]


def load_ddinter(raw_dir: Path = DDINTER_DIR) -> pd.DataFrame:
    """All rows of the 14 DDInter files (`drug_a`, `drug_b`, `level`), duplicates kept."""
    files = ddinter_files(raw_dir)
    missing = [f.name for f in files if not f.exists()]
    if missing:
        raise FileNotFoundError(f"{raw_dir} lacks DDInter files {missing}")
    rows = pd.concat([pd.read_csv(f, dtype=str, keep_default_na=False) for f in files], ignore_index=True)
    rows = rows.rename(columns={"Drug_A": "drug_a", "Drug_B": "drug_b", "Level": "level"})
    bad = set(rows.level) - set(LEVELS)
    if bad:
        raise ValueError(f"DDInter levels {sorted(bad)} not in {LEVELS}")
    return rows[["drug_a", "drug_b", "level"]]


@dataclass(frozen=True)
class Clinical:
    """DDInter at the level of Decagon drugs."""

    drugs: frozenset[str]  # STITCH IDs some DDInter name reaches: the checkable drugs
    levels: pd.Series  # highest DDInter level by unordered STITCH pair (lower ID first)
    names: pd.DataFrame  # each DDInter name: rxcuis, stitch_ids (";"-joined, may be ""), cause

    def label(self, a, b) -> pd.DataFrame:
        """`clinically_documented` and `ddinter_level` ("" if not listed) of pairs (`a[i]`, `b[i]`)."""
        a, b = np.asarray(a, dtype=object), np.asarray(b, dtype=object)
        key = pd.MultiIndex.from_arrays([np.where(a < b, a, b), np.where(a < b, b, a)])
        level = self.levels.reindex(key).fillna("").to_numpy(dtype=object)
        checkable = np.array([x in self.drugs and y in self.drugs for x, y in zip(a, b)], dtype=bool)
        state = np.where(np.isin(level, list(GRADED)), DOCUMENTED, np.where(checkable, NOT_DOCUMENTED, NOT_CHECKABLE))
        return pd.DataFrame({"clinically_documented": state, "ddinter_level": level})


def map_ddinter(rows: pd.DataFrame, resolve: ResolveName, stitch_of: dict[str, str]) -> Clinical:
    """DDInter rows → pairs of Decagon drugs, highest level kept.

    `resolve` is called once per distinct name; `stitch_of` maps an `IN` RxCUI to the
    STITCH ID it joins (`decagon_table.DrugTable.stitch_of`). A name reaching several
    `IN`s (a combination) stands for each of them; a pair of one Decagon drug with
    itself is dropped.
    """
    names = sorted(set(rows.drug_a) | set(rows.drug_b))
    records = []
    for n in names:
        res = resolve(n)
        rxcuis = [r for r, _ in res.ingredients]
        stitch = sorted({stitch_of[r] for r in rxcuis if r in stitch_of})
        records.append((n, ";".join(rxcuis), ";".join(stitch), res.cause or "", stitch))
    table = pd.DataFrame(records, columns=["name", "rxcuis", "stitch_ids", "cause", "_stitch"])
    stitch_of_name = dict(zip(table.name, table._stitch))

    rank = {lv: i for i, lv in enumerate(LEVELS)}
    p = rows.assign(sa=rows.drug_a.map(stitch_of_name), sb=rows.drug_b.map(stitch_of_name), lv=rows.level.map(rank))
    p = p.explode("sa").explode("sb").dropna(subset=["sa", "sb"])
    p = p[p.sa != p.sb]
    lo, hi = np.where(p.sa < p.sb, p.sa, p.sb), np.where(p.sa < p.sb, p.sb, p.sa)
    best = pd.Series(p.lv.to_numpy(), index=pd.MultiIndex.from_arrays([lo, hi])).groupby(level=[0, 1]).max()
    levels = best.map(dict(enumerate(LEVELS)))
    drugs = frozenset(s for v in table._stitch for s in v)
    return Clinical(drugs, levels, table.drop(columns="_stitch"))


def ddinter(
    stitch_of: dict[str, str],
    raw_dir: Path = DDINTER_DIR,
    cache_dir: Path = lookup.CACHE_DIR,
    fetch: lookup.Fetch = lookup.http_fetch,
) -> Clinical:
    """DDInter mapped onto the Decagon drugs, names resolved over a disk-cached RxNav."""
    rx = lookup.Cache(Path(cache_dir) / DDINTER_RXNAV_CACHE, lookup.RXNAV_INTERVAL, fetch)
    try:
        return map_ddinter(load_ddinter(raw_dir), Resolver(rx).name, stitch_of)
    finally:
        rx.save()


def _known(flag) -> np.ndarray:
    return np.where(np.asarray(flag, dtype=bool), DOCUMENTED, NOT_DOCUMENTED)


def _candidate(df: pd.DataFrame) -> np.ndarray:
    return ((df.known_to_decagon == NOT_DOCUMENTED) & (df.clinically_documented == NOT_DOCUMENTED)).to_numpy()


def label_pairs(pairs: pd.DataFrame, table: ScoreTable, clinical: Clinical) -> pd.DataFrame:
    """Covered pairs (`unit_a`, `unit_b`) with their two labels (`LABEL_COLUMNS`)."""
    rows = table.rows(pairs.unit_a, pairs.unit_b)
    out = pairs.reset_index(drop=True).copy()
    out["known_to_decagon"] = _known(table.pairs.labelled.to_numpy()[rows])
    lab = clinical.label(out.unit_a, out.unit_b)
    out["clinically_documented"] = lab.clinically_documented.to_numpy()
    out["ddinter_level"] = lab.ddinter_level.to_numpy()
    out["candidate"] = _candidate(out)
    return out


def label_triples(triples: pd.DataFrame, table: ScoreTable, clinical: Clinical) -> pd.DataFrame:
    """Triples (`unit_a`, `unit_b`, `cui`) with their two labels; the clinical one is
    the pair's, marked in `clinical_scope`."""
    rows = table.rows(triples.unit_a, triples.unit_b)
    col = pd.Series(np.arange(len(table.side_effects)), index=table.side_effects.cui)
    cols = col.reindex(triples.cui).to_numpy()
    if np.isnan(cols.astype(float)).any():
        raise KeyError("triples name side effects not in side_effects.csv")
    out = triples.reset_index(drop=True).copy()
    out["known_to_decagon"] = _known(np.asarray(table.labels[rows, cols.astype(np.int64)]).ravel())
    lab = clinical.label(out.unit_a, out.unit_b)
    out["clinically_documented"] = lab.clinically_documented.to_numpy()
    out["ddinter_level"] = lab.ddinter_level.to_numpy()
    out["clinical_scope"] = PAIR_PART
    out["candidate"] = _candidate(out)
    return out


def _grid(df: pd.DataFrame, value: str) -> pd.Series:
    """Σ `value` by (known to Decagon, clinically documented), every cell present."""
    full = pd.MultiIndex.from_product([DECAGON_STATES, STATES], names=["known_to_decagon", "clinically_documented"])
    return df.groupby(["known_to_decagon", "clinically_documented"])[value].sum().reindex(full, fill_value=0)


def _known_sum_p(pairs: pd.DataFrame, table: ScoreTable) -> np.ndarray:
    """Σ_r p_r over each pair's labelled side effects."""
    rows = table.rows(pairs.unit_a, pairs.unit_b)
    labels = table.labels[rows]
    p = np.asarray(table.scores[rows], dtype=np.float64)
    return np.asarray(labels.multiply(p).sum(axis=1)).ravel()


def triple_parts(pairs: pd.DataFrame, table: ScoreTable) -> pd.DataFrame:
    """All covered triples, folded per pair into their labelled and unlabelled part:
    two rows per labelled pair (`burden`, `n` triples), labels as in `label_triples`."""
    rows = table.rows(pairs.unit_a, pairs.unit_b)
    k = len(table.side_effects)
    known_n = np.asarray(table.labels[rows].sum(axis=1)).ravel()
    known_b = pairs.prevalence.to_numpy() * _known_sum_p(pairs, table)
    base = pairs[["unit_a", "unit_b", "clinically_documented"]].reset_index(drop=True)
    return pd.concat([
        base.assign(known_to_decagon=DOCUMENTED, burden=known_b, n=known_n),
        base.assign(known_to_decagon=NOT_DOCUMENTED, burden=pairs.burden.to_numpy() - known_b, n=k - known_n),
    ], ignore_index=True)


def burden_shares(pairs: pd.DataFrame, table: ScoreTable) -> pd.DataFrame:
    """The headline 2×2 (`SHARE_COLUMNS`): shares (%) of total burden over all covered
    pairs and over all covered triples, no count floor. `pairs`: `label_pairs` of
    every covered pair of `rankings.pair_burden`."""
    out = []
    for level, df in [("pairs", pairs), ("triples", triple_parts(pairs, table))]:
        g = _grid(df, "burden")
        out.append(g.rename("burden").reset_index().assign(level=level, burden_pct=100 * g.to_numpy() / g.sum()))
    return pd.concat(out, ignore_index=True)[SHARE_COLUMNS]


def weight_shares(pairs: pd.DataFrame, weight: str = "prevalence") -> pd.DataFrame:
    """The pair 2×2 as shares (%) of Σ `weight`, e.g. covered pair weight (prevalence)."""
    g = _grid(pairs, weight)
    return g.rename(weight).reset_index().assign(pct=100 * g.to_numpy() / g.sum())


def top_k(pairs: pd.DataFrame, triples: pd.DataFrame, table: ScoreTable, k: int = TOP_K) -> pd.DataFrame:
    """The 2×2 as counts for the top `k` ranked pairs and triples (`TOP_COLUMNS`),
    next to the base rate over all covered pairs and triples.

    `pairs`: `label_pairs` of every covered pair; `triples`: `label_triples` of the
    ranked triples of `rankings.triple_burden`.
    """
    top_pairs = pairs[pairs["rank"].notna() & (pairs["rank"] <= k)].assign(n=1)
    top_triples = triples[triples["rank"] <= k].assign(n=1)
    sets = [
        ("pairs", top_pairs, pairs.assign(n=1)),
        ("triples", top_triples, triple_parts(pairs, table)),
    ]
    out = []
    for level, top, base in sets:
        t, b = _grid(top, "n"), _grid(base, "n")
        out.append(pd.DataFrame({
            "top_n": t.astype(int), "top_pct": 100 * t / max(t.sum(), 1),
            "base_n": b.astype(int), "base_pct": 100 * b / max(b.sum(), 1),
        }).reset_index().assign(level=level))
    return pd.concat(out, ignore_index=True)[TOP_COLUMNS]


def candidates(pairs: pd.DataFrame, min_participants: int = CANDIDATE_MIN, precise: int = RANK_MIN) -> pd.DataFrame:
    """Novelty candidate pairs with ≥ `min_participants`, highest pair burden first;
    `imprecise` under `precise` participants."""
    c = pairs[pairs.candidate & (pairs.participants >= min_participants)]
    c = c.sort_values("burden", ascending=False, kind="stable").reset_index(drop=True)
    return c.assign(imprecise=c.participants < precise)[CANDIDATE_COLUMNS]


def notes(pairs: pd.DataFrame, clinical: Clinical) -> pd.DataFrame:
    """Footnote numbers: DDInter "Unknown" and checkable shares (%) of total burden and
    of covered pair weight; DDInter names reaching an `IN` and a Decagon drug."""
    def pct(mask, col):
        return 100 * pairs[col][mask].sum() / pairs[col].sum()

    unknown = (pairs.ddinter_level == UNKNOWN).to_numpy()
    checkable = (pairs.clinically_documented != NOT_CHECKABLE).to_numpy()
    rows = [
        ("ddinter_unknown_burden_pct", pct(unknown, "burden")),
        ("ddinter_unknown_weight_pct", pct(unknown, "prevalence")),
        ("checkable_burden_pct", pct(checkable, "burden")),
        ("checkable_weight_pct", pct(checkable, "prevalence")),
        ("ddinter_names", len(clinical.names)),
        ("ddinter_names_reaching_an_in", int((clinical.names.rxcuis != "").sum())),
        ("ddinter_names_reaching_decagon", int((clinical.names.stitch_ids != "").sum())),
        ("decagon_drugs_in_ddinter", len(clinical.drugs)),
    ]
    return pd.DataFrame(rows, columns=["what", "value"])
