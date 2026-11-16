"""Positive controls: does the model score known harmful pairs high, and does the
burden ranking reach them?

Spec: `docs/spec-burden.md` §11.1 (with §8 for burden and §7 for the minimum count).

The controls are the committed list `data/positive_controls.csv` (#60), written down
before any score was looked up: one row per drug or side effect of a Beers 2023 entry.
An entry's **control pairs** are every `a` × `b` drug, or any two distinct `member`
drugs; drug rows with no `id` are not in the 645 and form no pair. A control is one
(entry, pair): the same pair under two entries is two controls, as their side effects
differ.

- **Model view** (`model_view`): for each control and each of its entry's side effects
  r, the percentile of `p_r(a, b)` among all pairs' scores for r: the share of pairs,
  in %, scoring strictly below it. A tie counts against the control. **In the top
  10%** means at most a tenth of all pairs score as high or higher, which is the same
  as a percentile ≥ 90. **Recovered** = in the top 10% for at least one side effect.
  This view reads only the score table, so it doesn't depend on who takes which pair.
- **Not checkable:** a control whose entry maps to no side effect of the score table.
  It has no percentile and is neither recovered nor not recovered.
- **Burden view** (`burden_view`): the pair's overall rank from `rankings.pair_burden`
  (only pairs with ≥ 30 participants are ranked) and whether it is in the top 50.
- **Labelled in Decagon** is marked, as S-regime scores are partly expected to be high:
  `labelled` for the pair (any side effect), `labelled_mapped` for one of the entry's
  side effects.

The report (`entries`) gives "k of n recovered" per entry and over all controls, with no
overall pass/fail; `unchecked` lists what couldn't be checked and why.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

from sief_druco.burden import decagon_table as dt
from sief_druco.burden import nhanes, rq2
from sief_druco.burden import prevalence as prv
from sief_druco.burden import rankings as rk
from sief_druco.burden import scores as sc
from sief_druco.burden import survey

CONTROLS_PATH = dt.REPO_DIR / "data" / "positive_controls.csv"
OUT_DIR = rq2.OUT_DIR.parent / "controls"
TOP_PCT = 10  # recovered = in the top 10% of all pairs for a side effect
TOP = 50  # burden view: in the top 50 ranked pairs

DRUG_KINDS = ("a", "b", "member")
OUT_OF_SCOPE = "out of scope"
NOT_IN_645 = "not in the 645"
NO_PAIR = "no pair of drugs in the 645"
NO_SIDE_EFFECT = "no side effect in the score table"

PAIR_COLUMNS = ["entry", "label", "unit_a", "unit_b", "name_a", "name_b"]
TRIPLE_COLUMNS = [*PAIR_COLUMNS, "cui", "side_effect", "p", "percentile", "top", "labelled"]
CONTROL_COLUMNS = [
    *PAIR_COLUMNS, "labelled", "labelled_mapped", "checkable",
    "best_side_effect", "best_p", "best_percentile", "recovered",
]
ENTRY_COLUMNS = ["entry", "label", "controls", "checkable", "recovered", "recovered_text", "labelled", "ranked", "top50"]
UNCHECKED_COLUMNS = ["entry", "label", "name", "reason"]


def read(path: Path = CONTROLS_PATH) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False).replace("", pd.NA).astype({"entry": int})


def _in_table(controls: pd.DataFrame, table: sc.ScoreTable) -> pd.DataFrame:
    """The drug and side-effect rows whose `id` is in the score table."""
    drugs = set(table.pairs.stitch_a) | set(table.pairs.stitch_b)
    is_drug = controls.kind.isin(DRUG_KINDS) & controls.id.isin(drugs)
    is_effect = (controls.kind == "side_effect") & controls.id.isin(set(table.side_effects.cui))
    return controls[is_drug | is_effect]


def _pairs(drugs: pd.DataFrame) -> pd.DataFrame:
    """Every control (`PAIR_COLUMNS`) of the drug rows in the score table: an entry's
    `a` × `b` pairs, or any two distinct `member` drugs; unordered, lower STITCH ID first."""
    rows = []
    for entry, g in drugs.groupby("entry", sort=True):
        name = dict(zip(g.id, g.name))
        a, b, m = (g.id[g.kind == k].unique() for k in DRUG_KINDS)
        cand = [(x, y) for x in a for y in b] + [(x, y) for i, x in enumerate(m) for y in m[i + 1:]]
        seen = set()
        for x, y in cand:
            x, y = min(x, y), max(x, y)
            if x != y and (x, y) not in seen:
                seen.add((x, y))
                rows.append((entry, g.label.iloc[0], x, y, name[x], name[y]))
    return pd.DataFrame(rows, columns=PAIR_COLUMNS)


def model_view(controls: pd.DataFrame, table: sc.ScoreTable) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Each control's percentile per side effect of its entry (`TRIPLE_COLUMNS`) and per
    control (`CONTROL_COLUMNS`, its best side effect). Reads only the score table."""
    known = _in_table(controls, table)
    ctrl = _pairs(known[known.kind.isin(DRUG_KINDS)])
    effects = known[known.kind == "side_effect"].drop_duplicates(["entry", "id"])
    col_of = pd.Series(np.arange(len(table.side_effects)), index=table.side_effects.cui)
    cols = np.unique(col_of[effects.id].to_numpy())

    ctrl_rows = table.rows(ctrl.unit_a, ctrl.unit_b) if len(ctrl) else np.zeros(0, np.int64)
    triples = ctrl.assign(row=ctrl_rows).merge(
        effects[["entry", "id", "name"]].rename(columns={"id": "cui", "name": "side_effect"}), on="entry"
    )
    n = len(table.pairs)
    if len(triples):
        j = np.searchsorted(cols, col_of[triples.cui].to_numpy())
        col_scores = np.asarray(table.scores[:, cols])  # [pairs, mapped side effects], one pass
        p = col_scores[triples.row.to_numpy(), j]
        below = np.empty(len(triples), np.int64)
        for k in range(len(cols)):
            at = j == k
            below[at] = np.searchsorted(np.sort(col_scores[:, k]), p[at], side="left")
        labels = table.labels[triples.row.to_numpy(), col_of[triples.cui].to_numpy()]
    else:
        p, below, labels = np.zeros(0, np.float32), np.zeros(0, np.int64), np.zeros(0)
    triples = triples.assign(
        p=p.astype(np.float64),
        percentile=100 * below / n,
        top=(n - below) * 100 <= TOP_PCT * n,  # at most 10% of pairs score ≥ p
        labelled=np.asarray(labels).astype(bool),
    )

    best = triples.sort_values(["percentile", "p"], ascending=False, kind="stable").groupby(
        ["entry", "unit_a", "unit_b"], sort=False
    )
    best = pd.DataFrame({
        "best_side_effect": best.side_effect.first(),
        "best_p": best.p.first(),
        "best_percentile": best.percentile.first(),
        "recovered": best.top.any(),
        "labelled_mapped": best.labelled.any(),
    })
    out = ctrl.assign(labelled=table.pairs.labelled.to_numpy()[ctrl_rows].astype(bool))
    out = out.join(best, on=["entry", "unit_a", "unit_b"])
    out["checkable"] = out.best_percentile.notna()
    out["recovered"] = out.recovered.astype("boolean")  # NA when not checkable
    out["labelled_mapped"] = out.labelled_mapped.astype("boolean")
    return triples[TRIPLE_COLUMNS], out[CONTROL_COLUMNS]


def burden_view(control_df: pd.DataFrame, pair_burden: pd.DataFrame) -> pd.DataFrame:
    """The controls with their pair's participants, burden and overall rank from
    `rankings.pair_burden` (ranked only with ≥ 30 participants), and whether that rank is
    in the top 50. A pair nobody takes has 0 participants and no rank."""
    pb = pair_burden.set_index(["unit_a", "unit_b"])[["participants", "burden", "rank"]]
    out = control_df.join(pb, on=["unit_a", "unit_b"])
    out["participants"] = out.participants.fillna(0).astype(int)
    out["rank"] = out["rank"].astype("Int64")
    out["top50"] = (out["rank"] <= TOP).fillna(False).astype(bool)
    return out


def _k_of_n(g: pd.DataFrame) -> dict:
    checkable = int(g.checkable.sum())
    recovered = int(g.recovered.fillna(False).sum())
    return {
        "controls": len(g),
        "checkable": checkable,
        "recovered": recovered,
        "recovered_text": f"{recovered} of {checkable} recovered",
        "labelled": int(g.labelled.sum()),
        "ranked": int(g["rank"].notna().sum()),
        "top50": int(g.top50.sum()),
    }


def entries(controls: pd.DataFrame, report: pd.DataFrame) -> pd.DataFrame:
    """"k of n recovered" per Beers entry with a control pair, then over all controls
    (`entry` "all"), from a `burden_view` report. n counts checkable controls only. No
    pass/fail."""
    label = controls.drop_duplicates("entry").set_index("entry").label
    rows = [{"entry": e, "label": label[e], **_k_of_n(g)} for e, g in report.groupby("entry", sort=True)]
    rows.append({"entry": "all", "label": "all controls", **_k_of_n(report)})
    return pd.DataFrame(rows, columns=ENTRY_COLUMNS)


def unchecked(controls: pd.DataFrame, table: sc.ScoreTable, report: pd.DataFrame) -> pd.DataFrame:
    """What the report can't check (`UNCHECKED_COLUMNS`): out-of-scope entries, drugs not
    in the 645, entries left with no control pair, and entries with no side effect in
    the score table (their controls are not checkable)."""
    known = _in_table(controls, table)
    rows = []
    for entry, g in controls.groupby("entry", sort=True):
        label = g.label.iloc[0]
        if (g.kind == "out_of_scope").any():
            rows.append((entry, label, pd.NA, OUT_OF_SCOPE))
            continue
        drugs = g[g.kind.isin(DRUG_KINDS)]
        for name in drugs.name[~drugs.index.isin(known.index)]:
            rows.append((entry, label, name, NOT_IN_645))
        mine = report[report.entry == entry]
        if mine.empty:
            rows.append((entry, label, pd.NA, NO_PAIR))
        elif not mine.checkable.any():
            rows.append((entry, label, pd.NA, NO_SIDE_EFFECT))
    return pd.DataFrame(rows, columns=UNCHECKED_COLUMNS)


def run(
    handoff: Path, controls_path: Path = CONTROLS_PATH, out_dir: Path | None = None, **exposure_args
) -> dict[str, pd.DataFrame]:
    """The positive-control report on the main NHANES pool, written to
    `local/data/burden/controls/<hand-off name>/`."""
    handoff = Path(handoff)
    table = sc.read(handoff)
    controls = read(controls_path)
    triples, report = model_view(controls, table)

    e = rq2.exposure(nhanes.MAIN_POOL, **exposure_args)
    design = survey.Design.from_persons(e.persons)
    pb = rk.pair_burden(prv.pair_prevalence(e.persons, e.pairs, design=design), table)
    report = burden_view(report, pb)
    out = {
        "entries": entries(controls, report),
        "controls": report,
        "triples": triples,
        "unchecked": unchecked(controls, table, report),
    }
    out_dir = Path(out_dir or OUT_DIR / handoff.name)
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, df in out.items():
        df.to_csv(out_dir / f"{name}.csv", index=False)
    return out


def main(argv) -> None:
    if len(argv) != 1:
        sys.exit("usage: python -m sief_druco.burden.controls <hand-off folder>")
    out = run(Path(argv[0]))
    cols = [
        "entry", "name_a", "name_b", "labelled", "best_side_effect", "best_percentile", "recovered",
        "participants", "rank", "top50",
    ]
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.precision", 2):
        print(f"\n== entries\n{out['entries'].to_string(index=False)}")
        shown = out["controls"].sort_values(["rank", "best_percentile"], ascending=[True, False])
        print(f"\n== controls, ranked first\n{shown[cols].head(30).to_string(index=False)}")
        print(f"\n== unchecked\n{out['unchecked'].to_string(index=False)}")


if __name__ == "__main__":
    main(sys.argv[1:])
