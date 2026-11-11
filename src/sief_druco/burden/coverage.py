"""Coverage of real co-medication by the 645 drugs.

Spec: `docs/spec-burden.md` §6 (the tables) and §4 (the two failure reasons).

The main table has three columns adding to 100%, *covered* / *not in the model* /
*mapping failed*, and four rows, each unweighted and weighted:

- **ingredients:** distinct drugs (`drug` keys of `person_drugs`). NHANES's
  unknown / refused / don't know rows are no identifiable drug and are left out
  here, as in the pilot. Weighted, a drug counts by the pooled weight of the
  people taking it.
- **mentions:** source records × ingredients (a combination counts once per
  ingredient), weighted by the person's pooled weight.
- **participants** with ≥ 2 drugs (distinct `drug` keys, so two aliases of one
  Decagon drug are two drugs): *covered* if all their drugs are, *mapping failed*
  if any failed, else *not in the model*.
- **person-pairs:** `person_drugs.pairs`, with their status.

Coverage is descriptive: point estimates, no intervals.

A source drug can also fail on the **Decagon side**: it reached an `IN` that is a
Decagon drug which itself never reached RxNorm, so the join can't see it.
`mark_decagon_failures` turns such rows into *mapping failed* with the cause
`DECAGON_UNMATCHED`; the supplement splits failures by side and cause.
"""

import pandas as pd

from sief_druco.burden.person_drugs import (
    COVERED, MAPPING_FAILED, NOT_IN_MODEL, PERSON, STATUSES, UNCODED_PREFIX,
)

DECAGON_UNMATCHED = "decagon_unmatched"
SOURCE, DECAGON = "source", "decagon"

INGREDIENTS, MENTIONS, PARTICIPANTS, PERSON_PAIRS = "ingredients", "mentions", "participants", "person_pairs"
LEVELS = (INGREDIENTS, MENTIONS, PARTICIPANTS, PERSON_PAIRS)
TABLE_COLUMNS = ["level", "weighted", "total", *STATUSES]


def mark_decagon_failures(person_drugs: pd.DataFrame, rxcuis) -> pd.DataFrame:
    """Rows whose `IN` is in `rxcuis` (Decagon drugs that never reached RxNorm)
    become *mapping failed*, cause `DECAGON_UNMATCHED`."""
    t = person_drugs.copy()
    hit = t.rxcui.isin(set(rxcuis)) & (t.status == NOT_IN_MODEL)
    t.loc[hit, "status"] = MAPPING_FAILED
    t.loc[hit, "cause"] = DECAGON_UNMATCHED
    return t


def person_weights(persons: pd.DataFrame, table: pd.DataFrame) -> pd.Series:
    """Pooled weight of each row's person, aligned to `table`."""
    w = persons.set_index(PERSON).weight
    keys = pd.MultiIndex.from_frame(table[PERSON])
    out = w.reindex(keys)
    if out.isna().any():
        raise ValueError("rows for persons missing from the persons table")
    return pd.Series(out.to_numpy(), index=table.index)


def participant_status(person_drugs: pd.DataFrame) -> pd.DataFrame:
    """Per person: `n_drugs` (distinct `drug` keys) and the status of their drug list."""
    t = person_drugs.drop_duplicates([*PERSON, "drug"])
    g = t.groupby(PERSON).status
    out = pd.DataFrame({
        "n_drugs": g.size(),
        "any_failed": g.agg(lambda s: (s == MAPPING_FAILED).any()),
        "all_covered": g.agg(lambda s: (s == COVERED).all()),
    })
    out["status"] = NOT_IN_MODEL
    out.loc[out.all_covered, "status"] = COVERED
    out.loc[out.any_failed, "status"] = MAPPING_FAILED
    return out[["n_drugs", "status"]].reset_index()


def _shares(level: str, status: pd.Series, amount: pd.Series, weight: pd.Series) -> list[list]:
    rows = []
    for weighted, a in ((False, amount), (True, amount * weight)):
        by = a.groupby(status.to_numpy()).sum().reindex(STATUSES, fill_value=0)
        total = by.sum()
        rows.append([level, weighted, float(total), *(100 * by / total).tolist()])
    return rows


def table(persons: pd.DataFrame, person_drugs: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    """The main table (`TABLE_COLUMNS`): per level and weighting, the total (a count,
    or a pooled weight sum) and each status's share in %."""
    rows = []

    drugs = person_drugs.drop_duplicates([*PERSON, "drug"])
    drugs = drugs[~drugs.drug.str.startswith(UNCODED_PREFIX)]
    w = person_weights(persons, drugs)
    per_drug = pd.DataFrame({"status": drugs.status, "w": w}).groupby(drugs.drug.to_numpy())
    per_drug = per_drug.agg(status=("status", "first"), w=("w", "sum"))
    one = pd.Series(1.0, index=per_drug.index)
    # Weighted: each drug counts by its users' weight (`amount` 1 × the weight sum).
    rows += _shares(INGREDIENTS, per_drug.status, one, per_drug.w)

    rows += _shares(
        MENTIONS, person_drugs.status, person_drugs.mentions.astype(float),
        person_weights(persons, person_drugs),
    )

    people = participant_status(person_drugs)
    people = people[people.n_drugs >= 2].reset_index(drop=True)
    rows += _shares(
        PARTICIPANTS, people.status, pd.Series(1.0, index=people.index), person_weights(persons, people)
    )

    rows += _shares(
        PERSON_PAIRS, pairs.status, pd.Series(1.0, index=pairs.index), person_weights(persons, pairs)
    )
    return pd.DataFrame(rows, columns=TABLE_COLUMNS)


def failures(persons: pd.DataFrame, person_drugs: pd.DataFrame) -> pd.DataFrame:
    """The supplement: *mapping failed* by side and cause.

    Per (side, cause): distinct drugs (each uncoded row is its own), mentions, and
    the share of all mentions in %, unweighted and weighted.
    """
    w = person_weights(persons, person_drugs)
    m = person_drugs.mentions.astype(float)
    all_m, all_wm = m.sum(), (m * w).sum()
    t = person_drugs[person_drugs.status == MAPPING_FAILED].assign(m=m, wm=m * w)
    t = t.assign(side=(t.cause == DECAGON_UNMATCHED).map({True: DECAGON, False: SOURCE}))
    out = t.groupby(["side", "cause"]).agg(drugs=("drug", "nunique"), mentions=("m", "sum"), wm=("wm", "sum"))
    out["mentions_pct"] = 100 * out.mentions / all_m
    out["weighted_pct"] = 100 * out.wm / all_wm
    out["mentions"] = out.mentions.astype(int)
    return out.drop(columns="wm").reset_index()


def moved(before: pd.DataFrame, after: pd.DataFrame, alias_rxcuis=()) -> dict:
    """What the hand fixes and aliases moved, in mentions.

    `before` / `after`: `person_drugs` of the same source rows without and with the
    overrides. Returns `changed` (mentions per (status before, status after) that
    differ) and `via_alias` (covered mentions that join through an alias `IN`).
    """
    keys = [*PERSON, "round", "drug"]
    b = before[[*keys, "status", "mentions"]]
    a = after[[*keys, "status", "rxcui", "mentions"]]
    j = b.merge(a, on=keys, suffixes=("_before", "_after"), validate="one_to_one")
    if len(j) != len(b) or len(j) != len(a) or not (j.mentions_before == j.mentions_after).all():
        raise ValueError("before and after don't hold the same source rows")
    diff = j[j.status_before != j.status_after]
    changed = diff.groupby(["status_before", "status_after"]).mentions_after.sum().astype(int).to_dict()
    via_alias = int(after.mentions[(after.status == COVERED) & after.rxcui.isin(set(alias_rxcuis))].sum())
    return {"changed": changed, "via_alias": via_alias, "mentions": int(after.mentions.sum())}


def counts(persons: pd.DataFrame, person_drugs: pd.DataFrame, rx_unknown=(7, 9)) -> dict:
    """Head counts for the text: participants, drug users, ≥ 2 drugs, `RXDUSE` 7 / 9."""
    people = participant_status(person_drugs)
    return {
        "participants": len(persons),
        "with_drugs": len(people),
        "with_2_drugs": int((people.n_drugs >= 2).sum()),
        "rx_use_unknown": int(persons.rx_use.isin(rx_unknown).sum()),
    }

