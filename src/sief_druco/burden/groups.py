"""Burden by patient group (RQ4): age bands and drug-count bands.

Spec: `docs/spec-burden.md` §9 (with §8 for burden and §7 for the minimum count).

A person's burden is Σ over their covered pairs of Σ_r `p_r`: expected side-effect
signals across their pairs, a drug counted again in each pair it is part of. It is no
patient risk (aggregating a person's pairs is drug-set work, #21).

Bands are survey **domains** of the whole design, never subsets of the data:

- age (`age`, top-coded at 80): 0–17, 18–44, 45–64, 65+;
- drug count: 2, 3–4, 5–9, 10+, plus the 0–1 group, which holds no pair and appears
  only as a population share. The count is a person's distinct `drug` keys of
  `person_drugs`: every `IN`, those outside the 645 included, plus one per drug that
  never reached an `IN`. A combination pill is two ingredients, not one product.

Per band (`bands`): (a) its share of total burden, relative to everyone; (b) burden per
person per 100,000 band members with a 95% Taylor interval (`p` fixed), next to the
mean number of covered pairs per person. Per band (`band_tops`): (c) its top pairs by
band burden = band prevalence × Σ_r `p_r`, ranked only with ≥ 30 participants within
the band. Band prevalence and band burden are shares of the band's own population.
"""

import numpy as np
import pandas as pd
from scipy import stats

from sief_druco.burden import prevalence as prv
from sief_druco.burden import rankings as rk
from sief_druco.burden import survey
from sief_druco.burden.person_drugs import COVERED, PERSON
from sief_druco.burden.prevalence import PER, RANK_MIN
from sief_druco.burden.scores import ScoreTable

# Band label → lowest value; each band runs up to the next one's lowest value.
AGE_BANDS = {"0–17": 0, "18–44": 18, "45–64": 45, "65+": 65}
DRUG_BANDS = {"0–1": 0, "2": 2, "3–4": 3, "5–9": 5, "10+": 10}
NO_PAIRS = "0–1"  # the drug-count group shown only as a population share
TOP = 10

FEWER = f"fewer than {TOP} rankable pairs"
TOO_FEW = "too few participants to rank"

BAND_COLUMNS = [
    "band_type", "band", "persons", "population_pct", "burden_pct",
    "burden_per_person", "ci_low", "ci_high", "covered_pairs_per_person", "rankable_pairs", "note",
]
TOP_COLUMNS = [
    "band_type", "band", "unit_a", "unit_b", "participants", "prevalence", "ci_low", "ci_high",
    "sum_p", "burden", "band_rank", "overall_rank", "labelled",
]


def _band(values, edges: dict[str, int]) -> np.ndarray:
    lows = np.array(list(edges.values()))
    if np.any(np.asarray(values) < lows[0]):
        raise ValueError(f"values below the lowest band {lows[0]}")
    return np.array(list(edges))[np.searchsorted(lows, values, side="right") - 1]


def age_band(age) -> np.ndarray:
    return _band(age, AGE_BANDS)


def drug_band(count) -> np.ndarray:
    return _band(count, DRUG_BANDS)


def drug_count(persons: pd.DataFrame, person_drugs: pd.DataFrame) -> np.ndarray:
    """Each person's distinct `drug` keys (spec §9), 0 for persons with none."""
    n = person_drugs.drop_duplicates([*PERSON, "drug"]).groupby(PERSON).size()
    return n.reindex(pd.MultiIndex.from_frame(persons[PERSON]), fill_value=0).to_numpy()


def band_domains(persons: pd.DataFrame, person_drugs: pd.DataFrame) -> dict[tuple[str, str], np.ndarray]:
    """(band type, band) → bool [persons], age bands first, each in band order."""
    by_type = {
        "age": (age_band(persons.age.to_numpy()), AGE_BANDS),
        "drugs": (drug_band(drug_count(persons, person_drugs)), DRUG_BANDS),
    }
    return {(t, b): labels == b for t, (labels, edges) in by_type.items() for b in edges}


def _covered(persons, pairs) -> tuple:
    """Persons × covered pairs 0/1, and the covered pairs (`unit_a`, `unit_b`)."""
    y, distinct = prv.pair_matrix(persons, pairs)
    covered = (distinct.status == COVERED).to_numpy()
    return y[:, np.flatnonzero(covered)], distinct[covered].reset_index(drop=True)


def person_burden(persons: pd.DataFrame, pairs: pd.DataFrame, table: ScoreTable) -> tuple[np.ndarray, np.ndarray]:
    """Each person's burden (Σ over their covered pairs of Σ_r p_r) and covered pair count."""
    y, distinct = _covered(persons, pairs)
    rows = table.rows(distinct.unit_a, distinct.unit_b)
    sum_p = np.asarray(table.scores[rows], dtype=np.float64).sum(axis=1)
    return y @ sum_p, np.asarray(y.sum(axis=1)).astype(int)


def bands(
    persons: pd.DataFrame,
    person_drugs: pd.DataFrame,
    pairs: pd.DataFrame,
    table: ScoreTable,
    design=None,
    min_participants: int = RANK_MIN,
    top: int = TOP,
    level: float = 0.95,
) -> pd.DataFrame:
    """(a) and (b) of §9 for every band (`BAND_COLUMNS`), with each band's count of
    rankable pairs and the note its top list gets."""
    design = design or survey.Design.from_persons(persons)
    burden, n_pairs = person_burden(persons, pairs, table)
    y, _ = _covered(persons, pairs)  # for each pair's participants within a band
    w = design.weight
    total = w @ burden
    rows = []
    for (band_type, band), dom in band_domains(persons, person_drugs).items():
        row = dict(
            band_type=band_type, band=band, persons=int(dom.sum()),
            population_pct=100 * w[dom].sum() / w.sum(),
            burden_pct=100 * (w * burden)[dom].sum() / total,
        )
        if (band_type, band) == ("drugs", NO_PAIRS):
            rows.append(row | dict(note=""))
            continue
        est = survey.mean(design, burden, domain=dom)
        half = stats.t.ppf(0.5 + level / 2, est.df) * est.se
        rankable = int((np.asarray(y.T @ dom.astype(float)) >= min_participants).sum())
        rows.append(row | dict(
            burden_per_person=PER * est.value,
            ci_low=PER * (est.value - half),
            ci_high=PER * (est.value + half),
            covered_pairs_per_person=survey.mean(design, n_pairs, domain=dom).value,
            rankable_pairs=rankable,
            note=TOO_FEW if rankable == 0 else FEWER if rankable < top else "",
        ))
    out = pd.DataFrame(rows, columns=BAND_COLUMNS)
    out["rankable_pairs"] = out.rankable_pairs.astype("Int64")
    return out


def band_tops(
    persons: pd.DataFrame,
    person_drugs: pd.DataFrame,
    pairs: pd.DataFrame,
    table: ScoreTable,
    design=None,
    min_participants: int = RANK_MIN,
    top: int = TOP,
) -> pd.DataFrame:
    """(c) of §9: each band's top pairs by band burden (`TOP_COLUMNS`), ranked within
    the band, with the overall rank of §8 and "known to Decagon" (`labelled`) next to
    them. Prevalence and burden are per 100,000 band members."""
    design = design or survey.Design.from_persons(persons)
    covered = pairs[pairs.status == COVERED]
    overall = rk.pair_burden(prv.pair_prevalence(persons, covered, design=design), table, min_participants)
    overall_rank = overall.set_index(["unit_a", "unit_b"])["rank"]
    tables = []
    for (band_type, band), dom in band_domains(persons, person_drugs).items():
        if (band_type, band) == ("drugs", NO_PAIRS):
            continue
        in_band = covered[_in(covered, persons[dom])]
        if in_band.empty:
            continue
        prev = prv.pair_prevalence(persons, in_band, domain=dom, design=design)
        b = rk.pair_burden(prev, table, min_participants)
        b = b[b["rank"].notna()].head(top)
        idx = pd.MultiIndex.from_frame(b[["unit_a", "unit_b"]])
        tables.append(b.assign(
            band_type=band_type, band=band, band_rank=b["rank"],
            overall_rank=overall_rank.reindex(idx).to_numpy(),
        ))
    if not tables:
        return pd.DataFrame(columns=TOP_COLUMNS)
    out = pd.concat(tables, ignore_index=True)
    out["overall_rank"] = out.overall_rank.astype("Int64")
    return out[TOP_COLUMNS]


def _in(pairs: pd.DataFrame, persons: pd.DataFrame) -> np.ndarray:
    """Person-pairs held by one of `persons`."""
    keys = pd.MultiIndex.from_frame(persons[PERSON])
    return pd.MultiIndex.from_frame(pairs[PERSON]).isin(keys)


def person_pair_shares(persons: pd.DataFrame, person_drugs: pd.DataFrame, pairs: pd.DataFrame) -> pd.Series:
    """Each band's share (%) of the unweighted covered person-pairs, as #25's pilot counted them."""
    covered = pairs[pairs.status == COVERED]
    return pd.Series({
        key: 100 * _in(covered, persons[dom]).sum() / len(covered)
        for key, dom in band_domains(persons, person_drugs).items()
    })
