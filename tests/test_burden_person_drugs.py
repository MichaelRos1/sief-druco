import pandas as pd
import pytest

from sief_druco.burden import ingredients as ing
from sief_druco.burden import lookup, nhanes
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden.ingredients import Resolution

# --- Synthetic drugs: code -> resolution, and IN -> STITCH ID ---

MPH, DEXMPH, METF, LISIN, ESCIT, ROSU = "6901", "1091", "6809", "29046", "321988", "301542"
STITCH_OF = {MPH: "CID000004158", DEXMPH: "CID000004158", METF: "CID000004091", LISIN: "CID000003961"}


def _found(*ins):
    return Resolution(tuple((i, f"in{i}") for i in ins), None, ing.MMSL)


def _failed(cause):
    return Resolution((), cause, "")


RESOLUTIONS = {
    "d_mph": _found(MPH),
    "d_dexmph": _found(DEXMPH),
    "d_metf": _found(METF),
    "d_metf_er": _found(METF),
    "d_lisin": _found(LISIN),
    "d_combo": _found(METF, LISIN),  # a fixed-dose combination
    "d_escit": _found(ESCIT),  # not in the model
    "d_rosu": _found(ROSU),  # not in the model
    "c00049": _failed(ing.ACH_CODE),
    "d_retired": _failed(ing.RETIRED_CODE),
    "d_unknown": _failed(ing.NAME_NOT_FOUND),
    "": _failed(ing.NO_CODE),
}


def resolve(code, name):
    return RESOLUTIONS[code]


def rows(*people):
    """`(seqn, [codes])` or `(seqn, round, [codes])` -> drug rows."""
    out = []
    for person in people:
        seqn, rnd, codes = person if len(person) == 3 else (person[0], pd.NA, person[1])
        for code in codes:
            name = "55555" if code == "" else code.upper()
            out.append(dict(source="NHANES", cycle="J", seqn=seqn, round=rnd, code=code, name=name))
    return pd.DataFrame(out)


def build(*people):
    return pdr.build(rows(*people), resolve, STITCH_OF)


def pairs_of(*people):
    return pdr.pairs(build(*people))


def pair_set(pairs):
    return {(r.seqn, r.unit_a, r.unit_b, r.status) for r in pairs.itertuples()}


# --- person_drugs ---


def test_each_row_has_a_stitch_id_or_a_reason():
    t = build((1, ["d_mph", "d_escit", "c00049"]))
    covered = t[t.status == pdr.COVERED]
    assert covered.stitch_id.notna().all() and covered.cause.isna().all()
    assert t.loc[t.status == pdr.NOT_IN_MODEL, ["stitch_id", "cause"]].isna().all().all()
    failed = t[t.status == pdr.MAPPING_FAILED]
    assert failed.stitch_id.isna().all() and failed.rxcui.isna().all()
    assert failed.cause.tolist() == [ing.ACH_CODE]


def test_one_row_per_person_and_ingredient_counting_mentions():
    t = build((1, ["d_metf", "d_metf_er", "d_combo"]))
    assert sorted(t.rxcui) == sorted([METF, LISIN])
    assert t.set_index("rxcui").mentions.to_dict() == {METF: 3, LISIN: 1}


def test_a_failed_drug_keeps_its_own_row_and_each_uncoded_row_is_one_drug():
    t = build((1, ["c00049", "c00049", "", "", "d_unknown"]))
    assert len(t) == 4  # one c-code row (2 mentions), two uncoded rows, one unknown code
    assert t[t.code == "c00049"].mentions.tolist() == [2]


def test_each_failure_cause_is_kept_apart():
    t = build((1, ["", "c00049", "d_retired", "d_unknown", "d_escit"]))
    by_cause = t.groupby("cause", dropna=False).size().to_dict()
    assert {k: v for k, v in by_cause.items() if isinstance(k, str)} == {
        ing.NO_CODE: 1, ing.ACH_CODE: 1, ing.RETIRED_CODE: 1, ing.NAME_NOT_FOUND: 1,
    }
    assert set(t.status) == {pdr.MAPPING_FAILED, pdr.NOT_IN_MODEL}
    assert t[t.status == pdr.NOT_IN_MODEL].rxcui.tolist() == [ESCIT]


def test_the_full_drug_list_stays():
    t = build((1, ["d_mph", "d_escit", "", "d_unknown"]))
    assert len(t) == 4 and set(t.status) == set(pdr.STATUSES)


def test_missing_columns_are_refused():
    with pytest.raises(ValueError, match="lack columns"):
        pdr.build(rows((1, ["d_mph"])).drop(columns="round"), resolve, STITCH_OF)


# --- pairs: one synthetic person per rule of §5 ---


def test_aliases_of_one_stitch_id_dont_form_a_pair():
    assert pairs_of((1, ["d_mph", "d_dexmph"])).empty
    p = pairs_of((1, ["d_mph", "d_dexmph", "d_metf"]))
    assert pair_set(p) == {(1, "CID000004091", "CID000004158", pdr.COVERED)}


def test_a_fixed_dose_combination_forms_a_pair():
    assert pair_set(pairs_of((1, ["d_combo"]))) == {(1, "CID000003961", "CID000004091", pdr.COVERED)}


def test_the_same_ingredient_twice_counts_once():
    p = pairs_of((1, ["d_metf", "d_metf_er", "d_lisin", "d_combo"]))
    assert len(p) == 1


def test_no_drug_pairs_with_itself():
    assert pairs_of((1, ["d_metf", "d_metf_er"])).empty


def test_failed_plus_not_in_model_is_mapping_failed():
    p = pairs_of((1, ["d_escit", "d_unknown"]))
    assert p.status.tolist() == [pdr.MAPPING_FAILED]


def test_pair_status_covered_not_in_model_and_failed():
    p = pairs_of((1, ["d_metf", "d_lisin", "d_escit", "d_rosu", "c00049"]))
    counts = p.status.value_counts().to_dict()
    # 5 drugs, 10 pairs: 1 covered; 2 covered × 2 out + 1 out × out = 5 not in the model;
    # 4 with the failed code.
    assert counts == {pdr.COVERED: 1, pdr.NOT_IN_MODEL: 5, pdr.MAPPING_FAILED: 4}


def test_covered_pairs_are_ordered_by_stitch_id():
    p = pairs_of((1, ["d_mph", "d_metf"]))
    assert (p.unit_a < p.unit_b).all() and p.unit_a.iloc[0] == "CID000004091"


def test_pairs_stay_within_a_person():
    assert pairs_of((1, ["d_metf"]), (2, ["d_lisin"])).empty


def test_drugs_in_different_rounds_dont_pair():
    assert pairs_of((1, 1, ["d_metf"]), (1, 2, ["d_lisin"])).empty


def test_a_pair_seen_in_two_rounds_counts_once():
    p = pairs_of((1, 1, ["d_metf", "d_lisin"]), (1, 2, ["d_combo"]))
    assert len(p) == 1


def test_people_without_two_drugs_have_no_pairs():
    assert pairs_of((1, ["d_metf"])).columns.tolist() == pdr.PAIR_COLUMNS


# --- NHANES 2017–2018 against the pilot (raw files and lookup cache needed) ---

CACHE = lookup.CACHE_DIR / pdr.NHANES_RXNAV_CACHE

needs_local = pytest.mark.skipif(
    not (nhanes.RAW_DIR.is_dir() and CACHE.exists()),
    reason="raw NHANES files or the NHANES RxNav cache not in local/",
)


def _offline(url, form):
    raise OSError(f"tests make no network calls: {url}")


@pytest.fixture(scope="module")
def nhanes_2017(tmp_path_factory):
    from sief_druco.burden import decagon_table

    no_overrides = tmp_path_factory.mktemp("overrides") / "none.csv"
    table = decagon_table.build(overrides_path=no_overrides, fetch=_offline)
    drugs = pdr.build_nhanes(["J"], table.stitch_of, fetch=_offline)
    persons = nhanes.load_persons(["J"])
    return drugs, pdr.pairs(drugs), persons


@needs_local
def test_2017_18_matches_the_pilot(nhanes_2017):
    drugs, pairs, persons = nhanes_2017
    w = persons.set_index("seqn").weight

    rows = nhanes.load_drug_rows(["J"])
    assert len(rows) == 14_434

    # "99.4% of mentions reach an IN" (spec §4) is the #12 count: 14,214 of the
    # 14,300 coded rows by the first MMSL lookup alone. The fallbacks add 35.
    coded = rows[rows.code != ""]
    reaching = set(drugs.code[drugs.rxcui.notna()])
    assert coded.code.isin(reaching).sum() == 14_249
    assert coded.code.isin(reaching).mean() > 0.994

    # Ingredient mentions (pilot §6: 15,560, of which 13,064 covered and 185 failed
    # on the NHANES side). The pilot's 3 Decagon-side failures (risedronate, found
    # by its name check) are *not in the model* here: the name check is #55's.
    m = drugs.mentions.groupby(drugs.status).sum()
    assert m.to_dict() == {pdr.COVERED: 13_064, pdr.NOT_IN_MODEL: 2_308 + 3, pdr.MAPPING_FAILED: 185}
    assert drugs[drugs.status == pdr.MAPPING_FAILED].cause.notna().all()

    # Person-pairs (pilot §6: 42,246; 29,977 covered, 69.2% of the weight). The 17
    # risedronate pairs move from mapping failed to not in the model. De-duplicating
    # by STITCH ID merges 28 distinct covered pairs (pilot 8,884) and no person-pair.
    assert pairs.status.value_counts().to_dict() == {
        pdr.COVERED: 29_977, pdr.NOT_IN_MODEL: 11_963 + 17, pdr.MAPPING_FAILED: 306 - 17,
    }
    pw = pairs.seqn.map(w)
    assert pw[pairs.status == pdr.COVERED].sum() / pw.sum() == pytest.approx(0.692, abs=0.0005)
    covered = pairs[pairs.status == pdr.COVERED]
    assert len(covered.drop_duplicates(["unit_a", "unit_b"])) == 8_884 - 28
