import json
from urllib.parse import quote

import pandas as pd
import pytest

from sief_druco import decagon
from sief_druco.burden import decagon_table as dt
from sief_druco.burden import handfix as hf
from sief_druco.burden import ingredients as ing
from sief_druco.burden import lookup
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden.ingredients import Resolution
from sief_druco.burden.lookup import PUG, RXNAV, Cache


def offline(url, form):
    raise AssertionError(f"network call in a test: {url}")


def _cache(tmp_path, name, answers):
    """A cache that answers from `answers` (url -> JSON) and fails on anything else."""
    path = tmp_path / name
    path.write_text(json.dumps(answers))
    return Cache(path, 0.0, offline)


# --- Source units and the name check ---


def _person_drugs():
    res = {
        "d_rosu": Resolution((("301542", "rosuvastatin"),), None, ing.MMSL),
        "d_li": Resolution((("6448", "lithium"),), None, ing.MMSL),
        "d_metf": Resolution((("6809", "metformin"),), None, ing.MMSL),
        "d_pen": Resolution((), ing.NAME_NOT_FOUND, ""),
        "": Resolution((), ing.NO_CODE, ""),
    }
    rows = pd.DataFrame(
        [
            ["NHANES", "I", 1, pd.NA, "d_rosu", "ROSUVASTATIN"],
            ["NHANES", "I", 2, pd.NA, "d_rosu", "ROSUVASTATIN"],
            ["NHANES", "I", 2, pd.NA, "d_li", "LITHIUM"],
            ["NHANES", "I", 2, pd.NA, "d_metf", "METFORMIN"],
            ["NHANES", "P_", 3, pd.NA, "d_pen", "PENICILLIN"],
            ["NHANES", "P_", 3, pd.NA, "", "99999"],
        ],
        columns=pdr.ROW_COLUMNS,
    )
    return pdr.build(rows, lambda code, name: res[code], {"6809": "CID000004091"})


def test_source_units_are_the_uncovered_named_units():
    units = hf.source_units(_person_drugs())
    got = {(u.kind, u.unit_name): (u.mentions, u.persons) for u in units.itertuples()}
    # Covered metformin and the uncoded row are left out.
    assert got == {("in", "rosuvastatin"): (2, 2), ("in", "lithium"): (1, 1), ("failed", "PENICILLIN"): (1, 1)}


def test_name_check_exact_then_loose_on_whole_words():
    units = hf.source_units(_person_drugs())
    names = pd.DataFrame(
        [
            ["CID000002349", "penicillins"],
            ["CID000002349", "penicillin"],
            ["CID000028112", "lithium aluminum hydride"],
            ["CID004475485", "rosuvastatin ep impurity l"],
            ["CID000000001", "lithiumate"],  # not a whole word
        ],
        columns=["stitch_id", "name"],
    )
    hits = hf.name_check(units, names)
    got = {(h.unit_name, h.check, h.stitch_id) for h in hits.itertuples()}
    assert got == {
        ("PENICILLIN", "exact", "CID000002349"),
        ("lithium", "loose", "CID000028112"),
        ("rosuvastatin", "loose", "CID004475485"),
    }


def test_failed_names_get_no_loose_check():
    units = pd.DataFrame([["NHANES", "failed", "UREA", None, 1, 1]], columns=hf.UNIT_COLUMNS)
    names = pd.DataFrame([["CID000000001", "urea cream"]], columns=["stitch_id", "name"])
    assert hf.name_check(units, names).empty


# --- Comparing structures ---

DEC = {"CID": 838, "MolecularFormula": "C9H13NO3", "InChIKey": "UCTWMZQNUQWSLP-UHFFFAOYSA-N"}


def test_compare_prefers_the_same_connectivity():
    refs = [
        {"CID": 1, "MolecularFormula": "C9H13NO3", "InChIKey": "XXXXXXXXXXXXXX-UHFFFAOYSA-N"},
        {"CID": 5816, "MolecularFormula": "C9H13NO3", "InChIKey": "UCTWMZQNUQWSLP-VIFPVBQESA-N"},
    ]
    got = hf.compare(DEC, refs)
    assert got["ref_cid"] == 5816
    assert (got["same_formula"], got["same_connectivity"], got["same_inchikey"]) == (True, True, False)


def test_compare_without_a_record():
    got = hf.compare(DEC, [])
    assert got["ref_cid"] is None and not got["same_formula"]


# --- The evidence sheet ---


def _props(rxcui, name):
    return {
        f"{RXNAV}/rxcui/{rxcui}/properties.json": {"properties": {"rxcui": rxcui, "name": name, "tty": "IN"}},
        f"{RXNAV}/rxcui/{rxcui}/historystatus.json": {"rxcuiStatusHistory": {"metaData": {"status": "Active"}}},
    }


def test_evidence_sheet_rows_and_candidates(tmp_path):
    automatic = pd.DataFrame(
        [
            ["CID000000838", 838, "327361", "racepinephrine", dt.SAME_CONNECTIVITY, False, None],
            ["CID000003003", 3003, "1514", "betamethasone", dt.SAME_CONNECTIVITY, False, None],
            ["CID000003003", 3003, "3264", "dexamethasone", dt.SAME_CONNECTIVITY, False, None],
            ["CID000041774", 41774, None, None, None, False, "approximateTerm hit not confirmed"],
        ],
        columns=dt.AUTO_COLUMNS,
    )
    table = dt.DrugTable(
        automatic, pd.DataFrame(), pd.DataFrame([["CID000041774", "x", "automatic"]],
                                                columns=["stitch_id", "reason", "source"]), {}, pd.DataFrame()
    )
    hits = pd.DataFrame(
        [["NHANES", "in", "epinephrine", "3992", 3, 3, "loose", "CID000000838", "(+-)-adrenaline"]],
        columns=hf.HIT_COLUMNS,
    )
    sider = {838: "epinephrine", 3003: "dexamethasone", 41774: "acarbose"}
    pc = _cache(tmp_path, "pc.json", {
        Cache.key(f"{PUG}/compound/cid/property/Title,MolecularFormula,InChIKey/JSON",
                  {"cid": "838,3003,41774"}): {"PropertyTable": {"Properties": [
            DEC | {"Title": "(+-)-Adrenaline"},
            {"CID": 3003, "Title": "Dexamethasone", "MolecularFormula": "C22H29FO5",
             "InChIKey": "UREBDLICKHMUKA-CXSFZGCWSA-N"},
            {"CID": 41774, "Title": "Acarbose", "MolecularFormula": "C25H43NO18",
             "InChIKey": "XUFXOAAUWZOOIT-UHFFFAOYSA-N"},
        ]}},
        **{f"{PUG}/compound/name/{quote(n)}/property/MolecularFormula,InChIKey/JSON": {"_status": 404}
           for n in ["racepinephrine", "betamethasone", "dexamethasone", "acarbose"]},
        f"{PUG}/compound/name/epinephrine/property/MolecularFormula,InChIKey/JSON": {"PropertyTable": {
            "Properties": [{"CID": 5816, "MolecularFormula": "C9H13NO3",
                            "InChIKey": "UCTWMZQNUQWSLP-VIFPVBQESA-N"}]}},
    })
    rx = _cache(tmp_path, "rx.json", {
        **_props("327361", "racepinephrine"), **_props("1514", "betamethasone"),
        **_props("3264", "dexamethasone"), **_props("3992", "epinephrine"), **_props("149", "acarbose"),
        f"{RXNAV}/rxcui.json?name=acarbose&search=2": {"idGroup": {"rxnormId": ["149"]}},
        f"{RXNAV}/rxcui/149/related.json?tty=IN": {"relatedGroup": {"conceptGroup": [
            {"conceptProperties": [{"rxcui": "149", "name": "acarbose"}]}]}},
    })
    counts = pd.Series({"CID000000838": 50, "CID000003003": 400, "CID000041774": 20})
    sheet = hf.evidence_sheet(table, hits, sider, counts, pc, rx)

    got = {(r.stitch_id, r.candidate_rxcui): (r.why, r.candidate_from) for r in sheet.itertuples()}
    assert got == {
        ("CID000000838", "327361"): ("name_loose", "automatic"),
        ("CID000000838", "3992"): ("name_loose", "name_hit"),
        ("CID000003003", "1514"): ("multi_in", "automatic"),
        ("CID000003003", "3264"): ("multi_in", "automatic"),
        ("CID000041774", "149"): ("unmatched", "sider_name"),
    }
    epi = sheet[sheet.candidate_rxcui == "3992"].iloc[0]
    assert epi.same_connectivity and not epi.same_inchikey and epi.ref_cid == 5816
    assert "3 mentions" in epi.hit_units
    assert list(sheet.columns) == hf.SHEET_COLUMNS


# --- The committed override file ---

LOCAL = decagon.RAW_DIR.is_dir() and (lookup.CACHE_DIR / dt.PUBCHEM_CACHE).exists()
needs_local = pytest.mark.skipif(not LOCAL, reason="raw Decagon files or the lookup cache not in local/")


def test_the_committed_overrides_parse():
    raw = pd.read_csv(dt.OVERRIDES_PATH, dtype=str, keep_default_na=False)
    overrides = dt.load_overrides(raw.stitch_id)
    assert len(overrides) == len(raw) > 0
    assert (overrides.evidence != "").all() and (overrides.reason != "").all()


@pytest.fixture(scope="module")
def with_overrides():
    return dt.build(fetch=offline)


@needs_local
def test_every_hand_fix_drug_has_an_override(with_overrides):
    auto = with_overrides.automatic
    n_ins = auto.dropna(subset=["rxcui"]).groupby("stitch_id").rxcui.nunique()
    unmatched = set(auto.stitch_id) - set(n_ins.index)
    needed = unmatched | set(n_ins[n_ins > 1].index)
    assert len(unmatched) == 34 and len(needed) == 34 + 17
    overridden = set(with_overrides.ingredients.query("source == 'override'").stitch_id)
    overridden |= set(with_overrides.unmatched.query("source == 'override'").stitch_id)
    assert needed <= overridden
    assert (with_overrides.ingredients.relation != dt.UNCHECKED).all()


@needs_local
def test_no_epimers_or_diastereomers_are_joined(with_overrides):
    ins = with_overrides.ingredients.groupby("stitch_id").rxcui.agg(set)
    for pair in [{"1514", "3264"}, {"9068", "9071"}]:  # betamethasone / dexamethasone, quinidine / quinine
        assert not any(pair <= s for s in ins), pair
