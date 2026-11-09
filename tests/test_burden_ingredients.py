import json
from urllib.parse import quote

import pytest

from sief_druco.burden import ingredients as ing
from sief_druco.burden import lookup
from sief_druco.burden.lookup import RXNAV

# --- A small fake RxNav, so nothing here touches the network ---

# MMSL_CODE lookups: code -> RxCUIs, without and with allsrc=1.
MMSL = {"d00001": ["10"], "d00002": ["20"], "d00004": ["40"]}
MMSL_ALLSRC = {"d00003": ["30"], "d00005": ["50"]}
# Name lookups (search=2): name as written -> RxCUIs.
NAMES = {
    "ALPHA": ["11"],
    "BETA": ["12"],
    "ALPHA-BETA": [],
    "GAMMA": ["13"],
    "DELTA-OMEGA": ["14"],
}
# related.json?tty=IN. 10 is a two-ingredient product; 20 is a PIN; 30 and 50 are
# retired concepts with no IN (30 is remapped, 50 is not).
RELATED_IN = {
    "10": [("1", "alpha"), ("2", "beta")],
    "20": [("3", "gamma")],
    "40": [("4", "delta")],
    "11": [("1", "alpha")],
    "12": [("2", "beta")],
    "13": [("3", "gamma")],
    "14": [("5", "delta-omega")],
    "31": [("6", "epsilon")],
}
REMAPPED = {"30": ["31"]}


def _answer(url):
    if "idtype=MMSL_CODE" in url:
        code = url.split("&id=GN")[1].split("&")[0]
        table = MMSL_ALLSRC if "allsrc=1" in url else MMSL
        return {"idGroup": {"rxnormId": table[code]} if code in table else {}}
    if "rxcui.json?name=" in url:
        for name, rxcuis in NAMES.items():
            if url == f"{RXNAV}/rxcui.json?name={quote(name)}&search=2":
                return {"idGroup": {"rxnormId": rxcuis} if rxcuis else {}}
        return {"idGroup": {}}
    rxcui = url.split("/rxcui/")[1].split("/")[0]
    if url.endswith("related.json?tty=IN"):
        props = [{"rxcui": i, "name": n} for i, n in RELATED_IN.get(rxcui, [])]
        return {"relatedGroup": {"conceptGroup": [{"tty": "IN", "conceptProperties": props}]}}
    if url.endswith("historystatus.json"):
        remapped = [{"remappedRxCui": r} for r in REMAPPED.get(rxcui, [])]
        return {"rxcuiStatusHistory": {"metaData": {"status": "NotCurrent"},
                                       "derivedConcepts": {"remappedConcept": remapped}}}
    raise AssertionError(f"unexpected request {url}")


@pytest.fixture
def resolver(tmp_path):
    calls = []

    def fetch(url, form):
        calls.append(url)
        return 200, json.dumps(_answer(url))

    r = ing.Resolver(lookup.Cache(tmp_path / "rxnav.json", fetch=fetch))
    r.calls = calls
    return r


# --- Multum codes (NHANES) ---


def test_a_combination_product_reaches_each_of_its_ingredients(resolver):
    res = resolver.multum("d00001", "ALPHA; BETA")
    assert res.ingredients == (("1", "alpha"), ("2", "beta"))
    assert res.cause is None and res.step == ing.MMSL


def test_a_pin_lifts_to_its_ingredient(resolver):
    assert resolver.multum("d00002", "GAMMA HCL").ingredients == (("3", "gamma"),)


def test_allsrc_and_the_history_remap_are_tried_before_the_name(resolver):
    res = resolver.multum("d00003", "UNKNOWN NAME")
    assert res.ingredients == (("6", "epsilon"),) and res.step == ing.MMSL_ALLSRC


def test_the_name_parts_are_tried_when_the_code_reaches_nothing(resolver):
    res = resolver.multum("d09999", "ALPHA; GAMMA")
    assert res.ingredients == (("1", "alpha"), ("3", "gamma")) and res.step == ing.NAME


@pytest.mark.parametrize("name", ["55555", "77777", "99999"])
def test_the_special_name_codes_have_no_code(resolver, name):
    res = resolver.multum("", name)
    assert res.ingredients == () and res.cause == ing.NO_CODE
    assert resolver.calls == []


@pytest.mark.parametrize("code", ["a54115", "c00049", "h00015"])
def test_a_c_h_codes_are_not_looked_up(resolver, code):
    res = resolver.multum(code, "ALPHA")
    assert res.ingredients == () and res.cause == ing.ACH_CODE
    assert resolver.calls == []


def test_a_retired_code_with_no_remap_is_its_own_cause(resolver):
    assert resolver.multum("d00005", "NOWHERE").cause == ing.RETIRED_CODE


def test_a_code_and_name_rxnorm_doesnt_know_are_name_not_found(resolver):
    assert resolver.multum("d08888", "NOWHERE").cause == ing.NAME_NOT_FOUND


def test_a_name_match_for_only_some_parts_fails_as_a_whole(resolver):
    res = resolver.multum("d08888", "ALPHA; NOWHERE")
    assert res.ingredients == () and res.cause == ing.NAME_NOT_FOUND


def test_each_cause_is_one_of_the_four():
    assert set(ing.CAUSES) == {ing.NO_CODE, ing.ACH_CODE, ing.RETIRED_CODE, ing.NAME_NOT_FOUND}


# --- Names (MEPS, DDInter) ---


def test_a_name_is_tried_whole_before_it_is_split(resolver):
    res = resolver.name("DELTA-OMEGA", separators=("-",))
    assert res.ingredients == (("5", "delta-omega"),)


def test_a_combination_name_splits_into_ingredients(resolver):
    res = resolver.name("ALPHA-BETA", separators=("-",))
    assert res.ingredients == (("1", "alpha"), ("2", "beta"))


def test_without_separators_a_name_is_not_split(resolver):
    assert resolver.name("ALPHA-BETA").cause == ing.NAME_NOT_FOUND


def test_an_empty_name_is_not_found_without_a_request(resolver):
    assert resolver.name("  ").cause == ing.NAME_NOT_FOUND
    assert resolver.calls == []


def test_answers_come_from_the_cache_on_a_second_resolver(resolver, tmp_path):
    resolver.multum("d00001", "ALPHA; BETA")
    resolver.rx.save()

    def offline(url, form):
        raise AssertionError(f"request sent with a full cache: {url}")

    again = ing.Resolver(lookup.Cache(tmp_path / "rxnav.json", fetch=offline))
    assert again.multum("d00001", "ALPHA; BETA").ingredients == (("1", "alpha"), ("2", "beta"))
