import json

import pandas as pd
import pytest

from sief_druco import decagon
from sief_druco.burden import decagon_table as dt
from sief_druco.burden import lookup

# --- A small fake PubChem / RxNav, so nothing here touches the network ---

HC, HC_BUTYRATE, MPH = "CID000003640", "CID000003642", "CID000004158"
UNKNOWN, PLAIN = "CID000009999", "CID000000100"
IDS = [HC, HC_BUTYRATE, MPH, UNKNOWN, PLAIN]
PAIRS = pd.Series({HC: 396, HC_BUTYRATE: 26, MPH: 300, UNKNOWN: 50, PLAIN: 10})

TITLES = {
    3640: "Hydrocortisone",
    3642: "Hydrocortisone butyrate",
    4158: "Methylphenidate",
    9999: "Mystery",
    100: "Plainol",
}
SYNONYMS = {
    3640: ["Hydrocortisone", "UNII-WI4X0X7BPJ"],
    3642: ["Hydrocortisone butyrate", "UNII-05RMF7YPWN"],
    # One record whose UNIIs reach both enantiomers' INs.
    4158: ["Methylphenidate", "UNII-207ZZ9QZ49", "UNII-M32RH9MFGP"],
    9999: ["Mystery"],
    100: ["Plainol", "UNII-PLAIN00001"],
}
UNII_RXCUI = {
    "WI4X0X7BPJ": "5492",
    "05RMF7YPWN": "1000",
    "207ZZ9QZ49": "6901",
    "M32RH9MFGP": "1091",
    "PLAIN00001": "777",
}
# related.json?tty=IN: a PIN (1000) lifts to its IN.
RELATED_IN = {
    "5492": [("5492", "hydrocortisone")],
    "1000": [("5492", "hydrocortisone")],
    "6901": [("6901", "methylphenidate")],
    "1091": [("1091", "dexmethylphenidate")],
    "777": [("777", "plainol")],
}
PROPERTIES = {
    "5492": ("hydrocortisone", "IN"),
    "6901": ("methylphenidate", "IN"),
    "1091": ("dexmethylphenidate", "IN"),
    "777": ("plainol", "IN"),
    "4242": ("mysterine", "IN"),
    "1000": ("hydrocortisone butyrate", "PIN"),
}


def _answer(url, form):
    if "property/Title" in url:
        cids = [c for c in map(int, form["cid"].split(",")) if c in TITLES]
        return {"PropertyTable": {"Properties": [{"CID": c, "Title": TITLES[c]} for c in cids]}}
    if "cid/synonyms" in url:
        cids = [c for c in map(int, form["cid"].split(",")) if c in SYNONYMS]
        return {"InformationList": {"Information": [{"CID": c, "Synonym": SYNONYMS[c]} for c in cids]}}
    if "idtype=UNII_CODE" in url:
        rxcui = UNII_RXCUI.get(url.rsplit("=", 1)[1])
        return {"idGroup": {"rxnormId": [rxcui]} if rxcui else {}}
    if "related.json?tty=IN" in url:
        rxcui = url.split("/rxcui/")[1].split("/")[0]
        props = [{"rxcui": i, "name": n} for i, n in RELATED_IN.get(rxcui, [])]
        return {"relatedGroup": {"conceptGroup": [{"tty": "IN", "conceptProperties": props}]}}
    if url.endswith("/properties.json"):
        rxcui = url.split("/rxcui/")[1].split("/")[0]
        if rxcui not in PROPERTIES:
            return {}
        name, tty = PROPERTIES[rxcui]
        return {"properties": {"rxcui": rxcui, "name": name, "tty": tty}}
    if "rxcui.json?name=" in url:
        return {"idGroup": {}}
    if "approximateTerm" in url:
        return {"approximateGroup": {}}
    return None  # parent / same-connectivity / SMILES / PUG-View: PubChem's 404


class FakeApi:
    def __init__(self):
        self.calls = []

    def __call__(self, url, form):
        self.calls.append(url)
        body = _answer(url, form)
        return (404, "PUGREST.NotFound") if body is None else (200, json.dumps(body))


def offline(url, form):
    raise AssertionError(f"network call: {url}")


def caches(fetch, cache_dir=None):
    if cache_dir is None:
        return lookup.Cache("unused-pc.json", fetch=fetch), lookup.Cache("unused-rx.json", fetch=fetch)
    return (
        lookup.Cache(cache_dir / dt.PUBCHEM_CACHE, fetch=fetch),
        lookup.Cache(cache_dir / dt.RXNAV_CACHE, fetch=fetch),
    )


@pytest.fixture
def auto():
    return dt.automatic_table(IDS, *caches(FakeApi()))


def write_overrides(path, rows):
    pd.DataFrame(rows, columns=dt.OVERRIDE_COLUMNS).to_csv(path, index=False)
    return path


def resolve(auto, overrides=None):
    _, rx = caches(FakeApi())
    if overrides is None:
        overrides = pd.DataFrame(columns=dt.OVERRIDE_COLUMNS)
    return dt.resolve(auto, overrides, PAIRS, rx)


def row(stitch_id, main_in="", aliases="", verdict="matched", reason="hand check", evidence="SIDER name"):
    return [stitch_id, main_in, aliases, verdict, reason, evidence]


# --- Automatic table ---


def test_every_drug_is_matched_or_unmatched_with_a_reason(auto):
    assert list(auto.columns) == dt.AUTO_COLUMNS
    assert sorted(set(auto.stitch_id)) == sorted(IDS)
    per_drug = auto.groupby("stitch_id").rxcui.agg(lambda c: c.notna().sum())
    assert per_drug.to_dict() == {HC: 1, HC_BUTYRATE: 1, MPH: 2, UNKNOWN: 0, PLAIN: 1}
    lost = auto[auto.rxcui.isna()]
    assert lost.stitch_id.tolist() == [UNKNOWN]
    assert "no UNII" in lost.reason.iloc[0]
    assert auto.reason[auto.rxcui.notna()].isna().all()


def test_cid_is_the_stitch_number_without_leading_zeros(auto):
    assert auto.set_index("stitch_id").cid.groupby(level=0).first()[HC] == 3640
    with pytest.raises(ValueError):
        dt.stitch_cid("CIDm00003640")


def test_a_salt_record_reaches_the_ingredient_through_its_pin(auto):
    butyrate = auto[auto.stitch_id == HC_BUTYRATE].iloc[0]
    assert (butyrate.rxcui, butyrate["name"], butyrate.match_method) == ("5492", "hydrocortisone", dt.EXACT)


# --- Resolution without overrides ---


def test_hydrocortisone_goes_to_the_drug_with_the_most_labelled_pairs(auto):
    table = resolve(auto)
    assert table.stitch_of["5492"] == HC
    footnote = table.collisions.set_index("rxcui").loc["5492"]
    assert footnote.stitch_id == HC
    assert footnote.labelled_pairs == 396
    assert footnote.others == f"{HC_BUTYRATE} (26)"
    # The ester still keeps hydrocortisone as its own main IN.
    assert table.ingredients.query("stitch_id == @HC_BUTYRATE").role.tolist() == ["main"]


def test_the_collision_follows_the_pair_counts_not_the_id(auto):
    flipped = PAIRS.copy()
    flipped[HC], flipped[HC_BUTYRATE] = 26, 396
    _, rx = caches(FakeApi())
    table = dt.resolve(auto, pd.DataFrame(columns=dt.OVERRIDE_COLUMNS), flipped, rx)
    assert table.stitch_of["5492"] == HC_BUTYRATE


def test_a_drug_with_several_ins_keeps_them_as_unchecked_aliases(auto):
    mph = resolve(auto).ingredients.query("stitch_id == @MPH").set_index("name")
    # The main IN is the one the drug's own names confirm, not the lowest RxCUI (1091).
    assert mph.loc["methylphenidate"].role == "main"
    assert pd.isna(mph.loc["methylphenidate"].relation)
    assert (mph.loc["dexmethylphenidate"].role, mph.loc["dexmethylphenidate"].relation) == ("alias", dt.UNCHECKED)


def test_every_reached_in_points_to_a_drug(auto):
    table = resolve(auto)
    assert set(table.stitch_of) == {"5492", "6901", "1091", "777"}
    assert table.stitch_of["1091"] == MPH
    assert table.unmatched.to_dict("records") == [
        {"stitch_id": UNKNOWN, "reason": auto.set_index("stitch_id").reason[UNKNOWN], "source": "automatic"}
    ]


def test_a_missing_pair_count_raises(auto):
    _, rx = caches(FakeApi())
    with pytest.raises(KeyError):
        dt.resolve(auto, pd.DataFrame(columns=dt.OVERRIDE_COLUMNS), PAIRS.drop(PLAIN), rx)


# --- Overrides ---


def test_a_missing_override_file_means_no_overrides(tmp_path):
    overrides = dt.load_overrides(IDS, tmp_path / "absent.csv")
    assert overrides.empty
    assert list(overrides.columns) == dt.OVERRIDE_COLUMNS


def test_overrides_match_drop_and_unmatch(auto, tmp_path):
    path = write_overrides(
        tmp_path / "o.csv",
        [
            row(UNKNOWN, "4242", "777:acid/salt", reason="SIDER name mysterine"),
            row(MPH, "6901", verdict="dropped_in", reason="enantiomers sold as separate drugs"),
            row(PLAIN, verdict="unmatched", reason="title names a different compound"),
        ],
    )
    table = resolve(auto, dt.load_overrides(IDS, path))
    ing = table.ingredients.set_index(["stitch_id", "rxcui"])
    assert ing.loc[(UNKNOWN, "4242")][["name", "role", "source"]].tolist() == ["mysterine", "main", "override"]
    assert ing.loc[(UNKNOWN, "777")][["role", "relation"]].tolist() == ["alias", "acid/salt"]
    assert table.ingredients.query("stitch_id == @MPH").rxcui.tolist() == ["6901"]
    assert "1091" not in table.stitch_of
    assert PLAIN not in set(table.ingredients.stitch_id)
    plain = table.unmatched.set_index("stitch_id").loc[PLAIN]
    assert plain.tolist() == ["title names a different compound", "override"]
    assert table.stitch_of["777"] == UNKNOWN  # now only the alias claims it


@pytest.mark.parametrize("aliases", ["777:enantiomer", "777:", "777", "777:Acid/salt"])
def test_an_alias_with_an_unknown_tag_raises(tmp_path, aliases):
    path = write_overrides(tmp_path / "o.csv", [row(UNKNOWN, "4242", aliases)])
    with pytest.raises(ValueError, match="alias"):
        dt.load_overrides(IDS, path)


def test_an_override_outside_the_decagon_drugs_raises(tmp_path):
    path = write_overrides(tmp_path / "o.csv", [row("CID000123456", "4242")])
    with pytest.raises(ValueError, match="outside"):
        dt.load_overrides(IDS, path)


@pytest.mark.parametrize(
    "rows",
    [
        [row(UNKNOWN, "4242", verdict="maybe")],
        [row(UNKNOWN, "4242"), row(UNKNOWN, "777")],
        [row(UNKNOWN, "")],
        [row(UNKNOWN, "mysterine")],
        [row(UNKNOWN, "4242", "4242:prodrug")],
        [row(UNKNOWN, "4242", "777:prodrug;777:acid/salt")],
        [row(UNKNOWN, "4242", verdict="unmatched")],
        [row(UNKNOWN, verdict="unmatched", reason="")],
    ],
    ids=[
        "verdict",
        "repeated",
        "no-main",
        "main-by-name",
        "main-as-alias",
        "alias-twice",
        "unmatched-with-main",
        "unmatched-no-reason",
    ],
)
def test_malformed_override_rows_raise(tmp_path, rows):
    with pytest.raises(ValueError):
        dt.load_overrides(IDS, write_overrides(tmp_path / "o.csv", rows))


def test_an_override_in_that_is_not_an_rxnorm_ingredient_raises(auto, tmp_path):
    for main in ["1000", "31337"]:  # a PIN, and an RxCUI RxNav doesn't know
        path = write_overrides(tmp_path / "o.csv", [row(UNKNOWN, main)])
        with pytest.raises(ValueError, match="not an RxNorm IN"):
            resolve(auto, dt.load_overrides(IDS, path))


# --- Cache ---


def test_a_second_run_with_a_full_cache_makes_no_network_call(tmp_path):
    overrides = write_overrides(tmp_path / "o.csv", [row(UNKNOWN, "4242", "777:acid/salt")])
    api = FakeApi()
    first = dt.build(IDS, PAIRS, overrides, tmp_path / "cache", api)
    assert api.calls  # the first run did ask
    second = dt.build(IDS, PAIRS, overrides, tmp_path / "cache", offline)
    pd.testing.assert_frame_equal(first.automatic, second.automatic)
    pd.testing.assert_frame_equal(first.ingredients, second.ingredients)
    assert first.stitch_of == second.stitch_of


def test_the_cache_stores_errors_and_retries_server_failures(tmp_path, monkeypatch):
    monkeypatch.setattr(lookup.time, "sleep", lambda s: None)
    answers = iter([(503, "busy"), (200, '{"ok": 1}')])
    cache = lookup.Cache(tmp_path / "c.json", fetch=lambda url, form: next(answers))
    assert cache.get("https://x/a") == {"ok": 1}
    cache = lookup.Cache(tmp_path / "c.json", fetch=lambda url, form: (404, "not found"))
    assert cache.get("https://x/b")["_status"] == 404
    cache.save()
    reread = lookup.Cache(tmp_path / "c.json", fetch=offline)
    assert reread.get("https://x/b")["_status"] == 404


def test_a_request_that_keeps_failing_is_not_cached(tmp_path, monkeypatch):
    monkeypatch.setattr(lookup.time, "sleep", lambda s: None)
    cache = lookup.Cache(tmp_path / "c.json", fetch=lambda url, form: (503, "busy"))
    assert cache.get("https://x/a") == {"_status": "failed"}
    assert "https://x/a" not in cache.data


# --- The real 645 (raw Decagon files and the lookup cache needed) ---

needs_local = pytest.mark.skipif(
    not (decagon.RAW_DIR.is_dir() and (lookup.CACHE_DIR / dt.PUBCHEM_CACHE).exists()),
    reason="raw Decagon files or the PubChem / RxNav cache not in local/",
)


@pytest.fixture(scope="module")
def real():
    # Offline: every lookup must come from the cache.
    return dt.build(overrides_path=lookup.CACHE_DIR / "no-overrides.csv", fetch=offline)


@needs_local
def test_without_overrides_the_pilot_numbers_come_back(real):
    per_drug = real.automatic.groupby("stitch_id").rxcui.agg(lambda c: c.notna().sum())
    assert len(per_drug) == 645
    assert (per_drug > 0).sum() == 611
    assert (per_drug == 0).sum() == 34
    assert (per_drug > 1).sum() == 17
    methods = real.automatic.drop_duplicates("stitch_id").match_method.value_counts()
    assert methods.to_dict() == {
        dt.EXACT: 321,
        dt.SAME_CONNECTIVITY: 234,
        dt.SAME_CONNECTIVITY_SMILES: 30,
        dt.NAME: 18,
        dt.PARENT: 7,
        dt.APPROXIMATE: 1,
    }
    assert len(real.stitch_of) == 630
    assert len(real.unmatched) == 34


@needs_local
def test_real_hydrocortisone_resolves_to_cid_3640(real):
    assert real.stitch_of["5492"] == "CID000003640"
    counts = dt.labelled_pair_counts(decagon.load())
    assert (counts["CID000003640"], counts["CID000003642"]) == (396, 26)
    assert real.collisions.rxcui.tolist() == ["5492"]
