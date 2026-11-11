"""The evidence sheet for the hand fixes to the Decagon mapping (#55).

Spec: `docs/spec-burden.md` §4 (hand fixes, name check) and §11.1 (lithium,
hydrocodone). Pilot: `research/nhanes-coverage-pilot.md` §3.3 and §7.

The sheet covers three kinds of Decagon drugs:

- `unmatched`: the drugs that reached no RxNorm `IN`;
- `multi_in`: the drugs that reached more than one `IN`;
- `name_exact` / `name_loose`: the hits of the **name check**. Source units that
  are not covered (an `IN` outside the 645, or a failed source name) are compared
  with the Decagon drugs' names: exact (case-insensitive, `IN` names and failed
  names), then loose (an `IN` name as a whole word inside a Decagon name).

Decagon names are those of pilot §3.3 (PubChem title, the first 30 synonyms of the
exact and parent CID, the `IN` names reached; for an unmatched drug the first 10
synonyms of each same-connectivity CID) plus the SIDER 4.1 name of the same flat
CID.

Each sheet row is one Decagon drug × one candidate `IN`: an `IN` the automatic
route reached, an `IN` the SIDER name reaches, or the source `IN` of a name hit.
Its evidence: the SIDER name, PubChem formula and InChIKey of the Decagon CID next
to those of the PubChem record the `IN` name resolves to, and the `IN`'s RxNav
type and status. Same connectivity = the first 14 characters of the InChIKeys are
equal (same molecule up to stereochemistry).

The verdicts are proposed by hand from this sheet and written to
`data/decagon_overrides.csv`; nothing here writes that file. The sheet itself goes
to `local/` (`python -m sief_druco.burden.handfix`).
"""

import re
from pathlib import Path
from urllib.parse import quote

import pandas as pd

from sief_druco.burden import decagon_table as dt
from sief_druco.burden import lookup, nhanes
from sief_druco.burden import person_drugs as pdr
from sief_druco.burden.ingredients import NO_CODE, Resolver
from sief_druco.burden.lookup import PUG, RXNAV, Cache

SIDER_NAMES = lookup.LOCAL_DIR / "data" / "raw" / "sider" / "drug_names.tsv"
MEPS_DIR = lookup.LOCAL_DIR / "data" / "raw" / "meps"
OUT_DIR = lookup.LOCAL_DIR / "data" / "handfix"

MEPS_FILES = {2015: "h178a", 2016: "h188a", 2017: "h197a"}  # §11.2: 2015-2017 only
MEPS_RXNAV_CACHE = "rxnav_meps.json"
MEPS_SEPARATORS = ("-",)  # §4: MEPS combinations split at '-'

UNMATCHED, MULTI_IN, NAME_EXACT, NAME_LOOSE = "unmatched", "multi_in", "name_exact", "name_loose"
FROM_AUTO, FROM_SIDER, FROM_NAME_HIT = "automatic", "sider_name", "name_hit"

UNIT_COLUMNS = ["source", "kind", "unit_name", "rxcui", "mentions", "persons"]
HIT_COLUMNS = [*UNIT_COLUMNS, "check", "stitch_id", "decagon_name"]
SHEET_COLUMNS = [
    "stitch_id", "why", "title", "sider_name", "formula", "inchikey", "auto_ins", "auto_reason",
    "candidate_rxcui", "candidate_name", "candidate_from", "candidate_tty", "candidate_status",
    "ref_cid", "ref_formula", "ref_inchikey", "same_formula", "same_connectivity", "same_inchikey",
    "hit_units", "labelled_pairs",
]

_N_SYNONYMS, _N_SAME_SYNONYMS = 30, 10
_CHUNK = 100


# --- Names ---


def sider_names(path: Path = SIDER_NAMES) -> dict[int, str]:
    """SIDER 4.1 `drug_names.tsv`: flat STITCH ID (`CID1` + CID) → name, by CID."""
    t = pd.read_csv(path, sep="\t", header=None, names=["flat", "name"], dtype=str)
    return {int(f[4:]): n for f, n in zip(t.flat, t.name)}


def decagon_names(table: dt.DrugTable, sider: dict[int, str], pubchem: Cache, rxnav: Cache) -> pd.DataFrame:
    """`stitch_id`, `name` (lower case): the names the name check compares against."""
    route = dt._Route(pubchem, rxnav)
    ids = sorted(table.automatic.stitch_id.unique())
    cid_of = {s: dt.stitch_cid(s) for s in ids}
    titles = route.titles(cid_of.values())
    route.fetch_synonyms(cid_of.values())
    unmatched = set(table.unmatched.stitch_id)
    reached = table.automatic.dropna(subset=["name"]).groupby("stitch_id").name.agg(list)

    rows = []
    for s, cid in cid_of.items():
        names = [titles.get(cid, ""), *route.synonyms.get(cid, [])[:_N_SYNONYMS]]
        parents = [p for p in route.related_cids(cid, "parent") if p != cid][:1]
        if parents:
            route.fetch_synonyms(parents)
            names += route.synonyms.get(parents[0], [])[:_N_SYNONYMS]
        if s in unmatched:
            same = [x for x in route.related_cids(parents[0] if parents else cid, "same_connectivity")
                    if x != cid]
            route.fetch_synonyms(same)
            for x in same:
                names += route.synonyms.get(x, [])[:_N_SAME_SYNONYMS]
        names += reached.get(s, [])
        names.append(sider.get(cid, ""))
        rows += [(s, n) for n in dict.fromkeys(n.strip().lower() for n in names) if n]
    return pd.DataFrame(rows, columns=["stitch_id", "name"])


# --- Source units and the name check ---


def source_units(person_drugs: pd.DataFrame) -> pd.DataFrame:
    """Uncovered source units with mentions and persons (`UNIT_COLUMNS`).

    `kind` = `in` (reached an `IN` outside the 645; `unit_name` is the `IN` name) or
    `failed` (reached no `IN`; `unit_name` is the source name). Uncoded NHANES rows
    and MEPS missing-value codes (negative numbers) have no name to check.
    """
    t = person_drugs[person_drugs.status != pdr.COVERED]
    t = t[(t.cause != NO_CODE) & ~t.name.fillna("").str.match(r"^-\d+$")]
    t = t.assign(
        kind=t.status.map({pdr.NOT_IN_MODEL: "in", pdr.MAPPING_FAILED: "failed"}),
        unit_name=t.in_name.where(t.status == pdr.NOT_IN_MODEL, t.name).fillna("").str.strip(),
        person=t.cycle.astype(str) + "/" + t.seqn.astype(str),
    )
    t = t[t.unit_name != ""]
    out = t.groupby(["source", "kind", "unit_name"], dropna=False).agg(
        rxcui=("rxcui", "first"), mentions=("mentions", "sum"), persons=("person", "nunique")
    )
    return out.reset_index()[UNIT_COLUMNS].sort_values(["source", "mentions"], ascending=[True, False])


def name_check(units: pd.DataFrame, names: pd.DataFrame) -> pd.DataFrame:
    """Exact hits (any unit) and loose hits (`IN` units only) against Decagon names."""
    by_name = names.groupby("name").stitch_id.agg(sorted)
    hits = []
    for u in units.itertuples(index=False):
        key = u.unit_name.lower()
        exact = by_name.get(key, [])
        hits += [(*u, "exact", s, key) for s in exact]
        if u.kind != "in":
            continue
        word = re.compile(rf"(?<![a-z0-9]){re.escape(key)}(?![a-z0-9])")
        loose = names[names.name.str.contains(word) & (names.name != key) & ~names.stitch_id.isin(exact)]
        for s, group in loose.groupby("stitch_id"):
            hits.append((*u, "loose", s, "; ".join(sorted(group.name)[:3])))
    return pd.DataFrame(hits, columns=HIT_COLUMNS)


# --- MEPS drug rows (names only; the MEPS tables are #62) ---


def load_meps_drug_rows(years=tuple(MEPS_FILES), raw_dir: Path = MEPS_DIR) -> pd.DataFrame:
    """PMED records with a positive person weight, as `person_drugs.build` takes them."""
    tables = []
    for year in years:
        f = MEPS_FILES[year]
        d = pd.read_sas(Path(raw_dir) / f / f"{f}.ssp", format="xport", encoding="latin-1")
        d = d[d[f"PERWT{year % 100}F"] > 0]
        tables.append(pd.DataFrame({
            "source": "MEPS",
            "cycle": str(year),
            "seqn": d.DUPERSID.astype(str),
            "round": d.PURCHRD.astype("int64"),
            "code": "",
            "name": d.RXDRGNAM.fillna("").astype(str).str.strip(),
        }))
    return pd.concat(tables, ignore_index=True)


# --- Evidence ---


class _Evidence:
    def __init__(self, pubchem: Cache, rxnav: Cache):
        self.pc, self.rx = pubchem, rxnav

    def cid_properties(self, cids) -> dict[int, dict]:
        out, cids = {}, sorted(set(cids))
        for i in range(0, len(cids), _CHUNK):
            chunk = ",".join(map(str, cids[i : i + _CHUNK]))
            answer = self.pc.get(f"{PUG}/compound/cid/property/Title,MolecularFormula,InChIKey/JSON",
                                 {"cid": chunk})
            for p in (answer.get("PropertyTable") or {}).get("Properties") or []:
                out[p["CID"]] = p
        return out

    def name_records(self, name: str) -> list[dict]:
        answer = self.pc.get(f"{PUG}/compound/name/{quote(name)}/property/MolecularFormula,InChIKey/JSON")
        return (answer.get("PropertyTable") or {}).get("Properties") or []

    def rx_status(self, rxcui: str) -> tuple[str, str, str]:
        """(name, tty, status) of an RxCUI."""
        props = lookup.rx_properties(self.rx, rxcui)
        hist = self.rx.get(f"{RXNAV}/rxcui/{rxcui}/historystatus.json")
        status = ((hist.get("rxcuiStatusHistory") or {}).get("metaData") or {}).get("status", "")
        return props.get("name", ""), props.get("tty", ""), status

    def name_ins(self, name: str) -> list[tuple[str, str]]:
        answer = self.rx.get(f"{RXNAV}/rxcui.json?name={quote(name)}&search=2")
        return [h for rc in lookup.rx_ids(answer) for h in lookup.rx_ingredients(self.rx, rc)]


def compare(decagon: dict, refs: list[dict]) -> dict:
    """The PubChem record of an `IN` name closest to the Decagon CID, and the matches.

    Prefers a record with the same InChIKey, then the same connectivity, then the
    same formula, else the first record.
    """
    key, formula = decagon.get("InChIKey", ""), decagon.get("MolecularFormula", "")

    def score(r):
        k = r.get("InChIKey", "")
        return (bool(key) and k == key, bool(key) and k[:14] == key[:14], r.get("MolecularFormula") == formula)

    if not refs:
        return {"ref_cid": None, "ref_formula": "", "ref_inchikey": "",
                "same_formula": False, "same_connectivity": False, "same_inchikey": False}
    best = max(refs, key=score)  # max keeps the first of equal scores
    same_key, same_conn, same_formula = score(best)
    return {"ref_cid": best.get("CID"), "ref_formula": best.get("MolecularFormula", ""),
            "ref_inchikey": best.get("InChIKey", ""), "same_formula": same_formula,
            "same_connectivity": same_conn, "same_inchikey": same_key}


def evidence_sheet(
    table: dt.DrugTable,
    hits: pd.DataFrame,
    sider: dict[int, str],
    pair_counts: pd.Series,
    pubchem: Cache,
    rxnav: Cache,
) -> pd.DataFrame:
    """One row per hand-fix drug × candidate `IN` (`SHEET_COLUMNS`)."""
    ev = _Evidence(pubchem, rxnav)
    auto = table.automatic
    reached = auto.dropna(subset=["rxcui"])
    n_ins = reached.groupby("stitch_id").rxcui.nunique()
    why = {s: [UNMATCHED] for s in table.unmatched.stitch_id}
    for s in n_ins[n_ins > 1].index:
        why.setdefault(s, []).append(MULTI_IN)
    for h in hits.itertuples(index=False):
        tag = NAME_EXACT if h.check == "exact" else NAME_LOOSE
        if tag not in why.setdefault(h.stitch_id, []):
            why[h.stitch_id].append(tag)

    props = ev.cid_properties(dt.stitch_cid(s) for s in why)
    rows = []
    for s in sorted(why):
        cid = dt.stitch_cid(s)
        p = props.get(cid, {})
        mine = reached[reached.stitch_id == s]
        reason = auto.loc[auto.stitch_id == s, "reason"].dropna()
        drug_hits = hits[hits.stitch_id == s]
        candidates: dict[str, list[str]] = {}
        for rxcui in mine.rxcui:
            candidates.setdefault(rxcui, []).append(FROM_AUTO)
        if s in table.unmatched.stitch_id.values and sider.get(cid):
            for rxcui, _ in ev.name_ins(sider[cid]):
                candidates.setdefault(rxcui, []).append(FROM_SIDER)
        for rxcui in drug_hits.rxcui.dropna():
            candidates.setdefault(rxcui, []).append(FROM_NAME_HIT)
        base = {
            "stitch_id": s, "why": ", ".join(why[s]), "title": p.get("Title", ""),
            "sider_name": sider.get(cid, ""), "formula": p.get("MolecularFormula", ""),
            "inchikey": p.get("InChIKey", ""),
            "auto_ins": "; ".join(f"{r.rxcui} {r.name}" for r in mine.itertuples()),
            "auto_reason": reason.iloc[0] if len(reason) else "",
            "labelled_pairs": int(pair_counts.get(s, 0)),
        }
        if not candidates:
            rows.append(base | {"hit_units": _hit_units(drug_hits)} | compare(p, []))
        for rxcui, origin in sorted(candidates.items(), key=lambda kv: int(kv[0])):
            name, tty, status = ev.rx_status(rxcui)
            unit_hits = drug_hits[drug_hits.rxcui == rxcui] if FROM_NAME_HIT in origin else drug_hits.iloc[:0]
            rows.append(base | {
                "candidate_rxcui": rxcui, "candidate_name": name, "candidate_from": "+".join(dict.fromkeys(origin)),
                "candidate_tty": tty, "candidate_status": status,
                "hit_units": _hit_units(unit_hits if len(unit_hits) else drug_hits[drug_hits.rxcui.isna()]),
            } | compare(p, ev.name_records(name) if name else []))
    return pd.DataFrame(rows).reindex(columns=SHEET_COLUMNS)


def _hit_units(hits: pd.DataFrame) -> str:
    return "; ".join(
        f"{h.source} {h.check} '{h.unit_name}' ({h.mentions} mentions)" for h in hits.itertuples(index=False)
    )


# --- The run ---


def main(out_dir: Path = OUT_DIR, cache_dir: Path = lookup.CACHE_DIR) -> None:
    from sief_druco import decagon

    d = decagon.load()
    counts = dt.labelled_pair_counts(d)
    # The automatic table, before any override: the sheet is about what it reached.
    table = dt.build(d.drug_ids, counts, overrides_path=Path(out_dir) / "no-overrides.csv", cache_dir=cache_dir)
    pc = Cache(Path(cache_dir) / dt.PUBCHEM_CACHE, lookup.PUBCHEM_INTERVAL)
    rx = Cache(Path(cache_dir) / dt.RXNAV_CACHE, lookup.RXNAV_INTERVAL)
    meps_rx = Cache(Path(cache_dir) / MEPS_RXNAV_CACHE, lookup.RXNAV_INTERVAL)
    try:
        sider = sider_names()
        names = decagon_names(table, sider, pc, rx)
        nhanes_drugs = pdr.build_nhanes(nhanes.MAIN_POOL, table.stitch_of, cache_dir=cache_dir)
        resolver = Resolver(meps_rx)
        meps_drugs = pdr.build(
            load_meps_drug_rows(), lambda _, name: resolver.name(name, MEPS_SEPARATORS), table.stitch_of
        )
        units = pd.concat([source_units(nhanes_drugs), source_units(meps_drugs)], ignore_index=True)
        hits = name_check(units, names)
        sheet = evidence_sheet(table, hits, sider, counts, pc, rx)
    finally:
        for c in (pc, rx, meps_rx):
            c.save()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    units.to_csv(out_dir / "uncovered_units.csv", index=False)
    hits.to_csv(out_dir / "name_hits.csv", index=False)
    sheet.to_csv(out_dir / "evidence.csv", index=False)
    print(f"{sheet.stitch_id.nunique()} drugs, {len(sheet)} rows; {len(hits)} name hits -> {out_dir}")


if __name__ == "__main__":
    main()
