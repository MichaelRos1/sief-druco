"""The Decagon drug table: each of the 645 STITCH IDs → RxNorm ingredients (`IN`).

Spec: `docs/spec-burden.md` §4 (step 1 and the rules). Facts and the pilot this is
built from: `research/nhanes-stitch-mapping.md`, `research/nhanes-coverage-pilot.md`
§3.1 and §4.

Three layers:

1. `automatic_table`: the reverse route STITCH ID → PubChem CID → UNII → RxNorm
   `IN`, with the fallbacks of §4. One row per drug × `IN` reached; a drug that
   reached nothing has one row with no `rxcui` and a `reason`.
2. `load_overrides`: the hand fixes in `data/decagon_overrides.csv` (#55). A
   missing file means no overrides.
3. `resolve`: STITCH ID → one main `IN` + aliases, each tagged with its relation;
   and `IN` → STITCH ID, where an `IN` claimed by several drugs goes to the one with
   the most labelled pairs, and the collision is kept for a footnote.

Without an override, a drug that reached several `IN`s keeps all of them
provisionally (as the pilot did): the main one is the `IN` whose name the drug's
own PubChem names confirm (else the lowest RxCUI), the others are aliases tagged
`unchecked`. These are not alias relations under the rule of §4; #55 gives every
such drug an override, so none is left by the final run (#64).
"""

import re
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import quote

import numpy as np
import pandas as pd

from sief_druco.burden import lookup
from sief_druco.burden.lookup import PUG, PUG_VIEW, RXNAV, Cache

REPO_DIR = Path(__file__).resolve().parents[3]
OVERRIDES_PATH = REPO_DIR / "data" / "decagon_overrides.csv"

PUBCHEM_CACHE = "pubchem_decagon.json"
RXNAV_CACHE = "rxnav_decagon.json"

# Alias relations allowed by the alias rule of §4 (#26), and the tag of an extra
# `IN` the automatic route reached that no hand fix has ruled on yet.
RELATIONS = ("acid/salt", "prodrug", "racemate-enantiomer")
UNCHECKED = "unchecked"

VERDICTS = ("matched", "dropped_in", "unmatched")
OVERRIDE_COLUMNS = ["stitch_id", "main_in", "aliases", "verdict", "reason", "evidence"]

# match_method, in the order the steps are tried (pilot §3.1 steps 1-5).
EXACT, PARENT, SAME_CONNECTIVITY, SAME_CONNECTIVITY_SMILES, NAME, APPROXIMATE = (
    "exact_cid",
    "parent_cid",
    "same_connectivity",
    "same_connectivity_smiles",
    "name",
    "approximate_name_confirmed",
)

AUTO_COLUMNS = ["stitch_id", "cid", "rxcui", "name", "match_method", "name_confirmed", "reason"]

# A bare UNII among PubChem synonyms: 10 characters, at least one digit and one letter.
_BARE_UNII = re.compile(r"(?=.*[0-9])(?=.*[A-Z])[A-Z0-9]{10}")
_CHUNK = 100


def stitch_cid(stitch_id: str) -> int:
    """PubChem CID of a STITCH ID: strip `CID` and the leading zeros."""
    if not re.fullmatch(r"CID\d{9}", stitch_id):
        raise ValueError(f"not a STITCH ID: {stitch_id!r}")
    return int(stitch_id[3:])


def labelled_pair_counts(decagon) -> pd.Series:
    """Labelled pairs per drug (pairs with any of the 963 side effects), by STITCH ID."""
    labelled = np.asarray(decagon.labels.sum(axis=1)).ravel() > 0
    counts = np.bincount(decagon.pairs[labelled].ravel(), minlength=len(decagon.drug_ids))
    return pd.Series(counts, index=pd.Index(decagon.drug_ids, name="stitch_id"), name="labelled_pairs")


# --- 1. Automatic table ---


class _Route:
    """The lookups of one run; PubChem synonyms are kept across drugs."""

    def __init__(self, pubchem: Cache, rxnav: Cache):
        self.pc, self.rx = pubchem, rxnav
        self.synonyms: dict[int, list[str]] = {}

    # PubChem

    def fetch_synonyms(self, cids) -> None:
        cids = sorted(set(cids))
        for i in range(0, len(cids), _CHUNK):
            chunk = ",".join(map(str, cids[i : i + _CHUNK]))
            answer = self.pc.get(f"{PUG}/compound/cid/synonyms/JSON", {"cid": chunk})
            for info in (answer.get("InformationList") or {}).get("Information") or []:
                self.synonyms[info["CID"]] = info.get("Synonym") or []

    def titles(self, cids) -> dict[int, str]:
        out, cids = {}, sorted(set(cids))
        for i in range(0, len(cids), _CHUNK):
            chunk = ",".join(map(str, cids[i : i + _CHUNK]))
            answer = self.pc.get(f"{PUG}/compound/cid/property/Title/JSON", {"cid": chunk})
            for p in (answer.get("PropertyTable") or {}).get("Properties") or []:
                out[p["CID"]] = p.get("Title", "")
        return out

    def related_cids(self, cid: int, kind: str) -> list[int]:
        answer = self.pc.get(f"{PUG}/compound/cid/{cid}/cids/JSON?cids_type={kind}")
        return (answer.get("IdentifierList") or {}).get("CID") or []

    def same_connectivity_by_smiles(self, cid: int) -> list[int]:
        answer = self.pc.get(f"{PUG}/compound/cid/{cid}/property/IsomericSMILES/JSON")
        props = ((answer.get("PropertyTable") or {}).get("Properties") or [{}])[0]
        smiles = props.get("SMILES") or props.get("IsomericSMILES")
        if not smiles:
            return []
        answer = self.pc.get(
            f"{PUG}/compound/fastidentity/smiles/cids/JSON?identity_type=same_connectivity",
            {"smiles": smiles},
        )
        return (answer.get("IdentifierList") or {}).get("CID") or []

    def pugview_uniis(self, cid: int) -> set[str]:
        answer = self.pc.get(f"{PUG_VIEW}/data/compound/{cid}/JSON?heading=UNII")
        found = set()

        def walk(x):
            if isinstance(x, dict):
                if isinstance(x.get("String"), str):
                    found.add(x["String"].strip())
                for v in x.values():
                    walk(v)
            elif isinstance(x, list):
                for v in x:
                    walk(v)

        walk(answer.get("Record", {}))
        return found

    # RxNav

    def unii_ingredients(self, unii: str) -> list[tuple[str, str]]:
        answer = self.rx.get(f"{RXNAV}/rxcui.json?idtype=UNII_CODE&id={unii}")
        return [hit for rc in lookup.rx_ids(answer) for hit in lookup.rx_ingredients(self.rx, rc)]

    def name_ingredients(self, name: str) -> list[tuple[str, str]]:
        answer = self.rx.get(f"{RXNAV}/rxcui.json?name={quote(name)}&search=2")
        return [hit for rc in lookup.rx_ids(answer) for hit in lookup.rx_ingredients(self.rx, rc)]

    def approximate_ingredients(self, name: str) -> list[tuple[str, str]]:
        answer = self.rx.get(f"{RXNAV}/approximateTerm.json?term={quote(name)}&maxEntries=1&option=1")
        candidates = (answer.get("approximateGroup") or {}).get("candidate") or []
        return [
            hit
            for c in candidates
            if str(c.get("rank")) == "1"
            for hit in lookup.rx_ingredients(self.rx, c["rxcui"])
        ]

    # Steps

    def via_uniis(self, cids) -> tuple[list[tuple[str, str]], list[str]]:
        """`IN`s reached from the UNIIs of `cids`, and every UNII seen.

        `UNII-…` synonyms are taken as they are. A bare 10-character code reaching
        RxNorm counts only if PUG-View lists it as that CID's UNII (pilot §3.1).
        """
        hits, seen = [], []
        for cid in cids:
            syns = self.synonyms.get(cid, [])
            prefixed = sorted({s[5:] for s in syns if s.startswith("UNII-")})
            for unii in prefixed:
                seen.append(unii)
                hits += self.unii_ingredients(unii)
            for unii in sorted({s for s in syns if _BARE_UNII.fullmatch(s)} - set(prefixed)):
                found = self.unii_ingredients(unii)
                if found and unii in self.pugview_uniis(cid):
                    seen.append(unii)
                    hits += found
        return hits, seen

    def match(self, cid: int, title: str) -> tuple[list[tuple[str, str]], str, str, set[str]]:
        """`(IN hits, match_method, reason if unmatched, the drug's own names)`."""
        own = [title, *self.synonyms.get(cid, [])]
        hits, seen = self.via_uniis([cid])
        if hits:
            return hits, EXACT, "", _lower(own)

        parents = self.related_cids(cid, "parent")
        parent = parents[0] if parents else None
        if parent is not None and parent != cid:
            self.fetch_synonyms([parent])
            own += self.synonyms.get(parent, [])
            hits, more = self.via_uniis([parent])
            seen += more
            if hits:
                return hits, PARENT, "", _lower(own)

        same = [x for x in self.related_cids(parent or cid, "same_connectivity") if x not in (cid, parent)]
        if same:
            self.fetch_synonyms(same)
            hits, more = self.via_uniis(same)
            seen += more
            if hits:
                return hits, SAME_CONNECTIVITY, "", _lower(own)

        # Records with a structure but no parent and no relations in PubChem (#16).
        if parent is None:
            same = [x for x in self.same_connectivity_by_smiles(cid) if x != cid]
            if same:
                self.fetch_synonyms(same)
                hits, more = self.via_uniis(same)
                seen += more
                if hits:
                    return hits, SAME_CONNECTIVITY_SMILES, "", _lower(own)

        if title:
            hits = self.name_ingredients(title)
            if hits:
                return hits, NAME, "", _lower(own)
            # Rank-1 approximate hits on IUPAC titles are often wrong: keep one only
            # if the drug's own title or synonyms name the `IN`.
            approx = self.approximate_ingredients(title)
            confirmed = [h for h in approx if h[1].lower() in _lower([title, *self.synonyms.get(cid, [])])]
            if confirmed:
                return confirmed, APPROXIMATE, "", _lower(own)
            if approx:
                names = ", ".join(sorted({n for _, n in approx}))
                return [], "", f"approximateTerm hit not confirmed by name ({names})", _lower(own)
        if seen:
            uniis = ", ".join(sorted(set(seen)))
            return [], "", f"UNII found ({uniis}) but no RxNorm IN; title not found by name", _lower(own)
        reason = "no UNII on the exact, parent or same-connectivity CIDs; title not found by name"
        return [], "", reason, _lower(own)


def _lower(names) -> set[str]:
    return {n.lower() for n in names if n}


def automatic_table(stitch_ids, pubchem: Cache, rxnav: Cache) -> pd.DataFrame:
    """One row per drug × `IN` reached; an unmatched drug gets one row with a `reason`.

    `name_confirmed`: the `IN`'s name is among the drug's own PubChem title and
    synonyms (exact CID, and parent CID when that was looked up).
    """
    ids = sorted(stitch_ids)
    cid_of = {s: stitch_cid(s) for s in ids}
    route = _Route(pubchem, rxnav)
    title = route.titles(cid_of.values())
    route.fetch_synonyms(cid_of.values())
    rows = []
    for s in ids:
        cid = cid_of[s]
        hits, method, reason, own = route.match(cid, title.get(cid, ""))
        ins = dict(sorted(dict(hits).items(), key=lambda kv: int(kv[0])))
        if not ins:
            rows.append([s, cid, None, None, None, False, reason])
        for rxcui, name in ins.items():
            rows.append([s, cid, rxcui, name, method, name.lower() in own, None])
    return pd.DataFrame(rows, columns=AUTO_COLUMNS)


# --- 2. Overrides ---


def _parse_aliases(text: str, stitch_id: str) -> list[tuple[str, str]]:
    out = []
    for item in filter(None, (x.strip() for x in text.split(";"))):
        rxcui, sep, relation = (x.strip() for x in item.partition(":"))
        if not sep or not rxcui.isdigit():
            raise ValueError(f"{stitch_id}: alias {item!r} is not '<rxcui>:<relation>'")
        if relation not in RELATIONS:
            raise ValueError(
                f"{stitch_id}: alias {rxcui} has unknown relation {relation!r}; use one of {RELATIONS}"
            )
        out.append((rxcui, relation))
    return out


def load_overrides(stitch_ids, path: Path | None = OVERRIDES_PATH) -> pd.DataFrame:
    """The hand fixes, checked; `aliases` parsed into `[(rxcui, relation), ...]`.

    File format (`data/decagon_overrides.csv`): `main_in` is an RxCUI; `aliases` is
    `<rxcui>:<relation>` items joined by `;`, each relation one of `RELATIONS`;
    `verdict` is `matched`, `dropped_in` (some reached `IN`s dropped) or
    `unmatched` (no `main_in`, no aliases, a `reason`). `path` None, or a missing
    file, means no overrides.
    """
    if path is None or not Path(path).exists():
        return pd.DataFrame(columns=OVERRIDE_COLUMNS)
    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = set(OVERRIDE_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"{path.name} lacks columns {sorted(missing)}")
    df = df[OVERRIDE_COLUMNS].apply(lambda col: col.str.strip())

    known = set(stitch_ids)
    outside = sorted(set(df["stitch_id"]) - known)
    if outside:
        raise ValueError(f"overrides for STITCH IDs outside the {len(known)} Decagon drugs: {outside}")
    repeated = sorted(df["stitch_id"][df["stitch_id"].duplicated()])
    if repeated:
        raise ValueError(f"STITCH IDs with more than one override row: {repeated}")

    parsed = []
    for row in df.itertuples(index=False):
        s = row.stitch_id
        if row.verdict not in VERDICTS:
            raise ValueError(f"{s}: unknown verdict {row.verdict!r}; use one of {VERDICTS}")
        aliases = _parse_aliases(row.aliases, s)
        if row.verdict == "unmatched":
            if row.main_in or aliases:
                raise ValueError(f"{s}: an unmatched drug has no main_in and no aliases")
            if not row.reason:
                raise ValueError(f"{s}: an unmatched drug needs a reason")
        else:
            if not row.main_in.isdigit():
                raise ValueError(f"{s}: verdict {row.verdict} needs main_in as an RxCUI, got {row.main_in!r}")
            alias_ids = [a for a, _ in aliases]
            if row.main_in in alias_ids or len(set(alias_ids)) != len(alias_ids):
                raise ValueError(f"{s}: an IN appears twice among main_in and aliases")
        parsed.append(row._replace(aliases=aliases))
    return pd.DataFrame(parsed, columns=OVERRIDE_COLUMNS)


# --- 3. Resolution ---


@dataclass
class DrugTable:
    automatic: pd.DataFrame  # automatic_table(), before overrides
    # One row per drug × IN: stitch_id, rxcui, name, role ("main" / "alias"),
    # relation (None for the main IN; else one of RELATIONS or UNCHECKED),
    # source ("automatic" / "override").
    ingredients: pd.DataFrame
    unmatched: pd.DataFrame  # stitch_id, reason, source
    stitch_of: dict[str, str]  # IN rxcui -> the STITCH ID it joins
    # INs claimed by more than one drug: rxcui, name, stitch_id (chosen),
    # labelled_pairs, others ("<stitch_id> (<pairs>)" joined by "; "), for footnotes.
    collisions: pd.DataFrame


def _override_rows(o, rx: Cache) -> list[list]:
    rows = []
    for rxcui, role, relation in [(o.main_in, "main", None)] + [(a, "alias", r) for a, r in o.aliases]:
        props = lookup.rx_properties(rx, rxcui)
        if props.get("tty") != "IN":
            raise ValueError(
                f"{o.stitch_id}: override RxCUI {rxcui} is not an RxNorm IN ({props or 'unknown'})"
            )
        rows.append([o.stitch_id, rxcui, props["name"], role, relation, "override"])
    return rows


def _automatic_rows(group: pd.DataFrame) -> list[list]:
    group = group.sort_values("rxcui", key=lambda c: c.astype(int))
    confirmed = group[group["name_confirmed"]]
    main = (confirmed if len(confirmed) else group).iloc[0]["rxcui"]
    return [
        [r.stitch_id, r.rxcui, r.name, "main", None, "automatic"]
        if r.rxcui == main
        else [r.stitch_id, r.rxcui, r.name, "alias", UNCHECKED, "automatic"]
        for r in group.itertuples(index=False)
    ]


def resolve(
    automatic: pd.DataFrame, overrides: pd.DataFrame, pair_counts: pd.Series, rxnav: Cache
) -> DrugTable:
    """Apply the overrides and decide which STITCH ID each `IN` joins.

    `pair_counts`: labelled pairs per STITCH ID (`labelled_pair_counts`). An `IN`
    claimed by several drugs goes to the one with the most labelled pairs; a tie
    goes to the lower STITCH ID.
    """
    by_override = {o.stitch_id: o for o in overrides.itertuples(index=False)}
    rows, unmatched = [], []
    for s, group in automatic.groupby("stitch_id", sort=True):
        o = by_override.get(s)
        if o is not None and o.verdict == "unmatched":
            unmatched.append([s, o.reason, "override"])
        elif o is not None:
            rows += _override_rows(o, rxnav)
        elif group["rxcui"].isna().all():
            unmatched.append([s, group["reason"].iloc[0], "automatic"])
        else:
            rows += _automatic_rows(group)
    ingredients = pd.DataFrame(rows, columns=["stitch_id", "rxcui", "name", "role", "relation", "source"])
    unmatched = pd.DataFrame(unmatched, columns=["stitch_id", "reason", "source"])

    pairs = ingredients["stitch_id"].map(pair_counts)
    if pairs.isna().any():
        missing = sorted(set(ingredients["stitch_id"][pairs.isna()]))
        raise KeyError(f"no labelled-pair count for {missing}")
    ranked = ingredients.assign(labelled_pairs=pairs.astype(int)).sort_values(
        ["rxcui", "labelled_pairs", "stitch_id"], ascending=[True, False, True]
    )
    stitch_of, collisions = {}, []
    for rxcui, claims in ranked.groupby("rxcui", sort=True):
        chosen = claims.iloc[0]
        stitch_of[rxcui] = chosen["stitch_id"]
        if len(claims) > 1:
            others = "; ".join(f"{r.stitch_id} ({r.labelled_pairs})" for r in claims.iloc[1:].itertuples())
            collisions.append([rxcui, chosen["name"], chosen["stitch_id"], chosen["labelled_pairs"], others])
    collisions = pd.DataFrame(collisions, columns=["rxcui", "name", "stitch_id", "labelled_pairs", "others"])
    return DrugTable(automatic, ingredients, unmatched, stitch_of, collisions)


def build(
    stitch_ids=None,
    pair_counts: pd.Series | None = None,
    overrides_path: Path | None = OVERRIDES_PATH,
    cache_dir: Path = lookup.CACHE_DIR,
    fetch: lookup.Fetch = lookup.http_fetch,
) -> DrugTable:
    """The Decagon drug table; STITCH IDs and pair counts default to Decagon's own."""
    if stitch_ids is None or pair_counts is None:
        from sief_druco import decagon

        d = decagon.load()
        stitch_ids = d.drug_ids if stitch_ids is None else stitch_ids
        pair_counts = labelled_pair_counts(d) if pair_counts is None else pair_counts
    stitch_ids = [str(s) for s in stitch_ids]
    pc = Cache(Path(cache_dir) / PUBCHEM_CACHE, lookup.PUBCHEM_INTERVAL, fetch)
    rx = Cache(Path(cache_dir) / RXNAV_CACHE, lookup.RXNAV_INTERVAL, fetch)
    try:
        automatic = automatic_table(stitch_ids, pc, rx)
        overrides = load_overrides(stitch_ids, overrides_path)
        return resolve(automatic, overrides, pair_counts, rx)
    finally:
        pc.save()
        rx.save()
