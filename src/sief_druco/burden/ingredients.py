"""Source drugs → RxNorm ingredients (`IN`).

Spec: `docs/spec-burden.md` §4 (step 2 and the rules). The route is the pilot's
(`research/nhanes-coverage-pilot.md` §3.2, §5).

- **Multum codes** (NHANES `RXDDRGID`): `GN` + code as an `MMSL_CODE`, then the
  same with `allsrc=1`; if that reaches no `IN`, the `RXDDRUG` name split at `;`.
- **Names** (MEPS `RXDRGNAM`, DDInter): the whole name, then its parts split at the
  caller's separators.

Every RxCUI goes to its `IN`s through `related.json?tty=IN`, so a PIN lifts to its
`IN` and a multi-ingredient product splits into its `IN`s. An RxCUI with no `IN` (a
retired concept) is followed through its `historystatus` remap. A name split into
parts reaches its `IN`s only if every part does.

No source-side hand fixes (§4): a drug that fails keeps one of the four `CAUSES`.
"""

from dataclasses import dataclass
from urllib.parse import quote

from sief_druco.burden import lookup
from sief_druco.burden.lookup import RXNAV, Cache

# Why a source drug reached no IN.
NO_CODE = "no_code"  # NHANES 55555 / 77777 / 99999: unknown, refused, don't know
ACH_CODE = "a_c_h_code"  # not a Multum drug code: a / h (not in Multum), c (a class)
RETIRED_CODE = "retired_code"  # RxNorm knows the code or name, but only as a retired concept
NAME_NOT_FOUND = "name_not_found"  # neither the code nor the name (or one of its parts) found
CAUSES = (NO_CODE, ACH_CODE, RETIRED_CODE, NAME_NOT_FOUND)

UNCODED_NAMES = frozenset({"55555", "77777", "99999"})

# How the INs were reached.
MMSL, MMSL_ALLSRC, NAME = "mmsl_code", "mmsl_code_allsrc", "name"


@dataclass(frozen=True)
class Resolution:
    ingredients: tuple[tuple[str, str], ...]  # (IN rxcui, IN name), by rxcui
    cause: str | None  # one of CAUSES when no IN was reached, else None
    step: str  # MMSL, MMSL_ALLSRC or NAME; "" when no IN was reached


def _sorted(hits) -> tuple[tuple[str, str], ...]:
    return tuple(sorted(dict(hits).items(), key=lambda kv: int(kv[0])))


class Resolver:
    """Resolves source drugs over one disk-cached RxNav (`lookup.Cache`)."""

    def __init__(self, rxnav: Cache):
        self.rx = rxnav

    def to_ingredients(self, rxcui: str) -> list[tuple[str, str]]:
        """`IN`s of an RxCUI; a retired one goes through its remapped concepts."""
        hits = lookup.rx_ingredients(self.rx, rxcui)
        if hits:
            return hits
        answer = self.rx.get(f"{RXNAV}/rxcui/{rxcui}/historystatus.json")
        derived = (answer.get("rxcuiStatusHistory") or {}).get("derivedConcepts") or {}
        return [
            hit
            for m in derived.get("remappedConcept") or []
            for hit in lookup.rx_ingredients(self.rx, m["remappedRxCui"])
        ]

    def _name_rxcuis(self, name: str) -> list[str]:
        return lookup.rx_ids(self.rx.get(f"{RXNAV}/rxcui.json?name={quote(name)}&search=2"))

    def _parts(self, parts: list[str]) -> list[tuple[str, str]]:
        """`IN` hits of all `parts`, or none unless every part reached an `IN`."""
        hits = []
        for part in parts:
            found = [h for rc in self._name_rxcuis(part) for h in self.to_ingredients(rc)]
            if not found:
                return []
            hits += found
        return hits

    def multum(self, code: str, name: str) -> Resolution:
        """An NHANES drug row: `RXDDRGID` and `RXDDRUG`."""
        code, name = (code or "").strip(), (name or "").strip()
        if not code or name in UNCODED_NAMES:
            return Resolution((), NO_CODE, "")
        if not code.startswith("d"):
            return Resolution((), ACH_CODE, "")

        url = f"{RXNAV}/rxcui.json?idtype=MMSL_CODE&id=GN{code}"
        rxcuis, step = lookup.rx_ids(self.rx.get(url)), MMSL
        if not rxcuis:
            rxcuis, step = lookup.rx_ids(self.rx.get(url + "&allsrc=1")), MMSL_ALLSRC
        hits = [h for rc in rxcuis for h in self.to_ingredients(rc)]
        if hits:
            return Resolution(_sorted(hits), None, step)

        hits = self._parts([p.strip() for p in name.split(";") if p.strip()])
        if hits:
            return Resolution(_sorted(hits), None, NAME)
        return Resolution((), RETIRED_CODE if rxcuis else NAME_NOT_FOUND, "")

    def name(self, name: str, separators: tuple[str, ...] = ()) -> Resolution:
        """A drug name: whole, then split at every one of `separators` (e.g. `("-",)`).

        `RETIRED_CODE` if the whole name is known to RxNorm only as retired concepts;
        else `NAME_NOT_FOUND`.
        """
        name = (name or "").strip()
        if not name:
            return Resolution((), NAME_NOT_FOUND, "")
        whole = self._name_rxcuis(name)
        hits = [h for rc in whole for h in self.to_ingredients(rc)]
        if not hits and separators:
            parts = [name]
            for sep in separators:
                parts = [q.strip() for p in parts for q in p.split(sep)]
            parts = [p for p in parts if p]
            if len(parts) > 1:
                hits = self._parts(parts)
        if hits:
            return Resolution(_sorted(hits), None, NAME)
        return Resolution((), RETIRED_CODE if whole else NAME_NOT_FOUND, "")
