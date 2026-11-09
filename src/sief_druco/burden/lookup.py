"""Disk-cached RxNav and PubChem lookups.

Spec: `docs/spec-burden.md` §4 ("Cache every RxNav / PubChem lookup to disk").

One JSON file per cache under `local/data/cache/`, keyed by the request (the URL,
plus the form data of a POST). Every answer is stored, HTTP errors too (as
`{"_status": code}`), so a second run with a full cache sends no request. A
request that keeps failing (no connection, or 429 / 5xx on every attempt) raises
`LookupFailed` and is not stored; what was fetched so far is kept, and the next
run retries it.

`fetch` is injectable, so tests run without the network.
"""

import json
import logging
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Callable
from pathlib import Path

LOCAL_DIR = Path(__file__).resolve().parents[3] / "local"
CACHE_DIR = LOCAL_DIR / "data" / "cache"

RXNAV = "https://rxnav.nlm.nih.gov/REST"
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUG_VIEW = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view"

# Seconds between requests: PubChem allows 5 per second, RxNav 20.
PUBCHEM_INTERVAL = 0.25
RXNAV_INTERVAL = 0.06

RETRY_STATUS = {429, 500, 502, 503, 504}
ATTEMPTS = 5
SAVE_EVERY = 50
USER_AGENT = "sief-druco (research; github.com/MichaelRos1/sief-druco)"

log = logging.getLogger(__name__)

# fetch(url, form) -> (HTTP status, body). `form` is None for a GET.
Fetch = Callable[[str, dict | None], tuple[int, str]]


class LookupFailed(RuntimeError):
    """A request got no usable answer after every attempt."""


def http_fetch(url: str, form: dict | None) -> tuple[int, str]:
    data = None if form is None else urllib.parse.urlencode(form).encode()
    req = urllib.request.Request(url, data=data, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.status, r.read().decode()
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode(errors="replace")


class Cache:
    """Answers from `path` if present, else from `fetch`, rate-limited and stored."""

    def __init__(self, path: Path, min_interval: float = 0.0, fetch: Fetch = http_fetch):
        self.path = Path(path)
        self.data = json.loads(self.path.read_text()) if self.path.exists() else {}
        self.min_interval = min_interval
        self.fetch = fetch
        self._last = 0.0
        self._unsaved = 0

    @staticmethod
    def key(url: str, form: dict | None = None) -> str:
        return url if form is None else url + " POST " + json.dumps(form, sort_keys=True)

    def get(self, url: str, form: dict | None = None) -> dict:
        key = self.key(url, form)
        if key in self.data:
            return self.data[key]
        for attempt in range(ATTEMPTS):
            wait = self.min_interval - (time.monotonic() - self._last)
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()
            try:
                status, body = self.fetch(url, form)
            except OSError as e:  # URLError, timeouts, dropped connections
                log.warning("%s: %r (attempt %d)", key, e, attempt + 1)
                time.sleep(2 * (attempt + 1))
                continue
            if status in RETRY_STATUS:
                log.warning("%s: HTTP %d (attempt %d)", key, status, attempt + 1)
                time.sleep(3 * (attempt + 1))
                continue
            if status == 200:
                try:
                    value = json.loads(body)
                except ValueError:
                    value = {"_text": body}
            else:
                value = {"_status": status, "_text": body[:300]}
            self.data[key] = value
            self._unsaved += 1
            if self._unsaved >= SAVE_EVERY:
                self.save()
            return value
        # Raise rather than return an empty answer: the route would read it as
        # "nothing found" and list the drug as unmatched with a false reason.
        raise LookupFailed(f"{key}: no answer after {ATTEMPTS} attempts; rerun to retry")

    def save(self) -> None:
        if not self._unsaved and self.path.exists():
            return
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp = self.path.with_suffix(".tmp")
        tmp.write_text(json.dumps(self.data))
        tmp.replace(self.path)
        self._unsaved = 0


# --- RxNav answers ---


def rx_ids(answer: dict) -> list[str]:
    """RxCUIs of an `rxcui.json` answer."""
    return (answer.get("idGroup") or {}).get("rxnormId") or []


def rx_ingredients(rx: Cache, rxcui: str) -> list[tuple[str, str]]:
    """`(IN rxcui, IN name)` for an RxCUI, via `related.json?tty=IN`."""
    answer = rx.get(f"{RXNAV}/rxcui/{rxcui}/related.json?tty=IN")
    return [
        (p["rxcui"], p["name"])
        for g in (answer.get("relatedGroup") or {}).get("conceptGroup") or []
        for p in g.get("conceptProperties") or []
    ]


def rx_properties(rx: Cache, rxcui: str) -> dict:
    """`{"rxcui", "name", "tty", ...}` of an RxCUI, or `{}` if RxNav doesn't know it."""
    answer = rx.get(f"{RXNAV}/rxcui/{rxcui}/properties.json")
    return answer.get("properties") or {}
