"""The Decagon dataset: parsed, filtered and indexed once, then cached in `local/`.

Spec: `docs/spec-model.md` §1 (pair order of the hand-off) and §2 (Data).
Facts about the files: `research/decagon-files.md`.

Every ID set is sorted before numbering, so indices are deterministic. Pairs and
PPI edges are stored once, in one direction; reverse edges are added after
splitting, never here.
"""

import hashlib
from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

LOCAL_DIR = Path(__file__).resolve().parents[2] / "local"
RAW_DIR = LOCAL_DIR / "data" / "raw" / "decagon"
CACHE_PATH = LOCAL_DIR / "data" / "processed" / "decagon.npz"

FILES = {
    "combo": "bio-decagon-combo.csv",
    "mono": "bio-decagon-mono.csv",
    "ppi": "bio-decagon-ppi.csv",
    "targets": "bio-decagon-targets.csv",
}

N_DRUGS = 645
N_PAIRS = N_DRUGS * (N_DRUGS - 1) // 2
MIN_PAIRS_PER_SIDE_EFFECT = 500

# Bump when the parsed layout changes, so stale caches are rebuilt.
CACHE_VERSION = 1


def all_pairs(n: int = N_DRUGS) -> np.ndarray:
    """The canonical pair order: every `(a, b)` with `a < b`, row by row.

    Row i is `pair_id` i. Every per-pair array and the hand-off files use this order.
    """
    a, b = np.triu_indices(n, k=1)
    return np.stack([a, b], axis=1).astype(np.int64)


def pair_id(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row of the unordered pair `{a, b}` in `all_pairs()`; drug order doesn't matter."""
    a, b = np.minimum(a, b).astype(np.int64), np.maximum(a, b).astype(np.int64)
    if np.any(a == b):
        raise ValueError("a drug pair needs two distinct drugs")
    if np.any(a < 0) or np.any(b >= N_DRUGS):
        raise ValueError(f"drug indices must be in 0..{N_DRUGS - 1}")
    return a * (2 * N_DRUGS - a - 1) // 2 + (b - a - 1)


def pair_of(ids: np.ndarray) -> np.ndarray:
    """`(a, b)` rows for the given pair IDs, `a < b`."""
    ids = np.asarray(ids)
    if np.any(ids < 0) or np.any(ids >= N_PAIRS):
        raise ValueError(f"pair IDs must be in 0..{N_PAIRS - 1}")
    return all_pairs()[ids]


def verify_checksums(raw_dir: Path, names: list[str]) -> None:
    """Check each file against `SHA256SUMS` in `raw_dir`; raise on any mismatch."""
    expected = {}
    for line in (raw_dir / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split()
        expected[name] = digest
    for name in names:
        h = hashlib.sha256()
        with open(raw_dir / name, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        if h.hexdigest() != expected.get(name):
            raise ValueError(f"checksum mismatch for {name} in {raw_dir}")


@dataclass
class Decagon:
    # Index → ID (sorted); position in the array is the index.
    drug_ids: np.ndarray  # STITCH IDs, [645]
    protein_ids: np.ndarray  # Entrez IDs, [19089]
    side_effect_ids: np.ndarray  # combo CUIs after the ≥ 500 filter, [963]
    side_effect_names: np.ndarray
    mono_ids: np.ndarray  # mono CUIs, [10184]; disjoint from side_effect_ids
    mono_names: np.ndarray

    pairs: np.ndarray  # [207690, 2], all_pairs(); row = pair_id
    labels: sp.csr_matrix  # 0/1 [207690, 963], rows in pair order
    mono: sp.csr_matrix  # 0/1 [645, 10184], drug × mono CUI
    targets: np.ndarray  # [18690, 2] (drug, protein)
    ppi: np.ndarray  # [715612, 2] (protein, protein), each edge once

    def drug_index(self, stitch_ids) -> np.ndarray:
        return _lookup(self.drug_ids, stitch_ids)

    def protein_index(self, entrez_ids) -> np.ndarray:
        return _lookup(self.protein_ids, entrez_ids)

    def side_effect_index(self, cuis) -> np.ndarray:
        return _lookup(self.side_effect_ids, cuis)


def _lookup(sorted_ids: np.ndarray, ids) -> np.ndarray:
    # No cast to sorted_ids.dtype: a fixed-width string cast would cut off extra characters.
    ids = np.asarray(ids)
    idx = np.searchsorted(sorted_ids, ids)
    idx = np.minimum(idx, len(sorted_ids) - 1)
    if not np.array_equal(sorted_ids[idx], ids):
        raise KeyError(f"unknown IDs: {ids[sorted_ids[idx] != ids][:5]}")
    return idx


def load(raw_dir: Path = RAW_DIR, cache_path: Path = CACHE_PATH) -> Decagon:
    """The parsed dataset, from the cache if present, else parsed from the raw CSVs."""
    if cache_path.exists():
        with np.load(cache_path) as z:
            if int(z["version"]) == CACHE_VERSION:
                return _from_arrays(z)
    data = parse(raw_dir)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache_path, version=CACHE_VERSION, **_to_arrays(data))
    return data


def parse(raw_dir: Path = RAW_DIR) -> Decagon:
    verify_checksums(raw_dir, list(FILES.values()))
    read = lambda key: pd.read_csv(raw_dir / FILES[key], dtype=str)
    combo, mono, ppi, targets = read("combo"), read("mono"), read("ppi"), read("targets")

    drug_ids = np.unique(np.concatenate([combo["STITCH 1"], combo["STITCH 2"]]))
    assert len(drug_ids) == N_DRUGS
    targets = targets[targets["STITCH"].isin(drug_ids)]

    # Side effects: keep those with ≥ 500 distinct unordered pairs.
    a = np.searchsorted(drug_ids, combo["STITCH 1"].to_numpy())
    b = np.searchsorted(drug_ids, combo["STITCH 2"].to_numpy())
    triples = pd.DataFrame(
        {"pair": pair_id(a, b), "cui": combo["Polypharmacy Side Effect"]}
    ).drop_duplicates()
    pairs_per_cui = triples.groupby("cui")["pair"].nunique()
    kept = pairs_per_cui[pairs_per_cui >= MIN_PAIRS_PER_SIDE_EFFECT].index
    side_effect_ids = np.sort(kept.to_numpy(dtype=str))
    triples = triples[triples["cui"].isin(side_effect_ids)]
    names = combo.drop_duplicates("Polypharmacy Side Effect").set_index(
        "Polypharmacy Side Effect"
    )["Side Effect Name"]

    labels = sp.csr_matrix(
        (
            np.ones(len(triples), dtype=np.uint8),
            (triples["pair"], np.searchsorted(side_effect_ids, triples["cui"])),
        ),
        shape=(N_PAIRS, len(side_effect_ids)),
    )

    # Proteins: every PPI protein plus target proteins missing from the PPI graph.
    ppi_ids = ppi[["Gene 1", "Gene 2"]].to_numpy().astype(np.int64)
    target_genes = targets["Gene"].to_numpy().astype(np.int64)
    protein_ids = np.unique(np.concatenate([ppi_ids.ravel(), target_genes]))

    mono = mono.drop_duplicates(["STITCH", "Individual Side Effect"])
    mono_ids = np.unique(mono["Individual Side Effect"].to_numpy(dtype=str))
    mono_names = mono.drop_duplicates("Individual Side Effect").set_index(
        "Individual Side Effect"
    )["Side Effect Name"]
    mono_matrix = sp.csr_matrix(
        (
            np.ones(len(mono), dtype=np.uint8),
            (
                np.searchsorted(drug_ids, mono["STITCH"]),
                np.searchsorted(mono_ids, mono["Individual Side Effect"]),
            ),
        ),
        shape=(N_DRUGS, len(mono_ids)),
    )

    return Decagon(
        drug_ids=drug_ids.astype(str),
        protein_ids=protein_ids,
        side_effect_ids=side_effect_ids,
        side_effect_names=names[side_effect_ids].to_numpy(dtype=str),
        mono_ids=mono_ids,
        mono_names=mono_names[mono_ids].to_numpy(dtype=str),
        pairs=all_pairs(),
        labels=labels,
        mono=mono_matrix,
        targets=np.stack(
            [
                np.searchsorted(drug_ids, targets["STITCH"]),
                np.searchsorted(protein_ids, target_genes),
            ],
            axis=1,
        ),
        ppi=np.searchsorted(protein_ids, ppi_ids),
    )


# --- Cache layout: dense arrays as-is, sparse matrices as their CSR parts. ---


def _to_arrays(data: Decagon) -> dict[str, np.ndarray]:
    out = {}
    for f in fields(data):
        x = getattr(data, f.name)
        if sp.issparse(x):
            out |= {
                f"{f.name}.data": x.data,
                f"{f.name}.indices": x.indices,
                f"{f.name}.indptr": x.indptr,
                f"{f.name}.shape": np.array(x.shape),
            }
        else:
            out[f.name] = x
    return out


def _from_arrays(z) -> Decagon:
    kwargs = {}
    for f in fields(Decagon):
        if f.name in z:
            kwargs[f.name] = z[f.name]
        else:
            parts = [z[f"{f.name}.{p}"] for p in ("data", "indices", "indptr")]
            kwargs[f.name] = sp.csr_matrix(tuple(parts), shape=tuple(z[f"{f.name}.shape"]))
    return Decagon(**kwargs)
